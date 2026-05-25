library(tidyverse)
library(plotly)
library(htmlwidgets)

coefs <- read_csv(
  file.path("data", "ManuscriptSummaries", "beta_coefficients.csv"),
  show_col_types = FALSE
) |>
  janitor::clean_names()

# --------------------------------------------------------------------------- #
# Variable metadata
# --------------------------------------------------------------------------- #

bio_meta <- tibble::tribble(
  ~variable,  ~group,           ~label,
  "bio_01",   "Temperature",    "Mean Annual Temp",
  "bio_02",   "Variability",    "Mean Diurnal Range",
  "bio_03",   "Variability",    "Isothermality",
  "bio_04",   "Variability",    "Temp Seasonality",
  "bio_05",   "Temperature",    "Max Temp Warmest Mo.",
  "bio_06",   "Cold",           "Min Temp Coldest Mo.",
  "bio_07",   "Variability",    "Temp Annual Range",
  "bio_08",   "Temperature",    "Mean Temp Wettest Qtr",
  "bio_10",   "Temperature",    "Mean Temp Warmest Qtr",
  "bio_11",   "Cold",           "Mean Temp Coldest Qtr",
  "bio_12",   "Precipitation",  "Annual Precip",
  "bio_15",   "Precipitation",  "Precip Seasonality",
  "bio_16",   "Precipitation",  "Precip Wettest Qtr",
  "bio_17",   "Precipitation",  "Precip Driest Qtr",
  "bio_18",   "Precipitation",  "Precip Warmest Qtr"
)

# Expected sign relationships for the continental interior
# Each row: if var_a has sign X, var_b should have sign Y (same/opposite)
expected_pairs <- tibble::tribble(
  ~var_a,    ~var_b,    ~expected,   ~rationale,
  # same-sign groups
  "bio_01",  "bio_05",  "same",      "warmth metrics",
  "bio_01",  "bio_08",  "same",      "warmth metrics; spring wet season here",
  "bio_01",  "bio_10",  "same",      "warmth metrics",
  "bio_05",  "bio_08",  "same",      "warmth metrics",
  "bio_02",  "bio_04",  "same",      "temperature variability",
  "bio_02",  "bio_07",  "same",      "temperature variability",
  "bio_04",  "bio_07",  "same",      "temperature variability"
)

# --------------------------------------------------------------------------- #
# Classify signs and detect violations
# --------------------------------------------------------------------------- #

coefs_signed <- coefs |>
  left_join(bio_meta, by = "variable") |>
  mutate(
    sign_class = case_when(
      q2_5 > 0 & q97_5 > 0 ~ "positive",
      q2_5 < 0 & q97_5 < 0 ~ "negative",
      TRUE                  ~ "ambiguous"
    )
  )

# For each species, join pairs and check for violations
violations <- coefs_signed |>
  select(species, variable, sign_class) |>
  inner_join(expected_pairs, by = c("variable" = "var_a"), relationship = "many-to-many") |>
  inner_join(
    coefs_signed |> select(species, variable, sign_b = sign_class),
    by = c("species", "var_b" = "variable")
  ) |>
  filter(sign_class != "ambiguous", sign_b != "ambiguous") |>
  mutate(
    violated = case_when(
      expected == "same"     & sign_class != sign_b ~ TRUE,
      expected == "opposite" & sign_class == sign_b ~ TRUE,
      TRUE                                          ~ FALSE
    )
  )

# --------------------------------------------------------------------------- #
# Heatmap: species x variable, fill = estimate, mark ambiguous cells
# --------------------------------------------------------------------------- #

# Variable order: group then label
var_order <- bio_meta |>
  arrange(group, variable) |>
  pull(label)

sp_order <- sort(unique(coefs_signed$species))

# Sign-preserving sqrt compresses large outliers and opens up the near-zero
# region so the midpoint reads clearly; tooltip still shows raw β
plot_data <- coefs_signed |>
  mutate(
    label        = factor(label,   levels = var_order),
    species      = factor(species, levels = sp_order),
    estimate_t   = sign(estimate) * sqrt(abs(estimate)),
    tooltip      = glue::glue(
      "<b>{species}</b><br>
       {label}<br>
       β = {round(estimate, 3)} [{round(q2_5, 3)}, {round(q97_5, 3)}]<br>
       Sign: {sign_class}"
    )
  )

# Ambiguous cells for overlay
ambig <- plot_data |> filter(sign_class == "ambiguous")

# Violation cells for overlay — tooltip describes the conflicting pair
violation_cells <- violations |>
  filter(violated) |>
  rename(var_a = variable) |>
  left_join(bio_meta |> select(variable, label_a = label), by = c("var_a" = "variable")) |>
  left_join(bio_meta |> select(variable, label_b = label), by = c("var_b" = "variable")) |>
  mutate(
    tooltip_v = glue::glue(
      "<b>Unexpected combination</b><br>
       {label_a}: {sign_class}<br>
       {label_b}: {sign_b}<br>
       Expected: {expected} sign<br>
       ({rationale})"
    )
  ) |>
  rename(label = label_a) |>
  distinct(species, label, .keep_all = TRUE)

p <- plot_ly() |>
  add_heatmap(
    data        = plot_data,
    x           = ~species,
    y           = ~label,
    z           = ~estimate_t,
    text        = ~tooltip,
    hoverinfo   = "text",
    colorscale  = list(c(0, "#d6604d"), c(0.5, "white"), c(1, "#2166ac")),
    zmid        = 0,
    colorbar    = list(title = "√|β| (signed)"),
    xgap = 1, ygap = 1
  ) |>
  add_markers(
    data        = ambig,
    x           = ~species,
    y           = ~label,
    text        = ~tooltip,
    hoverinfo   = "text",
    marker      = list(symbol = "x", size = 5, color = "grey50"),
    showlegend  = FALSE
  ) |>
  add_markers(
    data        = violation_cells,
    x           = ~species,
    y           = ~label,
    text        = ~tooltip_v,
    hoverinfo   = "text",
    marker      = list(symbol = "circle-open", size = 12,
                       color  = "black", line = list(width = 2)),
    showlegend  = FALSE
  ) |>
  layout(
    title  = "Visual Display of Variables in Species Distribution Models",
    xaxis     = list(title = "", tickangle = -45, categoryorder = "array",
                     categoryarray = sp_order),
    yaxis     = list(title = "", categoryorder = "array",
                     categoryarray = rev(var_order)),
    margin    = list(l = 140, b = 160)
  )

saveWidget(p, "figures/coef_sign_explorer.html", selfcontained = TRUE)