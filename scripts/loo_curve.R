library(tidyverse)
library(plotly)
library(htmlwidgets)

loo_dat <- read_csv(
  file.path("data", "ManuscriptSummaries", "loo_roc.csv"),
  show_col_types = FALSE
) |>
  janitor::clean_names()

colnames(loo_dat)

fig <- plot_ly(
  loo_dat,
  x = ~fpr,
  y = ~tpr,
  split = ~species,
  type = 'scatter',
  mode = 'lines',
  text = ~species,
  hoverinfo = 'text',
  visible = "legendonly",
  line = list(color = 'steelblue')
) |>
  layout(
  title = 'Leave-One-Out ROC Curve',
    xaxis = list(
      title = "False Positive Rate", 
      range = c(0, 1)
    ),
    yaxis = list(
      title = "True Positive Rate", 
      range = c(0, 1),
      scaleanchor = "x", 
      scaleratio = 1
    )
  )


fig <- plot_ly(
  loo_dat,
  x = ~fpr,
  y = ~tpr,
  split = ~species,
  type = 'scatter',
  mode = 'lines',
  text = ~species,
  hoverinfo = 'text',
  line = list(color = 'steelblue')
) |>
  layout(
    paper_bgcolor = "rgba(0,0,0,0)",
    plot_bgcolor  = "rgba(0,0,0,0)",
    xaxis = list(
      title = "False Positive Rate",
      range = c(0, 1)#,
   #   gridcolor = "rgba(128,128,128,0.3)", 
   #   zerolinecolor = "rgba(128,128,128,0.5)"
    ),
    yaxis = list(
      title = "True Positive Rate", 
      range = c(0, 1)#,
    #  gridcolor = "rgba(128,128,128,0.3)", 
    #  zerolinecolor = "rgba(128,128,128,0.5)"
    ),
    font = list(color = "#ddd"),
    updatemenus = list(
      list(
        type = "buttons",
        direction = "right",
        x = 0.6, xanchor = "left",
        y = 0.15, yanchor = "top",
        buttons = list(
          list(method = "restyle", args = list("visible", TRUE),         label = "Show all"),
          list(method = "restyle", args = list("visible", "legendonly"), label = "Hide all")
        )
      )
    )
  )
fig

saveWidget(fig, "figures/loo_curve.html", selfcontained = TRUE)
