library(sf)
library(tidyverse)
library(tigris)
library(rnaturalearth)

cook = counties(state = 'Illinois') |>
  filter(NAME == 'Cook') |>
  select(geometry)

canada = rnaturalearth::ne_states(country = 'Canada') |>
  select(geometry)
us_states = rnaturalearth::ne_states(country = 'United States of America')|>
  select(geometry)
states = bind_rows(canada, us_states) |>
  sf::st_simplify() |>
  st_make_valid()


bb_ll = rbind(
  cbind(rep(-170, 50), seq(24, 85, length.out = 50)),   # western edge
  cbind(rep(-48,  50), seq(24, 85, length.out = 50)),   # eastern edge
  cbind(seq(-170, -48, length.out = 50), rep(24, 50)),  # southern edge
  cbind(seq(-170, -48, length.out = 50), rep(85, 50))   # northern edge
) |>
  sf::st_multipoint() |>
  sf::st_sfc(crs = 4326) |>
  sf::st_convex_hull()

states = sf::st_intersection(states, bb_ll)

ggplot() + 
  geom_sf(data = states) + 
  geom_sf(data = bb_ll, fill = NA)

lakes <- rnaturalearth::ne_download(
  scale = 10, type = "lakes", category = "physical", returnclass = "sf"
) |>
  filter(name_alt == 'Great Lakes' | is.na(name_alt)) |>
  dplyr::filter(
          name %in% c('Lake Winnipeg', 'Lake Manitoba', 'Lake Winnipegosis') |
          name_alt == 'Great Lakes' 
      ) |>
  st_transform(st_crs(cook)) |>
  select(geometry) |>
  sf::st_make_valid()

cook <- sf::st_difference(cook, st_union(lakes))


### import cluster results - current state
lyr_names = sf::st_layers(file.path('data', 'all_cluster_geometries.gpkg'))

spec = lapply(lyr_names$name, \(x)
  sf::st_read(
    file.path('data', 'all_cluster_geometries.gpkg'),
     layer = x, 
     quiet = TRUE) |>
    dplyr::mutate(species = x)
) |>
  bind_rows() |>
  sf::st_make_valid() |>
  sf::st_transform("ESRI:102008")

cook = st_transform(cook, st_crs(spec)) |>
  sf::st_make_valid() |>
  st_union()

## total area for each species
sp_union = spec |>
  summarize(
    geometry = st_union(geom), 
    .by = species
  ) 

total_areas = sp_union |>
  sf::st_transform("ESRI:102008")

## total number of clusters per species

n_clusts = spec |>
  sf::st_drop_geometry()|>
  dplyr::summarise(
    n_presence = dplyr::n(),
    .by = 'species'
  )

## identify the cluster(s) per species that overlap cook county
ints = spec[ unlist(st_intersects(cook, spec)), ]

## area of the cook county seed zone cluster(s)
total_areas_ha = as.numeric(units::set_units(st_area(ints), 'ha'))

## area of the overlapping portions of each seed zone with cook county

overlaps = st_intersection(cook, spec)
as.numeric(units::set_units(st_area(overlaps), 'ha'))

## create domain for map... use MT/WY border to Maine, south to southern border of arkansas



## plot all cook county overlapping seed zones (very light grey. )


ggplot(data = spec) + 
  geom_sf()

