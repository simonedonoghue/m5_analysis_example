# In the follwoing script, I load the shapefiles for the 2011 census tracts and 
# M5 stops, clean and merge them. I subsequently compute the distance for 
# each tract to the nearest m5 stop and assign treatment/control status, as 
# well as the timing of treatment. The script yields three maps:
# MAP 1: Classification of tracts based on treatment status
# MAP 2: Classification of tracts based on distance from nearest M5 stopte
# MAP 3: Staggered timing of treatment for treated tracts

# BEGINNING OF SCRIPT

# 0. SET UP -------------------------------------------------------------------

# NB: remember to install missing packages
library(sf)
library(tidyverse)
library(tmap)
library(units)
library(here)

tmap_mode("plot")

# Create paths
data_path <- here("data")
output_path <- here("output")

# Set distance parameters
treat_radius      <- set_units(500, "m")
donut_inner       <- set_units(800, "m")
donut_outer       <- set_units(2000, "m")
study_area_radius <- set_units(3000, "m")

# Target projected CRS
target_crs <- 3035   # ETRS89 / LAEA Europe


# 1. LOAD AND PREPARE CENSUS TRACTS -------------------------------------------

lombardy <- st_read(
  here(data_path,
       "R03_11_WGS84 sez censimento lombardia",
       "R03_11_WGS84.shp"),
  quiet = TRUE
)

# Validate geometries
lombardy <- st_make_valid(lombardy)

# Create  identifiers
lombardy <- lombardy %>%
  mutate(
    tract_id   = as.character(SEZ2011),
    cod_comune = PRO_COM
  )

# Transform to projected CRS
lombardy <- st_transform(lombardy, target_crs)

stopifnot(st_is_valid(lombardy))


# 2. LOAD AND PREPARE M5 STATIONS ---------------------------------------------

metro <- st_read(
  here(data_path, "tpl_metrofermate", "tpl_metrofermate.shp"),
  quiet = TRUE
)

# Filter M5 stations
m5 <- metro %>%
  filter(str_detect(linee, "5")) %>%
  st_make_valid() %>%
  st_transform(target_crs)

stopifnot(st_crs(m5) == st_crs(lombardy))


# 3. COMPUTE DISTANCES --------------------------------------------------------

# Create centroids
centroids <- st_centroid(lombardy)

# Efficient nearest station computation
nearest_idx <- st_nearest_feature(centroids, m5)

lombardy$dist_m5 <- st_distance(
  centroids,
  m5[nearest_idx, ],
  by_element = TRUE
)

# Ensure units are metres
units(lombardy$dist_m5)

# 4. DEFINE STUDY AREA AND TREATMENT ------------------------------------------

# Study area (tracts within 3000m)
study_area <- lombardy %>%
  filter(dist_m5 <= study_area_radius)

# Treatment definition (T is 0-500m, C is 800-2000m)
study_area <- study_area %>%
  mutate(
    treated = dist_m5 <= treat_radius,
    
    # Donut control reduces spatial spillovers
    control_donut = dist_m5 > donut_inner &
      dist_m5 <= donut_outer,
    
    treatment_status = case_when(
      treated ~ "Treated (≤500m)",
      control_donut ~ "Control (800–2000m)",
      TRUE ~ "Excluded"
    )
  )

table(study_area$treatment_status)

# 5. MAP 1: TREATMENT AND CONTROL AREAS --------------------------------------

map_treatment <- 
  tm_shape(study_area) +
  tm_polygons(
    col = "treatment_status",
    palette = c(
      "Treated (≤500m)" = "darkred",
      "Control (800–2000m)" = "steelblue",
      "Excluded" = "grey85"
    ),
    title = "Treatment status"
  ) +
  tm_shape(m5) +
  tm_dots(size = 0.05, col = "black") +
  tm_layout(
    legend.outside = TRUE,
    frame = FALSE,
    main.title = "Spatial Definition of Treatment and Control Areas (M5)",
    main.title.size = 1.2
  )

map_treatment

# 6. ADD TREATMENT TIMING ----------------------------------------------------

# Opening year by station
m5_opening <- tibble(
  nome = c(
    "BIGNAMI", "PONALE", "BICOCCA", "CA GRANDA", "ISTRIA", "MARCHE", "ZARA",
    "ISOLA", "GARIBALDI FS", "MONUMENTALE",
    "CENISIO", "GERUSALEMME", "DOMODOSSOLA", "PORTELLO",
    "LOTTO M5", "SEGESTA", "SAN SIRO IPPODROMO",
    "SAN SIRO STADIO", "TRE TORRI"
  ),
  opening_year = c(
    2013, 2013, 2013, 2013, 2013, 2013, 2013,
    2014, 2014, 2014,
    2015, 2015, 2015, 2015,
    2015, 2015, 2015,
    2015, 2015
  )
)

m5 <- m5 %>%
  left_join(m5_opening, by = "nome")

# Recompute nearest station for study area only
centroids_sa <- st_centroid(study_area)

nearest_idx_sa <- st_nearest_feature(centroids_sa, m5)

study_area$first_treat_year <- m5$opening_year[nearest_idx_sa]

# Controls are never treated
study_area <- study_area %>%
  mutate(
    first_treat_year = ifelse(treated, first_treat_year, NA)
  )

table(study_area$first_treat_year, useNA = "ifany")

# 7. MAP 2: DISTANCE GRADIENT ------------------------------------------------

map_distance <- 
  tm_shape(study_area) +
  tm_polygons(
    col = "dist_m5",
    style = "quantile",
    palette = "-viridis",   # reversed (often visually softer)
    alpha = 0.7,            # transparency softens intensity
    title = "Distance to nearest M5 station (m)"
  ) +
  tm_shape(m5) +
  tm_dots(size = 0.05, col = "black") +
  tm_layout(
    legend.outside = TRUE,
    frame = FALSE,
    main.title = "Distance to M5 Stations",
    main.title.size = 1.2
  )

map_distance

# 8. MAP3: STAGGERED TREATMENT TIMING ----------------------------------------

map_timing <- 
  tm_shape(study_area) +
  tm_polygons(
    col = "first_treat_year",
    palette = "viridis",
    title = "Treatment year",
    colorNA = "grey90",
    textNA = "Never treated"
  ) +
  tm_shape(m5) +
  tm_dots(size = 0.05, col = "black") +
  tm_layout(
    legend.outside = TRUE,
    frame = FALSE,
    main.title = "Staggered Exposure to the M5 Metro Line",
    main.title.size = 1.2
  )

map_timing

# 9. EXPORT MAPS -------------------------------------------------------------

st_write(
  study_area,
  here(output_path, "study_area_m5_spatial_base.shp"),
  delete_dsn = TRUE,
  quiet = TRUE
)

st_write(
  study_area,
  here(output_path, "study_area_m5_spatial_timing.shp"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# END OF SCRIPT.
