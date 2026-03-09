# BBS-syntopy project outline ---------------------------------------------
# Matěj Tvarůžka

# Here I focus on building a conceptual tool for my thesis,
# which will analyse local co-occurence (syntopy) on two spatial scales (between and withing transects) using data from Breeding Bird Survey.

### Dataset:
# routes length: 25.4 miles (= ~ 39.42893 km)
# stops distance: 0.5 mile (= ~ 804.672 m)
# 50 stops in each transect
# 3 minute point count at each stop
# record in radius of 0.25 mile (= ~ 402.336 m)
# number of transects: 4100+

# Installing bbsAssistant package from GitHub forked repository -----------------------------
# Original repository here: https://github.com/trashbirdecology/bbsAssistant

remotes::install_github("matejtvar/bbsAssistant")


# Libraries ---------------------------------------------------------------

library(bbsAssistant) # for downloading, importing, and munging the official releases of the North American Breeding Bird Survey (BBS) data
ls("package:bbsAssistant") # view functions and data in package bbsAssistant



# ============================================================
# 0. Libraries
# ============================================================

library(tidyverse)
library(sf)


# ============================================================
# 1. Load BBS data and subset to 2018
# ============================================================

bbs <- grab_bbs_data()

bbs_2018 <- bbs$observations %>%
  filter(Year == 2018)

# ============================================================
# 2. Prepare route metadata (points)
# ============================================================

routes_sf <- bbs$routes %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# ============================================================
# 3. Join observations to routes (attribute join)
# ============================================================

bbs_joined <- routes_sf %>%
  left_join(
    bbs_2018,
    by = c("RTENO", "CountryNum", "StateNum", "Route")
  )

# ============================================================
# 4. Compute species presence per route
# ============================================================

stop_cols <- paste0("Stop", 1:50)

bbs_presence <- bbs_joined %>%
  mutate(
    present = rowSums(across(all_of(stop_cols)), na.rm = TRUE) > 0,
    present = as.integer(present)
  ) %>%
  filter(present == 1)

# ============================================================
# 5. Compute species richness per route
# ============================================================

route_richness <- bbs_presence %>%
  group_by(RTENO) %>%
  summarise(
    species_richness = n_distinct(AOU),
    .groups = "drop"
  ) %>%
  st_drop_geometry()

# ============================================================
# 6. Create route metadata table (one row per route)
# ============================================================

route_meta <- bbs_presence %>%
  distinct(
    RTENO,
    RouteName,
    StateNum,
    Route,
    CountryNum,
    geometry
  )

# ============================================================
# 7. Join richness back to route metadata
# ============================================================

routes_final <- route_meta %>%
  left_join(route_richness, by = "RTENO")


# ============================================================
# 10. Final selection and sanity checks
# ============================================================

routes_final <- routes_final %>%
  select(
    RTENO,
    RouteName,
    StateNum,
    Route,
    CountryNum,
    species_richness,
    geometry
  )
head(routes_final)
# Quick checks
summary(routes_final$species_richness)
nrow(routes_final)

routes_lines <- read_sf("routes_geometry/vy474dv5024.shp")
head(routes_lines)
routes_lines$rteno
routeroutes_linesroutes_lines <- routes_lines %>%
  mutate(RTENO_clean = as.character(RTENO))
plot(routes_lines)

# Function to check if a column is unique
# Explore which columns are unique
sapply(routes_lines, function(x) length(unique(x)))
sapply(routes_final, function(x) length(unique(x)))


routes_final <- routes_final %>%
  left_join(
    routes_lines %>% select(RTENO_clean, geometry),
    by = "RTENO_clean",
    suffix = c("_pt", "")
  ) %>%
  st_drop_geometry() %>%
  st_as_sf(crs = 4326)


# Ensure both layers use the same CRS
routes_lines <- st_transform(routes_lines, 4326)
routes_final <- st_transform(routes_final, 4326)


# Spatial join: add point attributes to lines
routes_spatial_join <- st_join(
  routes_lines,
  routes_final,
  join = st_nearest_feature # finds the nearest point for each line
)
routes_spatial_join <- routes_spatial_join |> 
  select(RTENO, RouteName, CountryNum, StateNum, Route, species_richness, length, geometry)

library(ggplot2)

ggplot() +
  geom_sf(data = routes_spatial_join, color = "green", size = 1)
# ============================================================
# 11. Write output for Python
# ============================================================

st_write(
  routes_final,
  "bbs_routes_2018_route_level.gpkg",
  layer = "routes",
  delete_dsn = TRUE
)
