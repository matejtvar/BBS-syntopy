
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



# Step 1: Downloading & tiding data --------------------------------------------------------

library(tidyverse) # for data transformations
library(sf) # for spatial data manipulations
library(here) # for selecting the file path within the project

# Downloading the data in the form of list
bbs <- grab_bbs_data()
names(bbs)
head(bbs$species_list,10)
head(bbs$observations)
head(bbs$routes)
str(bbs$observations)

# Convert route coordinates to an sf object
routes_sf <- bbs$routes |> 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) # WGS84

bbs_data <- bbs$observations |> 
  filter(Year %in% c(2017, 2018, 2019)) |> 
  pivot_longer(
    cols = starts_with("Stop"),
    names_to = "Stop_Number",
    values_to = "Count"
  ) |> 
  mutate(
    Stop_PA = if_else(!is.na(Count) & Count > 0, 1L, 0L)
  )

# Join routes with observations
bbs_data <- routes_sf |> 
  left_join(bbs_data, by = c("RTENO", "CountryNum", "StateNum", "Route"))

bbs_data <- bbs_data |> 
  filter(Active == 1) |> 
  select(RTENO, Year, AOU, geometry, Stop_PA, everything(), 
         -RouteName, -(RouteDataID:RPID), -(Active:RouteTypeDetailID), -Count, -(CountryNum:Route))
bbs_data
bbs_data_wider <- bbs_data |> 
  pivot_wider(
    names_from  = Stop_Number,
    values_from = Stop_PA
  )
bbs_data_wider
