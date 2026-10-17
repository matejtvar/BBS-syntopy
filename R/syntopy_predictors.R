# ---- BBS Syntopy Analysis ----
# Author: Matěj Tvarůžka
# Description: Data analysis of syntopy predictors with BBS dataset
# GitHub repository: https://github.com/matejtvar/BBS-syntopy

here::i_am("R/syntopy_predictors.R")
library(here)
here() # project path

# Libraries
library(dplyr)
library(sf)
library(purrr)
library(assertthat)
library(testthat)
library(usethis)
library(tidyr)
library(ggplot2)
library(targets)
library(cluster)
library(tibble)
library(readxl)

# 1. Syntopy on route scale ----------------------------------------------

# Selecting year
year <- 2019

dat_year <- dat |>
  dplyr::filter(Year == year) |>
  sf::st_as_sf(crs = 4326)

# Sour Fisher's Hypergeometric route syntopy function
source(here::here("R/functions/calculate_route_syntopy.R"))

# Map calculate_route_syntopy across all sympatric sister pairs
route_syntopy <- map2_dfr(
  sympatric_pairs$sp1,
  sympatric_pairs$sp2,
  \(s1, s2) {
    calculate_route_syntopy(
      sp1 = s1,
      sp2 = s2,
      range_list = ranges,
      bbs_dat = dat_year,
      season_filter = c("breeding", "resident")
    )
  }
)

# Bind the resulting hypergeometric counts (N, m1, m2, k, empirical_log_OR) directly
route_data <- bind_cols(sympatric_pairs, route_syntopy)

# Filter out pairs with no routes sampled (N == 0) for modeling
route_data <- route_data |>
  filter(N > 0, !is.na(empirical_log_OR))
# Note: for 2019 and 2017 there is just 1 pair with N=0, while for 2018 there are 2 pairs

# 2. Syntopy on stop scale ------------------------------------------------
stop_data <- route_data |>
  dplyr::filter(k > 0) |>
  dplyr::select(sp1, sp2)

# Source the function
source(here::here("R/functions/calculate_stop_syntopy.R"))

# Map across sympatric pairs (using your filtered route_data or sympatric_pairs)
stop_syntopy_results <- map2_dfr(
  stop_data$sp1,
  stop_data$sp2,
  \(s1, s2) {
    calculate_stop_syntopy(
      sp1 = s1,
      sp2 = s2,
      bbs_dat = dat_year
    )
  }
)

# Bind the stop metrics back into your dataset
stop_data <- bind_cols(stop_data, stop_syntopy_results)

# 3. Range symmetry and other range metrics -------------------------------

# Source function
source(here::here("R/functions/calculate_symmetry.R"))

# Calculate range symmetry

range_metrics <- sympatric_pairs |>
  dplyr::mutate(
    symmetry = map2_dbl(sp1, sp2, \(sp1, sp2) {
      calculate_symmetry(
        sp1 = sp1, sp2 = sp2,
        range_list = ranges, season_filter = c("breeding", "resident"))
    })
  )

# Source function
source(here::here("R/functions/calculate_cen_dist.R"))

# Calculate centroid distance of the ranges for every pair in sister_pairs

range_metrics <- range_metrics |>
  mutate(
    centroid_dist_km = map2_dbl(sp1, sp2, \(s1, s2) {
      calculate_centroid_distance(
        sp1 = s1,
        sp2 = s2,
        range_list = ranges,
        season_filter = c("breeding", "resident"),
        unit = "km"
      )
    })
  )

# Visual exploration of predictors
par(mfrow = c(2,2))
hist(range_metrics$sympatry, breaks = 12, main="Degree of Sympatry", xlab="Range overlap (%)", col="darkseagreen3")
hist(range_metrics$symmetry, breaks = 12, main="Degree of Range Symmetry", xlab="Range symmetry (from 0 to 0.5)", col="darkseagreen4")
hist(range_metrics$age_myr, breaks = 12, main="Evolutionary age", xlab="Evolutionary age (myr)", col="steelblue")
hist(range_metrics$centroid_dist_km, breaks = 12, main="Distance Between Species Range Centroids", xlab="Geographical Distance (Km)", col="steelblue2")
par(mfrow = c(1,1))

# We have to standardise predictors prior to modelling
range_metrics <- range_metrics |>
  dplyr::mutate(sympatry = (sympatry - mean(sympatry)) / sd(sympatry)) |>
  dplyr::mutate(symmetry = (symmetry - mean(symmetry)) / sd(symmetry)) |>
  dplyr::mutate(age = (age_myr - mean(age_myr)) / sd(age_myr)) |>
  dplyr::mutate(cen_dist = (centroid_dist_km - mean(centroid_dist_km)) / sd(centroid_dist_km)) |>
  dplyr::select(-age_myr, -centroid_dist_km)

# 4. Trait and niche similarity predictors --------------------------------

# Calculate Gower distance
# gower_dist <- cluster::daisy(avoniche, metric = "gower")
# # Convert to a standard square matrix for easy lookup
# gower_matrix <- as.matrix(gower_dist)
# gawdis package?

# 5. Environmental predictors ---------------------------------------------
# list2env(c(data_import, syntopy_predictors), envir = environment())


list(
  route_data = route_data,
  stop_data = stop_data,
  range_metrics = range_metrics
)
