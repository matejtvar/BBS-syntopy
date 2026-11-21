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
library(CooccurrenceAffinity)

# 1. Syntopy on route scale ----------------------------------------------

# Selecting year
year <- 2019

dat_year <- dat |>
  dplyr::filter(Year == year) |>
  sf::st_as_sf(crs = 4326)

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

route_data <- bind_cols(sympatric_pairs, route_syntopy)

# Filter out pairs with no routes sampled (N == 0) for modeling
route_data <- route_data |>
  filter(N > 0)
# Note: for 2019 and 2017 there is just 1 pair with N=0, while for 2018 there are 2 pairs

# Calculate empirical co-occurrence indicies - alpha_mle is not working!
# route_data_indicies <- route_data |>
#   dplyr::filter(m1 > 0 & m2 > 0) |>
#   dplyr::rowwise() |>
#   dplyr::mutate(
#     # affinity2by2 returns a list containing alpha_mle, jaccard, sorensen, simpson
#     res = list(CooccurrenceAffinity::affinity2by2(x = k, marg = c(m1, m2, N))),
#     empirical_log_OR = res$alpha_mle,
#     jaccard          = res$jaccard,
#     sorensen         = res$sorensen,
#     simpson          = res$simpson
#   ) |>
#   ungroup() |>
#   dplyr::select(-res)

# 2. Syntopy on stop scale ------------------------------------------------
syntopic_routes <- route_data |>
  dplyr::filter(k > 0) |> # only pairs co-occurring at some routes can be used
  dplyr::select(sp1, sp2)

# Source the function
source(here::here("R/functions/calculate_stop_syntopy.R"))

# Map across sympatric pairs
stop_syntopy <- map2_dfr(
  syntopic_routes$sp1,
  syntopic_routes$sp2,
  \(s1, s2) {
    calculate_stop_syntopy(
      sp1 = s1,
      sp2 = s2,
      bbs_dat = dat_year
    )
  }
)

# Calculate empirical co-occurrence indicies - alpha_mle is not working!
# stop_data_indicies <- stop_syntopy |>
#   dplyr::filter(m1 > 0 & m2 > 0) |>
#   dplyr::rowwise() |>
#   dplyr::mutate(
#     # affinity2by2 returns a list containing alpha_mle, jaccard, sorensen, simpson
#     res = list(CooccurrenceAffinity::affinity2by2(x = k, marg = c(m1, m2, N))),
#     empirical_log_OR = res$alpha_mle,
#     jaccard          = res$jaccard,
#     sorensen         = res$sorensen,
#     simpson          = res$simpson
#   ) |>
#   ungroup() |>
#   dplyr::select(-res)

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

# 4. Environmental predictors ---------------------------------------------
# list2env(c(data_import, syntopy_predictors), envir = environment())
ndvi_yr <- ndvi_sf |>
  dplyr::filter(Year == year)

stop_data_ndvi <- stop_data |>
  dplyr::left_join(ndvi_yr, by = c("RTENO"))

# 5. Trait and niche similarity predictors --------------------------------

# Calculate Gower distance
# gower_dist <- cluster::daisy(avoniche, metric = "gower")
# # Convert to a standard square matrix for easy lookup
# gower_matrix <- as.matrix(gower_dist)
# gawdis package?

list(
  route_data = route_data,
  stop_data = stop_syntopy,
  range_metrics = range_metrics
)
