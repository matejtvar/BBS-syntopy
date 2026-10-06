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

# 1. Calculate degree of sympatry for each pair -------------------------------------

# create empty vector
sister_pairs$sympatry <- vector(mode = "numeric", length = nrow(sister_pairs))

# Assign species scientific names as list element names
names(ranges) <- sister_names

# Verify that indexing by species name now works
ranges[["Passerina ciris"]]

# Check if all element names match their internal scientific_name
all(names(ranges) == sapply(ranges, function(x) x$scientific_name[1]))

# Calculate sympatry for every pair in sister_pairs
# using breeding and resident ranges

sister_pairs <- sister_pairs |>
  dplyr::mutate(
    sympatry = map2_dbl(sp1, sp2, \(sp1, sp2) {
      calculate_sympatry(
        sp1 = sp1, sp2 = sp2,
        range_list = ranges, season_filter = c("breeding", "resident"))
    })
  )

# Check with map
ranges_sf <- dplyr::bind_rows(ranges)

ranges_sf |>
  dplyr::filter(scientific_name %in% c("Piranga flava", "Piranga rubra")) |>
  dplyr::filter(season == "breeding") |>
  ggplot() +
  geom_sf(aes(fill = scientific_name), alpha = 0.5, color = "black") +
  theme_minimal() +
  scale_fill_manual(
    values = c("Piranga flava" = "firebrick", "Piranga rubra" = "steelblue")
  ) +
  labs(
    title = "Test Dataset Range Overlap",
    fill = "Species"
  )

# 2. Check the number of valid pairs for analysis -------------------------

# Filter sister pairs for sympatry threshold 5 %
sympatric_pairs <- sister_pairs |>
  dplyr::filter(sympatry > 5)
nrow(sympatric_pairs)
# There are 59 sympatric sister pairs in the data over the 3 selected years

# Extract unique species-route-year occurrences from BBS data
bbs_occurrences <- dat |>
  sf::st_drop_geometry() |>
  dplyr::select(RTENO, Year, Scientific_Name) |>
  dplyr::distinct()

# Check for presence of each species in each year
# (Species present in a given year if observed on >= 1 route)
species_by_year <- bbs_occurrences |>
  dplyr::group_by(Year, Scientific_Name) |>
  dplyr::summarise(detected = TRUE, .groups = "drop")

# Determine co-presence for each sympatric pair by year
pairs_by_year <- sympatric_pairs |>
  dplyr::cross_join(tibble::tibble(Year = c(2017, 2018, 2019))) |>
  dplyr::left_join(
    species_by_year,
    by = c("sp1" = "Scientific_Name", "Year" = "Year")
  ) |>
  dplyr::rename(sp1_present = detected) |>
  dplyr::left_join(
    species_by_year,
    by = c("sp2" = "Scientific_Name", "Year" = "Year")
  ) |>
  dplyr::rename(sp2_present = detected) |>
  dplyr::mutate(
    sp1_present = tidyr::replace_na(sp1_present, FALSE),
    sp2_present = tidyr::replace_na(sp2_present, FALSE),
    both_present = sp1_present & sp2_present
  )

# Count of sympatric pairs detected per year
annual_summary <- pairs_by_year |>
  dplyr::group_by(Year) |>
  dplyr::summarise(
    total_sympatric_pairs = dplyr::n(),
    sp1_detected = sum(sp1_present),
    sp2_detected = sum(sp2_present),
    both_detected_in_bbs = sum(both_present)
  )

# Find routes where both sp1 and sp2 were observed in the same year
route_syntopy_by_year <- sympatric_pairs |>
  dplyr::inner_join(bbs_occurrences, by = c("sp1" = "Scientific_Name"), relationship = "many-to-many") |>
  dplyr::inner_join(bbs_occurrences, by = c("sp2" = "Scientific_Name", "RTENO" = "RTENO", "Year" = "Year")) |>
  dplyr::group_by(Year, sp1, sp2) |>
  dplyr::summarise(n_shared_routes = dplyr::n(), .groups = "drop")

# Summary count of sister pairs sharing at least one route per year
route_syntopy_summary <- route_syntopy_by_year |>
  dplyr::group_by(Year) |>
  dplyr::summarise(pairs_sharing_routes = dplyr::n_distinct(sp1, sp2))

