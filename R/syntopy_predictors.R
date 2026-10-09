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

here::here(load("data/data_import.RData"))

# 1. Calculate degree of sympatry for each pair -------------------------------------

# Assign species scientific names as list element names
names(ranges) <- sister_names

# Verify that indexing by species name now works
ranges[["Passerina ciris"]]

# Check if all element names match their internal scientific_name
all(names(ranges) == sapply(ranges, function(x) x$scientific_name[1]))

# Source function
source(here::here("R/functions/calculate_sympatry.R"))

# Calculate sympatry for every pair in sister_pairs

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
  dplyr::inner_join(bbs_occurrences, by = c("sp1" = "Scientific_Name"),relationship = "many-to-many") |>
  dplyr::inner_join(bbs_occurrences, by = c("sp2" = "Scientific_Name", "RTENO" = "RTENO", "Year" = "Year")) |>
  dplyr::group_by(Year, sp1, sp2) |>
  dplyr::summarise(n_shared_routes = dplyr::n(), .groups = "drop")

# Summary count of sister pairs sharing at least one route per year
route_syntopy_summary <- route_syntopy_by_year |>
  dplyr::group_by(Year) |>
  dplyr::summarise(pairs_sharing_routes = dplyr::n_distinct(sp1, sp2))
# There are around 50 sister pairs sharing at least one route for each year
# We need to select all 59 sympatric pairs (despite their absence of shared routes),
# because we are interested in the degree of syntopy in the sympatric zone

# 3. Range symmetry and other range metrics -------------------------------

# Source function
source(here::here("R/functions/calculate_symmetry.R"))

# Calculate range symmetry

sympatric_pairs <- sympatric_pairs |>
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

sympatric_pairs <- sympatric_pairs |>
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
hist(sympatric_pairs$sympatry, breaks = 12, main="Degree of Sympatry", xlab="Range overlap (%)", col="darkseagreen3")
hist(sympatric_pairs$symmetry, breaks = 12, main="Degree of Range Symmetry", xlab="Range symmetry (%)", col="darkseagreen4")
hist(sympatric_pairs$age_myr, breaks = 12, main="Evolutionary age", xlab="Evolutionary age (myr)", col="steelblue")
hist(sympatric_pairs$centroid_dist_km, breaks = 12, main="Distance Between Species Range Centroids", xlab="Geographical Distance (Km)", col="brown")

# 4. Syntopy on route scale ----------------------------------------------



# Save objects needed for downstream scripts
# save(
#   dat,
#   sympatric_pairs,
#   sister_names,
#   ranges,
#   file = here::here("data/syntopy_pred_dat.RData")
# )
