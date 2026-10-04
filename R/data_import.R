# ---- BBS Syntopy Analysis ----
# Author: Matěj Tvarůžka
# Description: Data analysis of syntopy predictors with BBS dataset
# GitHub repository: https://github.com/matejtvar/BBS-syntopy

here::i_am("R/data_import.R")
library(here)
here() # project path

# To install from GitHub:
#remotes::install_github("ebird/ebirdst")
#remotes::install_github("trashbirdecology/bbsAssistant")

# Libraries
library(clootl)
library(dplyr)
library(bbsAssistant)
library(ape)
library(ebirdst)
library(diverge)
library(sf)
library(purrr)

testing
# 1. Preparing BBS data ---------------------------------------------------

# Download BBS data
# bbs <- grab_bbs_data()  Using bbsAssistent package
# load BBS dataset
# Merge BBS dataframes and clean the data

here::here(load("data/bbs_data/bbs_dataset.RData")) # Release 2019

# Define target years once
target_years <- c(2017, 2018, 2019)

# select methodologicaly suitable routes for chosen years
weather <- bbs$weather |>
  dplyr::filter(RunType == 1, Year %in% target_years) |>
  dplyr::select(RTENO, Year)
message(sprintf("Filtered weather: %d -> %d rows", nrow(bbs$weather), nrow(weather)))

# extract routes geometries and convert them to a spatial object
routes_sf <- bbs$routes |>
  dplyr::select(RTENO, Latitude, Longitude) |>
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

message(sprintf("Routes loaded: %d spatial features", nrow(routes_sf)))

# select only routes with valid weather conditions
valid_routes_years <- weather |>
  dplyr::inner_join(routes_sf, by = "RTENO") |>
  sf::st_as_sf()

# Validate
if (nrow(valid_routes_years) == 0) stop("No valid routes found after filtering!")
if (anyNA(sf::st_geometry(valid_routes_years))) warning("Some geometries are NA")

# filter only passerines
species_list <- bbs$species_list |>
  dplyr::filter(ORDER == "Passeriformes")|>
  dplyr::select(Scientific_Name, AOU)

# filter observations for the target years and add scientific names
observations <-bbs$observations |>
  dplyr::filter(Year %in% target_years) |>
  dplyr::inner_join(species_list, by = "AOU")  |>
  dplyr::select(-RouteDataID, -CountryNum, -StateNum, -Route, -RPID, -AOU)

# final join: Combine observations with valid spatial routes using RTENO and Year
dat <- observations |>
  dplyr::inner_join(valid_routes_years, by = c("RTENO", "Year")) |>
  dplyr::relocate(RTENO, Year, Scientific_Name, geometry)

# 2. Select sister taxa from the phylogeny -----------------------------------

# extract vector of species names
species_names <- ebirdst::ebirdst_runs |>
  dplyr::filter(scientific_name %in% dat$Scientific_Name) |>
  dplyr::pull(scientific_name)
length(species_names) # 284 species are present in both BBS and ebirdst

# check the phylogenetic tree
tree <- clootl::extractTree(species = species_names,
                         taxonomy_year = 2023,
                         version = "1.5",
                         data_path= here::here("data/AvesDataLite-main"))

plot(tree, type = "fan", cex = 0.3, tip.color = "darkblue")

# extract sister species pairs
sister_pairs <- diverge::extract_sisters(tree)
dim(sister_pairs) # 96 species pairs

# extract vector of sister species names
sp1 <- sister_pairs$sp1
sp2 <- sister_pairs$sp2
sister_names <- c(sp1,sp2)
length(sister_names)
# We need to get range data for 192 species

# 3. Get the species ranges from ebirdst ----------------------------------

# An access key is required to download eBird Status and Trends data
# set_ebirdst_access_key("") at the website https://ebird.org/st/request
# If you want to download range data, use the function below
# by selecting the code and pressing ctrl + shift + c to uncomment

# purrr::walk(sister_names, function(species) {
#   tryCatch(
#     {
#       ebirdst::ebirdst_download_status(
#         species = species,
#         path = ebirdst::ebirdst_data_dir(),
#         download_abundance = FALSE,
#         download_occurrence = FALSE,
#         download_count = FALSE,
#         download_ranges = TRUE,
#         download_regional = FALSE,
#         download_pis = FALSE,
#         download_ppms = FALSE,
#         download_all = FALSE,
#         pattern = NULL,
#         dry_run = FALSE,
#         force = FALSE,
#         show_progress = TRUE
#       )
#       message(sprintf("✓ Downloaded ranges for %s", species))
#     },
#     error = function(e) {
#       warning(sprintf("✗ Failed to download ranges for %s: %s", species, e$message))
#     }
#   )
# }, .progress = TRUE)

# Make a list of ranges
ranges <- purrr::map(sister_names, function(species) {
  tryCatch(
    load_ranges(species = species, resolution = "27km", smoothed = TRUE, path = ebirdst_data_dir()),
    error = function(e) {
      warning(sprintf("Failed to load ranges for %s: %s", species, e$message))
      NULL
    }
  )
}, .progress = TRUE)

# Check for failures
failed <- which(sapply(ranges, is.null))
if (length(failed) > 0) {
  warning(sprintf("%d species failed to load ranges", length(failed)))
}

# Save objects needed for downstream scripts
save(
  dat,
  sister_pairs,
  sister_names,
  ranges,
  file = here::here("data/imported_bbs_data.RData")
)
