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
library(readxl)

# 1. Preparing BBS data ---------------------------------------------------

# Download BBS data
# bbs <- grab_bbs_data()  Using bbsAssistent package

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

# plot(tree, type = "fan", cex = 0.3, tip.color = "darkblue")

# extract sister species pairs
sister_pairs <- diverge::extract_sisters(tree)
dim(sister_pairs) # 96 species pairs

# extract sister species age
source(here::here("R/functions/extract_age.R"))
sister_pairs <- extract_sister_ages(sister_pairs = sister_pairs, phylo_tree = tree)

# 3. Get the species ranges from ebirdst ----------------------------------

# extract vector of sister species names
sp1 <- sister_pairs$sp1
sp2 <- sister_pairs$sp2
sister_names <- c(sp1,sp2)
length(sister_names)
# We need to get range data for 192 species

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

# Assign species scientific names as list element names
names(ranges) <- sister_names

# Verify that indexing by species name now works
ranges[["Passerina ciris"]]

# Check if all element names match their internal scientific_name
all(names(ranges) == sapply(ranges, function(x) x$scientific_name[1]))

# 4. Calculate degree of sympatry -------------------------------------

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

species_1 <- "Piranga flava"
species_2 <- "Piranga rubra"

ranges_sf |>
  dplyr::filter(scientific_name %in% c(species_1, species_2)) |>
  dplyr::filter(season == "breeding") |>
  ggplot() +
  geom_sf(aes(fill = scientific_name), alpha = 0.5, color = "black") +
  theme_minimal() +
  scale_fill_manual(
    values = stats::setNames(c("firebrick", "steelblue"), c(species_1, species_2))
  ) +
  labs(
    title = "Test Dataset Range Overlap",
    fill = "Species"
  )

# 5. Check the number of valid pairs for analysis -------------------------

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


# 6. Import Avoniche and Avonet -------------------------------------------

# Avoniche data
# Available here: https://datadryad.org/dataset/doi:10.5061/dryad.1zcrjdg56

avoniche <- read.csv(here::here("data/Avoniche2_eBird2021.csv"))

avoniche <- avoniche |>
  dplyr::filter(Species2 %in% sister_names) |>
  dplyr::select(-Order2, -Avibase_id)

# Avonet data
# Available here: https://figshare.com/s/b990722d72a26b5bfead

avonet <- readxl::read_excel("data/AVONET Supplementary dataset 1.xlsx",
                             sheet = "AVONET2_eBird")
avonet <- avonet |>
  dplyr::filter(Species2 %in% sister_names) |>
  dplyr::select(-Order2, -Avibase.ID2)


# End list

list(
  dat                   = dat,
  sister_pairs          = sister_pairs,
  sympatric_pairs       = sympatric_pairs,
  ranges                = ranges,
  sister_names          = sister_names,
  avoniche              = avoniche,
  avonet                = avonet
)
