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


# 1. Preparing BBS data ---------------------------------------------------

# Download BBS data
# bbs <- grab_bbs_data()  Using bbsAssistent package
# load BBS dataset
# Merge BBS dataframes and clean the data

here::here(load("data/bbs_data/bbs_dataset.RData")) # Release 2019

# select methodologicaly suitable routes for chosen years
weather <- bbs$weather |>
  dplyr::filter(RunType == 1, Year %in% c(2017, 2018, 2019)) |>
  dplyr::select(RTENO, Year)

# extract routes geometries and convert them to a spatial object
routes_sf <- bbs$routes |>
  dplyr::select(RTENO, Latitude, Longitude) |>
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# select only routes with valid weather conditions
valid_routes_years <- weather |>
  dplyr::inner_join(routes_sf, by = "RTENO") |>
  sf::st_as_sf()

# filter only passerines
species_list <- bbs$species_list |>
  dplyr::filter(ORDER == "Passeriformes")|>
  dplyr::select(Scientific_Name, AOU)

# filter observations for the target years and add scientific names
observations <-bbs$observations |>
  dplyr::filter(Year %in% c(2017, 2018, 2019)) |>
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


# 3. Get the species ranges from ebirdst ----------------------------------

# An access key is required to download eBird Status and Trends data
# set_ebirdst_access_key("")
# Download species ranges
# for( i in sister_names){
#   ebirdst_download_status(species = i, path = ebirdst_data_dir(),
#                           download_abundance = FALSE,
#                           download_occurrence = FALSE,
#                           download_count = FALSE,
#                           download_ranges = TRUE,
#                           download_regional = FALSE,
#                           download_pis = FALSE,
#                           download_ppms = FALSE,
#                           download_all = FALSE,
#                           pattern = NULL,
#                           dry_run = FALSE,
#                           force = FALSE,
#                           show_progress = TRUE)}
# Make a list of ranges
ranges <- lapply(sister_names, function(species) {
  load_ranges(
    species = species,
    resolution = "27km",
    smoothed = TRUE,
    path = ebirdst_data_dir()
  )
})

# Save objects needed for downstream scripts
save(
  dat,
  sister_pairs,
  sister_names,
  ranges,
  file = here::here("data/imported_bbs_data.RData")
)
