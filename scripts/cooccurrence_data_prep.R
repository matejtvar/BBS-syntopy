
# 1. Sister pairs extraction from phylogeny -----------------------------------

# Libraries
require(diverge)
require(clootl)
require(tidyverse)

# Selecting taxonomy from eBird using clootl package
taxonomy <- taxonomyGet(taxonomy_year = 2023)
# Selecting passerine family
names_vec <- taxonomy |> 
  filter(ORDER1 == "Passeriformes") |> 
  select(SCI_NAME) |> 
  pull(SCI_NAME)

# ==== Names check
# Ensure names are clean (remove trailing/leading whitespace, ensure the taxonomy matching)
names_vec <- trimws(names_vec)
valid_sci_names <- taxonomy$SCI_NAME
final_names <- names_vec[names_vec %in% valid_sci_names]

# Check how many were matched
message(paste("Matched", length(final_names), "out of", length(names_vec), "species."))

# Extract the tree
tree <- extractTree(species = final_names, taxonomy_year = 2023)

# Extract sister species pairs using extract_sisters() function from diverge package
sister_pairs <- extract_sisters(tree)

# Inspect the results
head(sister_pairs)

# With the sister pairs of passerines extracted, we can look which of those are presented in eBird S&T and BBS data

# 2. Identification of available species data - ebirdst & BBS ----------------------

# ==== eBird S&T
library(ebirdst)

# Get the list of species that actually have eBird S&T models
available_ebird_models <- ebirdst_runs |> 
  select(species_code, scientific_name) |> 
  mutate(sci_name_underscore = gsub(" ", "_", scientific_name))

# Create a lookup vector for quick checking
ebird_species_list <- available_ebird_models$sci_name_underscore

sister_pairs_ebirdst <- sister_pairs |> 
  # Check if Species 1 is in eBird S&T
  filter(sp1 %in% ebird_species_list) |> 
  # Check if Species 2 is in eBird S&T
  filter(sp2 %in% ebird_species_list)


message("Final N of sister pairs with eBird models: ", nrow(sister_pairs_ebirdst))
# So we have about 262 sister pairs of passerines from ebirdst!
# Now we can look at the BBS data and the overlap with ebirdst

# ==== BBS
# Selecting species list from 2019 released dataset → downloaded with bbsAssistent() package
bbs_species_list <- bbs$species_list |> 
  select(Scientific_Name) |> 
  mutate(sci_name_underscore = gsub(" ", "_", Scientific_Name))

bbs_species_list_vec <- bbs_species_list$sci_name_underscore

sister_pairs_bbs <- sister_pairs_ebirdst |> 
  # Check if Species 1 is in BBS
  filter(sp1 %in% bbs_species_list_vec) |> 
  # Check if Species 2 is in BBS
  filter(sp2 %in% bbs_species_list_vec)

message("Final N of sister pairs: ", nrow(sister_pairs_bbs))
check <- unique(c(sister_pairs_bbs$sp1, sister_pairs_bbs$sp2))
glimpse(check) # each column has a unique set of species
sister_pairs_total <- sister_pairs_bbs

# The max count of sister pairs for BBS and ebird S&T data overlap is 54!
# Could there be a potential taxonomical mismatches and should I filter with rather with species codes?

# 3. Selecting ranges + range overlap calculation ----------------------------
library(sf)
library(terra)
library(ebirdst)
library(tidyverse)


# Get a unique list of all species involved in sister pairs
spp_to_download <- unique(c(sister_pairs_total$sp1, sister_pairs_total$sp2))

# Remove underscores to match eBird's scientific name format
spp_names_clean <- gsub("_", " ", spp_to_download)

# Create a lookup table using the internal ebirdst taxonomy
# This ensures we get the exact 6-letter code needed for load_ranges()
ebird_lookup <- ebirdst_runs |>  
  filter(scientific_name %in% spp_names_clean) %>% 
  select(species_code, scientific_name)

# Extract just the codes
codes_to_download <- ebird_lookup$species_code

# An access key is required to download eBird Status and Trends data
set_ebirdst_access_key("56fnsc0aep49")

# Download the data (ranges only to save space/time)
for (code in codes_to_download) {
  ebirdst_download_status(
    species = code,
    path = ebirdst_data_dir(),
    download_abundance = FALSE,
    download_occurrence = FALSE,
    download_count = FALSE,
    download_ranges = TRUE,
    download_regional = FALSE,
    download_pis = FALSE,
    download_ppms = FALSE,
    download_all = FALSE,
    pattern = NULL,
    dry_run = FALSE,
    force = FALSE,
    show_progress = TRUE
  )
}

# Load the breeding ranges (BBS data represents the breeding season)
# Using 'low resolution' (lr) is usually sufficient for range overlap
ranges <- lapply(codes_to_download, function(code) {
  load_ranges(species = code,
              resolution = c("9km", "27km"),
              smoothed = TRUE,
              path = ebirdst_data_dir()
  )
})



# Get the codes of species that were actually successfully loaded
available_codes <- names(ranges)

# Create a logical check for sister pairs
# We keep only pairs where BOTH sp1 and sp2 have a corresponding range in the list
sister_pairs_check <- sister_pairs_total |> 
  filter(sp1_code %in% available_codes & sp2_code %in% available_codes)

message(paste("Original pairs:", nrow(sister_pairs_total)))
message(paste("Pairs with available range data:", nrow(sister_pairs_check)))