
# 1. Sister pairs extraction from phylogeny -----------------------------------

# Libraries
require(diverge)
require(clootl)
require(tidyverse)
require(ape)

# Taxonomy from eBird
taxonomy <- taxonomyGet(taxonomy_year = 2023)
glimpse(taxonomy)

names_vec <- taxonomy |> 
  filter(ORDER1 == "Passeriformes") |> 
  select(SCI_NAME) |> 
  pull(SCI_NAME)

glimpse(names_vec)

# ==== Names check
# 1. Ensure your names are clean (remove trailing/leading whitespace)
names_vec <- trimws(names_vec)

# 2. Get the official list of scientific names from the 2021 taxonomy
# This ensures we only ask for species that clootl definitely recognizes
valid_sci_names <- taxonomy$SCI_NAME

# 3. Filter your vector (Remeš's defensive approach)
# This removes any species that might have changed names in the 2021 version
final_names <- names_vec[names_vec %in% valid_sci_names]

# Check how many were matched
message(paste("Matched", length(final_names), "out of", length(names_vec), "species."))

# 4. Extract the tree
# Since it's a vector of scientific names, clootl handles it automatically
tree <- extractTree(species = final_names, taxonomy_year = 2023)

# 5. Extract sister species pairs (as per Remeš MS)
# This is the 'syntopy' preparation step from the original script
sister_pairs <- extract_sisters(tree)
# Inspect the results
head(sister_pairs)

# With the sister pairs of passerines extracted, we can look which of those are presented in ebirdst data

# 2. Identification of available species data - ebirdst ----------------------

library(ebirdst)

# Get the list of species that actually have eBird S&T models
# These are the species you can actually run load_ranges() on.
available_ebird_models <- ebirdst_runs |> 
  select(species_code, scientific_name) |> 
  mutate(sci_name_underscore = gsub(" ", "_", scientific_name))

# Create a lookup vector for quick checking
ebird_species_list <- available_ebird_models$sci_name_underscore

sister_pairs_filtered <- sister_pairs %>%
  # 1. Check if Species 1 is in eBird S&T
  filter(sp1 %in% ebird_species_list) %>%
  # 2. Check if Species 2 is in eBird S&T
  filter(sp2 %in% ebird_species_list)


message("Final N of sister pairs with eBird models: ", nrow(sister_pairs_filtered))
# So we have about 262 sister pairs of passerines from ebirdst!
# Now we can look at the BBS data and the overlap with ebirdst


# 3. BBS & ebirdst sp overlap ---------------------------------------------
bbs_species_list <- bbs$species_list |> 
  select(Scientific_Name) |> 
  mutate(sci_name_underscore = gsub(" ", "_", Scientific_Name))

bbs_species_list_vec <- bbs_species_list$sci_name_underscore

final_species_count <- sister_pairs_filtered |> 
  # 1. Check if Species 1 is in BBS
  filter(sp1 %in% bbs_species_list_vec) |> 
  # 2. Check if Species 2 is in BBS
  filter(sp2 %in% bbs_species_list_vec)

message("Final N of sister pairs: ", nrow(final_species_count))
check <- unique(c(final_species_count$sp1, final_species_count$sp2))
glimpse(check)

# The max count of sister pairs for BBS and ebird S&T data overlap is 54!
# But I have to check the potential taxonomical mismatches and filter with codes instead!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 4. Selecting ranges + range overlap calculation ----------------------------
library(sf)
library(terra)
library(ebirdst)
library(tidyverse)


# Get a unique list of all species involved in sister pairs
spp_to_download <- unique(c(sister_pairs$sp1, sister_pairs$sp2))

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

ebirdst_data_dir()

# Load the breeding ranges (BBS data represents the breeding season)
# Using 'low resolution' (lr) is usually sufficient for range overlap
warbler_ranges <- lapply(codes_to_download, function(code) {
  load_ranges(species = code,
              resolution = c("9km", "27km"),
              smoothed = TRUE,
              path = ebirdst_data_dir()
  )
})

# Get the codes of species that were actually successfully loaded
available_codes <- names(warbler_ranges)

# Create a logical check for your sister pairs
# We keep only pairs where BOTH sp1 and sp2 have a corresponding range in the list
sister_pairs_ready <- sister_pairs |> 
  filter(sp1_code %in% available_codes & sp2_code %in% available_codes)

message(paste("Original pairs:", nrow(sister_pairs_final)))
message(paste("Pairs with available range data:", nrow(sister_pairs_ready)))


glimpse(ebird_lookup)

# 1. Clean names for joining (convert underscores to spaces to match scientific_name)
sister_pairs_mapped <- sister_pairs %>%
  mutate(sp1_clean = gsub("_", " ", sp1),
         sp2_clean = gsub("_", " ", sp2))

# 2. Join the lookup table to get the codes for sp1 and sp2
sister_pairs_final <- sister_pairs_mapped %>%
  left_join(ebird_lookup, by = c("sp1_clean" = "scientific_name")) %>%
  rename(sp1_code = species_code) %>%
  left_join(ebird_lookup, by = c("sp2_clean" = "scientific_name")) %>%
  rename(sp2_code = species_code)

# 3. Filter: Only keep rows where both species have valid codes
sister_pairs_ready <- sister_pairs_final %>%
  filter(!is.na(sp1_code) & !is.na(sp2_code))

message("Pairs retained after mapping: ", nrow(sister_pairs_ready))

# Subset to breeding season
# Subset using the column name directly
sp_range <- subset(warbler_ranges[["ovenbi"]])

# Plot the range
plot(sp_range, col = "darkgreen", main = "Ovenbird Breeding Range")

warbler_ranges

library(leaflet)
library(sf)

# Convert your breeding range to sf
# Replace 'ovenbi_range' with your subsetted SpatVector
altyel1_sf <- st_as_sf(subset(warbler_ranges[[1]], warbler_ranges[[1]]$season == "breeding"))

# Convert both species in a pair
sp1_sf <- st_as_sf(subset(warbler_ranges[[code1]], season == "breeding"))
sp2_sf <- st_as_sf(subset(warbler_ranges[[code2]], season == "breeding"))

# Create interactive map
leaflet() %>%
  addTiles() %>% # Adds base map (OpenStreetMap)
  addPolygons(data = sp1_sf, color = "red", weight = 1, fillOpacity = 0.5, group = "Species 1") %>%
  addPolygons(data = sp2_sf, color = "blue", weight = 1, fillOpacity = 0.5, group = "Species 2") %>%
  addLayersControl(overlayGroups = c("Species 1", "Species 2"), options = layersControlOptions(collapsed = FALSE))

# Function to generate map for any pair
plot_pair_map <- function(code1, code2) {
  s1 <- st_as_sf(subset(warbler_ranges[[code1]], season == "breeding"))
  s2 <- st_as_sf(subset(warbler_ranges[[code2]], season == "breeding"))
  
  leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(data = s1, color = "red", fillOpacity = 0.4) %>%
    addPolygons(data = s2, color = "blue", fillOpacity = 0.4)
}


"C:\Users\Matěj Tvarůžka\AppData\Roaming\R\data\R\ebirdst\2023\elwwar1"