# 1. Sister pairs extraction from phylogeny -----------------------------------

# Libraries
require(diverge)
require(clootl)
require(tidyverse)
require(terra)
require(sf)
require(ebirdst)

# Get the taxonomy from clootl
taxonomy <- taxonomyGet(taxonomy_year = tax_year)
# Selecting the vector of species codes
codes_vec <- taxonomy |> 
  filter(ORDER1 == "Passeriformes") |> 
  select(SPECIES_CODE) |> 
  pull(SPECIES_CODE)
clean_codes <- gsub("[0-9]", "", codes_vec)

# Extract the tree
tree_sci_names <- extractTree(species = sci_names_vec, taxonomy_year = 2023)
is.ultrametric(tree_sci_names) # False
# Extract sister species pairs using extract_sisters() function from diverge package
sci_names_sister <- extract_sisters(tree_sci_names)
# tree_codes <- extractTree(species = clean_codes, taxonomy_year = 2023)
# codes_sister <- extract_sisters(tree_codes)
# Selecting with species codes doesn't work

tree_extraction <- function(target_order, tax_year) {
  # Extracting phylogenetic tree for given order and taxonomy year
  message("Fetching taxonomy and building tree for: ", target_order)
  library(clootl)
  library(ape)
  tax_full <- taxonomyGet(taxonomy_year = tax_year)
  sci_names <- tax_full |> 
    filter(ORDER1 == target_order) |> 
    select(SCI_NAME) |> 
    pull(SCI_NAME)
  
  tree <- extractTree(species = sci_names, taxonomy_year = tax_year)
  tree_info <- list(
    phylo_tree = tree,
    ultrametric = is.ultrametric(tree),
    metadata = paste("Order:", target_order, "Year:", tax_year)
  )
  
  return(tree_info)
}
tree <- tree_extraction("Passeriformes", 2023)
sister_pairs <- extract_sisters(tree[[1]])
# 2. Identification of available species data - ebirdst & BBS ----------------------
# With the sister pairs of passerines extracted, we can look which of those are presented in eBird S&T and BBS data

extract_underscored_names <- function(df, col_name) {
  #Extract underscored species names from a data frame
  raw_names <- df |>  pull({{col_name}})
  
  # Step 2: Clean the names
  # GAP 2: Use gsub to replace " " with "_"
  clean_names <- gsub(" ", "_", raw_names)
  
  # Step 3: Return unique values only
  return(unique(clean_names))
}

filter_available_pairs <- function(pairs_df, lookup_vec) {
  # Filter pairs based on a lookup vector
  filtered_df <- pairs_df |> 
    filter(sp1 %in% lookup_vec) |> 
    filter(sp2 %in% lookup_vec)
  
  # Reporting is part of good programming!
  message("Input pairs: ", nrow(pairs_df), " -> Output pairs: ", nrow(filtered_df))
  
  return(filtered_df)
}

# Create Lookup Vectors
# ==== eBird S&T
ebird_available <- extract_underscored_names(ebirdst_runs, scientific_name)
# ==== BBS
# Selecting species list from 2019 released dataset → downloaded with bbsAssistent() package
bbs_available   <- extract_underscored_names(bbs$species_list, Scientific_Name)
# Filter the shared pairs between the datasets
ebirdst_pairs <- filter_available_pairs(sister_pairs, ebird_available)
bbs_ebirds_pairs <- filter_available_pairs(ebirdst_pairs, bbs_available)

# 3. Selecting ranges using ebirdst ----------------------------

# Create a lookup table using the internal ebirdst taxonomy
all_sisters_vec <- c(bbs_ebirds_pairs$sp1, bbs_ebirds_pairs$sp2)
spp_names_clean <- gsub("_", " ", unique(all_sisters_vec))
ebird_lookup <- ebirdst_runs |>  
  filter(scientific_name %in% spp_names_clean) |> 
  select(species_code, scientific_name)
# Extract codes vector
codes_to_download <- ebird_lookup$species_code

# An access key is required to download eBird Status and Trends data
# set_ebirdst_access_key("56fnsc0aep49")

#Download the data (ranges only to save space/time)
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
ranges <- lapply(codes_to_download, function(code) {
  load_ranges(species = code,
              resolution = c("9km", "27km"),
              smoothed = TRUE,
              path = ebirdst_data_dir()
  )
})

class(ranges)
ranges[[1]]

library(sf)
library(dplyr)

# 1. Corrected function to handle an individual sf object from the list
get_breeding_poly <- function(sf_item) {
  # Ensure we are working with the sf object inside the list element
  sf_item |> 
    filter(season %in% c("breeding", "resident")) |> 
    st_union() 
}

poly1 <- get_breeding_poly(ranges[[39]]) # Setophaga pensylvanica
poly2 <- get_breeding_poly(ranges[[107]]) # Setophaga petechia
# Find where they overlap
overlap_poly <- st_intersection(poly1, poly2)

library(ggspatial)

ggplot() +
  # 1. Add the basemap layer (zoomed to your data)
  # 'cartolight' or 'osm' are good neutral choices
  annotation_map_tile(type = "cartolight", zoomin = 0) + 
  
  # 2. Add your species layers
  geom_sf(data = poly1, fill = "#377eb8", alpha = 0.5, color = "black", size = 0.1) +
  geom_sf(data = poly2, fill = "#e41a1c", alpha = 0.5, color = "black", size = 0.1) +
  geom_sf(data = overlap_poly, fill = "purple", alpha = 0.7) +
  
  # 3. Ensure the coordinate system matches the tiles (usually Web Mercator)
  annotation_north_arrow(location = "bl", which_north = "true") +
  theme_minimal() +
  coord_sf()

# 4. Range overlap calculation -----------------------------------------------

names(ranges) <- codes_to_download
calculate_overlap <- function(code1, code2, range_list) {
  # 1. Extract polygons
  s1 <- range_list[[code1]]
  s2 <- range_list[[code2]]
  
  # 2. Filter for Breeding OR Resident
  s1_b <- s1[s1$season %in% c("breeding", "resident"), ]
  s2_b <- s2[s2$season %in% c("breeding", "resident"), ]
  
  # 3. If either is missing a range, return 0
  if (nrow(s1_b) == 0 | nrow(s2_b) == 0) return(0)
  
  # 4. Calculate areas (st_area returns units, so we convert to numeric)
  area1 <- as.numeric(st_area(s1_b))
  area2 <- as.numeric(st_area(s2_b))
  
  # 5. Find intersection
  # st_make_valid handles geometry errors often found in complex range maps
  inter <- st_intersection(st_make_valid(s1_b), st_make_valid(s2_b))
  
  if (nrow(inter) == 0) {
    return(0)
  } else {
    area_inter <- sum(as.numeric(st_area(inter)))
    # Index: Intersection / smaller of the two ranges
    return(area_inter / min(area1, area2))
  }
}

# Convert names to match scientific_name in ebird_lookup
sister_pairs <- sister_pairs_bbs |>
  mutate(sp1_clean = gsub("_", " ", sp1),
         sp2_clean = gsub("_", " ", sp2)) |>
  left_join(ebird_lookup, by = c("sp1_clean" = "scientific_name")) |>
  rename(code1 = species_code) |>
  left_join(ebird_lookup, by = c("sp2_clean" = "scientific_name")) |>
  rename(code2 = species_code)

# Run the calculation across all 54 pairs
sister_pairs$range_overlap <- sapply(1:nrow(sister_pairs), function(i) {
  calculate_overlap(
    sister_pairs$code1[i], 
    sister_pairs$code2[i], 
    ranges
  )
})
### It gives this warning:  attribute variables are assumed to be spatially constant throughout all geometries

ggplot(sister_pairs, aes(x = range_overlap)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  theme_minimal() +
  labs(
    title = "Distribution of Range Overlap in 54 Passerine Sister Pairs",
    x = "Symmetric Overlap Index (0 = Allopatric, 1 = Completely Overlapping)",
    y = "Number of Pairs"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1))

sum(sister_pairs$range_overlap > 0, na.rm = TRUE)
# 45 pairs have some degree of range overlap
sister_pairs$sympatry <- sister_pairs$range_overlap >= 0.05
sum(sister_pairs$sympatry == T)
# 34 pairs are sympatric for 5 % treshold

# 5. Inspecting N of BBS routes in sympatric zone -------------------------

# Convert BBS routes to a spatial object (sf)
bbs_routes_sf <- st_as_sf(bbs$routes, 
                          coords = c("Longitude", "Latitude"), 
                          crs = 4326)
# Filter data for 2018
observations <- bbs$observations |> 
  filter(Year == 2018)

count_routes_sympatry <- function(code1, code2, range_list, routes_sf) {
  # Function to count number of routes in sympatric range
  
  # Get breeding ranges
  r1 <- range_list[[code1]][range_list[[code1]]$season %in% c("breeding", "resident"), ]
  r2 <- range_list[[code2]][range_list[[code2]]$season %in% c("breeding", "resident"), ]
  
  if(nrow(r1) == 0 | nrow(r2) == 0) return(0)
  
  # Find the intersection polygon
  overlap_poly <- st_intersection(st_make_valid(r1), st_make_valid(r2))
  
  if(nrow(overlap_poly) == 0) return(0)
  
  # Count routes that fall inside this polygon
  routes_inside <- st_intersects(routes_sf, overlap_poly, sparse = FALSE)
  return(sum(routes_inside))
}

# lookup for AOU codes using sister_pairs names
aou_lookup <- bbs$species_list |> 
  mutate(sci_name_underscore = gsub(" ", "_", Scientific_Name)) |> 
  select(AOU, sci_name_underscore)

# Add AOU codes to sister_pairs
sister_pairs <- sister_pairs |> 
  left_join(aou_lookup, by = c("sp1" = "sci_name_underscore")) |> 
  rename(aou1 = AOU) |> 
  left_join(aou_lookup, by = c("sp2" = "sci_name_underscore")) |> 
  rename(aou2 = AOU) |> 
  select(-sp1_clean, sp2_clean)
names(sister_pairs)

count_shared_routes <- function(a1, a2, obs_df) {
  # Function to count shared BBS routes
  
  # Find routes where Sp1 was seen
  r_sp1 <- obs_df$RouteDataID[obs_df$AOU == a1]
  # Find routes where Sp2 was seen
  r_sp2 <- obs_df$RouteDataID[obs_df$AOU == a2]
  
  # Intersection of Route IDs
  shared_routes <- intersect(r_sp1, r_sp2)
  return(length(unique(shared_routes)))
}

data <- sister_pairs |>
  filter(sympatry == T) |> 
  mutate(
    sympatry_routes = map2_dbl(code1, code2, ~count_routes_sympatry(.x, .y, ranges, bbs_routes_sf)),
    bbs_shared_routes = map2_dbl(aou1, aou2, ~count_shared_routes(.x, .y, observations))
  )
summary(data)
# There is 34 sympatric sister pairs with min. 3 BBS routes
ggplot(data, aes(x = range_overlap)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  theme_minimal() +
  labs(
    title = "Distribution of Range Overlap in 34 Passerine Sister Pairs",
    x = "Degree of sympatry",
    y = "Number of Pairs"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1))

# Save the object
saveRDS(data, "sister_pairs_data.rds")
library(readr)
write_csv(data, "sister_pairs_results.csv")
# To load it back in:
# data <- readRDS("sister_pairs_data.rds")

# 6. Pres/Abs matrix ---------------------------------------------

sister_aou <- unique(c(data$aou1, data$aou2))

# 2. Create the Route-Level Matrix
route_matrix <- observations |> 
  filter(AOU %in% sister_aou) |> 
  # Group by Route and Species to see if they were EVER there
  group_by(RouteDataID, AOU) |> 
  summarise(Present = 1, .groups = "drop") |> 
  pivot_wider(names_from = AOU, values_from = Present, values_fill = 0)
head(route_matrix)

sympatric_routes_id <- function(code1, code2, range_list, routes_sf) {
  # Function to get the IDs of routes within the overlap zone
  r1 <- range_list[[code1]][range_list[[code1]]$season %in% c("breeding", "resident"), ]
  r2 <- range_list[[code2]][range_list[[code2]]$season %in% c("breeding", "resident"), ]
  
  if(nrow(r1) == 0 | nrow(r2) == 0) return(NULL)
  
  overlap_poly <- st_intersection(st_make_valid(r1), st_make_valid(r2))
  if(nrow(overlap_poly) == 0) return(NULL)
  
  # Logical vector of routes inside the polygon
  is_inside <- st_intersects(routes_sf, overlap_poly, sparse = FALSE)
  
  # Return the Route identifiers (adjust column name based on your bbs$routes structure)
  return(routes_sf$RouteDataID[is_inside])
}

