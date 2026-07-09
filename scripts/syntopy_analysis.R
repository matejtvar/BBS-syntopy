# ---- BBS Syntopy Analysis ----
# Author: Matěj Tvarůžka
# Description: Data analysis of BBS dataset for my master's thesis

# Check for missing packages and install them
cran_pkgs <- c("tidyverse", "clootl", "ape", "diverge", "terra", "sf")
is_installed <- cran_pkgs %in% rownames(installed.packages())
if(any(is_installed == FALSE)){
  install.packages(cran_pkgs[!is_installed])
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# To install from GitHub:
#remotes::install_github("bbsBayes/bbsBayes")
#remotes::install_github("ebird/ebirdst")
#remotes::install_github("trashbirdecology/bbsAssistant")

# Libraries
require(clootl)
require(tidyverse)
require(bbsAssistant)
require(ape)
require(ebirdst)
require(diverge)
require(terra)
require(sf)
require(ggspatial)

# Download BBS data
# bbs <- grab_bbs_data()  Using bbsAssistent package
# bbs <- fetch_bbs_data() Using bbsBayes

# Step 1: Select Species Pairs -----------------------------------
# load BBS dataset
# extract vectore of species names
# select main groups for the phylogeny
# extract the tree for selected groups
# extract sister species pairs

## Version A: selecting only several families based on some criteria(?)
names <- bbs$species_list |> 
  filter(ORDER == "Passeriformes") |> 
  group_by(Family) %>% 
  summarise(nrow(Species), n_species=n()) |> 
  arrange(desc(n_species)) |> 
  filter(n_species > 10)
families <- names$Family
which(names$n_species > 20)
sum(names[1:nrow(names),2])
# If I will select 12 most species rich passerine families in BBS
# I will end up with 297 species

# BBS names vector
bbs_names_vec <- bbs$species_list |>
  filter(Family %in% families) |> 
  pull(Scientific_Name)

# Intersection of BBS and eBirdst names vectors
names_vec <- ebirdst_runs |> 
  filter(scientific_name %in% bbs_names_vec) |> 
  pull(scientific_name)

# Passerine taxonomy names vector
taxonomy <- taxonomyGet(taxonomy_year = 2024)
tax_passerines <- taxonomy[taxonomy$ORDER1 == "Passeriformes",]
length(tax_passerines$SCI_NAME)
tax_names_vec <- tax_passerines$SCI_NAME

# Get the trees
extractTree(species = ebirst_names_vec, taxonomy_year = 2023)
treeA <- extractTree(species = names_vec, taxonomy_year = 2025, version = "1.6", force = TRUE)
plot(treeA, type = "fan", cex = 0.3, tip.color = "darkblue")
extract_sisters(treeA)

## Version B: Selecting all possible species
all_names <- bbs$species_list |> 
  filter(ORDER == "Passeriformes") |>
  pull(Scientific_Name, AOU)
length(all_names)

all_ebirdst <- ebirdst_runs |> 
  filter(scientific_name %in% all_names) |>
  pull(scientific_name)
length(all_ebirdst)
# The total number of BBS species is 352 and 298 of them have ebirdst data

# Filter data for 2018
data_2018 <- bbs$observations |> 
  filter(Year == 2018) 
names_2018 <- data_2018$AOU
names_2018 <- bbs$species_list |> 
  filter(ORDER == "Passeriformes") |>
  filter(AOU %in% names_2018) |> 
  pull(Scientific_Name)
length(names_2018)

names <- ebirdst_runs |> 
  filter(scientific_name %in% names_2018) |>
  pull(scientific_name)
length(names) # 280 species

tree_2018 <- extractTree(species = names,
                         taxonomy_year = 2023,
                         version = "1.5",
                         data_path="/home/matejtvar/Documents/02_Research/BBS-syntopy/AvesDataLite-main")
plot(tree_2018, type = "fan", cex = 0.3, tip.color = "darkblue")
all_sisters <- extract_sisters(tree_2018)
dim(all_sisters) # 94 species pairs
sp1 <- all_sisters$sp1
sp2 <- all_sisters$sp2
sister_names <- c(sp1,sp2) 
length(sister_names)
# Matching table of codes and names
ebird_lookup <- ebirdst_runs |>
  filter(scientific_name %in% sister_names) |> 
  select(species_code, scientific_name)

# An access key is required to download eBird Status and Trends data
# set_ebirdst_access_key("")
# Download species ranges
for( i in sister_names){
  ebirdst_download_status(species = i, path = ebirdst_data_dir(),
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
                          show_progress = TRUE)
}
# Make a list of ranges
ranges <- lapply(sister_names, function(species) {
  load_ranges(
    species = species,
    resolution = "27km",
    smoothed = TRUE,
    path = ebirdst_data_dir()
  )
})
extract_codes <- sapply(ranges, function(df){
  df$species_code[1]
})
names(ranges) <- extract_codes
head(ranges, 1)
head(all_sisters)

# 2. Calculate Range Overlap ----------------------------------------------

# Pairs dataframe
pairs_data <- all_sisters |> 
  left_join(ebird_lookup, by = c("sp1" = "scientific_name")) |>
  rename(code1 = species_code) |>
  left_join(ebird_lookup, by = c("sp2" = "scientific_name")) |>
  rename(code2 = species_code)
head(pairs_data)

# Run sympatry calculation across all pairs
pairs_data$range_overlap <- sapply(1:nrow(pairs_data), function(i) {
  calculate_overlap(
    pairs_data$code1[i],
    pairs_data$code2[i],
    ranges
  )
})
# Gives warning at the end:
# attribute variables are assumed to be spatially constant throughout all geometries

sympatric_pairs <- pairs_data |> 
  filter(range_overlap >= 5)
nrow(sympatric_pairs)

ggplot(sympatric_pairs, aes(x = range_overlap)) +
  geom_histogram(fill = "steelblue", color = "white", boundary = 0) +
  theme_minimal() +
  labs(
    title = "Distribution of Range Overlap in 56 Passerine Sister Pairs",
    x = "Symmetric Overlap Index (0 = Allopatric, 1 = Completely Overlapping)",
    y = "Number of Pairs"
  )

# Run symmetry calculation across all pairs
pairs_data$range_symmetry <- sapply(1:nrow(pairs_data), function(i) {
  calculate_symmetry(
    pairs_data$code1[i],
    pairs_data$code2[i],
    ranges
  )
})

# 3. Inspecting N of BBS routes in sympatric zone -------------------------

# Convert BBS routes to a spatial object (sf)
bbs_routes_sf <- st_as_sf(bbs$routes,
                          coords = c("Longitude", "Latitude"),
                          crs = 4326
)

# Filter data for 2018
observations <- bbs$observations |>
  filter(Year == 2018)

count_routes_sympatry <- function(code1, code2, range_list, routes_sf) {
  # Function to count number of routes in sympatric range
  
  # Get breeding ranges
  r1 <- range_list[[code1]][range_list[[code1]]$season %in% c("breeding", "resident"), ]
  r2 <- range_list[[code2]][range_list[[code2]]$season %in% c("breeding", "resident"), ]
  
  if (nrow(r1) == 0 | nrow(r2) == 0) {
    return(0)
  }
  
  # Find the intersection polygon
  overlap_poly <- st_intersection(st_make_valid(r1), st_make_valid(r2))
  
  if (nrow(overlap_poly) == 0) {
    return(0)
  }
  
  # Count routes that fall inside this polygon
  routes_inside <- st_intersects(routes_sf, overlap_poly, sparse = FALSE)
  return(sum(routes_inside))
}

# lookup for AOU codes using sister_pairs names
aou_lookup <- bbs$species_list |>
  select(AOU, Scientific_Name)

# Add AOU codes to sister_pairs
pairs_data <- pairs_data |>
  left_join(aou_lookup, by = c("sp1" = "Scientific_Name")) |>
  rename(aou1 = AOU) |>
  left_join(aou_lookup, by = c("sp2" = "Scientific_Name")) |>
  rename(aou2 = AOU)
names(pairs_data)

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

data <- pairs_data |>
  filter(range_overlap >= 5) |>
  mutate(
    sympatry_routes = map2_dbl(code1, code2, ~ count_routes_sympatry(.x, .y, ranges, bbs_routes_sf)),
    bbs_shared_routes = map2_dbl(aou1, aou2, ~ count_shared_routes(.x, .y, observations))
  )
summary(data)
nrow(data)
View(data)
hist(data$bbs_shared_routes)

# Save the object
saveRDS(data, "syntopy_data_2018.rds")
write_csv(data, "syntopy_data_2018.csv")
# data <- readRDS("sister_pairs_data.rds")