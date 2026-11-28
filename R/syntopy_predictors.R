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
library(gawdis)
library(FD)
library(factoextra)
library(gridExtra)
library(corrplot)

# 1. Syntopy at route scale ----------------------------------------

# Selecting year
year <- 2019

dat_year <- dat |>
  dplyr::filter(Year == year) |>
  sf::st_as_sf(crs = 4326)

source(here::here("R/functions/calculate_route_counts.R"))

# Map calculate_route_counts across all sympatric sister pairs
route_counts <- map2_dfr(
  sympatric_pairs$sp1,
  sympatric_pairs$sp2,
  \(s1, s2) {
    calculate_route_counts(
      sp1 = s1,
      sp2 = s2,
      range_list = ranges,
      bbs_dat = dat_year,
      season_filter = c("breeding", "resident")
    )
  }
)

route_data <- bind_cols(sympatric_pairs, route_counts)

# Filter out pairs with no routes sampled (N == 0) for modeling
route_data <- route_data |>
  filter(N > 0)
# Note: for 2019 and 2017 there is just 1 pair with N=0, while for 2018 there are 2 pairs

# Calculate co-occurrence indicies for route counts
calc_log_OR <- function(N, m1, m2, k){
  a <- k
  b <- m1 - k
  c <- m2 - k
  d <- N - m1 - m2 + k
  # Haldane-Anscombe 0.5 continuity correction for empirical log-odds ratio
  empirical_log_or <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  return(empirical_log_or)
}

route_data_indicies <- route_data |>
  dplyr::mutate(log_OR = purrr::pmap_dbl(list(N, m1, m2, k), calc_log_OR))

# Calculate affinity using CooccurrenceAffinity package
aff <- purrr::pmap(route_data[c("N", "m1", "m2", "k")], \(N, m1, m2, k) CooccurrenceAffinity::affinity2by2(x = k, marg = c(m1, m2, N)))
# Unlist into a matrix
aff_mat <- do.call(rbind, lapply(aff, \(x) unlist(x[c("est", "jaccard", "sorensen", "simpson")])))
route_data_indicies <- cbind(route_data_indicies, as.data.frame(aff_mat))
route_data_indicies <- route_data_indicies |>
  dplyr::rename(alpha_mle = est)

# 2. Syntopy at stop scale ------------------------------------------------
syntopic_routes <- route_data |>
  dplyr::filter(k > 0) |> # only pairs co-occurring at some routes can be used
  dplyr::select(sp1, sp2)

# Source the function
source(here::here("R/functions/calculate_stop_counts.R"))

# Map across sympatric pairs
stop_counts <- map2_dfr(
  syntopic_routes$sp1,
  syntopic_routes$sp2,
  \(s1, s2) {
    calculate_stop_counts(
      sp1 = s1,
      sp2 = s2,
      bbs_dat = dat_year
    )
  }
)
# Calculate co-occurrence indicies for stop counts
stop_data_indicies <- stop_counts |>
  dplyr::filter(m1 > 0 & m2 > 0) |>
  dplyr::mutate(log_OR = purrr::pmap_dbl(list(N, m1, m2, k), calc_log_OR))

# Calculate affinity using CooccurrenceAffinity package
aff <- purrr::pmap(stop_data_indicies[c("N", "m1", "m2", "k")], \(N, m1, m2, k) CooccurrenceAffinity::affinity2by2(x = k, marg = c(m1, m2, N)))
# Unlist into a matrix
aff_mat <- do.call(rbind, lapply(aff, \(x) unlist(x[c("est", "jaccard", "sorensen", "simpson")])))
stop_data_indicies <- cbind(stop_data_indicies, as.data.frame(aff_mat))
stop_data_indicies <- stop_data_indicies |>
  dplyr::rename(alpha_mle = est)

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
  dplyr::select(RTENO, Year, NDVI_filled) |>
  dplyr::filter(Year == year)

stop_data <- stop_counts |>
  dplyr::left_join(ndvi_yr, by = c("RTENO"))

# 5. Traits and niche similarity --------------------------------

# Calculate Gower distance for Avoniche traits
avoniche_names <- avoniche$Species2

dietary_categories <- avoniche |>
  dplyr::select(-Family2, -Inferred, -MainDiet, -MainNiche) |>
  dplyr::select(Species2, Pa, Pt, Ne, Se, Fr, In, Ap, Vt, Ca) |>
  tibble::column_to_rownames(var = "Species2")

foraging_strategy <- avoniche |>
  dplyr::select(-Family2, -Inferred, -MainDiet, -MainNiche, -Pa, -Pt, -Ne, -Se, -Fr, -In, -Ap, -Vt, -Ca) |>
  tibble::column_to_rownames(var = "Species2")

gd_diet <- FD::gowdis(dietary_categories)
gd_fs <-  FD::gowdis(foraging_strategy)

diet_matrix <- as.matrix(gd_diet)
colnames(diet_matrix)[1:192] <- paste(avoniche_names)
rownames(diet_matrix)[1:192] <- paste(avoniche_names)

fs_matrix <- as.matrix(gd_fs)
colnames(fs_matrix)[1:192] <- paste(avoniche_names)
rownames(fs_matrix)[1:192] <- paste(avoniche_names)

# Calculate PCA difference for selected traits

# Gower habitaty

traits <- c("Beak.Width", "Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Depth")

m <- avonet  |>
  select(Species2, all_of(traits)) |>
  pivot_longer(
    cols = all_of(traits),
    names_to = "trait",
    values_to = "value"
  ) |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "white", na.rm = TRUE) +
  facet_wrap(~ trait, scales = "free") +
  theme_minimal() +
  labs(
    title = "Distribution of Selected Traits Across Sister Species",
    x = "Trait Measurement",
    y = "Count"
  )
# m

traits_df <- avonet |>
  dplyr::select(Species2, all_of(traits), pair_id, role)
summary(traits_df)

res.pca <- prcomp(log(traits_df[, traits]))
print(res.pca)
summary(res.pca)

# Visualise eigenvalues
eig.val<-get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, col.var="lightblue")
# Results for variables
var <- get_pca_var(res.pca)
var
head(var$cos2) # quality of representation of variables
# corrplot(var$cos2, is.corr=FALSE)
# fviz_pca_var(res.pca,
#              col.var = "cos2", # Color by the quality of representation
#              gradient.cols = c("darkorchid4", "gold", "darkorange"),
#              repel = TRUE
# )


# Contributions of variables to PC1
a<-fviz_contrib(res.pca, choice = "var", axes = 1)
# Contributions of variables to PC2
b<-fviz_contrib(res.pca, choice = "var", axes = 2)
grid.arrange(a,b, ncol=2, top='Contribution of the variables to the first two PCs')

# PCA results for species pairs
# add PCA3
# body mass = difference in logs of mass
sister_pca_diff <- cbind(traits_df[, c("pair_id", "Species2")], res.pca$x[, 1:2]) |>
  group_by(pair_id) |>
  summarise(sp1 = first(Species2),
            sp2 = last(Species2),
            PC1_diff = abs(diff(PC1)),
            PC2_diff = abs(diff(PC2)),
            PCA_Euclidean_Dist = sqrt(diff(PC1)^2 + diff(PC2)^2),  # funkce dist vegan?
            .groups = "drop")

list(
  route_data = route_data,
  stop_data = stop_data,
  range_metrics = range_metrics,
  route_data_indicies = route_data_indicies,
  stop_data_indicies = stop_data_indicies
)
