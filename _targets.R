library(targets)

tar_option_set(
  packages = c("dplyr", "sf", "ebirdst", "clootl", "diverge", "purrr", "ggplot2", "tidyr", "ape", "here"),
  format = "rds"
)

list(
  # 1. Track the raw data file, external directories, and custom function files as dependencies
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(taxonomy_dir, "data/AvesDataLite-main", format = "file"),
  tar_target(ebirdst_dir, ebirdst::ebirdst_data_dir(), format = "file"),
  tar_target(fn_extract_age, "R/functions/extract_age.R", format = "file"),
  tar_target(fn_sympatry, "R/functions/calculate_sympatry.R", format = "file"),
  tar_target(fn_symmetry, "R/functions/calculate_symmetry.R", format = "file"),
  tar_target(fn_cen_dist, "R/functions/calculate_cen_dist.R", format = "file"),
  tar_target(fn_route_syntopy, "R/functions/calculate_route_syntopy.R", format = "file"),

  # 2. Data Import Target: Sources R/data_import.R
  tar_target(
    data_import,
    {
      # Mention function, dataset, and directory dependencies so targets watches them for changes
      dummy_deps <- c(bbs_file, taxonomy_dir, ebirdst_dir, fn_extract_age)

      # Execute your data import script directly
      source(here::here("R/data_import.R"))

      # Return the objects created by your script to cache them
      list(
        dat = dat,
        sister_pairs = sister_pairs,
        ranges = ranges,
        sister_names = sister_names
      )
    }
  ),

  # 3. Syntopy Predictors Target: Sources R/syntopy_predictors.R
  tar_target(
    syntopy_predictors,
    {
      # Mention script and helper function dependencies so targets watches them
      dummy_deps <- c(fn_sympatry, fn_symmetry, fn_cen_dist, fn_route_syntopy)

      dat <- data_import$dat
      sister_pairs <- data_import$sister_pairs
      ranges <- data_import$ranges

      source(here::here("R/syntopy_predictors.R"))

      sympatric_pairs
    }
  )
)
