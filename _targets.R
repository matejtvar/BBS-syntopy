library(targets)

tar_option_set(
  packages = c(
    "dplyr", "sf", "ebirdst", "clootl", "diverge", "purrr",
    "ggplot2", "tidyr", "ape", "here", "readxl", "tibble", "cluster"
  ),
  format = "rds"
)

list(
  # 1. Track raw input files and directories as file targets
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(taxonomy_dir, "data/AvesDataLite-main", format = "file"),
  tar_target(ebirdst_dir, ebirdst::ebirdst_data_dir(), format = "file"),

  # 2. Data Import Target: Sources R/data_import.R
  tar_target(
    data_import,
    {
      file_deps <- c(bbs_file, taxonomy_dir, ebirdst_dir)

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
      # Extract objects from data_import
      dat <- data_import$dat
      sister_pairs <- data_import$sister_pairs
      sister_names <- data_import$sister_names
      ranges <- data_import$ranges


      source(here::here("R/syntopy_predictors.R"))

      sympatric_pairs
    }
  )
)
