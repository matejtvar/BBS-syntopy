library(targets)

tar_option_set(
  packages = c(
    "dplyr", "sf", "ebirdst", "clootl", "diverge", "purrr",
    "ggplot2", "tidyr", "ape", "here", "readxl", "tibble", "cluster"
  ),
  format = "rds"
)
# tar_source()

list(
  # 1. Track raw input files and directories as file targets
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(taxonomy_dir, "data/AvesDataLite-main", format = "file"),
  tar_target(ebirdst_dir, ebirdst::ebirdst_data_dir(), format = "file"),
  tar_target(avoniche_file, "data/Avoniche2_eBird2021.csv", format = "file"),
  tar_target(avonet_file, "data/AVONET Supplementary dataset 1.xlsx", format = "file"),

  # 2. Data Import & Initial Filtering Target: Sources R/data_import.R
  tar_target(
    data_import,
    {
      file_deps <- c(bbs_file, taxonomy_dir, ebirdst_dir, avoniche_file, avonet_file)

      # Execute data import script locally and capture its returned list of objects
      source(here::here("R/data_import.R"), local = TRUE)$value
    }
  ),

  # 3. Syntopy Predictors Target: Sources R/syntopy_predictors.R
  tar_target(
    syntopy_predictors,
    {
      # Unpack all objects from data_import directly into this target's execution environment
      # (Makes dat, sympatric_pairs, ranges, gower_matrix, etc., directly available)
      list2env(data_import, envir = environment())

      # Execute syntopy predictors script and return its final output list
      source(here::here("R/syntopy_predictors.R"), local = TRUE)$value
    }
  )
)
