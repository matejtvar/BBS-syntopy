library(targets)

tar_option_set(
  packages = c(
    # Core & spatial
    "dplyr", "sf", "tidyr", "purrr", "tibble", "here", "readxl", "ggplot2",
    # Phylogeny, ranges & traits
    "ebirdst", "clootl", "diverge", "ape", "phytools", "cluster",
    "CooccurrenceAffinity", "gawdis", "FD",
    # Visualization & modeling
    "factoextra", "gridExtra", "corrplot", "brms"
  ),
  format = "rds"
)

list(
  # 1. Track raw input files and directories as file targets
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(taxonomy_dir, "data/AvesDataLite-main", format = "file"),
  tar_target(ebirdst_dir, ebirdst::ebirdst_data_dir(), format = "file"),
  tar_target(avoniche_file, "data/Avoniche2_eBird2021.csv", format = "file"),
  tar_target(avonet_file, "data/AVONET Supplementary dataset 1.xlsx", format = "file"),

  # 2. Track helper function scripts for each pipeline stage
  tar_target(
    import_funcs,
    c(
      "R/functions/extract_age.R",
      "R/functions/calculate_sympatry.R"
    ),
    format = "file"
  ),

  tar_target(
    predictor_funcs,
    c(
      "R/functions/calculate_route_counts.R",
      "R/functions/calculate_stop_counts.R",
      "R/functions/calculate_symmetry.R",
      "R/functions/calculate_cen_dist.R"
    ),
    format = "file"
  ),

  # 3. Data Import & Initial Filtering Target
  tar_target(
    data_import,
    {
      file_deps <- c(bbs_file, taxonomy_dir, ebirdst_dir, avoniche_file, avonet_file)
      func_deps <- import_funcs

      source(here::here("R/data_import.R"), local = TRUE)$value
    }
  ),

  # 4. Syntopy Predictors Target
  tar_target(
    syntopy_predictors,
    {
      func_deps <- predictor_funcs
      list2env(data_import, envir = environment())

      source(here::here("R/syntopy_predictors.R"), local = TRUE)$value
    }
  ),

  # 5. Fit Models Target
  tar_target(
    fit_models,
    {
      list2env(syntopy_predictors, envir = environment())

      source(here::here("R/fit_models.R"), local = TRUE)$value
    }
  )
)
