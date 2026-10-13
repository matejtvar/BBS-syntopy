library(targets)

tar_option_set(
  packages = c("dplyr", "sf", "ebirdst", "clootl", "diverge", "purrr", "ggplot2", "tidyr", "ape"),
  format = "rds"
)

list(
  # Track raw data files and script files
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(script_import, "R/data_import.R", format = "file"),
  tar_target(script_predictors, "R/syntopy_predictors.R", format = "file"),

  # 1. Run data_import.R script directly
  tar_target(
    data_import_step,
    {
      # Accessing script_import ensures changes in R/data_import.R trigger a re-run
      source(script_import)
      list(dat = dat, sister_pairs = sister_pairs, sister_names = sister_names, ranges = ranges)
    }
  ),

  # 2. Run syntopy_predictors.R script directly
  tar_target(
    syntopy_predictors_step,
    {
      # Force dependency on step 1 and script changes
      data_import_outputs <- data_import_step
      source(script_predictors)
      sympatric_pairs
    }
  )
)
