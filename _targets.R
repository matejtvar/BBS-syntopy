library(targets)

# 1. Global pipeline settings
tar_option_set(
  packages = c("dplyr", "sf", "ebirdst", "clootl", "diverge", "purrr", "ggplot", "tidyr"),
  format = "rds"
)

# 2. Source main workflow scripts AND subfolder functions
# Option A: If functions/ is inside R/ (e.g. R/functions/)
tar_source("R")

# Option B: If functions/ is at the root project directory level
# tar_source(c("R", "functions"))


# 3. Define the Targets pipeline
list(
  # Track raw files
  tar_target(bbs_file, "data/bbs_data/bbs_dataset.RData", format = "file"),
  tar_target(tree_path, "data/AvesDataLite-main", format = "file"),

  # 1. Data Import
  tar_target(bbs_data, load_bbs_data(bbs_file)),
  tar_target(sister_names, get_sister_names(bbs_data, tree_path)),
  tar_target(ranges, load_species_ranges(sister_names)),

  # 2. Syntopy Predictors (uses calculate_symmetry() and calculate_sympatry())
  tar_target(
    syntopy_data,
    compute_syntopy_predictors(ranges, bbs_data) # Function defined in syntopy_predictors.R
  ),

  # 3. Modelling
  tar_target(models, run_models(syntopy_data)), # Defined in modelling.R

  # 4. Final Analysis
  tar_target(final_outputs, run_final_analysis(models)) # Defined in final_analysis.R
)
