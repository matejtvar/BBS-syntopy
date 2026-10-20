# BBS-syntopy

Here I focus on building a conceptual tool for my thesis, which will analyse local co-occcurence (syntopy) on two spatial scales (between and withing transects) using data from Breeding Bird Survey.

# Thesis Implementation Roadmap

## Current Status & Key Scripts

- **Data Import & Cleaning:** [`R/data_import.R`](R/data_import.R)
- **Syntopy & Predictor Extraction:** [`R/syntopy_predictors.R`](R/syntopy_predictors.R)
- **Core Functions:**
  - Route Syntopy: [`R/functions/calculate_route_syntopy.R`](R/functions/calculate_route_syntopy.R)
  - Stop Syntopy: [`R/functions/calculate_stop_syntopy.R`](R/functions/calculate_stop_syntopy.R)
  - Sympatry Calculation: [`R/functions/calculate_sympatry.R`](R/functions/calculate_sympatry.R)
  - Range Symmetry Calculation: [`R/functions/calculate_symmetry.R`](R/functions/calculate_symmetry.R)
  - Extract evolutionary age: [`R/functions/extract_age.R`](R/functions/extract_age.R)
- **Pipeline Orchestration:** [`_targets.R`](_targets.R)
