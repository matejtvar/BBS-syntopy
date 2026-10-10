# BBS-syntopy
Here I focus on building a conceptual tool for my thesis, which will analyse local co-occcurence (syntopy) on two spatial scales (between and withing transects) using data from Breeding Bird Survey.


BBS-syntopy/
├── _targets.R
├── R/
│   ├── 01_data_import.R          # Functions for reading & cleaning BBS + eBird + phylogeny
│   ├── 02_syntopy_predictors.R   # Functions for calculating spatial & phylogenetic predictors
│   └── functions/
│       ├── extract_age.R         # extract_sister_ages()
│       ├── calculate_sympatry.R  # calculate_sympatry()
│       ├── calculate_symmetry.R # calculate_symmetry()
│       └── calculate_cen_dist.R # calculate_centroid_distance()
