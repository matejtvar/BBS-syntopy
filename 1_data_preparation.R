
# BBS-syntopy project outline ---------------------------------------------
# Matěj Tvarůžka

# Here I focus on building a conceptual tool for my thesis,
# which will analyse local co-occurence (syntopy) on two spatial scales (between and withing transects) using data from Breeding Bird Survey.

### Basic Info ###

# The BBS is a long-term, large-scale, international avian monitoring program initiated in 1966 to track the status and trends of North American bird populations.
# 


### Dataset:
# routes length: 25.4 miles (= ~ 39.42893 km)
# stops distance: 0.5 mile (= ~ 804.672 m)
# 3 minute point count at each stop
# record in radius of 0.25 mile (= ~ 402.336 m)
# number of transects: 4100+


# Installing bbsAssistant package from GitHub forked repository -----------------------------
# Original repository here: https://github.com/trashbirdecology/bbsAssistant

remotes::install_github("matejtvar/bbsAssistant")


# Libraries ---------------------------------------------------------------

library(bbsAssistant) # for downloading, importing, and munging the official releases of the North American Breeding Bird Survey (BBS) data
ls("package:bbsAssistant") # view functions and data in package bbsAssistant

library(tidyverse) # for data transformations
library(sf) # for spatial data manipulations
library(here) # for selecting the file path within the project

# Downloading data --------------------------------------------------------

bbs <- grab_bbs_data()
names(bbs)
View(bbs)
