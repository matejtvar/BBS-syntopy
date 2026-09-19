# ---- BBS Syntopy Analysis ----
# Author: Matěj Tvarůžka
# Description: Data analysis of syntopy predictors with BBS dataset
# GitHub repository: https://github.com/matejtvar/BBS-syntopy

here::i_am("R/data_import.R")
library(here)
here() # project path

# To install from GitHub:
#remotes::install_github("ebird/ebirdst")
#remotes::install_github("trashbirdecology/bbsAssistant")

# Libraries
library(clootl)
library(dplyr)
library(bbsAssistant)
library(ape)
library(ebirdst)
library(diverge)
library(terra)
library(sf)

