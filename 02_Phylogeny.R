
### Here I work with the phylogeny of Clements 2021 to extract species pairs for my diploma thesis

##
# 1. Installing clootl package and getting the tree examples -----------------------

install.packages("devtools")
library(devtools)  
install_github("eliotmiller/clootl")
library(clootl)

ex1 <- extractTree(species=c("Turdus migratorius",
                             "Setophaga dominica",
                             "Setophaga ruticilla",
                             "Sitta canadensis"))

plot(ex1)
cites <- getCitations(ex1)
cites

### Getting the 2021 taxonomy
ex2 <- extractTree(species=c("Turdus migratorius",
                             "Setophaga dominica",
                             "Setophaga ruticilla",
                             "Sitta canadensis"),
                   taxonomy_year = 2021)
plot(ex2)

# 2. Getting avonet data -----------------------------------------------------
install.packages("phytools")
install.packages("readxl")
library("phytools")
library(readxl)

# Define the file URL and destination
file_url <- "https://figshare.com/ndownloader/files/34480856"
destfile <- tempfile(fileext = ".xlsx")

# Download the file
download.file(file_url, destfile, mode = "wb")

# Read the sheet that corresponds to eBird taxonomy
dat <- as.data.frame(read_excel(destfile, sheet = "AVONET2_eBird"))

# Create a column with underscores for simplicity later
dat$underscores <- sub(" ", "_", dat$Species2)

# Take a random sample of 50 species
spp <- sample(dat$Species2, 50)

# Extract a tree for these species
prunedTree <- extractTree(species=spp, label_type="scientific", taxonomy_year=2021, version="1.5")


prunedDat <- dat[dat$Species2 %in% spp,]

# Pull a vector of traits out, log transform for normality
x <- log(prunedDat$Mass)
names(x) <- prunedDat$underscores

# Plot log body mass across the phylogeny
contMap(prunedTree, x, fsize=0.5)


# 3. Selecting Parulidae --------------------------------------------------

# 1) Getting taxonomy
tax2021 <- taxonomyGet(2021)

# 2) Getting the full tree
tree <- extractTree(taxonomy_year = 2021, version = "1.5")

# Identify all species belonging to Parulidae
# Note: Use exact family name as it appears in the taxonomy
parulidae_species <- tax2021$SCI_NAME[tax2021$FAMILY == "Parulidae (New World Warblers)"]

# 3) Extract the tree using these specific species
parulidae_tree <- extractTree(species = parulidae_species, taxonomy_year = 2021)

# Plot to verify
plot(parulidae_tree, show.tip.label = FALSE)

