
# Sister pairs extractio from phylogeny -----------------------------------

# Libraries
require(diverge)
require(clootl)
require(tidyverse)
require(ape)

taxonomy <- taxonomyGet(taxonomy_year = 2021)
glimpse(taxonomy)

warblers <- taxonomy |> 
  filter(FAMILY == "Parulidae (New World Warblers)")
glimpse(warblers)

warb_names_vec <- warblers |> 
  select(SCI_NAME) |> 
  pull(SCI_NAME)
glimpse(warb_names_vec)

# 1. Ensure your names are clean (remove trailing/leading whitespace)
warb_names_vec <- trimws(warb_names_vec)

# 2. Get the official list of scientific names from the 2021 taxonomy
# This ensures we only ask for species that clootl definitely recognizes
taxonomy <- taxonomyGet(taxonomy_year = 2021)
valid_sci_names <- taxonomy$SCI_NAME

# 3. Filter your vector (Remeš's defensive approach)
# This removes any species that might have changed names in the 2021 version
final_names <- warb_names_vec[warb_names_vec %in% valid_sci_names]

# Check how many were matched
message(paste("Matched", length(final_names), "out of", length(warb_names_vec), "warblers."))

# 4. Extract the tree
# Since it's a vector of scientific names, clootl handles it automatically
tree <- extractTree(species = final_names, taxonomy_year = 2021)
# Circular layout
plot(tree, type = "fan", cex = 0.5, tip.color = "darkblue")


# 5. Extract sister species pairs (as per Remeš MS)
# This is the 'syntopy' preparation step from the original script
sister_pairs <- extract_sisters(tree)

# Inspect the results
head(sister_pairs)





?extractTree
?extract_sisters
