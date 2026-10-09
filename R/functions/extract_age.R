library(ape)
library(dplyr)
library(purrr)

#' Extract Evolutionary Divergence Age for Sister Pairs
#'
#' @param sister_pairs Data frame containing sister pairs with columns `sp1` and `sp2`
#' @param phylo_tree Phylogenetic tree object of class `phylo`
#'
#' @return Data frame `sister_pairs` with an added `age_myr` column

extract_sister_ages <- function(sister_pairs, phylo_tree) {

  # Ensure tree is ultrametric
  if (!ape::is.ultrametric(phylo_tree)) {
    warning("Tree is not strictly ultrametric. Adjusting with chronos().")
    phylo_tree <- ape::chronos(phylo_tree)
  }

  # Calculate depth from root to all nodes/tips
  node_depths <- ape::node.depth.edgelength(phylo_tree)
  total_tree_depth <- max(node_depths)

  # Calculate divergence age for each sister pair
  ages <- purrr::map2_dbl(
    sister_pairs$sp1,
    sister_pairs$sp2,
    function(sp1, sp2) {

      # Check if both tips exist in tree
      if (!(sp1 %in% phylo_tree$tip.label) || !(sp2 %in% phylo_tree$tip.label)) {
        return(NA_real_)
      }

      # Find Most Recent Common Ancestor (MRCA) node index
      mrca_node <- ape::getMRCA(phylo_tree, c(sp1, sp2))

      # Divergence age = Total height minus MRCA node height
      mrca_height <- node_depths[mrca_node]
      divergence_age <- total_tree_depth - mrca_height

      return(as.numeric(divergence_age))
    }
  )

  sister_pairs <- sister_pairs |>
    dplyr::mutate(age_myr = ages)

  return(sister_pairs)
}
