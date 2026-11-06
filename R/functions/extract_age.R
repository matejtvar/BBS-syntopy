library(ape)
library(dplyr)
library(purrr)
library(phytools)

#' Extract Evolutionary Divergence Age for Sister Pairs via Cophenetic Distance
#'
#' @param sister_pairs Data frame containing sister pairs with columns `sp1` and `sp2`
#' @param phylo_tree Phylogenetic tree object of class `phylo`
#'
#' @return Data frame `sister_pairs` with an added `age_myr` column

extract_sister_ages <- function(sister_pairs, phylo_tree) {

  coph_matrix <- ape::cophenetic.phylo(phylo_tree)

  ages <- purrr::map2_dbl(
    sister_pairs$sp1,
    sister_pairs$sp2,
    function(sp1, sp2) {
      # Check if both species exist in the tree tips
      if (!(sp1 %in% rownames(coph_matrix)) || !(sp2 %in% rownames(coph_matrix))) {
        return(NA_real_)
      }

      # Extract divergence age
      div_age <- coph_matrix[sp1, sp2] / 2
      return(as.numeric(div_age))
    }
  )

  # 4. Bind back to dataframe
  sister_pairs <- sister_pairs |>
    dplyr::mutate(age_myr = ages)

  return(sister_pairs)
}
