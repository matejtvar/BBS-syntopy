#' Calculate Range Symmetry Index Between Two Species
#'
#' @param sp1 Scientific name of species 1 (character)
#' @param sp2 Scientific name of species 2 (character)
#' @param range_list Named list of `sf` objects (eBird Status & Trends format)
#' @param season_filter String specifying season to evaluate (default: "breeding")
#'
#' @return Numeric value representing the range symmetry index
#' @values 0: ranges are vastly different in size, 0.5: ranges are equal in size
#' @export
calculate_symmetry <- function(sp1, sp2, range_list, season_filter = "breeding") {

  # --- Input Assertions ---
  assertthat::assert_that(
    is.character(sp1) && length(sp1) == 1 && nchar(sp1) > 0,
    msg = "Argument 'sp1' must be a non-empty character string."
  )
  assertthat::assert_that(
    is.character(sp2) && length(sp2) == 1 && nchar(sp2) > 0,
    msg = "Argument 'sp2' must be a non-empty character string."
  )
  assertthat::assert_that(
    is.list(range_list) && !is.null(names(range_list)),
    msg = "Argument 'range_list' must be a named list."
  )
  assertthat::assert_that(
    sp1 %in% names(range_list),
    msg = paste0("Species '", sp1, "' is missing from 'range_list'.")
  )
  assertthat::assert_that(
    sp2 %in% names(range_list),
    msg = paste0("Species '", sp2, "' is missing from 'range_list'.")
  )

  # Extract sf objects
  r1_raw <- range_list[[sp1]]
  r2_raw <- range_list[[sp2]]

  assertthat::assert_that(
    inherits(r1_raw, "sf"),
    msg = paste0("Range entry for '", sp1, "' is not an sf object.")
  )
  assertthat::assert_that(
    inherits(r2_raw, "sf"),
    msg = paste0("Range entry for '", sp2, "' is not an sf object.")
  )

  # Filter target season (e.g. "breeding" or "resident") if column exists
  r1 <- if ("season" %in% names(r1_raw)) {
    dplyr::filter(r1_raw, season %in% season_filter)
  } else {
    r1_raw
  }

  r2 <- if ("season" %in% names(r2_raw)) {
    dplyr::filter(r2_raw, season %in% season_filter)
  } else {
    r2_raw
  }

  # Check if geometries exist for chosen season
  if (nrow(r1) == 0 || nrow(r2) == 0) {
    warning(paste("No range polygon found for season:", season_filter))
    return(0)
  }

  # Combine multi-row seasonal polygons into a unified geometry
  geom1 <- sf::st_union(r1)
  geom2 <- sf::st_union(r2)

  # --- Surface Area Calculations ---
  areaSP1 <- sum(sf::st_area(geom1))
  areaSP2 <- sum(sf::st_area(geom2))

  if (as.numeric(areaSP1) == 0 || as.numeric(areaSP2) == 0) {
    return(0)
  }

  # --- Index Calculation ---
  index <- as.numeric(min(areaSP1, areaSP2) / (areaSP1 + areaSP2))
  return(min(max(index, 0), 0.5))
}
