#' Calculate Geographic Distance Between Species Range Centroids
#'
#' @param sp1 Scientific name of species 1 (character)
#' @param sp2 Scientific name of species 2 (character)
#' @param range_list Named list of `sf` objects
#' @param season_filter Vector specifying seasons to evaluate (default: c("breeding", "resident"))
#' @param unit Measurement unit for output distance (e.g., "km", "m"). Default: "km"
#' @return Numeric distance between centroids in the specified unit
#' @export
calculate_centroid_distance <- function(sp1, sp2, range_list,
                                        season_filter = c("breeding", "resident"),
                                        unit = "km") {

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

  # Filter seasons if column exists
  r1 <- if ("season" %in% names(r1_raw)) dplyr::filter(r1_raw, season %in% season_filter) else r1_raw
  r2 <- if ("season" %in% names(r2_raw)) dplyr::filter(r2_raw, season %in% season_filter) else r2_raw

  if (nrow(r1) == 0 || nrow(r2) == 0) {
    warning(paste("No range polygon found for species in season:", paste(season_filter, collapse = ", ")))
    return(NA_real_)
  }

  # Union seasonal geometries into a single geometry per species
  geom1 <- sf::st_union(r1)
  geom2 <- sf::st_union(r2)

  # Calculate centroids
  # suppressWarnings handles longitude/latitude centroid warnings from s2
  centroid1 <- suppressWarnings(sf::st_centroid(geom1))
  centroid2 <- suppressWarnings(sf::st_centroid(geom2))

  # Calculate great-circle / geodesic distance
  dist_raw <- sf::st_distance(centroid1, centroid2)

  # Convert units (e.g., meters to kilometers)
  dist_converted <- units::set_units(dist_raw, unit, mode = "standard")

  return(as.numeric(dist_converted))
}
