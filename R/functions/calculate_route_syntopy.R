library(sf)
library(dplyr)
library(assertthat)

#' Calculate Route-Level Counts for Fisher's Hypergeometric Model
#'
#' Extracts N, m1, m2, and k within the sympatric overlap zone for a sister pair,
#' preparing inputs for Bayesian Fisher's Non-Central Hypergeometric modeling.
#'
#' @param sp1 Scientific name of species 1
#' @param sp2 Scientific name of species 2
#' @param range_list Named list of sf range objects
#' @param bbs_dat Spatial sf object of BBS observations (RTENO, Year, Scientific_Name, geometry)
#' @param season_filter Seasons to include (default: c("breeding", "resident"))
#'
#' @return A 1-row tibble with N, m1, m2, k, and empirical_log_OR
calculate_route_syntopy <- function(sp1, sp2, range_list, bbs_dat, season_filter = c("breeding", "resident")) {

  # --- Assertions ---
  assert_that(is.character(sp1) && length(sp1) == 1 && nchar(sp1) > 0,
              msg = "sp1 must be a non-empty character string.")
  assert_that(is.character(sp2) && length(sp2) == 1 && nchar(sp2) > 0,
              msg = "sp2 must be a non-empty character string.")
  assert_that(is.list(range_list) && !is.null(names(range_list)),
              msg = "range_list must be a named list.")
  assert_that(sp1 %in% names(range_list), msg = paste("Species missing from range_list:", sp1))
  assert_that(sp2 %in% names(range_list), msg = paste("Species missing from range_list:", sp2))
  assert_that(inherits(bbs_dat, "sf"), msg = "bbs_dat must be an sf object.")
  assert_that(all(c("RTENO", "Year", "Scientific_Name") %in% names(bbs_dat)),
              msg = "bbs_dat must contain RTENO, Year, and Scientific_Name columns.")

  # 1. Extract ranges and filter seasons
  r1_raw <- range_list[[sp1]]
  r2_raw <- range_list[[sp2]]

  r1 <- if ("season" %in% names(r1_raw)) dplyr::filter(r1_raw, season %in% season_filter) else r1_raw
  r2 <- if ("season" %in% names(r2_raw)) dplyr::filter(r2_raw, season %in% season_filter) else r2_raw

  if (nrow(r1) == 0 || nrow(r2) == 0) {
    return(tibble::tibble(N = 0, m1 = 0, m2 = 0, k = 0, empirical_log_OR = NA_real_))
  }

  # 2. Compute spatial overlap polygon
  g1 <- sf::st_make_valid(sf::st_union(r1))
  g2 <- sf::st_make_valid(sf::st_union(r2))

  overlap_poly <- suppressWarnings(sf::st_intersection(g1, g2))

  if (length(overlap_poly) == 0 || sf::st_is_empty(overlap_poly)) {
    return(tibble::tibble(N = 0, m1 = 0, m2 = 0, k = 0, empirical_log_OR = NA_real_))
  }

  # Match CRS
  if (sf::st_crs(bbs_dat) != sf::st_crs(overlap_poly)) {
    bbs_dat <- sf::st_transform(bbs_dat, sf::st_crs(overlap_poly))
  }

  # 3. Filter routes inside spatial overlap polygon
  routes_unique <- bbs_dat |>
    dplyr::select(RTENO, geometry) |>
    dplyr::distinct(RTENO, .keep_all = TRUE)

  routes_in_overlap <- suppressWarnings(
    routes_unique[sf::st_intersects(routes_unique, overlap_poly, sparse = FALSE)[, 1], ]
  )

  # Total sampled routes in overlap (N)
  N_routes <- nrow(routes_in_overlap)

  if (N_routes == 0) {
    return(tibble::tibble(N = 0, m1 = 0, m2 = 0, k = 0, empirical_log_OR = NA_real_))
  }

  # 4. Filter observations strictly on routes within the overlap polygon
  bbs_overlap_obs <- bbs_dat |>
    sf::st_drop_geometry() |>
    dplyr::filter(RTENO %in% routes_in_overlap$RTENO, Scientific_Name %in% c(sp1, sp2)) |>
    dplyr::select(RTENO, Scientific_Name) |>
    dplyr::distinct()

  # 5. Calculate presence counts across unique routes in overlap
  sp1_routes <- bbs_overlap_obs |> dplyr::filter(Scientific_Name == sp1) |> dplyr::pull(RTENO)
  sp2_routes <- bbs_overlap_obs |> dplyr::filter(Scientific_Name == sp2) |> dplyr::pull(RTENO)

  m1 <- length(sp1_routes)                     # Routes with Species 1
  m2 <- length(sp2_routes)                     # Routes with Species 2
  k  <- length(intersect(sp1_routes, sp2_routes)) # Routes with BOTH (Co-occurrence)

  # Contingency matrix entries for empirical log OR calculation:
  # a = k (both present)
  # b = m1 - k (only sp1)
  # c = m2 - k (only sp2)
  # d = N - m1 - m2 + k (both absent)
  a <- k
  b <- m1 - k
  c <- m2 - k
  d <- N_routes - m1 - m2 + k

  # Haldane-Anscombe 0.5 continuity correction for empirical log-odds ratio
  empirical_log_or <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))

  return(tibble::tibble(
    N = N_routes,
    m1 = m1,
    m2 = m2,
    k = k,
    empirical_log_OR = empirical_log_or
  ))
}
