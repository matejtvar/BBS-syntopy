library(sf)
library(dplyr)
library(tidyr)
library(assertthat)

#' Calculate Route-Specific Stop-Level Counts for Micro-Scale Syntopy Modeling
#'
#' Extracts stop-level N, m1, m2, and k across 50 BBS stops separately for each route
#' where at least one of the sister species was observed.
#'
#' @param sp1 Scientific name of species 1
#' @param sp2 Scientific name of species 2
#' @param bbs_dat Spatial sf object of BBS observations containing RTENO, Year,
#'   Scientific_Name, geometry, and Stop1 through Stop50 columns.
#'
#' @return A tibble with one row per route containing RTENO, N, m1, m2, k, and empirical_log_OR
calculate_stop_syntopy <- function(sp1, sp2, bbs_dat) {

  # --- Assertions ---
  assert_that(is.character(sp1) && length(sp1) == 1 && nchar(sp1) > 0,
              msg = "sp1 must be a non-empty character string.")
  assert_that(is.character(sp2) && length(sp2) == 1 && nchar(sp2) > 0,
              msg = "sp2 must be a non-empty character string.")
  assert_that(inherits(bbs_dat, "sf"), msg = "bbs_dat must be an sf object.")
  assert_that(all(c("RTENO", "Year", "Scientific_Name") %in% names(bbs_dat)),
              msg = "bbs_dat must contain RTENO, Year, and Scientific_Name columns.")
  assert_that(any(grepl("^Stop", names(bbs_dat))),
              msg = "bbs_dat must contain Stop columns (e.g., Stop1 to Stop50).")

  # 1. Identify routes where at least one of the target species was observed
  relevant_routes <- bbs_dat |>
    sf::st_drop_geometry() |>
    dplyr::filter(Scientific_Name %in% c(sp1, sp2)) |>
    dplyr::pull(RTENO) |>
    unique()

  if (length(relevant_routes) == 0) {
    return(tibble::tibble(
      RTENO = character(0),
      N = integer(0),
      m1 = integer(0),
      m2 = integer(0),
      k = integer(0),
      empirical_log_OR = numeric(0)
    ))
  }

  # 2. Create a complete grid of Route x Stop for the relevant routes (50 stops per route)
  stop_names <- paste0("Stop", 1:50)
  base_grid <- tidyr::expand_grid(
    RTENO = relevant_routes,
    Stop_Name = stop_names
  )

  # 3. Extract observations for sp1, reshape, and aggregate by route-stop
  obs_sp1 <- bbs_dat |>
    sf::st_drop_geometry() |>
    dplyr::filter(RTENO %in% relevant_routes, Scientific_Name == sp1) |>
    tidyr::pivot_longer(
      cols = dplyr::any_of(stop_names),
      names_to = "Stop_Name",
      values_to = "Count_sp1"
    ) |>
    dplyr::group_by(RTENO, Stop_Name) |>
    dplyr::summarise(Count_sp1 = sum(Count_sp1, na.rm = TRUE), .groups = "drop")

  # 4. Extract observations for sp2, reshape, and aggregate by route-stop
  obs_sp2 <- bbs_dat |>
    sf::st_drop_geometry() |>
    dplyr::filter(RTENO %in% relevant_routes, Scientific_Name == sp2) |>
    tidyr::pivot_longer(
      cols = dplyr::any_of(stop_names),
      names_to = "Stop_Name",
      values_to = "Count_sp2"
    ) |>
    dplyr::group_by(RTENO, Stop_Name) |>
    dplyr::summarise(Count_sp2 = sum(Count_sp2, na.rm = TRUE), .groups = "drop")

  # 5. Merge with base grid and evaluate presence/absence at each stop
  stop_comparison <- base_grid |>
    dplyr::left_join(obs_sp1, by = c("RTENO", "Stop_Name")) |>
    dplyr::left_join(obs_sp2, by = c("RTENO", "Stop_Name")) |>
    dplyr::mutate(
      sp1_pres = !is.na(Count_sp1) & Count_sp1 > 0,
      sp2_pres = !is.na(Count_sp2) & Count_sp2 > 0
    )

  # 6. Calculate stop-level metrics independently for each route and attach species names
  route_summary <- stop_comparison |>
    dplyr::group_by(RTENO) |>
    dplyr::summarise(
      N = dplyr::n(),
      m1 = sum(sp1_pres),
      m2 = sum(sp2_pres),
      k = sum(sp1_pres & sp2_pres),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      sp1 = sp1,
      sp2 = sp2,
      a = k,
      b = m1 - k,
      c = m2 - k,
      d = N - m1 - m2 + k,
      empirical_log_OR = log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
    ) |>
    dplyr::select(sp1, sp2, RTENO, N, m1, m2, k, empirical_log_OR)

  return(route_summary)
}
