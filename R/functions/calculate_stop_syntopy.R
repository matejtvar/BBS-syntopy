library(sf)
library(dplyr)
library(tidyr)
library(assertthat)

#' Calculate Stop-Level Counts for Micro-Scale Syntopy Modeling
#'
#' Extracts stop-level N, m1, m2, and k across 50 BBS stops on routes
#' where at least one of the sister species was observed, preparing inputs
#' for micro-scale co-occurrence modeling.
#'
#' @param sp1 Scientific name of species 1
#' @param sp2 Scientific name of species 2
#' @param bbs_dat Spatial sf object of BBS observations containing RTENO, Year,
#'   Scientific_Name, geometry, and Stop1 through Stop50 columns.
#'
#' @return A 1-row tibble with N, m1, m2, k, and empirical_log_OR
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
    return(tibble::tibble(N = 0, m1 = 0, m2 = 0, k = 0, empirical_log_OR = NA_real_))
  }

  # 2. Create a complete grid of Route x Stop for the relevant routes (50 stops per route)
  stop_names <- paste0("Stop", 1:50)
  base_grid <- tidyr::expand_grid(
    RTENO = relevant_routes,
    Stop_Name = stop_names
  )

  # 3. Extract observations for sp1, reshape, and aggregate by route-stop to prevent duplicates
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

  # 4. Extract observations for sp2, reshape, and aggregate by route-stop to prevent duplicates
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

  # 5. Merge with base grid cleanly (1-to-1 matching ensured via group_by/summarise)
  stop_comparison <- base_grid |>
    dplyr::left_join(obs_sp1, by = c("RTENO", "Stop_Name")) |>
    dplyr::left_join(obs_sp2, by = c("RTENO", "Stop_Name")) |>
    dplyr::mutate(
      sp1_pres = !is.na(Count_sp1) & Count_sp1 > 0,
      sp2_pres = !is.na(Count_sp2) & Count_sp2 > 0
    )

  # 6. Calculate stop-level contingency counts
  stop_n  <- nrow(stop_comparison)
  stop_m1 <- sum(stop_comparison$sp1_pres)
  stop_m2 <- sum(stop_comparison$sp2_pres)
  stop_k  <- sum(stop_comparison$sp1_pres & stop_comparison$sp2_pres)

  # Contingency matrix entries:
  a <- stop_k
  b <- stop_m1 - stop_k
  c <- stop_m2 - stop_k
  d <- stop_n - stop_m1 - stop_m2 + stop_k

  # Haldane-Anscombe 0.5 continuity correction for stop-level log odds ratio
  stop_log_or <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))

  return(tibble::tibble(
    N = stop_n,
    m1 = stop_m1,
    m2 = stop_m2,
    k = stop_k,
    empirical_log_OR = stop_log_or
  ))
}
