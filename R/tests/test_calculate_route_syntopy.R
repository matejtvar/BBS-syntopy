library(testthat)
library(sf)
library(dplyr)
library(purrr)

# Source the custom function
source(here::here("R/functions/calculate_route_syntopy.R"))

# 1. Create Artificial Ranges (Polygons) -----------------------------------
poly_a_coords <- matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol = 2, byrow = TRUE)
poly_b_coords <- matrix(c(5,0, 15,0, 15,10, 5,10, 5,0), ncol = 2, byrow = TRUE) # Overlaps x=5 to x=10 with A
poly_c_coords <- matrix(c(20,20, 30,20, 30,30, 20,30, 20,20), ncol = 2, byrow = TRUE) # No overlap

poly_a <- sf::st_make_valid(sf::st_polygon(list(poly_a_coords)))
poly_b <- sf::st_make_valid(sf::st_polygon(list(poly_b_coords)))
poly_c <- sf::st_make_valid(sf::st_polygon(list(poly_c_coords)))

mock_ranges <- list(
  "Species A" = sf::st_sf(season = "breeding", geometry = sf::st_sfc(poly_a, crs = 4326)),
  "Species B" = sf::st_sf(season = "breeding", geometry = sf::st_sfc(poly_b, crs = 4326)),
  "Species C" = sf::st_sf(season = "breeding", geometry = sf::st_sfc(poly_c, crs = 4326))
)

# 2. Create Artificial BBS Routes inside/outside Overlap Zone --------------
pt1 <- sf::st_point(c(6, 5))  # Overlap zone: Sp A + Sp B present (both)
pt2 <- sf::st_point(c(7, 5))  # Overlap zone: Sp A only
pt3 <- sf::st_point(c(8, 5))  # Overlap zone: Sp B only
pt4 <- sf::st_point(c(9, 5))  # Overlap zone: Neither species detected
pt5 <- sf::st_point(c(2, 5))  # Outside overlap (Species A only range)

mock_bbs <- sf::st_sf(
  RTENO = c("R1", "R1", "R2", "R3", "R4", "R5"),
  Year = rep(2019, 6),
  Scientific_Name = c("Species A", "Species B", "Species A", "Species B", "None", "Species A"),
  geometry = sf::st_sfc(pt1, pt1, pt2, pt3, pt4, pt5, crs = 4326)
)

mock_pairs <- tibble::tibble(sp1 = "Species A", sp2 = "Species B")

# 3. Unit Tests -----------------------------------------------------------

test_that("calculate_route_syntopy extracts correct Hypergeometric parameters", {

  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE), add = TRUE)

  res <- calculate_route_syntopy(
    sp1 = "Species A",
    sp2 = "Species B",
    range_list = mock_ranges,
    bbs_dat = mock_bbs,
    season_filter = "breeding"
  )

  # Total overlap routes (R1, R2, R3, R4) -> N = 4
  expect_equal(res$N, 4)

  # Species A present on R1, R2 -> m1 = 2
  expect_equal(res$m1, 2)

  # Species B present on R1, R3 -> m2 = 2
  expect_equal(res$m2, 2)

  # Both co-occur on R1 -> k = 1
  expect_equal(res$k, 1)

  # Check empirical log-odds ratio output is a valid numeric scalar
  expect_true(is.numeric(res$empirical_log_OR))
})

test_that("calculate_route_syntopy returns 0s for non-overlapping ranges", {

  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE), add = TRUE)

  res <- calculate_route_syntopy(
    sp1 = "Species A",
    sp2 = "Species C",
    range_list = mock_ranges,
    bbs_dat = mock_bbs,
    season_filter = "breeding"
  )

  expect_equal(res$N, 0)
  expect_equal(res$k, 0)
  expect_true(is.na(res$empirical_log_OR))
})
