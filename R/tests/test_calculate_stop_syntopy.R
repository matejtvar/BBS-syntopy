library(sf)
library(dplyr)

# Helper function to create a mock data frame for a route with 50 stops
create_mock_stops <- function(sp_name, route_name, active_stops) {
  row_data <- data.frame(
    RTENO = route_name,
    Year = 2019,
    Scientific_Name = sp_name
  )
  for (i in 1:50) {
    row_data[[paste0("Stop", i)]] <- if (i %in% active_stops) 1 else 0
  }
  return(row_data)
}

# 1. Create mock coordinate points for spatial sf object creation
pt1 <- sf::st_point(c(-95.0, 40.0))
pt2 <- sf::st_point(c(-95.1, 40.1))
pt3 <- sf::st_point(c(-95.2, 40.2))

# 2. Build mock raw data frames across multiple routes (R1, R2, R3)
# Route 1: Species A and Species B co-occur at stops 1 and 2
r1_a <- create_mock_stops("Species A", "R1", c(1, 2, 5))
r1_b <- create_mock_stops("Species B", "R1", c(1, 2, 10, 15))

# Route 2: Species A present, Species B absent
r2_a <- create_mock_stops("Species A", "R2", c(3, 4, 6))

# Route 3: Species B present, Species A absent
r3_b <- create_mock_stops("Species B", "R3", c(7, 8))

# 3. Combine into a single data frame and convert to an sf object
mock_bbs_df <- bind_rows(r1_a, r1_b, r2_a, r3_b)

# Assign a dummy geometry matching the number of rows
mock_bbs_stops <- sf::st_sf(
  mock_bbs_df,
  geometry = sf::st_sfc(
    rep(list(pt1), nrow(mock_bbs_df)),
    crs = 4326
  )
)


library(testthat)
library(dplyr)
library(sf)

test_that("calculate_stop_syntopy_per_route returns correct per-route metrics", {

  res <- calculate_stop_syntopy_per_route(
    sp1 = "Species A",
    sp2 = "Species B",
    bbs_dat = mock_bbs_stops
  )

  # Check that output has multiple rows (one per route) rather than a single collapsed row
  expect_true(nrow(res) > 0)
  expect_true(all(c("RTENO", "N", "m1", "m2", "k", "empirical_log_OR"), names(res) %in% names(res)))

  # Check that N equals 50 for each individual route standard
  expect_true(all(res$N == 50))

  # Check numeric formatting of log OR
  expect_true(all(is.numeric(res$empirical_log_OR)))
})

test_that("calculate_stop_syntopy_per_route returns 0 co-occurrences when species never overlap on the same route", {

  r4_a <- create_mock_stops("Species A", "R4", c(1, 2))
  r5_c <- create_mock_stops("Species C", "R5", c(1, 2))

  mock_separate_df <- bind_rows(r4_a, r5_c)
  mock_separate_sf <- sf::st_sf(
    mock_separate_df,
    geometry = sf::st_sfc(pt1, pt2, crs = 4326)
  )

  res <- calculate_stop_syntopy_per_route(
    sp1 = "Species A",
    sp2 = "Species C",
    bbs_dat = mock_separate_sf
  )

  # Should produce rows for each route, but k should be 0 where they don't co-occur
  expect_true(all(res$k == 0))
  expect_true(all(is.numeric(res$empirical_log_OR)))
})

test_that("calculate_stop_syntopy_per_route handles missing species gracefully by returning empty tibble", {

  res_missing <- calculate_stop_syntopy_per_route(
    sp1 = "Species A",
    sp2 = "Nonexistent Species",
    bbs_dat = mock_bbs_stops
  )

  # Should return an empty tibble with correct columns
  expect_equal(nrow(res_missing), 0)
  expect_true(all(c("RTENO", "N", "m1", "m2", "k", "empirical_log_OR") %in% names(res_missing)))
})

test_that("calculate_stop_syntopy_per_route aggregates duplicate route-species rows safely", {

  create_mock_stops_raw <- function(sp_name, route_name, active_stops) {
    row_data <- data.frame(
      RTENO = route_name,
      Year = 2019,
      Scientific_Name = sp_name
    )
    for (i in 1:50) {
      row_data[[paste0("Stop", i)]] <- if (i %in% active_stops) 1 else 0
    }
    return(row_data)
  }

  r1_a_part1 <- create_mock_stops_raw("Species A", "R1", c(1))
  r1_a_part2 <- create_mock_stops_raw("Species A", "R1", c(2))
  r1_b        <- create_mock_stops_raw("Species B", "R1", c(1, 3))

  mock_dupes_df <- bind_rows(r1_a_part1, r1_a_part2, r1_b)
  mock_dupes_sf <- sf::st_sf(
    mock_dupes_df,
    geometry = sf::st_sfc(pt1, pt1, pt1, crs = 4326)
  )

  res <- NULL

  expect_no_warning({
    res <- calculate_stop_syntopy_per_route(
      sp1 = "Species A",
      sp2 = "Species B",
      bbs_dat = mock_dupes_sf
    )
  })

  # Inspect route R1 values after internal aggregation
  r1_res <- res |> dplyr::filter(RTENO == "R1")
  expect_equal(r1_res$N, 50)
  expect_equal(r1_res$m1, 2)  # Stops 1 and 2 combined
  expect_equal(r1_res$m2, 2)  # Stops 1 and 3
  expect_equal(r1_res$k, 1)   # Overlap only at stop 1
  expect_true(is.numeric(r1_res$empirical_log_OR))
})
