library(testthat)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)

# Source the custom stop syntopy function
source(here::here("R/functions/calculate_stop_syntopy.R"))

# 1. Setup Mock Data for Standard Tests -----------------------------------
pt1 <- sf::st_point(c(6, 5))
pt2 <- sf::st_point(c(7, 5))
pt3 <- sf::st_point(c(8, 5))

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

r1_a <- create_mock_stops("Species A", "R1", c(1, 2, 3, 4))
r1_b <- create_mock_stops("Species B", "R1", c(1, 2, 5))
r2_a <- create_mock_stops("Species A", "R2", c(10))
r3_b <- create_mock_stops("Species B", "R3", c(10))

mock_bbs_stops_df <- bind_rows(r1_a, r1_b, r2_a, r3_b)

mock_bbs_stops <- sf::st_sf(
  mock_bbs_stops_df,
  geometry = sf::st_sfc(
    pt1, pt1, pt2, pt3,
    crs = 4326
  )
)


# 2. Unit Tests -----------------------------------------------------------

test_that("calculate_stop_syntopy extracts correct stop-level Hypergeometric parameters", {

  res <- calculate_stop_syntopy(
    sp1 = "Species A",
    sp2 = "Species B",
    bbs_dat = mock_bbs_stops
  )

  expect_equal(res$N, 150)
  expect_equal(res$m1, 5)
  expect_equal(res$m2, 4)
  expect_equal(res$k, 2)
  expect_true(is.numeric(res$empirical_log_OR))
})

test_that("calculate_stop_syntopy returns 0s when species never appear together", {

  r4_a <- create_mock_stops("Species A", "R4", c(1, 2))
  r5_c <- create_mock_stops("Species C", "R5", c(1, 2))

  mock_separate_df <- bind_rows(r4_a, r5_c)
  mock_separate_sf <- sf::st_sf(
    mock_separate_df,
    geometry = sf::st_sfc(pt1, pt2, crs = 4326)
  )

  res <- calculate_stop_syntopy(
    sp1 = "Species A",
    sp2 = "Species C",
    bbs_dat = mock_separate_sf
  )

  expect_equal(res$k, 0)
  expect_true(is.numeric(res$empirical_log_OR))
})

test_that("calculate_stop_syntopy handles missing species gracefully", {

  res_missing <- calculate_stop_syntopy(
    sp1 = "Species A",
    sp2 = "Nonexistent Species",
    bbs_dat = mock_bbs_stops
  )

  expect_equal(res_missing$N, 0)
  expect_equal(res_missing$m2, 0)
  expect_equal(res_missing$k, 0)
  expect_true(is.na(res_missing$empirical_log_OR))
})

test_that("calculate_stop_syntopy aggregates duplicate rows and runs without warnings", {

  # Helper to create raw stop rows for duplicate tests
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
  r1_b       <- create_mock_stops_raw("Species B", "R1", c(1, 3))

  mock_dupes_df <- bind_rows(r1_a_part1, r1_a_part2, r1_b)
  mock_dupes_sf <- sf::st_sf(
    mock_dupes_df,
    geometry = sf::st_sfc(pt1, pt1, pt1, crs = 4326)
  )

  res <- NULL

  expect_no_warning({
    res <- calculate_stop_syntopy(
      sp1 = "Species A",
      sp2 = "Species B",
      bbs_dat = mock_dupes_sf
    )
  })

  expect_equal(res$N, 50)
  expect_equal(res$m1, 2)
  expect_equal(res$m2, 2)
  expect_equal(res$k, 1)
  expect_true(is.numeric(res$empirical_log_OR))
})
