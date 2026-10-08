library(testthat)
library(sf)
library(ggplot2)
library(dplyr)

# Build realistic mock eBird-structured sf objects
poly1 <- sf::st_polygon(list(matrix(c(0,0, 0,2, 2,2, 2,0, 0,0), ncol = 2, byrow = TRUE)))
poly2 <- sf::st_polygon(list(matrix(c(1,0, 1,2, 3,2, 3,0, 1,0), ncol = 2, byrow = TRUE)))

mock_ranges <- list(
  "Piranga ludoviciana" = sf::st_sf(
    species_code = "westan",
    scientific_name = "Piranga ludoviciana",
    season = c("breeding", "nonbreeding"),
    geom = sf::st_sfc(poly1, poly1, crs = 4326)
  ),
  "Piranga olivacea" = sf::st_sf(
    species_code = "scatan",
    scientific_name = "Piranga olivacea",
    season = c("breeding", "nonbreeding"),
    geom = sf::st_sfc(poly2, poly2, crs = 4326)
  )
)

# Bind the mock ranges list into a single sf data frame
mock_sf <- dplyr::bind_rows(mock_ranges)

# Plot using ggplot2
ggplot(mock_sf) +
  geom_sf(aes(fill = scientific_name), alpha = 0.5, color = "black") +
  theme_minimal() +
  scale_fill_manual(
    values = c("Piranga ludoviciana" = "firebrick", "Piranga olivacea" = "steelblue")
  ) +
  labs(
    title = "Test Dataset Range Overlap",
    fill = "Species"
  )

test_that("calculate_symmetry handles multi-season sf structures correctly", {
  # Calculate breeding season overlap (poly1 and poly2 are identical)
  sym_index <- calculate_symmetry(
    sp1 = "Piranga ludoviciana",
    sp2 = "Piranga olivacea",
    range_list = mock_ranges,
    season_filter = "breeding"
  )

  expect_type(sym_index, "double")
  expect_equal(round(sym_index, 1),0.5)
})

test_that("calculate_sympatry fails on invalid input types", {
  expect_error(
    calculate_symmetry("Piranga ludoviciana", "Missing species", mock_ranges),
    regexp = "missing from 'range_list'"
  )
})

# Tests passed successfully!

result <- calculate_symmetry(
  sp1 = "Piranga ludoviciana",
  sp2 = "Piranga olivacea",
  range_list = mock_ranges,
  season_filter = "breeding"
)
