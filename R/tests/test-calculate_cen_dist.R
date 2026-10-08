library(ggplot2)
library(sf)
library(dplyr)

plot_centroid_distance <- function(sp1, sp2, range_list, season_filter = c("breeding", "resident")) {

  # 1. Safely extract raw sf objects
  r1_raw <- range_list[[sp1]]
  r2_raw <- range_list[[sp2]]

  # 2. Filter season if column exists, otherwise use raw
  r1_sub <- if ("season" %in% names(r1_raw)) dplyr::filter(r1_raw, season %in% season_filter) else r1_raw
  r2_sub <- if ("season" %in% names(r2_raw)) dplyr::filter(r2_raw, season %in% season_filter) else r2_raw

  # Assert non-empty geometries
  if (nrow(r1_sub) == 0 || nrow(r2_sub) == 0) {
    stop(paste("No ranges found for season(s):", paste(season_filter, collapse = ", ")))
  }

  # 3. Union geometries
  g1 <- sf::st_union(r1_sub)
  g2 <- sf::st_union(r2_sub)

  # 4. Calculate centroids
  c1 <- suppressWarnings(sf::st_centroid(g1))
  c2 <- suppressWarnings(sf::st_centroid(g2))

  # 5. Build line from raw coordinates (avoids sfg/sfc class mismatch)
  coords1 <- sf::st_coordinates(c1)[1, 1:2]
  coords2 <- sf::st_coordinates(c2)[1, 1:2]

  line_sfg <- sf::st_linestring(rbind(coords1, coords2))
  line_sf  <- sf::st_sf(geometry = sf::st_sfc(line_sfg, crs = sf::st_crs(r1_raw)))

  # Distance readout
  dist_km <- round(as.numeric(units::set_units(sf::st_distance(c1, c2), "km")), 1)

  # 6. Assemble spatial data frames for plotting
  ranges_sf <- sf::st_sf(
    species = c(sp1, sp2),
    geometry = c(g1, g2)
  )

  centroids_sf <- sf::st_sf(
    species = c(sp1, sp2),
    geometry = c(c1, c2)
  )

  # 7. Render plot
  ggplot() +
    geom_sf(data = ranges_sf, aes(fill = species), alpha = 0.35, color = NA) +
    geom_sf(data = line_sf, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_sf(data = centroids_sf, aes(color = species), size = 4) +
    scale_fill_manual(values = c("firebrick", "steelblue")) +
    scale_color_manual(values = c("darkred", "navy")) +
    theme_minimal() +
    labs(
      title = paste("Range Centroid Distance:", sp1, "vs.", sp2),
      subtitle = paste0("Geodesic Centroid Distance: ", dist_km, " km"),
      fill = "Species Range",
      color = "Range Centroid",
      x = "Longitude",
      y = "Latitude"
    )
}

# Plot Western Tanager vs Scarlet Tanager
plot_centroid_distance("Piranga ludoviciana", "Piranga olivacea", ranges)
