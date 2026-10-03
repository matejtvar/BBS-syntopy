calculate_sympatry <- function(code1, code2, range_list, target_crs = 5070) {
  # Calculation of sympatry
  # --- SAFETY CHECK 1: Are the codes valid? ---
  if (is.na(code1) | is.na(code2)) return(NA)
  
  # --- SAFETY CHECK 2: Do these codes exist in the range list? ---
  if (is.null(range_list[[code1]]) | is.null(range_list[[code2]])) {
    message("Warning: Range data missing for ", code1, " or ", code2)
    return(NA)
  }
  # Transform to Equal Area CRS 5070 (Albers Equal Area)
  s1 <- st_transform(range_list[[code1]], target_crs)
  s2 <- st_transform(range_list[[code2]], target_crs)
  
  # Seasonal filtering
  s1_b <- s1[s1$season %in% c("breeding", "resident"), ]
  s2_b <- s2[s2$season %in% c("breeding", "resident"), ]
  
  if (nrow(s1_b) == 0 | nrow(s2_b) == 0) {
    return(0)
  }
  
  # Total range area calculation
  area1 <- sum(as.numeric(st_area(s1_b)))
  area2 <- sum(as.numeric(st_area(s2_b)))
  
  # Apply st_make_valid before intersecting
  inter_geom <- st_intersection(st_make_valid(s1_b), st_make_valid(s2_b))
  if (nrow(inter_geom) == 0) {
    return(0)
  } else {
    area_inter <- sum(as.numeric(st_area(inter_geom)))
    
    # Sympatry Index Formula
    index <- (area_inter/min(area1, area2)*100)
    return(index)
  }
}