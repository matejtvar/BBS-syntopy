calculate_symmetry <- function(code1, code2, range_list, target_crs = 5070) {
  # Calculation of range symmetry
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
  
  # Sympatry Index Formula
  index <- (min(area1, area2)/sum(area1, area2))
  return(index)
}