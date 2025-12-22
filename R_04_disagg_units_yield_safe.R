# 04_disagg_units_yield_safe.R
#
# Purpose:
#   Assign “NASS-consistent” yields to harvest units, preserving:
#     - county total area ~ NASS area
#     - county mean yield ~ NASS yield
#   while applying floors/ceilings and beta weighting.
#
# Main API:
#   disagg_units_yield_safe(
#     units_sf,
#     nass_area_ha,
#     nass_yield_kg_ha,
#     county,
#     year,
#     beta,
#     to_crs,
#     min_frac,
#     max_frac,
#     min_abs
#   )
#
# Returns:
#   - list(units_with_yield, check)



disagg_units_yield_safe <- function(units_sf,
                                    nass_area_ha,
                                    nass_yield_kg_ha,
                                    county,
                                    year,
                                    beta = NULL,
                                    to_crs = "EPSG:4326",
                                    min_frac = min_y_frac,   # min yield = 5% of county yield
                                    max_frac = max_y_frac,    # max yield = 300% of county yield
                                    min_abs  = 80      # also enforce at least 50 kg/ha absolute
) {
  stopifnot(inherits(units_sf, "sf"))
  needed_cols <- c("unit_id", "unit_area_ha", "unit_w_sum")
  if (!all(needed_cols %in% names(units_sf))) {
    stop("units_sf must have columns: unit_id, unit_area_ha, unit_w_sum")
  }
  
  #-------------------------------------------------
  # 1. Scale area so total == NASS harvested area
  #-------------------------------------------------
  area_total_raw <- sum(units_sf$unit_area_ha, na.rm = TRUE)
  if (!is.finite(area_total_raw) || area_total_raw <= 0) {
    stop("Bad area total in units_sf")
  }
  
  scaleA <- nass_area_ha / area_total_raw
  units_sf$area_ha_scaled <- units_sf$unit_area_ha * scaleA
  # sum(area_ha_scaled) ~= nass_area_ha
  
  #-------------------------------------------------
  # 2. Weight strength (analogous to w^beta)
  #-------------------------------------------------
  units_sf$w_beta <- (units_sf$unit_w_sum)^beta
  units_sf$w_beta[!is.finite(units_sf$w_beta) | units_sf$w_beta < 0] <- 0
  
  # weighted mean of w_beta across units
  wb_mean <- sum(units_sf$w_beta * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  # raw yield allocation by weight ratio
  yield_raw <- nass_yield_kg_ha * (units_sf$w_beta / wb_mean)
  
  #-------------------------------------------------
  # 3. Safeguards: floor + ceiling
  #-------------------------------------------------
  # floor and ceiling in kg/ha based on county yield
  y_floor_frac <- nass_yield_kg_ha * min_frac   # e.g. 5% of county mean
  y_floor_abs  <- min_abs                       # e.g. at least 50 kg/ha
  y_floor      <- max(y_floor_frac, y_floor_abs)
  
  y_ceiling    <- nass_yield_kg_ha * max_frac   # e.g. 3x county yield
  
  yield_clamped <- pmin(pmax(yield_raw, y_floor), y_ceiling)
  
  #-------------------------------------------------
  # 4. Recenter so county mean ~ NASS yield
  #-------------------------------------------------
  # After clamping, mean may have drifted.
  # Compute area-weighted mean of clamped yields:
  aw_mean_clamped <- sum(yield_clamped * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  # Scale all yields so weighted mean ~ nass_yield_kg_ha again
  # (This preserves relative differences but restores county average.)
  adj_factor <- if (is.finite(aw_mean_clamped) && aw_mean_clamped > 0) {
    nass_yield_kg_ha / aw_mean_clamped
  } else {
    1
  }
  
  yield_adj <- yield_clamped * adj_factor
  
  # OPTIONAL final tiny floor just to kill any resurrected 1 kg/ha edge cases
  # (This might introduce a *very* small bias away from the exact county mean,
  #  but protects against pathological 0 clusters.)
  yield_final <- pmax(yield_adj, y_floor)
  
  #-------------------------------------------------
  # 5. Production per unit
  #-------------------------------------------------
  prod_final <- yield_final * units_sf$area_ha_scaled
  
  units_sf$yield_kg_ha <- yield_final
  units_sf$prod_kg     <- prod_final
  
  #-------------------------------------------------
  # 6. Attach lon/lat centroid for export
  #-------------------------------------------------
  if (!is.null(to_crs)) {
    units_ll <- sf::st_transform(units_sf, to_crs)
  } else {
    units_ll <- units_sf
  }
  cents <- sf::st_centroid(units_ll)
  xy    <- sf::st_coordinates(cents)
  units_sf$lon <- xy[,1]
  units_sf$lat <- xy[,2]
  
  # metadata
  units_sf$county <- county
  units_sf$year   <- as.integer(year)
  
  #-------------------------------------------------
  # 7. QA check
  #-------------------------------------------------
  total_area_scaled <- sum(units_sf$area_ha_scaled, na.rm = TRUE)
  total_prod_kg     <- sum(units_sf$prod_kg,        na.rm = TRUE)
  mean_yield_check  <- total_prod_kg / total_area_scaled
  
  check <- data.frame(
    county = county,
    year   = year,
    nass_area_ha               = nass_area_ha,
    units_total_area_scaled_ha = total_area_scaled,
    nass_yield_kg_ha           = nass_yield_kg_ha,
    area_weighted_mean_yield   = mean_yield_check,
    expected_prod_kg           = nass_area_ha * nass_yield_kg_ha,
    units_total_prod_kg        = total_prod_kg,
    y_floor_used_kg_ha         = y_floor,
    y_ceiling_used_kg_ha       = y_ceiling,
    adj_factor_applied         = adj_factor
  )
  
  list(
    units_with_yield = units_sf %>%
      dplyr::select(
        unit_id,
        county,
        year,
        lon, lat,
        area_ha_scaled,
        yield_kg_ha,
        prod_kg,
        unit_area_ha,
        unit_w_sum,
        geometry
      ),
    check = check
  )
}
