# 04_disagg_units_yield.R
#
# Purpose:
#   Disaggregate NASS county mean yield to units (harvest units / clusters)
#   using beta-weighted unit strengths, while enforcing crop-specific
#   absolute caps and optional relative bounds.
#
# Main API:
#   disagg_units_yield_safe()
#
# Assumptions:
#   - units_sf has columns: unit_id, unit_area_ha, unit_w_sum, geometry
#   - nass_area_ha and nass_yield_kg_ha are in ha and kg/ha
#   - beta >= 0; higher beta amplifies differences between high- and low-weight units.


# --------------------------------------------------------------------
# Internal helper: default caps by crop
# --------------------------------------------------------------------

.get_crop_caps <- function(crop) {
  crop_low <- tolower(crop %||% "")
  
  if (crop_low %in% c("cotton", "lint", "upland cotton")) {
    # kg/ha
    list(y_min = 100, y_max = 4000, min_abs = 80)
    
  } else if (crop_low %in% c("winter wheat", "wheat")) {
    list(y_min = 500, y_max = 9000, min_abs = 400)
    
  } else if (crop_low %in% c("corn", "maize")) {
    list(y_min = 700, y_max = 16000, min_abs = 600)
    
  } else {
    # generic, relatively wide caps
    list(y_min = 100, y_max = 10000, min_abs = 50)
  }
}

# small infix helper (you can also define this in a utils file)
`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------------------------------------------------
# Main function
# --------------------------------------------------------------------

#' Disaggregate NASS yield to units with crop-specific safeguards
#'
#' @param units_sf sf; must contain `unit_id`, `unit_area_ha`, `unit_w_sum`, `geometry`
#' @param nass_area_ha numeric; NASS harvested area in ha
#' @param nass_yield_kg_ha numeric; NASS mean yield (kg/ha)
#' @param county character; county name (for metadata / diagnostics)
#' @param year integer; year (for metadata / diagnostics)
#' @param beta numeric; exponent applied to unit weight sums
#' @param crop character; crop name used to select default caps (e.g. "cotton")
#' @param to_crs character; target CRS for lon/lat export (default "EPSG:4326")
#' @param min_frac numeric or NULL; optional minimum fraction of county mean
#' @param max_frac numeric or NULL; optional maximum fraction of county mean
#' @param y_min_crop numeric or NULL; optional crop-specific min yield (kg/ha)
#' @param y_max_crop numeric or NULL; optional crop-specific max yield (kg/ha)
#' @param min_abs numeric or NULL; optional legacy absolute minimum yield (kg/ha)
#'
#' @return list with:
#'   - units_with_yield: sf with per-unit yields, production, etc.
#'   - check: data.frame of diagnostics (area, mean yield, deviation from NASS)
#' @export
disagg_units_yield_safe <- function(units_sf,
                                    nass_area_ha,
                                    nass_yield_kg_ha,
                                    county,
                                    year,
                                    beta,
                                    crop   = NULL,
                                    to_crs = "EPSG:4326",
                                    # optional relative bounds (can be NULL)
                                    min_frac = NULL,
                                    max_frac = NULL,
                                    # optional crop-specific caps;
                                    # if NULL they fall back to defaults
                                    y_min_crop = NULL,
                                    y_max_crop = NULL,
                                    min_abs    = NULL) {
  stopifnot(inherits(units_sf, "sf"))
  needed_cols <- c("unit_id", "unit_area_ha", "unit_w_sum")
  if (!all(needed_cols %in% names(units_sf))) {
    stop("units_sf must have columns: ", paste(needed_cols, collapse = ", "))
  }
  
  # --------------------------------------------------------------
  # 0) Resolve crop-specific defaults
  # --------------------------------------------------------------
  caps <- .get_crop_caps(crop %||% "generic")
  
  if (is.null(y_min_crop)) y_min_crop <- caps$y_min
  if (is.null(y_max_crop)) y_max_crop <- caps$y_max
  if (is.null(min_abs))    min_abs    <- caps$min_abs
  
  # --------------------------------------------------------------
  # 1) Scale unit areas so total matches NASS harvested area
  # --------------------------------------------------------------
  area_total_raw <- sum(units_sf$unit_area_ha, na.rm = TRUE)
  if (!is.finite(area_total_raw) || area_total_raw <= 0) {
    stop("Bad area total in units_sf")
  }
  
  scaleA <- nass_area_ha / area_total_raw
  units_sf$area_ha_scaled <- units_sf$unit_area_ha * scaleA
  
  # --------------------------------------------------------------
  # 2) Beta-weighted allocation of yield
  # --------------------------------------------------------------
  units_sf$w_beta <- (units_sf$unit_w_sum)^beta
  units_sf$w_beta[!is.finite(units_sf$w_beta) | units_sf$w_beta < 0] <- 0
  
  wb_mean <- sum(units_sf$w_beta * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  yield_raw <- nass_yield_kg_ha * (units_sf$w_beta / wb_mean)
  
  # --------------------------------------------------------------
  # 3) Safeguards: crop caps + optional relative bounds
  # --------------------------------------------------------------
  # If relative bounds are NULL, they are effectively inactive.
  y_floor_rel <- if (is.null(min_frac)) -Inf else nass_yield_kg_ha * min_frac
  y_ceiling_rel <- if (is.null(max_frac)) Inf else nass_yield_kg_ha * max_frac
  
  # final floor = max(relative floor, legacy min_abs, crop minimum)
  y_floor <- max(y_min_crop, min_abs, y_floor_rel, na.rm = TRUE)
  
  # final ceiling = min(relative ceiling, crop maximum)
  y_ceiling <- min(y_ceiling_rel, y_max_crop, na.rm = TRUE)
  
  # Defensive guard: if floor > ceiling, relax floor
  if (is.finite(y_floor) && is.finite(y_ceiling) && y_floor > y_ceiling) {
    warning("For ", county, " ", year,
            ": floor > ceiling; relaxing floor to ceiling.")
    y_floor <- y_ceiling
  }
  
  yield_clamped <- pmin(pmax(yield_raw, y_floor), y_ceiling)
  
  # --------------------------------------------------------------
  # 4) Recenter so area-weighted mean ~ NASS mean (if possible)
  # --------------------------------------------------------------
  aw_mean_clamped <- sum(yield_clamped * units_sf$area_ha_scaled, na.rm = TRUE) /
    sum(units_sf$area_ha_scaled, na.rm = TRUE)
  
  adj_factor <- if (is.finite(aw_mean_clamped) && aw_mean_clamped > 0) {
    nass_yield_kg_ha / aw_mean_clamped
  } else {
    1
  }
  
  yield_adj <- yield_clamped * adj_factor
  
  # Final clipping to [y_floor, y_ceiling] again (in case recenter overshot)
  yield_final <- pmin(pmax(yield_adj, y_floor), y_ceiling)
  
  # --------------------------------------------------------------
  # 5) Production, centroids, metadata
  # --------------------------------------------------------------
  prod_final <- yield_final * units_sf$area_ha_scaled
  
  units_sf$yield_kg_ha <- yield_final
  units_sf$prod_kg     <- prod_final
  
  if (!is.null(to_crs)) {
    units_ll <- sf::st_transform(units_sf, to_crs)
  } else {
    units_ll <- units_sf
  }
  cents <- sf::st_centroid(units_ll)
  xy    <- sf::st_coordinates(cents)
  units_sf$lon <- xy[, 1]
  units_sf$lat <- xy[, 2]
  
  units_sf$county <- county
  units_sf$year   <- as.integer(year)
  
  # --------------------------------------------------------------
  # 6) QA diagnostics
  # --------------------------------------------------------------
  total_area_scaled <- sum(units_sf$area_ha_scaled, na.rm = TRUE)
  total_prod_kg     <- sum(units_sf$prod_kg,        na.rm = TRUE)
  mean_yield_check  <- total_prod_kg / total_area_scaled
  
  nass_dev_pct <- 100 * (mean_yield_check - nass_yield_kg_ha) / nass_yield_kg_ha
  
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
    adj_factor_applied         = adj_factor,
    nass_deviation_pct         = nass_dev_pct
  )
  
  list(
    units_with_yield = units_sf |>
      dplyr::select(
        .data$unit_id,
        .data$county,
        .data$year,
        .data$lon, .data$lat,
        .data$area_ha_scaled,
        .data$yield_kg_ha,
        .data$prod_kg,
        .data$unit_area_ha,
        .data$unit_w_sum,
        .data$geometry
      ),
    check = check
  )
}


