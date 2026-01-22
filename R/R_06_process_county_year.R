# 06_process_county_year.R
#
# High-level wrapper for a single county-year:
#   1) take weight raster
#   2) segment into harvest units
#   3) select harvested subset & renormalize weights
#   4) disaggregate NASS yield to units
#   5) cluster into production zones
#   6) merge excess floor clusters
#   7) return clusters + QA bundle
#
# Assumes:
#   - w_norm_orig is an NDVI/ET-based weight raster (sum ~ 1 over eligible crop)
#   - All helpers are in the package namespace.

#' Process a single county-year (NAYD pipeline)
#'
#' @param w_norm_orig SpatRaster; county/year weight raster (any CRS, will be
#'   reprojected to EPSG:5070 internally).
#' @param county_geom sf; county polygon geometry (any CRS).
#' @param max_field_ha_thresh numeric; max area (ha) for a “natural” field
#'   before tiling into pseudo-units.
#' @param tile_target_ha numeric; target area (ha) for fishnet tiles in large
#'   fields (e.g. ~200 ha).
#' @param min_field_ha numeric; minimum allowed harvest-unit area (ha).
#' @param max_field_ha_final numeric; maximum allowed harvest-unit area (ha)
#'   after tiling.
#' @param cluster_dist_m numeric; cut height (meters) for hierarchical
#'   clustering of unit centroids (e.g. 10000 for ~10 km).
#' @param nass_area_ha numeric; NASS harvested area (ha) for this county-year.
#' @param nass_yield_kg_ha numeric; NASS county mean yield (kg/ha).
#' @param county_name character; county name.
#' @param crop_name character; crop name (used by disagg caps).
#' @param year_val integer; year.
#' @param beta numeric; exponent applied to unit-level weights.
#' @param min_frac,max_frac optional numeric; relative min/max yield fractions
#'   vs NASS mean. If \code{NULL}, the disagg function will rely on crop caps.
#' @param y_min_crop,y_max_crop,min_abs optional numeric; crop-specific caps.
#'   If \code{NULL}, defaults are taken from \code{.get_crop_caps(crop_name)}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{prep}{list with reprojected weight raster and capped weights}
#'     \item{seg}{segmentation objects (units before harvest selection)}
#'     \item{harv}{selected harvest units}
#'     \item{dis}{units with disaggregated yield + QA check}
#'     \item{clusters}{final merged clusters (sf + table)}
#'     \item{info}{county/year metadata}
#'   }
#' @export
process_county_year <- function(
    w_norm_orig,
    county_geom,
    max_field_ha_thresh,
    tile_target_ha,
    min_field_ha,
    max_field_ha_final,
    cluster_dist_m,
    nass_area_ha,
    nass_yield_kg_ha,
    county_name,
    crop_name,
    year_val,
    beta,
    min_frac  = NULL,
    max_frac  = NULL,
    y_min_crop = NULL,
    y_max_crop = NULL,
    min_abs    = NULL
) {
  # ---- basic checks ----
  stopifnot(inherits(w_norm_orig, "SpatRaster"))
  stopifnot(inherits(county_geom, "sf"))
  if (!is.numeric(nass_area_ha) || nass_area_ha <= 0)
    stop("nass_area_ha must be > 0")
  if (!is.numeric(nass_yield_kg_ha) || nass_yield_kg_ha <= 0)
    stop("nass_yield_kg_ha must be > 0")
  if (!is.numeric(beta))
    stop("beta must be numeric")
  if (!is.numeric(cluster_dist_m) || cluster_dist_m <= 0)
    stop("cluster_dist_m must be > 0")
  
  year_val <- as.integer(year_val)
  
  # ---- 1) Reproject raster + county to EPSG:5070 and crop/mask ----
  w_norm_5070 <- terra::project(w_norm_orig, "EPSG:5070", method = "bilinear")
  county_5070 <- sf::st_transform(county_geom, crs = terra::crs(w_norm_5070))
  
  w_norm_crop <- terra::crop(w_norm_5070, terra::vect(county_5070))
  w_norm_crop <- terra::mask(w_norm_crop, terra::vect(county_5070))
  
  if (all(is.na(terra::values(w_norm_crop, mat = FALSE)))) {
    stop("Weight raster has no non-NA cells after cropping/masking to county.")
  }
  
  # ---- 2) Segment fields + tile giants into harvest units ----
  seg <- segment_fields_diag(
    w_smooth           = w_norm_crop,
    max_field_ha_thresh = max_field_ha_thresh,
    tile_target_ha      = tile_target_ha,
    min_field_ha        = min_field_ha,
    max_field_ha_final  = max_field_ha_final
  )
  
  all_units_sf <- seg$all_units_sf
  if (nrow(all_units_sf) == 0L) {
    stop("No units produced by segmentation for county ", county_name,
         " year ", year_val, ". Check weights / thresholds.")
  }
  
  # ---- 3) Select harvested units to match NASS area (weight-based) ----
  harv <- build_harvest_units(
    all_units_sf = all_units_sf,
    w_smooth     = w_norm_crop,
    target_ha    = nass_area_ha
  )
  
  # Optionally prune by cumulative weight (e.g. keep 95% of w_cap)
  harvest_units_sf <- prune_harvest_units_by_weight(
    harvest_units_sf = harv$harvest_units_sf,
    keep_cumw        = 0.95,
    min_units        = 2
  )
  
  # dynamic floor for *cluster* merging (not for disagg caps)
  min_abs_dynamic <- max(80, 0.2 * nass_yield_kg_ha)
  
  # ---- 4) Disaggregate NASS yield to harvested units ----
  dis <- disagg_units_yield_safe(
    units_sf         = harvest_units_sf,
    nass_area_ha     = nass_area_ha,
    nass_yield_kg_ha = nass_yield_kg_ha,
    county           = county_name,
    year             = year_val,
    beta             = beta,
    crop             = crop_name,
    to_crs           = "EPSG:4326",
    min_frac         = min_frac,
    max_frac         = max_frac,
    y_min_crop       = y_min_crop,
    y_max_crop       = y_max_crop,
    min_abs          = min_abs
  )
  
  units_with_yield <- dis$units_with_yield
  qa_check         <- dis$check
  
  # ---- 5) Cluster units into production zones ----
  cl <- cluster_units_5km(
    units_with_yield_sf = units_with_yield,
    cluster_dist_m      = cluster_dist_m
  )
  
  eps <- 1e-6
  floors_before <- sum(
    round(cl$cluster_sf$cluster_yield_kg_ha, 6) <=
      round(min_abs_dynamic, 6) + eps,
    na.rm = TRUE
  )
  
  # ---- 6) Merge excess floor clusters ----
  merged <- merge_excess_floor_clusters(
    clusters_sf     = cl$cluster_sf,
    min_abs         = min_abs_dynamic,
    keep_floor      = 1,           # keep only one floor cluster
    prefer_nonfloor = TRUE,
    county          = county_name,
    year            = year_val,
    add_lonlat      = TRUE
  )
  
  clusters_sf_final <- merged$cluster_sf
  
  clusters_table_out <- merged$cluster_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      .data$county, .data$year,
      .data$lon, .data$lat,
      .data$cluster_area_ha,
      .data$cluster_prod_kg,
      .data$cluster_yield_kg_ha,
      .data$n_units
    ) |>
    dplyr::mutate(
      nass_area_ha      = nass_area_ha,
      nass_yield_kg_ha  = nass_yield_kg_ha
    )
  
  floors_after <- sum(
    round(clusters_sf_final$cluster_yield_kg_ha, 6) <=
      round(min_abs_dynamic, 6) + eps,
    na.rm = TRUE
  )
  
  # ---- 7) QA logging (area + area-weighted mean yield) ----
  total_cluster_area <- sum(clusters_table_out$cluster_area_ha, na.rm = TRUE)
  total_cluster_prod <- sum(clusters_table_out$cluster_prod_kg, na.rm = TRUE)
  mean_yield_clusters <- if (total_cluster_area > 0) {
    total_cluster_prod / total_cluster_area
  } else {
    NA_real_
  }
  
  cat(sprintf(
    "[%s %d] NASS: area=%.0f ha | yield=%.1f kg/ha  | clusters=%d (floors: %d→%d) || clusters: area=%.0f ha, AW-mean=%.1f kg/ha\n",
    county_name, year_val,
    nass_area_ha,
    nass_yield_kg_ha,
    nrow(clusters_sf_final),
    floors_before, floors_after,
    total_cluster_area,
    mean_yield_clusters
  ))
  
  # ---- 8) Return bundle ----
  list(
    prep = list(
      w_smooth = w_norm_crop,
      w_cap    = harv$w_cap
    ),
    seg = list(
      all_units_sf     = seg$all_units_sf,
      smallmed_poly_sf = seg$smallmed_poly_sf,
      pseudo_sf2       = seg$pseudo_sf2
    ),
    harv = list(
      harvest_units_sf = harvest_units_sf,
      sel_units_full   = harv$sel_units_full
    ),
    dis = list(
      units_with_yield = units_with_yield,
      check            = qa_check
    ),
    clusters = list(
      clusters_sf    = clusters_sf_final,
      clusters_table = clusters_table_out
    ),
    info = list(
      county = county_name,
      crop   = crop_name,
      year   = year_val
    )
  )
}

