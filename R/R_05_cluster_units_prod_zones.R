# 05_cluster_units_5km.R
#
# Purpose:
#   Cluster harvest units (with yield) into ~5–10 km production zones using
#   complete-linkage clustering on centroids.
#
# Main API:
#   cluster_units_5km(units_with_yield_sf, cluster_dist_m = 10000)
#
# Returns:
#   list(cluster_sf, cluster_table)


#' Cluster harvest units into ~5–10 km production zones
#'
#' @param units_with_yield_sf sf; polygons with unit-level yield and area.
#'   Must contain columns:
#'   \itemize{
#'     \item \code{unit_id}
#'     \item \code{county}
#'     \item \code{year}
#'     \item \code{area_ha_scaled}
#'     \item \code{prod_kg}
#'     \item \code{yield_kg_ha} (optional but recommended)
#'   }
#' @param cluster_dist_m numeric; maximum within-cluster linkage distance in
#'   meters (e.g. 10000 for ~10 km clusters).
#'
#' @return A list with:
#'   \describe{
#'     \item{cluster_sf}{sf MULTIPOLYGONS with cluster_id, area, production, yield, n_units, lon/lat}
#'     \item{cluster_table}{cluster attributes as a plain data.frame (no geometry)}
#'   }
#' @export
cluster_units_5km <- function(units_with_yield_sf,
                              cluster_dist_m = 10000) {
  stopifnot(inherits(units_with_yield_sf, "sf"))
  req <- c("unit_id", "county", "year",
           "area_ha_scaled", "prod_kg")
  miss <- setdiff(req, names(units_with_yield_sf))
  if (length(miss)) {
    stop("units_with_yield_sf is missing required columns: ",
         paste(miss, collapse = ", "))
  }
  
  # Optional yield column: if absent, we reconstruct from area+prod
  if (!"yield_kg_ha" %in% names(units_with_yield_sf)) {
    units_with_yield_sf$yield_kg_ha <-
      units_with_yield_sf$prod_kg / units_with_yield_sf$area_ha_scaled
  }
  
  # ---- 1. Transform to an equal-area CRS (meters) ----
  units_m <- sf::st_transform(units_with_yield_sf, 5070)  # EPSG:5070
  
  # 2. Centroids (POINT in meters)
  cent_m <- sf::st_centroid(units_m)
  n <- nrow(cent_m)
  
  # === Short-circuit: n < 2 ===
  if (n < 2L) {
    # Single “cluster” with everything
    cent_ll <- sf::st_transform(cent_m, 4326)
    xy      <- sf::st_coordinates(cent_ll)
    
    units_m$cluster_id <- 1L
    
    cluster_area_ha <- sum(units_with_yield_sf$area_ha_scaled, na.rm = TRUE)
    cluster_prod_kg <- sum(units_with_yield_sf$prod_kg,        na.rm = TRUE)
    cluster_yield_kg_ha <- if (cluster_area_ha > 0) {
      cluster_prod_kg / cluster_area_ha
    } else {
      NA_real_
    }
    n_units <- n
    
    cluster_poly <- sf::st_union(units_m) |> sf::st_make_valid()
    cluster_sf <- sf::st_sf(
      cluster_id        = 1L,
      cluster_area_ha   = cluster_area_ha,
      cluster_prod_kg   = cluster_prod_kg,
      cluster_yield_kg_ha = cluster_yield_kg_ha,
      n_units           = n_units,
      geometry          = cluster_poly
    )
    
    cluster_cent <- sf::st_centroid(cluster_poly) |>
      sf::st_transform(4326)
    xy2 <- sf::st_coordinates(cluster_cent)
    
    cluster_tbl <- tibble::tibble(
      cluster_id        = 1L,
      lon               = xy2[, 1],
      lat               = xy2[, 2],
      cluster_area_ha   = cluster_area_ha,
      cluster_prod_kg   = cluster_prod_kg,
      cluster_yield_kg_ha = cluster_yield_kg_ha,
      n_units           = n_units
    )
    
    return(list(cluster_sf = cluster_sf,
                cluster_table = cluster_tbl))
  }
  
  # ---- 3. Centroid coordinates (matrix of X,Y in meters) ----
  cent_xy <- sf::st_coordinates(cent_m)
  
  # 4. Distance matrix and hierarchical clustering
  dist_mat <- stats::dist(cent_xy)          # class "dist"
  hc       <- stats::hclust(dist_mat, method = "complete")
  cluster_id <- stats::cutree(hc, h = cluster_dist_m)
  
  # Attach cluster ID to each unit
  units_m$cluster_id <- cluster_id
  
  # ---- 5. Summarise to cluster polygons ----
  cluster_sf <- units_m |>
    dplyr::group_by(.data$cluster_id, .data$county, .data$year) |>
    dplyr::summarise(
      cluster_area_ha   = sum(.data$area_ha_scaled, na.rm = TRUE),
      cluster_prod_kg   = sum(.data$prod_kg,        na.rm = TRUE),
      cluster_yield_kg_ha =
        cluster_prod_kg / cluster_area_ha,
      n_units           = dplyr::n(),
      geometry          = sf::st_union(.data$geometry),
      .groups           = "drop"
    )
  
  # ---- 6. Cluster centroids in lon/lat ----
  cluster_centroid_m  <- sf::st_centroid(cluster_sf)
  cluster_centroid_ll <- sf::st_transform(cluster_centroid_m, 4326)
  coords_ll <- sf::st_coordinates(cluster_centroid_ll)
  
  cluster_sf$lon <- coords_ll[, 1]
  cluster_sf$lat <- coords_ll[, 2]
  
  # final table (no geometry)
  cluster_table <- cluster_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      .data$cluster_id,
      .data$county,
      .data$year,
      .data$lon, .data$lat,
      .data$cluster_area_ha,
      .data$cluster_yield_kg_ha,
      .data$cluster_prod_kg,
      .data$n_units
    )
  
  list(
    cluster_sf    = cluster_sf,
    cluster_table = cluster_table
  )
}


