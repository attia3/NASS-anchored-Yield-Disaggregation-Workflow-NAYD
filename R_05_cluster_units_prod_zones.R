# 05_cluster_units_5km.R
#
# Purpose:
#   Cluster harvest units (with yield) into ~5–10 km production zones using
#   complete-linkage clustering on centroids.
#
# Main API:
#   cluster_units_5km(units_with_yield_sf, cluster_dist_m)
#   merge_excess_floor_clusters(clusters_sf, min_abs, keep_floor, ...)
#
# Returns:
#   - list(cluster_sf, cluster_table)


cluster_units_5km <- function(units_with_yield_sf,
                              cluster_dist_m = cluster_dist_m) {
  
  # 1. Make sure we're in EPSG:5070 (meters)
  units_m <- st_transform(units_with_yield_sf, 5070)
  
  # 2. Centroids (POINT in meters)
  cent_m <- st_centroid(units_m)
  
  n <- nrow(cent_m)
  
  # === Short-circuit: n < 2 ===
  if (n < 2) {
    # build a single cluster 1: aggregate equals the unit itself
    # centroid in WGS84 for lon/lat
    cent_ll <- sf::st_transform(cent_m, 4326)
    xy      <- sf::st_coordinates(cent_ll)
    out_sf  <- units_m
    out_sf$cluster_id <- 1L
    
    # aggregate to cluster-level (one row)
    cluster_area_ha   <- sum(units_with_yield_sf$area_ha_scaled, na.rm = TRUE)
    cluster_prod_kg   <- sum(units_with_yield_sf$prod_kg,        na.rm = TRUE)
    cluster_yield_kg_ha <- if (cluster_area_ha > 0) cluster_prod_kg / cluster_area_ha else NA_real_
    n_units <- n
    
    cluster_poly <- sf::st_union(units_m) |> sf::st_make_valid()
    cluster_sf <- sf::st_sf(
      cluster_id = 1L,
      cluster_area_ha = cluster_area_ha,
      cluster_prod_kg = cluster_prod_kg,
      cluster_yield_kg_ha = cluster_yield_kg_ha,
      n_units = n_units,
      geometry = cluster_poly
    )
    
    # table with lon/lat of centroid
    cluster_cent <- sf::st_centroid(cluster_poly) |> sf::st_transform(4326)
    xy2 <- sf::st_coordinates(cluster_cent)
    cluster_tbl <- tibble::tibble(
      cluster_id = 1L,
      lon = xy2[,1], lat = xy2[,2],
      cluster_area_ha = cluster_area_ha,
      cluster_prod_kg = cluster_prod_kg,
      cluster_yield_kg_ha = cluster_yield_kg_ha,
      n_units = n_units
    )
    
    return(list(cluster_sf = cluster_sf, cluster_table = cluster_tbl))
  }
  
  # 3. Extract centroid coordinates (matrix of X,Y in meters)
  cent_xy <- st_coordinates(cent_m)
  
  # 4. Distance matrix (in meters)
  dist_mat <- dist(cent_xy)  # class "dist"
  
  # hierarchical clustering using complete linkage
  hc <- hclust(dist_mat, method = "complete")
  
  # cut the tree at 20 km (20,000 meters) # changed to 10 km max dis for more realistic clusters
  cluster_id <- cutree(hc, h = cluster_dist_m)
  
  # attach cluster ID to each unit
  units_m$cluster_id <- cluster_id
  
  
  # group and summarise
  cluster_sf <- units_m %>%
    group_by(cluster_id, county, year) %>%
    summarise(
      cluster_area_ha      = sum(area_ha_scaled, na.rm = TRUE),
      cluster_prod_kg      = sum(prod_kg,        na.rm = TRUE),
      # area-weighted mean yield for the block:
      cluster_yield_kg_ha  = cluster_prod_kg / cluster_area_ha,
      n_units              = dplyr::n(),
      # dissolve all polygons in this cluster so we have one MULTIPOLYGON
      geometry             = st_union(geometry),
      .groups = "drop"
    )
  
  # get centroid of each aggregated block
  cluster_centroid_m  <- st_centroid(cluster_sf)
  
  # convert centroid to lon/lat for export/ML
  cluster_centroid_ll <- st_transform(cluster_centroid_m, 4326)
  coords_ll <- st_coordinates(cluster_centroid_ll)
  
  cluster_sf$lon <- coords_ll[,1]
  cluster_sf$lat <- coords_ll[,2]
  
  # final table you can export / model on (drop geometry if you want CSV)
  cluster_table <- cluster_sf %>%
    st_drop_geometry() %>%
    select(
      cluster_id,
      county,
      year,
      lon, lat,
      cluster_area_ha,
      cluster_yield_kg_ha,
      cluster_prod_kg,
      n_units
    )
  
  list( 
    cluster_sf = cluster_sf,
    cluster_table = cluster_table
    )
}

# Merge excess "floor" clusters (those at min_abs) into nearby clusters
# Keeps 'keep_floor' floor clusters (largest area), merges the rest into the
# nearest cluster (prefer non-floor; else nearest kept floor). Totals preserved.
merge_excess_floor_clusters <- function(
    clusters_sf,
    min_abs,
    keep_floor = 1,
    prefer_nonfloor = TRUE,
    county,
    year,
    add_lonlat = TRUE
) {
  stopifnot(inherits(clusters_sf, "sf"))
  req <- c("cluster_id","cluster_area_ha","cluster_prod_kg","cluster_yield_kg_ha","n_units")
  miss <- setdiff(req, names(clusters_sf))
  if (length(miss)) stop("cluster_sf missing: ", paste(miss, collapse = ", "))
  
  # --- neutralize AGR to avoid sf’s column-index assertions
  suppressWarnings(sf::st_agr(clusters_sf) <- sf::NA_agr_)
  
  eps <- 1e-6
  is_floor <- clusters_sf$cluster_yield_kg_ha <= (min_abs + eps)
  n_floor  <- sum(is_floor, na.rm = TRUE)
  
  # If nothing to merge, just enrich and return
  if (n_floor <= keep_floor) {
    out_sf <- clusters_sf
  } else {
    floor_sf    <- clusters_sf[is_floor, ]
    nonfloor_sf <- clusters_sf[!is_floor, ]
    
    # keep largest 'keep_floor' floors by area
    floor_sf <- dplyr::arrange(floor_sf, dplyr::desc(.data$cluster_area_ha))
    floor_keep  <- floor_sf[seq_len(keep_floor), ]
    floor_extra <- if (nrow(floor_sf) > keep_floor) floor_sf[-seq_len(keep_floor), ] else floor_sf[0,]
    
    # pool to merge into
    target_pool <- if (prefer_nonfloor && nrow(nonfloor_sf) > 0) nonfloor_sf else floor_keep
    
    # mutable copy
    upd_pool <- target_pool
    
    # centroids for distance
    target_cents <- sf::st_centroid(upd_pool)
    
    for (i in seq_len(nrow(floor_extra))) {
      f_row  <- floor_extra[i, , drop = FALSE]
      f_cent <- sf::st_centroid(f_row)
      
      # nearest target
      dmat <- sf::st_distance(f_cent, target_cents)
      j <- which.min(as.numeric(dmat))
      
      # merge attributes
      merged_area <- upd_pool$cluster_area_ha[j] + f_row$cluster_area_ha
      merged_prod <- upd_pool$cluster_prod_kg[j] + f_row$cluster_prod_kg
      merged_yld  <- if (merged_area > 0) merged_prod / merged_area else upd_pool$cluster_yield_kg_ha[j]
      merged_geom <- suppressWarnings(sf::st_make_valid(sf::st_union(sf::st_geometry(upd_pool[j,]), sf::st_geometry(f_row))))
      
      upd_pool$cluster_area_ha[j]   <- merged_area
      upd_pool$cluster_prod_kg[j]   <- merged_prod
      upd_pool$cluster_yield_kg_ha[j]<- merged_yld
      sf::st_geometry(upd_pool)[j]  <- merged_geom
    }
    
    # rebuild out set
    if (prefer_nonfloor && nrow(nonfloor_sf) > 0) {
      # replace target pool with updated
      out_sf <- dplyr::bind_rows(
        upd_pool,
        floor_keep  # unchanged keepers if target_pool was nonfloor
      )
      # remove original target_pool rows to avoid duplicates
      # simplest: distinct by geometry+id (light dedupe)
      out_sf <- dplyr::distinct(out_sf, cluster_id, .keep_all = TRUE)
    } else {
      # all floors: we merged extras into floor_keep
      out_sf <- upd_pool
    }
  }
  
  # clean + order
  out_sf <- sf::st_make_valid(out_sf)
  out_sf <- dplyr::arrange(out_sf, dplyr::desc(.data$cluster_yield_kg_ha))
  
  # attach county/year
  out_sf$county <- county
  out_sf$year   <- as.integer(year)
  
  # centroids (WGS84) if requested
  if (isTRUE(add_lonlat)) {
    cent_ll <- out_sf |>
      sf::st_centroid() |>
      sf::st_transform(4326)
    xy <- sf::st_coordinates(cent_ll)
    out_sf$lon <- xy[,1]
    out_sf$lat <- xy[,2]
  }
  
  # table to save
  out_tbl <- out_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(
      county, year, lon, lat,
      cluster_area_ha, cluster_prod_kg, cluster_yield_kg_ha
    )
  
  list(cluster_sf = out_sf, cluster_table = out_tbl)
}


