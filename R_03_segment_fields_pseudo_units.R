# 03_segment_fields_pseudo_units.R
#
# Purpose:
#   Convert a continuous weight raster (w_smooth) into a set of “harvest units”
#   by:
#     1) segmentation into contiguous patches,
#     2) tiling very large patches with a fishnet,
#     3) enforcing min/max area per unit.
#
# Main API:
#   segment_fields_diag(
#     w_smooth,
#     max_field_ha_thresh,
#     tile_target_ha,
#     min_field_ha,
#     max_field_ha_final
#   )
#
# Returns:
#   - list(patch_id, cell_ha, smallmed_poly_sf, pseudo_sf2, all_units_sf)

make_fishnet <- function(poly_sf, tile_size_m) {
  bb <- st_bbox(poly_sf)
  
  xs <- seq(bb["xmin"], bb["xmax"], by = tile_size_m)
  ys <- seq(bb["ymin"], bb["ymax"], by = tile_size_m)
  
  cells <- vector("list", (length(xs)-1)*(length(ys)-1))
  k <- 1
  for (i in seq_len(length(xs)-1)) {
    for (j in seq_len(length(ys)-1)) {
      x1 <- xs[i]; x2 <- xs[i+1]
      y1 <- ys[j]; y2 <- ys[j+1]
      cells[[k]] <- st_polygon(list(rbind(
        c(x1,y1),
        c(x2,y1),
        c(x2,y2),
        c(x1,y2),
        c(x1,y1)
      )))
      k <- k + 1
    }
  }
  
  grid <- st_sfc(cells, crs = st_crs(poly_sf))
  st_sf(tile_id = seq_along(grid), geometry = grid)
}

summarize_tile <- function(tile_sf, w_smooth, cell_ha) {
  # tile_sf: 1-row sf polygon
  v <- terra::vect(tile_sf)
  
  wv <- terra::values( terra::mask(w_smooth, v) )
  av <- terra::values( terra::mask(cell_ha,  v) )
  
  tile_area_ha <- sum(av, na.rm = TRUE)
  tile_w_sum   <- sum(wv, na.rm = TRUE)
  
  if (!is.finite(tile_area_ha) || tile_area_ha <= 0) return(NULL)
  
  list(area_ha = tile_area_ha,
       w_sum   = tile_w_sum)
}

segment_fields_diag <- function(w_smooth,
                                max_field_ha_thresh = NULL,
                                tile_target_ha      = NULL,
                                min_field_ha        = NULL,
                                max_field_ha_final  = NULL) {
  
  # 0. Base area of support from w_smooth -------------------------------
  # Assumes w_smooth is in an equal-area CRS (e.g., EPSG:5070)
  # and is NA outside the NDVI-eligible crop mask.
  
  # support mask: non-NA cells in w_smooth
  support_mask <- !is.na(w_smooth)
  support_mask <- terra::ifel(support_mask, 1, NA)
  
  cell_ha_raw <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha_sup <- terra::mask(cell_ha_raw, support_mask)
  
  support_area_ha <- terra::global(cell_ha_sup, "sum", na.rm = TRUE)[1, 1]
  cat(">> support (from mask) area (ha):", round(support_area_ha, 1), "\n")
  
  # 1. Segment into patches ---------------------------------------------
  patch_id <- terra::patches(support_mask, directions = 4)
  
  # area from patch_id BEFORE sieve
  cell_ha <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha <- terra::mask(cell_ha, patch_id)
  area_patches_before <- terra::global(cell_ha, "sum", na.rm = TRUE)[1, 1]
  cat(">> area in patches before sieve (ha):", round(area_patches_before, 1), "\n")
  
  # Remove microscopic specks (<50 cells)
  patch_id <- terra::sieve(patch_id, 50)
  
  cell_ha <- terra::cellSize(w_smooth, unit = "ha")
  cell_ha <- terra::mask(cell_ha, patch_id)
  area_patches_after <- terra::global(cell_ha, "sum", na.rm = TRUE)[1, 1]
  cat(">> area in patches after sieve (ha): ", round(area_patches_after, 1), "\n")
  
  # 2. Zonal stats per patch --------------------------------------------
  df_area <- as.data.frame(
    terra::zonal(cell_ha, patch_id, fun = "sum", na.rm = TRUE)
  ) |>
    dplyr::rename(patch = patches, area_ha = area)
  
  df_wsum <- as.data.frame(
    terra::zonal(w_smooth, patch_id, fun = "sum", na.rm = TRUE)
  ) |>
    dplyr::rename(patch = patches, w_sum = w_norm)
  
  fields <- df_area |>
    dplyr::inner_join(df_wsum, by = "patch") |>
    dplyr::filter(
      is.finite(area_ha), area_ha > 0,
      is.finite(w_sum)    # removed w_sum > 0 condition
    ) |>
    dplyr::mutate(
      w_density = w_sum / area_ha
    )
  
  cat(">> total area in fields (ha):",
      round(sum(fields$area_ha, na.rm = TRUE), 1), "\n")
  
  # 3. Split fields into small/medium vs giant --------------------------
  fields_smallmed <- fields %>% dplyr::filter(area_ha <= max_field_ha_thresh)
  fields_giant    <- fields %>% dplyr::filter(area_ha >  max_field_ha_thresh)
  
  cat(">> small/med fields area (ha):",
      round(sum(fields_smallmed$area_ha, na.rm = TRUE), 1), "\n")
  cat(">> giant fields area (ha):   ",
      round(sum(fields_giant$area_ha,    na.rm = TRUE), 1), "\n")
  
  # keep only small/med in this mask
  keep_ids <- fields_smallmed$patch
  
  smallmed_mask <- patch_id
  smallmed_mask[ !(patch_id[] %in% keep_ids) ] <- NA
  
  smallmed_poly <- terra::as.polygons(smallmed_mask, dissolve = TRUE, minArea = 0)
  smallmed_poly <- smallmed_poly[!is.na(smallmed_poly$patches), ]
  smallmed_poly_sf <- sf::st_as_sf(smallmed_poly)
  
  smallmed_poly_sf <- smallmed_poly_sf |>
    dplyr::left_join(
      fields_smallmed %>% dplyr::select(patch, area_ha, w_sum),
      by = c("patches" = "patch")
    ) |>
    dplyr::rename(
      unit_id      = patches,
      unit_area_ha = area_ha,
      unit_w_sum   = w_sum
    )
  
  smallmed_poly_sf <- sf::st_set_crs(
    smallmed_poly_sf,
    sf::st_crs(terra::crs(w_smooth))
  )
  
  cat(">> area in smallmed_poly_sf (ha):",
      round(sum(smallmed_poly_sf$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  # 4. Tilings for giant fields -----------------------------------------
  tile_side_m <- sqrt(tile_target_ha * 10000)  # ha -> m² -> side length
  pseudo_list <- list()
  
  for (p in fields_giant$patch) {
    
    this_patch_mask <- terra::ifel(patch_id == p, 1, NA)
    this_poly       <- terra::as.polygons(this_patch_mask, dissolve = TRUE, minArea = 0)
    
    id_col <- names(this_poly)[1]
    this_poly <- this_poly[!is.na(this_poly[[id_col]]), ]
    if (nrow(this_poly) == 0) {
      cat("patch", p, " -> no polygon rows after mask\n")
      next
    }
    
    this_poly_sf <- sf::st_as_sf(this_poly)
    this_poly_sf <- sf::st_set_crs(this_poly_sf, sf::st_crs(terra::crs(patch_id)))
    this_poly_sf <- sf::st_make_valid(this_poly_sf)
    
    grid_sf <- make_fishnet(this_poly_sf, tile_side_m)
    grid_sf <- sf::st_make_valid(grid_sf)
    
    tiles_clipped <- suppressWarnings(sf::st_intersection(grid_sf, this_poly_sf))
    if (nrow(tiles_clipped) == 0) {
      cat("patch", p, " -> tiles_clipped empty after intersection\n")
      next
    }
    
    for (i in seq_len(nrow(tiles_clipped))) {
      tile_i <- tiles_clipped[i, , drop = FALSE]
      
      s <- summarize_tile(tile_i, w_smooth, cell_ha)
      if (is.null(s)) next
      
      pseudo_list[[length(pseudo_list) + 1]] <- sf::st_sf(
        patch_parent = p,
        unit_id      = paste0(p, "_", i),
        unit_area_ha = s$area_ha,
        unit_w_sum   = s$w_sum,
        geometry     = sf::st_geometry(tile_i),
        crs          = sf::st_crs(this_poly_sf)
      )
    }
    
    cat("patch", p, ": tiles added, pseudo_list length =",
        length(pseudo_list), "\n")
  }
  
  if (length(pseudo_list) > 0) {
    pseudo_sf2 <- do.call(rbind, pseudo_list)
  } else {
    pseudo_sf2 <- sf::st_sf(
      patch_parent = character(0),
      unit_id      = character(0),
      unit_area_ha = numeric(0),
      unit_w_sum   = numeric(0),
      geometry     = sf::st_sfc(crs = sf::st_crs(terra::crs(patch_id)))
    )
  }
  
  pseudo_sf2 <- sf::st_set_crs(pseudo_sf2, sf::st_crs(smallmed_poly_sf))
  
  cat(">> area in pseudo_sf2 (ha):",
      round(sum(pseudo_sf2$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  smallmed_poly_sf$unit_id <- as.character(smallmed_poly_sf$unit_id)
  pseudo_sf2$unit_id       <- as.character(pseudo_sf2$unit_id)
  
  # 5. Combine and filter units ------------------------------------------
  all_units0 <- dplyr::bind_rows(
    smallmed_poly_sf[, c("unit_id","unit_area_ha","unit_w_sum","geometry")],
    pseudo_sf2[,        c("unit_id","unit_area_ha","unit_w_sum","geometry")]
  ) %>%
    dplyr::mutate(
      unit_w_density = unit_w_sum / unit_area_ha
    )
  
  cat(">> area in all_units0 (before size filter, ha):",
      round(sum(all_units0$unit_area_ha, na.rm = TRUE), 1), "\n")
  
  all_units_sf <- all_units0 %>%
    dplyr::filter(
      unit_area_ha >= min_field_ha,
      unit_area_ha <= max_field_ha_final
    ) %>%
    dplyr::arrange(dplyr::desc(unit_w_sum)) %>%
    dplyr::mutate(
      cum_area = cumsum(unit_area_ha)
    )
  
  cat(">> area in all_units_sf (after size filter, ha):",
      round(sum(all_units_sf$unit_area_ha, na.rm = TRUE), 1), "\n")
  cat(">> max cum_area (ha):",
      round(max(all_units_sf$cum_area, na.rm = TRUE), 1), "\n")
  
  list(
    patch_id         = patch_id,
    cell_ha          = cell_ha,
    smallmed_poly_sf = smallmed_poly_sf,
    pseudo_sf2       = pseudo_sf2,
    all_units_sf     = all_units_sf
  )
}
  
# Build renormalized harvest weights and push them back to polygons
# Inputs:
#   all_units_sf : sf with cols [unit_id, unit_area_ha, unit_w_sum, geometry]
#   w_smooth     : SpatRaster (same CRS/extent as units) – smoothed cotton weight
#   target_ha    : numeric – NASS harvested hectares target for this county-year
# Returns:
#   list(harvest_units_sf = sf with unit_w_sum replaced by capped, renormalized weights,
#        w_cap            = SpatRaster of capped, normalized weights, sum ~ 1)
build_harvest_units <- function(all_units_sf, w_smooth, target_ha) {
  stopifnot(inherits(all_units_sf, "sf"))
  req <- c("unit_id","unit_area_ha","unit_w_sum","geometry")
  if (!all(req %in% names(all_units_sf))) {
    stop("all_units_sf must contain: ", paste(req, collapse=", "))
  }
  if (!inherits(w_smooth, "SpatRaster")) stop("w_smooth must be a SpatRaster")
  if (!is.numeric(target_ha) || target_ha <= 0) stop("target_ha must be > 0")
  
  # Rank by cotton strength and compute cumulative area
  harvest_cands <- all_units_sf |>
    dplyr::arrange(dplyr::desc(.data$unit_w_sum)) |>
    dplyr::mutate(cum_area = cumsum(.data$unit_area_ha))
  
  A_total  <- sum(harvest_cands$unit_area_ha, na.rm = TRUE)
  A_target <- min(target_ha, A_total)
  
  # Find the index where cumulative area first reaches target
  k <- which(harvest_cands$cum_area >= A_target)[1]
  
  # Full units to include (ignore partial top-off unit for polygon simplicity)
  sel_units_full <- if (is.na(k)) {
    harvest_cands$unit_id
  } else if (k > 1) {
    harvest_cands$unit_id[seq_len(k-1)]
  } else {
    character(0)
  }
  
  # If nothing selected yet (tiny weights), pick the top unit to avoid empty set
  if (length(sel_units_full) == 0 && nrow(harvest_cands) > 0) {
    sel_units_full <- harvest_cands$unit_id[1]
  }
  
  harvest_units_sf <- all_units_sf |>
    dplyr::filter(.data$unit_id %in% sel_units_full)
  
  if (nrow(harvest_units_sf) == 0) {
    stop("No units selected for harvest subset; check inputs/target_ha.")
  }
  
  # Rasterize selected polygons -> mask -> normalize to sum ~ 1
  sel_vect <- terra::vect(harvest_units_sf)  # CRS must match w_smooth
  sel_mask <- w_smooth * NA
  sel_mask <- terra::rasterize(sel_vect, sel_mask, field = 1, touches = TRUE)
  
  w_cap <- terra::mask(w_smooth, sel_mask)
  s <- as.numeric(terra::global(w_cap, "sum", na.rm = TRUE))
  if (!is.na(s) && s > 0) w_cap <- w_cap / s
  
  # Push w_cap back to polygons: sum of pixel weights per unit
  cap_extract <- terra::extract(w_cap, sel_vect, fun = "sum", na.rm = TRUE)
  # terra::extract returns ID column + "sum"
  harvest_units_sf$unit_w_sum_cap <- cap_extract[, "w_norm"]
  
  # Renormalize polygon weights to sum exactly 1
  total_cap <- sum(harvest_units_sf$unit_w_sum_cap, na.rm = TRUE)
  if (!is.finite(total_cap) || total_cap <= 0) {
    stop("Capped weight per-unit sums to 0; inspect w_cap/selection.")
  }
  harvest_units_sf$unit_w_sum_cap <- harvest_units_sf$unit_w_sum_cap / total_cap
  # Replace unit_w_sum with the capped/renormalized version
  harvest_units_sf <- harvest_units_sf |>
    dplyr::mutate(unit_w_sum = unit_w_sum_cap) |>
    dplyr::select(-unit_w_sum_cap)
  list(
    harvest_units_sf = harvest_units_sf,
    sel_units_full = sel_units_full,
    w_cap            = w_cap
  )
}

#Keep only the top share of w_cap-weight (default 95%), ensure >= 2 units
prune_harvest_units_by_weight <- function(harvest_units_sf,
                                          keep_cumw = 0.95,
                                          min_units = 2) {
  stopifnot("unit_w_sum" %in% names(harvest_units_sf))
  df <- harvest_units_sf |>
    dplyr::arrange(dplyr::desc(.data$unit_w_sum)) |>
    dplyr::mutate(cumw = cumsum(.data$unit_w_sum))

  # indices to keep until cumw >= keep_cumw
  k <- which(df$cumw >= keep_cumw)[1]
  if (is.na(k)) k <- nrow(df)
  k <- max(k, min_units)  # ensure at least min_units are kept

  kept <- df[seq_len(k), ]
  # renormalize unit_w_sum to 1
  kept$unit_w_sum <- kept$unit_w_sum / sum(kept$unit_w_sum, na.rm = TRUE)
  kept
}

  
  
  
  
  
  
  
  
  
  

  
  