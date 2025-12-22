

compute_weights_batch <- function(
    var_test,
    tx_shp,          # sf with county polygons; must have column "COUNTY"
    nass_df,         # data.frame with columns: county, year, area_planted_ha, area_harvested_ha
    years,
    crop_name,
    crop_code,
    sensor        = "LANDSAT",
    season_start,    # e.g. "04-15"
    season_end,      # e.g. "09-15"
    CRS_TARGET       # e.g. "EPSG:5070"
) {
  years <- as.integer(years)
  
  for (yr in years) {
    cat("\n=== Year", yr, "===\n")
    
    # 1) CDL mosaic for this year
    cdl_fp_year <- cdl_for_processing(yr)
    if (is.na(cdl_fp_year)) {
      cat("Year", yr, ": no CDL mosaic\n")
      next
    }
    cdl_TX <- terra::rast(cdl_fp_year)
    
    # 2) Counties that actually have trials in this year
    sub_yr <- var_test[var_test$year == yr, ]
    if (nrow(sub_yr) == 0L) {
      cat("  - No variety data in", yr, "→ skip year\n")
      next
    }
    
    cnties <- unique(sub_yr$county)
    
    for (c in cnties) {
      cname <- norm_county(c)
      
      # county geometry from tx_shp
      county_sf <- tx_shp[tx_shp$COUNTY == c, ]
      if (nrow(county_sf) == 0L) {
        cat(sprintf("  - SKIP %s %d: county geometry not found in tx_shp\n", cname, yr))
        next
      }
      
      # cached weight filename
      out_fp <- .cache_path("weight", cname, yr, norm_county(crop_name), sensor)
      if (file.exists(out_fp)) {
        cat(sprintf("  ✓ SKIP %s_%s_%d_%s_%s already exists\n",
                    "weight", cname, yr, norm_county(crop_name), sensor))
        next
      }
      
      # --- lookup NASS harvested area for this county-year ---
      row_i <- nass_df[nass_df$county == cname & nass_df$year == yr, ]
      
      if (nrow(row_i) == 0L) {
        cat(sprintf("  - SKIP %s %d: no NASS row\n", cname, yr))
        next
      }
      
      # keep only row with largest harvested area
      row_i <- row_i[which.max(row_i$area_harvested_ha), , drop = FALSE]
      
      NASS_planed_area_ha  <- row_i$area_planted_ha
      NASS_harvest_area_ha <- row_i$area_harvested_ha
      
      # --- skip if no valid harvested area ---
      if (!is.finite(NASS_harvest_area_ha) ||
          is.na(NASS_harvest_area_ha)      ||
          NASS_harvest_area_ha <= 0) {
        cat(sprintf("  - SKIP %s %d: no valid NASS harvested area\n", cname, yr))
        next
      }
      
      # --- compute weights ---
      wt <- tryCatch(
        get_weights_for_county_year(
          county_sf            = county_sf,
          county_name          = cname,
          sensor               = sensor,
          year                 = yr,
          season_start         = season_start,
          season_end           = season_end,
          cdl_rast             = cdl_TX,
          CRS_TARGET           = CRS_TARGET,
          crop_code            = crop_code,
          crop                 = crop_name,
          a_ndvi               = 0.70,
          a_et                 = 0.30,
          NASS_planed_area_ha  = NASS_planed_area_ha,
          NASS_harvest_area_ha = NASS_harvest_area_ha,
          lower_thresh         = 0.7,
          upper_thresh         = 2.50,   # max ratio vs NASS area
          area_buffer_frac     = 0.15,
          pre2008_mask_strategy = "auto"
        ),
        error = function(e) {
          cat(sprintf("  - ERROR %s %d: %s\n", cname, yr, conditionMessage(e)))
          NULL
        }
      )
      
      if (is.null(wt) || !("weight" %in% names(wt))) {
        cat(sprintf("  - SKIP %s %d: weights NULL\n", cname, yr))
        next
      }
      
      wsum <- as.numeric(terra::global(wt$weight, "sum", na.rm = TRUE))
      cat(sprintf("  ✓ %s %d (weights cached, sum = %.6f, NASS_harvest_area_ha = %.1f)\n",
                  cname, yr, wsum, NASS_harvest_area_ha))
    }
  }
  
  invisible(TRUE)
}


compute_weights_batch(
  var_test   = var_test,
  tx_shp     = test_counties_sf2,
  nass_df    = nass_df,
  years      = 2018:2024,
  crop_name  = crop_name,
  crop_code  = crop_code,
  sensor     = "LANDSAT",
  season_start = season_start,
  season_end   = season_end,
  CRS_TARGET   = CRS_TARGET
)
