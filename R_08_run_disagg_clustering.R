# 08_run_disagg_clustering.R
#
# Main driver script:
#  1) load NASS tables, county shapefiles, and variety-testing data
#  2) loop over years × counties × beta values
#  3) run process_county_year() to disaggregate NASS yields
#  4) validate segmented yields against variety trials and plot 1:1
#
# NOTE:
#  - This script assumes that weight rasters for each county-year
#    have already been created by get_weights_for_county_year()
#    and cached via .cache_path("weight", ...).

# ------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------
library(sf)
library(terra)
library(dplyr)
library(data.table)
library(stringr)
library(purrr)
library(ggplot2)
library(ggspatial)
library(patchwork)
library(rlang)
library(viridisLite)
library(viridis)
library(cowplot)
library(readxl)
library(tidyr)
library(tibble)

# ------------------------------------------------------------------
# 2. Source project functions (use relative paths in the repo)
# ------------------------------------------------------------------
source("R/03_segment_fields_pseudo_units.R")
source("R/04_disagg_units_yield_safe.R")
source("R/05_cluster_units_5km.R")
source("R/06_process_county_year.R")
source("R/02_get_weights_county_year.R")
source("R/01_ndvi_et_season.R")
# optional: source("R/07_compute_weights_batch.R")

# ------------------------------------------------------------------
# 3. Paths & basic settings
#    (for the repo, keep these relative to the project root)
# ------------------------------------------------------------------

# Root dirs (adapt if you use a different structure)
data_dir  <- "data"
cache_dir <- "cache"     # where .cache_path() writes
out_dir   <- file.path(cache_dir, "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Weight directory (used only for sanity checks if needed)
weight_dir <- file.path(cache_dir, "weight")

# NASS tables (choose the crop; wheat example here)
# -------------------------------------------------
# Cotton (example)
# nass_csv <- file.path(data_dir, "TX_CO_PLAN_HARV_Y.csv")

# Winter wheat (example)
nass_csv <- file.path(data_dir, "TX_WW_PLAN_HARV_Y.csv")

nass_df <- read.csv(nass_csv, stringsAsFactors = FALSE)

# Shapefiles
tx_shp_path  <- file.path(data_dir, "TX_counties.shp")
thp_shp_path <- file.path(data_dir, "TXhighPlains_counties.shp")

TX_shp  <- st_read(tx_shp_path, quiet = TRUE)  |> st_make_valid() |> st_transform(4326)
THP_shp <- st_read(thp_shp_path, quiet = TRUE) |> st_make_valid() |> st_transform(4326)

THP_counties <- THP_shp$COUNTY

# Variety-testing data (wheat example; cotton version commented)
# --------------------------------------------------------------
# Cotton:
# var_test <- read.csv(file.path(data_dir, "RACE_validation_Cotton.csv"),
#                      stringsAsFactors = FALSE)

# Wheat:
var_test <- read_excel(file.path(data_dir, "wheat_var_test.xlsx"))  # handle encoding in Excel

# clean county names (remove non-breaking spaces etc.)
var_test$county <- trimws(gsub("\u00A0", " ", enc2utf8(var_test$county)))

# THP core counties (optional, used for filtering)
thp_core_path <- file.path(data_dir, "THP_core_counties.gpkg")
core_counties <- sf::st_read(thp_core_path, quiet = TRUE) |>
  dplyr::mutate(COUNTY = norm_county(COUNTY))

# ------------------------------------------------------------------
# 4. Select counties & years to process
# ------------------------------------------------------------------

# Counties present in variety trials (excluding 2024 if desired)
# NOTE: adapt filter if you want to include/exclude years
years <- sort(unique(var_test$year))

test_counties <- unique(var_test$county)                   # all years
# test_counties <- unique(var_test$county[var_test$year != 2024])  # example filter

test_counties_sf <- TX_shp[TX_shp$COUNTY %in% test_counties, ]

# Subset NASS table to relevant counties/years
nass_table <- nass_df[
  nass_df$county %in% norm_county(test_counties) &
    nass_df$year   %in% years,
]

# ------------------------------------------------------------------
# 5. Disaggregation + clustering loop
# ------------------------------------------------------------------

results_list <- list()
idx <- 1

# beta values for testing
betas <- c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8)

# You may have multiple crops; define here for .cache_path():
crop_name <- "winter_wheat"  # or "cotton"
CRS_TARGET <- "EPSG:5070"    # ensure this matches the rest of your project

for (yr in years) {
  sub_yr <- var_test[var_test$year == yr, ]
  cnties <- unique(sub_yr$county)
  
  for (c in cnties) {
    county_i <- c
    year_i   <- yr
    
    # --- NASS area and yield for this county-year ---
    nass_rows <- nass_table[
      nass_table$county == norm_county(county_i) &
        nass_table$year   == year_i,
    ]
    
    if (nrow(nass_rows) == 0L) {
      message("[SKIP] No NASS row for ", county_i, " ", year_i)
      next
    }
    
    # If multiple entries exist, take the one with max harvested area
    nass_rows <- nass_rows[which.max(nass_rows$area_harvested_ha), , drop = FALSE]
    
    nass_area_ha_i      <- nass_rows$area_harvested_ha
    nass_yield_kg_ha_i  <- nass_rows$yield_kg_ha
    
    # county geometry (EPSG:4326 here; process_county_year will handle CRS)
    county_geom_i <- test_counties_sf[test_counties_sf$COUNTY == county_i, ]
    if (nrow(county_geom_i) == 0L) {
      message("[SKIP] No county geometry for ", county_i)
      next
    }
    
    # weight raster path
    w_path_i <- .cache_path(
      "weight",
      norm_county(county_i),
      year_i,
      crop_name = crop_name,
      sensor    = "LANDSAT"
    )
    
    if (!file.exists(w_path_i)) {
      message("[SKIP] No weight raster for ", county_i, " ", year_i,
              " at ", w_path_i)
      next
    }
    
    w_norm_orig_i <- terra::rast(w_path_i)
    
    for (b in betas) {
      message("[INFO] ", county_i, " ", year_i, " - beta = ", b)
      
      out <- process_county_year(
        w_norm_orig        = w_norm_orig_i,
        county_geom        = county_geom_i,
        max_field_ha_thresh = 300,
        tile_target_ha     = 200,
        min_field_ha       = 5,
        max_field_ha_final = 500,
        cluster_dist_m     = 10000,
        nass_area_ha       = nass_area_ha_i,
        nass_yield_kg_ha   = nass_yield_kg_ha_i,
        county_name        = county_i,
        year_val           = year_i,
        beta               = b,
        min_y_frac         = 0.01,
        max_y_frac         = 15,
        min_abs_dynamic    = 50
      )
      
      cl_tab <- out$clusters$clusters_table %>%
        dplyr::mutate(
          beta   = b,
          county = county_i,
          year   = year_i
        )
      
      results_list[[idx]] <- cl_tab
      idx <- idx + 1
    }
  }
}

# stack all county-year cluster rows into one big DF
seg_df_ww <- dplyr::bind_rows(results_list)

# ------------------------------------------------------------------
# 6. Save segmented yields
# ------------------------------------------------------------------
seg_out_csv <- file.path(out_dir, "segmented_wheat.csv")
write.csv(seg_df_ww, seg_out_csv, row.names = FALSE)

# (optional) re-read to confirm
# seg_df_ww <- read.csv(seg_out_csv, stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# 7. Validation: variety trials vs segmented yields
# ------------------------------------------------------------------

# seg_df_ww must contain: county, year, lon, lat, cluster_yield_kg_ha, beta
# Make sure those columns exist in seg_df_ww before proceeding.

# Convert segmented clusters (for all betas) later, but for matching we
# work beta-by-beta.

# Convert variety trials to a single yield column if needed
if ("lint_yield_lb_ac" %in% names(var_test)) {
  var_test <- var_test %>%
    dplyr::mutate(var_yield_kg_ha = lint_yield_lb_ac * 1.12085)
}

if ("Yield_bu_ac" %in% names(var_test)) {
  var_test <- var_test %>%
    dplyr::mutate(var_yield_kg_ha = Yield_bu_ac * 67.25)
}

# add point id (to join later)
var_test <- var_test %>%
  mutate(point_id = dplyr::row_number())

# Variety trials as sf (WGS84)
var_sf <- st_as_sf(
  var_test,
  coords = c("longitude", "latitude"),
  crs    = 4326
)

# Reproject to metric CRS (e.g. EPSG:5070)
var_sf_p <- st_transform(var_sf, 5070)

# Radii to test (km → m)
max_radius_m <- c(20000, 30000, 40000)

# Matching function for a single point and radius
match_one_point <- function(i, max_radius_m, seg_sf_p) {
  obs_y <- var_sf_p$var_yield_kg_ha[i]
  yr_i  <- var_sf_p$year[i]
  cty_i <- var_sf_p$county[i]
  
  # restrict candidate segments: same year + same county
  seg_idx_pool <- which(seg_sf_p$year == yr_i &
                          seg_sf_p$county == cty_i)
  
  if (length(seg_idx_pool) == 0L) {
    return(tibble(
      point_id   = i,
      radius_m   = max_radius_m,
      match_idx  = NA_integer_,
      dist_m     = NA_real_,
      y_seg_match = NA_real_
    ))
  }
  
  # distances only to candidate segments
  dists <- st_distance(var_sf_p[i, ], seg_sf_p[seg_idx_pool, , drop = FALSE])
  dists <- as.numeric(dists)  # meters
  
  res_list <- vector("list", length(max_radius_m))
  
  for (k in seq_along(max_radius_m)) {
    r <- max_radius_m[k]
    
    cand_local <- which(dists <= r)
    
    if (length(cand_local) == 0L) {
      res_list[[k]] <- tibble(
        point_id   = i,
        radius_m   = r,
        match_idx  = NA_integer_,
        dist_m     = NA_real_,
        y_seg_match = NA_real_
      )
    } else {
      seg_idx_cand <- seg_idx_pool[cand_local]
      y_seg        <- seg_sf_p$cluster_yield_kg_ha[seg_idx_cand]
      
      # among candidates, pick segment with yield closest to observed
      j_loc  <- which.min(abs(y_seg - obs_y))
      seg_idx <- seg_idx_cand[j_loc]
      
      res_list[[k]] <- tibble(
        point_id   = i,
        radius_m   = r,
        match_idx  = seg_idx,
        dist_m     = dists[cand_local[j_loc]],
        y_seg_match = seg_sf_p$cluster_yield_kg_ha[seg_idx]
      )
    }
  }
  
  bind_rows(res_list)
}

# ------------------------------------------------------------------
# 8. Loop over betas and radii, compute matches
# ------------------------------------------------------------------

res_list <- list()
idx <- 1

for (b in betas) {
  seg_df_beta <- seg_df_ww %>% filter(beta == b)
  
  if (nrow(seg_df_beta) == 0L) {
    warning("No segments for beta = ", b)
    next
  }
  
  seg_sf <- st_as_sf(seg_df_beta,
                     coords = c("lon", "lat"),
                     crs    = 4326)
  seg_sf_p <- st_transform(seg_sf, 5070)
  
  for (r in max_radius_m) {
    matches_tbl <- purrr::map_dfr(
      seq_len(nrow(var_sf_p)),
      ~ match_one_point(.x, max_radius_m = r, seg_sf_p = seg_sf_p)
    )
    
    res_list[[idx]] <- matches_tbl %>%
      mutate(
        beta      = b,
        radius_km = r / 1000
      )
    idx <- idx + 1
  }
}

all_matches <- bind_rows(res_list)

df_all <- var_test %>%
  left_join(all_matches, by = "point_id")

# Keep only valid yield pairs
df_plot_all <- df_all %>%
  filter(!is.na(var_yield_kg_ha),
         !is.na(y_seg_match))

# ------------------------------------------------------------------
# 9. Goodness-of-fit stats per (beta, radius_km)
# ------------------------------------------------------------------

stats_df <- df_plot_all %>%
  group_by(beta, radius_km) %>%
  summarise(
    n        = n(),
    r2       = cor(var_yield_kg_ha, y_seg_match)^2,
    rmse     = sqrt(mean((y_seg_match - var_yield_kg_ha)^2)),
    nrmse    = rmse / (max(var_yield_kg_ha) - min(var_yield_kg_ha)) * 100,
    mean_obs = mean(var_yield_kg_ha),
    d_index  = 1 - sum((y_seg_match - mean_obs)^2) /
      sum((abs(y_seg_match - mean_obs) + abs(var_yield_kg_ha - mean_obs))^2),
    .groups  = "drop"
  ) %>%
  mutate(
    stats_lab = sprintf("RMSE = %.1f kg/ha\nnRMSE = %.1f%%\nd = %.2f",
                        rmse, nrmse, d_index),
    x_pos = 100,
    y_pos = 9000
  )

df_plot_all <- df_plot_all %>%
  mutate(
    beta_fac   = factor(beta, levels = betas,
                        labels = paste0("β = ", betas)),
    radius_fac = factor(radius_km, levels = c(20, 30, 40),
                        labels = c("20 km", "30 km", "40 km"))
  )

stats_df <- stats_df %>%
  mutate(
    beta_fac   = factor(beta, levels = betas,
                        labels = paste0("β = ", betas)),
    radius_fac = factor(radius_km, levels = c(20, 30, 40),
                        labels = c("20 km", "30 km", "40 km"))
  )

# ------------------------------------------------------------------
# 10. Plot: 1:1 comparisons in a beta × radius panel
# ------------------------------------------------------------------

p_all <- ggplot(df_plot_all,
                aes(x = var_yield_kg_ha, y = y_seg_match)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed(xlim = c(0, 9000), ylim = c(0, 9000)) +
  geom_text(
    data = stats_df,
    aes(x = x_pos, y = y_pos, label = stats_lab),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size  = 3.5,
    color = "blue",
    fontface = "plain"
  ) +
  facet_grid(radius_fac ~ beta_fac) +
  labs(
    x = "Variety trial yield (kg/ha)",
    y = "NASS-anchored cluster yield (kg/ha)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

p_all


