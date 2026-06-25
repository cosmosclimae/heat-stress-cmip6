# ============================================================
# Figure_Supp.r  — Standalone script to generate ALL SI outputs
# For: TAC paper (Heat stress + NEX-GDDP-CMIP6)
#
# OUTPUTS (in dir_out):
#   - SI_FigS1_range_MEAN_WBGT_Tw.png        (ensmax - ensmin; mean indices; pseudoHist 1991–2020)
#   - SI_FigS3_peak_SSP126.png              (Fig3 remake for SSP1-2.6; present vs 2100)
#   - SI_FigS3_peak_SSP245.png              (Fig3 remake for SSP2-4.5; present vs 2100)
#   - SI_FigS4_chronic_SSP126.png           (Fig4 remake for SSP1-2.6; 2050 & 2100; Tw35 + WBGT32)
#   - SI_FigS4_chronic_SSP245.png           (Fig4 remake for SSP2-4.5; 2050 & 2100; Tw35 + WBGT32)
#   - SI_TableS1_population_ranges.csv      (exposure table with ensmean/min/max, millions)
#   - SI_TableS1_population_ranges.tex      (LaTeX table body; plug into your tabularx)
#
# NOTES
# - You MUST edit the PATHS section to match your local structure.
# - Assumes NetCDF has monthly layers with "time" variable in CF units.
# - Chronic definition = >= 30 days/yr (same as paper).
# - Population files: user indicated pattern like "SSP1_2050_corr.tif" (note _corr).
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(cowplot)
  library(grid)
  library(lubridate)
  library(dplyr)
  library(tidyr)
})

# ----------------------------
# PATHS (EDIT THESE)
# ----------------------------
dir_out <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_SI"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# --- precomputed climatologies (from precompute_climatos.sh) ---
DIR_CLIM <- "C:/Data/HEAT_STRESS_MONTHLY_V2/CLIMATOS"
clim_tag <- function(ncfile) sub("\\.nc$", "", basename(ncfile))
read_mean_clim <- function(ncfile, var)
  rast(file.path(DIR_CLIM, sprintf("MEAN_%s__%s.nc", clim_tag(ncfile), var)))
read_days_clim <- function(ncfile, var)
  rast(file.path(DIR_CLIM, sprintf("DAYS_%s__%s.nc", clim_tag(ncfile), var)))

# Present (pseudoHist) ensemble files (needed for SI Fig S1 and "Present" for S3)
f_pres_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
f_pres_ensmax  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmax_pseudoHist_1991-2020.nc"
f_pres_ensmin  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmin_pseudoHist_1991-2020.nc"

# Future ensemble-mean files (SSP126/SSP245)
f_126_2050_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2041-2070/HEATSTRESS_ensmean_SSP126_2041-2070.nc"
f_126_2100_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2071-2100/HEATSTRESS_ensmean_SSP126_2071-2100.nc"
f_245_2050_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2041-2070/HEATSTRESS_ensmean_SSP245_2041-2070.nc"
f_245_2100_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2071-2100/HEATSTRESS_ensmean_SSP245_2071-2100.nc"

# OPTIONAL: if you have ensmin/ensmax per scenario/period for population ranges.
# If you don't, set to NULL and the script will fall back to model-range not available.
f_126_2050_ensmin <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2041-2070/HEATSTRESS_ensmin_SSP126_2041-2070.nc"
f_126_2050_ensmax <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2041-2070/HEATSTRESS_ensmax_SSP126_2041-2070.nc"
f_126_2100_ensmin <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2071-2100/HEATSTRESS_ensmin_SSP126_2071-2100.nc"
f_126_2100_ensmax <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2071-2100/HEATSTRESS_ensmax_SSP126_2071-2100.nc"

f_245_2050_ensmin <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2041-2070/HEATSTRESS_ensmin_SSP245_2041-2070.nc"
f_245_2050_ensmax <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2041-2070/HEATSTRESS_ensmax_SSP245_2041-2070.nc"
f_245_2100_ensmin <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2071-2100/HEATSTRESS_ensmin_SSP245_2071-2100.nc"
f_245_2100_ensmax <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2071-2100/HEATSTRESS_ensmax_SSP245_2071-2100.nc"

# If you also want SSP585 in Table S1, add these (and pop files below)
f_585_2050_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2041-2070/HEATSTRESS_ensmean_SSP585_2041-2070.nc"
f_585_2100_ensmean <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmean_SSP585_2071-2100.nc"
f_585_2050_ensmin  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2041-2070/HEATSTRESS_ensmin_SSP585_2041-2070.nc"
f_585_2050_ensmax  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2041-2070/HEATSTRESS_ensmax_SSP585_2041-2070.nc"
f_585_2100_ensmin  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmin_SSP585_2071-2100.nc"
f_585_2100_ensmax  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmax_SSP585_2071-2100.nc"

# Population rasters (0.25°) — user note: *_corr.tif
dir_pop <- "Z:/DATA_SOC/Pixel/Population"  # EDIT
pop_files <- list(
  "SSP1-2.6" = list(
    "2041-2070" = file.path(dir_pop, "SSP1/SSP1_2050_corr.tif"),
    "2071-2100" = file.path(dir_pop, "SSP1/SSP1_2100_corr.tif")
  ),
  "SSP2-4.5" = list(
    "2041-2070" = file.path(dir_pop, "SSP2/SSP2_2050_corr.tif"),
    "2071-2100" = file.path(dir_pop, "SSP2/SSP2_2100_corr.tif")
  ),
  "SSP5-8.5" = list(
    "2041-2070" = file.path(dir_pop, "SSP5/SSP5_2050_corr.tif"),
    "2071-2100" = file.path(dir_pop, "SSP5/SSP5_2100_corr.tif")
  )
)

# ----------------------------
# SETTINGS
# ----------------------------
CHRONIC_DAYS <- 30

# Mean indices for SI Fig S1
v_wbgt_mean <- "wbgt_shade_mean"
v_tw_mean   <- "tw_mean"

# Peak indices for SI Fig S3
vars_peak <- c("wbgt_shade_peak", "tw_peak")
nice_peak <- c(wbgt_shade_peak="WBGT", tw_peak="Tw")

# Exceedance-days variables (monthly sums) for chronic maps + exposure
v_tw35_days   <- "ndays_tw_peak_ge35"
v_wbgt32_days <- "ndays_wbgt_shade_peak_ge32"

# ----------------------------
# BASEMAP + CONTINENTS
# ----------------------------
world_ll <- ne_countries(scale = "medium", returnclass = "sf")
world_54030 <- st_transform(world_ll, 54030)

# continents dissolve (NaturalEarth "continent")
cont_ll <- world_ll |> st_make_valid() |> select(continent, geometry) |> filter(!is.na(continent))
cont_ll <- cont_ll |> group_by(continent) |> summarise(geometry = st_union(geometry), .groups="drop")
cont_54030 <- st_transform(cont_ll, 54030)

# Map to your table regions
continent_map <- c(
  "Africa" = "Africa",
  "Asia" = "Asia",
  "Europe" = "Europe",
  "North America" = "North America",
  "Oceania" = "Oceania",
  "South America" = "South America"
)

# ============================================================
# HELPERS
# ============================================================

stop_if_missing <- function(f) {
  if (is.null(f)) return(invisible(TRUE))
  if (!file.exists(f)) stop("Missing file: ", f)
  invisible(TRUE)
}

read_var <- function(ncfile, var) {
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0, nrow(r) > 0, ncol(r) > 0)
  r
}

to_m180_180 <- function(r) {
  if (xmax(r) > 180) rotate(r) else r
}

nc_time_ymd <- function(ncfile) {
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  t <- ncvar_get(nc, "time")
  u <- ncatt_get(nc, "time", "units")$value
  origin <- sub(".*since\\s+", "", u)
  origin <- strsplit(origin, " ")[[1]][1]
  d <- as.Date(t, origin = origin)
  data.frame(
    idx   = seq_along(d),
    date  = d,
    year  = as.integer(format(d, "%Y")),
    month = as.integer(format(d, "%m")),
    ndays = lubridate::days_in_month(d)
  )
}

clim_weighted_monthly_mean <- function(ncfile, var) {
  read_mean_clim(ncfile, var) |> to_m180_180()
}

# monthly sums -> yearly sums -> mean across years -> now reads DAYS_*.nc
mean_annual_days <- function(ncfile, var) {
  read_days_clim(ncfile, var) |> to_m180_180()
}

align_to <- function(src, tgt, method) {
  src <- to_m180_180(src)
  tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap="out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method = method)
}

to_df <- function(r, nm="value") {
  as.data.frame(r, xy=TRUE, na.rm=TRUE) |> setNames(c("x","y",nm))
}

# plot helpers
style_theme <- function(show_legend=FALSE){
  theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=14, face="bold"),
      legend.text  = element_text(size=13),
      legend.key.height = unit(4.2, "mm"),
      legend.key.width  = unit(42, "mm"),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

scale_range <- function(vmax, legend_title){
  scale_fill_viridis_c(option="cividis", limits=c(0, vmax), oob=squish, name=legend_title)
}

extract_legend <- function(p){
  cowplot::get_legend(p + theme(plot.title = element_blank()))
}

make_range_map <- function(r, title, vmax, legend_title, show_legend=FALSE, world=world_54030) {
  rr <- project(r, "ESRI:54030")
  df <- to_df(rr, "value")
  resx <- res(rr)[1]; resy <- res(rr)[2]
  ggplot() +
    geom_tile(data=df, aes(x, y, fill=value), width=resx, height=resy) +
    geom_sf(data=world, fill=NA, color="grey20", linewidth=0.2) +
    coord_sf(crs=54030, expand=FALSE) +
    scale_range(vmax, legend_title) +
    labs(title=title) +
    style_theme(show_legend)
}

coldhot_cols <- c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")

guide_steps <- function() {
  guide_coloursteps(
    title.position = "top",
    title.hjust = 0.5,
    barwidth  = unit(42, "mm"),
    barheight = unit(4.4, "mm"),
    ticks = TRUE,
    even.steps = TRUE
  )
}

breaks_pretty <- function(lims, n=4) pretty(lims, n=n)

make_steps_map <- function(r, title, lims, breaks, legend_title, show_legend=FALSE, world=world_54030) {
  rr <- project(r, "ESRI:54030")
  df <- to_df(rr, "val")
  resx <- res(rr)[1]; resy <- res(rr)[2]
  ggplot() +
    geom_tile(data=df, aes(x, y, fill=val), width=resx, height=resy) +
    geom_sf(data=world, fill=NA, color="grey20", linewidth=0.2) +
    coord_sf(crs = 54030, expand = FALSE) +
    scale_fill_stepsn(
      colors = coldhot_cols,
      limits = lims,
      breaks = breaks,
      oob = squish,
      name = legend_title,
      guide = guide_steps()
    ) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=13, face="bold"),
      legend.text  = element_text(size=12),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

row_header <- function(txt){
  wrap_elements(full = textGrob(txt, gp=gpar(fontface="bold", fontsize=18))) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
}

# Chronic binary helper
bin01 <- function(r, thr = CHRONIC_DAYS) {
  classify(r, rcl = matrix(c(-Inf, thr, 0,
                            thr,  Inf, 1), ncol=3, byrow=TRUE))
}

cols_bin <- c("0"="grey92", "1"="#d7191c")

style_theme_bin <- function(show_legend=FALSE){
  theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=13, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=14, face="bold"),
      legend.text  = element_text(size=13),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width  = unit(18, "mm"),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

make_bin_map <- function(r_bin01, title, legend_title=NULL, show_legend=FALSE) {
  rr <- project(r_bin01, "ESRI:54030")
  df <- to_df(rr, "z")
  resx <- res(rr)[1]; resy <- res(rr)[2]

  if (is.null(legend_title)) legend_title <- paste0("Chronic exposure (≥", CHRONIC_DAYS, " days/yr)")

  ggplot() +
    geom_tile(data=df, aes(x, y, fill=factor(z)), width=resx, height=resy) +
    geom_sf(data=world_54030, fill=NA, color="grey20", linewidth=0.2) +
    coord_sf(crs=54030, expand=FALSE) +
    scale_fill_manual(
      values=cols_bin,
      breaks=c("1","0"),
      labels=c(paste0("≥", CHRONIC_DAYS, " days/yr"), paste0("<", CHRONIC_DAYS, " days/yr")),
      name=legend_title
    ) +
    labs(title=title) +
    style_theme_bin(show_legend)
}

# ============================================================
# SI FIG S1 — Ensemble range on MEAN WBGT & Tw (ensmax - ensmin)
# ============================================================

message("SI Fig S1 — range on mean WBGT/Tw (ensmax - ensmin)")

stop_if_missing(f_pres_ensmax)
stop_if_missing(f_pres_ensmin)

wbgt_max <- clim_weighted_monthly_mean(f_pres_ensmax, v_wbgt_mean)
wbgt_min <- clim_weighted_monthly_mean(f_pres_ensmin, v_wbgt_mean)
tw_max   <- clim_weighted_monthly_mean(f_pres_ensmax, v_tw_mean)
tw_min   <- clim_weighted_monthly_mean(f_pres_ensmin, v_tw_mean)

wbgt_rng <- wbgt_max - wbgt_min
tw_rng   <- tw_max - tw_min

vmax_rng <- max(
  quantile(values(wbgt_rng, mat=FALSE), 0.995, na.rm=TRUE),
  quantile(values(tw_rng,   mat=FALSE), 0.995, na.rm=TRUE),
  na.rm=TRUE
)

pS1a <- make_range_map(wbgt_rng, "WBGT mean", vmax_rng, "ensmax − ensmin (°C)", show_legend=FALSE)
pS1b <- make_range_map(tw_rng,   "Tw mean",   vmax_rng, "ensmax − ensmin (°C)", show_legend=FALSE)

legS1 <- extract_legend(make_range_map(wbgt_rng, NULL, vmax_rng, "ensmax − ensmin (°C)", show_legend=TRUE))

figS1 <- (pS1a / pS1b) / wrap_elements(full = legS1) +
  plot_layout(heights = c(1, 1, 0.12))

ggsave(file.path(dir_out, "SI_FigS1_range_MEAN_WBGT_Tw.png"), figS1, width=10.5, height=9.2, dpi=300)


# ============================================================
# SI FIG S3 — Fig3 remake for SSP126 / SSP245 (PEAK indices)
# ============================================================

make_figS3_peak <- function(f_present, f_future, scen_label, out_png){
  stop_if_missing(f_present); stop_if_missing(f_future)

  pres_peak <- setNames(lapply(vars_peak, \(v) clim_weighted_monthly_mean(f_present, v)), vars_peak)
  fut_peak  <- setNames(lapply(vars_peak, \(v) clim_weighted_monthly_mean(f_future,  v)), vars_peak)

  fut_peak_a <- mapply(align_to, fut_peak, pres_peak, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)

  lims_peak <- lapply(vars_peak, function(v){
    a <- values(pres_peak[[v]], mat=FALSE); a <- a[is.finite(a)]
    b <- values(fut_peak_a[[v]], mat=FALSE); b <- b[is.finite(b)]
    vv <- c(a,b)
    as.numeric(quantile(vv, probs=c(0.005, 0.995), na.rm=TRUE))
  })
  names(lims_peak) <- vars_peak
  breaks_peak <- lapply(vars_peak, function(v) breaks_pretty(lims_peak[[v]], n=4))
  names(breaks_peak) <- vars_peak

  maps_pres <- lapply(vars_peak, function(v){
    make_steps_map(
      pres_peak[[v]], title = nice_peak[[v]],
      lims = lims_peak[[v]], breaks = breaks_peak[[v]],
      legend_title = paste0(nice_peak[[v]], " (°C)"),
      show_legend = FALSE
    )
  })
  maps_fut <- lapply(vars_peak, function(v){
    make_steps_map(
      fut_peak_a[[v]], title = nice_peak[[v]],
      lims = lims_peak[[v]], breaks = breaks_peak[[v]],
      legend_title = paste0(nice_peak[[v]], " (°C)"),
      show_legend = FALSE
    )
  })

  legends <- lapply(vars_peak, function(v){
    p_leg <- make_steps_map(
      fut_peak_a[[v]], title=NULL,
      lims=lims_peak[[v]], breaks=breaks_peak[[v]],
      legend_title=paste0(nice_peak[[v]], " (°C)"),
      show_legend=TRUE
    ) + theme(plot.title = element_blank())
    cowplot::get_legend(p_leg)
  })

  legend_row <- wrap_elements(full = legends[[1]]) | wrap_elements(full = legends[[2]])

  row1 <- (maps_pres[[1]] | maps_pres[[2]])
  row2 <- (maps_fut[[1]]  | maps_fut[[2]])

  fig <- row_header("Present (1991–2020)") /
    row1 /
    row_header(paste0(scen_label, " (2071–2100)")) /
    row2 /
    legend_row +
    plot_layout(heights = c(0.08, 1, 0.08, 1, 0.18))

  ggsave(file.path(dir_out, out_png), fig, width=10.5, height=8.8, dpi=300)
  invisible(fig)
}

message("SI Fig S3 — peaks SSP126/SSP245")
make_figS3_peak(f_pres_ensmean, f_126_2100_ensmean, "SSP1-2.6", "SI_FigS3_peak_SSP126.png")
make_figS3_peak(f_pres_ensmean, f_245_2100_ensmean, "SSP2-4.5", "SI_FigS3_peak_SSP245.png")

# ============================================================
# SI FIG S4 — Fig4 remake for SSP126 / SSP245 (CHRONIC exceedance)
# ============================================================

make_figS4_chronic <- function(f_2050, f_2100, scen_label, out_png){
  stop_if_missing(f_2050); stop_if_missing(f_2100)

  # target grid = present days grid to keep stable layout
  tgt <- mean_annual_days(f_pres_ensmean, v_tw35_days)

  # compute days
  tw50   <- mean_annual_days(f_2050, v_tw35_days)   |> align_to(tgt, method="near")
  wb50   <- mean_annual_days(f_2050, v_wbgt32_days) |> align_to(tgt, method="near")
  tw00   <- mean_annual_days(f_2100, v_tw35_days)   |> align_to(tgt, method="near")
  wb00   <- mean_annual_days(f_2100, v_wbgt32_days) |> align_to(tgt, method="near")

  # binary
  tw50b <- bin01(tw50)
  wb50b <- bin01(wb50)
  tw00b <- bin01(tw00)
  wb00b <- bin01(wb00)

  p1 <- make_bin_map(tw50b, "Tw ≥ 35°C")
  p2 <- make_bin_map(wb50b, "WBGT ≥ 32°C")
  p3 <- make_bin_map(tw00b, "Tw ≥ 35°C")
  p4 <- make_bin_map(wb00b, "WBGT ≥ 32°C")

  leg <- extract_legend(make_bin_map(tw50b, NULL, show_legend=TRUE))
  legend_row <- wrap_elements(full = leg)

  fig <- row_header(paste0("2050 (2041–2070) — ", scen_label)) /
    (p1 | p2) /
    row_header(paste0("2100 (2071–2100) — ", scen_label)) /
    (p3 | p4) /
    legend_row +
    plot_layout(heights = c(0.07, 1, 0.07, 1, 0.12))

  ggsave(file.path(dir_out, out_png), fig, width=10.5, height=9.2, dpi=300)
  invisible(fig)
}

message("SI Fig S4 — chronic SSP126/SSP245")
make_figS4_chronic(f_126_2050_ensmean, f_126_2100_ensmean, "SSP1-2.6", "SI_FigS4_chronic_SSP126.png")
make_figS4_chronic(f_245_2050_ensmean, f_245_2100_ensmean, "SSP2-4.5", "SI_FigS4_chronic_SSP245.png")

# ============================================================
# SI TABLE S1 — Population exposure with ensmean/ensmin/ensmax
# ============================================================

# Prepare a continent raster on target grid
make_continent_raster <- function(tgt_latlon){
  # build in 54030 for good area behavior in plotting,
  # but for zonal sums we just need consistent IDs on the same grid.
  # We'll rasterize in lat/lon directly to match your NEX/POP grid.
  cont <- cont_ll
  cont_vect <- vect(cont)  # terra vector (lat/lon)
  # rasterize with continent codes (as character -> we will map later)
  # terra rasterize needs numeric field; create an integer ID
  cont$cid <- as.integer(factor(cont$continent))
  cont_vect <- vect(cont)
  cr <- rast(tgt_latlon); values(cr) <- NA
  r_cont <- rasterize(cont_vect, cr, field="cid", touches=TRUE)
  # store lookup
  lut <- data.frame(cid=unique(cont$cid), continent=unique(cont$continent))
  list(r_cont=r_cont, lut=lut)
}

sum_pop_mask_global <- function(pop_r, mask01){
  pop_r <- to_m180_180(pop_r)
  mask01 <- to_m180_180(mask01)
  pop_r <- align_to(pop_r, mask01, method="bilinear")
  v <- values(pop_r * mask01, mat=FALSE)
  sum(v, na.rm=TRUE)
}

sum_pop_mask_by_continent <- function(pop_r, mask01, r_cont, lut){
  pop_r <- to_m180_180(pop_r)
  mask01 <- to_m180_180(mask01)
  pop_r <- align_to(pop_r, mask01, method="bilinear")
  r_cont <- align_to(r_cont, mask01, method="near")

  x <- pop_r * mask01
  z <- zonal(x, r_cont, fun="sum", na.rm=TRUE)
  colnames(z) <- c("cid","pop")

  z <- left_join(z, lut, by="cid") |> select(continent, pop)
  z
}

exposure_from_nc_pop <- function(ncfile, popfile, var_days){
  stop_if_missing(ncfile)
  stop_if_missing(popfile)

  pop <- rast(popfile) |> to_m180_180()

  days <- mean_annual_days(ncfile, var_days)
  chronic <- bin01(days)  # 0/1

  # continent raster built on chronic grid
  cc <- make_continent_raster(chronic)
  r_cont <- cc$r_cont
  lut <- cc$lut

  g <- sum_pop_mask_global(pop, chronic) / 1e6
  z <- sum_pop_mask_by_continent(pop, chronic, r_cont, lut) |> mutate(pop = pop/1e6)

  list(global=g, by_continent=z)
}

# Scenario registry for Table S1
scenarios <- list(
  "SSP1-2.6" = list(
    "2041-2070" = list(mean=f_126_2050_ensmean, min=f_126_2050_ensmin, max=f_126_2050_ensmax),
    "2071-2100" = list(mean=f_126_2100_ensmean, min=f_126_2100_ensmin, max=f_126_2100_ensmax)
  ),
  "SSP2-4.5" = list(
    "2041-2070" = list(mean=f_245_2050_ensmean, min=f_245_2050_ensmin, max=f_245_2050_ensmax),
    "2071-2100" = list(mean=f_245_2100_ensmean, min=f_245_2100_ensmin, max=f_245_2100_ensmax)
  ),
  "SSP5-8.5" = list(
    "2041-2070" = list(mean=f_585_2050_ensmean, min=f_585_2050_ensmin, max=f_585_2050_ensmax),
    "2071-2100" = list(mean=f_585_2100_ensmean, min=f_585_2100_ensmin, max=f_585_2100_ensmax)
  )
)

# If you don't want SSP5-8.5 in SI Table S1, comment it out above.

# Compute Table S1
message("SI Table S1 — population exposure ranges (ensmean/min/max)")

regions_order <- c("Global","Africa","Asia","Europe","North America","Oceania","South America")

compute_block <- function(scen_label, period_label, nc_triplet, popfile, var_days, metric_name){
  # returns rows for Global + continents, with ensmean/min/max if available
  get_one <- function(tag, ncfile){
    if (is.null(ncfile) || !file.exists(ncfile)) return(NULL)
    ex <- exposure_from_nc_pop(ncfile, popfile, var_days)
    # Global
    out_g <- data.frame(region="Global", value=ex$global, tag=tag)
    # Continents
    out_c <- ex$by_continent |>
      filter(continent %in% names(continent_map)) |>
      mutate(region = continent_map[continent]) |>
      transmute(region, value=pop, tag=tag)
    bind_rows(out_g, out_c)
  }

  a_mean <- get_one("ensmean", nc_triplet$mean)
  a_min  <- get_one("ensmin",  nc_triplet$min)
  a_max  <- get_one("ensmax",  nc_triplet$max)

  df <- bind_rows(a_mean, a_min, a_max)
  if (nrow(df) == 0) return(NULL)

  df |>
    mutate(
      scenario = scen_label,
      period   = period_label,
      metric   = metric_name
    )
}

all_rows <- list()

for (sc in names(scenarios)) {
  for (per in names(scenarios[[sc]])) {
    popfile <- pop_files[[sc]][[per]]
    # WBGT32
    all_rows[[length(all_rows)+1]] <- compute_block(
      sc, per, scenarios[[sc]][[per]], popfile,
      v_wbgt32_days, "WBGT$_{\\geq32}$"
    )
    # Tw35
    all_rows[[length(all_rows)+1]] <- compute_block(
      sc, per, scenarios[[sc]][[per]], popfile,
      v_tw35_days, "Tw$_{\\geq35}$"
    )
  }
}

tab <- bind_rows(all_rows) |>
  filter(!is.na(region)) |>
  mutate(region = factor(region, levels=regions_order))

if (nrow(tab) == 0) {
  warning("No rows in Table S1. Check NC/pop file paths, and especially ensmin/ensmax scenario files.")
} else {
  # Pivot to ensmean/min/max columns
  tab_w <- tab |>
    tidyr::pivot_wider(names_from = tag, values_from = value) |>
    arrange(metric, region, scenario, period)

  # Write CSV
  write.csv(tab_w, file.path(dir_out, "SI_TableS1_population_ranges.csv"), row.names = FALSE)

  # Write a LaTeX-friendly table body (you can wrap it in your tabularx)
  # Columns: Region | Metric | Scenario | Period | ensmean | ensmin | ensmax
  tex_lines <- tab_w |>
    mutate(
      ensmean = ifelse(is.na(ensmean), "", sprintf("%.1f", ensmean)),
      ensmin  = ifelse(is.na(ensmin),  "", sprintf("%.1f", ensmin)),
      ensmax  = ifelse(is.na(ensmax),  "", sprintf("%.1f", ensmax))
    ) |>
    transmute(
      line = paste0(
        as.character(region), " & ",
        metric, " & ",
        scenario, " & ",
        period, " & ",
        ensmean, " & ",
        ensmin, " & ",
        ensmax, " \\\\"
      )
    )

  writeLines(tex_lines$line, con = file.path(dir_out, "SI_TableS1_population_ranges.tex"))
}

message("DONE — SI outputs written to: ", dir_out)

