# ============================================================
# FIG1 + SI1 + SI2 — "Figure 2 style"
# Legends extracted and placed in a dedicated bottom row
# ============================================================

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


# ----------------------------
# Paths
# ----------------------------
dir_ens <- "C:/Data/HEAT_STRESS_MONTHLY/ENSEMBLE/pseudoHist/1991-2020"
dir_era <- "C:/Data/HEAT_STRESS_MONTHLY/ERA5/packs"
dir_fig <- "C:/Data/HEAT_STRESS_MONTHLY/Figures"
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

f_ens <- file.path(dir_ens, "HEATSTRESS_ensmean_pseudoHist_1991-2020.nc")
f_max <- file.path(dir_ens, "HEATSTRESS_ensmax_pseudoHist_1991-2020.nc")
f_min <- file.path(dir_ens, "HEATSTRESS_ensmin_pseudoHist_1991-2020.nc")
f_era <- file.path(dir_era, "HEATSTRESS_mon_ERA5_gr025_1991-2020.nc")

vars_mean <- c("wbgt_mean", "tw_mean", "hi_mean", "humidex_mean")
vars_days <- c(wbgt30="ndays_wbgt_ge30", wbgt32="ndays_wbgt_ge32", tw35="ndays_tw_ge35")

nice_mean <- c(wbgt_mean="WBGT", tw_mean="Tw", hi_mean="Heat Index", humidex_mean="Humidex")

# ----------------------------
# Basemap
# ----------------------------
world_ll    <- ne_countries(scale = "medium", returnclass = "sf")
world_54030 <- st_transform(world_ll, 54030)  # ESRI:54030

# ----------------------------
# Helpers: read + time
# ----------------------------
read_var <- function(ncfile, var) {
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0, nrow(r) > 0, ncol(r) > 0)
  r
}

nc_time_df <- function(ncfile) {
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  t <- ncvar_get(nc, "time")
  u <- ncatt_get(nc, "time", "units")$value
  origin <- sub(".*since\\s+", "", u)
  origin <- strsplit(origin, " ")[[1]][1]
  d <- as.Date(t, origin = origin)
  data.frame(idx = seq_along(d), year = as.integer(format(d, "%Y")))
}

to_m180_180 <- function(r) {
  if (xmax(r) > 180) rotate(r) else r
}

# ----------------------------
# Climatologies (1 layer)
# ----------------------------
clim_mean <- function(ncfile, var) {
  r <- read_var(ncfile, var)
  mean(r, na.rm = TRUE) |> to_m180_180()
}

mean_annual_days <- function(ncfile, var) {
  r  <- read_var(ncfile, var)
  td <- nc_time_df(ncfile)
  stopifnot(nlyr(r) == nrow(td))
  yrs <- unique(td$year)
  annual <- lapply(yrs, function(y) {
    idx <- td$idx[td$year == y]
    sum(r[[idx]], na.rm = TRUE)
  })
  mean(rast(annual), na.rm = TRUE) |> to_m180_180()
}

# ----------------------------
# Align: src -> tgt
# ----------------------------
align_to <- function(src, tgt, method) {
  src <- to_m180_180(src)
  tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap="out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method = method)
}

# ----------------------------
# Robust limits
# ----------------------------
sym_lim_q <- function(r, q = 0.995) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) == 0) return(1)
  as.numeric(quantile(abs(v), probs = q, na.rm = TRUE))
}

to_df <- function(r, nm="value") {
  as.data.frame(r, xy=TRUE, na.rm=TRUE) |> setNames(c("x","y",nm))
}

# ----------------------------
# Common style (titles + legends)
# ----------------------------
style_theme <- function(show_legend=FALSE){
  theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=14, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=9, face="bold"),
      legend.text  = element_text(size=8),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width  = unit(38, "mm"),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

# Bias palette (diverging)
scale_bias <- function(lim, legend_title){
  scale_fill_gradient2(
    low="#2166ac", mid="white", high="#b2182b",
    midpoint=0, limits=c(-lim, lim), oob=squish, name=legend_title
  )
}

# Range palette (sequential, readable)
scale_range <- function(vmax, legend_title){
  scale_fill_viridis_c(option="cividis", limits=c(0, vmax), oob=squish, name=legend_title)
}

# ----------------------------
# Map builders (Fig2-style: optional legend)
# ----------------------------
make_bias_map <- function(r, title, lim, legend_title, show_legend=FALSE, world=world_54030) {
  rr <- project(r, "ESRI:54030")
  df <- to_df(rr, "bias")
  resx <- res(rr)[1]; resy <- res(rr)[2]

  ggplot() +
    geom_tile(data=df, aes(x, y, fill=bias), width=resx, height=resy) +
    geom_sf(data=world, fill=NA, color="grey20", linewidth=0.2) +
    coord_sf(crs=54030, expand=FALSE) +
    scale_bias(lim, legend_title) +
    labs(title=title) +
    style_theme(show_legend)
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

extract_legend <- function(p){
  cowplot::get_legend(p + theme(plot.title = element_blank()))
}

# ============================================================
# 1) Compute climatologies (1 layer)
# ============================================================
ens_mean <- setNames(lapply(vars_mean, \(v) clim_mean(f_ens, v)), vars_mean)
era_mean <- setNames(lapply(vars_mean, \(v) clim_mean(f_era, v)), vars_mean)

ens_days <- lapply(vars_days, \(v) mean_annual_days(f_ens, v))
era_days <- lapply(vars_days, \(v) mean_annual_days(f_era, v))

ensmax_days <- lapply(vars_days, \(v) mean_annual_days(f_max, v))
ensmin_days <- lapply(vars_days, \(v) mean_annual_days(f_min, v))

# ============================================================
# 2) Align ERA5 to ENS (ENS target grid)
# ============================================================
era_mean_a <- mapply(align_to, era_mean, ens_mean, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)
era_days_a <- mapply(align_to, era_days, ens_days, MoreArgs=list(method="near"),     SIMPLIFY=FALSE)

# ============================================================
# 3) Bias + range
# ============================================================
bias_mean  <- mapply(`-`, ens_mean, era_mean_a, SIMPLIFY=FALSE)
bias_days  <- mapply(`-`, ens_days, era_days_a, SIMPLIFY=FALSE)
range_days <- mapply(`-`, ensmax_days, ensmin_days, SIMPLIFY=FALSE)

# ============================================================
# FIG 1 — Bias on exceedance days (3 panels + legend row)
# ============================================================
lim1 <- max(sym_lim_q(bias_days$wbgt30), sym_lim_q(bias_days$wbgt32), sym_lim_q(bias_days$tw35))
leg_title1 <- "ENS − ERA5 (days/yr)"

p1 <- make_bias_map(bias_days$wbgt30, "WBGT ≥ 30°C", lim1, leg_title1, show_legend=FALSE)
p2 <- make_bias_map(bias_days$wbgt32, "WBGT ≥ 32°C", lim1, leg_title1, show_legend=FALSE)
p3 <- make_bias_map(bias_days$tw35,   "Tw ≥ 35°C",   lim1, leg_title1, show_legend=FALSE)

row_maps1 <- (p1 | p2 | p3)
leg1 <- extract_legend(make_bias_map(bias_days$wbgt30, NULL, lim1, leg_title1, show_legend=TRUE))
fig1 <- row_maps1 / wrap_elements(full = leg1) + plot_layout(heights = c(1, 0.12))

ggsave(file.path(dir_fig, "Fig1_bias_thresholds.png"), fig1, width=12, height=4.6, dpi=300)

# ============================================================
# SI1 — Bias on annual mean indices (2×2 + legend row) [robust]
# ============================================================
limS1 <- max(sym_lim_q(bias_mean$wbgt_mean), sym_lim_q(bias_mean$tw_mean),
             sym_lim_q(bias_mean$hi_mean),   sym_lim_q(bias_mean$humidex_mean))
leg_titleS1 <- "ENS − ERA5"

s1 <- make_bias_map(bias_mean$wbgt_mean,    nice_mean[["wbgt_mean"]],    limS1, leg_titleS1, FALSE)
s2 <- make_bias_map(bias_mean$tw_mean,      nice_mean[["tw_mean"]],      limS1, leg_titleS1, FALSE)
s3 <- make_bias_map(bias_mean$hi_mean,      nice_mean[["hi_mean"]],      limS1, leg_titleS1, FALSE)
s4 <- make_bias_map(bias_mean$humidex_mean, nice_mean[["humidex_mean"]], limS1, leg_titleS1, FALSE)

# 2×2 grid of maps (patchwork is fine here)
mapsS1 <- (s1 | s2) / (s3 | s4)

# extract legend (single)
legS1 <- extract_legend(make_bias_map(bias_mean$wbgt_mean, NULL, limS1, leg_titleS1, TRUE))

# Combine using cowplot (rock-solid for plot + legend layouts)
si1 <- cowplot::plot_grid(
  mapsS1,
  legS1,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggsave(file.path(dir_fig, "SI1_bias_means_4panel.png"),
       si1, width=10.5, height=7.6, dpi=300)

# ============================================================
# SI2 — Ensemble range (max−min) on exceedance days (3 panels + legend row)
# ============================================================
vmaxS2 <- max(
  quantile(values(range_days$wbgt30, mat=FALSE), 0.995, na.rm=TRUE),
  quantile(values(range_days$wbgt32, mat=FALSE), 0.995, na.rm=TRUE),
  quantile(values(range_days$tw35,   mat=FALSE), 0.995, na.rm=TRUE),
  na.rm=TRUE
)
leg_titleS2 <- "max − min (days/yr)"

r1 <- make_range_map(range_days$wbgt30, "WBGT ≥ 30°C", vmaxS2, leg_titleS2, FALSE)
r2 <- make_range_map(range_days$wbgt32, "WBGT ≥ 32°C", vmaxS2, leg_titleS2, FALSE)
r3 <- make_range_map(range_days$tw35,   "Tw ≥ 35°C",   vmaxS2, leg_titleS2, FALSE)

row_maps2 <- (r1 | r2 | r3)
legS2 <- extract_legend(make_range_map(range_days$wbgt30, NULL, vmaxS2, leg_titleS2, TRUE))
si2 <- row_maps2 / wrap_elements(full = legS2) + plot_layout(heights = c(1, 0.12))

ggsave(file.path(dir_fig, "SI2_range_thresholds.png"), si2, width=12, height=4.6, dpi=300)

message("DONE -> ", dir_fig)





# ============================================================
# FIGURE 2 — 2×4 maps + row titles + 1 legend per indicator
# Present (1991–2020) vs SSP5-8.5 (2071–2100)
# ============================================================

# ----------------------------
# PATHS (EDIT)
# ----------------------------
dir_fig <- "C:/Data/HEAT_STRESS_MONTHLY/Figures"
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# Your NEX-GDDP-CMIP6 ensemble files (EDIT)
f_pres <- "C:/Data/HEAT_STRESS_MONTHLY/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
f_fut  <- "C:/Data/HEAT_STRESS_MONTHLY/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmean_SSP585_2071-2100.nc"

vars_mean <- c("wbgt_mean", "tw_mean", "hi_mean", "humidex_mean")
nice_names <- c(wbgt_mean="WBGT", tw_mean="Tw", hi_mean="Heat Index", humidex_mean="Humidex")

# ----------------------------
# Helpers
# ----------------------------
read_var <- function(ncfile, var) {
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0, nrow(r) > 0, ncol(r) > 0)
  r
}

to_m180_180 <- function(r) {
  if (xmax(r) > 180) rotate(r) else r
}

nc_time_ymd <- function(ncfile) {
  nc <- nc_open(ncfile)
  on.exit(nc_close(nc))
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

# Weighted climatological mean for monthly means (1 layer out)
clim_weighted_monthly_mean <- function(ncfile, var) {
  r <- read_var(ncfile, var) |> to_m180_180()
  td <- nc_time_ymd(ncfile)
  stopifnot(nlyr(r) == nrow(td))
  w <- td$ndays

  out <- app(r, fun = function(...) {
    v <- c(...)
    ok <- is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * w[ok]) / sum(w[ok])
  })
  out
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

# ----------------------------
# World outlines (coasts/continents)
# ----------------------------
world_ll   <- rnaturalearth::ne_countries(scale="medium", returnclass="sf")
world_54030 <- sf::st_transform(world_ll, 54030)  # ESRI:54030

# ----------------------------
# Color scale: cold → hot (binned)
# ----------------------------
coldhot_cols <- c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")

guide_steps <- function() {
  guide_coloursteps(
    title.position = "top",
    title.hjust = 0.5,
    barwidth  = unit(34, "mm"),   # compact so 4 legends fit
    barheight = unit(3.2, "mm"),
    ticks = TRUE,
    even.steps = TRUE
  )
}

breaks_pretty <- function(lims, n=4) pretty(lims, n=n)

# Unit labels: keep °C for WBGT/Tw/HI; Humidex often treated as index (drop °C)
unit_lab <- function(v) if (v == "humidex_mean") "" else " (°C)"

# ----------------------------
# Plot (map only; legend optional)
# ----------------------------
make_map <- function(r, title, lims, breaks, legend_title, show_legend=FALSE) {
  rr <- project(r, "ESRI:54030")
  df <- to_df(rr, "val")
  resx <- res(rr)[1]; resy <- res(rr)[2]

  p <- ggplot() +
    geom_tile(data=df, aes(x, y, fill=val), width=resx, height=resy) +
    geom_sf(data=world_54030, fill=NA, color="grey20", linewidth=0.2) +
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
      plot.title = element_text(hjust=0.5, face="bold", size=14, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=9, face="bold"),
      legend.text  = element_text(size=8),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )

  p
}

# ----------------------------
# Row header grob (tight)
# ----------------------------
row_header <- function(txt){
  wrap_elements(full = textGrob(txt, gp=gpar(fontface="bold", fontsize=15))) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
}

# ============================================================
# 1) Compute weighted climatological means
# ============================================================
pres_mean <- setNames(lapply(vars_mean, \(v) clim_weighted_monthly_mean(f_pres, v)), vars_mean)
fut_mean  <- setNames(lapply(vars_mean, \(v) clim_weighted_monthly_mean(f_fut,  v)), vars_mean)

# Align future to present grid (bilinear for continuous means)
fut_mean_a <- mapply(align_to, fut_mean, pres_mean, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)

# Per-index limits across both periods (robust 0.5–99.5%)
lims_idx <- lapply(vars_mean, function(v){
  a <- values(pres_mean[[v]], mat=FALSE); a <- a[is.finite(a)]
  b <- values(fut_mean_a[[v]], mat=FALSE); b <- b[is.finite(b)]
  vv <- c(a,b)
  qs <- quantile(vv, probs=c(0.005, 0.995), na.rm=TRUE)
  as.numeric(qs)
})
names(lims_idx) <- vars_mean

# Per-index breaks (few ticks)
breaks_idx <- lapply(vars_mean, function(v) breaks_pretty(lims_idx[[v]], n=4))
names(breaks_idx) <- vars_mean

# ============================================================
# 2) Build maps (no legends) + extract legends separately
# ============================================================
maps_pres <- lapply(vars_mean, function(v){
  make_map(
    pres_mean[[v]],
    title = nice_names[[v]],
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_names[[v]], unit_lab(v)),
    show_legend = FALSE
  )
})

maps_fut <- lapply(vars_mean, function(v){
  make_map(
    fut_mean_a[[v]],
    title = nice_names[[v]],
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_names[[v]], unit_lab(v)),
    show_legend = FALSE
  )
})

# Legend extraction from a dummy plot per indicator (same scale)
legends <- lapply(vars_mean, function(v){
  p_leg <- make_map(
    fut_mean_a[[v]],
    title = NULL,
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_names[[v]], unit_lab(v)),
    show_legend = TRUE
  ) + theme(plot.title = element_blank())
  cowplot::get_legend(p_leg)
})

# Convert legend grobs into patchwork elements (keeps alignment)
legend_row <- wrap_elements(full = legends[[1]]) |
              wrap_elements(full = legends[[2]]) |
              wrap_elements(full = legends[[3]]) |
              wrap_elements(full = legends[[4]])

# ============================================================
# 3) Assemble figure: header + 2×4 maps + legend row
# ============================================================
row1 <- (maps_pres[[1]] | maps_pres[[2]] | maps_pres[[3]] | maps_pres[[4]])
row2 <- (maps_fut[[1]]  | maps_fut[[2]]  | maps_fut[[3]]  | maps_fut[[4]])

fig2 <- row_header("Present (1991–2020)") /
        row1 /
        row_header("SSP5-8.5 (2071–2100)") /
        row2 /
        legend_row +
        plot_layout(heights = c(0.06, 1, 0.06, 1, 0.16))

ggsave(file.path(dir_fig, "Fig2_means_2x4_SSP585.png"),
       fig2, width=16, height=8.8, dpi=300)

message("DONE: ", file.path(dir_fig, "Fig2_means_2x4_SSP585.png"))
