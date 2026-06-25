# ============================================================
# FIG 1–4 — FULL SCRIPT (V2 variables)
# - Fig1: ERA5 validation (bias on exceedance days; PEAK ndays vars)
# - SI1: bias on MEAN indices
# - SI2: ensemble range on exceedance days (max-min)
# - Fig2: MEAN climatologies (2×4, present vs end-century)
# - Fig3: PEAK climatologies (mean of daily peaks; 2×2, present vs end-century)
# - Fig4: Chronic exceedance (>=30 days/yr) + co-occurrence (present / 2050 / 2100)
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
# PATHS (EDIT)
# ----------------------------
dir_ens <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020"
dir_era <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ERA5/packs"
dir_fig <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_NEW"
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# --- precomputed climatologies (from precompute_climatos.sh) ---
# Heavy time-aggregation is done once in CDO; R just reads single-layer fields.
DIR_CLIM <- "C:/Data/HEAT_STRESS_MONTHLY_V2/CLIMATOS"
clim_tag <- function(ncfile) sub("\\.nc$", "", basename(ncfile))
read_mean_clim <- function(ncfile, var)
  rast(file.path(DIR_CLIM, sprintf("MEAN_%s__%s.nc", clim_tag(ncfile), var)))
read_days_clim <- function(ncfile, var)
  rast(file.path(DIR_CLIM, sprintf("DAYS_%s__%s.nc", clim_tag(ncfile), var)))

# ENS: present validation + range
f_ens <- file.path(dir_ens, "HEATSTRESS_ensmean_pseudoHist_1991-2020.nc")
f_max <- file.path(dir_ens, "HEATSTRESS_ensmax_pseudoHist_1991-2020.nc")
f_min <- file.path(dir_ens, "HEATSTRESS_ensmin_pseudoHist_1991-2020.nc")

# ERA5 validation
f_era <- file.path(dir_era, "HEATSTRESS_mon_ERA5_gr025_1991-2020.nc")

# Fig2/3/4 ENSMEAN files (EDIT)
f_pres <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
f_fut  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmean_SSP585_2071-2100.nc"

# Fig4 periods (EDIT)
f_2050 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2041-2070/HEATSTRESS_ensmean_SSP585_2041-2070.nc"
f_2100 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmean_SSP585_2071-2100.nc"

# ----------------------------
# VARIABLES (FIXED)
# ----------------------------
# Mean indices (Fig2, SI1)
vars_mean <- c("wbgt_shade_mean", "tw_mean", "hi_mean", "humidex_mean")
nice_mean <- c(
  wbgt_shade_mean="WBGT",
  tw_mean="Tw",
  hi_mean="Heat Index",
  humidex_mean="Humidex"
)

# Peak indices (Fig3 only)
vars_peak <- c("wbgt_shade_peak", "tw_peak")
nice_peak <- c(
  wbgt_shade_peak="WBGT",
  tw_peak="Tw"
)

# Exceedance days (Fig1, SI2, Fig4) — PEAK ndays
vars_days <- c(
  wbgt30 = "ndays_wbgt_shade_peak_ge30",
  wbgt32 = "ndays_wbgt_shade_peak_ge32",
  tw35   = "ndays_tw_peak_ge35"
)

# Chronic definition for Fig4
CHRONIC_DAYS <- 30

# ----------------------------
# BASEMAP
# ----------------------------
world_ll    <- ne_countries(scale = "medium", returnclass = "sf")
world_54030 <- st_transform(world_ll, 54030)  # ESRI:54030

# ============================================================
# HELPERS
# ============================================================

read_var <- function(ncfile, var) {
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0, nrow(r) > 0, ncol(r) > 0)
  r
}

to_m180_180 <- function(r) {
  if (xmax(r) > 180) rotate(r) else r
}

# time -> Date + year (+ month + ndays for weighted climatology)
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

# simple mean over time -> now reads precomputed MEAN_*.nc (timmean in CDO)
clim_mean <- function(ncfile, var) {
  read_mean_clim(ncfile, var) |> to_m180_180()
}

# weighted climatological mean -> now reads precomputed MEAN_*.nc
# (days-of-month weighting dropped: <0.1 C, and cancels in the change row)
clim_weighted_monthly_mean <- function(ncfile, var) {
  read_mean_clim(ncfile, var) |> to_m180_180()
}

# mean annual exceedance days:
# - input is monthly sums (cell_methods "time: sum") over months
# - we sum months per year then average across years -> 1 layer
# mean annual exceedance days -> now reads precomputed DAYS_*.nc
# (timmean -yearsum done in CDO)
mean_annual_days <- function(ncfile, var) {
  read_days_clim(ncfile, var) |> to_m180_180()
}

# Align src -> tgt
align_to <- function(src, tgt, method) {
  src <- to_m180_180(src)
  tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap="out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method = method)
}

# Robust symmetric limit for bias
sym_lim_q <- function(r, q = 0.995) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) == 0) return(1)
  as.numeric(quantile(abs(v), probs = q, na.rm = TRUE))
}

to_df <- function(r, nm="value") {
  as.data.frame(r, xy=TRUE, na.rm=TRUE) |> setNames(c("x","y",nm))
}

# ============================================================
# COMMON PLOT STYLE (Fig2 style: legends separated)
# ============================================================

style_theme <- function(show_legend=FALSE){
  theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=16, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=13, face="bold"),
      legend.text  = element_text(size=12),
      legend.key.height = unit(4.2, "mm"),
      legend.key.width  = unit(42, "mm"),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

scale_bias <- function(lim, legend_title){
  scale_fill_gradient2(
    low="#2166ac", mid="white", high="#b2182b",
    midpoint=0, limits=c(-lim, lim), oob=squish, name=legend_title
  )
}

scale_range <- function(vmax, legend_title){
  scale_fill_viridis_c(option="cividis", limits=c(0, vmax), oob=squish, name=legend_title)
}

extract_legend <- function(p){
  cowplot::get_legend(p + theme(plot.title = element_blank()))
}

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

# For Fig2/Fig3 (binned cold->hot)
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

# ============================================================
# PART A — FIG1 + SI1 + SI2 (validation + range)
# ============================================================

# 1) Compute present climatologies (1 layer)
ens_mean <- setNames(lapply(vars_mean, \(v) clim_mean(f_ens, v)), vars_mean)
era_mean <- setNames(lapply(vars_mean, \(v) clim_mean(f_era, v)), vars_mean)

ens_days <- setNames(lapply(vars_days, \(v) mean_annual_days(f_ens, v)), names(vars_days))
era_days <- setNames(lapply(vars_days, \(v) mean_annual_days(f_era, v)), names(vars_days))

ensmax_days <- setNames(lapply(vars_days, \(v) mean_annual_days(f_max, v)), names(vars_days))
ensmin_days <- setNames(lapply(vars_days, \(v) mean_annual_days(f_min, v)), names(vars_days))

# 2) Align ERA5 to ENS grid
era_mean_a <- mapply(align_to, era_mean, ens_mean, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)
era_days_a <- mapply(align_to, era_days, ens_days, MoreArgs=list(method="near"),     SIMPLIFY=FALSE)

# 3) Bias + range
bias_mean  <- mapply(`-`, ens_mean, era_mean_a, SIMPLIFY=FALSE)
bias_days  <- mapply(`-`, ens_days, era_days_a, SIMPLIFY=FALSE)
range_days <- mapply(`-`, ensmax_days, ensmin_days, SIMPLIFY=FALSE)

# ----------------------------
# FIG1 — Bias on exceedance days (3 panels + legend row)
# ----------------------------
lim1 <- max(sym_lim_q(bias_days$wbgt30), sym_lim_q(bias_days$wbgt32), sym_lim_q(bias_days$tw35))
leg_title1 <- "ENS − ERA5 (days/yr)"

p1 <- make_bias_map(bias_days$wbgt30, "WBGT ≥ 30°C", lim1, leg_title1, show_legend=FALSE)
p2 <- make_bias_map(bias_days$wbgt32, "WBGT ≥ 32°C", lim1, leg_title1, show_legend=FALSE)
p3 <- make_bias_map(bias_days$tw35,   "Tw ≥ 35°C",   lim1, leg_title1, show_legend=FALSE)

row_maps1 <- (p1 | p2 | p3)
leg1 <- extract_legend(make_bias_map(bias_days$wbgt30, NULL, lim1, leg_title1, show_legend=TRUE))
fig1 <- row_maps1 / wrap_elements(full = leg1) + plot_layout(heights = c(1, 0.12))

ggsave(file.path(dir_fig, "Fig1_bias_thresholds.png"), fig1, width=12, height=4.6, dpi=300)

# ----------------------------
# SI1 — Bias on MEAN indices (2×2 + legend row)
# ----------------------------
limS1 <- max(sym_lim_q(bias_mean$wbgt_shade_mean),
             sym_lim_q(bias_mean$tw_mean),
             sym_lim_q(bias_mean$hi_mean),
             sym_lim_q(bias_mean$humidex_mean))
leg_titleS1 <- "ENS − ERA5"

s1 <- make_bias_map(bias_mean$wbgt_shade_mean, nice_mean[["wbgt_shade_mean"]], limS1, leg_titleS1, FALSE)
s2 <- make_bias_map(bias_mean$tw_mean,         nice_mean[["tw_mean"]],         limS1, leg_titleS1, FALSE)
s3 <- make_bias_map(bias_mean$hi_mean,         nice_mean[["hi_mean"]],         limS1, leg_titleS1, FALSE)
s4 <- make_bias_map(bias_mean$humidex_mean,    nice_mean[["humidex_mean"]],    limS1, leg_titleS1, FALSE)

mapsS1 <- (s1 | s2) / (s3 | s4)
legS1  <- extract_legend(make_bias_map(bias_mean$wbgt_shade_mean, NULL, limS1, leg_titleS1, TRUE))

si1 <- cowplot::plot_grid(
  mapsS1, legS1,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggsave(file.path(dir_fig, "SI1_bias_means_4panel.png"), si1, width=10.5, height=7.6, dpi=300)

# ----------------------------
# SI2 — Ensemble range (max−min) on exceedance days (3 panels + legend row)
# ----------------------------
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

# ============================================================
# PART B — FIG2 (MEAN climatologies, 2×4)
# Present vs end-century (SSP5-8.5)
# ============================================================

# Weighted climatological means (monthly means -> 1 layer)
pres_mean <- setNames(lapply(vars_mean, \(v) clim_weighted_monthly_mean(f_pres, v)), vars_mean)
fut_mean  <- setNames(lapply(vars_mean, \(v) clim_weighted_monthly_mean(f_fut,  v)), vars_mean)

# Align future to present grid
fut_mean_a <- mapply(align_to, fut_mean, pres_mean, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)

# Limits across both periods (robust 0.5–99.5%)
lims_idx <- lapply(vars_mean, function(v){
  a <- values(pres_mean[[v]], mat=FALSE); a <- a[is.finite(a)]
  b <- values(fut_mean_a[[v]], mat=FALSE); b <- b[is.finite(b)]
  vv <- c(a,b)
  qs <- quantile(vv, probs=c(0.005, 0.995), na.rm=TRUE)
  as.numeric(qs)
})
names(lims_idx) <- vars_mean

breaks_idx <- lapply(vars_mean, function(v) breaks_pretty(lims_idx[[v]], n=4))
names(breaks_idx) <- vars_mean

unit_lab_mean <- function(v) if (v == "humidex_mean") "" else " (°C)"

maps_pres <- lapply(vars_mean, function(v){
  make_steps_map(
    pres_mean[[v]],
    title = nice_mean[[v]],
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_mean[[v]], unit_lab_mean(v)),
    show_legend = FALSE
  )
})

maps_fut <- lapply(vars_mean, function(v){
  make_steps_map(
    fut_mean_a[[v]],
    title = nice_mean[[v]],
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_mean[[v]], unit_lab_mean(v)),
    show_legend = FALSE
  )
})

legends <- lapply(vars_mean, function(v){
  p_leg <- make_steps_map(
    fut_mean_a[[v]],
    title = NULL,
    lims = lims_idx[[v]],
    breaks = breaks_idx[[v]],
    legend_title = paste0(nice_mean[[v]], unit_lab_mean(v)),
    show_legend = TRUE
  ) + theme(plot.title = element_blank())
  cowplot::get_legend(p_leg)
})

legend_row <- wrap_elements(full = legends[[1]]) |
              wrap_elements(full = legends[[2]]) |
              wrap_elements(full = legends[[3]]) |
              wrap_elements(full = legends[[4]])

# --- Projected change row (end-century minus present) ---
# Sequential, positive-only scale: warming is positive everywhere, so a diverging
# scale would waste half the colourbar. 0 -> per-index 99th percentile.
seq_warm <- c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929",
              "#ec7014","#cc4c02","#993404","#662506")
q_hi_chg <- function(r, q = 0.99){
  v <- values(r, mat = FALSE); v <- v[is.finite(v)]
  if (!length(v)) return(1); as.numeric(quantile(v, q, na.rm = TRUE))
}
make_change_map <- function(r, title, hi, breaks, legend_title, robust_file=NULL,
                            show_legend=FALSE, world=world_54030){
  rr <- project(r, "ESRI:54030")
  df <- to_df(rr, "v"); resx <- res(rr)[1]; resy <- res(rr)[2]
  p <- ggplot() +
    geom_tile(data=df, aes(x, y, fill=v), width=resx, height=resy)
  # stipple where models disagree on the sign of change (robust == 0);
  # no per-panel legend (a single shared stippling legend is built below)
  if (!is.null(robust_file) && file.exists(robust_file)) {
    rb  <- project(rast(robust_file)[[1]], "ESRI:54030", method="near")
    pts <- to_df(rb, "rob"); pts <- pts[pts$rob == 0, c("x","y")]
    if (nrow(pts) > 8000) pts <- pts[round(seq(1, nrow(pts), length.out = 8000)), ]
    if (nrow(pts)) p <- p + geom_point(data=pts, aes(x, y), size=0.05, alpha=0.5,
                                       colour="grey15", show.legend=FALSE)
  }
  p +
    geom_sf(data=world, fill=NA, color="grey20", linewidth=0.2) +
    coord_sf(crs=54030, expand=FALSE) +
    scale_fill_stepsn(colors=seq_warm, limits=c(0, hi), breaks=breaks, oob=squish,
                      name=legend_title, guide=guide_steps()) +
    labs(title=title) +
    theme_void() +
    theme(
      plot.title  = element_text(hjust=0.5, face="bold", size=16, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=13, face="bold"),
      legend.text  = element_text(size=12),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0)
    )
}

chg_mean   <- setNames(lapply(vars_mean, function(v) fut_mean_a[[v]] - pres_mean[[v]]), vars_mean)
hi_chg     <- setNames(lapply(vars_mean, function(v) ceiling(q_hi_chg(chg_mean[[v]], 0.99))), vars_mean)
breaks_chg <- setNames(lapply(vars_mean, function(v) breaks_pretty(c(0, hi_chg[[v]]), n=4)), vars_mean)

# inter-model agreement masks (from agreement_cdo.sh, MODE=permodel, SSP585).
# If absent, maps are drawn without stippling (graceful).
dir_agr  <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_SI"
robust_f <- function(v) file.path(dir_agr, sprintf("robust_%s_SSP585.nc", v))

maps_chg <- lapply(vars_mean, function(v){
  make_change_map(chg_mean[[v]], nice_mean[[v]], hi_chg[[v]], breaks_chg[[v]],
                  paste0("\u0394 ", nice_mean[[v]], unit_lab_mean(v)),
                  robust_file = robust_f(v), show_legend = FALSE)
})

# 4 horizontal binned legends (one per index), placed below the change row
legends_chg <- lapply(vars_mean, function(v){
  p_leg <- make_change_map(chg_mean[[v]], NULL, hi_chg[[v]], breaks_chg[[v]],
                           paste0("\u0394 ", nice_mean[[v]], unit_lab_mean(v)),
                           show_legend = TRUE) + theme(plot.title = element_blank())
  cowplot::get_legend(p_leg)
})
legend_row_chg <- wrap_elements(full = legends_chg[[1]]) |
                  wrap_elements(full = legends_chg[[2]]) |
                  wrap_elements(full = legends_chg[[3]]) |
                  wrap_elements(full = legends_chg[[4]])

# shared stippling legend entry — only shown if any disagreement exists
any_disagree <- any(vapply(vars_mean, function(v){
  f <- robust_f(v)
  if (!file.exists(f)) return(FALSE)
  any(values(rast(f)[[1]], mat = FALSE) == 0, na.rm = TRUE)
}, logical(1)))

.stipple_lab <- "models disagree on sign of change (< 4/5)"
stipple_leg <- cowplot::get_legend(
  ggplot(data.frame(x = 0, y = 0, grp = .stipple_lab)) +
    geom_point(aes(x, y, colour = grp), size = 1.6, alpha = 0.7) +
    scale_colour_manual(name = NULL, values = setNames("grey15", .stipple_lab)) +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
)

row1 <- (maps_pres[[1]] | maps_pres[[2]] | maps_pres[[3]] | maps_pres[[4]])
row2 <- (maps_fut[[1]]  | maps_fut[[2]]  | maps_fut[[3]]  | maps_fut[[4]])
row3 <- (maps_chg[[1]]  | maps_chg[[2]]  | maps_chg[[3]]  | maps_chg[[4]])

base_stack <- row_header("Present (1991–2020)") /
              row1 /
              row_header("SSP5-8.5 (2071–2100)") /
              row2 /
              legend_row /
              row_header("Projected change (2071\u20132100 \u2212 1991\u20132020)") /
              row3 /
              legend_row_chg

if (any_disagree) {
  fig2 <- base_stack / wrap_elements(full = stipple_leg) +
          plot_layout(heights = c(0.10, 1, 0.10, 1, 0.22, 0.14, 1, 0.22, 0.08))
} else {
  fig2 <- base_stack +
          plot_layout(heights = c(0.10, 1, 0.10, 1, 0.22, 0.14, 1, 0.22))
}

ggsave(file.path(dir_fig, "Fig2_means_2x4_SSP585.png"),
       fig2, width=16, height=14.5, dpi=300)

# ============================================================
# PART C — FIG3 (PEAK climatologies: mean of daily peaks, 2×2)
# Present vs end-century (SSP5-8.5)
# ============================================================

pres_peak <- setNames(lapply(vars_peak, \(v) clim_weighted_monthly_mean(f_pres, v)), vars_peak)
fut_peak  <- setNames(lapply(vars_peak, \(v) clim_weighted_monthly_mean(f_fut,  v)), vars_peak)

fut_peak_a <- mapply(align_to, fut_peak, pres_peak, MoreArgs=list(method="bilinear"), SIMPLIFY=FALSE)

lims_peak <- lapply(vars_peak, function(v){
  a <- values(pres_peak[[v]], mat=FALSE); a <- a[is.finite(a)]
  b <- values(fut_peak_a[[v]], mat=FALSE); b <- b[is.finite(b)]
  vv <- c(a,b)
  qs <- quantile(vv, probs=c(0.005, 0.995), na.rm=TRUE)
  as.numeric(qs)
})
names(lims_peak) <- vars_peak

breaks_peak <- lapply(vars_peak, function(v) breaks_pretty(lims_peak[[v]], n=4))
names(breaks_peak) <- vars_peak

unit_lab_peak <- function(v) " (°C)"

maps_pres_p <- lapply(vars_peak, function(v){
  make_steps_map(
    pres_peak[[v]],
    title = nice_peak[[v]],
    lims = lims_peak[[v]],
    breaks = breaks_peak[[v]],
    legend_title = paste0(nice_peak[[v]], unit_lab_peak(v)),
    show_legend = FALSE
  )
})

maps_fut_p <- lapply(vars_peak, function(v){
  make_steps_map(
    fut_peak_a[[v]],
    title = nice_peak[[v]],
    lims = lims_peak[[v]],
    breaks = breaks_peak[[v]],
    legend_title = paste0(nice_peak[[v]], unit_lab_peak(v)),
    show_legend = FALSE
  )
})

legends_p <- lapply(vars_peak, function(v){
  p_leg <- make_steps_map(
    fut_peak_a[[v]],
    title = NULL,
    lims = lims_peak[[v]],
    breaks = breaks_peak[[v]],
    legend_title = paste0(nice_peak[[v]], unit_lab_peak(v)),
    show_legend = TRUE
  ) + theme(plot.title = element_blank())
  cowplot::get_legend(p_leg)
})

legend_row_p <- wrap_elements(full = legends_p[[1]]) |
                wrap_elements(full = legends_p[[2]])

row1p <- (maps_pres_p[[1]] | maps_pres_p[[2]])
row2p <- (maps_fut_p[[1]]  | maps_fut_p[[2]])

fig3 <- row_header("Present (1991–2020)") /
        row1p /
        row_header("SSP5-8.5 (2071–2100)") /
        row2p /
        legend_row_p +
        plot_layout(heights = c(0.08, 1, 0.08, 1, 0.18))

ggsave(file.path(dir_fig, "Fig3_peaks_2x2_SSP585.png"),
       fig3, width=10.5, height=8.8, dpi=300)

# ============================================================
# PART D — FIG4 (CHRONIC exceedance days + co-occurrence, 3×3)
# Present / 2050 / 2100 (SSP5-8.5)
# ============================================================
# --- Fig4 (SSP5-8.5 only)
f_pres_ens <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
f_2050_585 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2041-2070/HEATSTRESS_ensmean_SSP585_2041-2070.nc"
f_2100_585 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP585/2071-2100/HEATSTRESS_ensmean_SSP585_2071-2100.nc"

# --- Fig5 (co-occurrence, SSP126/245/585; 2050 & 2100)
f_2050_126 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2041-2070/HEATSTRESS_ensmean_SSP126_2041-2070.nc"
f_2100_126 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP126/2071-2100/HEATSTRESS_ensmean_SSP126_2071-2100.nc"

f_2050_245 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2041-2070/HEATSTRESS_ensmean_SSP245_2041-2070.nc"
f_2100_245 <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/SSP245/2071-2100/HEATSTRESS_ensmean_SSP245_2071-2100.nc"

f_2050_585 <- f_2050_585
f_2100_585 <- f_2100_585

# ----------------------------
# VARIABLES (V2 names)
# ----------------------------
v_tw35   <- "ndays_tw_peak_ge35"
v_wbgt32 <- "ndays_wbgt_shade_peak_ge32"

CHRONIC_DAYS <- 30

# ----------------------------
# BASEMAP
# ----------------------------
world_ll    <- ne_countries(scale="medium", returnclass="sf")
world_54030 <- st_transform(world_ll, 54030)  # ESRI:54030

# ============================================================
# HELPERS
# ============================================================

read_var <- function(ncfile, var) {
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0)
  r
}

to_m180_180 <- function(r) {
  if (xmax(r) > 180) rotate(r) else r
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

# monthly sums -> yearly sums -> mean across years (1 layer)
# monthly sums -> yearly sums -> mean across years -> now reads DAYS_*.nc
mean_annual_days <- function(ncfile, var) {
  read_days_clim(ncfile, var) |> to_m180_180()
}

align_to <- function(src, tgt, method="near") {
  src <- to_m180_180(src)
  tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap="out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method=method)
}

to_df <- function(r, nm="value") {
  as.data.frame(r, xy=TRUE, na.rm=TRUE) |> setNames(c("x","y",nm))
}

bin01 <- function(r, thr = CHRONIC_DAYS) {
  # returns 0 if <thr, 1 if >=thr
  classify(r, rcl = matrix(c(-Inf, thr, 0,
                            thr, Inf, 1), ncol=3, byrow=TRUE))
}

# ============================================================
# PLOTTING HELPERS (binary 0/1 red maps)
# ============================================================

cols_bin <- c("0"="grey92", "1"="#d7191c")

style_theme_bin <- function(show_legend=FALSE){
  theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=13, margin=margin(b=2)),
      plot.margin = margin(0,0,0,0),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size=12, face="bold"),
      legend.text  = element_text(size=11),
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

row_header <- function(txt){
  wrap_elements(full = textGrob(txt, gp=gpar(fontface="bold", fontsize=18))) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
}

extract_legend <- function(p){
  cowplot::get_legend(p + theme(plot.title = element_blank()))
}

# ============================================================
# FIG4 — Tw & WBGT chronic exposure (SSP5-8.5), 3 periods × 2 cols
# ============================================================

# 1) Compute mean annual days per period
tw_pres   <- mean_annual_days(f_pres_ens, v_tw35)
wbgt_pres <- mean_annual_days(f_pres_ens, v_wbgt32)

tw_2050   <- mean_annual_days(f_2050_585, v_tw35)
wbgt_2050 <- mean_annual_days(f_2050_585, v_wbgt32)

tw_2100   <- mean_annual_days(f_2100_585, v_tw35)
wbgt_2100 <- mean_annual_days(f_2100_585, v_wbgt32)

# 2) Align to present grid (near)
tw_2050_a   <- align_to(tw_2050,   tw_pres, method="near")
wbgt_2050_a <- align_to(wbgt_2050, tw_pres, method="near")

tw_2100_a   <- align_to(tw_2100,   tw_pres, method="near")
wbgt_2100_a <- align_to(wbgt_2100, tw_pres, method="near")

# 3) Binary chronic masks
twP  <- bin01(tw_pres)
wbP  <- bin01(wbgt_pres)

tw50 <- bin01(tw_2050_a)
wb50 <- bin01(wbgt_2050_a)

tw00 <- bin01(tw_2100_a)
wb00 <- bin01(wbgt_2100_a)

# 4) Build panels
p11 <- make_bin_map(twP,  "Tw ≥ 35°C")
p12 <- make_bin_map(wbP,  "WBGT ≥ 32°C")

p21 <- make_bin_map(tw50, "Tw ≥ 35°C")
p22 <- make_bin_map(wb50, "WBGT ≥ 32°C")

p31 <- make_bin_map(tw00, "Tw ≥ 35°C")
p32 <- make_bin_map(wb00, "WBGT ≥ 32°C")

row1 <- (p11 | p12)
row2 <- (p21 | p22)
row3 <- (p31 | p32)

# Legend (single)
leg4 <- extract_legend(make_bin_map(twP, NULL, show_legend=TRUE))
legend_row4 <- wrap_elements(full = leg4)

fig4 <- row_header("Present (1991–2020)") /
        row1 /
        row_header("2050 (2041–2070) — SSP5-8.5") /
        row2 /
        row_header("2100 (2071–2100) — SSP5-8.5") /
        row3 /
        legend_row4 +
        plot_layout(heights = c(0.06, 1, 0.06, 1, 0.06, 1, 0.12))

ggsave(file.path(dir_fig, "Fig4_chronic_Tw35_WBGT32_SSP585.png"),
       fig4, width=10.5, height=12.0, dpi=300)

message("DONE Fig4 -> ", file.path(dir_fig, "Fig4_chronic_Tw35_WBGT32_SSP585.png"))

# ============================================================
# FIG5 — Co-occurrence chronic exposure (Tw35 & WBGT32), 2 periods × 3 SSP
# Binary 0/1 maps (red = exposed)
# ============================================================

# Helper: compute co-occurrence binary map for a file
coocc_bin_from_file <- function(ncfile, var_tw=v_tw35, var_wbgt=v_wbgt32, tgt=NULL) {
  tw   <- mean_annual_days(ncfile, var_tw)
  wbgt <- mean_annual_days(ncfile, var_wbgt)

  if (!is.null(tgt)) {
    tw   <- align_to(tw,   tgt, method="near")
    wbgt <- align_to(wbgt, tgt, method="near")
  }

  twb <- bin01(tw) * bin01(wbgt)   # co-occurrence 1 only if both chronic
  twb
}

# Use present Tw grid as a common target for all maps (stable layout)
tgt_grid <- tw_pres

# 2050
co_2050_126 <- coocc_bin_from_file(f_2050_126, tgt=tgt_grid)
co_2050_245 <- coocc_bin_from_file(f_2050_245, tgt=tgt_grid)
co_2050_585 <- coocc_bin_from_file(f_2050_585, tgt=tgt_grid)

# 2100
co_2100_126 <- coocc_bin_from_file(f_2100_126, tgt=tgt_grid)
co_2100_245 <- coocc_bin_from_file(f_2100_245, tgt=tgt_grid)
co_2100_585 <- coocc_bin_from_file(f_2100_585, tgt=tgt_grid)

# Panels: 2 rows × 3 cols
q11 <- make_bin_map(co_2050_126, "SSP1-2.6", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)
q12 <- make_bin_map(co_2050_245, "SSP2-4.5", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)
q13 <- make_bin_map(co_2050_585, "SSP5-8.5", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)

q21 <- make_bin_map(co_2100_126, "SSP1-2.6", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)
q22 <- make_bin_map(co_2100_245, "SSP2-4.5", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)
q23 <- make_bin_map(co_2100_585, "SSP5-8.5", legend_title="Co-occurrence chronic exposure", show_legend=FALSE)

row2050 <- (q11 | q12 | q13)
row2100 <- (q21 | q22 | q23)

# Legend (single)
leg5 <- extract_legend(make_bin_map(co_2050_585, NULL, legend_title="Co-occurrence chronic exposure", show_legend=TRUE))
legend_row5 <- wrap_elements(full = leg5)

fig5 <- row_header("2050 (2041–2070)") /
        row2050 /
        row_header("2100 (2071–2100)") /
        row2100 /
        legend_row5 +
        plot_layout(heights = c(0.08, 1, 0.08, 1, 0.14))

ggsave(file.path(dir_fig, "Fig5_coocc_chronic_2050_2100_3SSP.png"),
       fig5, width=12.5, height=8.8, dpi=300)

message("DONE Fig5 -> ", file.path(dir_fig, "Fig5_coocc_chronic_2050_2100_3SSP.png"))
# ============================================================
# PART E — CHRONIC-THRESHOLD SENSITIVITY (reviewer R2 #1)
# Exposed land area vs the persistence threshold (15 / 30 / 45 days/yr),
# per hazard x scenario x period. One summary figure + CSV.
# Re-thresholds the precomputed DAYS_*.nc fields -> cheap.
# ============================================================
message("PART E — chronic-threshold sensitivity (15/30/45 d/yr)")

thr_set <- c(15, 30, 45)
haz     <- setNames(c(v_tw35, v_wbgt32), c("Tw \u2265 35\u00B0C", "WBGT \u2265 32\u00B0C"))

# scenario x period x file (Present is scenario-agnostic baseline)
scen_files <- list(
  list(scen="SSP1-2.6", period="2041\u20132070", f=f_2050_126),
  list(scen="SSP2-4.5", period="2041\u20132070", f=f_2050_245),
  list(scen="SSP5-8.5", period="2041\u20132070", f=f_2050_585),
  list(scen="SSP1-2.6", period="2071\u20132100", f=f_2100_126),
  list(scen="SSP2-4.5", period="2071\u20132100", f=f_2100_245),
  list(scen="SSP5-8.5", period="2071\u20132100", f=f_2100_585)
)

# per-cell area (km^2) on the reference grid
.ref     <- mean_annual_days(f_pres_ens, v_tw35)
area_km2 <- cellSize(.ref, unit="km")

exposed_Mkm2 <- function(days_r, thr){
  m <- days_r >= thr                       # 1 where chronic, 0 else, NA on ocean
  as.numeric(global(m * area_km2, "sum", na.rm=TRUE)) / 1e6
}

rows <- list()
for (h in seq_along(haz)){
  hv <- haz[[h]]; hname <- names(haz)[h]
  for (sf in scen_files){
    d <- mean_annual_days(sf$f, hv) |> align_to(.ref, method="near")
    for (thr in thr_set){
      rows[[length(rows)+1]] <- data.frame(
        hazard    = hname,
        scenario  = sf$scen,
        period    = sf$period,
        threshold = thr,
        area_Mkm2 = exposed_Mkm2(d, thr)
      )
    }
  }
}
sens <- do.call(rbind, rows)
write.csv(sens, file.path(dir_fig, "chronic_threshold_sensitivity.csv"), row.names=FALSE)

fig_sens <- ggplot(sens, aes(x=factor(threshold), y=area_Mkm2, fill=scenario)) +
  geom_col(position=position_dodge(width=0.78), width=0.72) +
  facet_grid(hazard ~ period, scales="free_y") +
  scale_fill_manual(values=c("SSP1-2.6"="#4575b4", "SSP2-4.5"="#fdae61", "SSP5-8.5"="#d73027")) +
  labs(x="Chronic persistence threshold (days yr\u207B\u00B9)",
       y="Exposed land area (million km\u00B2)",
       fill=NULL,
       title="Sensitivity of chronic humid-heat exposure to the persistence threshold") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        legend.position="bottom",
        strip.text=element_text(face="bold"),
        panel.grid.minor=element_blank())

ggsave(file.path(dir_fig, "Fig_threshold_sensitivity.png"),
       fig_sens, width=10, height=7, dpi=300)

message("DONE PART E -> Fig_threshold_sensitivity.png + chronic_threshold_sensitivity.csv")
