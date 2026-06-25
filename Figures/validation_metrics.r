# ============================================================
# validation_metrics.r
# Quantitative validation of heat-stress indices against ERA5 (1991-2020)
# Produces, for each index/threshold field:
#   - mean bias (ENS - ERA5)        [area-weighted, cos-lat]
#   - RMSE                          [area-weighted]
#   - MAE                           [area-weighted]
#   - spatial Pearson correlation   [area-weighted]
#   - N cells used
# Outputs:
#   ERA5_validation_metrics.csv
#   ERA5_validation_metrics.tex   (table body; plug into your tabular)
#
# Answers reviewer R2 #3 (replace qualitative ERA5 comparison by statistics).
# Uses the SAME conventions (paths, var names, helpers) as figures.r.
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
  library(lubridate)
})

# ----------------------------
# PATHS (EDIT)
# ----------------------------
dir_ens <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020"
dir_era <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ERA5/packs"
dir_out <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Validation"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

f_ens <- file.path(dir_ens, "HEATSTRESS_ensmean_pseudoHist_1991-2020.nc")
f_era <- file.path(dir_era, "HEATSTRESS_mon_ERA5_gr025_1991-2020.nc")

# Restrict statistics to land cells (recommended for human heat stress).
# Ocean WBGT/Tw values inflate correlations and are not human-relevant.
LAND_ONLY <- TRUE

# ----------------------------
# VARIABLES (consistent with figures.r)
# ----------------------------
# Mean indices (climatological means)
vars_mean  <- c(wbgt_shade_mean = "WBGT",
                tw_mean         = "Tw",
                hi_mean         = "Heat Index",
                humidex_mean    = "Humidex")
units_mean <- c(wbgt_shade_mean = "\u00B0C",
                tw_mean         = "\u00B0C",
                hi_mean         = "\u00B0C",
                humidex_mean    = "\u2013")   # Humidex is dimensionless

# Exceedance-day fields (mean annual days)
vars_days <- c(ndays_wbgt_shade_peak_ge30 = "WBGT \u2265 30\u00B0C",
               ndays_wbgt_shade_peak_ge32 = "WBGT \u2265 32\u00B0C",
               ndays_tw_peak_ge35         = "Tw \u2265 35\u00B0C")

# ============================================================
# HELPERS
# ============================================================
read_var <- function(ncfile, var){
  r <- rast(ncfile, subds = var)
  stopifnot(inherits(r, "SpatRaster"), nlyr(r) > 0)
  r
}
to_m180_180 <- function(r) if (xmax(r) > 180) rotate(r) else r

nc_time_ymd <- function(ncfile){
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  t <- ncvar_get(nc, "time")
  u <- ncatt_get(nc, "time", "units")$value
  origin <- sub(".*since\\s+", "", u)
  origin <- strsplit(origin, " ")[[1]][1]
  d <- as.Date(t, origin = origin)
  data.frame(idx = seq_along(d),
             year = as.integer(format(d, "%Y")),
             ndays = lubridate::days_in_month(d))
}

clim_weighted_monthly_mean <- function(ncfile, var){
  r <- read_var(ncfile, var) |> to_m180_180()
  td <- nc_time_ymd(ncfile); stopifnot(nlyr(r) == nrow(td))
  w <- td$ndays
  app(r, fun = function(...){
    v <- c(...); ok <- is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * w[ok]) / sum(w[ok])
  })
}

mean_annual_days <- function(ncfile, var){
  r  <- read_var(ncfile, var) |> to_m180_180()
  td <- nc_time_ymd(ncfile); stopifnot(nlyr(r) == nrow(td))
  yrs <- unique(td$year)
  annual <- lapply(yrs, function(y){
    idx <- td$idx[td$year == y]
    sum(r[[idx]], na.rm = TRUE)
  })
  mean(rast(annual), na.rm = TRUE)
}

align_to <- function(src, tgt, method){
  src <- to_m180_180(src); tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap = "out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method = method)
}

# weighted statistics
wmean <- function(x, w) sum(x * w) / sum(w)
wcor  <- function(x, y, w){
  mx <- wmean(x, w); my <- wmean(y, w)
  covxy <- wmean((x - mx) * (y - my), w)
  vx <- wmean((x - mx)^2, w); vy <- wmean((y - my)^2, w)
  covxy / sqrt(vx * vy)
}

# land mask (1 on land, NA elsewhere) on a target grid
land_mask_on <- function(tgt){
  suppressPackageStartupMessages({ library(rnaturalearth); library(sf) })
  land <- ne_countries(scale = "medium", returnclass = "sf") |> st_make_valid()
  rasterize(vect(land), tgt, field = 1, background = NA, touches = TRUE)
}

# core comparison: ENS field vs ERA5 field (ERA5 aligned to ENS grid)
metrics_field <- function(ens_r, era_r, lmask = NULL, method = "bilinear"){
  era_a <- align_to(era_r, ens_r, method = method)
  wlat  <- init(ens_r, "y"); wlat <- cos(wlat * pi / 180)   # area weight

  ve <- values(ens_r, mat = FALSE)
  va <- values(era_a, mat = FALSE)
  vw <- values(wlat,  mat = FALSE)

  ok <- is.finite(ve) & is.finite(va) & is.finite(vw)
  if (!is.null(lmask)){
    vm <- values(lmask, mat = FALSE)
    ok <- ok & is.finite(vm)
  }
  ve <- ve[ok]; va <- va[ok]; vw <- vw[ok]
  d  <- ve - va

  data.frame(
    bias = wmean(d, vw),
    rmse = sqrt(wmean(d^2, vw)),
    mae  = wmean(abs(d), vw),
    r    = wcor(ve, va, vw),
    n    = length(ve)
  )
}

# ============================================================
# RUN
# ============================================================
# build land mask once on the ENS climatology grid
lmask <- NULL
if (LAND_ONLY){
  g0 <- clim_weighted_monthly_mean(f_ens, names(vars_mean)[1])
  lmask <- land_mask_on(g0)
}

rows <- list()

# mean indices (bilinear regrid for continuous fields)
for (v in names(vars_mean)){
  message("Validating mean index: ", vars_mean[[v]])
  ens <- clim_weighted_monthly_mean(f_ens, v)
  era <- clim_weighted_monthly_mean(f_era, v)
  m <- metrics_field(ens, era, lmask, method = "bilinear")
  m$variable <- vars_mean[[v]]; m$unit <- units_mean[[v]]
  rows[[length(rows) + 1]] <- m
}

# exceedance-day fields (nearest regrid for count-like fields)
for (v in names(vars_days)){
  message("Validating exceedance field: ", vars_days[[v]])
  ens <- mean_annual_days(f_ens, v)
  era <- mean_annual_days(f_era, v)
  m <- metrics_field(ens, era, lmask, method = "near")
  m$variable <- vars_days[[v]]; m$unit <- "days/yr"
  rows[[length(rows) + 1]] <- m
}

tab <- do.call(rbind, rows)
tab <- tab[, c("variable", "unit", "bias", "rmse", "mae", "r", "n")]
tab$bias <- round(tab$bias, 2)
tab$rmse <- round(tab$rmse, 2)
tab$mae  <- round(tab$mae, 2)
tab$r    <- round(tab$r, 3)

# CSV (human-readable, UTF-8)
write.csv(tab, file.path(dir_out, "ERA5_validation_metrics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# LaTeX body (escape special glyphs)
esc_tex <- function(s){
  s <- gsub("\u00B0", "$^{\\circ}$", s, fixed = TRUE)  # degree
  s <- gsub("\u2265", "$\\geq$",     s, fixed = TRUE)   # >=
  s <- gsub("\u2013", "--",          s, fixed = TRUE)   # en dash
  s
}
tex <- apply(tab, 1, function(z) paste0(
  esc_tex(z["variable"]), " & ", esc_tex(z["unit"]), " & ",
  z["bias"], " & ", z["rmse"], " & ", z["mae"], " & ", z["r"], " \\\\"
))
writeLines(tex, file.path(dir_out, "ERA5_validation_metrics.tex"))

message("DONE -> ", file.path(dir_out, "ERA5_validation_metrics.csv"))
print(tab)

# Suggested table header for the manuscript:
#   Variable & Unit & Bias & RMSE & MAE & Spatial r \\
# (statistics area-weighted by cos(latitude); land cells only if LAND_ONLY=TRUE)
