# ============================================================
# intermodel_agreement.r
# Inter-model agreement on the projected change of heat-stress indices.
# AR6-style map: ensemble-mean change (colour) + stippling where the
# models do NOT agree on the sign of change (low robustness).
#
# Output (in dir_fig): SI_agreement_<SCENARIO>_<PERIOD>.png
# Answers reviewer R2 minor #2 (inter-model agreement maps).
#
# IMPORTANT: this REQUIRES per-model files. The ensemble mean cannot give
# agreement (the spread is already collapsed). A conservative fallback that
# uses only ensmin/ensmax is provided at the very bottom (USE_FALLBACK).
# ============================================================

suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(ggplot2); library(patchwork)
  library(scales); library(sf); library(rnaturalearth); library(rnaturalearthdata)
  library(cowplot); library(grid); library(lubridate)
})

# ----------------------------
# PATHS / SETTINGS (EDIT)
# ----------------------------
dir_fig <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_SI"
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# The 5 CMIP6 models — EDIT to match your filenames exactly.
MODELS <- c("MODEL1", "MODEL2", "MODEL3", "MODEL4", "MODEL5")

# Per-model file builders — EDIT to match YOUR per-model layout.
# These must point to the per-model monthly files (one per model), with the
# same internal variables as the ensemble files.
f_model_present <- function(model){
  sprintf("C:/Data/HEAT_STRESS_MONTHLY_V2/PERMODEL/pseudoHist/1991-2020/%s/HEATSTRESS_%s_pseudoHist_1991-2020.nc",
          model, model)
}
f_model_future <- function(model, scenario, period){
  sprintf("C:/Data/HEAT_STRESS_MONTHLY_V2/PERMODEL/%s/%s/%s/HEATSTRESS_%s_%s_%s.nc",
          scenario, period, model, model, scenario, period)
}

# Headline change to map (edit freely)
SCENARIO   <- "SSP585"
PERIOD     <- "2071-2100"
SCEN_LABEL <- "SSP5-8.5 (2071\u20132100)"

# Indices to assess (mean climatologies)
vars <- c(wbgt_shade_mean = "WBGT", tw_mean = "Tw")

# Agreement rule
N_MODELS      <- length(MODELS)
K_AGREE       <- 4              # robust if >= K_AGREE of N agree on the sign (AR6-like ~80%)
STIPPLE_WHERE <- "disagree"     # "disagree" -> dots where NOT robust (AR6 simple method)
                                # "agree"    -> dots where robust
MAX_STIPPLE   <- 8000           # cap number of stipple points (readability)

# ============================================================
# HELPERS
# ============================================================
read_var <- function(ncfile, var){
  r <- rast(ncfile, subds = var); stopifnot(nlyr(r) > 0); r
}
to_m180_180 <- function(r) if (xmax(r) > 180) rotate(r) else r

nc_time_ymd <- function(ncfile){
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  t <- ncvar_get(nc, "time"); u <- ncatt_get(nc, "time", "units")$value
  origin <- sub(".*since\\s+", "", u); origin <- strsplit(origin, " ")[[1]][1]
  d <- as.Date(t, origin = origin)
  data.frame(idx = seq_along(d), ndays = lubridate::days_in_month(d))
}

clim_weighted_monthly_mean <- function(ncfile, var){
  r <- read_var(ncfile, var) |> to_m180_180()
  td <- nc_time_ymd(ncfile); stopifnot(nlyr(r) == nrow(td)); w <- td$ndays
  app(r, fun = function(...){
    v <- c(...); ok <- is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * w[ok]) / sum(w[ok])
  })
}

align_to <- function(src, tgt, method = "bilinear"){
  src <- to_m180_180(src); tgt <- to_m180_180(tgt)
  src <- crop(src, ext(tgt), snap = "out")
  if (!same.crs(src, tgt)) src <- project(src, tgt)
  resample(src, tgt, method = method)
}

to_df <- function(r, nm = "value") as.data.frame(r, xy = TRUE, na.rm = TRUE) |> setNames(c("x", "y", nm))

sym_lim_q <- function(r, q = 0.995){
  v <- values(r, mat = FALSE); v <- v[is.finite(v)]
  if (!length(v)) return(1)
  as.numeric(quantile(abs(v), q, na.rm = TRUE))
}

world_ll    <- ne_countries(scale = "medium", returnclass = "sf")
world_54030 <- st_transform(world_ll, 54030)

# ============================================================
# CORE: per-model delta -> ensemble-mean change + agreement mask
# ============================================================
agreement_field <- function(var){
  # present per model, align all to model-1 present grid
  pres <- lapply(MODELS, function(m) clim_weighted_monthly_mean(f_model_present(m), var))
  tgt  <- pres[[1]]
  pres <- lapply(seq_along(MODELS), function(i)
            if (i == 1) pres[[i]] else align_to(pres[[i]], tgt))
  fut  <- lapply(MODELS, function(m)
            align_to(clim_weighted_monthly_mean(f_model_future(m, SCENARIO, PERIOD), var), tgt))

  deltas <- lapply(seq_along(MODELS), function(i) fut[[i]] - pres[[i]])
  dstack <- rast(deltas)

  mean_delta <- mean(dstack, na.rm = TRUE)
  nvalid <- sum(!is.na(dstack))                 # models available per cell
  npos   <- sum(dstack > 0, na.rm = TRUE)       # models with positive change
  nneg   <- nvalid - npos
  robust <- (npos >= K_AGREE) | (nneg >= K_AGREE)
  robust <- robust * 1                          # logical -> 0/1 numeric

  list(mean_delta = mean_delta, robust = robust)
}

# ============================================================
# PLOT
# ============================================================
make_agreement_map <- function(var, title){
  af <- agreement_field(var)
  md <- project(af$mean_delta, "ESRI:54030")
  rb <- project(af$robust,     "ESRI:54030", method = "near")
  lim <- sym_lim_q(af$mean_delta)

  df   <- to_df(md, "delta")
  resx <- res(md)[1]; resy <- res(md)[2]

  # stipple points
  mask_val <- if (STIPPLE_WHERE == "disagree") 0 else 1
  pts <- to_df(rb, "rob")
  pts <- pts[pts$rob == mask_val, c("x", "y")]
  if (nrow(pts) > MAX_STIPPLE) pts <- pts[round(seq(1, nrow(pts), length.out = MAX_STIPPLE)), ]

  ggplot() +
    geom_tile(data = df, aes(x, y, fill = delta), width = resx, height = resy) +
    geom_point(data = pts, aes(x, y), size = 0.05, alpha = 0.5, colour = "grey15") +
    geom_sf(data = world_54030, fill = NA, colour = "grey20", linewidth = 0.2) +
    coord_sf(crs = 54030, expand = FALSE) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim), oob = squish,
                         name = paste0("\u0394 ", title, " (\u00B0C)")) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 2)),
      legend.position = "bottom",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text  = element_text(size = 12),
      legend.key.height = unit(4.2, "mm"),
      legend.key.width  = unit(42, "mm")
    )
}

p1 <- make_agreement_map("wbgt_shade_mean", "WBGT")
p2 <- make_agreement_map("tw_mean", "Tw")

cap <- if (STIPPLE_WHERE == "disagree")
  sprintf("Stippling: grid cells where fewer than %d/%d models agree on the sign of change.", K_AGREE, N_MODELS) else
  sprintf("Stippling: grid cells where at least %d/%d models agree on the sign of change.", K_AGREE, N_MODELS)

fig <- (p1 | p2) +
  plot_annotation(
    title   = paste0("Inter-model agreement on projected change \u2014 ", SCEN_LABEL),
    caption = cap,
    theme   = theme(plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
                    plot.caption = element_text(hjust = 0.5, size = 11))
  )

ggsave(file.path(dir_fig, sprintf("SI_agreement_%s_%s.png", SCENARIO, PERIOD)),
       fig, width = 14, height = 6.4, dpi = 300)

message("DONE -> ", file.path(dir_fig, sprintf("SI_agreement_%s_%s.png", SCENARIO, PERIOD)))


# ============================================================
# OPTIONAL FALLBACK — conservative robustness from ensmin/ensmax ONLY
# (no per-model files needed)
# ------------------------------------------------------------
# Rationale: ensmin/ensmax bracket all members at each cell. Therefore:
#   - if ensmin(future) - ensmax(present) > 0  => ALL models increase (robust+)
#   - if ensmax(future) - ensmin(present) < 0  => ALL models decrease (robust-)
# This is a SUFFICIENT (conservative) condition: it UNDER-estimates agreement
# and cannot express "4/5". Use only if you cannot access the members.
# Set USE_FALLBACK <- TRUE and fill the four paths per index.
# ============================================================
USE_FALLBACK <- FALSE

if (USE_FALLBACK){
  f_pres_ensmin <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmin_pseudoHist_1991-2020.nc"
  f_pres_ensmax <- "C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmax_pseudoHist_1991-2020.nc"
  f_fut_ensmin  <- sprintf("C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/%s/%s/HEATSTRESS_ensmin_%s_%s.nc",
                           SCENARIO, PERIOD, SCENARIO, PERIOD)
  f_fut_ensmax  <- sprintf("C:/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE/%s/%s/HEATSTRESS_ensmax_%s_%s.nc",
                           SCENARIO, PERIOD, SCENARIO, PERIOD)

  robust_minmax <- function(var){
    pmin <- clim_weighted_monthly_mean(f_pres_ensmin, var)
    pmax <- clim_weighted_monthly_mean(f_pres_ensmax, var)
    fmin <- align_to(clim_weighted_monthly_mean(f_fut_ensmin, var), pmin)
    fmax <- align_to(clim_weighted_monthly_mean(f_fut_ensmax, var), pmin)
    pmax <- align_to(pmax, pmin)
    inc <- (fmin - pmax) > 0     # guaranteed all-models increase
    dec <- (fmax - pmin) < 0     # guaranteed all-models decrease
    (inc | dec) * 1              # 1 = unanimous sign (lower bound on agreement)
  }

  make_fallback_map <- function(var, title){
    rb <- project(robust_minmax(var), "ESRI:54030", method = "near")
    df <- to_df(rb, "rob")
    resx <- res(rb)[1]; resy <- res(rb)[2]
    ggplot() +
      geom_tile(data = df, aes(x, y, fill = factor(rob)), width = resx, height = resy) +
      geom_sf(data = world_54030, fill = NA, colour = "grey20", linewidth = 0.2) +
      coord_sf(crs = 54030, expand = FALSE) +
      scale_fill_manual(values = c("0" = "grey92", "1" = "#1a9850"),
                        breaks = c("1", "0"),
                        labels = c("Unanimous sign", "Not unanimous"),
                        name = "Conservative agreement") +
      labs(title = title) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.position = "bottom",
            legend.title = element_text(size = 13, face = "bold"),
            legend.text  = element_text(size = 12))
  }

  fb <- (make_fallback_map("wbgt_shade_mean", "WBGT") |
         make_fallback_map("tw_mean", "Tw")) +
    plot_annotation(
      title = paste0("Conservative inter-model agreement (min/max bound) \u2014 ", SCEN_LABEL),
      caption = "Green = all 5 models share the sign of change (lower bound; underestimates agreement).",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                    plot.caption = element_text(hjust = 0.5, size = 11)))

  ggsave(file.path(dir_fig, sprintf("SI_agreement_FALLBACK_%s_%s.png", SCENARIO, PERIOD)),
         fb, width = 14, height = 6.4, dpi = 300)
  message("DONE (fallback) -> conservative agreement map")
}
