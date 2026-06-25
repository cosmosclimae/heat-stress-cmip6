# ============================================================
# make_maps_from_cdo.r  — LIGHT R: only draws maps from CDO outputs.
#   validation_cdo.sh -> bias_<VAR>.nc
#   agreement_cdo.sh  -> meandelta_<VAR>_<SCEN>.nc, robust_<VAR>_<SCEN>.nc
# ============================================================

suppressPackageStartupMessages({
  library(terra); library(ggplot2); library(patchwork); library(scales)
  library(sf); library(rnaturalearth); library(rnaturalearthdata); library(grid)
})

# ---------------- PATHS / SETTINGS (EDIT) ----------------
dir_val <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Validation"
dir_agr <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_SI"
dir_fig <- "C:/Data/HEAT_STRESS_MONTHLY_V2/Figures_SI"
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

DO_BIAS_MAP  <- TRUE
DO_AGREEMENT <- TRUE
PERIOD_LABEL <- "2071\u20132100"
MAX_STIPPLE  <- 8000
HI_OVERRIDE  <- NA      # set e.g. 6 to force the shared upper limit; NA = auto (q99)

# scenarios to show (file tag -> nice label). Must match agreement_cdo.sh outputs.
scenarios <- c(SSP245 = "SSP2-4.5", SSP585 = "SSP5-8.5")
# indices (file tag -> nice label)
vars_idx  <- c(wbgt_shade_mean = "WBGT", tw_mean = "Tw")
# bias fields to map (validation)
vars_bias <- c(wbgt_shade_mean = "WBGT", tw_mean = "Tw",
               hi_mean = "Heat Index", humidex_mean = "Humidex")

# ---------------- HELPERS ----------------
world_54030 <- st_transform(ne_countries(scale = "medium", returnclass = "sf"), 54030)
to_df  <- function(r, nm = "value") as.data.frame(r, xy = TRUE, na.rm = TRUE) |> setNames(c("x","y",nm))
read1  <- function(f) rast(f)[[1]]
sym_lim_q <- function(r, q = 0.995){ v <- values(r, mat=FALSE); v <- v[is.finite(v)]
  if (!length(v)) return(1); as.numeric(quantile(abs(v), q, na.rm=TRUE)) }
q_hi <- function(r, q = 0.99){ v <- values(r, mat=FALSE); v <- v[is.finite(v)]
  if (!length(v)) return(1); as.numeric(quantile(v, q, na.rm=TRUE)) }

base_theme <- function() theme_void() +
  theme(plot.title   = element_text(hjust=0.5, face="bold", size=15, margin=margin(b=2)),
        legend.position = "bottom",
        legend.title = element_text(size=13, face="bold"),
        legend.text  = element_text(size=12),
        legend.key.height = unit(4.2,"mm"), legend.key.width = unit(46,"mm"))

# warm sequential palette (white -> yellow -> red -> dark red)
warm_cols <- c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929",
               "#ec7014","#cc4c02","#993404","#662506")

# ---------------- AGREEMENT (scenario x index, shared sequential scale) ----------------
if (DO_AGREEMENT){
  cells <- list()
  for (sc in names(scenarios)) for (v in names(vars_idx)){
    fmd <- file.path(dir_agr, sprintf("meandelta_%s_%s.nc", v, sc))
    frb <- file.path(dir_agr, sprintf("robust_%s_%s.nc",    v, sc))
    if (file.exists(fmd) && file.exists(frb))
      cells[[paste(sc, v)]] <- list(sc = sc, v = v, md = fmd, rb = frb)
    else message("skip (missing): ", sc, " ", v)
  }

  if (length(cells)){
    HI <- if (!is.na(HI_OVERRIDE)) HI_OVERRIDE else
          ceiling(max(vapply(cells, function(z) q_hi(read1(z$md), 0.99), numeric(1))))

    seq_map <- function(z){
      md <- project(read1(z$md), "ESRI:54030")
      rb <- project(read1(z$rb), "ESRI:54030", method = "near")
      df <- to_df(md, "delta"); resx <- res(md)[1]; resy <- res(md)[2]
      pts <- to_df(rb, "rob"); pts <- pts[pts$rob == 0, c("x","y")]   # stipple where NOT robust
      if (nrow(pts) > MAX_STIPPLE) pts <- pts[round(seq(1, nrow(pts), length.out = MAX_STIPPLE)), ]
      ggplot() +
        geom_tile(data = df, aes(x, y, fill = delta), width = resx, height = resy) +
        (if (nrow(pts)) geom_point(data = pts, aes(x, y), size = 0.05, alpha = 0.5, colour = "grey15") else NULL) +
        geom_sf(data = world_54030, fill = NA, colour = "grey20", linewidth = 0.2) +
        coord_sf(crs = 54030, expand = FALSE) +
        scale_fill_gradientn(colours = warm_cols, limits = c(0, HI), oob = squish,
                             name = "\u0394 (\u00B0C)") +
        labs(title = paste0(scenarios[[z$sc]], " \u2014 ", vars_idx[[z$v]])) +
        base_theme()
    }

    plots  <- lapply(cells, seq_map)
    n_scen <- length(scenarios); n_var <- length(vars_idx)
    fig <- wrap_plots(plots, ncol = n_var, guides = "collect") +
      plot_annotation(
        title   = paste0("Inter-model agreement on projected change (", PERIOD_LABEL, ")"),
        caption = "Stippling: cells where fewer than 4/5 models agree on the sign of change.",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                      plot.caption = element_text(hjust = 0.5, size = 11))) &
      theme(legend.position = "bottom")

    ggsave(file.path(dir_fig, "SI_agreement_map.png"), fig,
           width = 7 * n_var, height = 3.4 * n_scen + 1.2, dpi = 300)
    message("DONE -> SI_agreement_map.png  (shared scale 0..", HI, " degC)")
  }
}

# ---------------- BIAS MAPS (validation; diverging, unchanged) ----------------
if (DO_BIAS_MAP){
  div_map <- function(r, title, legend_title){
    rr <- project(r, "ESRI:54030"); lim <- sym_lim_q(r)
    df <- to_df(rr, "v"); resx <- res(rr)[1]; resy <- res(rr)[2]
    ggplot() +
      geom_tile(data = df, aes(x, y, fill = v), width = resx, height = resy) +
      geom_sf(data = world_54030, fill = NA, colour = "grey20", linewidth = 0.2) +
      coord_sf(crs = 54030, expand = FALSE) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                           midpoint = 0, limits = c(-lim, lim), oob = squish, name = legend_title) +
      labs(title = title) + base_theme()
  }
  ps <- lapply(names(vars_bias), function(v){
    f <- file.path(dir_val, paste0("bias_", v, ".nc"))
    if (!file.exists(f)) { message("skip (missing): ", f); return(NULL) }
    ut <- if (v == "humidex_mean") "" else " (\u00B0C)"
    div_map(read1(f), vars_bias[[v]], paste0("ENS \u2212 ERA5", ut))
  })
  ps <- Filter(Negate(is.null), ps)
  if (length(ps)){
    fig_bias <- wrap_plots(ps, ncol = 2) +
      plot_annotation(title = "Bias vs ERA5 (1991\u20132020)",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))
    ggsave(file.path(dir_fig, "Fig_bias_vs_ERA5.png"), fig_bias, width = 12.5, height = 8.4, dpi = 300)
    message("DONE -> Fig_bias_vs_ERA5.png")
  }
}
