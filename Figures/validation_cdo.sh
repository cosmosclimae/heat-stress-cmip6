#!/usr/bin/env bash
# ============================================================
# validation_cdo.sh
# Fast ERA5 validation metrics with CDO (no R for the heavy part).
# For each index/threshold field, computes (area-weighted by CDO):
#   - bias  = fldmean(ENS - ERA5)
#   - RMSE  = sqrt(fldmean((ENS - ERA5)^2))
#   - MAE   = fldmean(|ENS - ERA5|)
#   - r     = fldcor(ENS, ERA5)            (spatial Pearson correlation)
# Writes: ERA5_validation_metrics.csv
#         bias_<VAR>.nc                    (bias field, for the map in R)
#
# Answers reviewer R2 #3. fldmean/fldcor use true grid-cell-area weights.
#
# Usage:
#   ./validation_cdo.sh
# (edit the PATHS block below first)
# ============================================================
set -euo pipefail

# ---------------- PATHS (EDIT) ----------------
# WSL mounts Windows C: at /mnt/c (your figures.r uses C:/Data/...).
ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
F_ENS="$ROOT/ENSEMBLE/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
F_ERA="$ROOT/ERA5/packs/HEATSTRESS_mon_ERA5_gr025_1991-2020.nc"
OUTDIR="$ROOT/Validation"
LAND_ONLY=1                 # 1 = restrict stats to land (recommended), 0 = global
# Optional: a land-fraction file (CMIP6 sftlf, %), same grid as ENS. If empty,
# a topography>0 land mask is built automatically.
SFTLF=""                    # e.g. "$ROOT/fixed/sftlf_ENSgrid.nc"  (var: sftlf)
# Conditional relative bias for day-count fields is computed only where the event
# is climatologically frequent in ERA5 (>= FREQ_MIN days/yr). Indices get "NA".
FREQ_MIN=30
# ----------------------------------------------

mkdir -p "$OUTDIR"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# Variables: mean indices (continuous -> remapbil) and day-counts (-> remapnn)
MEAN_VARS=(wbgt_shade_mean tw_mean hi_mean humidex_mean)
MEAN_UNIT=("degC" "degC" "degC" "-")
DAYS_VARS=(ndays_wbgt_shade_peak_ge30 ndays_wbgt_shade_peak_ge32 ndays_tw_peak_ge35)

# ---------------- land mask (once, on ENS grid) ----------------
LANDMASK=""
if [ "$LAND_ONLY" -eq 1 ]; then
  LANDMASK="$TMP/land.nc"
  if [ -n "$SFTLF" ]; then
    # sftlf in % -> land where > 50%
    cdo -s -f nc4 -O gtc,50 "$SFTLF" "$LANDMASK"
  else
    # topography > 0 as land proxy (misses below-sea-level land; fine for global stats)
    cdo -s -f nc4 -O remapbil,"$F_ENS" -topo "$TMP/topo_global.nc"
    cdo -s -f nc4 -O gtc,0 "$TMP/topo_global.nc" "$LANDMASK"
  fi
fi

apply_mask () {  # $1 in  $2 out
  if [ "$LAND_ONLY" -eq 1 ]; then
    cdo -s -f nc4 -O ifthen "$LANDMASK" "$1" "$2"
  else
    cp "$1" "$2"
  fi
}

num () { cdo -s -f nc4 -O outputf,%.6g "$1" | tr -d '[:space:]'; }  # 1-value nc -> scalar

CSV="$OUTDIR/ERA5_validation_metrics.csv"
echo "variable,unit,bias,rmse,mae,r,cond_rel_bias_pct" > "$CSV"

process () {   # $1 var  $2 unit  $3 climcmd("mean"|"days")  $4 remap("remapbil"|"remapnn")
  local VAR="$1" UNIT="$2" KIND="$3" REMAP="$4"
  echo "[INFO] $VAR ($KIND)"

  local ens_clim="$TMP/ens_${VAR}.nc" era_clim="$TMP/era_${VAR}.nc" era_rg="$TMP/era_${VAR}_rg.nc"

  if [ "$KIND" = "mean" ]; then
    cdo -s -f nc4 -O timmean -selname,"$VAR" "$F_ENS" "$ens_clim"
    cdo -s -f nc4 -O timmean -selname,"$VAR" "$F_ERA" "$era_clim"
  else
    cdo -s -f nc4 -O timmean -yearsum -selname,"$VAR" "$F_ENS" "$ens_clim"
    cdo -s -f nc4 -O timmean -yearsum -selname,"$VAR" "$F_ERA" "$era_clim"
  fi

  # regrid ERA5 onto ENS grid
  cdo -s -f nc4 -O "$REMAP","$ens_clim" "$era_clim" "$era_rg"

  # mask both fields identically
  local ens_m="$TMP/ensm_${VAR}.nc" era_m="$TMP/eram_${VAR}.nc"
  apply_mask "$ens_clim" "$ens_m"
  apply_mask "$era_rg"  "$era_m"

  # bias field (kept for the map)
  local bias="$OUTDIR/bias_${VAR}.nc"
  cdo -s -f nc4 -O sub "$ens_m" "$era_m" "$bias"

  # scalar metrics (absolute; two-step writes, robust on Git-Bash/WSL)
  local b rmse mae c crel
  cdo -s -f nc4 -O fldmean "$bias" "$TMP/m_b.nc";              b=$(num "$TMP/m_b.nc")
  cdo -s -f nc4 -O sqrt -fldmean -sqr "$bias" "$TMP/m_r.nc";   rmse=$(num "$TMP/m_r.nc")
  cdo -s -f nc4 -O fldmean -abs "$bias" "$TMP/m_m.nc";         mae=$(num "$TMP/m_m.nc")
  cdo -s -f nc4 -O fldcor "$ens_m" "$era_m" "$TMP/m_c.nc";     c=$(num "$TMP/m_c.nc")

  # conditional relative bias (day-count fields only): regional NMB =
  # mean(bias) / mean(ERA5 count) over cells where ERA5 >= FREQ_MIN d/yr.
  # Meaningless for interval-scale indices (arbitrary zero) -> NA.
  if [ "$KIND" = "days" ]; then
    cdo -s -f nc4 -O gec,"$FREQ_MIN" "$era_m" "$TMP/m_fm.nc"          # frequent-event mask
    cdo -s -f nc4 -O fldmean -ifthen "$TMP/m_fm.nc" "$bias"  "$TMP/m_bb.nc"
    cdo -s -f nc4 -O fldmean -ifthen "$TMP/m_fm.nc" "$era_m" "$TMP/m_ee.nc"
    local bb ee; bb=$(num "$TMP/m_bb.nc"); ee=$(num "$TMP/m_ee.nc")
    crel=$(awk -v b="$bb" -v e="$ee" 'BEGIN{
        ab=(b<0?-b:b); ae=(e<0?-e:e);
        if (!(e==e) || !(b==b) || ae<1e-6 || ae>1e6 || ab>1e6) print "NA";
        else printf "%.1f", 100*b/e }')
  else
    crel="NA"
  fi

  echo "${VAR},${UNIT},${b},${rmse},${mae},${c},${crel}" >> "$CSV"
}

for i in "${!MEAN_VARS[@]}"; do
  process "${MEAN_VARS[$i]}" "${MEAN_UNIT[$i]}" "mean" "remapbil"
done
for v in "${DAYS_VARS[@]}"; do
  process "$v" "days/yr" "days" "remapnn"
done

echo "[OK] -> $CSV"
cat "$CSV"
