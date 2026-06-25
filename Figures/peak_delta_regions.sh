#!/usr/bin/env bash
# Delta of mean PEAK indices (2071-2100 SSP585 minus present) by hotspot region.
# Reports both the regional max and the more robust 99th percentile.
set -uo pipefail
ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
ENS="$ROOT/ENSEMBLE"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# present-day ensemble-mean file (adapt if your naming differs)
PRES="$ENS/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
FUT="$ENS/SSP585/2071-2100/HEATSTRESS_ensmean_ssp585_2071-2100.nc"

# regions as  name:latS,latN,lonW,lonE
REGIONS=(
  "S_America:-20,5,-75,-45"
  "Congo:-5,5,10,30"
  "S_Asia:20,30,70,90"
  "SE_Asia:0,20,95,120"
  "Sahel:10,20,-10,30"
  "ArabianGulf:24,30,48,56"
)

num(){ cdo -s outputf,%.2f "$1" 2>/dev/null | tr -d '[:space:]'; }

for V in wbgt_shade_peak tw_peak; do
  echo "================= $V : Delta peak (2071-2100 SSP585 - present) ================="
  D="$TMP/d_$V.nc"
  cdo -s -f nc4 -O sub \
      -timmean -selname,"$V" "$FUT" \
      -timmean -selname,"$V" "$PRES" "$D" 2>/dev/null

  cdo -s -f nc4 -O fldmax "$D" "$TMP/gmax.nc" 2>/dev/null
  cdo -s -f nc4 -O fldpctl,99 "$D" "$TMP/gp99.nc" 2>/dev/null
  printf "  GLOBAL    max=%s  p99=%s  (deg C)\n" "$(num "$TMP/gmax.nc")" "$(num "$TMP/gp99.nc")"
  echo   "  per-region  (max / p99, deg C):"
  for R in "${REGIONS[@]}"; do
    name="${R%%:*}"; box="${R##*:}"
    # box is latS,latN,lonW,lonE -> sellonlatbox wants lonW,lonE,latS,latN
    sel=$(echo "$box" | awk -F, '{print $3","$4","$1","$2}')
    cdo -s -f nc4 -O fldmax    -sellonlatbox,"$sel" "$D" "$TMP/rmax.nc" 2>/dev/null
    cdo -s -f nc4 -O fldpctl,99 -sellonlatbox,"$sel" "$D" "$TMP/rp99.nc" 2>/dev/null
    printf "    %-12s %5s / %5s\n" "$name" "$(num "$TMP/rmax.nc")" "$(num "$TMP/rp99.nc")"
  done
  echo ""
done
echo "[done] Use p99 (more robust) for the manuscript; max can hit isolated cells."
