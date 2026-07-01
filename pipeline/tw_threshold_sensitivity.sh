#!/usr/bin/env bash
# Sensitivity to the PHYSIOLOGICAL Tw threshold (30 / 32 / 35 C): exposed land
# area (million km^2) with chronic exceedance (>=30 days/yr), by scenario/period.
# Uses ensemble-mean monthly ndays_tw_peak_ge{30,32,35}, collapsed to a 12-month
# climatology, summed to days/yr, masked at >=30 d/yr, area via gridarea.
set -uo pipefail
ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
ENS="$ROOT/ENSEMBLE"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
num(){ cdo -s outputf,%.2f "$1" 2>/dev/null | tr -d '[:space:]'; }

THRS=(30 32 35)
CHRONIC=30
declare -a SCN=(SSP126 SSP245 SSP585)
declare -a PER=(2041-2070 2071-2100)

# gridarea (km^2) from any one field, built once
REF=$(ls "$ENS/SSP585/2071-2100"/HEATSTRESS_ensmean_*2071-2100.nc | head -1)
cdo -s -f nc4 -O gridarea -seltimestep,1 -selname,ndays_tw_peak_ge35 "$REF" "$TMP/area_m2.nc" 2>/dev/null

printf "%-8s %-12s %8s %8s %8s\n" "scen" "period" "Tw>=30" "Tw>=32" "Tw>=35"
for s in "${SCN[@]}"; do
  for p in "${PER[@]}"; do
    f=$(ls "$ENS/$s/$p"/HEATSTRESS_ensmean_*"$p".nc 2>/dev/null | head -1)
    [ -z "$f" ] && { printf "%-8s %-12s   (file not found)\n" "$s" "$p"; continue; }
    line=$(printf "%-8s %-12s" "$s" "$p")
    for t in "${THRS[@]}"; do
      V="ndays_tw_peak_ge${t}"
      # climatology -> days/yr -> chronic mask -> area (km^2) -> sum -> Mkm^2
      cdo -s -f nc4 -O timsum -ymonmean -selname,"$V" "$f" "$TMP/ann.nc" 2>/dev/null
      cdo -s -f nc4 -O gec,$CHRONIC "$TMP/ann.nc" "$TMP/mask.nc" 2>/dev/null
      cdo -s -f nc4 -O fldsum -mul "$TMP/mask.nc" "$TMP/area_m2.nc" "$TMP/a.nc" 2>/dev/null
      km2=$(num "$TMP/a.nc")
      mkm2=$(awk -v x="$km2" 'BEGIN{printf "%.2f", x/1e6/1e6}')   # m2 -> km2(/1e6) -> Mkm2(/1e6)
      line+=$(printf " %8s" "$mkm2")
    done
    echo "$line"
  done
done
echo ""
echo "[done] Area in million km^2. Lower Tw thresholds -> larger chronic area."
