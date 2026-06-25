#!/usr/bin/env bash
# Seasonal DURATION of chronic humid heat from MONTHLY ndays fields.
# FIX: collapse the 30 years to a 12-month climatology (ymonmean) FIRST,
# so all statistics are per-year averages, not summed over 30 years.
set -uo pipefail
ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
ENS="$ROOT/ENSEMBLE"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
num(){ cdo -s outputf,%.2f "$1" 2>/dev/null | tr -d '[:space:]'; }

declare -A VAR=( [Tw35]="ndays_tw_peak_ge35" [WBGT32]="ndays_wbgt_shade_peak_ge32" )
FUT=$(ls "$ENS/SSP585/2071-2100"/HEATSTRESS_ensmean_*2071-2100.nc 2>/dev/null | head -1)
echo "using file: $FUT"
echo "time steps in file: $(cdo -s ntime "$FUT" 2>/dev/null)"; echo ""

HALF_MONTH=15; CHRONIC=30

for H in Tw35 WBGT32; do
  V="${VAR[$H]}"
  echo "================= $H  ($V) ================="

  # 1) 12-month climatology of monthly exceedance-days (mean seasonal cycle)
  cdo -s -f nc4 -O ymonmean -selname,"$V" "$FUT" "$TMP/clim.nc" 2>/dev/null
  echo "   clim time steps (should be 12): $(cdo -s ntime "$TMP/clim.nc" 2>/dev/null)"

  # 2) annual exceedance-days = sum of the 12 climatological months
  cdo -s -f nc4 -O timsum "$TMP/clim.nc" "$TMP/ann.nc" 2>/dev/null

  # 3) chronic mask (>=30 days/yr) -> 1 on chronic cells, missing elsewhere
  cdo -s -f nc4 -O gec,$CHRONIC "$TMP/ann.nc" "$TMP/mask01.nc" 2>/dev/null
  cdo -s -f nc4 -O setrtomiss,-0.5,0.5 "$TMP/mask01.nc" "$TMP/mask.nc" 2>/dev/null
  # area (cells) in mask, for area-weighted means restricted to chronic zone
  frac=$(num <(cdo -s -f nc4 -O fldmean "$TMP/mask01.nc" "$TMP/fr.nc" 2>/dev/null; echo "$TMP/fr.nc") 2>/dev/null)

  # (a) exposed months/yr: # climatological months with >= HALF_MONTH days
  cdo -s -f nc4 -O timsum -gec,$HALF_MONTH "$TMP/clim.nc" "$TMP/months.nc" 2>/dev/null
  cdo -s -f nc4 -O fldmean -mul "$TMP/months.nc" "$TMP/mask.nc" "$TMP/mo.nc" 2>/dev/null

  # (b) annual exceedance-days within chronic cells: mean & max
  cdo -s -f nc4 -O fldmean -mul "$TMP/ann.nc" "$TMP/mask.nc" "$TMP/dm.nc" 2>/dev/null
  cdo -s -f nc4 -O fldmax  -mul "$TMP/ann.nc" "$TMP/mask.nc" "$TMP/dx.nc" 2>/dev/null

  # (c) most-exposed single climatological month (days): mean & max
  cdo -s -f nc4 -O timmax "$TMP/clim.nc" "$TMP/pk.nc" 2>/dev/null
  cdo -s -f nc4 -O fldmean -mul "$TMP/pk.nc" "$TMP/mask.nc" "$TMP/pm.nc" 2>/dev/null
  cdo -s -f nc4 -O fldmax  -mul "$TMP/pk.nc" "$TMP/mask.nc" "$TMP/px.nc" 2>/dev/null

  printf "  exposed months/yr (>=%d d/month) : mean = %s  / 12\n" "$HALF_MONTH" "$(num "$TMP/mo.nc")"
  printf "  annual exceedance-days           : mean = %s   max = %s\n" "$(num "$TMP/dm.nc")" "$(num "$TMP/dx.nc")"
  printf "  most-exposed month (days)        : mean = %s   max = %s\n" "$(num "$TMP/pm.nc")" "$(num "$TMP/px.nc")"
  echo ""
done
echo "[done] Seasonal/cumulative durations from monthly counts (NOT consecutive-day spells)."
