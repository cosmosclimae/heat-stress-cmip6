#!/usr/bin/env bash
# ============================================================
# precompute_climatos.sh  (generic)
# Precompute every single-layer climatology that figures.r AND Figure_Supp.r
# may need, in CDO, so R only reads + plots.
#
# For each ensemble file (ensmean/ensmin/ensmax, all scenarios/periods) and ERA5:
#   MEAN_<basename>__<var>.nc = timmean -selname,var          (mean & peak indices)
#   DAYS_<basename>__<var>.nc = timmean -yearsum -selname,var (exceedance-day fields)
# Variables absent from a file are skipped (e.g. peak indices not in ERA5).
#
# The R helpers rebuild these names from basename(ncfile)+var.
# NOTE: MEAN uses plain timmean (days-of-month weighting negligible; cancels in deltas).
# ============================================================
set -uo pipefail   # not -e: missing vars are handled per-call

ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
OUT="$ROOT/CLIMATOS"
ENS="$ROOT/ENSEMBLE"
ERA="$ROOT/ERA5/packs/HEATSTRESS_mon_ERA5_gr025_1991-2020.nc"
mkdir -p "$OUT"

VARS_MEAN="wbgt_shade_mean tw_mean hi_mean humidex_mean wbgt_shade_peak tw_peak"
VARS_DAYS="ndays_wbgt_shade_peak_ge30 ndays_wbgt_shade_peak_ge32 ndays_tw_peak_ge35"

process_file () {
  local f="$1"; [ -f "$f" ] || return 0
  local base; base="$(basename "$f" .nc)"
  local avail; avail="$(cdo -s showname "$f" 2>/dev/null | tr ' ' '\n')"

  for v in $VARS_MEAN; do
    printf '%s\n' "$avail" | grep -qx "$v" || continue
    local out="$OUT/MEAN_${base}__${v}.nc"
    [ -f "$out" ] && continue
    cdo -s -f nc4 -O timmean -selname,"$v" "$f" "$out" \
      && echo "  MEAN $base :: $v" || echo "  [fail] MEAN $base :: $v"
  done

  for v in $VARS_DAYS; do
    printf '%s\n' "$avail" | grep -qx "$v" || continue
    local out="$OUT/DAYS_${base}__${v}.nc"
    [ -f "$out" ] && continue
    cdo -s -f nc4 -O timmean -yearsum -selname,"$v" "$f" "$out" \
      && echo "  DAYS $base :: $v" || echo "  [fail] DAYS $base :: $v"
  done
}

echo "[INFO] scanning ensemble files under $ENS ..."
shopt -s nullglob
n=0
for f in "$ENS"/*/*/HEATSTRESS_ens*.nc; do process_file "$f"; n=$((n+1)); done
echo "[INFO] ERA5 ..."
process_file "$ERA"

echo "[OK] processed $n ensemble files + ERA5"
echo "[OK] climatologies -> $OUT  ($(ls -1 "$OUT" | wc -l) files)"
