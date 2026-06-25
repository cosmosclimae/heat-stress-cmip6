#!/usr/bin/env bash
# ============================================================
# verify_agreement.sh  — sanity-check the "100% agreement" result.
# For each index, independently of the robust-mask logic:
#   (1) prints the 5 per-model mean deltas  -> must be DISTINCT (real 5 models)
#   (2) ensmin of the 5 deltas, then its GLOBAL MINIMUM over land:
#       if > 0  => every model warms at every land cell => agreement is genuine
#   (3) counts land cells where at least one model is non-positive (ensmin<=0)
#       i.e. the only cells where disagreement on sign is even possible.
# ============================================================
set -uo pipefail

ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
SCEN="ssp585"; PERIOD="2071-2100"
MODELS=(IPSL-CM6A-LR MPI-ESM1-2-HR MRI-ESM2-0 GFDL-ESM4 UKESM1-0-LL)
VARS=(wbgt_shade_mean tw_mean hi_mean humidex_mean)

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
f_pres(){ echo "$ROOT/$1/packs/HEATSTRESS_pack_$1_pseudoHist_1991-2020.nc"; }
f_fut(){  echo "$ROOT/$1/packs/HEATSTRESS_pack_$1_${SCEN}_${PERIOD}.nc"; }
num(){ cdo -s outputf,%.5g "$1" | tr -d '[:space:]'; }

for VAR in "${VARS[@]}"; do
  echo "==================== $VAR ($SCEN $PERIOD) ===================="
  deltas=(); i=0
  echo "  per-model mean delta (must differ across models):"
  for m in "${MODELS[@]}"; do
    pc="$TMP/p_$i.nc"; fc="$TMP/f_$i.nc"; dl="$TMP/d_$i.nc"
    cdo -s -f nc4 -O timmean -selname,"$VAR" "$(f_pres "$m")" "$pc"
    cdo -s -f nc4 -O timmean -selname,"$VAR" "$(f_fut  "$m")" "$fc"
    cdo -s -f nc4 -O sub "$fc" "$pc" "$dl"
    cdo -s -f nc4 -O fldmean "$dl" "$TMP/fm.nc"; md=$(num "$TMP/fm.nc")
    printf "    %-16s  mean Δ = %s\n" "$m" "$md"
    deltas+=("$dl"); i=$((i+1))
  done

  # cell-wise minimum delta across the 5 models
  cdo -s -f nc4 -O ensmin "${deltas[@]}" "$TMP/dmin.nc"
  cdo -s -f nc4 -O ensmax "${deltas[@]}" "$TMP/dmax.nc"
  gmin=$(cdo -s outputf,%.5g -fldmin "$TMP/dmin.nc" | tr -d '[:space:]')
  gmax=$(cdo -s outputf,%.5g -fldmax "$TMP/dmax.nc" | tr -d '[:space:]')

  # area-weighted fraction where the least-warming model is <= 0 (0 => none)
  negsum=$(cdo -s outputf,%.5g -fldsum -lec,0 "$TMP/dmin.nc" | tr -d '[:space:]')

  echo "  --------------------------------------------------------"
  echo "  global MIN of per-cell minimum Δ : $gmin   (> 0  => ALL models warm at ALL land cells)"
  echo "  global MAX of per-cell maximum Δ : $gmax"
  echo "  weighted sum where min-model Δ <= 0 (0 => no cell can disagree): $negsum"
  echo ""
done
echo "[done] If every 'global MIN' > 0 and the 5 means differ, the 100% agreement is real."
