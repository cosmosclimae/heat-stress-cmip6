#!/usr/bin/env bash
# Area-weighted % of land where ALL 5 models warm (npos=5) vs >=4/5 (robust).
# fldmean of a 0/1 mask IS the area-weighted fraction -> no gridarea needed.
# HDF5 '_QuantizeBitRound...' messages are a harmless netCDF version-mismatch
# warning (files written with quantization); silenced on reads.
set -uo pipefail
ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
MODELS=(IPSL-CM6A-LR MPI-ESM1-2-HR MRI-ESM2-0 GFDL-ESM4 UKESM1-0-LL)
SCEN="ssp585"; PERIOD="2071-2100"
VARS=(wbgt_shade_mean tw_mean hi_mean humidex_mean)
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
num(){ cdo -s outputf,%.6g "$1" 2>/dev/null | tr -d '[:space:]'; }

printf "%-16s  all5(npos=5)   >=4/5(robust)\n" "index"
for VAR in "${VARS[@]}"; do
  pos=()
  for m in "${MODELS[@]}"; do
    P="$ROOT/$m/packs/HEATSTRESS_pack_${m}_pseudoHist_1991-2020.nc"
    F="$ROOT/$m/packs/HEATSTRESS_pack_${m}_${SCEN}_${PERIOD}.nc"
    cdo -s -f nc4 -O sub -timmean -selname,"$VAR" "$F" -timmean -selname,"$VAR" "$P" "$TMP/d.nc" 2>/dev/null
    cdo -s -f nc4 -O gtc,0 "$TMP/d.nc" "$TMP/pos_$m.nc" 2>/dev/null
    pos+=("$TMP/pos_$m.nc")
  done
  cdo -s -f nc4 -O enssum "${pos[@]}" "$TMP/npos.nc" 2>/dev/null
  cdo -s -f nc4 -O fldmean -eqc,5 "$TMP/npos.nc" "$TMP/f5.nc" 2>/dev/null; f5=$(num "$TMP/f5.nc")
  cdo -s -f nc4 -O fldmean -gec,4 "$TMP/npos.nc" "$TMP/f4.nc" 2>/dev/null; f4=$(num "$TMP/f4.nc")
  awk -v v="$VAR" -v f5="$f5" -v f4="$f4" 'BEGIN{
    if (f5+0!=f5 || f4+0!=f4) { printf "%-16s  (read failed)\n", v }
    else printf "%-16s  %8.3f%%     %8.3f%%\n", v, 100*f5, 100*f4
  }'
done
