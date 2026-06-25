#!/usr/bin/env bash
# ============================================================
# agreement_cdo.sh  — inter-model agreement on projected change (CDO only).
# Loops over SCENARIOS; outputs are scenario-tagged so they don't clobber.
#   meandelta_<VAR>_<SCEN>.nc   ensemble-mean change (future - present)
#   robust_<VAR>_<SCEN>.nc      1 = >=K/N models share the sign, else 0
# R only draws (make_maps_from_cdo.r).
#
# WSL: Windows C: is /mnt/c. Per-model packs use LOWERCASE scenario (ssp585);
# ensemble files use UPPERCASE (SSP585). Both handled.
# ============================================================
set -euo pipefail

MODE="permodel"    # "permodel" (rigorous, uses members) or "minmax" (ensemble envelope)

ROOT="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2"
ENS="$ROOT/ENSEMBLE"
OUTDIR="$ROOT/Figures_SI"

SCENARIOS=(SSP245 SSP585)   # EDIT: add SSP126 if you want the three
PERIOD="2071-2100"
N=5
K=4
VARS=(wbgt_shade_mean tw_mean hi_mean humidex_mean)
MODELS=(IPSL-CM6A-LR MPI-ESM1-2-HR MRI-ESM2-0 GFDL-ESM4 UKESM1-0-LL)

mkdir -p "$OUTDIR"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
NNEG=$((N - K))

clim () { cdo -s -f nc4 -O timmean -selname,"$2" "$1" "$3"; }   # add -yearsum for day counts

f_pres () { echo "$ROOT/$1/packs/HEATSTRESS_pack_$1_pseudoHist_1991-2020.nc"; }      # $1 model
f_fut  () { echo "$ROOT/$1/packs/HEATSTRESS_pack_$1_${2}_${PERIOD}.nc"; }            # $1 model $2 scen_lc

agree_permodel () {   # $1 SCEN(uppercase)  $2 SCEN_LC
  local SCEN="$1" SCEN_LC="$2"
  for m in "${MODELS[@]}"; do
    for f in "$(f_pres "$m")" "$(f_fut "$m" "$SCEN_LC")"; do
      [ -f "$f" ] || { echo "[ERROR] missing: $f" >&2; exit 2; }
    done
  done
  for VAR in "${VARS[@]}"; do
    echo "[INFO] permodel: $VAR ($SCEN $PERIOD)"
    local deltas=() posmasks=() i=0
    for m in "${MODELS[@]}"; do
      local pc="$TMP/p_${i}.nc" fc="$TMP/f_${i}.nc" dl="$TMP/d_${i}.nc" pm="$TMP/pos_${i}.nc"
      clim "$(f_pres "$m")" "$VAR" "$pc"
      clim "$(f_fut "$m" "$SCEN_LC")" "$VAR" "$fc"
      cdo -s -f nc4 -O sub "$fc" "$pc" "$dl"
      cdo -s -f nc4 -O gtc,0 "$dl" "$pm"
      deltas+=("$dl"); posmasks+=("$pm"); i=$((i+1))
    done
    cdo -s -f nc4 -O ensmean "${deltas[@]}" "$OUTDIR/meandelta_${VAR}_${SCEN}.nc"
    cdo -s -f nc4 -O enssum  "${posmasks[@]}" "$TMP/npos.nc"
    cdo -s -f nc4 -O gec,"$K"    "$TMP/npos.nc" "$TMP/apos.nc"
    cdo -s -f nc4 -O lec,"$NNEG" "$TMP/npos.nc" "$TMP/aneg.nc"
    cdo -s -f nc4 -O add "$TMP/apos.nc" "$TMP/aneg.nc" "$OUTDIR/robust_${VAR}_${SCEN}.nc"
    echo "  -> meandelta_${VAR}_${SCEN}.nc / robust_${VAR}_${SCEN}.nc"
  done
}

agree_minmax () {   # $1 SCEN(uppercase)
  local SCEN="$1"
  local PM="$ENS/pseudoHist/1991-2020/HEATSTRESS_ensmean_pseudoHist_1991-2020.nc"
  local Pm="$ENS/pseudoHist/1991-2020/HEATSTRESS_ensmin_pseudoHist_1991-2020.nc"
  local PM2="$ENS/pseudoHist/1991-2020/HEATSTRESS_ensmax_pseudoHist_1991-2020.nc"
  local FM="$ENS/$SCEN/$PERIOD/HEATSTRESS_ensmean_${SCEN}_${PERIOD}.nc"
  local Fm="$ENS/$SCEN/$PERIOD/HEATSTRESS_ensmin_${SCEN}_${PERIOD}.nc"
  local FM2="$ENS/$SCEN/$PERIOD/HEATSTRESS_ensmax_${SCEN}_${PERIOD}.nc"
  for f in "$PM" "$Pm" "$PM2" "$FM" "$Fm" "$FM2"; do
    [ -f "$f" ] || { echo "[ERROR] missing: $f (future ensmin/ensmax needed for minmax mode)" >&2; exit 2; }
  done
  for VAR in "${VARS[@]}"; do
    echo "[INFO] minmax: $VAR ($SCEN $PERIOD)"
    clim "$PM" "$VAR" "$TMP/pmean.nc"; clim "$FM" "$VAR" "$TMP/fmean.nc"
    clim "$Pm" "$VAR" "$TMP/pmin.nc";  clim "$PM2" "$VAR" "$TMP/pmax.nc"
    clim "$Fm" "$VAR" "$TMP/fmin.nc";  clim "$FM2" "$VAR" "$TMP/fmax.nc"
    cdo -s -f nc4 -O sub "$TMP/fmean.nc" "$TMP/pmean.nc" "$OUTDIR/meandelta_${VAR}_${SCEN}.nc"
    cdo -s -f nc4 -O gtc,0 -sub "$TMP/fmin.nc" "$TMP/pmax.nc" "$TMP/inc.nc"
    cdo -s -f nc4 -O ltc,0 -sub "$TMP/fmax.nc" "$TMP/pmin.nc" "$TMP/dec.nc"
    cdo -s -f nc4 -O add "$TMP/inc.nc" "$TMP/dec.nc" "$OUTDIR/robust_${VAR}_${SCEN}.nc"
    echo "  -> meandelta_${VAR}_${SCEN}.nc / robust_${VAR}_${SCEN}.nc"
  done
}

for SCEN in "${SCENARIOS[@]}"; do
  SCEN_LC="$(echo "$SCEN" | tr 'A-Z' 'a-z')"
  if [ "$MODE" = "permodel" ]; then agree_permodel "$SCEN" "$SCEN_LC"
  elif [ "$MODE" = "minmax" ]; then agree_minmax "$SCEN"
  else echo "[ERROR] MODE must be permodel or minmax" >&2; exit 1; fi
done
echo "[OK] agreement fields (>=${K}/${N}) -> $OUTDIR"
