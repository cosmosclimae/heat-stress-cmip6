#!/usr/bin/env bash
set -euo pipefail

# Usage: $0 IN_DIR OUT_DIR
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 IN_DIR OUT_DIR" >&2
  exit 1
fi

IN_DIR="$1"
OUT_ROOT="$2"
mkdir -p "${OUT_ROOT}"

# Inputs (single big files)
TAS="${IN_DIR}/tas_ERA5_day_gr025_1979-2020.nc"
TASMAX="${IN_DIR}/tasmax_ERA5_day_gr025_1979-2020.nc"
HURS="${IN_DIR}/hurs_ERA5_day_gr025_1979-2020.nc"

# Outputs (as before)
OUT_MON_ALL="${OUT_ROOT}/HEATSTRESS_mon_ERA5_gr025_1979-2020.nc"
OUT_MON_9120="${OUT_ROOT}/HEATSTRESS_mon_ERA5_gr025_1991-2020.nc"

# Working subdir
BYYEAR_DIR="${OUT_ROOT}/byyear"
mkdir -p "${BYYEAR_DIR}"

# If final target exists, skip everything
if [ -f "${OUT_MON_9120}" ]; then
  echo "[INFO] Exists ${OUT_MON_9120} -> skip"
  exit 0
fi

# -----------------------------
# Helper: compute one year
# -----------------------------
compute_year() {
  local Y="$1"
  local OUT_Y="${BYYEAR_DIR}/HEATSTRESS_mon_ERA5_gr025_${Y}.nc"

  if [ -f "${OUT_Y}" ]; then
    echo "[INFO] Year ${Y} exists -> skip"
    return 0
  fi

  echo "[INFO] Processing year ${Y}"

  # Extract year slices to avoid loading the full 1979-2020 stack into memory
  local TAS_Y="$(mktemp --suffix=.nc)"
  local TASMAX_Y="$(mktemp --suffix=.nc)"
  local HURS_Y="$(mktemp --suffix=.nc)"
  local TMP_DAY="$(mktemp --suffix=.nc)"
  local TMP_MON="$(mktemp --suffix=.nc)"

  # Select year (fast + stable)
  cdo -L -z zip_5 selyear,${Y} "${TAS}"    "${TAS_Y}"
  cdo -L -z zip_5 selyear,${Y} "${TASMAX}" "${TASMAX_Y}"
  cdo -L -z zip_5 selyear,${Y} "${HURS}"   "${HURS_Y}"

  # Daily indices
# Daily indices
cdo -L -z zip_5 \
  -expr,"\
    # =========================
    # Base variables
    # =========================
    tmean = tas - 273.15; \
    tmax  = tasmax - 273.15; \
    rhm   = max(min(hurs,99),1); \
 
    # =========================
    # Wet-bulb temperature (Stull 2011)
    # =========================
    tw_mean = tmean*atan(0.151977*sqrt(rhm+8.313659)) \
              + atan(tmean+rhm) \
              - atan(rhm-1.676331) \
              + 0.00391838*pow(rhm,1.5)*atan(0.023101*rhm) \
              - 4.686035; \

    tw_peak = tmax*atan(0.151977*sqrt(rhm+8.313659)) \
              + atan(tmax+rhm) \
              - atan(rhm-1.676331) \
              + 0.00391838*pow(rhm,1.5)*atan(0.023101*rhm) \
              - 4.686035; \

    # =========================
    # WBGT (SHADE / no solar load)
    # =========================
    wbgt_shade_mean = 0.7*tw_mean + 0.3*tmean; \
    wbgt_shade_peak = 0.7*tw_peak + 0.3*tmax; \

    # =========================
    # Dew point from (T, RH) and Humidex
    # =========================
    a = 17.625; b = 243.04; \

    gamma_m = (a*tmean)/(b+tmean) + log(rhm/100.); \
    tdew_m  = (b*gamma_m)/(a-gamma_m); \
    e_dew_m = 6.112*exp((17.67*tdew_m)/(tdew_m+243.5)); \
    humidex_mean = tmean + 0.5555*(e_dew_m - 10.); \

    gamma_x = (a*tmax)/(b+tmax) + log(rhm/100.); \
    tdew_x  = (b*gamma_x)/(a-gamma_x); \
    e_dew_x = 6.112*exp((17.67*tdew_x)/(tdew_x+243.5)); \
    humidex_peak = tmax + 0.5555*(e_dew_x - 10.); \

    # =========================
    # Heat Index (Rothfusz) — mean + peak (for consistency)
    # =========================
    Tm_F = tmean*9.0/5.0 + 32.0; \
    Tx_F = tmax *9.0/5.0 + 32.0; \

    HI_F_raw_m = -42.379 \
                 + 2.04901523*Tm_F \
                 + 10.14333127*rhm \
                 - 0.22475541*Tm_F*rhm \
                 - 0.00683783*Tm_F*Tm_F \
                 - 0.05481717*rhm*rhm \
                 + 0.00122874*Tm_F*Tm_F*rhm \
                 + 0.00085282*Tm_F*rhm*rhm \
                 - 0.00000199*Tm_F*Tm_F*rhm*rhm; \
    HI_F_m = (Tm_F >= 80.0 && rhm >= 40.0) ? HI_F_raw_m : Tm_F; \
    hi_mean = (HI_F_m - 32.0)*5.0/9.0; \

    HI_F_raw_x = -42.379 \
                 + 2.04901523*Tx_F \
                 + 10.14333127*rhm \
                 - 0.22475541*Tx_F*rhm \
                 - 0.00683783*Tx_F*Tx_F \
                 - 0.05481717*rhm*rhm \
                 + 0.00122874*Tx_F*Tx_F*rhm \
                 + 0.00085282*Tx_F*rhm*rhm \
                 - 0.00000199*Tx_F*Tx_F*rhm*rhm; \
    HI_F_x = (Tx_F >= 80.0 && rhm >= 40.0) ? HI_F_raw_x : Tx_F; \
    hi_peak = (HI_F_x - 32.0)*5.0/9.0; \

    # =========================
    # Threshold flags (use PEAK for extremes / Fig3)
    # =========================
    f_tw_peak_ge35           = tw_peak           >= 35; \
    f_wbgt_shade_peak_ge30   = wbgt_shade_peak   >= 30; \
    f_wbgt_shade_peak_ge32   = wbgt_shade_peak   >= 32; \
  " \
  -merge "${TAS_Y}" "${TASMAX_Y}" "${HURS_Y}" \
    "${TMP_DAY}"

  # Monthly aggregation: means + counts (days above thresholds)
cdo -L -z zip_5 \
   -merge \
    -monmean -selname,tw_mean,wbgt_shade_mean,hi_mean,humidex_mean "${TMP_DAY}" \
    -monmean -selname,tw_peak,wbgt_shade_peak,hi_peak,humidex_peak "${TMP_DAY}" \
    -monsum  -selname,f_tw_peak_ge35,f_wbgt_shade_peak_ge30,f_wbgt_shade_peak_ge32,f_wbgt_rad_peak_ge30,f_wbgt_rad_peak_ge32 "${TMP_DAY}" \
  "${TMP_MON}"

  # Rename variables (same convention as your original script)
  cdo -L -z zip_5 \
    -chname,\
f_tw_peak_ge35,ndays_tw_peak_ge35,\
f_wbgt_shade_peak_ge30,ndays_wbgt_shade_peak_ge30,\
f_wbgt_shade_peak_ge32,ndays_wbgt_shade_peak_ge32,\
f_wbgt_rad_peak_ge30,ndays_wbgt_rad_peak_ge30,\
f_wbgt_rad_peak_ge32,ndays_wbgt_rad_peak_ge32 \
    "${TMP_MON}" "${OUT_Y}"

  rm -f "${TAS_Y}" "${TASMAX_Y}" "${HURS_Y}" "${TMP_DAY}" "${TMP_MON}"
}

# -----------------------------
# Run per year
# -----------------------------
for Y in $(seq 1991 2020); do
  compute_year "${Y}"
done

# -----------------------------
# Concat all years (ordered) + subset 1991-2020
# -----------------------------
echo "[INFO] Concatenating yearly files -> ${OUT_MON_ALL}"

# Build ordered list to avoid glob-order surprises
YEAR_FILES=()
for Y in $(seq 1991 2020); do
  YEAR_FILES+=("${BYYEAR_DIR}/HEATSTRESS_mon_ERA5_gr025_${Y}.nc")
done

cdo -L -z zip_5 mergetime "${YEAR_FILES[@]}" "${OUT_MON_ALL}"

echo "[INFO] Selecting 1991-2020 -> ${OUT_MON_9120}"
cdo -L -z zip_5 seldate,1991-01-01,2020-12-31 "${OUT_MON_ALL}" "${OUT_MON_9120}"

echo "[OK] ERA5 monthly ready: ${OUT_MON_9120}"
