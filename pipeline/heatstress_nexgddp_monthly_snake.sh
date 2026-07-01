#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################
# Heat-stress monthly (CDO) - NEX-GDDP-CMIP6
# Snakemake-compatible version with stable outputs.
#
# Inputs (daily):
#   IN_ROOT/[model]/[scenario]/[var]/
#     [var]_day_[model]_[scen]_<member+grid>_<year>_v2.0.nc
#
# Uses only: tas (K), tasmax (K), hurs (% 0..100)
#
# Scenario time windows (ENFORCED):
#   historical: 1991-01-01 .. 2014-12-31
#   ssp126    : 2015-01-01 .. 2100-12-31
#   ssp245    : 2020-01-01 .. 2100-12-31
#   ssp585    : 2020-01-01 .. 2100-12-31
#
# Outputs (monthly, final sliced):
#   OUT_ROOT/MODEL/SCEN/HEATSTRESS_mon_MODEL_SCEN_START-END.nc
############################################

# Usage: $0 IN_ROOT OUT_ROOT MODEL SCEN
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 IN_ROOT OUT_ROOT MODEL SCEN" >&2
  exit 1
fi

IN_ROOT="$1"
OUT_ROOT="$2"
MODEL="$3"
SCEN="$4"
VERSION_TAG="v2.0"

case "${SCEN}" in
  historical)
    START="1991-01-01"; END="2014-12-31"
    ;;
  ssp126)
    START="2015-01-01"; END="2100-12-31"
    ;;
  ssp245|ssp585)
    START="2020-01-01"; END="2100-12-31"
    ;;
  *)
    echo "[ERROR] Unsupported scenario: ${SCEN}" >&2
    exit 2
    ;;
esac

IN_DIR="${IN_ROOT}/${MODEL}/${SCEN}"
TASMAX_DIR="${IN_DIR}/tasmax"
TAS_DIR="${IN_DIR}/tas"
HURS_DIR="${IN_DIR}/hurs"
SFCWIND_DIR="${IN_DIR}/sfcWind"
RSDS_DIR="${IN_DIR}/rsds"

[ -d "${SFCWIND_DIR}" ] || { echo "[ERROR] Missing ${SFCWIND_DIR}" >&2; exit 3; }
[ -d "${RSDS_DIR}" ]    || { echo "[ERROR] Missing ${RSDS_DIR}" >&2; exit 3; }
[ -d "${TASMAX_DIR}" ] || { echo "[ERROR] Missing ${TASMAX_DIR}" >&2; exit 3; }
[ -d "${TAS_DIR}" ]    || { echo "[ERROR] Missing ${TAS_DIR}" >&2; exit 3; }
[ -d "${HURS_DIR}" ]   || { echo "[ERROR] Missing ${HURS_DIR}" >&2; exit 3; }

OUT_DIR="${OUT_ROOT}/${MODEL}/${SCEN}"
mkdir -p "${OUT_DIR}"

FINAL_OUT="${OUT_DIR}/HEATSTRESS_mon_${MODEL}_${SCEN}_${START:0:4}-${END:0:4}.nc"
if [ -f "${FINAL_OUT}" ]; then
  echo "[INFO] Exists ${FINAL_OUT} -> skip"
  exit 0
fi

echo "==== Heat-stress monthly: model=${MODEL} scen=${SCEN} window=${START}..${END} ===="

YEARLY_MON_FILES=()

for TASMAX_FILE in "${TASMAX_DIR}/tasmax_day_${MODEL}_${SCEN}_"*"_${VERSION_TAG}.nc"; do
  BASENAME="$(basename "${TASMAX_FILE}")"
  PREFIX="tasmax_day_${MODEL}_${SCEN}_"
  SUFFIX="${BASENAME#${PREFIX}}"

  TAS_FILE="${TAS_DIR}/tas_day_${MODEL}_${SCEN}_${SUFFIX}"
  HURS_FILE="${HURS_DIR}/hurs_day_${MODEL}_${SCEN}_${SUFFIX}"
  SFCWIND_FILE="${SFCWIND_DIR}/sfcWind_day_${MODEL}_${SCEN}_${SUFFIX}"
  RSDS_FILE="${RSDS_DIR}/rsds_day_${MODEL}_${SCEN}_${SUFFIX}"

  if [ ! -f "${TAS_FILE}" ] || [ ! -f "${HURS_FILE}" ]; then
    echo "  [WARN] Missing matching tas/hurs for ${BASENAME} -> skip" >&2
    continue
  fi

  YEAR="$(echo "${BASENAME}" | sed -n 's/.*_\([0-9]\{4\}\)_v2\.0\.nc$/\1/p')"
  [ -n "${YEAR}" ] || { echo "  [WARN] No year parsed from ${BASENAME} -> skip" >&2; continue; }

  OUT_MON_YEAR="${OUT_DIR}/HEATSTRESS_mon_${MODEL}_${SCEN}_${YEAR}.nc"
  YEARLY_MON_FILES+=("${OUT_MON_YEAR}")

  if [ -f "${OUT_MON_YEAR}" ]; then
    continue
  fi

  TMP_DAY="$(mktemp --suffix=.nc)"
  TMP_MON="$(mktemp --suffix=.nc)"

cdo -L -z zip_5 \
  -expr,"\
    # =========================
    # Base variables
    # =========================
    tmean = tas - 273.15; \
    tmax  = tasmax - 273.15; \
    rhm   = max(min(hurs,99),1); \
    wind  = max(sfcWind, 0.1); \
    sw    = max(rsds, 0.0); \

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
    # Optional: radiation+wind proxy (SENSITIVITY ONLY)
    # =========================
    k_sw = 0.006; \
    damp = 1.0/(1.0 + 0.25*wind); \
    t_sun = tmax + k_sw*sw*damp; \
    wbgt_rad_proxy_peak = 0.7*tw_peak + 0.3*t_sun; \

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
    f_tw_peak_ge32           = tw_peak           >= 32; \
    f_tw_peak_ge30           = tw_peak           >= 30; \
    f_wbgt_shade_peak_ge30   = wbgt_shade_peak   >= 30; \
    f_wbgt_shade_peak_ge32   = wbgt_shade_peak   >= 32; \
    f_wbgt_rad_peak_ge30     = wbgt_rad_proxy_peak >= 30; \
    f_wbgt_rad_peak_ge32     = wbgt_rad_proxy_peak >= 32; \
  " \
  -merge "${TAS_FILE}" "${TASMAX_FILE}" "${HURS_FILE}" "${SFCWIND_FILE}" "${RSDS_FILE}" \
  "${TMP_DAY}"

cdo -L -z zip_5 \
  -merge \
    -monmean -selname,tw_mean,wbgt_shade_mean,hi_mean,humidex_mean "${TMP_DAY}" \
    -monmean -selname,tw_peak,wbgt_shade_peak,hi_peak,humidex_peak "${TMP_DAY}" \
    -monsum  -selname,f_tw_peak_ge30,f_tw_peak_ge32,f_tw_peak_ge35,f_wbgt_shade_peak_ge30,f_wbgt_shade_peak_ge32,f_wbgt_rad_peak_ge30,f_wbgt_rad_peak_ge32 "${TMP_DAY}" \
  "${TMP_MON}"


  cdo -L -z zip_5 \
  -chname,\
f_tw_peak_ge35,ndays_tw_peak_ge35,\
f_wbgt_shade_peak_ge30,ndays_wbgt_shade_peak_ge30,\
f_wbgt_shade_peak_ge32,ndays_wbgt_shade_peak_ge32,\
f_wbgt_rad_peak_ge30,ndays_wbgt_rad_peak_ge30,\
f_wbgt_rad_peak_ge32,ndays_wbgt_rad_peak_ge32 \
  "${TMP_MON}" "${OUT_MON_YEAR}"

  rm -f "${TMP_DAY}" "${TMP_MON}"
done

if [ "${#YEARLY_MON_FILES[@]}" -eq 0 ]; then
  echo "[ERROR] No yearly outputs created for ${MODEL}/${SCEN}" >&2
  exit 4
fi

TMP_MERGED="$(mktemp --suffix=.nc)"
cdo -L -z zip_5 mergetime "${YEARLY_MON_FILES[@]}" "${TMP_MERGED}"

cdo -L -z zip_5 seldate,"${START}","${END}" "${TMP_MERGED}" "${FINAL_OUT}"

rm -f "${TMP_MERGED}" "${YEARLY_MON_FILES[@]}"

echo "[OK] Wrote ${FINAL_OUT}"
