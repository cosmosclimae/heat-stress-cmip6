#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################
# Heat-stress monthly (CDO) - NEX-GDDP-CMIP6
# Snakemake-compatible version with stable outputs.
#
# Inputs (daily):
#   /mnt/e/CMIP6/NEX-GDDP/[model]/[scenario]/[var]/
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
#   ${OUT_ROOT}/${MODEL}/${SCEN}/HEATSTRESS_mon_${MODEL}_${SCEN}_${START}-${END}.nc
# with START/END fixed per scenario above.
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

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 MODEL SCEN" >&2
  exit 1
fi

MODEL="$1"
SCEN="$2"

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

[ -d "${TASMAX_DIR}" ] || { echo "[ERROR] Missing ${TASMAX_DIR}" >&2; exit 3; }
[ -d "${TAS_DIR}" ] || { echo "[ERROR] Missing ${TAS_DIR}" >&2; exit 3; }
[ -d "${HURS_DIR}" ] || { echo "[ERROR] Missing ${HURS_DIR}" >&2; exit 3; }

OUT_DIR="${OUT_ROOT}/${MODEL}/${SCEN}"
mkdir -p "${OUT_DIR}"

FINAL_OUT="${OUT_DIR}/HEATSTRESS_mon_${MODEL}_${SCEN}_${START:0:4}-${END:0:4}.nc"
if [ -f "${FINAL_OUT}" ]; then
  echo "[INFO] Exists ${FINAL_OUT} -> skip"
  exit 0
fi

echo "==== Heat-stress monthly: model=${MODEL} scen=${SCEN} window=${START}..${END} ===="

YEARLY_MON_FILES=()

# Loop tasmax yearly files
for TASMAX_FILE in "${TASMAX_DIR}/tasmax_day_${MODEL}_${SCEN}_"*"_${VERSION_TAG}.nc"; do
  BASENAME="$(basename "${TASMAX_FILE}")"
  PREFIX="tasmax_day_${MODEL}_${SCEN}_"
  SUFFIX="${BASENAME#${PREFIX}}"

  TAS_FILE="${TAS_DIR}/tas_day_${MODEL}_${SCEN}_${SUFFIX}"
  HURS_FILE="${HURS_DIR}/hurs_day_${MODEL}_${SCEN}_${SUFFIX}"

  if [ ! -f "${TAS_FILE}" ] || [ ! -f "${HURS_FILE}" ]; then
    echo "  [WARN] Missing matching tas/hurs for ${BASENAME} -> skip" >&2
    continue
  fi

  # Year extraction (robust)
  YEAR="$(echo "${BASENAME}" | grep -oE '_[0-9]{4}_'"${VERSION_TAG//./\\.}"'\.nc$' | tr -dc '0-9')"
  [ -n "${YEAR}" ] || { echo "  [WARN] No year parsed from ${BASENAME} -> skip" >&2; continue; }

  OUT_MON_YEAR="${OUT_DIR}/heatstress_mon_${MODEL}_${SCEN}_${YEAR}.nc"
  YEARLY_MON_FILES+=("${OUT_MON_YEAR}")

  if [ -f "${OUT_MON_YEAR}" ]; then
    continue
  fi

  TMP_DAY="$(mktemp --suffix=.nc)"
  TMP_MON="$(mktemp --suffix=.nc)"

  # Daily indices + flags
  cdo -L -z zip_5 \
    -expr,"\
      tmean = tas - 273.15; \
      tmax  = tasmax - 273.15; \
      rhm   = max(min(hurs,99),1); \

      tw = tmean*atan(0.151977*sqrt(rhm+8.313659)) + atan(tmean+rhm) - atan(rhm-1.676331) \
           + 0.00391838*pow(rhm,1.5)*atan(0.023101*rhm) - 4.686035; \

      wbgt = 0.7*tw + 0.3*tmean; \

      a = 17.625; b = 243.04; \
      gamma = (a*tmean)/(b+tmean) + log(rhm/100.); \
      tdew = (b*gamma)/(a-gamma); \
      e_dew = 6.112*exp((17.67*tdew)/(tdew+243.5)); \
      humidex = tmean + 0.5555*(e_dew - 10.); \

      T_F = tmax*9.0/5.0 + 32.0; \
      HI_F_raw = -42.379 \
                 + 2.04901523*T_F \
                 + 10.14333127*rhm \
                 - 0.22475541*T_F*rhm \
                 - 0.00683783*T_F*T_F \
                 - 0.05481717*rhm*rhm \
                 + 0.00122874*T_F*T_F*rhm \
                 + 0.00085282*T_F*rhm*rhm \
                 - 0.00000199*T_F*T_F*rhm*rhm; \
      HI_F = (T_F >= 80.0 && rhm >= 40.0) ? HI_F_raw : T_F; \
      hi = (HI_F - 32.0)*5.0/9.0; \

      f_tw_ge35   = tw   >= 35; \
      f_wbgt_ge30 = wbgt >= 30; \
      f_wbgt_ge32 = wbgt >= 32; \
    " \
    -merge "${TAS_FILE}" "${TASMAX_FILE}" "${HURS_FILE}" \
    "${TMP_DAY}"

  # Monthly aggregate
  cdo -L -z zip_5 \
    -merge \
      -monmean -selname,tw,wbgt,hi,humidex "${TMP_DAY}" \
      -monsum  -selname,f_tw_ge35,f_wbgt_ge30,f_wbgt_ge32 "${TMP_DAY}" \
    "${TMP_MON}"

  # Rename + write yearly monthly
  cdo -L -z zip_5 \
    -chname,tw,tw_mean,wbgt,wbgt_mean,hi,hi_mean,humidex,humidex_mean,\
f_tw_ge35,ndays_tw_ge35,f_wbgt_ge30,ndays_wbgt_ge30,f_wbgt_ge32,ndays_wbgt_ge32 \
    "${TMP_MON}" "${OUT_MON_YEAR}"

  rm -f "${TMP_DAY}" "${TMP_MON}"
done

if [ "${#YEARLY_MON_FILES[@]}" -eq 0 ]; then
  echo "[ERROR] No yearly outputs created for ${MODEL}/${SCEN}" >&2
  exit 4
fi

# Merge time then slice exact window
TMP_MERGED="$(mktemp --suffix=.nc)"
cdo -L -z zip_5 mergetime "${YEARLY_MON_FILES[@]}" "${TMP_MERGED}"

# Slice to scenario window (this enforces your dates)
cdo -L -z zip_5 seldate,"${START}","${END}" "${TMP_MERGED}" "${FINAL_OUT}"

rm -f "${TMP_MERGED}"

echo "[OK] Wrote ${FINAL_OUT}"
