#!/usr/bin/env bash
set -euo pipefail

# Usage: $0 IN_DIR OUT_DIR
if [ "$#" -lt 4 ]; then
  echo "Usage: $0 NEX_IN OUT MODEL SCEN" >&2
  exit 1
fi

NEX_IN="$1"
OUT="$2"
MODEL="$3"
SCEN="$4"

IN_DIR="${NEX_IN}/${MODEL}/${SCEN}"
OUT_ROOT="${OUT}/${MODEL}/${SCEN}"

mkdir -p "${OUT_ROOT}"

TAS="${IN_DIR}/tas_ERA5_day_gr025_1979-2020.nc"
TASMAX="${IN_DIR}/tasmax_ERA5_day_gr025_1979-2020.nc"
HURS="${IN_DIR}/hurs_ERA5_day_gr025_1979-2020.nc"

OUT_MON_ALL="${OUT_ROOT}/HEATSTRESS_mon_ERA5_gr025_1979-2020.nc"
OUT_MON_9120="${OUT_ROOT}/HEATSTRESS_mon_ERA5_gr025_1991-2020.nc"

mkdir -p "${OUT_ROOT}"


if [ -f "${OUT_MON_9120}" ]; then
  echo "[INFO] Exists ${OUT_MON_9120} -> skip"
  exit 0
fi

TMP_DAY="$(mktemp --suffix=.nc)"
TMP_MON="$(mktemp --suffix=.nc)"

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
               + 2.04901523*T_F + 10.14333127*rhm \
               - 0.22475541*T_F*rhm \
               - 0.00683783*T_F*T_F - 0.05481717*rhm*rhm \
               + 0.00122874*T_F*T_F*rhm + 0.00085282*T_F*rhm*rhm \
               - 0.00000199*T_F*T_F*rhm*rhm; \
    HI_F = (T_F >= 80.0 && rhm >= 40.0) ? HI_F_raw : T_F; \
    hi = (HI_F - 32.0)*5.0/9.0; \
    f_tw_ge35 = tw >= 35; \
    f_wbgt_ge30 = wbgt >= 30; \
    f_wbgt_ge32 = wbgt >= 32; \
  " \
  -merge "${TAS}" "${TASMAX}" "${HURS}" \
  "${TMP_DAY}"

cdo -L -z zip_5 \
  -merge \
    -monmean -selname,tw,wbgt,hi,humidex "${TMP_DAY}" \
    -monsum  -selname,f_tw_ge35,f_wbgt_ge30,f_wbgt_ge32 "${TMP_DAY}" \
  "${TMP_MON}"

cdo -L -z zip_5 \
  -chname,tw,tw_mean,wbgt,wbgt_mean,hi,hi_mean,humidex,humidex_mean,\
f_tw_ge35,ndays_tw_ge35,f_wbgt_ge30,ndays_wbgt_ge30,f_wbgt_ge32,ndays_wbgt_ge32 \
  "${TMP_MON}" "${OUT_MON_ALL}"

cdo -L -z zip_5 seldate,1991-01-01,2020-12-31 "${OUT_MON_ALL}" "${OUT_MON_9120}"

rm -f "${TMP_DAY}" "${TMP_MON}"
echo "[OK] ERA5 monthly ready: ${OUT_MON_9120}"
