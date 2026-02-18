#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

BASE="/mnt/e/CMIP6/NEX-GDDP"
VAR="hurs"
SCEN="ssp126"
OUTDIR="${BASE}/_DERIVED/${VAR}"
TMPDIR="/tmp/hurs_work"
mkdir -p "${OUTDIR}" "${TMPDIR}"

models=(
  "GFDL-ESM4"
  "MPI-ESM1-2-HR"
  "MRI-ESM2-0"
  "IPSL-CM6A-LR"
  "UKESM1-0-LL"
)

# --- helper: collect files whose name contains .YYYY. or _YYYY_ or YYYY-
# Adjust patterns if your filenames differ.
collect_year_files() {
  local dir="$1"; shift
  local -a years=("$@")
  local -a files=()
  local y

  for y in "${years[@]}"; do
    # Try multiple common patterns
    files+=( "${dir}"/*"${y}"*.nc )
  done

  # Filter only existing files (nullglob already handles missing globs)
  if (( ${#files[@]} == 0 )); then
    echo "ERROR: No files matched in ${dir} for years: ${years[*]}" >&2
    echo "Check filename pattern. Example listing:" >&2
    ls -1 "${dir}"/*.nc 2>/dev/null | head -n 10 >&2 || true
    exit 2
  fi

  # Return list on stdout, one per line
  printf "%s\n" "${files[@]}"
}

for m in "${models[@]}"; do
  echo "======================================"
  echo "MODEL: ${m}"
  echo "======================================"

  HIST_DIR="${BASE}/${m}/historical/${VAR}"
  SSP_DIR="${BASE}/${m}/${SCEN}/${VAR}"

  HIST_OUT="${TMPDIR}/${VAR}_${m}_hist_1991-2014.nc"
  SSP_OUT="${TMPDIR}/${VAR}_${m}_ssp_2015-2020.nc"
  MERGE_OUT="${TMPDIR}/${VAR}_${m}_1991-2020.nc"
  FINAL_OUT="${OUTDIR}/${VAR}_mean_${m}_1991-2020.nc"

  # Build explicit file lists by year (file-based selection)
  mapfile -t HIST_FILES < <(collect_year_files "${HIST_DIR}" $(seq 1991 2014))
  mapfile -t SSP_FILES  < <(collect_year_files "${SSP_DIR}"  2015 2016 2017 2018 2019 2020)

  echo "Hist files: ${#HIST_FILES[@]}"
  echo "SSP  files: ${#SSP_FILES[@]}"

  # 1) merge hist 1991-2014
  cdo -L -mergetime "${HIST_FILES[@]}" "${HIST_OUT}"

  # 2) merge ssp 2015-2020
  cdo -L -mergetime "${SSP_FILES[@]}" "${SSP_OUT}"

  # 3) merge both
  cdo -L -mergetime "${HIST_OUT}" "${SSP_OUT}" "${MERGE_OUT}"

  # 4) mean 1991-2020
  cdo -L -timmean "${MERGE_OUT}" "${FINAL_OUT}"

  echo "DONE -> ${FINAL_OUT}"
done

echo "ALL MODELS DONE."

