#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash heatstress_ensemble_stats.sh OUT_MEAN.nc OUT_MIN.nc OUT_MAX.nc IN1.nc IN2.nc ...
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 OUT_MEAN OUT_MIN OUT_MAX IN1 IN2 [IN3 ...]" >&2
  exit 1
fi

OUT_MEAN="$1"
OUT_MIN="$2"
OUT_MAX="$3"
shift 3

mkdir -p "$(dirname "$OUT_MEAN")"

# Quick sanity: fail fast if an input is missing
for f in "$@"; do
  [ -f "$f" ] || { echo "[ERROR] Missing input: $f" >&2; exit 2; }
done

# Ensemble statistics
cdo -L -z zip_5 ensmean "$@" "$OUT_MEAN"
cdo -L -z zip_5 ensmin  "$@" "$OUT_MIN"
cdo -L -z zip_5 ensmax  "$@" "$OUT_MAX"

echo "[OK] Wrote:"
echo "  mean: $OUT_MEAN"
echo "  min : $OUT_MIN"
echo "  max : $OUT_MAX"
