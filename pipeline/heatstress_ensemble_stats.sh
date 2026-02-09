#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 5 ]; then
  echo "Usage: $0 OUT_MEAN OUT_MIN OUT_MAX IN1 IN2 [IN3 ...]" >&2
  exit 1
fi

OUT_MEAN="$1"
OUT_MIN="$2"
OUT_MAX="$3"
shift 3

mkdir -p "$(dirname "$OUT_MEAN")"

echo "[INFO] Checking inputs..."
for f in "$@"; do
  if [ ! -f "$f" ]; then
    echo "[ERROR] Missing input: $f" >&2
    exit 2
  fi
  echo "  - $f"
  # Quick check readability
  cdo -s sinfo "$f" >/dev/null || { echo "[ERROR] Unreadable/unsupported: $f" >&2; exit 3; }
done

# Optional: enforce identical grid/time/var structure by copying to temporary nc4c
# This avoids CDO choking on minor netCDF encoding differences across models.
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "[INFO] Normalising inputs to netCDF4 classic model..."
NORM=()
i=0
for f in "$@"; do
  i=$((i+1))
  nf="${TMPDIR}/in_${i}.nc"
  # copy keeps data but normalises encoding; also drops weird compression/chunking issues
  cdo -L -f nc4c copy "$f" "$nf"
  NORM+=("$nf")
done

echo "[INFO] Computing ensemble stats..."
cdo -L -z zip_5 ensmean "${NORM[@]}" "$OUT_MEAN"
cdo -L -z zip_5 ensmin  "${NORM[@]}" "$OUT_MIN"
cdo -L -z zip_5 ensmax  "${NORM[@]}" "$OUT_MAX"

echo "[OK] Wrote:"
echo "  mean: $OUT_MEAN"
echo "  min : $OUT_MIN"
echo "  max : $OUT_MAX"
