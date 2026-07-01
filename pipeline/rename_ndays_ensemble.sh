#!/usr/bin/env bash
# Rename f_tw_peak_ge30/ge32 -> ndays_tw_peak_ge30/ge32 in ENSEMBLE files only.
# In-place, idempotent: only renames if the source name is actually present.
set -uo pipefail
ENS="/mnt/c/Data/HEAT_STRESS_MONTHLY_V2/ENSEMBLE"

# pairs: oldname newname
PAIRS=(
  "f_tw_peak_ge30 ndays_tw_peak_ge30"
  "f_tw_peak_ge32 ndays_tw_peak_ge32"
)

# all ENSEMBLE netCDF files (ensmean/ensmin/ensmax, all scenarios/periods)
mapfile -t FILES < <(find "$ENS" -type f -name '*.nc' | sort)
echo "found ${#FILES[@]} ENSEMBLE files"; echo ""

for f in "${FILES[@]}"; do
  names=$(cdo -s showname "$f" 2>/dev/null | tr ' ' '\n')
  # build chname argument only for pairs whose OLD name exists and NEW does not
  args=""
  for p in "${PAIRS[@]}"; do
    old="${p%% *}"; new="${p##* }"
    if grep -qx "$old" <<<"$names" && ! grep -qx "$new" <<<"$names"; then
      args+="${old},${new},"
    fi
  done
  if [ -z "$args" ]; then
    echo "skip  $(basename "$f")  (nothing to rename)"
    continue
  fi
  args="${args%,}"   # strip trailing comma
  tmp="${f%.nc}.renaming.nc"
  if cdo -s -O chname,"$args" "$f" "$tmp" 2>/dev/null; then
    mv "$tmp" "$f"
    echo "DONE  $(basename "$f")  [$args]"
  else
    rm -f "$tmp"
    echo "FAIL  $(basename "$f")  (chname error — left untouched)"
  fi
done
echo ""
echo "[done] Verify with: cdo -s showname <a file> | tr ' ' '\\n' | grep tw_peak_ge"
