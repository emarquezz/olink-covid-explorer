#!/usr/bin/env bash
# Downloads the four CSV files from the Gisby et al. (2021) Olink COVID-19 study.
# Repo: https://github.com/Palash63/longitudinal_olink_proteomics
# URLs verified via GitHub API on 2026-09-14: branch = main, files in data/.

set -euo pipefail

BASE="https://raw.githubusercontent.com/Palash63/longitudinal_olink_proteomics/main/data"

PLASMA_NPX="${BASE}/plasma_npx_level.csv"
PLASMA_META="${BASE}/plasma_sample_level.csv"
SERUM_NPX="${BASE}/serum_npx_level.csv"
SERUM_META="${BASE}/serum_sample_level.csv"

OUT="data"
mkdir -p "$OUT"

echo "Downloading to $OUT/ ..."
curl -fL "$PLASMA_NPX"  -o "$OUT/plasma_npx_level.csv"
curl -fL "$PLASMA_META" -o "$OUT/plasma_sample_level.csv"
curl -fL "$SERUM_NPX"   -o "$OUT/serum_npx_level.csv"
curl -fL "$SERUM_META"  -o "$OUT/serum_sample_level.csv"

echo "Done. Files:"
ls -lh "$OUT"/*.csv
