#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p /tmp/mpl_cache rendered
export MPLCONFIGDIR=/tmp/mpl_cache
FN="${FIGNIMBLE_BIN:-$HOME/.local/bin/fignimble}"

echo "=== Running Fignimble Nature preflight checks ==="
"$FN" check graphical_abstract.py fig1.py fig2.py fig3.py fig4.py || true

echo "=== Exporting vector PDFs and inspection PNGs ==="
for f in graphical_abstract fig1 fig2 fig3 fig4; do
  "$FN" export "${f}.py" "rendered/${f}.pdf"
  "$FN" export "${f}.py" "rendered/${f}.png" --dpi 220
done
ls -lh rendered/
