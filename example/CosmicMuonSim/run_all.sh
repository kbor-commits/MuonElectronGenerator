#!/usr/bin/env bash
# Build CosmicMuonSim, generate muons with each recipe, and plot the distributions.
# Usage: ./run_all.sh [nEvents]          (default 1000000)
# Geant4 must be set up first, e.g.:  source <geant4-install>/bin/geant4.sh
set -euo pipefail
N=${1:-1000000}
HERE=$(cd "$(dirname "$0")" && pwd)
command -v geant4-config >/dev/null || { echo "Source geant4.sh first"; exit 1; }

mkdir -p "$HERE/build" && cd "$HERE/build"
cmake "$HERE" -DCMAKE_BUILD_TYPE=Release > cmake.log
make -j4 > make.log

mkdir -p "$HERE/output" && cd "$HERE/output"
for r in Reyna Blanpied Chris Jokisch Haruo; do
  echo "== $r: $N events"
  "$HERE/build/cosmic_muons" "$r" "$N" "muons_$r" > "log_$r.txt" 2>&1
done
python3 "$HERE/plot_muons.py" "$HERE/output" "$HERE/plots" Reyna
