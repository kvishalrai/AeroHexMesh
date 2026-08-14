#!/bin/bash
# Runs the pyHyp extrusion end to end on the bundled ONERA M6 wing example
# (m6_small.fmt, a Plot3D surface mesh) with ../generate_volume_mesh.py. Good
# for checking your pyHyp build is set up correctly before pointing it at
# your own surface mesh.
#
# Usage: ./generate_mesh.sh
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Lmod's init script (sourced below) turns on `set -e` as a side effect, and
# `module load`'s own exit status is unreliable even on success - so bracket
# both in an explicit +e/-e regardless of what state we're left in.
set +e
source /etc/profile.d/modules.sh 2>/dev/null
eval "${AEROHEXMESH_PYHYP_MODULE_SETUP:-module load anaconda/2023.07 && source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate pyhyp-env}"
set -e

python3 ../generate_volume_mesh.py \
    --config config.json \
    --input-file m6_small.fmt \
    --file-type PLOT3D \
    --output-file volumeMesh.xyz \
    |& tee logMeshGeneration.txt
