#!/bin/bash
# Regenerates the C-grid linear mesh from the examples/ fixture (naca0012_sharp.p3d/.nmf)
# through run_pipeline.py, using flow_config_cgrd.json (currently idim=179, jmax=65,
# porder=4, matching the examples/ fixture's dimensions).
#
# Usage: ./generate_cgrd_mesh.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Lmod's init script (sourced below) turns on `set -e` as a side effect, and
# `module load`'s own exit status is unreliable even on success - so bracket
# both in an explicit +e/-e regardless of what state we're left in.
set +e
source /etc/profile.d/modules.sh 2>/dev/null
module load anaconda >/dev/null 2>&1
set -e

WORK_DIR="runs/env_0_cgrd"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
cp -v examples/naca0012_sharp.p3d "$WORK_DIR/airfoil.p3d"
cp -v examples/naca0012_sharp.nmf "$WORK_DIR/airfoil.nmf"

python3 run_pipeline.py \
    --work-dir "$WORK_DIR" \
    --config flow_config_cgrd.json \
    --airfoil-file airfoil
