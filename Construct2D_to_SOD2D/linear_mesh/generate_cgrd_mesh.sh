#!/bin/bash
# Runs steps 1-3 of ../README.md's pipeline end to end, on the bundled
# C-grid example mesh (examples/naca0012_sharp.p3d/.nmf -- a C-grid wraps
# the airfoil and peels off downstream into a wake region, shaped like
# the letter C; see the top-level README's "big picture" section if
# that's unfamiliar). Good for checking your environment/build is set up
# correctly before running on your own airfoil.
#
# flow_config_cgrd.json's settings (idim=179, jmax=65, porder=4) match
# this specific example mesh's own dimensions -- if you copy this script
# to run a different mesh, you'll likely need a different config too.
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

# run_pipeline.py writes all of its output (extruded grid, Gmsh mesh,
# partitioned SOD2D mesh) into one self-contained work directory, so
# different mesh variants never collide. It also expects its two input
# files to already be sitting there, named exactly {airfoil-file}.p3d/
# .nmf -- hence the copy step below.
WORK_DIR="runs/env_0_cgrd"

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
cp -v examples/naca0012_sharp.p3d "$WORK_DIR/airfoil.p3d"
cp -v examples/naca0012_sharp.nmf "$WORK_DIR/airfoil.nmf"

python3 run_pipeline.py \
    --work-dir "$WORK_DIR" \
    --config flow_config_cgrd.json \
    --airfoil-file airfoil
