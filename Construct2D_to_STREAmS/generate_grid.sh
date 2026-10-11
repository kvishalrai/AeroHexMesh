#!/bin/bash
# Construct2D -> STREAmS-2 grid in one step.
#
# Usage: ./generate_grid.sh <case_dir> <airfoil.dat>
#
# <case_dir> must contain grid_options.in (Construct2D's batch options, loaded automatically from
# the working directory). Runs Construct2D (GRID, SMTH, QUIT) there on <airfoil.dat> (path relative
# to <case_dir> or absolute), which writes <name>.p3d/.nmf/_stats.p3d with <name> from grid_options.in,
# then converts <name>.p3d to <case_dir>/grid2d.xyz with p3d_to_streams.py.
#
# Example: ./generate_grid.sh examples/naca0012_re23k naca0012.dat
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
C2D_DIR="$SCRIPT_DIR/../Construct2D"
CASE_DIR="$(cd "$1" && pwd)"
AIRFOIL="$2"
if [ -z "$AIRFOIL" ] || [ ! -f "$CASE_DIR/grid_options.in" ]; then
    echo "usage: $0 <case_dir with grid_options.in> <airfoil.dat>"; exit 1
fi

# (re)build construct2d from scratch if it is missing or older than any source file
# (its Makefile does not track module dependencies, so incremental builds can fail)
if [ ! -x "$C2D_DIR/construct2d" ] || [ -n "$(find "$C2D_DIR/src" -newer "$C2D_DIR/construct2d")" ]; then
    echo "Building construct2d..."
    (cd "$C2D_DIR" && cp -n Makefile_Linux_MacOSX Makefile && rm -f *.o *.mod construct2d && make >build.log 2>&1)         || { echo "construct2d build failed, see $C2D_DIR/build.log"; exit 1; }
fi

cd "$CASE_DIR"
NAME=$(sed -n "s/^ *name *= *'\([^']*\)'.*/\1/p" grid_options.in)
printf 'GRID\nSMTH\nQUIT\n' | "$C2D_DIR/construct2d" "$AIRFOIL" > construct2d.log
grep -A12 "Grid quality information" construct2d.log || true
python3 "$SCRIPT_DIR/p3d_to_streams.py" "$NAME.p3d" grid2d.xyz
