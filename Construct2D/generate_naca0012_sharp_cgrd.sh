#!/bin/bash
# Example: drive Construct2D non-interactively to regenerate
# naca0012_sharp.p3d/naca0012_sharp.nmf/naca0012_sharp_stats.p3d (C-grid,
# 179x65) from sample_airfoils/naca0012_sharp.dat.
#
# See generate_naca0012_ogrd.sh for how the grid_options.in auto-load and
# GRID/SMTH/QUIT menu sequence works -- same mechanism, C-grid topology.
#
# Usage: ./generate_naca0012_sharp_cgrd.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ ! -x ./construct2d ]; then
    echo "Building construct2d..."
    cp -n Makefile_Linux_MacOSX Makefile
    make
fi

# For a non-OGRD topology, imax = nsrf + 2*nwke (the wake-cut adds nwke
# points on each side); nsrf=129, nwke=25 -> imax=179. jmax=65. These
# reproduce this repo's existing naca0012_sharp.p3d exactly (179x65), the
# dimensions the rest of the pipeline is built/tested around.
cat > grid_options.in <<'EOF'
&SOPT
  nsrf = 129
  radi = 150.0
  nwke = 25
  fdst = 1.0
  fwkl = 1.0
  fwki = 10.0
/
&VOPT
  name = 'naca0012_sharp'
  jmax = 65
  slvr = 'HYPR'
  topo = 'CGRD'
  ypls = 20.0
  recd = 1000000.0
  stp1 = 1000
  stp2 = 20
  nrmt = 1
  nrmb = 1
  alfa = 1.0
  epsi = 15.0
  epse = 0.0
  funi = 0.0
  asmt = 20
/
&OOPT
  gdim = 2
  npln = 2
  dpln = 1.0
/
EOF

printf 'GRID\nSMTH\nQUIT\n' | ./construct2d sample_airfoils/naca0012_sharp.dat

rm -f grid_options.in

echo
echo "Wrote naca0012_sharp.p3d, naca0012_sharp.nmf, naca0012_sharp_stats.p3d in $SCRIPT_DIR"
