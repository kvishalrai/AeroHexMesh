#!/bin/bash
# Example: drive Construct2D non-interactively to regenerate naca0012.p3d/
# naca0012.nmf/naca0012_stats.p3d (O-grid, 128x65) from sample_airfoils/naca0012.dat.
#
# Construct2D is an interactive menu program, but it auto-loads settings
# from a fixed-name file called grid_options.in in the working directory if
# one is present (see src/menu.f90:set_defaults) -- that covers every
# surface/volume/output option non-interactively. The only menu interaction
# left is driving the main loop itself: GRID -> SMTH (generate the smoothed
# grid) -> QUIT, piped via stdin.
#
# Usage: ./generate_naca0012_ogrd.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ ! -x ./construct2d ]; then
    echo "Building construct2d..."
    cp -n Makefile_Linux_MacOSX Makefile
    make
fi

# nsrf directly sets the O-grid's i-dimension (topo='OGRD' -> imax=nsrf);
# jmax matches. These reproduce this repo's existing naca0012.p3d exactly
# (128x65), the dimensions the rest of the pipeline is built/tested around.
cat > grid_options.in <<'EOF'
&SOPT
  nsrf = 128
  radi = 15.0
  nwke = 50
  fdst = 1.0
  fwkl = 1.0
  fwki = 10.0
/
&VOPT
  name = 'naca0012'
  jmax = 65
  slvr = 'HYPR'
  topo = 'OGRD'
  ypls = 0.9
  recd = 1000000.0
  stp1 = 1000
  stp2 = 20
  nrmt = 1
  nrmb = 1
  alfa = 1.0
  epsi = 5.0
  epse = 0.0
  funi = 0.2
  asmt = 10
/
&OOPT
  gdim = 2
  npln = 2
  dpln = 1.0
/
EOF

printf 'GRID\nSMTH\nQUIT\n' | ./construct2d sample_airfoils/naca0012.dat

rm -f grid_options.in

echo
echo "Wrote naca0012.p3d, naca0012.nmf, naca0012_stats.p3d in $SCRIPT_DIR"
