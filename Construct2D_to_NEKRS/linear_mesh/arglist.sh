#!/bin/bash
# This script is STEP 3 of the pipeline in ../README.md: it takes the
# jagged-walled 2D mesh from step 1 and turns it into a 3D mesh, by
# repeating ("extruding") the 2D cross-section along a straight line and
# telling the solver the flow repeats identically at every cross-section
# (spanwise-periodic). It also fills in the matching mesh-size settings in
# the NekRS case that will run this mesh (step 4).
#
# Why three separate tool calls just to build one mesh? Nek5000 ships a
# handful of small, old command-line programs, each one only able to
# read/write ONE specific file format:
#   - .re2  = the modern, compact, BINARY mesh format
#   - .rea  = the old, human-readable, PLAIN-TEXT mesh format
# re2torea converts binary->text, n2to3 (the only one of the three that
# actually knows how to extrude a mesh to 3D) only understands the text
# format, and reatore2 converts text->binary again at the end. So the
# path is: binary -> text -> (extrude) -> text -> binary. None of the
# files in between (bootstrap.rea, 2d.rea, 3d.rea, ...) matter on their
# own -- they're thrown away as soon as this script finishes; only the
# final 3D mesh gets written where you asked for it.
#
# NUMBER_ELEMENTS_Z (how many mesh elements span the wing) can't be read
# from the 2D mesh -- it's a free choice of spanwise resolution, not a
# property of the airfoil shape. You give it here, ONCE, as an argument,
# and this script passes the same value both to the extrusion tool AND
# into the NekRS case's settings, so the two can never quietly disagree.
#
# Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]
#   source_2d.re2  the 2D mesh to extrude (gmsh2nek's output from step 1,
#                  e.g. smooth_2D/naca_gen_spline_info.re2)
#   nmf_file       the airfoil's .nmf file (only used to look up/patch
#                  case settings below, not for the extrusion itself)
#   output_re2     where to write the resulting 3D mesh
#   nlevels_z      how many mesh elements along the span (default 3)
#   step2_usr      optional: step 4's .usr file, to fill in NUMBER_ELEMENTS_*
#   step2_par      optional: step 4's .par file, to fill in number_elements_*
set -e

SRC_RE2="${1:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
NMF_FILE="${2:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
OUT_RE2="${3:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
NLEVELS_Z="${4:-3}"
STEP2_USR="$5"
STEP2_PAR="$6"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Turn every path argument into an absolute path right away. Below, this
# script cd's into a scratch directory to do its work; if we waited and
# resolved these paths later, a relative path typed on the command line
# would be looked up relative to the wrong directory.
SRC_RE2="$(cd "$(dirname "$SRC_RE2")" && pwd)/$(basename "$SRC_RE2")"
NMF_FILE="$(cd "$(dirname "$NMF_FILE")" && pwd)/$(basename "$NMF_FILE")"
OUT_DIR="$(cd "$(dirname "$OUT_RE2")" && pwd)"
OUT_NAME="$(basename "$OUT_RE2")"
[ -n "$STEP2_USR" ] && STEP2_USR="$(cd "$(dirname "$STEP2_USR")" && pwd)/$(basename "$STEP2_USR")"
[ -n "$STEP2_PAR" ] && STEP2_PAR="$(cd "$(dirname "$STEP2_PAR")" && pwd)/$(basename "$STEP2_PAR")"

for tool in re2torea reatore2 n2to3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not found on PATH -- build it via Nek5000's tools/maketools" >&2
        exit 1
    fi
done

# Do all the work in a throwaway scratch directory, so this script never
# litters the case folder with intermediate files -- only the final mesh
# (copied out at the very end) and whatever case files we patch survive.
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# ==========================================================================
# Part A: convert the 2D mesh from binary (.re2) to plain-text (.rea) form,
# because that's the only form the next tool (n2to3) can read.
# ==========================================================================
#
# The conversion tool, re2torea, needs a full, valid plain-text .rea file
# alongside the binary .re2 -- not because it reads the mesh from that
# text file (it doesn't: the real geometry always comes from the binary
# .re2), but because a .rea file also carries some bookkeeping sections
# (run parameters, boundary-condition-type flags) that only exist in text
# form. base.rea and tail.rea (in this same directory) are exactly that
# bookkeeping, split into a "header" and "footer" with a gap in between
# for the actual mesh data. Since re2torea gets the real mesh from the
# binary file anyway, whatever we put in that gap is irrelevant -- it
# scans straight past it and throws it away -- so a single placeholder
# line is enough.
printf 'placeholder mesh data - discarded by re2torea scanout\n' > _ph.txt
cat "$SCRIPT_DIR/base.rea" _ph.txt "$SCRIPT_DIR/tail.rea" > bootstrap.rea
cp "$SRC_RE2" bootstrap.re2

re2torea <<EOF
bootstrap
2d
EOF

# ==========================================================================
# Part B: extrude the 2D mesh into a 3D mesh, periodic along the span.
# ==========================================================================
#
# n2to3 is interactive -- each line below answers one of its prompts, in
# order. What each answer means:
n2to3 <<EOF
2d
3d
0
$NLEVELS_Z
0
4.0
1
no
P
EOF
#   2d          input mesh (the file written by Part A, without ".rea")
#   3d          output mesh basename
#   0           output format: 0 = plain text (we still need one more
#               conversion step below, so stay in text form for now)
#   $NLEVELS_Z  how many mesh elements along the span
#   0           span starts at z = 0
#   4.0         span ends at z = 4.0 (in the mesh's own length units)
#   1           uniform spacing between spanwise elements
#   no          not a CEM (electromagnetics) case -- always "no" here
#   P           periodic boundary condition along the span: the flow at
#               one end of the span is forced to exactly match the flow
#               at the other end, which is what makes "repeat the 2D
#               cross-section forever" a valid simplification in the
#               first place

# ==========================================================================
# Part C: convert the finished 3D mesh from text back to Nek5000's normal
# binary form, since that's what NekRS actually expects to load.
# ==========================================================================
reatore2 <<EOF
3d
final
EOF

cp final.re2 "$OUT_DIR/$OUT_NAME"
echo "Wrote $OUT_DIR/$OUT_NAME"

# ==========================================================================
# Finally, make sure the NekRS case that will use this mesh has matching
# settings (element counts) -- see compute_case_params.py for how these
# numbers get worked out and why they can't just be guessed.
# ==========================================================================
PATCH_ARGS=()
[ -n "$STEP2_USR" ] && PATCH_ARGS+=(--patch-step2-usr "$STEP2_USR")
[ -n "$STEP2_PAR" ] && PATCH_ARGS+=(--patch-step2-par "$STEP2_PAR")
python3 "$SCRIPT_DIR/compute_case_params.py" "$NMF_FILE" --nz "$NLEVELS_Z" "${PATCH_ARGS[@]}"
