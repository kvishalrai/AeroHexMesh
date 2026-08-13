#!/bin/bash
# Build a 3D Nek5000 mesh (periodic in the spanwise/Z direction) from a 2D
# .re2 mesh (e.g. smooth_2D/naca_gen_spline_info.re2), via Nek5000's
# re2torea/n2to3 tools, then report -- and optionally patch -- the
# example_nekrs-style case parameters that go with the mesh just built.
# Everything in between is a throwaway intermediate (re2torea needs a
# companion .rea -- base.rea/tail.rea provide its header/footer, see the
# comment below -- and n2to3 itself works in the old ASCII .rea format);
# only the final 3D mesh gets written where you asked for it.
#
# NUMBER_ELEMENTS_Z isn't derivable from the 2D mesh (it's a free spanwise-
# resolution choice, not topology) but is set ONCE here instead of
# separately in the n2to3 call and in whatever case files consume the
# result.
#
# Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]
#   source_2d.re2  the 2D linear mesh to extrude (gmsh2nek's output, e.g.
#                  smooth_2D/naca_gen_spline_info.re2)
#   nmf_file       native 2D Construct2D .nmf the mesh came from (only used
#                  to report/patch case parameters, not the extrusion itself)
#   output_re2     where to write the resulting 3D mesh
#   nlevels_z      number of spanwise elements for n2to3 (default 3)
#   step2_usr      optional example_nekrs-style .usr to patch NUMBER_ELEMENTS_*
#   step2_par      optional example_nekrs-style .par to patch number_elements_*
set -e

SRC_RE2="${1:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
NMF_FILE="${2:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
OUT_RE2="${3:?Usage: ./arglist.sh <source_2d.re2> <nmf_file> <output_re2> [nlevels_z] [step2_usr] [step2_par]}"
NLEVELS_Z="${4:-3}"
STEP2_USR="$5"
STEP2_PAR="$6"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_RE2="$(cd "$(dirname "$SRC_RE2")" && pwd)/$(basename "$SRC_RE2")"
NMF_FILE="$(cd "$(dirname "$NMF_FILE")" && pwd)/$(basename "$NMF_FILE")"
OUT_DIR="$(cd "$(dirname "$OUT_RE2")" && pwd)"
OUT_NAME="$(basename "$OUT_RE2")"
# Resolve these now, before cd-ing into the temp work dir below --
# compute_case_params.py runs from there, not from wherever this script
# was invoked.
[ -n "$STEP2_USR" ] && STEP2_USR="$(cd "$(dirname "$STEP2_USR")" && pwd)/$(basename "$STEP2_USR")"
[ -n "$STEP2_PAR" ] && STEP2_PAR="$(cd "$(dirname "$STEP2_PAR")" && pwd)/$(basename "$STEP2_PAR")"

for tool in re2torea reatore2 n2to3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not found on PATH -- build it via Nek5000's tools/maketools" >&2
        exit 1
    fi
done

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# --- 2D .re2 -> ASCII .rea, via re2torea ---------------------------------
# re2torea needs a full, valid companion .rea alongside the .re2: it reads
# the PARAMETERS/LOGICAL header from it, then rewinds and copies everything
# from PRESOLVE onward as the footer, with the real geometry coming from
# the binary .re2 in between. base.rea/tail.rea already provide exactly
# that header/footer; whatever sits between them here is only scanned past
# and discarded (re2torea's own scanout), so a placeholder line is fine.
printf 'placeholder mesh data - discarded by re2torea scanout\n' > _ph.txt
cat "$SCRIPT_DIR/base.rea" _ph.txt "$SCRIPT_DIR/tail.rea" > bootstrap.rea
cp "$SRC_RE2" bootstrap.re2

re2torea <<EOF
bootstrap
2d
EOF

# --- 2D extrusion to 3D (periodic Z), via n2to3 --------------------------
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

# --- 3D ASCII .rea -> final .re2, via reatore2 ---------------------------
reatore2 <<EOF
3d
final
EOF

cp final.re2 "$OUT_DIR/$OUT_NAME"
echo "Wrote $OUT_DIR/$OUT_NAME"

# --- report / patch case parameters --------------------------------------
PATCH_ARGS=()
[ -n "$STEP2_USR" ] && PATCH_ARGS+=(--patch-step2-usr "$STEP2_USR")
[ -n "$STEP2_PAR" ] && PATCH_ARGS+=(--patch-step2-par "$STEP2_PAR")
python3 "$SCRIPT_DIR/compute_case_params.py" "$NMF_FILE" --nz "$NLEVELS_Z" "${PATCH_ARGS[@]}"
