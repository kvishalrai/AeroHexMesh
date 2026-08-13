#!/usr/bin/env python3
"""
Compute the ring-width / wall-range Nek5000 case parameters directly from
a native 2D Construct2D .nmf file -- works for both O-grid and C-grid
without needing a --mesh-type flag, since the wall's element range is
already exactly what the NMF's VISCOUS entry encodes (f1==3, excluding the
C-grid's ONE_TO_ONE wake-cut-fold pairing of the same face). No need to run
the mesh through gmsh/gmsh2nek and inspect the output just to find these
numbers.

Optionally patches the values directly into case files -- but the two case
*types* need genuinely different things patched, not just the same values
in different files, so they're separate flags rather than one generic
"patch this .usr" option:

  --patch-step1-usr FILE   a smooth_2D-style .usr (Fischer's spline-fit
                            case). Patches NUMBER_ELEMENTS_X (ring width),
                            WALL_START_ELEMENT, WALL_ELEMENT_COUNT. Its
                            own NUMBER_ELEMENTS_Y means something entirely
                            different here (how many near-wall RINGS to
                            actually spline-correct, not the mesh's total
                            radial ring count) -- pass --smoothing-levels
                            explicitly to also patch that; it is NOT
                            derived from the mesh and has no safe default,
                            silently defaulting it to the full radial
                            count would smooth every ring instead of just
                            the near-wall one(s).

  --patch-step2-usr FILE   an example_nekrs-style .usr (applies the
                            precomputed spline correction to the real 3D
                            mesh). Patches NUMBER_ELEMENTS_X and
                            NUMBER_ELEMENTS_Y (here genuinely the mesh's
                            full radial ring count, auto-derived). Has no
                            WALL_START_ELEMENT/WALL_ELEMENT_COUNT to patch
                            (get_smooth_data works on every element, wall
                            or not). Also patches NUMBER_ELEMENTS_Z if
                            --nz is given.

  --patch-step2-par FILE   the matching .par's [CASEDATA] section --
                            same values as --patch-step2-usr.

Usage:
    python3 compute_case_params.py mesh.nmf
    python3 compute_case_params.py mesh.nmf --smoothing-levels 1 \\
        --patch-step1-usr smooth_2D/naca_gen_spline_info.usr
    python3 compute_case_params.py mesh.nmf --nz 3 \\
        --patch-step2-usr example_nekrs/naca.usr \\
        --patch-step2-par example_nekrs/naca.par
"""

import argparse
import re
import sys

from p3d_to_gmsh_nek import read_2d_nmf


def compute_wall_params(nmf_path):
    idim, jmax, boundaries = read_2d_nmf(nmf_path)

    wall = [
        b for b in boundaries
        if b['f1'] == 3 and b['name'].upper() not in ('ONE_TO_ONE', 'ONE-TO-ONE')
    ]
    if len(wall) != 1:
        raise ValueError(
            f"Expected exactly one non-ONE_TO_ONE f1==3 (VISCOUS) boundary "
            f"in {nmf_path}, found {len(wall)}"
        )
    wb = wall[0]

    return {
        'ring_width': idim - 1,               # NUMBER_ELEMENTS_X in both usr styles
        'radial_rings': jmax - 1,              # NUMBER_ELEMENTS_Y, step2 meaning only
        'wall_start_element': wb['s1'],        # 1-based, matches Fortran element numbering
        'wall_element_count': wb['e1'] - wb['s1'],
    }


def _patch_defines(path, replacements):
    """replacements: dict of {DEFINE_NAME: value}. Matches the single
    token after '#define NAME', whether it's currently a literal or a
    symbolic reference (e.g. smooth_2D's default
    '#define WALL_ELEMENT_COUNT  NUMBER_ELEMENTS_X')."""
    with open(path) as f:
        text = f.read()
    for name, value in replacements.items():
        pattern = re.compile(rf'(#define\s+{name}\s+)\S+')
        text, n = pattern.subn(rf'\g<1>{value}', text)
        if n == 0:
            print(f"  WARNING: no '#define {name}' in {path}, not patched", file=sys.stderr)
    with open(path, 'w') as f:
        f.write(text)
    print(f"Patched {path}: " + ", ".join(f"{k}={v}" for k, v in replacements.items()))


def _patch_par(path, replacements):
    """replacements: dict of {par_key: value}, e.g. {'number_elements_x': 230}."""
    with open(path) as f:
        text = f.read()
    for key, value in replacements.items():
        pattern = re.compile(rf'(?mi)^({re.escape(key)}\s*=\s*)\S+')
        text, n = pattern.subn(rf'\g<1>{value}', text)
        if n == 0:
            print(f"  WARNING: no '{key} =' in {path}, not patched", file=sys.stderr)
    with open(path, 'w') as f:
        f.write(text)
    print(f"Patched {path}: " + ", ".join(f"{k}={v}" for k, v in replacements.items()))


def main():
    parser = argparse.ArgumentParser(
        description="Compute (and optionally patch in) the ring/wall element "
                    "parameters for a Nek5000/NekRS airfoil case, directly "
                    "from the native 2D Construct2D .nmf file (O-grid or C-grid).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("nmf_file", help="Native 2D Neutral Map File (.nmf)")
    parser.add_argument(
        "--nz", type=int,
        help="NUMBER_ELEMENTS_Z (spanwise element count) for --patch-step2-* -- "
             "not derivable from the .nmf, it's your own n2to3 extrusion choice",
    )
    parser.add_argument(
        "--smoothing-levels", type=int,
        help="NUMBER_ELEMENTS_Y for --patch-step1-usr -- how many near-wall "
             "rings smooth_geom should spline-correct. Not derivable from the "
             "mesh; typically 1 (just the wall-adjacent ring).",
    )
    parser.add_argument("--patch-step1-usr", nargs='+', default=[], metavar="FILE",
                         help="smooth_2D-style .usr file(s) (Fischer's spline-fit case)")
    parser.add_argument("--patch-step2-usr", nargs='+', default=[], metavar="FILE",
                         help="example_nekrs-style .usr file(s) (applies the correction to the 3D mesh)")
    parser.add_argument("--patch-step2-par", nargs='+', default=[], metavar="FILE",
                         help="example_nekrs-style .par file(s) ([CASEDATA] section)")
    args = parser.parse_args()

    p = compute_wall_params(args.nmf_file)

    print(f"Computed from {args.nmf_file}:")
    print(f"  ring width (NUMBER_ELEMENTS_X)      = {p['ring_width']}")
    print(f"  radial rings (step2 NUMBER_ELEMENTS_Y) = {p['radial_rings']}")
    print(f"  WALL_START_ELEMENT                  = {p['wall_start_element']}")
    print(f"  WALL_ELEMENT_COUNT                  = {p['wall_element_count']}")
    if p['wall_element_count'] == p['ring_width'] and p['wall_start_element'] == 1:
        print("  -> wall spans the whole ring: this is an O-grid-shaped mesh")
    else:
        print("  -> wall is a sub-range of the ring: this is a C-grid-shaped mesh")

    for usr in args.patch_step1_usr:
        step1 = {
            'NUMBER_ELEMENTS_X': p['ring_width'],
            'WALL_START_ELEMENT': p['wall_start_element'],
            'WALL_ELEMENT_COUNT': p['wall_element_count'],
        }
        if args.smoothing_levels is not None:
            step1['NUMBER_ELEMENTS_Y'] = args.smoothing_levels
        else:
            print(
                f"  NOTE: --smoothing-levels not given, leaving {usr}'s "
                f"NUMBER_ELEMENTS_Y (near-wall rings to smooth) untouched",
                file=sys.stderr,
            )
        _patch_defines(usr, step1)

    for usr in args.patch_step2_usr:
        step2 = {
            'NUMBER_ELEMENTS_X': p['ring_width'],
            'NUMBER_ELEMENTS_Y': p['radial_rings'],
        }
        if args.nz is not None:
            step2['NUMBER_ELEMENTS_Z'] = args.nz
        _patch_defines(usr, step2)

    for par in args.patch_step2_par:
        step2_par = {
            'number_elements_x': p['ring_width'],
            'number_elements_y': p['radial_rings'],
        }
        if args.nz is not None:
            step2_par['number_elements_z'] = args.nz
        _patch_par(par, step2_par)


if __name__ == "__main__":
    main()
