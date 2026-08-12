"""Build an interactive before/after visualization of wall-spline smoothing
node movement.

Reads the wall_coords_before/after_rank*.dat diagnostic dump
(MeshElasticitySolver.f90:dump_wall_coords) and the compact wall_spline.dat
corner+coefficient table (wall_spline.py) used for that same run, and
produces a single self-contained HTML page: raw (before) node positions, the
fitted spline curve and its corner points, and the smoothed (after)
positions, with a slider to exaggerate the (otherwise sub-percent-of-chord)
connecting line for visibility, pan/zoom, and hover tooltips.

Usage:
    python3 plot_wall_migration.py \\
        --run-dir /path/to/high_order_mesh_test_new \\
        --spline-table naca_o_wall_spline.dat \\
        --output wall_migration.html \\
        --mesh-label "O-grid mesh (naca0012, idim=128)" \\
        --minq-before 0.025872 --minq-after 0.025866
"""

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np

TEMPLATE_PATH = Path(__file__).parent / "wall_migration_template.html"


def load_wall_coords(run_dir, tag):
    """Read and merge every wall_coords_<tag>_rank*.dat in run_dir, sorted by
    global id. Returns (ids, xyz) as numpy arrays."""
    paths = sorted(glob.glob(str(Path(run_dir) / f"wall_coords_{tag}_rank*.dat")))
    if not paths:
        raise FileNotFoundError(f"No wall_coords_{tag}_rank*.dat found in {run_dir}")

    ids, xyz = [], []
    for path in paths:
        with open(path) as f:
            next(f)  # header line
            for line in f:
                parts = line.split()
                ids.append(int(parts[0]))
                xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])

    ids = np.array(ids)
    xyz = np.array(xyz)
    order = np.argsort(ids)
    return ids[order], xyz[order]


def build_raw(run_dir):
    """Match before/after wall nodes by global id and return the RAW rows
    ([x0,y0,x1,y1,disp,gid]) the template's JS expects, plus the raw
    displacement-magnitude array for reporting."""
    ids_b, xyz_b = load_wall_coords(run_dir, "before")
    ids_a, xyz_a = load_wall_coords(run_dir, "after")
    if not np.array_equal(ids_b, ids_a):
        raise ValueError(
            "before/after global-id sets differ -- mismatched dump files? "
            f"({len(ids_b)} vs {len(ids_a)} nodes)"
        )

    disp = np.linalg.norm((xyz_a - xyz_b)[:, :2], axis=1)
    raw = [
        [
            float(xyz_b[i, 0]), float(xyz_b[i, 1]),
            float(xyz_a[i, 0]), float(xyz_a[i, 1]),
            float(disp[i]), int(ids_b[i]),
        ]
        for i in range(len(ids_b))
    ]
    return raw, disp


def load_spline_table(path):
    """Parse wall_spline.py's compact corner+coefficient table:
        line 1: n_corners
        n_corners rows: x y s
        (n_corners-1) rows: ax bx cx dx ay by cy dy
    """
    with open(path) as f:
        n = int(f.readline())
        corners = [tuple(map(float, f.readline().split())) for _ in range(n)]
        segs = [list(map(float, f.readline().split())) for _ in range(n - 1)]
    return corners, segs


def sample_curve(corners, segs, samples_per_segment=20):
    """Dense (x,y) sampling of the fitted spline, evaluating each segment's
    own closed-form cubic -- the same evaluation MeshElasticitySolver.f90
    does at runtime, just sampled densely here for a smooth plotted curve."""
    curve = []
    for i, (ax, bx, cx, dx, ay, by, cy, dy) in enumerate(segs):
        s0, s1 = corners[i][2], corners[i + 1][2]
        for t in np.linspace(0.0, s1 - s0, samples_per_segment):
            x = ax + t * (bx + t * (cx + t * dx))
            y = ay + t * (by + t * (cy + t * dy))
            curve.append([float(x), float(y)])
    return curve


def render_html(template_path, raw, curve, corners, title, mesh_label, minq_before, minq_after):
    with open(template_path) as f:
        html = f.read()

    minq_note = ""
    if minq_before is not None and minq_after is not None:
        minq_note = (
            f" minQ held at {minq_before:.6g} before and {minq_after:.6g} "
            f"after (no collapse)."
        )

    html = (
        html
        .replace("__DATA__", json.dumps(raw))
        .replace("__CURVE__", json.dumps(curve))
        .replace("__CORNERS__", json.dumps(corners))
        .replace("__MESH_LABEL__", mesh_label)
        .replace("__MINQ_NOTE__", minq_note)
        .replace("__N_CORNERS__", str(len(corners)))
        .replace("__TITLE__", title)
    )
    return html


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--run-dir", required=True,
                    help="Directory with wall_coords_before/after_rank*.dat")
    p.add_argument("--spline-table", required=True,
                    help="Path to the compact wall_spline.dat table used for that run")
    p.add_argument("--output", required=True, help="Output HTML path")
    p.add_argument("--title", default="Wall Node Migration")
    p.add_argument("--mesh-label", default="airfoil wall mesh",
                    help='e.g. "O-grid mesh (naca0012, idim=128)"')
    p.add_argument("--minq-before", type=float, default=None,
                    help="minQ before elasticity, for the header note (optional)")
    p.add_argument("--minq-after", type=float, default=None,
                    help="minQ after elasticity, for the header note (optional)")
    p.add_argument("--template", default=str(TEMPLATE_PATH))
    args = p.parse_args()

    raw, disp = build_raw(args.run_dir)
    corners_full, segs = load_spline_table(args.spline_table)
    curve = sample_curve(corners_full, segs)
    corner_pts = [[x, y] for x, y, s in corners_full]

    html = render_html(
        args.template, raw, curve, corner_pts,
        args.title, args.mesh_label, args.minq_before, args.minq_after,
    )

    Path(args.output).write_text(html)

    print(f"Wrote {args.output} ({len(html):,} bytes)")
    print(f"  wall nodes:  {len(raw)}")
    print(f"  max |disp|:  {disp.max():.6e}")
    print(f"  mean |disp|: {disp.mean():.6e}")
    print(f"  corners:     {len(corner_pts)}  (curve samples: {len(curve)})")


if __name__ == "__main__":
    main()
