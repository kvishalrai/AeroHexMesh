#!/usr/bin/env python3
"""
Convert a Construct2D 2D C-grid (.p3d, ASCII plot3d) to a STREAmS-2 grid2d.xyz.

  python3 p3d_to_streams.py naca0012_re23k.p3d grid2d.xyz [--ref old_grid2d.xyz]

grid2d.xyz is what STREAmS-2 reads for grid_type = 6 (curvilinear, grid_dim = 2): a stream
binary with int32 nx, ny, 1 followed by float64 x(i,j) and y(i,j), i fastest.
The grid is oriented as in STREAmS's airfoil example: i runs lower wake -> pressure side ->
leading edge -> suction side -> upper wake, j runs wall -> far field. The two copies of the
wake cut (i and nx+1-i on j = 1) are made bit-identical. Wall spacing, surface spacing,
far-field radius and stretching ratios are printed for checking.
"""
import argparse
import numpy as np

def read_p3d(fname):
    with open(fname) as f:
        nx, ny = map(int, f.readline().split())
        v = np.array(f.read().split(), dtype=float)
    x = v[:nx * ny].reshape(ny, nx).T
    y = v[nx * ny:2 * nx * ny].reshape(ny, nx).T
    return x, y

def read_xyz(fname):
    with open(fname, 'rb') as f:
        nx, ny, _ = np.fromfile(f, np.int32, 3)
        r = np.fromfile(f, np.float64, 2 * nx * ny)
    return r[:nx * ny].reshape(ny, nx).T, r[nx * ny:].reshape(ny, nx).T

def info(x, y, tag):
    nx, ny = x.shape
    i_le = np.argmin(x[:, 0])
    jac = (np.gradient(x, axis=0) * np.gradient(y, axis=1)
           - np.gradient(x, axis=1) * np.gradient(y, axis=0))
    print(f'{tag}: {nx}x{ny}, i=1 at ({x[0,0]:.3f},{y[0,0]:.3f}), LE i={i_le+1}, '
          f'y[i=nx/4,j=1]={y[nx//4,0]:+.4f}, jac sign {np.sign(jac.min())}/{np.sign(jac.max())}')

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('p3d', help='Construct2D grid (.p3d)')
    ap.add_argument('out', help='STREAmS grid file, normally grid2d.xyz')
    ap.add_argument('--ref', help='existing grid2d.xyz to print for comparison')
    args = ap.parse_args()

    x, y = read_p3d(args.p3d)
    nx, ny = x.shape
    if args.ref:
        info(*read_xyz(args.ref), 'reference')
    info(x, y, 'input')
    # wake cut length: points where i and nx+1-i coincide on j=1
    nw = 0
    while np.hypot(x[nw, 0] - x[nx-1-nw, 0], y[nw, 0] - y[nx-1-nw, 0]) < 1e-9:
        nw += 1
    ite = nw                                # 1-based index of the TE node
    # lower surface must come first (pressure side for aoa > 0)
    if y[ite + (nx - 2*ite)//4, 0] > 0:
        x, y = x[::-1], y[::-1]
    # j=1 at wall already; enforce exact symmetry of the wake cut
    for i in range(nw):
        k = nx - 1 - i
        xm, ym = 0.5*(x[i, 0] + x[k, 0]), 0.5*(y[i, 0] + y[k, 0])
        x[i, 0] = x[k, 0] = xm; y[i, 0] = y[k, 0] = ym
    info(x, y, 'output')
    xs, ys = x[ite-1:nx-ite+1, 0], y[ite-1:nx-ite+1, 0]
    print(f'TE node i={ite} at ({x[ite-1,0]:.12f},{y[ite-1,0]:.3e}); wall pts {nx-2*ite+2}; '
          f'chord x [{xs.min():.3e}, {xs.max():.6f}], thickness {ys.max()-ys.min():.5f}')
    dn = np.hypot(x[:, 1]-x[:, 0], y[:, 1]-y[:, 0])[ite-1:nx-ite+1]
    ds = np.hypot(np.diff(xs), np.diff(ys))
    r = np.hypot(x[:, -1]-0.5, y[:, -1])
    print(f'wall-normal dn1 min/max {dn.min():.2e}/{dn.max():.2e}; surface ds min/max '
          f'{ds.min():.2e}/{ds.max():.2e}; farfield r {r.min():.1f}-{r.max():.1f}')
    gi = np.hypot(np.diff(x, axis=0), np.diff(y, axis=0)); gj = np.hypot(np.diff(x, axis=1), np.diff(y, axis=1))
    print(f'max stretching ratio i {np.max(gi[1:]/gi[:-1]):.3f}, j {np.max(gj[:,1:]/gj[:,:-1]):.3f}')
    with open(args.out, 'wb') as f:
        np.array([nx, ny, 1], np.int32).tofile(f)
        x.T.astype(np.float64).tofile(f); y.T.astype(np.float64).tofile(f)
    print('wrote', args.out)

if __name__ == '__main__':
    main()
