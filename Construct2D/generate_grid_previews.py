#!/usr/bin/env python3
"""Regenerate the grid preview images used in README.md, using PostPycess's
own read_grid()/plot_grid() (see postpycess.py) rather than reimplementing
Plot3D parsing or grid plotting.

Usage: python3 generate_grid_previews.py
Requires: numpy, matplotlib (headless/Agg backend, no display needed).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from postpycess import read_grid, plot_grid

OUT_DIR = 'doc/images'


def save_preview(p3d_file, out_name, title, xlim=None, ylim=None, color='#2a6f97'):
    imax, jmax, kmax, x, y, threed = read_grid(p3d_file)
    plot_grid(x, y, plaincolor=color)
    fig = plt.gcf()
    ax = fig.gca()
    ax.set_title(title)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    fig.set_size_inches(6, 5)
    fig.savefig(f'{OUT_DIR}/{out_name}', dpi=150, bbox_inches='tight',
                facecolor='white')
    plt.close(fig)
    print(f'Wrote {OUT_DIR}/{out_name}')


if __name__ == '__main__':
    import os
    os.makedirs(OUT_DIR, exist_ok=True)

    save_preview('naca0012.p3d', 'ogrd_full.png',
                  'O-grid: full domain (NACA 0012)')
    save_preview('naca0012.p3d', 'ogrd_closeup.png',
                  'O-grid: near-wall clustering (NACA 0012)',
                  xlim=(-0.15, 1.15), ylim=(-0.35, 0.35))
    save_preview('naca0012_sharp.p3d', 'cgrd_full.png',
                  'C-grid: full domain (NACA 0012, sharp TE)',
                  color='#bc4749')
    save_preview('naca0012_sharp.p3d', 'cgrd_closeup.png',
                  'C-grid: near-wall clustering (NACA 0012, sharp TE)',
                  xlim=(-0.15, 1.15), ylim=(-0.35, 0.35), color='#bc4749')
