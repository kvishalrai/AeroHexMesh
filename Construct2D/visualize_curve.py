#!/usr/bin/env python3

#  InspectCut -- an interactive viewer for the xcut-pipeline's curve .dat
#  files (the single closed curve fed to Construct2D, before meshing).
#  Same interaction style as postpycess.py, this project's own CFD-grid
#  postprocessor: a text menu loop, matplotlib windows that block until
#  closed, then back to the menu.
#
#  Unlike postpycess.py this doesn't need a mesh -- it only needs a
#  plain "header line + x y pairs" .dat curve (step1's own output), or
#  a header-less "x y per line" raw digitized-coordinate file (this
#  project's other established convention); both are auto-detected. If
#  a same-named <file>.meta.json sidecar exists next to it (written by
#  step1's own write_meta), every "i_..." entry in it -- a single int
#  (one point) or a 2-element [start, end] (an inclusive 1-based index
#  range) -- is picked up automatically as a NAMED SEGMENT, with no
#  hardcoded knowledge of which pipeline (single-body xcut, or
#  multi-element branch-cut) wrote it. Those segments drive the
#  colored-overview and zoom-to-segment views, and label points in the
#  click-to-inspect view.
#
#  Usage:
#      python3 visualize_curve.py
#      (then follow the prompts, same as postpycess.py)

import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmobj


################################################################################
#
# Main program
#
################################################################################
def main():

  if sys.version_info[0] < 3:
    print('Error: this tool requires Python 3. Run it with python3 instead.')
    return

  loadnew = True

# Load a new curve every time this loop is entered
  while loadnew:

    loadnew = False

    x, y, header, segments, fname = load_curve()
    if x is None:
      continue

#   Default program options
    markersize = 3
    zoompad = 0.6
    colormap = 'tab20'
    current_options = [markersize, zoompad, colormap]

    menudone = False
    while not menudone:

      choice = main_menu(fname, len(x), segments)

      if choice == '1':
        plot_plain(x, y, fname, markersize)
      elif choice == '2':
        if segments:
          plot_segments(x, y, segments, fname, markersize, colormap)
        else:
          print('No named segments available for this file (no .meta.json '
                'sidecar found next to it, or it had no usable "i_..." '
                'entries).\n')
      elif choice == '3':
        if segments:
          zoom_to_segment(x, y, segments, fname, markersize, zoompad)
        else:
          print('No named segments available to zoom to.\n')
      elif choice == '4':
        click_inspect(x, y, segments, fname, markersize)
      elif choice == '5':
        zoom_to_region(x, y, fname, markersize, zoompad)
      elif choice == 'O' or choice == 'o':
        new_options = change_options(current_options)
        markersize = new_options[0]
        zoompad = new_options[1]
        colormap = new_options[2]
      elif choice == 'L' or choice == 'l':
        menudone = True
        loadnew = True
      elif choice == 'Q' or choice == 'q':
        menudone = True
      else:
        print('Error: option ' + choice + ' not recognized.\n')


################################################################################
#
# Checks if a string can be converted to int
#
################################################################################
def is_int(string):

  try:
    int(string)
    return True
  except ValueError:
    return False


################################################################################
#
# Checks if a string can be converted to float
#
################################################################################
def is_float(string):

  try:
    float(string)
    return True
  except ValueError:
    return False


################################################################################
#
# Prompts for a curve .dat file, reads it, and picks up its .meta.json
# sidecar (if any) as a list of named segments.
#
################################################################################
def load_curve():

  print('\nInspectCut -- interactive viewer for xcut-pipeline curve files')
  fname = input('\nEnter path to a curve .dat file (blank to cancel): ').strip()
  print('')
  if not fname:
    return None, None, None, None, None

  if not os.path.isfile(fname):
    print('Error: can\'t read the file ' + fname + '\n')
    return None, None, None, None, None

  x, y, header = read_curve_dat(fname)
  print('Read ' + str(len(x)) + ' points from ' + fname)
  if header:
    print('Header: ' + header)

  segments = None
  metapath = os.path.splitext(fname)[0] + '.meta.json'
  if os.path.isfile(metapath):
    try:
      with open(metapath) as f:
        meta = json.load(f)
      segments = extract_segments(meta, len(x))
      print('Read ' + str(len(segments)) + ' named segment(s) from ' + metapath)
    except Exception as e:
      print('Warning: could not parse ' + metapath + ' (' + str(e) + ')')
  print('')

  return x, y, header, segments, fname


################################################################################
#
# Reads a curve .dat file. Auto-detects whether the first line is data
# (two numbers -- this project's header-less raw-digitization convention)
# or a header line to skip (this project's step1-output convention).
#
################################################################################
def read_curve_dat(fname):

  with open(fname) as f:
    lines = [l for l in f if l.strip()]

  first_tokens = lines[0].split()
  first_is_data = False
  if len(first_tokens) >= 2 and is_float(first_tokens[0]) and is_float(first_tokens[1]):
    first_is_data = True

  header = None
  start = 0
  if not first_is_data:
    header = lines[0].strip()
    start = 1

  xs = []
  ys = []
  for l in lines[start:]:
    tokens = l.split()
    xs.append(float(tokens[0]))
    ys.append(float(tokens[1]))

  return np.array(xs), np.array(ys), header


################################################################################
#
# Scans a step1-style meta dict for "i_..." entries that look like a
# 1-based point index (a bare int) or an inclusive 1-based index range
# (a 2-element list/tuple), both within [1, npts]. Returns a list of
# dicts {name, i0, i1} sorted by i0, with no hardcoded knowledge of any
# particular pipeline's own field names.
#
################################################################################
def extract_segments(meta, npts):

  segments = []
  for key, val in meta.items():
    if not key.startswith('i_'):
      continue
    if isinstance(val, bool):
      continue
    if isinstance(val, int) and 1 <= val <= npts:
      segments.append(dict(name=key, i0=val, i1=val))
    elif (isinstance(val, (list, tuple)) and len(val) == 2 and
          all(isinstance(v, int) and not isinstance(v, bool) for v in val) and
          1 <= val[0] <= npts and 1 <= val[1] <= npts and val[0] <= val[1]):
      segments.append(dict(name=key, i0=val[0], i1=val[1]))

  segments.sort(key=lambda s: s['i0'])
  return segments


################################################################################
#
# Returns the (comma-joined) name(s) of every segment covering 1-based
# index i1 (a point can legitimately belong to more than one recorded
# segment, e.g. a shared closure/revisit point), or None if uncovered.
#
################################################################################
def segment_names_for_index(segments, i1):

  names = [s['name'] for s in segments if s['i0'] <= i1 <= s['i1']]
  return ', '.join(names) if names else None


################################################################################
#
# Main menu
#
################################################################################
def main_menu(fname, npts, segments):

  print('Curve: ' + os.path.basename(fname) + '  (' + str(npts) + ' points)')
  if segments:
    print(str(len(segments)) + ' named segment(s) available (from .meta.json)')
  else:
    print('No .meta.json sidecar (or no usable segments) -- plain views only')

  print('\nSelect an option:\n')
  print(' 1) Plot the whole curve (plain)')
  print(' 2) Plot the whole curve, colored by named segment')
  print(' 3) Zoom in on one named segment')
  print(' 4) Click-to-inspect (identify point indices interactively)')
  print(' 5) Zoom to a region (by curve position, index range, or x-range)')
  print(' O) Change options')
  print(' L) Load a different curve file')
  print(' Q) Quit')

  choice = input('\nInput: ')
  print('')

  return choice


################################################################################
#
# Change program options
#
################################################################################
def change_options(current_options):

  markersize, zoompad, colormap = current_options
  cmaplist = ['tab20', 'tab10', 'nipy_spectral', 'gist_rainbow', 'jet']

  optdone = False
  while not optdone:

    print('Select option to change:\n')
    print(' 1) Marker size (current = ' + str(markersize) + ')')
    print(' 2) Zoom padding fraction around a segment (current = ' + str(zoompad) + ')')
    print(' 3) Colormap for segment coloring (current = ' + colormap + ')')
    print(' Q) Go back to the main menu')

    optselect = input('\nInput: ')
    print('')

    if optselect == '1':
      newval = input('Enter new marker size (or Q to cancel): ')
      if is_float(newval):
        markersize = float(newval)
      elif newval not in ('Q', 'q'):
        print('Error: not a number.\n')
    elif optselect == '2':
      newval = input('Enter new zoom padding fraction, e.g. 0.6 (or Q to cancel): ')
      if is_float(newval):
        zoompad = float(newval)
      elif newval not in ('Q', 'q'):
        print('Error: not a number.\n')
    elif optselect == '3':
      print('Available colormaps:\n')
      for i, c in enumerate(cmaplist):
        print(' ' + str(i + 1) + ') ' + c)
      sel = input('\nInput: ')
      if is_int(sel) and 1 <= int(sel) <= len(cmaplist):
        colormap = cmaplist[int(sel) - 1]
      elif sel not in ('Q', 'q'):
        print('Error: option ' + sel + ' not recognized.\n')
    elif optselect == 'Q' or optselect == 'q':
      optdone = True
    else:
      print('Error: option ' + optselect + ' not recognized.\n')

  print('')
  return [markersize, zoompad, colormap]


################################################################################
#
# Plain plot of the whole curve
#
################################################################################
def plot_plain(x, y, fname, markersize):

  print('Plotting curve ...\n')

  fig, ax = plt.subplots()
  ax.plot(x, y, '-o', color='black', ms=markersize, lw=1.0)
  ax.set_aspect('equal')
  ax.set_title(os.path.basename(fname) + '  (' + str(len(x)) + ' points)')
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  plt.show()


################################################################################
#
# Plot the whole curve with every named segment drawn in its own color.
# Any point not covered by a segment is drawn in plain gray so the full
# curve stays visible.
#
################################################################################
def plot_segments(x, y, segments, fname, markersize, colormap):

  print('Plotting curve colored by named segment ...\n')

  cmap = cmobj.get_cmap(colormap)
  ncolors = getattr(cmap, 'N', 20)

  fig, ax = plt.subplots()
  ax.plot(x, y, '-', color='0.8', lw=0.8, zorder=1)

  for idx, seg in enumerate(segments):
    i0, i1 = seg['i0'] - 1, seg['i1'] - 1
    color = cmap(idx % ncolors)
    if i0 == i1:
      ax.plot(x[i0], y[i0], 'o', color=color, ms=markersize + 3,
               zorder=3, label=seg['name'])
    else:
      ax.plot(x[i0:i1 + 1], y[i0:i1 + 1], '-o', color=color, ms=markersize,
               lw=1.5, zorder=2, label=seg['name'])

  ax.set_aspect('equal')
  ax.set_title(os.path.basename(fname) + ' -- colored by named segment')
  ax.set_xlabel('x')
  ax.set_ylabel('y')

  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.72, box.height])
  ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1.0), loc=2, borderaxespad=0.)

  plt.show()


################################################################################
#
# List the named segments and zoom in on the one the user picks, with a
# little context (a handful of points on each side) drawn alongside it.
#
################################################################################
def zoom_to_segment(x, y, segments, fname, markersize, zoompad):

  print('Named segments:\n')
  for idx, seg in enumerate(segments):
    rng = str(seg['i0']) if seg['i0'] == seg['i1'] else str(seg['i0']) + '-' + str(seg['i1'])
    print(' ' + str(idx + 1) + ') ' + seg['name'] + '  (i=' + rng + ')')
  print(' Q) cancel')

  sel = input('\nZoom to which segment? ')
  print('')

  if sel == 'Q' or sel == 'q':
    return
  if not is_int(sel) or not (1 <= int(sel) <= len(segments)):
    print('Error: not a valid segment number.\n')
    return

  seg = segments[int(sel) - 1]
  i0, i1 = seg['i0'] - 1, seg['i1'] - 1

  print('Plotting zoomed view of ' + seg['name'] + ' ...\n')
  zoom_to_index_range(x, y, fname, markersize, zoompad, i0, i1, seg['name'], npad=5)


################################################################################
#
# Zoom to different sections of the curve, without needing a .meta.json
# sidecar: a handful of quick presets that scrub along the curve by
# point index (front/quarter/middle/three-quarter/end -- works for any
# curve, since it never assumes where the "real" airfoil surface is
# versus a wake or bridge), or a custom index range or x-coordinate
# range typed in directly.
#
################################################################################
def zoom_to_region(x, y, fname, markersize, zoompad):

  npts = len(x)

  print('Zoom to a region:\n')
  print(' 1) Start of curve (first ~5%)')
  print(' 2) Quarter point')
  print(' 3) Middle of curve')
  print(' 4) Three-quarter point')
  print(' 5) End of curve (last ~5%)')
  print(' 6) Custom index range')
  print(' 7) Custom x-coordinate range')
  print(' Q) cancel')

  sel = input('\nInput: ')
  print('')

  presets = {'1': 0.0, '2': 0.25, '3': 0.5, '4': 0.75, '5': 1.0}

  if sel in presets:
    frac = presets[sel]
    center = int(round(frac * (npts - 1)))
    halfwidth = max(int(round(0.025 * npts)), 5)
    i0 = max(0, center - halfwidth)
    i1 = min(npts - 1, center + halfwidth)
    print('Plotting zoomed view near i=' + str(i0 + 1) + '-' + str(i1 + 1) + ' ...\n')
    zoom_to_index_range(x, y, fname, markersize, zoompad, i0, i1,
                         'i=' + str(i0 + 1) + '-' + str(i1 + 1), npad=5)

  elif sel == '6':
    i0s = input('Start index (1-based): ')
    i1s = input('End index (1-based): ')
    print('')
    if not is_int(i0s) or not is_int(i1s):
      print('Error: indices must be integers.\n')
      return
    i0, i1 = int(i0s), int(i1s)
    if i0 > i1:
      i0, i1 = i1, i0
    if not (1 <= i0 <= npts) or not (1 <= i1 <= npts):
      print('Error: indices must be within 1-' + str(npts) + '.\n')
      return
    print('Plotting zoomed view of i=' + str(i0) + '-' + str(i1) + ' ...\n')
    zoom_to_index_range(x, y, fname, markersize, zoompad, i0 - 1, i1 - 1,
                         'i=' + str(i0) + '-' + str(i1), npad=5)

  elif sel == '7':
    xmins = input('Minimum x: ')
    xmaxs = input('Maximum x: ')
    print('')
    if not is_float(xmins) or not is_float(xmaxs):
      print('Error: x-range must be numeric.\n')
      return
    xmin, xmax = float(xmins), float(xmaxs)
    if xmin > xmax:
      xmin, xmax = xmax, xmin
    mask = (x >= xmin) & (x <= xmax)
    if not mask.any():
      print('Error: no points found with x in [' + str(xmin) + ', ' + str(xmax) + '].\n')
      return
    print('Plotting zoomed view of x in [' + str(xmin) + ', ' + str(xmax) + '] ...\n')
    zoom_to_xrange(x, y, fname, markersize, zoompad, xmin, xmax)

  elif sel == 'Q' or sel == 'q':
    return
  else:
    print('Error: option ' + sel + ' not recognized.\n')


################################################################################
#
# Shared rendering for a zoom on a contiguous 0-based index range [i0,i1],
# with a little context (npad points on each side) drawn alongside it.
#
################################################################################
def zoom_to_index_range(x, y, fname, markersize, zoompad, i0, i1, label, npad=5):

  lo = max(0, i0 - npad)
  hi = min(len(x) - 1, i1 + npad)

  fig, ax = plt.subplots()
  ax.plot(x, y, '-', color='0.8', lw=0.8, zorder=1, label='full curve')
  ax.plot(x[lo:hi + 1], y[lo:hi + 1], '-o', color='0.4', ms=markersize,
           zorder=2, label='context')
  ax.plot(x[i0:i1 + 1], y[i0:i1 + 1], '-o', color='C1', ms=markersize + 2,
           zorder=3, label=label)

  xseg = x[i0:i1 + 1]
  yseg = y[i0:i1 + 1]
  xr = max(xseg.max() - xseg.min(), 1e-9)
  yr = max(yseg.max() - yseg.min(), 1e-9)
  xmid = 0.5 * (xseg.max() + xseg.min())
  ymid = 0.5 * (yseg.max() + yseg.min())
  halfspan = max(xr, yr) * (0.5 + zoompad)
  ax.set_xlim(xmid - halfspan, xmid + halfspan)
  ax.set_ylim(ymid - halfspan, ymid + halfspan)

  ax.set_aspect('equal')
  ax.set_title('Zoom: ' + label + '  -- ' + os.path.basename(fname))
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.legend(fontsize=8)

  plt.show()


################################################################################
#
# Zoom rendering for a (possibly non-contiguous) x-coordinate range --
# e.g. upper and lower surface can both cross the same x window, so
# this doesn't try to reduce it to one index range: it just plots the
# whole curve and frames the view on that x-band (with matching y
# padding), so every point in that band stays visible regardless of
# how many separate arcs of the curve pass through it.
#
################################################################################
def zoom_to_xrange(x, y, fname, markersize, zoompad, xmin, xmax):

  mask = (x >= xmin) & (x <= xmax)
  yseg = y[mask]

  fig, ax = plt.subplots()
  ax.plot(x, y, '-o', color='0.4', ms=markersize, lw=1.0, zorder=1)

  xr = max(xmax - xmin, 1e-9)
  yr = max(yseg.max() - yseg.min(), 1e-9)
  xmid = 0.5 * (xmin + xmax)
  ymid = 0.5 * (yseg.max() + yseg.min())
  halfspan = max(xr, yr) * (0.5 + zoompad)
  ax.set_xlim(xmid - halfspan, xmid + halfspan)
  ax.set_ylim(ymid - halfspan, ymid + halfspan)

  ax.set_aspect('equal')
  label = 'x in [' + '{:.6g}'.format(xmin) + ', ' + '{:.6g}'.format(xmax) + ']'
  ax.set_title('Zoom: ' + label + '  -- ' + os.path.basename(fname))
  ax.set_xlabel('x')
  ax.set_ylabel('y')

  plt.show()


################################################################################
#
# Click-to-inspect: click near any point on the curve to see its 1-based
# index (matching the p3d convention once meshed with NWKE=0), its (x,y)
# coordinate, and which named segment(s), if any, it belongs to.
#
################################################################################
def click_inspect(x, y, segments, fname, markersize):

  print('Click near any point to identify it. Close the window when done.\n')

  fig, ax = plt.subplots()
  ax.plot(x, y, '-o', color='0.4', ms=markersize, lw=1.0, picker=5)
  ax.set_aspect('equal')
  ax.set_title('Click a point to inspect it -- ' + os.path.basename(fname))
  ax.set_xlabel('x')
  ax.set_ylabel('y')

  annot = ax.annotate('', xy=(0, 0), xytext=(15, 15), textcoords='offset points',
                       bbox=dict(boxstyle='round', fc='lightyellow'),
                       arrowprops=dict(arrowstyle='->'))
  annot.set_visible(False)
  marker, = ax.plot([], [], 'o', color='red', ms=markersize + 5, zorder=5)

  def on_pick(event):
    ind = event.ind[0]
    i1 = ind + 1
    label = 'i=' + str(i1) + '\n(' + '{:.6g}'.format(x[ind]) + ', ' + '{:.6g}'.format(y[ind]) + ')'
    segname = segment_names_for_index(segments, i1) if segments else None
    if segname:
      label += '\n' + segname
    annot.xy = (x[ind], y[ind])
    annot.set_text(label)
    annot.set_visible(True)
    marker.set_data([x[ind]], [y[ind]])
    fig.canvas.draw_idle()
    print(label.replace('\n', '   '))

  fig.canvas.mpl_connect('pick_event', on_pick)
  plt.show()


### End of functions ###
################################################################################

# Run main program
if __name__ == '__main__':
  main()
