from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.small_cell.cake_plot
from matplotlib import pyplot as plt
import sys, math
from matplotlib.ticker import FuncFormatter
import numpy as np

help_str = """
Make a cake plot from DIALS spotfinder spots

A cake plot is the azimuthal angle of a spot on an image vs. its resolution.
Powder rings will appear as vertical stripes, with defects in geometry
causing them to appear wavy. A cake plot is also insensitive to badly masked
regions of the detector compared to a 1d radial average as the aziumuthal
angle of a spot isn't averaged into the 1d trace.

This script plots the results from cctbx.xfel.small_cell.cake_plot_prep.
Run that script first to generate cake.npy, which this script uses.

Usage:
cctbx.xfel.small_cell.cake_plot_prep
"""

def run(args):
  if "-h" in args or "--help" in args:
    print(help_str)
    return

  def resolution(x, pos):
    if x <= 0:
      return '-'
    return "%.3f"%(1/math.sqrt(x))
  formatter = FuncFormatter(resolution)
  plt.gca().xaxis.set_major_formatter(formatter)

  print('loading')
  with open('cake.npy', 'rb') as f:
    datasets = 0
    while True:
      try:
        x = np.load(f)
        y = np.load(f)
        d = np.load(f)
        azi = np.load(f)
        datasets += 1
        print("Dataset %d loaded"%datasets)
      except ValueError:
        break

      plt.figure(1)
      plt.xlabel("Resolution (A)"); plt.ylabel("Azimuthal angle relative to (0,1,0) (deg)"); plt.title("Azimuthal angle vs. resolution")
      plt.scatter((1/d**2), azi, s=.1, c='blue')

      plt.figure(2)
      plt.xlabel("mm"); plt.ylabel("mm"); plt.title("Spot position in lab space")
      plt.scatter(x,y,s=.1,c='blue')

  plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
