from __future__ import absolute_import, division, print_function

from mmtbx.validation.ramalyze import draw_ramachandran_plot

import math
from scipy import interpolate
import numpy as np
from copy import deepcopy

def calculate_indexes(x, y, xmin, x_step, ymin, y_step):
  i = math.floor((x-xmin) / x_step)
  j = math.floor((y-ymin) / y_step)
  return int(i), int(j)

class Grid2D(object):
  def __init__(self, data, xmin, xmax, ymin, ymax):
    # data - list of lists [[],[],...[]], data[x][y]
    self.g = deepcopy(data)
    self.xmax = xmax
    self.xmin = xmin
    self.ymax = ymax
    self.ymin = ymin
    n_x = len(data)
    n_y = len(data[0])
    self.x_step = (xmax-xmin)/n_x
    self.y_step = (ymax-ymin)/n_y
    self.interpolation_f = None

  def filter_out_less_populated_points(self, threshold=1):
    n = 0
    for x in range(len(self.g)):
      for y in range(len(self.g[0])):
        if self.g[x][y] <= threshold:
          if self.g[x][y] != 0:
            n += 1
          self.g[x][y] = 0
    return n

  @classmethod
  def make_Grid2d(cls, points, x_nbins, y_nbins,
      xmin=None, xmax=None, ymin=None, ymax=None, allow_outside=True):
    # points - list of tuples [(x,y), (x,y)...]
    if xmin is None:
      print ("xmin = ", min([x[0] for x in points]) )
      xmin =  min([x[0] for x in points])
    if xmax is None:
      print ("xmax = ", max([x[0] for x in points]) )
      xmax = max([x[0] for x in points])
    if ymin is None:
      print ("ymin = ", min([x[1] for x in points]) )
      ymin = min([x[1] for x in points])
    if ymax is None:
      print ("ymax = ", max([x[1] for x in points]) )
      ymax = max([x[1] for x in points])

    x_step = (xmax-xmin)/x_nbins
    y_step = (ymax-ymin)/y_nbins
    data = []
    for i in range(x_nbins):
      data.append([0]*y_nbins)
    for p in points:
      if (p[0]<xmin or p[0]>xmax) and allow_outside:
        continue
      if (p[1]<ymin or p[1]>ymax) and allow_outside:
        continue
      i,j = calculate_indexes(p[0], p[1], xmin, x_step, ymin, y_step)
      # print (p)
      data[i][j] += 1
    return cls(data, xmin, xmax, ymin, ymax)

  def set_interpolation_f(self, interpolation_type='linear'):
    x = []
    y = []
    for target_list, vmin, vmax, vstep in [(x, self.xmin, self.xmax, self.x_step),
        (y, self.ymin, self.ymax, self.y_step)]:
      current_tick = vmin - vstep/2
      while current_tick <= vmax + vstep/2 + 1e-6:
        target_list.append(current_tick)
        current_tick += vstep

    z = np.array(self.g)
    z = np.swapaxes(z, 0, 1)
    z = np.pad(z,  pad_width=1, mode='wrap')

    self.interpolation_f = interpolate.interp2d(x,y,z, kind=interpolation_type)

  def get_interpolated_score(self, x, y):
    if self.interpolation_f is None:
      raise RuntimeError('function is not set yet')
    return self.interpolation_f([x],[y])[0]

  def plot_distribution(self, fname, title="Default title"):
    npz = np.array(self.g)
    npz = np.swapaxes(npz, 0, 1)
    # if you want ramachandran-like scaling
    npz = npz ** 0.25
    npz.astype(float)
    p = draw_ramachandran_plot(
        points=[],
        rotarama_data=npz,
        position_type=0, # RAMA_GENERAL
        title=title,
        show_labels=False,
        markerfacecolor="white",
        markeredgecolor="black",
        show_filling=True,
        show_contours=False,
        markersize=10,
        point_style='bo')
    # fname = "4sc_%s_%s.png" % (k, r_type)
    print ("saving", fname)
    p.save_image(fname, dpi=300)
