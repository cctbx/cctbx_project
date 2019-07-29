from __future__ import absolute_import, division, print_function

from mmtbx.validation.ramalyze import draw_ramachandran_plot

import math
from scipy import interpolate
import numpy as np

def calculate_indexes(x, y, xmin, x_step, ymin, y_step):
  i = math.floor((x-xmin) / x_step)
  j = math.floor((y-ymin) / y_step)
  return int(i), int(j)

class Grid2D(object):
  def __init__(self, data, xmin, xmax, ymin, ymax):
    # data - list of lists [[],[],...[]], data[x][y]
    self.g = data
    self.xmax = xmax
    self.xmin = xmin
    self.ymax = ymax
    self.ymin = ymin
    n_x = len(data)
    n_y = len(data[0])


    self.x_step = (xmax-xmin)/n_x
    self.y_step = (ymax-ymin)/n_y

    # print("x y steps:", self.x_step, self.y_step)

    # assert 360 % self.phi_step == 0
    # assert 360 // self.phi_step % 2 == 0
    # assert 360 % self.psi_step == 0
    # assert 360 // self.psi_step % 2 == 0
    # self.n_phi_half = 360 // self.phi_step // 2
    # self.n_psi_half = 360 // self.psi_step // 2

    # for i in range(360 // self.phi_step):
    #   self.g.append([0]*(360 // self.psi_step))

    # self.n_points = 0
    # self.mean = None
    # self.std = None
    self.interpolation_f = self.set_interpolation_f()

  @classmethod
  def make_Grid2d(cls, points, x_nbins, y_nbins,
      xmin=None, xmax=None, ymin=None, ymax=None):
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
      i,j = calculate_indexes(p[0], p[1], xmin, x_step, ymin, y_step)
      data[i][j] += 1
    return cls(data, xmin, xmax, ymin, ymax)

  def set_interpolation_f(self):
    x = []
    y = []
    for target_list, vmin, vmax, vstep in [(x, self.xmin, self.xmax, self.x_step),
        (y, self.ymin, self.ymax, self.y_step)]:
      current_tick = vmin - vstep/2
      while current_tick <= vmax + vstep/2 + 1e-6:
        target_list.append(current_tick)
        current_tick += vstep

    z = []
    # print "x,y", x, y
    for i in range(len(self.g[0])+2):
      z.append([0] * (len(self.g)+2) )
    for i in range(len(z)):
      for j in range(len(z)):
        # figure out where to get value
        ii = i-1
        jj = j-1
        if i == 0:
          ii = len(self.g)-1
        if i == len(z) - 1:
          ii = 0
        if j == 0:
          jj = len(self.g[0])-1
        if j == len(z) - 1:
          jj = 0
        z[i][j] = self.g[jj][ii]
    # print ("len(x)", len(x))
    # print ("len(y)", len(y))
    # print ("len(z)", len(z), len(z[0]))
    # print ("x", x)
    # print ("y", y)
    self.interpolation_f = interpolate.interp2d(x,y,z, kind='linear')

  def get_interpolated_score(self, x, y):
    # i, j = self._calc_ij(x,y)
    # return [self.g[i][j]]
    if self.interpolation_f is None:
      self.set_interpolation_f()
    # int_sc = self.interpolation_f([point[0]], [point[1]])
    # print "Get interp score:", point, int_sc, self.g[i][j]
    return self.interpolation_f([x],[y])[0]
    # return [self.g[i][j], self.interpolation_f([point[0]], [point[1]])[0]]

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
