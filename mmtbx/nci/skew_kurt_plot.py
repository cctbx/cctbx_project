from __future__ import absolute_import, division, print_function

import libtbx.load_env
from libtbx import easy_pickle
from mmtbx.utils.grid2d import Grid2D
import numpy as np
import matplotlib.pyplot as plt
from libtbx.utils import Sorry
import os

from matplotlib.ticker import FormatStrFormatter

import math

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def rotate_translate(origin, point, angle, translate):
  x,y = rotate(origin, point, angle)
  return x+translate[0], y+translate[1]

def get_filling_data(data, x_nbins, y_nbins, xmin, xmax, ymin, ymax,
    threshold=2, interpolation='linear'):
  grid = Grid2D.make_Grid2d(data, x_nbins=x_nbins, y_nbins=y_nbins,
      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
  grid.set_interpolation_f(interpolation)
  filtered_out = grid.filter_out_less_populated_points(threshold=threshold)
  npz = np.array(grid.g)
  npz = np.swapaxes(npz, 0, 1)
  npz = npz ** 0.5
  for x in range(len(npz)):
    for y in range( int(len(npz[0]) / 2), len(npz[0])):
      npz[x][y] = - npz[x][y]
  npz = npz.astype(float)
  nmin = np.amin(npz)
  nmax = np.amax(npz)
  for x in range(len(npz)):
    for y in range(len(npz[0])):
      if npz[x][y] > 0:
        npz[x][y] = npz[x][y] / float(nmax)
      if npz[x][y] < 0:
        npz[x][y] = npz[x][y] / abs(float(nmin))
  grid = Grid2D(npz, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
  grid.set_interpolation_f(interpolation)
  n_steps = 1000
  smooth_npz = []
  for i in range(n_steps):
    smooth_npz.append([0]*n_steps)
  for i, x in enumerate(np.arange(xmin, xmax, (xmax-xmin) / n_steps)):
    for j, y in enumerate(np.arange(ymin, ymax, (ymax-ymin) / n_steps)):
      smooth_npz[i][j] = grid.get_interpolated_score(x,y)
  smooth_npz = np.array(smooth_npz)
  smooth_npz.astype(float)
  # easy_pickle.dump(obj=smooth_npz, file_name="smooth_npz.pkl")
  return smooth_npz

def make_figure(
    file_name,
    theta1_coords,
    Rha_coords,
    dot_size,
    type='all',
    colorblind_friendly=True,
    override_palette = {},
    resolution=300,
    do_rotate=False,
    do_y_log=False,
    ):
  """ Plots skew-kurtosis plot and saves it to .png file

  Args:
      file_name (str): file name for the produced figure without extension
      theta1_coords (iterabe): list of skew-kurtosis coordinates of theta1 [(),(),..]
      Rha_coords (iterable): list of skew-kurtosis coordinates of Rha [(),(),..]
      type (str, optional): Type of contours to plot. Defaults to 'all'. Allowed 'all', 'alpha', 'beta'.
      colorblind_friendly (bool, optional): Use color-blind friendly palette. Defaults to True.
      override_palette (dict, optional): keys and values to overide color_palette values.
      resolution (int, optional): Resolution of outputted figure. Defaults to 300.
  """
  assert type in ['all', 'alpha', 'beta']
  pkl_fn = libtbx.env.find_in_repositories(
      relative_path="mmtbx/nci")+'/filtered_raw_data.pkl'
  assert os.path.isfile(pkl_fn)
  db = easy_pickle.load(pkl_fn)
  # settings
  color_palette = {'colormap': 'PiYG',
      'contour_left':'green',
      'contour_right':'mediumvioletred',
      'contour_thick_left':1.5,
      'contour_thick_right':1.5,
      'theta_color':'green',
      'theta_contour':'black',
      'Rha_color': 'deeppink',
      'Rha_contour': 'black',
      }
  if colorblind_friendly:
    color_palette = {'colormap': 'cividis',
        'contour_left':'yellow',
        'contour_right':'darkblue',
        'contour_thick_left':1.5,
        'contour_thick_right':1.5,
        'theta_color':'yellow',
        'theta_contour':'black',
        'Rha_color': 'darkblue',
        'Rha_contour': 'white',
        }
  for k,v in override_palette.items():
    if k in color_palette:
      color_palette[k] = v
    else:
      raise Sorry("%s parameter is not in the used palette options." % k)

  contours = {'all': [-0.14, 0.14],
      'alpha': [-0.14, 0.14],
      'beta': [-0.14, 0.14]}

  aspect = 0.5

  xmin = -2
  xmax = 2
  ymin = 0
  ymax = 7

  ta = 0.3
  tx = 1.
  ty = 0.5
  if do_rotate:
    xmin = -1.5
    xmax =  1.5
  if do_y_log:
    ymin = 1
    aspect = 2.5

  data = []
  for x, y in zip(db[type]['t1o'][0], db['all']['t1o'][1]):
    if do_rotate: x,y = rotate_translate((0,0), (x,y), -ta, (-tx,-ty))
    data.append([x,y])
  for x, y in zip(db[type]['dhao'][0], db['all']['dhao'][1]):
    if do_rotate: x,y = rotate_translate((0,0), (x,y), ta, (tx,-ty))
    data.append([x ,y])
  blue_filling = get_filling_data(data, x_nbins=50, y_nbins=50,
      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
      threshold=1, interpolation='cubic')
  fig = plt.figure(figsize=(10,10),)
  ax = plt.subplot(111)

  if do_y_log:
    ax.set_xscale('linear')
    ax.set_yscale('log')
    yax = ax.axes.get_yaxis()
    yax.set_major_formatter(FormatStrFormatter('%d'))
    yax.set_minor_formatter(FormatStrFormatter('%d'))

  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  im = ax.imshow(
      blue_filling,
      origin="lower",
      cmap=plt.get_cmap(color_palette['colormap']),
      extent=[xmin, xmax, ymin, ymax],
      aspect=aspect,
      # interpolation='bicubic',
      # norm=norm,
      alpha=0.7
      )
  plt.contour( blue_filling, contours[type],
      origin="lower",
      linestyles=['solid', 'solid'],
      colors=[color_palette['contour_right'],color_palette['contour_left']],
      linewidths=[color_palette['contour_thick_right'], color_palette['contour_thick_left']],
      extent=[xmin, xmax, ymin, ymax],
      )

  for theta1_c in theta1_coords:
    x,y = theta1_c[0], theta1_c[1]
    if do_rotate: x,y = rotate_translate((0,0), (x,y), -ta, (-tx,-ty))
    ax.scatter([x], [y], s=dot_size, c=color_palette['theta_color'], edgecolor=color_palette['theta_contour'])
  for Rha_c in Rha_coords:
    x,y = Rha_c[0], Rha_c[1]
    if do_rotate: x,y = rotate_translate((0,0), (x,y), ta, (tx,-ty))
    ax.scatter([x], [y], s=dot_size, c=color_palette['Rha_color'], edgecolor=color_palette['Rha_contour'])
  if do_rotate:
    xax = ax.axes.get_xaxis()
    xax = xax.set_visible(False)
    yax = ax.axes.get_yaxis()
    yax = yax.set_visible(False)
  else:
    ax.set_xlabel('Skew')
    ax.set_ylabel('Kurtosis')
  fig.savefig("%s.png" % file_name, dpi=resolution)

if __name__ == '__main__':
  make_figure()
