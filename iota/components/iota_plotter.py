# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from past.builtins import range
from six.moves import range
from six.moves import zip

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 01/30/2019
Description : IOTA Plotter module. Exists to provide custom MatPlotLib plots
              that are dynamic, fast-updating, interactive, and keep up with
              local wxPython upgrades.
'''

import os
import wx
from collections import Counter
import numpy as np
import math

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# workaround to avoid unused import warning
assert Axes3D
assert cm
assert colors

from iota.components.iota_ui_base import IOTABaseFrame, IOTABasePanel

class PlotWindow(IOTABaseFrame):
  def __init__(self, parent, id, title, plot_panel=None):
    IOTABaseFrame.__init__(self, parent, id, title, size=(500, 500))

    self.initialize_toolbar()
    self.tb_btn_quit = self.add_tool(label='Quit',
                                     bitmap=('actions', 'exit'),
                                     shortHelp='Quit')
    self.tb_btn_save = self.add_tool(label='Save',
                                     bitmap=('actions', 'save_all'),
                                     shortHelp='Save image in various formats')
    self.realize_toolbar()

    self.Bind(wx.EVT_TOOL, self.onSave, self.tb_btn_save)
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)

    self.plot_panel = plot_panel

  def plot(self):
    if self.plot_panel:
      self.main_sizer.Add(self.plot_panel, 1, flag=wx.EXPAND)
      self.SetSize(self.plot_panel.canvas.GetSize())
      self.plot_panel.canvas.draw()
      self.Layout()

  def onSave(self, e):
    save_dlg = wx.FileDialog(self,
                             message="Save Image",
                             defaultDir=os.curdir,
                             defaultFile="*",
                             wildcard="*",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      script_filepath = save_dlg.GetPath()
      self.plot_panel.figure.savefig(script_filepath, format='pdf',
                                     bbox_inches=0)

  def onQuit(self, e):
    self.Close()

class Plotter(IOTABasePanel):
  ''' Generic Plotter (will plot anything given specific data) '''

  def __init__(self, parent,
               params=None,
               info=None,
               *args, **kwargs):
    IOTABasePanel.__init__(self, parent=parent, *args, **kwargs)
    self.info = info
    self.params = params
    # if viz_dir:

    self.font = {'fontfamily':'sans-serif', 'fontsize':12}

  def initialize_figure(self, figsize=(9, 9)):
    self.figure = Figure(figsize=figsize)
    self.canvas = FigureCanvas(self, -1, self.figure)
    self.main_sizer.Add(self.canvas, 1, flag=wx.EXPAND)

  def plot_res_histogram(self):

    # Get resolution values
    hres = zip(*self.info.stats['res']['lst'])[2]
    lres = zip(*self.info.stats['lres']['lst'])[2]

    # Plot figure
    gsp = gridspec.GridSpec(2, 1)
    hr = self.figure.add_subplot(gsp[0, :])
    hr_n, hr_bins, hr_patches = hr.hist(hres, 20, facecolor='b', alpha=0.75,
                                         histtype='stepfilled')
    hr_height = (np.max(hr_n) + 9) // 10 * 10
    hr.axis([np.min(hres), np.max(hres), 0, hr_height])
    reslim = 'High Resolution Limit ({})'.format(r'$\AA$')
    hr.set_xlabel(reslim, fontsize=15)
    hr.set_ylabel('No. of frames', fontsize=15)

    lr = self.figure.add_subplot(gsp[1, :])
    lr_n, lr_bins, lr_patches = lr.hist(lres, 20, facecolor='b', alpha=0.75,
                                         histtype='stepfilled')
    lr_height = (np.max(lr_n) + 9) // 10 * 10
    lr.axis([np.min(lres), np.max(lres), 0, lr_height])
    reslim = 'Low Resolution Limit ({})'.format(r'$\AA$')
    lr.set_xlabel(reslim, fontsize=15)
    lr.set_ylabel('No. of frames', fontsize=15)

    self.figure.set_tight_layout(tight=True)


  def plot_spotfinding_heatmap(self):

    hlist = self.info.stats['h']['lst']
    alist = self.info.stats['a']['lst']

    ch = max(hlist) - min(hlist) + 1
    ca = max(alist) - min(alist) + 1
    ints = zip(hlist, alist)
    ic = Counter(ints)

    hm_data = np.zeros((ch, ca))
    for i in ic.items():
      hm_data[i[0][0]-min(hlist), i[0][1]-min(alist)] = i[1]

    rows = range(min(hlist), max(hlist) + 1)
    cols = range(min(alist), max(alist) + 1)
    row_labels = [str(i) for i in rows]
    col_labels = [str(j) for j in cols]

    ax = self.figure.add_subplot(111)
    ax.pcolor(hm_data, cmap='Reds')

    ax.set_yticks(np.arange(len(rows))+.5, minor=False)
    ax.set_xticks(np.arange(len(cols))+.5, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    ax.set_xticklabels(col_labels, minor=False)
    ax.set_xlabel('Spot area')
    ax.set_ylabel('Spot height')

    ax.set_xlim(0, len(cols))
    ax.set_ylim(0, len(rows))

    # Annotate
    for y in range(hm_data.shape[0]):
        for x in range(hm_data.shape[1]):
            ax.text(x + 0.5, y + 0.5, '%3d' % hm_data[y, x],
                     horizontalalignment='center',
                     verticalalignment='center',
                     )

    self.figure.set_tight_layout(True)

  def plot_beam_xy(self, write_files=False, return_values=False, threeD=False):
    """ Plot beam center coordinates and a histogram of distances from the median
        of beam center coordinates to each set of coordinates. Superpose a
        predicted mis-indexing shift by L +/- 1 (calculated for each axis).
    """

    # Calculate beam center coordinates and distances
    beamX = zip(*self.info.stats['beamX']['lst'])[2]
    beamY = zip(*self.info.stats['beamY']['lst'])[2]
    beamZ = zip(*self.info.stats['distance']['lst'])[2]
    beamXY = zip(beamX, beamY)

    beam_dist = [math.hypot(i[0] - np.median(beamX), i[1] - np.median(beamY))
                 for i in beamXY]
    beam_dist_std = np.std(beam_dist)
    beamXYdist = zip(beamX, beamY, beam_dist)

    # Separate out outliers
    outliers = [i for i in beamXYdist if i[2] > 2 * beam_dist_std]
    clean = [i for i in beamXYdist if i[2] <= 2 * beam_dist_std]
    cbeamX = [i[0] for i in clean]
    cbeamY = [j[1] for j in clean]
    obeamX = [i[0] for i in outliers]
    obeamY = [j[1] for j in outliers]

    wavelength = self.info.stats['wavelength']['median']
    det_distance = self.info.stats['distance']['median']
    a = np.median([i[0] for i in self.info.cluster_iterable])
    b = np.median([i[1] for i in self.info.cluster_iterable])
    c = np.median([i[2] for i in self.info.cluster_iterable])

    # Calculate predicted L +/- 1 misindexing distance for each cell edge
    aD = det_distance * math.tan(2 * math.asin(wavelength / (2 * a)))
    bD = det_distance * math.tan(2 * math.asin(wavelength / (2 * b)))
    cD = det_distance * math.tan(2 * math.asin(wavelength / (2 * c)))


    # Plot figure
    if threeD:
      self.figure.set_size_inches(w=8, h=8)
      ax1 = self.figure.add_subplot(111, projection='3d')
      Axes3D.mouse_init(ax1)
    else:
      self.figure.set_size_inches(w=9, h=13)
      gsp = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
      ax1 = self.figure.add_subplot(gsp[0, :], aspect='equal')

    # Calculate axis limits of beam center scatter plot
    ax1_delta = np.ceil(np.max(beam_dist))
    xmax = round(np.median(beamX) + ax1_delta)
    xmin = round(np.median(beamX) - ax1_delta)
    ymax = round(np.median(beamY) + ax1_delta)
    ymin = round(np.median(beamY) - ax1_delta)
    zmax = round(np.ceil(self.info.stats['distance']['max']))
    zmin = round(np.floor(self.info.stats['distance']['min']))

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    if threeD:
      ax1.set_zlim(zmin, zmax)

    # Plot beam center scatter plot
    if threeD:
      ax1.scatter(beamX, beamY, beamZ, alpha=1, s=20, c='grey', lw=1)
      ax1.plot([self.info.stats['beamX']['median']],
               [self.info.stats['beamY']['median']],
               [self.info.stats['distance']['median']],
               markersize=8, marker='o', c='yellow', lw=2)
    else:
      ax1.scatter(cbeamX, cbeamY, alpha=1, s=20, c='grey', lw=1)
      ax1.scatter(obeamX, obeamY, alpha=1, s=20, c='red', lw=1)
      ax1.plot(np.median(beamX), np.median(beamY), markersize=8, marker='o',
               c='yellow', lw=2)

      # Plot projected mis-indexing limits for all three axes
      from matplotlib.patches import Circle
      circle_a = Circle((np.median(beamX), np.median(beamY)), radius=aD,
                            color='r', fill=False, clip_on=True)
      circle_b = Circle((np.median(beamX), np.median(beamY)), radius=bD,
                            color='g', fill=False, clip_on=True)
      circle_c = Circle((np.median(beamX), np.median(beamY)), radius=cD,
                            color='b', fill=False, clip_on=True)
      ax1.add_patch(circle_a)
      ax1.add_patch(circle_b)
      ax1.add_patch(circle_c)

    # Set labels
    ax1.set_xlabel('BeamX (mm)', fontsize=15)
    ax1.set_ylabel('BeamY (mm)', fontsize=15)
    if threeD:
      ax1.set_zlabel('Distance (mm)', fontsize=15)
      ax1.set_title('Beam XYZ Coordinates')
    else:
      ax1.set_title('Beam XY Coordinates')

    if not threeD:
      # Plot histogram of distances to each beam center from median
      ax2 = self.figure.add_subplot(gsp[1, :])
      ax2_n, ax2_bins, ax2_patches = ax2.hist(beam_dist, 20, facecolor='b',
                                              alpha=0.75, histtype='stepfilled')
      ax2_height = (np.max(ax2_n) + 9) // 10 * 10
      ax2.axis([0, np.max(beam_dist), 0, ax2_height])
      ax2.set_xlabel('Distance from median (mm)', fontsize=15)
      ax2.set_ylabel('No. of images', fontsize=15)

    self.figure.set_tight_layout(True)
