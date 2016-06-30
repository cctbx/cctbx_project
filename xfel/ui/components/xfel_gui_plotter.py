from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/30/2016
Last Changed: 06/30/2016
Description : XFEL UI Plots and Charts
'''

import wx

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import xfel.ui.components.xfel_gui_controls as gctr


class SingleBarPlot(gctr.CtrlBase):
  def __init__(self, parent,
               label='',
               label_size=(80, -1),
               label_style='normal',
               content_style='normal',
               gauge_size=(250, 15),
               button=False,
               button_label='View Stats',
               button_size=wx.DefaultSize,
               choice_box=True,
               choice_label='',
               choice_label_size=(120, -1),
               choice_size=(100, -1),
               choice_style='normal',
               choices=[],
               gauge_max=100):
    gctr.CtrlBase.__init__(self, parent=parent, label_style=label_style,
                           content_style=content_style)

    self.sizer = wx.FlexGridSizer(1, 6, 0, 5)
    self.sizer.AddGrowableCol(3)

    dpi = wx.ScreenDC().GetPPI()[0]
    figsize = (gauge_size[0] / dpi, gauge_size[1] / dpi)
    self.status_figure = Figure(figsize=figsize)
    self.status_figure.patch.set_alpha(0)
    self.ax = self.status_figure.add_subplot(111)
    self.canvas = FigureCanvas(self, -1, self.status_figure)

    if choice_box:
      self.bins = gctr.ChoiceCtrl(self,
                                  label=choice_label,
                                  label_size=choice_label_size,
                                  label_style=choice_style,
                                  ctrl_size=choice_size,
                                  choices=choices)

    self.txt_iso = wx.StaticText(self, label=label, size=label_size)
    self.txt_max = wx.StaticText(self, label='{:.2f}'.format(gauge_max))
    self.txt_min = wx.StaticText(self, label='0')
    self.sizer.Add(self.txt_iso, flag=wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.txt_min, flag=wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.canvas, 1, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.txt_max, flag=wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.bins, flag=wx.ALIGN_CENTER_VERTICAL)

    if button:
      self.btn = wx.Button(self, label=button_label, size=button_size)
      self.sizer.Add(self.btn, 1, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)

    self.SetSizer(self.sizer)


class PopUpCharts(object):
  ''' Class to generate chargs and graphs that will appear in separate
  windows when user requests them, e.g. unit cell histogram chart '''

  def __init__(self):
    pass

  def plot_uc_histogram(self, isoform, info):

    fig = plt.figure(figsize=(12, 9))
    gsp = GridSpec(2, 3)
    sub_a = fig.add_subplot(gsp[0])
    sub_b = fig.add_subplot(gsp[1])
    sub_c = fig.add_subplot(gsp[2])
    sub_alpha = fig.add_subplot(gsp[3])
    sub_beta = fig.add_subplot(gsp[4])
    sub_gamma = fig.add_subplot(gsp[5])
