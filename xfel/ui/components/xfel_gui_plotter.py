from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/30/2016
Last Changed: 06/30/2016
Description : XFEL UI Plots and Charts
'''

import wx
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import xfel.ui.components.xfel_gui_controls as gctr

class SingleBarPlot(gctr.CtrlBase):
  def __init__(self, parent,
               label='',
               label_size=wx.DefaultSize,
               label_style='normal',
               content_style='normal',
               gauge_size=(250, 15),
               button_label='View Stats',
               button_size=wx.DefaultSize,
               choice_box=True,
               choice_label='',
               choice_label_size=wx.DefaultSize,
               choice_size=(100, -1),
               choice_style='normal',
               choices=[]):
    gctr.CtrlBase.__init__(self, parent=parent, label_style=label_style,
                           content_style=content_style)

    self.sizer = wx.FlexGridSizer(1, 4, 0, 5)
    self.sizer.AddGrowableCol(1)
    self.dpi = wx.ScreenDC().GetPPI()[0]

    figsize = (gauge_size[0] / self.dpi, gauge_size[1] / self.dpi)
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
    self.sizer.Add(self.txt_iso, flag=wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.canvas, 1, flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
    self.sizer.Add(self.bins, flag=wx.ALIGN_CENTER_VERTICAL)
    self.btn = wx.Button(self, label=button_label, size=button_size)
    self.sizer.Add(self.btn, 1, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)

    self.SetSizer(self.sizer)

  def redraw_axes(self, value, goal, xmax):
    ''' Re-draw axes with latest values '''

    self.ax.clear()
    bar = self.ax.barh(0, value, height=1, align='center', color='#7570b3')

    xloc = xmax / 0.8
    label = '{:.2f}'.format(value)
    yloc = bar[0].get_y() + bar[0].get_height() / 2.0
    self.ax.text(xloc, yloc, label, horizontalalignment='right',
                 verticalalignment='center', weight='bold',
                 clip_on=False)

    self.ax.axvline(x=goal, lw=4, c='#d95f02')
    self.ax.set_xlim(xmax=xmax / 0.8)
    self.ax.axis('off')
    self.canvas.draw()
    self.Fit()


class NoBarPlot(gctr.CtrlBase):
  def __init__(self, parent,
               label='',
               label_size=wx.DefaultSize,
               label_style='bold'):
    gctr.CtrlBase.__init__(self, parent=parent, label_style=label_style,
                           content_style='bold')

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.info_sizer=wx.FlexGridSizer(1, 3, 0, 10)

    self.iso_txt = wx.StaticText(self, label=label, size=label_size)
    self.num_txt = wx.StaticText(self, label='0')
    self.end_txt = wx.StaticText(self, label='images integrated')
    self.iso_txt.SetFont(wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD))
    self.num_txt.SetFont(wx.Font(20, wx.DEFAULT, wx.NORMAL, wx.BOLD))
    self.end_txt.SetFont(wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD))

    self.info_sizer.Add(self.iso_txt, flag=wx.ALIGN_CENTER_VERTICAL)
    self.info_sizer.Add(self.num_txt, flag=wx.ALIGN_CENTER_VERTICAL)
    self.info_sizer.Add(self.end_txt, flag=wx.ALIGN_CENTER_VERTICAL)

    self.sizer.Add(self.info_sizer, flag=wx.ALIGN_CENTER)
    self.SetSizer(self.sizer)

  def update_number(self, number):
    self.num_txt.SetLabel(str(number))

class PopUpCharts(object):
  ''' Class to generate chargs and graphs that will appear in separate
  windows when user requests them, e.g. unit cell histogram chart '''

  def __init__(self):
    pass

  def plot_uc_histogram(self, info):

    # Initialize figure
    fig = plt.figure(figsize=(12, 10))
    gsp = GridSpec(2, 3)

    # Extract uc dimensions from info list
    a = [i[0] for i in info]
    b = [i[1] for i in info]
    c = [i[2] for i in info]
    alpha = [i[3] for i in info]
    beta = [i[4] for i in info]
    gamma = [i[5] for i in info]

    nbins = int(np.sqrt(len(info))) * 2

    fig.suptitle('Histogram of Unit Cell Dimensions ({} images)'
                 ''.format(len(info)), fontsize=18)

    sub_a = fig.add_subplot(gsp[0])
    sub_a.hist(a, nbins, normed=False, facecolor='#2c7fb8',
             alpha=0.75, histtype='stepfilled')
    sub_a.set_xlabel("a-edge ($\AA$)")
    sub_a.set_ylabel('Number of images')

    sub_b = fig.add_subplot(gsp[1], sharey=sub_a)
    sub_b.hist(b, nbins, normed=False, facecolor='#2c7fb8',
             alpha=0.75, histtype='stepfilled')
    sub_b.set_xlabel("b-edge ($\AA$)")
    plt.setp(sub_b.get_yticklabels(), visible=False)
    sub_b.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_b.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    sub_c = fig.add_subplot(gsp[2], sharey=sub_a)
    sub_c.hist(c, nbins, normed=False, facecolor='#2c7fb8',
             alpha=0.75, histtype='stepfilled')
    sub_c.set_xlabel("c-edge ($\AA$)")
    plt.setp(sub_c.get_yticklabels(), visible=False)
    sub_c.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_c.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    sub_alpha = fig.add_subplot(gsp[3])
    sub_alpha.hist(alpha, nbins, normed=False, facecolor='#7fcdbb',
                   alpha=0.75,  histtype='stepfilled')
    sub_alpha.set_xlabel(r'$\alpha (\circ)$')
    sub_alpha.set_ylabel('Number of images')

    sub_beta = fig.add_subplot(gsp[4], sharey=sub_alpha)
    sub_beta.hist(beta, nbins, normed=False, facecolor='#7fcdbb',
                  alpha=0.75, histtype='stepfilled')
    sub_beta.set_xlabel(r'$\beta (\circ)$')
    plt.setp(sub_beta.get_yticklabels(), visible=False)
    sub_beta.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_beta.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    sub_gamma = fig.add_subplot(gsp[5], sharey=sub_alpha)
    sub_gamma.hist(gamma, nbins, normed=False, facecolor='#7fcdbb',
                   alpha=0.75, histtype='stepfilled')
    sub_gamma.set_xlabel(r'$\gamma (\circ)$')
    plt.setp(sub_gamma.get_yticklabels(), visible=False)
    sub_gamma.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_gamma.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    gsp.update(wspace=0)
    plt.show()
