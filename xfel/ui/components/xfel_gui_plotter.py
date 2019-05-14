from __future__ import absolute_import, division, print_function
from six.moves import range

'''
Author      : Lyubimov, A.Y.
Created     : 06/30/2016
Last Changed: 06/30/2016
Description : XFEL UI Plots and Charts
'''

import wx
import numpy as np
from scitbx.array_family import flex

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import xfel.ui.components.xfel_gui_controls as gctr

class DoubleBarPlot(gctr.CtrlBase):
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

  def redraw_axes(self, valuea, valueb, goal, xmax, minimalist=False):
    ''' Re-draw axes with latest values '''

    self.ax.clear()
    bar = self.ax.barh(0, valueb, height=1, align='center', color='green')
    bar = self.ax.barh(0, valuea, height=1, align='center', color='#7570b3')

    xloc = xmax / 0.8
    label = '{:.1f}({:.1f})'.format(valueb, valuea)
    yloc = bar[0].get_y() + bar[0].get_height() / 2.0
    self.ax.text(xloc, yloc, label, horizontalalignment='right',
                 verticalalignment='center', weight='bold',
                 clip_on=False)

    if not minimalist:
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

  def __init__(self, interactive = True):
    import matplotlib.pyplot as plt
    self.plt=plt
    self.interactive = interactive

  def reject_outliers(self, data, iqr_ratio = 1.5):
    eps = 1e-6
    outliers = flex.bool(len(data), False)
    if iqr_ratio is None:
      return outliers
    from scitbx.math import five_number_summary
    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(data)
    #print "Five number summary: min %.1f, q1 %.1f, med %.1f, q3 %.1f, max %.1f"%(min_x, q1_x, med_x, q3_x, max_x)
    iqr_x = q3_x - q1_x
    cut_x = iqr_ratio * iqr_x
    outliers.set_selected(data > q3_x + cut_x + eps, True)
    outliers.set_selected(data < q1_x - cut_x - eps, True)
    #print "Rejecting", outliers.count(True), "out of", len(outliers)
    return outliers

  def plot_uc_histogram(self, info_list, legend_list, extra_title = None, xsize = 10, ysize = 10, high_vis = False, iqr_ratio = 1.5, ranges = None, title = None):
    """
    Plot a 3x3 grid of plots showing unit cell dimensions.
    @param info list of lists of dictionaries. The outer list groups seperate lists
    of cells to be plotted in the same graph, where each dictionary describes one cell.
    @param extra_title will be added to the title of the plot
    @param xsize if class initialized with not interacive, this is the x size of the
    plot to save in inches
    @param ysize as xsize
    @param iqr_ratio Inter-quartile range multiplier for rejecting outliers
    @param ranges Limits for the a, b and c axes. Tuple of 6 doubles, low then high in pairs for each.
    @return if not interactive, returns the path of the saved image, otherwise None
    """
    if ranges is not None:
      assert len(ranges) == 6
      alim = ranges[0:2]
      blim = ranges[2:4]
      clim = ranges[4:6]
    else:
      alim = blim = clim = None

    plot_ratio = max(min(xsize, ysize)/2.5, 3)
    if high_vis:
      text_ratio = plot_ratio*4
      separator = "\n"
    else:
      text_ratio = plot_ratio*3
      separator = "\n"

    # Initialize figure
    fig = plt.figure(figsize=(xsize, ysize))
    gsp = GridSpec(3, 4)
    legend_sub_a = fig.add_subplot(gsp[3])
    legend_sub_b = fig.add_subplot(gsp[7])
    legend_sub_c = fig.add_subplot(gsp[11])
    sub_ba = fig.add_subplot(gsp[0])
    sub_cb = fig.add_subplot(gsp[1])
    sub_ac = fig.add_subplot(gsp[2])
    sub_a = fig.add_subplot(gsp[4])
    sub_b = fig.add_subplot(gsp[5], sharey=sub_a)
    sub_c = fig.add_subplot(gsp[6], sharey=sub_a)
    sub_alpha = fig.add_subplot(gsp[8])
    sub_beta = fig.add_subplot(gsp[9], sharey=sub_alpha)
    sub_gamma = fig.add_subplot(gsp[10], sharey=sub_alpha)
    total = 0
    abc_hist_ylim = 0

    legend_sub_a.axis('off')
    legend_sub_b.axis('off')
    legend_sub_c.axis('off')

    for legend, info in zip(legend_list, info_list):
      if len(info) == 0:
        continue
      # Extract uc dimensions from info list
      a = flex.double([i['a'] for i in info])
      b = flex.double([i['b'] for i in info])
      c = flex.double([i['c'] for i in info])
      alpha = flex.double([i['alpha'] for i in info])
      beta = flex.double([i['beta'] for i in info])
      gamma = flex.double([i['gamma'] for i in info])
      if ranges is not None:
        sel = (a >= alim[0]) & (a <= alim[1]) & (b >= blim[0]) & (b <= blim[1]) & (c >= clim[0]) & (c <= clim[1])
        a = a.select(sel)
        b = b.select(sel)
        c = c.select(sel)
        alpha = alpha.select(sel)
        beta = beta.select(sel)
        gamma = gamma.select(sel)

      accepted = flex.bool(len(a), True)
      for d in [a, b, c, alpha, beta, gamma]:
        outliers = self.reject_outliers(d, iqr_ratio)
        accepted &= ~outliers

      a = a.select(accepted)
      b = b.select(accepted)
      c = c.select(accepted)
      alpha = alpha.select(accepted)
      beta = beta.select(accepted)
      gamma = gamma.select(accepted)

      total += len(a)
      nbins = int(np.sqrt(len(a))) * 2
      n_str = "N: %d "%len(a)

      from matplotlib.ticker import FormatStrFormatter
      hists = []
      for name, dimension, sub, lim in \
        [('a', a, sub_a, alim), ('b', b, sub_b, blim), ('c', c, sub_c, clim)]:
        stats = flex.mean_and_variance(dimension)
        mean = stats.mean()
        try:
          stddev = stats.unweighted_sample_standard_deviation()
        except RuntimeError:
          raise Exception("Not enough data to produce a histogram")
        varstr = "%.2f +/- %.2f"%(mean, stddev)
        if len(legend) > 0:
          dim_legend = legend + separator + varstr
        else:
          dim_legend = varstr
        if len(info_list) > 1 and name == "a":
          dim_legend = n_str + dim_legend
        hist = sub.hist(dimension, nbins, normed=False,
                 alpha=0.75, histtype='stepfilled', label = dim_legend, range = lim)
        sub.set_xlabel("%s-edge (%s $\AA$)"%(name, varstr)).set_fontsize(text_ratio)
        xloc = plt.MaxNLocator(5)
        if not high_vis:
          sub.xaxis.set_major_locator(xloc)
          sub.xaxis.set_major_formatter(FormatStrFormatter("%5.1f"))
        if name == 'a':
          sub.set_ylabel('Number of images').set_fontsize(text_ratio)
        else:
          self.plt.setp(sub.get_yticklabels(), visible=False)
        hists.append(hist)

      abc_hist_ylim = max(1.2*max([max(h[0]) for h in hists]), abc_hist_ylim)
      sub_a.set_ylim([0, abc_hist_ylim])

      for (n1, n2, d1, d2, lim1, lim2, sub) in \
        [('a', 'b', a, b, alim, blim, sub_ba),
         ('b', 'c', b, c, blim, clim, sub_cb),
         ('c', 'a', c, a, clim, alim, sub_ac)]:
        if len(info_list) == 1:
          sub.hist2d(d1, d2, bins=100, range=[lim1, lim2] if ranges is not None else None)
        else:
          sub.plot(d1.as_numpy_array(), d2.as_numpy_array(), '.')
          if ranges is not None:
            sub.set_xlim(lim1)
            sub.set_ylim(lim2)
        sub.set_xlabel("%s axis"%n1).set_fontsize(text_ratio)
        sub.set_ylabel("%s axis"%n2).set_fontsize(text_ratio)
        # plt.setp(sub.get_yticklabels(), visible=False)

      for ax in (sub_a, sub_b, sub_c, sub_alpha, sub_beta, sub_gamma):
        ax.tick_params(axis='both', which='both', left='off', right='off')
        ax.set_yticklabels([])
      for ax in (sub_ba, sub_cb, sub_ac):
        ax.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off')
        ax.set_xticklabels([])
        ax.set_yticklabels([])

      for (name, angle, sub) in \
        [(r'$\alpha$', alpha, sub_alpha),
         (r'$\beta$', beta, sub_beta),
         (r'$\gamma$', gamma, sub_gamma)]:
        sub.hist(angle, nbins, normed=False,
                 alpha=0.75, histtype='stepfilled')
        stats = flex.mean_and_variance(angle)
        mean = stats.mean()
        stddev = stats.unweighted_sample_standard_deviation()
        sub.set_xlabel(r'%s (%.2f +/- %.2f$^\circ$)' % (name, mean, stddev)).set_fontsize(text_ratio)
        xloc = plt.MaxNLocator(5)
        if not high_vis:
          sub.xaxis.set_major_locator(xloc)
          sub.xaxis.set_major_formatter(FormatStrFormatter("%5.1f"))
        if name == '\alpha':
          sub.set_ylabel('Number of images').set_fontsize(text_ratio)
        else:
          self.plt.setp(sub.get_yticklabels(), visible=False)

    sub_b.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_b.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_c.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_c.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_ba.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_ba.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_cb.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_cb.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_ac.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_ac.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_beta.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_beta.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    sub_gamma.xaxis.get_major_ticks()[0].label1.set_visible(False)
    sub_gamma.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    h, l = sub_a.get_legend_handles_labels()
    legend_sub_a.legend(h, l, fontsize=text_ratio)
    h, l = sub_b.get_legend_handles_labels()
    legend_sub_b.legend(h, l, fontsize=text_ratio)
    h, l = sub_c.get_legend_handles_labels()
    legend_sub_c.legend(h, l, fontsize=text_ratio)

    # if len(legend_list) > 0:
    #   import matplotlib.patches as mpatches
    #   rgb_alphas = [a_hist[2][i] for i in range(len(a_hist[2]))]
    #   assert len(legend_list) == len(rgb_alphas)
    #   patches = [mpatches.Patch(
    #     color=rgb_alphas[i].get_facecolor(),
    #     label=legend_list[i])
    #   for i in range(len(rgb_alphas))]
    #   fig.legend(patches, 'upper right')

    gsp.update(wspace=0)
    if title is None: title = "Unit cell distribution"
    fig.suptitle(title + " (%d xtals)" % total)

    if not self.interactive:
      fig.set_size_inches(xsize*1.05+.5, ysize*.95)
      fig.savefig("ucell_tmp.png", bbox_inches='tight', dpi=100)
      plt.close(fig)
      return "ucell_tmp.png"

  def plot_uc_3Dplot(self, info, iqr_ratio = 1.5):
    assert self.interactive

    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D # import dependency

    fig = self.plt.figure(figsize=(12, 10))
    # Extract uc dimensions from info list
    a = flex.double([i['a'] for i in info])
    b = flex.double([i['b'] for i in info])
    c = flex.double([i['c'] for i in info])
    alpha = flex.double([i['alpha'] for i in info])
    beta = flex.double([i['beta'] for i in info])
    gamma = flex.double([i['gamma'] for i in info])
    n_total = len(a)

    accepted = flex.bool(n_total, True)
    for d in [a, b, c, alpha, beta, gamma]:
      outliers = self.reject_outliers(d, iqr_ratio)
      accepted &= ~outliers

    a = a.select(accepted)
    b = b.select(accepted)
    c = c.select(accepted)

    AA = "a-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(a),
                                        flex.mean_and_variance(a).unweighted_sample_standard_deviation())
    BB = "b-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(b),
                                        flex.mean_and_variance(b).unweighted_sample_standard_deviation())
    CC = "c-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(c),
                                        flex.mean_and_variance(c).unweighted_sample_standard_deviation())


    subset = min(len(a),1000)

    flex.set_random_seed(123)
    rnd_sel = flex.random_double(len(a))<(subset/n_total)

    a = a.select(rnd_sel)
    b = b.select(rnd_sel)
    c = c.select(rnd_sel)

    fig.suptitle('{} randomly selected cells out of total {} images'
                 ''.format(len(a),n_total), fontsize=18)

    ax = fig.add_subplot(111, projection='3d')

    for ia in range(len(a)):
      ax.scatter(a[ia],b[ia],c[ia],c='r',marker='+')

    ax.set_xlabel(AA)
    ax.set_ylabel(BB)
    ax.set_zlabel(CC)
