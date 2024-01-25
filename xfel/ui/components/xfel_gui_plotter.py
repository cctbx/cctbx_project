from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip

'''
Author      : Lyubimov, A.Y.
Created     : 06/30/2016
Last Changed: 06/30/2016
Description : XFEL UI Plots and Charts
'''

import wx
import numpy as np
from scitbx.array_family import flex

import matplotlib as mpl
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


class CommonUnitCellKey(object):
  """Handle unit cell parameters when setting a common legend for histograms"""
  sep = '\n'
  symbols = ['a', 'b', 'c', r'$\alpha$', r'$\beta$', r'$\gamma$']
  units = 3 * [r'$\AA$'] + 3 * [r'$^\circ$']

  def __init__(self, name='', crystals=0):
    self.name = name
    self.crystals = crystals
    self.means = []
    self.stds = []

  @property
  def prefix(self):
    return self.name + (' ' if self.name else '') + "(N: %d)" % self.crystals

  @property
  def line_list(self):
    return ['%s: %.2f +/- %.2f%s' % (d, m, s, u) for d, m, s, u
            in zip(self.symbols, self.means, self.stds, self.units)]

  def lines(self, mask=6*(True,)):
    return self.sep.join(line for line, m in zip(self.line_list, mask) if m)

  @classmethod
  def common_lines_of(cls, uc_keys):
    line_lists = zip(*[uc_key.line_list for uc_key in uc_keys])
    return [ll.count(ll[0]) == len(ll) for ll in line_lists] # true if all same


class PopUpCharts(object):
  ''' Class to generate chargs and graphs that will appear in separate
  windows when user requests them, e.g. unit cell histogram chart '''

  def __init__(self, interactive = True, figure = None):
    import matplotlib.pyplot as plt
    self.plt=plt
    self.interactive = interactive
    self.figure = figure

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

  def plot_uc_histogram(self, info_list, legend_list, extra_title = None, xsize = 10, ysize = 10, high_vis = False, iqr_ratio = 1.5, ranges = None, angle_ranges = None, title = None, image_fname=None, hist_scale=None):
    """
    Plot a 3x3 grid of plots showing unit cell dimensions.
    @param info list of lists of dictionaries. The outer list groups seperate lists
    of cells to be plotted in the same graph, where each dictionary describes one cell.
    @param extra_title if given will be appended to the title of the plot
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

    if angle_ranges is not None:
      assert len(angle_ranges) == 6
      allim = angle_ranges[0:2]
      belim = angle_ranges[2:4]
      galim = angle_ranges[4:6]
    else:
      allim = belim = galim = None

    plot_ratio = max(min(xsize, ysize)/2.5, 3)
    text_ratio = plot_ratio * (4 if high_vis else 3)

    # Initialize figure
    if self.figure:
      fig = self.figure
    else:
      fig = plt.figure(figsize=(xsize, ysize))
    gsp = GridSpec(3, 4)
    sub_ba = fig.add_subplot(gsp[0, 0])
    sub_cb = fig.add_subplot(gsp[0, 1])
    sub_ac = fig.add_subplot(gsp[0, 2])
    sub_a = fig.add_subplot(gsp[1, 0], sharex=sub_ba)
    sub_b = fig.add_subplot(gsp[1, 1], sharex=sub_cb, sharey=sub_a)
    sub_c = fig.add_subplot(gsp[1, 2], sharex=sub_ac, sharey=sub_a)
    sub_alpha = fig.add_subplot(gsp[2, 0])
    sub_beta = fig.add_subplot(gsp[2, 1], sharey=sub_alpha)
    sub_gamma = fig.add_subplot(gsp[2, 2], sharey=sub_alpha)
    sub_key = fig.add_subplot(gsp[:, 3])
    total = 0
    abc_hist_ylim = 0
    legend_keys = []

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
        axis_sel = (a >= alim[0]) & (a <= alim[1]) & (b >= blim[0]) & (b <= blim[1]) & (c >= clim[0]) & (c <= clim[1])
      else:
        axis_sel = flex.bool(len(a), True)
      if angle_ranges is not None:
        angle_sel = (alpha >= allim[0]) & (alpha <= allim[1]) & (beta >= belim[0]) & (beta <= belim[1]) & (gamma >= galim[0]) & (gamma <= galim[1])
      else:
        angle_sel = flex.bool(len(a), True)
      sel = axis_sel & angle_sel
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
      hists = []
      legend_key = CommonUnitCellKey(name=legend, crystals=len(a))

      for (d, sub, lim) in [(a, sub_a, alim), (b, sub_b, blim), (c, sub_c, clim)]:
        stats = flex.mean_and_variance(d)
        mean = stats.mean()
        try:
          stddev = stats.unweighted_sample_standard_deviation()
        except RuntimeError:
          raise Exception("Not enough data to produce a histogram")
        legend_key.means.append(mean)
        legend_key.stds.append(stddev)
        hist = sub.hist(d, nbins, alpha=0.75, histtype='stepfilled',
                        label='placeholder-label', range=lim)
        hists.append(hist)
        if len(info_list) == 1:
          sub.set_xlabel(legend_key.line_list[-1]).set_fontsize(text_ratio)

      abc_hist_ylim = max(1.2*max([max(h[0]) for h in hists]), abc_hist_ylim)
      sub_a.set_ylim([0, abc_hist_ylim])

      for (n1, n2, d1, d2, lim1, lim2, sub) in \
        [('a', 'b', a, b, alim, blim, sub_ba),
         ('b', 'c', b, c, blim, clim, sub_cb),
         ('c', 'a', c, a, clim, alim, sub_ac)]:
        if len(info_list) == 1:
          hist_kwargs = {
            'norm': mpl.colors.LogNorm() if hist_scale == "log" else None,
            'range': [lim1, lim2] if ranges is not None else None}
          sub.hist2d(d1, d2, bins=100, **hist_kwargs)
        else:
          sub.plot(d1.as_numpy_array(), d2.as_numpy_array(), '.', alpha=0.1,
                   markeredgewidth=0, markersize=2)
          if ranges is not None:
            sub.set_xlim(lim1)
            sub.set_ylim(lim2)
        sub.set_xlabel("%s axis" % n1).set_fontsize(text_ratio)
        sub.set_ylabel("%s axis" % n2).set_fontsize(text_ratio)

      for (angle, sub, lim) in [(alpha, sub_alpha, allim), (beta, sub_beta, belim), (gamma, sub_gamma, galim)]:
        sub.hist(angle, nbins, alpha=0.75, histtype='stepfilled', range=lim)
        stats = flex.mean_and_variance(angle)
        mean = stats.mean()
        stddev = stats.unweighted_sample_standard_deviation()
        legend_key.means.append(mean)
        legend_key.stds.append(stddev)
        if len(info_list) == 1:
          sub.set_xlabel(legend_key.line_list[-1]).set_fontsize(text_ratio)
      legend_keys.append(legend_key)

    # Set up general subplot and legend information
    sub_a.set_ylabel('Number of images').set_fontsize(text_ratio)
    self.plt.setp(sub_b.get_yticklabels(), visible=False)
    self.plt.setp(sub_c.get_yticklabels(), visible=False)
    for sub in (sub_a, sub_b, sub_c):
      if not high_vis:
        sub.xaxis.set_major_locator(plt.MaxNLocator(4))

    sub_alpha.set_ylabel('Number of images').set_fontsize(text_ratio)
    self.plt.setp(sub_beta.get_yticklabels(), visible=False)
    self.plt.setp(sub_gamma.get_yticklabels(), visible=False)
    for sub in [sub_alpha, sub_beta, sub_gamma]:
      if not high_vis:
        sub.xaxis.set_major_locator(plt.MaxNLocator(4))

    for ax in (sub_a, sub_b, sub_c, sub_alpha, sub_beta, sub_gamma):
      ax.tick_params(axis='both', which='both', left='off', right='off')
      ax.set_yticklabels([])
    for ax in (sub_ba, sub_cb, sub_ac):
      ax.tick_params(axis='both', which='both', bottom='off', top='off',
                     left='off', right='off')
      plt.setp(ax.get_xticklabels(), visible=False)
      ax.set_yticklabels([])

    # Prepare common legend by using existing handles and CommonUnitCellKeys
    handles, _ = sub_a.get_legend_handles_labels()
    common_key_lines = CommonUnitCellKey.common_lines_of(legend_keys)
    if len(info_list) == 1:
      labels = [k.lines() for k in legend_keys]
    else:
      unique = [not common for common in common_key_lines]
      labels = [k.prefix + k.sep + k.lines(unique) for k in legend_keys]
      if any(common_key_lines):
        handles.append(mpl.lines.Line2D([0], [0], alpha=0))  # empty handle
        labels.append(legend_keys[0].lines(common_key_lines))
    sub_key.legend(handles, labels, fontsize=text_ratio, labelspacing=1, loc=6)
    sub_key.axis('off')

    gsp.update(wspace=0)
    title = "Unit cell distribution" if title is None else title
    title += " (%d xtals)" % total
    title += " %s" % extra_title if extra_title else ""
    fig.suptitle(title)

    if not self.interactive:
      image_fname = image_fname or "ucell_tmp.png"
      fig.set_size_inches(xsize*1.05+.5, ysize*.95)
      fig.savefig(image_fname, bbox_inches='tight', dpi=100)
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

    AA = r"a-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(a),
                                        flex.mean_and_variance(a).unweighted_sample_standard_deviation())
    BB = r"b-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(b),
                                        flex.mean_and_variance(b).unweighted_sample_standard_deviation())
    CC = r"c-edge (%.2f +/- %.2f $\AA$)" % (flex.mean(c),
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
