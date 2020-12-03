
"""
Base classes for visualization of MolProbity analysis using matplotlib.
"""

from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate
from six.moves import filterfalse
from six.moves import zip
from six.moves import range

class rotarama_plot_mixin(object):
  extent = [0, 360, 0, 360]
  def __init__(self):
    assert hasattr(self, "figure")
    self._points = []
    self._xyz = [] # only used by Phenix GUI (not offline plotting)
    self.plot = self.figure.add_subplot(111)
    self.plot.set_position([0.1, 0.1, 0.85, 0.85])

  def draw_plot(self,
                 stats,
                 title,
                 points=None,
                 show_labels=True,
                 colormap='jet',
                 contours=None,
                 xyz=None,
                 extent=None,
                 y_marks=None,
                 markerfacecolor="white",
                 markeredgecolor="black",
                 show_filling=True,
                 markersize=5,
                 point_style='bo'):
    # points = [(x,y,label, isoutlier(bool)), (), ...]
    import matplotlib.cm
    self._points = []
    self._xyz = []
    cm = getattr(matplotlib.cm, colormap)
    self.plot.clear()
    if (extent is None):
      extent = self.extent
    else :
      assert (len(extent) == 4)
    if show_filling:
      im = self.plot.imshow(stats, origin="lower", cmap=cm, extent=extent)
      # This code will draw bar with actual numbers that are represented by
      # color (sloppily over the axis values)
      # from mpl_toolkits.axes_grid1 import make_axes_locatable
      # divider = make_axes_locatable(self.plot)
      # cax = divider.append_axes('bottom', size='5%', pad=0.05)
      # self.figure.colorbar(im, cax=cax, orientation='horizontal')
    if (contours is not None):
      self.plot.contour(stats, contours,
        origin="lower",
        colors='k',
        extent=extent,
        zorder=9)
    if (y_marks is None):
      self.set_labels()
    else :
      self.set_labels(y_marks=y_marks)
    self.plot.set_title(title)
    if points is not None:
      if xyz is not None: assert (len(xyz) == len(points))
      out = list(filter(lambda x: x[3], points))
      out_columns = list(zip(*out))
      # ^^^^ is doing e.g. this:
      # >>> l = [(1,2), (3,4), (8,9)]
      # >>> zip(*l)
      # [(1, 3, 8), (2, 4, 9)]
      non_out = list(filterfalse(lambda x: x[3], points))
      non_out_columns = list(zip(*non_out))
      if len(out) > 0:
        self.plot.plot(tuple(out_columns[0]), tuple(out_columns[1]), point_style,
          markerfacecolor='red', markersize=markersize,
          markeredgecolor=markeredgecolor)
        if show_labels :
          for x, y, label, _ in out:
            self.plot.text(x, y, label, color='black')
      if len(non_out) > 0:
        # make non-outliers smaller and translucent if there are many points
        if len(non_out) > 2500:
          markersize = 3
        self.plot.plot(non_out_columns[0], non_out_columns[1], point_style,
          markerfacecolor=markerfacecolor, markersize=markersize,
          markeredgecolor=markeredgecolor)
    self.canvas.draw()

class residue_bin(slots_getstate_setstate):
  __slots__ = ["residues", "marks", "labels"]
  def __init__(self):
    self.residues = []
    self.marks = []
    self.labels = []

  def add_residue(self, residue):
    self.residues.append(residue)
    n_res = len(self.residues)
    if (n_res % 10 == 0):
      self.marks.append(n_res)
      if (residue is not None):
        self.labels.append(residue.residue_group_id_str())
      else :
        self.labels.append("")

  def add_empty(self, n):
    for i in range(n):
      self.add_residue(None)

  def n_res(self):
    return len(self.residues)

  def get_residue_range(self):
    bin_start = bin_end = None
    i = 0
    while (i < len(self.residues)):
      if (self.residues[i] is not None):
        bin_start = self.residues[i].residue_group_id_str()
        break
      i += 1
    i = len(self.residues) - 1
    while (i >= 0):
      if (self.residues[i] is not None):
        bin_end = self.residues[i].residue_group_id_str()
        break
      i -= 1
    return "%s - %s" % (bin_start, bin_end)

  def x_values(self):
    return list(range(len(self.residues)))

  def get_selected(self, index):
    return self.residues[index]

  def get_real_space_plot_values(self):
    import numpy
    y = []
    for residue in self.residues :
      if (residue is not None):
        y.append(residue.get_real_space_plot_values())
      else :
        y.append([numpy.NaN] * 4)
    return numpy.array(y).transpose()

  def get_outlier_plot_values(self):
    import numpy
    y = []
    for residue in self.residues :
      if (residue is not None):
        y.append(residue.get_outlier_plot_values())
      else :
        y.append([numpy.NaN] * 4)
    return numpy.array(y).transpose()

class residue_binner(object):
  def __init__(self, res_list, bin_size=100, one_chain_per_bin=False):
    self.bins = []
    last_chain = last_resseq = None
    for i, residue in enumerate(res_list):
      if (len(self.bins) == 0) or (self.bins[-1].n_res() == bin_size):
        self.bins.append(residue_bin())
      chain = residue.chain_id
      resseq = residue.resseq_as_int()
      if (last_chain is not None):
        # FIXME needs to take icode into account!
        if (chain != last_chain) or (resseq > (last_resseq + 10)):
          if ((chain != last_chain and one_chain_per_bin) or
              (self.bins[-1].n_res() > (bin_size - 20))):
            self.bins.append(residue_bin())
          else :
            self.bins[-1].add_empty(10)
        elif (resseq > (last_resseq + 1)) and (self.bins[-1].n_res() > 0):
          gap = resseq - (last_resseq + 1)
          i = 0
          while (i < gap) and (self.bins[-1].n_res() < bin_size):
            self.bins[-1].add_empty(1)
            i += 1
      self.bins[-1].add_residue(residue)
      last_chain = chain
      last_resseq = resseq

  def get_bin(self, i_bin):
    return self.bins[i_bin]

  def get_ranges(self):
    return [ bin.get_residue_range() for bin in self.bins ]

class multi_criterion_plot_mixin(object):
  def __init__(self, binner, y_limits):
    self.binner = binner
    self.y_limits = y_limits
    self.disabled = False

  def plot_range(self, i_bin):
    if (self.disabled) : return
    # TODO: fix y-ticks, x-ticks width residue ID, add legend
    self.figure.clear()
    bin = self.binner.get_bin(i_bin)
    self._current_bin = bin
    x = bin.x_values()
    # b_iso, cc, fmodel, two_fofc
    y = bin.get_real_space_plot_values()
    (b_iso, cc, fmodel, two_fofc) = bin.get_real_space_plot_values()
    # electron density (observed and calculated)
    rho_min, rho_max = self.y_limits["rho"]
    plot_density = not ( (rho_min is None) and (rho_max is None) )
    top_position = [0.075, 0.625, 0.825, 0.325]
    bot_position = [0.075, 0.075, 0.825, 0.55]
    if (plot_density):
      n_plots = 2
      p1 = self.figure.add_subplot(n_plots, 1, n_plots)
      p1.set_position(top_position)
      p1.plot(x, y[2], "-", linewidth=1, color="g")
      p1.plot(x, y[3], "-", linewidth=1, color="k")
      p1.set_title("Multi-criterion validation")
      p1.set_ylabel("Density", fontproperties=self.get_font("label"))
      p1.set_xlim(-1, 101)
      p1.set_ylim(rho_min, rho_max)
      p1.xaxis.set_major_formatter(self.null_fmt)
    else:
      n_plots = 1
    # CC
    if (plot_density):
      p2 = self.figure.add_subplot(n_plots, 1, n_plots - 1, sharex=p1)
      p2.set_position(bot_position)
    else:
      p2 = self.figure.add_subplot(n_plots, 1, n_plots)
    p2.plot(x, y[1], "-", linewidth=1, color="b")
    p2.set_ylabel("Local real-space CC", fontproperties=self.get_font("label"),
      color='b')
    p2.set_xlim(-1, 101)
    cc_min, cc_max = self.y_limits["cc"]
    p2.set_ylim(cc_min, cc_max)
    p2.get_yticklabels()[-1].set_visible(False)
    # B_iso
    p3 = p2.twinx()
    if (plot_density):
      p3.set_position(bot_position)
    p3.plot(x, y[0], "-", linewidth=1, color="r")
    p3.set_ylabel("B-factor", fontproperties=self.get_font("label"),
      color='r')
    b_min, b_max = self.y_limits["b"]
    p3.set_ylim(b_min, b_max)
    p3.get_yticklabels()[-1].set_visible(False)
    # labels
    p2.set_xlabel("Residue", fontproperties=self.get_font("label"))
    p2.set_xticks(bin.marks)
    p2.set_xticklabels(bin.labels, fontproperties=self.get_font("basic"))
    p2.get_yticklabels()[-1].set_visible(False)
    # rama, rota, cbeta, clash
    y2 = bin.get_outlier_plot_values()
    p2.plot(x, y2[0] * 0.99, "o")
    p2.plot(x, y2[1] * 0.98, "^")
    p2.plot(x, y2[2] * 0.97, "s")
    p2.plot(x, y2[3] * 0.96, "d")
    legend_lines = p2.lines + p3.lines
    legend_labels = ["CC", "Ramachandran", "Rotamer", "C-beta", "Bad clash",
                     "B-factor"]
    if (plot_density):
      legend_lines += p1.lines
      legend_labels.extend(["Fc", "2mFo-DFc"])
    self.figure.legend(legend_lines, legend_labels,
                       prop=self.get_font("legend"))
    self.canvas.draw()
