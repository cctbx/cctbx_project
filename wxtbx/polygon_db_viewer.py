
from wxtbx import plots
from mmtbx import polygon
from libtbx.utils import Sorry
from libtbx import group_args
import wx

STD_FLAGS = wx.ALL|wx.ALIGN_CENTER_VERTICAL

all_keys = polygon.keys_to_show + polygon.other_numerical_keys
all_captions = polygon.key_captions + polygon.other_captions

class ConfigFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.statusbar.SetStatusText("Reference: Urzhumtseva et al. (2009) "+
      "Acta Cryst. D65:297-300.")
    p = wx.Panel(self, -1)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr.Add(p, 1, wx.EXPAND)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    p.SetSizer(szr2)
    box1 = wx.StaticBox(p, -1, "Plot correlation between statistics")
    bszr = wx.StaticBoxSizer(box1, wx.VERTICAL)
    szr2.Add(bszr, 0, STD_FLAGS, 5)
    grid = wx.FlexGridSizer(cols=5)
    bszr.Add(grid, 1, wx.EXPAND)
    grid.Add(wx.StaticText(p, -1, "X axis:"), 0, STD_FLAGS, 5)
    self.x_chooser = wx.Choice(p, -1, choices=all_captions)
    grid.Add(self.x_chooser, 0, STD_FLAGS|wx.EXPAND, 5)
    grid.Add((10,1))
    grid.Add(wx.StaticText(p, -1, "Y axis:"), 0, STD_FLAGS, 5)
    self.y_chooser = wx.Choice(p, -1, choices=all_captions)
    grid.Add(self.y_chooser, 0, STD_FLAGS|wx.EXPAND, 5)
    grid.Add(wx.StaticText(p, -1, "Limits:"), 0, STD_FLAGS, 5)
    inner_szr1 = wx.BoxSizer(wx.HORIZONTAL)
    grid.Add(inner_szr1, 0)
    self.x_min = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr1.Add(self.x_min, 0, STD_FLAGS, 5)
    inner_szr1.Add(wx.StaticText(p, -1, "to"), 0, STD_FLAGS, 5)
    self.x_max = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr1.Add(self.x_max, 0, STD_FLAGS, 5)
    grid.Add((10,1))
    grid.Add(wx.StaticText(p, -1, "Limits:"), 0, STD_FLAGS, 5)
    inner_szr2 = wx.BoxSizer(wx.HORIZONTAL)
    grid.Add(inner_szr2, 0)
    self.y_min = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr2.Add(self.y_min, 0, STD_FLAGS, 5)
    inner_szr2.Add(wx.StaticText(p, -1, "to"), 0, STD_FLAGS, 5)
    self.y_max = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr2.Add(self.y_max, 0, STD_FLAGS, 5)
    plot_btn = wx.Button(p, -1, "Show plot")
    self.Bind(wx.EVT_BUTTON, self.OnPlotCorr, plot_btn)
    bszr.Add(plot_btn, 0, wx.ALIGN_RIGHT|wx.ALL, 5)
    # histogram controls
    box2 = wx.StaticBox(p, -1, "Plot histogram")
    bszr2 = wx.StaticBoxSizer(box2, wx.VERTICAL)
    szr2.Add(bszr2, 0, STD_FLAGS|wx.EXPAND, 5)
    grid2 = wx.FlexGridSizer(cols=2)
    bszr2.Add(grid2, 1, wx.EXPAND)
    grid2.Add(wx.StaticText(p, -1, "Statistic:"), 0, STD_FLAGS, 5)
    self.h_chooser = wx.Choice(p, -1, choices=all_captions)
    grid2.Add(self.h_chooser, 0, STD_FLAGS|wx.EXPAND, 5)
    grid2.Add(wx.StaticText(p, -1, "Limits:"), 0, STD_FLAGS, 5)
    inner_szr4 = wx.BoxSizer(wx.HORIZONTAL)
    grid2.Add(inner_szr4)
    self.h_min = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr4.Add(self.h_min, 0, STD_FLAGS, 5)
    inner_szr4.Add(wx.StaticText(p, -1, "to"), 0, STD_FLAGS, 5)
    self.h_max = wx.TextCtrl(p, -1, "", size=(80, -1))
    inner_szr4.Add(self.h_max, 0, STD_FLAGS, 5)
    grid2.Add(wx.StaticText(p, -1, "Resolution range:"), 0, STD_FLAGS, 5)
    inner_szr2 = wx.BoxSizer(wx.HORIZONTAL)
    grid2.Add(inner_szr2)
    self.d_min = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr2.Add(self.d_min, 0, STD_FLAGS, 5)
    inner_szr2.Add(wx.StaticText(p, -1, "to"), 0, STD_FLAGS, 5)
    self.d_max = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr2.Add(self.d_max, 0, STD_FLAGS, 5)
    grid2.Add(wx.StaticText(p, -1, "Number of bins:"), 0, STD_FLAGS, 5)
    inner_szr3 = wx.BoxSizer(wx.HORIZONTAL)
    grid2.Add(inner_szr3, 0, STD_FLAGS, 5)
    self.n_bins = wx.TextCtrl(p, -1, "20", size=(80,-1))
    inner_szr3.Add(self.n_bins, 0, STD_FLAGS, 5)
    inner_szr3.Add(wx.StaticText(p, -1, "Reference value:"), 0, STD_FLAGS, 5)
    self.ref_value = wx.TextCtrl(p, -1, "", size=(80,-1))
    inner_szr3.Add(self.ref_value, 0, STD_FLAGS, 5)
    hist_btn = wx.Button(p, -1, "Show histogram")
    self.Bind(wx.EVT_BUTTON, self.OnPlotHist, hist_btn)
    bszr2.Add(hist_btn, 0, wx.ALIGN_RIGHT|wx.ALL, 5)
    szr2.Layout()
    szr.Fit(p)
    self.Fit()
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self._db = None

  def get_db (self) :
    if (self._db is None) :
      self._db = polygon.load_db()
    return self._db

  def OnPlotCorr (self, evt) :
    x_key = all_keys[self.x_chooser.GetSelection()]
    y_key = all_keys[self.y_chooser.GetSelection()]
    if (x_key == y_key) :
      raise Sorry("X and Y selections are the same statistic!")
    db = self.get_db()
    x_stats = db[x_key]
    y_stats = db[y_key]
    limits = group_args(x_min=None, x_max=None, y_min=None, y_max=None)
    for lim in ["x_min", "x_max", "y_min", "y_max"] :
      lim_txt = getattr(self, lim).GetValue()
      if (lim_txt != "") :
        try :
          lim_val = float(lim_txt)
        except ValueError :
          raise Sorry("Invalid value for limit - must be a decimal number.")
        else :
          setattr(limits, lim, lim_val)
    x_values = []
    y_values = []
    for x_, y_ in zip(x_stats, y_stats) :
      try :
        x = float(x_)
        y = float(y_)
      except ValueError : pass
      else :
        if   (limits.x_min is not None) and (x < limits.x_min) : continue
        elif (limits.x_max is not None) and (x > limits.x_max) : continue
        elif (limits.y_min is not None) and (y < limits.y_min) : continue
        elif (limits.y_max is not None) and (x > limits.y_max) : continue
        else :
          x_values.append(x)
          y_values.append(y)
    assert (len(x_values) == len(y_values))
    if (len(x_values) == 0) :
      raise Sorry("No points within specified limits.")
    plt = CorrPlotFrame(self)
    plt.set_plot(x_values, y_values, self.x_chooser.GetStringSelection(),
      self.y_chooser.GetStringSelection())
    plt.Show()

  def OnPlotHist (self, evt) :
    h_key = all_keys[self.h_chooser.GetSelection()]
    db = self.get_db()
    h_values = db[h_key]
    limits = group_args(d_min=None, d_max=None, h_min=None, h_max=None)
    for lim in ["d_min", "d_max", "h_min", "h_max"] :
      lim_txt = getattr(self, lim).GetValue()
      if (lim_txt != "") :
        try :
          lim_val = float(lim_txt)
        except ValueError :
          raise Sorry("Invalid resolution '%s' - must be a decimal number." %
            lim_txt)
        else :
          setattr(limits, lim, lim_val)
    reference_value = None
    reference_value_txt = self.ref_value.GetValue()
    if (reference_value_txt != "") :
      try :
        reference_value = float(reference_value_txt)
      except ValueError :
        raise Sorry("Invalid reference value '%s' - must be a decimal number."%
          reference_value_txt)
    try :
      n_bins = float(self.n_bins.GetValue())
    except ValueError :
      raise Sorry("Number of bins must be a decimal number.")
    data = []
    print limits.h_min, limits.h_max
    for (d_, v_) in zip(db['high_resolution'], db[h_key]) :
      try :
        d_min = float(d_)
        value = float(v_)
      except ValueError : pass
      else :
        if   (limits.d_min is not None) and (d_min < limits.d_min) : continue
        elif (limits.d_max is not None) and (d_min > limits.d_max) : continue
        elif (limits.h_min is not None) and (value < limits.h_min) : continue
        elif (limits.h_max is not None) and (value > limits.h_max) : continue
        else :
          data.append(value)
    if (len(data) == 0) :
      raise Sorry("No points within specified limits.")
    plt = HistogramFrame(self)
    plt.show_histogram(
      data=data,
      n_bins=n_bins,
      reference_value=reference_value,
      xlabel=self.h_chooser.GetStringSelection())
    plt.Show()

  def OnShowTable (self, evt) :
    pass

  def OnClose (self, evt) :
    self.Destroy()

  def OnDestroy (self, evt) :
    parent = self.GetParent()
    if (parent is not None) :
      parent.plot_frame = None

class CorrPlot (plots.plot_container) :
  def set_plot (self, x, y, x_label, y_label) :
    from scitbx.array_family import flex
    self.figure.clear()
    ax = self.figure.add_subplot(111)
    ax.plot(x, y, '.')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, color="0.5")
    c = flex.linear_correlation(flex.double(x), flex.double(y))
    cc = c.coefficient()
    ax.set_title("%s vs. %s (%d structures, CC = %.3f)" % (x_label, y_label,
      len(x), cc))
    self.canvas.draw()
    #self.parent.statusbar.SetStatusText("Correlation coefficient (CC): %.3f" %
    #  cc)
    self.parent.Refresh()

class CorrPlotFrame (plots.plot_frame) :
  show_controls_default = False
  def create_plot_panel (self) :
    self.statusbar = self.CreateStatusBar()
    return CorrPlot(self, figure_size=(12,8))

  def set_plot (self, *args, **kwds) :
    self.plot_panel.set_plot(*args, **kwds)

class HistogramPlot (plots.histogram) :
  def show_histogram (self, data, n_bins, reference_value, xlabel) :
    from scitbx.array_family import flex
    mean = flex.mean(flex.double(data))
    p = plots.histogram.show_histogram(self,
      data=data,
      n_bins=n_bins,
      reference_value=reference_value,
      draw_now=False)
    p.axvline(mean, color='g', linewidth=2)
    p.get_axes().set_ylabel("Number of structures")
    p.get_axes().set_xlabel(xlabel)
    if (reference_value is not None) :
      labels = ["Reference value", "Mean (%.3f)" % mean]
    else :
      labels = ["Mean (%.3f)" % mean]
    p.set_title("Distribution of %s (%d structures)" % (xlabel, len(data)))
    self.figure.legend(p.lines, labels)
    self.canvas.draw()
    self.parent.Refresh()

class HistogramFrame (plots.plot_frame) :
  show_controls_default = False
  def create_plot_panel (self) :
    self.statusbar = self.CreateStatusBar()
    return HistogramPlot(self, figure_size=(12,8))

  def show_histogram (self, *args, **kwds) :
    self.plot_panel.show_histogram(*args, **kwds)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = ConfigFrame(None, -1, "PDB statistics from POLYGON database")
  frame.Show()
  app.MainLoop()
