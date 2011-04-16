
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
    szr2.Layout()
    szr.Fit(p)
    self.Fit()
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
    pass

  def OnShowTable (self, evt) :
    pass

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
    ax.set_title("%s vs. %s (CC = %.3f)" % (x_label, y_label, cc))
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

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = ConfigFrame(None, -1, "PDB statistics from POLYGON database")
  frame.Show()
  app.MainLoop()
