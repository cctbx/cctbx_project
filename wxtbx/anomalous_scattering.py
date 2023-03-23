
from __future__ import absolute_import, division, print_function
import wxtbx.plots
from wxtbx.phil_controls import strctrl, floatctrl
from libtbx.utils import Sorry
import wx
from six.moves import range

class AnomPlot(wxtbx.plots.plot_container):
  def __init__(self, *args, **kwds):
    wxtbx.plots.plot_container.__init__(self, *args, **kwds)
    self.canvas.mpl_connect('motion_notify_event', self.parent.OnHover)

  def show_plot(self,
      elements,
      range_type="wavelength",
      range_values=(0.75,2.5),
      n_points=100,
      include_fp=True,
      table="sasaki"):
    assert (range_type in ["wavelength", "energy"])
    assert (table in ["sasaki", "henke"])
    assert (n_points >= 10)
    x_min, x_max = range_values
    assert (x_max > x_min)
    if (range_type == "wavelength"):
      x_label = "Wavelength (Angstroms)"
    else :
      x_label= "Energy (eV)"
    table_module = None
    if (table == "sasaki"):
      from cctbx.eltbx import sasaki as table_module
    elif (table == "henke"):
      from cctbx.eltbx import henke as table_module
    assert (table_module is not None)
    increment = (x_max - x_min) / n_points
    x_points = [ x_min + increment*x for x in range(n_points+1) ]
    self.figure.clear()
    ax = self.figure.add_subplot(111)
    labels = []
    for element in elements :
      try :
        scatt_table = table_module.table(element.upper())
      except ValueError as e :
        raise Sorry(str(e))
      fp = []
      fdp = []
      for x in x_points :
        if (range_type == "wavelength"):
          scatt = scatt_table.at_angstrom(x)
        else :
          scatt = scatt_table.at_ev(x)
        fp.append(scatt.fp())
        fdp.append(scatt.fdp())
      fdp_line = ax.plot(x_points, fdp)[0]
      ax.plot(x_points, fp, linestyle='--', color=fdp_line.get_color())
      labels.append("%s f''" % scatt_table.label())
      labels.append("%s f'" % scatt_table.label())
    self.figure.legend(ax.lines, labels)
    ax.set_xlabel(x_label)
    ax.set_ylabel("e-")
    self.canvas.draw()

  def OnHover(self, mpl_event):
    (xdata, ydata) = (mpl_event.xdata, mpl_event.ydata)

class AnomPlotFrame(wxtbx.plots.plot_frame):
  def draw_top_panel(self):
    self.top_panel = ControlPanel(parent=self, style=wx.SUNKEN_BORDER)
    self.statusbar = self.CreateStatusBar()

  def create_plot_panel(self):
    return AnomPlot(parent=self, figure_size=(12,6))

  def update_plot(self, *args, **kwds):
    self.plot_panel.show_plot(*args, **kwds)
    self.Refresh()

  def OnHover(self, mpl_event):
    (xdata, ydata) = (mpl_event.xdata, mpl_event.ydata)
    if (xdata is None) or (ydata is None):
      self.statusbar.SetStatusText("")
    else :
      wavelength = energy = 0
      if (xdata < 10):
        wavelength = xdata
        if (xdata != 0):
          energy = 12398 / xdata
      else :
        energy = xdata
        if (xdata != 0):
          wavelength = 12398 / energy
      status_label = "Wavelength=%.3f A   Energy=%.3f eV  e-=%.3f" % (
        wavelength, energy, ydata)
      self.statusbar.SetStatusText(status_label)

class ControlPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr1 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr1)
    txt1 = wx.StaticText(self, -1, "Elements:")
    szr1.Add(txt1, 0, wx.ALL, 5)
    self.elem_ctrl = strctrl.StrCtrl(
      parent=self,
      name="Elements",
      size=(200,-1))
    self.elem_ctrl.SetOptional(False)
    szr1.Add(self.elem_ctrl, 0, wx.ALL, 5)
    szr1.Add(wx.StaticText(self, -1, "Table:"), 0, wx.ALL, 5)
    self.table_choice = wx.Choice(parent=self,
      choices=["Sasaki", "Henke"])
    szr1.Add(self.table_choice, 0, wx.ALL, 5)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr2)
    self.range_choice = wx.Choice(parent=self,
      choices=["Wavelength", "Energy"])
    szr2.Add(self.range_choice, 0, wx.ALL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnChangeRange, self.range_choice)
    self.range_min = floatctrl.FloatCtrl(
      parent=self, name="Range minimum")
    self.range_min.SetMin(0.5)
    self.range_min.SetMax(5.0)
    self.range_min.SetOptional(False)
    self.range_max = floatctrl.FloatCtrl(
      parent=self, name="Range maximum")
    self.range_max.SetMin(0.5)
    self.range_max.SetMax(5.0)
    self.range_max.SetOptional(False)
    szr2.Add(wx.StaticText(self, -1, "min:"), 0, wx.ALL, 5)
    szr2.Add(self.range_min, 0, wx.ALL, 5)
    szr2.Add(wx.StaticText(self, -1, "max:"), 0, wx.ALL, 5)
    szr2.Add(self.range_max, 0, wx.ALL, 5)
    btn = wx.Button(self, -1, "Redraw plot")
    self.Bind(wx.EVT_BUTTON, self.OnUpdatePlot, btn)
    szr2.Add(btn, 0, wx.ALL, 5)

  def OnChangeRange(self, evt):
    range_type = self.range_choice.GetStringSelection()
    if (range_type == "Wavelength"):
      self.range_min.SetMin(0.75)
      self.range_min.SetMax(5.0)
      self.range_max.SetMin(0.75)
      self.range_max.SetMax(5.0)
    else :
      self.range_min.SetMin(2500)
      self.range_min.SetMax(16000)
      self.range_max.SetMin(2500)
      self.range_max.SetMax(16000)

  def OnUpdatePlot(self, evt):
    elements_str = self.elem_ctrl.GetPhilValue()
    elements = elements_str.replace(",", " ").split()
    from cctbx.eltbx import chemical_elements
    allowed_list = chemical_elements.proper_upper_list()
    for elem in elements :
      if (len(elem) > 2):
        raise Sorry(
          "You must enter standard chemical element symbols (H, Zn, etc.).")
      if (not elem.upper() in allowed_list):
        raise Sorry("The symbol '%s' is not a recognized chemical element." %
          elem)
    range_min = self.range_min.GetPhilValue()
    range_max = self.range_max.GetPhilValue()
    range_type = self.range_choice.GetStringSelection().lower()
    table = self.table_choice.GetStringSelection().lower()
    if (range_max <= range_min):
      raise Sorry("The range maximum must be greater than the minimum")
    self.GetParent().update_plot(
      elements=elements,
      range_type=range_type,
      range_values=(range_min,range_max),
      table=table,
      include_fp=True)
