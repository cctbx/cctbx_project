
"""
Classes for display of MolProbity validation statistics for multi-model PDB
files, used in GUI for phenix.ensemble_refinement.
"""

from __future__ import division
from libtbx.utils import Sorry
from wxtbx import plots, app
from mmtbx.command_line import validation_summary
from wxGUI2 import Base
import wx
import sys

class ensemble_validation_plot (plots.histogram) :
  def show_plot (self,
      values,
      as_histogram=False,
      n_bins=20,
      reference_value=None,
      title=None) :
    if (as_histogram) :
      self.show_histogram(
        data=values,
        n_bins=n_bins,
        reference_value=reference_value,
        x_label="Model Score",
        y_label="Number of Models",
        title=title)
    else :
      x = range(1, len(values) + 1)
      self.figure.clear()
      p = self.figure.add_subplot(111)
      ax = p.plot(x, values, '-^', color=(0.0,0.5,1.0))
      if (reference_value is not None) :
        p.axhline(reference_value, color='red')
      p.set_xlabel("Model Number")
      p.set_ylabel("Score")
      if (title is not None) :
        p.set_title(title)
      self.canvas.draw()

class ensemble_validation_panel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Plot controls
    box1 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box1)
    txt1 = wx.StaticText(self, label="Show statistic:")
    style = wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL
    box1.Add(txt1, 0, style , 5)
    self.stats_menu = wx.Choice(self,
      choices=validation_summary.molprobity_stat_labels)
    self.Bind(wx.EVT_CHOICE, self.OnSelectPlot, self.stats_menu)
    box1.Add(self.stats_menu, 0, style, 5)
    self.hist_box = wx.CheckBox(self, label="Display as histogram")
    box1.Add(self.hist_box, 0, style, 5)
    self.Bind(wx.EVT_CHECKBOX, self.OnSelectPlot, self.hist_box)
    n_bins_txt = wx.StaticText(self, label="Number of bins")
    box1.AddSpacer(40)
    box1.Add(n_bins_txt, 0, style, 5)
    self.n_bins_ctrl = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER,
                                   name="n_bins", size=(60, -1))
    self.n_bins_ctrl.SetValue("20")
    box1.Add(self.n_bins_ctrl, 0, style, 5)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnSelectPlot, self.n_bins_ctrl)
    self.n_bins_ctrl.Disable()

    # Data display
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box2)
    min_label = Base.BoldText(self, 'Minimum:')
    self.min_value = Base.PlainText(self, '')
    max_label = Base.BoldText(self, 'Maximum:')
    self.max_value = Base.PlainText(self, '')
    mean_label = Base.BoldText(self, 'Mean:')
    self.mean_value = Base.PlainText(self, '')
    self.value_ctrls = (self.min_value, self.max_value, self.mean_value)
    style = wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL
    for label, ctrl in zip([min_label, max_label, mean_label],self.value_ctrls):
      box2.Add(label, 0, style, 5)
      box2.Add(ctrl, 0, style, 5)
      box2.AddSpacer(40)
    self.plot = ensemble_validation_plot(
      parent=self,
      transparent=False)
    sizer.Add(self.plot, 1, wx.EXPAND|wx.ALL, 0)
    self.ensemble = None

  def set_ensemble (self, ensemble) :
    assert (type(ensemble).__name__ == 'ensemble')
    self.ensemble = ensemble

  def OnSelectPlot (self, event) :

    # plot line graph or histogram
    n_bins = 20
    as_histogram = self.hist_box.GetValue()
    if (as_histogram):
      self.n_bins_ctrl.Enable()
      n_bins = self.n_bins_ctrl.GetValue()
      try:
        n_bins = int(n_bins)
      except:
        raise Sorry('Please enter an integer for the number of histogram bins.')
    else:
      self.n_bins_ctrl.Disable()

    # check that n_bins is reasonable
    min_n_bins = 1
    max_n_bins = 1000
    if (n_bins < min_n_bins):
      n_bins = min_n_bins
    if (n_bins > max_n_bins):
      n_bins = max_n_bins
    self.n_bins_ctrl.SetValue(str(n_bins))

    # get attribute name and values for data
    label = self.stats_menu.GetStringSelection()
    i_label = validation_summary.molprobity_stat_labels.index(label)
    stat = self.ensemble.__slots__[i_label]
    values = getattr(self.ensemble, stat)
    if (len(values) == 0) or (values.count(None) == len(values)) :
      self.plot.figure.clear()
      [ ctrl.SetLabel("") for ctrl in self.value_ctrls ]
      return
    mean = sum(values) / len(values)
    dist = [ min(values), max(values), mean ]
    [ ctrl.SetLabel("%.3f" % x) for x, ctrl in zip(dist, self.value_ctrls) ]
    self.plot.show_plot(
      values=values,
      as_histogram=as_histogram,
      n_bins=n_bins,
      reference_value=mean,
      title=label)
    self.Refresh()

if (__name__ == "__main__") :
  result = validation_summary.run(sys.argv[1:])
  if (type(result).__name__ != 'ensemble') :
    raise Sorry("Not an ensemble, graphics not available.")
  app = app.CCTBXApp(0)
  frame = wx.Frame(None, -1, "Ensemble validation")
  szr = wx.BoxSizer(wx.VERTICAL)
  panel = ensemble_validation_panel(frame)
  panel.set_ensemble(result)
  szr.Add(panel, 1, wx.EXPAND|wx.ALL)
  frame.SetSizer(szr)
  szr.Fit(panel)
  frame.Fit()
  panel.OnSelectPlot(None)
  frame.Show()
  app.MainLoop()
