
"""
Classes for display of MolProbity validation statistics for multi-model PDB
files, used in GUI for phenix.ensemble_refinement.
"""

from __future__ import division
from wxtbx import plots, app
from mmtbx.command_line import validation_summary
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
        x_label="Number of models",
        y_label="Model score",
        title=title)
    else :
      x = range(1, len(values) + 1)
      self.figure.clear()
      p = self.figure.add_subplot(111)
      ax = p.plot(x, values, '-^', color=(0.0,0.5,1.0))
      if (reference_value is not None) :
        p.axhline(reference_value, color='red')
      p.set_xlabel("Model number")
      p.set_ylabel("Score")
      if (title is not None) :
        p.set_title(title)
      self.canvas.draw()

class ensemble_validation_panel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)
    box1 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box1)
    txt1 = wx.StaticText(self, label="Show statistic:")
    box1.Add(txt1, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.stats_menu = wx.Choice(self,
      choices=validation_summary.ensemble.__slot_labels__)
    self.Bind(wx.EVT_CHOICE, self.OnSelectPlot, self.stats_menu)
    box1.Add(self.stats_menu, 0,
      wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.hist_box = wx.CheckBox(self, label="Display as histogram")
    box1.Add(self.hist_box, 0,
      wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHECKBOX, self.OnSelectPlot, self.hist_box)
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(box2)
    txt2 = wx.StaticText(self, label="Minimum:")
    txt3 = wx.StaticText(self, label="Maximum:")
    txt4 = wx.StaticText(self, label="Mean:")
    val1 = wx.TextCtrl(self, name="Minimum", size=(80,-1))
    val2 = wx.TextCtrl(self, name="Minimum", size=(80,-1))
    val3 = wx.TextCtrl(self, name="Minimum", size=(80,-1))
    self.value_ctrls = (val1, val2, val3)
    box2.Add(txt2, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(val1, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(txt3, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(val2, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(txt4, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(val3, 0, wx.LEFT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.plot = ensemble_validation_plot(
      parent=self,
      transparent=False)
    sizer.Add(self.plot, 1, wx.EXPAND|wx.ALL, 0)
    self.ensemble = None

  def set_ensemble (self, ensemble) :
    assert (type(result).__name__ == 'ensemble')
    self.ensemble = ensemble

  def OnSelectPlot (self, event) :
    as_histogram = self.hist_box.GetValue()
    label = self.stats_menu.GetStringSelection()
    i_label = validation_summary.ensemble.__slot_labels__.index(label)
    stat = validation_summary.ensemble.__slots__[i_label]
    values = getattr(self.ensemble, stat)
    if (len(values) == 0) or (values.count(None) == len(values)) :
      self.plot.figure.clear()
      [ ctrl.SetValue("") for ctrl in self.value_ctrls ]
      return
    mean = sum(values) / len(values)
    dist = [ min(values), max(values), mean ]
    [ ctrl.SetValue("%g" % x) for x, ctrl in zip(dist, self.value_ctrls) ]
    self.plot.show_plot(
      values=values,
      as_histogram=as_histogram,
      n_bins=20,
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
