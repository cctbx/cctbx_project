
"""
wxPython-specific code for displaying EM-Ringer plots.  The actual Matplotlib
plotting code lives in mmtbx.ringer.em_scoring.
"""

from __future__ import absolute_import, division, print_function
from wxtbx.phil_controls import floatctrl
from wxtbx import phil_controls
import wxtbx.plots
import wx

class peaks_plot(wxtbx.plots.plot_container):
  pass

class emringer_plot_frame(wxtbx.plots.plot_frame):
  """Base class for all frames enclosing plots."""
  def __init__(self, *args, **kwds):
    self._result = None
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)

  def SetResult(self, result):
    self._result = result
    self._load_result()

  def _load_result(self):
    raise NotImplementedError()

class peaks_plot_frame(emringer_plot_frame):
  def create_plot_panel(self):
    return peaks_plot(self)

  def draw_top_panel(self):
    wxtbx.plots.plot_frame.draw_top_panel(self)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self.top_panel.SetSizer(szr)
    label = wx.StaticText(self.top_panel, -1, "Density threshold:")
    szr.Add(label, 0, wx.ALL, 5)
    self.threshold_choice = wx.Choice(self.top_panel, -1, size=(160,-1))
    szr.Add(self.threshold_choice, 0, wx.ALL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnSetThreshold, self.threshold_choice)

  def _load_result(self):
    self.threshold_choice.SetItems(["%.3f"%x for x in self._result.thresholds])
    self._set_threshold(0)

  def OnSetThreshold(self, evt):
    self._set_threshold(evt.GetEventObject().GetSelection())

  def _set_threshold(self, i):
    assert (self._result is not None)
    self._result.draw_wx_peaks_plot(
      plot=self.plot_panel,
      i_threshold=i)
    self.plot_panel.canvas.draw()

class threshold_plot(wxtbx.plots.plot_container):
  pass

class threshold_plot_frame(emringer_plot_frame):
  def create_plot_panel(self):
    return threshold_plot(self)

  def _load_result(self):
    assert (self._result is not None)
    self._result.draw_wx_progression_plot(plot=self.plot_panel)
    self.plot_panel.canvas.draw()

class rolling_plot(wxtbx.plots.plot_container):
  pass

class rolling_plot_frame(emringer_plot_frame):
  def create_plot_panel(self):
    return rolling_plot(self)

  def draw_top_panel(self):
    wxtbx.plots.plot_frame.draw_top_panel(self)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self.top_panel.SetSizer(szr)
    label1 = wx.StaticText(self.top_panel, -1, "Chain ID:")
    szr.Add(label1, 0, wx.ALL, 5)
    self.chain_choice = wx.Choice(self.top_panel, -1, size=(160,-1))
    szr.Add(self.chain_choice, 0, wx.ALL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnSetChain, self.chain_choice)
    #label2 = wx.StaticText(self.top_panel, -1, "Density threshold:")
    #szr.Add(label2, 0, wx.ALL, 5)
    #self.threshold_entry = floatctrl.FloatCtrl(
    #  parent=self.top_panel,
    #  size=(80,-1),
    #  value=0,
    #  style=wx.TE_PROCESS_ENTER)
    #self.threshold_entry.SetMin(0)
    #self.threshold_entry.SetMax(1.0)
    #szr.Add(self.threshold_entry, 0, wx.ALL, 5)
    #self.Bind(phil_controls.EVT_PHIL_CONTROL, self.OnSetThreshold,
    #  self.threshold_entry)

  def _load_result(self):
    self.chain_choice.SetItems(self._result.chain_ids)
    self._redraw()

  def OnSetThreshold(self, evt):
    self._redraw()

  def OnSetChain(self, evt):
    self._redraw()

  def _redraw(self):
    chain_id = self.chain_choice.GetStringSelection()
    #threshold = self.threshold_entry.GetPhilValue()
    self._result.draw_wx_plot(
      plot=self.plot_panel,
      chain_id=chain_id)
#      threshold=threshold)
    self.plot_panel.canvas.draw()
