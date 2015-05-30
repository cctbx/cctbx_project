
"""
wxPython-specific code for displaying EM-Ringer plots.  The actual Matplotlib
plotting code lives in mmtbx.ringer.em_scoring.
"""

from __future__ import division
import wxtbx.plots
import wx

class peaks_plot (wxtbx.plots.plot_container) :
  pass

class emringer_plot_frame (wxtbx.plots.plot_frame) :
  """Base class for all frames enclosing plots."""
  def __init__ (self, *args, **kwds) :
    self._result = None
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)

  def SetResult (self, result) :
    self._result = result
    self._load_result()

  def _load_result (self) :
    raise NotImplementedError()

class peaks_plot_frame (emringer_plot_frame) :
  def create_plot_panel (self) :
    return peaks_plot(self)

  def draw_top_panel (self) :
    wxtbx.plots.plot_frame.draw_top_panel(self)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self.top_panel.SetSizer(szr)
    label = wx.StaticText(self.top_panel, -1, "Density threshold:")
    szr.Add(label, 0, wx.ALL, 5)
    self.threshold_choice = wx.Choice(self.top_panel, -1, size=(160,-1))
    szr.Add(self.threshold_choice, 0, wx.ALL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnSetThreshold, self.threshold_choice)

  def _load_result (self) :
    self.threshold_choice.SetItems(["%.3f"%x for x in self._result.thresholds])
    self._set_threshold(0)

  def OnSetThreshold (self, evt) :
    self._set_threshold(evt.GetEventObject().GetSelection())

  def _set_threshold (self, i) :
    assert (self._result is not None)
    self._result.draw_wx_peaks_plot(
      plot=self.plot_panel,
      i_threshold=i)
    self.plot_panel.canvas.draw()

class threshold_plot (wxtbx.plots.plot_container) :
  pass

class threshold_plot_frame (emringer_plot_frame) :
  def create_plot_panel (self) :
    return threshold_plot(self)

  def _load_result (self) :
    assert (self._result is not None)
    self._result.draw_wx_progression_plot(plot=self.plot_panel)
    self.plot_panel.canvas.draw()
