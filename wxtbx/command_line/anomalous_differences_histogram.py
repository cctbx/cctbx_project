# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function
import wxtbx.plots
from libtbx.utils import Sorry
import wx
import os
import sys

def display_file_info(file_name, obs_type="amplitudes", parent=None,
    n_bins=40, out=sys.stdout):
  from iotbx import file_reader
  hkl_in = file_reader.any_file(file_name, force_type="hkl")
  hkl_in.check_file_type("hkl")
  anom_data = []
  for array in hkl_in.file_server.miller_arrays :
    if (not array.anomalous_flag()) : continue
    if (array.is_xray_amplitude_array()) or (array.is_xray_intensity_array()):
      anom_data.append(array)
  if (len(anom_data) == 0):
    raise Sorry("No anomalous data arrays found.")
  x_label = "D_anom(F) / <F>"
  if (obs_type == "intensities"):
    x_label = "D_anom(I) / <I>"
  for array in anom_data :
    array_name = os.path.basename(file_name)+":"+array.info().label_string()
    frame = AnomHistPlotFrame(
      parent=parent,
      title="Bijvoet ratios for %s" % array_name,
      array=array,
      array_name=array_name)
    frame.Show()

class AnomHistPlotFrame(wxtbx.plots.plot_frame):
  def __init__(self, *args, **kwds):
    self._array = kwds.pop("array")
    self._array_name = kwds.pop("array_name")
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)
    self.recalculate()

  def create_plot_panel(self):
    return wxtbx.plots.histogram(parent=self, figure_size=(12,6))

  def draw_top_panel(self):
    self.top_panel = wx.Panel(parent=self, style=wx.RAISED_BORDER)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.top_panel.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr2, 1)
    self.meas_box = wx.CheckBox(self.top_panel, label="Measurable differences only")
    szr2.Add(self.meas_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.log_box = wx.CheckBox(self.top_panel, label="Log scale")
    szr2.Add(self.log_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr2.Add(wx.StaticText(self.top_panel, label="Data type:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.obs_type = wx.Choice(self.top_panel,
      choices=["amplitudes","intensities"])
    szr2.Add(self.obs_type, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr2.Add(wx.StaticText(self.top_panel, label="# of bins:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.n_bins = wx.SpinCtrl(self.top_panel, -1)
    self.n_bins.SetValue(40)
    szr2.Add(self.n_bins, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHECKBOX, lambda evt: self.recalculate(), self.meas_box)
    self.Bind(wx.EVT_CHECKBOX, lambda evt: self.recalculate(), self.log_box)
    self.Bind(wx.EVT_CHOICE, lambda evt: self.recalculate(), self.obs_type)
    self.Bind(wx.EVT_SPINCTRL, lambda evt: self.recalculate(), self.n_bins)

  def recalculate(self):
    from scitbx.array_family import flex
    measurable_only = self.meas_box.GetValue()
    obs_type = self.obs_type.GetStringSelection()
    n_bins = self.n_bins.GetValue()
    d_ano_rel = self._array.bijvoet_ratios(
      obs_type=obs_type,
      measurable_only=measurable_only)
    hist = flex.histogram(d_ano_rel, n_slots=n_bins)
    hist.show(f=sys.stdout)
    if (obs_type == "intensities"):
      x_label = "D_anom(I[hkl]) / I_mean[hkl]"
    else :
      x_label = "D_anom(F[hkl]) / F_mean[hkl]"
    self.show_histogram(
      data=list(d_ano_rel),
      n_bins=n_bins,
      x_label=x_label,
      y_label="# hkl",
      title="Bijvoet ratios for %s" % self._array_name,
      log_scale=self.log_box.GetValue())

  def show_histogram(self, *args, **kwds):
    self.plot_panel.show_histogram(*args, **kwds)
    self.Refresh()

def run(args=(), params=None, out=sys.stdout):
  import wxtbx.app
  app = wxtbx.app.CCTBXApp(0)
  if (len(args) == 0):
    from wxtbx.command_line.inspect_r_free_flags import ask_for_file
    file_name = ask_for_file(parent=None)
  else :
    file_name = args[0]
  display_file_info(file_name)
  app.MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])
