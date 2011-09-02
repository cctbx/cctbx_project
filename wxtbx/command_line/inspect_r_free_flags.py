# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import wxtbx.plots
from libtbx.utils import Sorry, Abort, Usage
from libtbx import group_args
import wx
import os
import sys

def run (args) :
  app = wx.App(0)
  if (len(args) == 0) :
    ask_for_file_name()
  else :
    file_name = args[0]
    if (not os.path.isfile(file_name)) :
      raise Usage("wxtbx.inspect_r_free_flags data.mtz")
    display_file_info(file_name)
  app.MainLoop()

def ask_for_file_name (parent=None) :
  import iotbx.file_reader
  file_name = wx.FileSelector(
    message="Select a reflections file",
    flags=wx.OPEN,
    wildcard=iotbx.file_reader.get_wildcard_strings(["hkl"]))
  display_file_info(file_name, parent)

def display_file_info (file_name, parent=None) :
  from iotbx import reflection_file_editor
  from iotbx import file_reader
  from iotbx import data_plots
  if (file_name == "") or (file_name is None) :
    raise Abort()
  hkl_in = file_reader.any_file(file_name, force_type="hkl")
  hkl_in.check_file_type("hkl")
  labels = []
  tables = []
  stats = []
  arrays_and_flags = hkl_in.file_server.get_r_free_flags(
    file_name=hkl_in.file_name,
    label=None,
    test_flag_value=None,
    disable_suitability_test=False,
    parameter_scope=None,
    return_all_valid_arrays=True)
  for array, test_flag_value in arrays_and_flags :
    labels.append(array.info().label_string())
    (n_bins, n_free, sse, accu) = reflection_file_editor.get_r_free_stats(
      miller_array=array,
      test_flag_value=test_flag_value)
    table = data_plots.table_data(
      title="Distribution of R-free test set",
      column_labels=["Reflection #", "Fraction counted"],
      graph_names=["Fraction free vs. resolution"],
      graph_labels=[["Reflection #", "ntest/total"]],
      graph_columns=[[0,1]],
      data=[accu.reflection_counts, accu.free_fractions])
    tables.append(table)
    stats.append(group_args(
      percent_free=n_free/array.indices().size(),
      completeness=array.completeness()))
  if (len(labels) == 0) :
    raise Sorry("No recognizable R-free flags found in this file.")
  frame = wx.Frame(parent, -1, "R-free flags info",
    style=wx.DEFAULT_FRAME_STYLE)
  szr = wx.BoxSizer(wx.VERTICAL)
  frame.SetSizer(szr)
  panel = RfreeInspector(frame, -1)
  panel.set_data(
    file_name=file_name,
    labels=labels,
    tables=tables,
    stats=stats,
    flags_and_values=[ (a.data(), f) for (a, f) in arrays_and_flags ])
  szr.Add(panel, 1, wx.ALL|wx.EXPAND)
  szr.Fit(panel)
  frame.Fit()
  frame.Show()

class RfreeInspector (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    grid = wx.FlexGridSizer(cols=2, rows=4)
    self.sizer.Add(grid, 0, wx.ALL, 5)
    grid.Add(wx.StaticText(self, -1, "File name:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.TextCtrl(self, -1, size=(400,-1),
      style=wx.TE_READONLY)
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid.Add(wx.StaticText(self, -1, "Column(s):"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.labels_choice = wx.Choice(self, -1, size=(200,-1))
    grid.Add(self.labels_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid.Add(wx.StaticText(self, -1, "Completeness:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.comp_txt = wx.StaticText(self, -1, "", size=(100,-1))
    grid.Add(self.comp_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid.Add(wx.StaticText(self, -1, "Overall % free:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.percent_txt = wx.StaticText(self, -1, "", size=(100,-1))
    grid.Add(self.percent_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt = wx.StaticText(self, -1,
"""The graph on the right shows the total nubmer of test set \
reflections as a function of resolution.  A set picked at random will have an \
approximately straight line; if the test set was picked in thin resolution \
shells to avoid bias caused by non-crystallographic symmetry, it will \
increase in sharp steps.  CCTBX/PHENIX can not currently extend a set picked \
in shells; we suggest creating a new set in such cases.""")
    txt.Wrap(780)
    self.sizer.Add(txt, 0, wx.ALL, 5)
    plot_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.sizer.Add(plot_sizer, 1, wx.ALL|wx.EXPAND)
    self.plot1 = RfreeFlagsPlot(self,
      figure_size=(5,5))
    plot_sizer.Add(self.plot1, 1, wx.ALL|wx.EXPAND)
    self.plot2 = wxtbx.plots.iotbx_data_plot_base(
      parent=self,
      tables=[],
      size=(400,400),
      show_data_points=False)
    plot_sizer.Add(self.plot2, 1, wx.ALL|wx.EXPAND)
    self.SetSizer(self.sizer)
    self.Bind(wx.EVT_CHOICE, self.OnChooseArray, self.labels_choice)
    self.stats = []

  def set_data (self, file_name, labels, tables, stats, flags_and_values) :
    self.stats = stats
    self.flags_and_values = flags_and_values
    self.file_txt.SetValue(file_name)
    self.labels_choice.SetItems(labels)
    self.labels_choice.SetSelection(0)
    self.plot2.set_tables(tables)
    self.OnChooseArray(None)

  def OnChooseArray (self, evt) :
    array_selection = self.labels_choice.GetSelection()
    flags, value = self.flags_and_values[array_selection]
    self.plot1.show_pie(flags, value)
    self.plot2.set_plot(graph_name="Fraction free vs. resolution",
      table_index=array_selection)
    self.comp_txt.SetLabel("%.2f%%" %
      (100 * self.stats[array_selection].completeness))
    self.percent_txt.SetLabel("%.2f%%" %
      (100 * self.stats[array_selection].percent_free))
    self.Refresh()

class RfreeFlagsPlot (wxtbx.plots.plot_container) :
  def show_pie (self, flags, test_flag_value) :
    values = sorted(list(set(flags)))
    n_values = len(values)
    labels = []
    colors = []
    counts = []
    for i, value in enumerate(values) :
      if (value == test_flag_value) :
        colors.append((1.,0.,0.))
      else :
        colors.append((0.5*(n_values-i)/n_values, 0.5*(n_values-i)/n_values,
          1.))
      labels.append(str(value))
      counts.append(flags.count(value))
    self.figure.clear()
    p = self.figure.add_subplot(111)
    p.set_position([0.1,0.1,0.8, 0.8])
    p.pie(counts, labels=labels, colors=colors)
    p.set_title("Flag values")
    self.canvas.draw()

if (__name__ == "__main__") :
 run(sys.argv[1:])
