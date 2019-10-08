
"""
Various tools and controls for displaying tabular data.
"""

from __future__ import absolute_import, division, print_function
import wxtbx.plots
import wx
import sys

class TableView(wx.ListCtrl):
  """
  Table display using wx.ListCtrl in combination with iotbx.data_plots.
  """
  def __init__(self, *args, **kwds):
    style = kwds.get("style", 0)
    if (not style & wx.LC_REPORT):
      style = style | wx.LC_REPORT
      kwds['style'] = style
    super(TableView, self).__init__(*args, **kwds)
    self._table = None

  def SetTable(self, table):
    self._table = table
    w, h = self.GetSize()
    rows = self._table.export_rows()
    labels = rows[0]
    col_width = ((w - 20) / len(labels))
    for i_lab, label in enumerate(labels):
      self.InsertColumn(i_lab, label)
      self.SetColumnWidth(i_lab, col_width)
    for row in rows[1:] :
      assert len(row) == len(labels), labels
      idx = self.InsertStringItem(sys.maxunicode, row[0])
      for i, cell in enumerate(row[1:]):
        self.SetStringItem(idx, i+1, cell)

  def OnViewGraphs(self, evt):
    """
    Open a wxtbx.plots.loggraph_frame window with the table.
    """
    assert hasattr(self._table, "get_graph")
    graph_frame = wxtbx.plots.loggraph(
      parent=self.GetParent(),
      title=self._table.title,
      tables=[self._table])
    graph_frame.Show()

# TESTING
if (__name__ == "__main__"):
  from iotbx import data_plots
  table = data_plots.table_data(
    title = "Resolution shell statistics",
    column_labels = ["1/resol^2", "Nrefl", "R-free", "FOM"],
    graph_names = ["R-free vs. resolution", "FOM vs. resolution"],
    graph_columns = [[0,2], [0,3]],
    data = [[0.02, 0.04, 0.06, 0.08, 0.10],
            [2004, 2084, 2037, 1949, 1783],
            [0.25, 0.23, 0.27, 0.28, 0.38],
            [0.89, 0.88, 0.83, 0.75, None]])
  app = wx.App(0)
  frame = wx.Frame(parent=None, title="Test frame", size=(600,280))
  panel = wx.Panel(parent=frame)
  ctrl = TableView(panel, pos=(10,10), size=(580,180))
  ctrl.SetTable(table)
  btn = wx.Button(parent=panel, label="Show graphs", pos=(10,210))
  frame.Bind(wx.EVT_BUTTON, ctrl.OnViewGraphs, btn)
  frame.Show()
  app.MainLoop()
