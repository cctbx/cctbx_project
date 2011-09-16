
from rstbx.viewer import controls, results_base
from libtbx import str_utils
from libtbx import easy_pickle
import wx
import os
import sys

columns = ["Solution", "Metric fit", "RMSD", "#spots", "Crystal system",
    "Unit cell", "Volume"]
column_sizes = [60, 80, 80, 60, 120, 250, 80]
column_alignments = [0, 1, 1, 1, 1, 0, 1]

class LatticeData (results_base.TableData) :
  def GetItemText (self, item, col) :
    n_items = self.GetItemCount()
    solution = self.table[item]
    if (col == 0) :
      return solution['counter']
    elif (col == 1) :
      return "%s%s" % (str_utils.format_value("%.4f",
        solution['max_angular_difference']), chr(161))
    elif (col == 2) :
      return str_utils.format_value("%.4f", solution['residual'])
    elif (col == 3) :
      return str(solution['count_GOOD'])
    elif (col == 4) :
      return "%s %s" % (solution['system'], solution['bravais'])
    elif (col == 5) :
      unit_cell = solution['orient'].unit_cell()
      return "%g %g %g %g %g %g" % unit_cell.parameters()
    elif (col == 6) :
      unit_cell = solution['orient'].unit_cell()
      return str(int(unit_cell.volume()))

class LatticeListCtrl (controls.ListBase) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['style'] = wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.LC_VIRTUAL
    controls.ListBase.__init__(self, *args, **kwds)
    for i, label in enumerate(columns) :
      self.InsertColumn(i, label, column_alignments[i])
      self.SetColumnWidth(i, column_sizes[i])
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick, self)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick, self)
    self.dataSource = results_base.EmptyData()
    self.RefreshAllItems()

  def SetIndexingResults (self, labelit_possible) :
    assert isinstance(labelit_possible, list)
    self.dataSource = LatticeData(labelit_possible)
    self.RefreshAllItems()

  def OnSelect (self, event) :
    pass

  def OnDeSelect (self, event) :
    pass

  def OnDoubleClick (self, event) :
    pass

  def OnRightClick (self, event) :
    pass

class IndexingPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.lattice_list = LatticeListCtrl(self, -1, size=(760, 200))
    txt1 = wx.StaticText(self, -1, "Indexing results:")
    f1 = txt1.GetFont()
    f1.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(f1)
    szr.Add(txt1, 0, wx.ALL|wx.EXPAND, 5)
    szr.Add(self.lattice_list, 0, wx.ALL|wx.EXPAND, 10)

  def __getattr__ (self, name) :
    return getattr(self.lattice_list, name)

if (__name__ == "__main__") :
  labelit_file = sys.argv[1]
  assert (os.path.isfile(labelit_file))
  labelit_result = easy_pickle.load(labelit_file)
  assert isinstance(labelit_result, list)
  app = wx.App(0)
  app.locale = wx.Locale(wx.LANGUAGE_ENGLISH)
  frame = wx.Frame(None, -1, "Lattice solutions", size=(800,320))
  panel = wx.Panel(frame)
  llist = LatticeListCtrl(panel, -1, size=(760, 200))
  llist.SetIndexingResults(labelit_result)
  frame.Show()
  app.MainLoop()
