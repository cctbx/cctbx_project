
from libtbx import str_utils
from libtbx import easy_pickle
import wx
import os
import sys

columns = ["Solution", "Metric fit", "RMSD", "#spots", "Crystal system",
    "Unit cell", "Volume"]
column_sizes = [60, 80, 80, 60, 120, 250, 80]
column_alignments = [0, 1, 1, 1, 1, 0, 1]

class LatticeData (object) :
  def __init__ (self, labelit_possible=()) :
    self.labelit_possible = labelit_possible

  def GetItemCount (self) :
    return len(self.labelit_possible)

  def GetItemImage (self, item) :
    return 0

  def GetItemText (self, item, col) :
    n_items = self.GetItemCount()
    solution = self.labelit_possible[item]
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

class LatticeListCtrl (wx.ListCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['style'] = wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.LC_VIRTUAL
    wx.ListCtrl.__init__(self, *args, **kwds)
    for i, label in enumerate(columns) :
      self.InsertColumn(i, label, column_alignments[i])
      self.SetColumnWidth(i, column_sizes[i])
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick, self)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick, self)
    self.dataSource = LatticeData()
    self.RefreshAllItems()

  def RefreshAllItems (self) :
    n_items = self.dataSource.GetItemCount()
    self.SetItemCount(n_items)
    if (n_items > 0) :
      self.RefreshItems(0, n_items - 1)

  def OnGetItemImage (self, item) :
    return self.dataSource.GetItemImage(item)

  def OnGetItemAttr (self, item) :
    pass

  def OnGetItemText (self, item, col) :
    return self.dataSource.GetItemText(item, col)

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
