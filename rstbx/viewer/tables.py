
import wx

class BinData (object) :
  def __init__ (self, table) :
    assert isinstance(table, list)
    self.table = table

  def GetItemCount (self) :
    return len(self.table)

  def GetItemImage (self, item) :
    return 0

  def GetItemText (self, item, col) :
    n_items = self.GetItemCount()
    assert (item < n_items) and (0 <= col <= 4)
    bin = self.table[item]
    if (col == 0) :
      return "%d" % bin.i_bin
    elif (col == 1) :
      return "%g - %g" % bin.d_max_min
    elif (col == 2) :
      return "%d / %d" % bin.completeness
    elif (col == 3) :
      return "%8.1f" % bin.mean_I
    else :
      return "%8.1f" % bin.mean_I_sigI

class ResolutionBinTable (wx.ListCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['style'] = wx.LC_REPORT|wx.LC_VIRTUAL|wx.LC_SINGLE_SEL
    kwds['size'] = (600,200)
    wx.ListCtrl.__init__(self, *args, **kwds)
    labels = ["Bin","Resolution","Completeness", "<I>", "<I/sig(I)>"]
    widths = [60,180,140,100,100]
    for i, label in enumerate(labels) :
      self.InsertColumn(i, label, wx.LIST_FORMAT_RIGHT)
      self.SetColumnWidth(i, widths[i])
    self.dataSource = BinData([])
    self.RefreshAllItems()

  def SetBins (self, table) :
    self.dataSource = BinData(table)
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

def exercise () :
  from libtbx import easy_pickle
  import os
  import sys
  indexing_file = sys.argv[1]
  p = easy_pickle.load(indexing_file)
  table = p['table_raw']
  assert os.path.isfile(indexing_file)
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Resolution bins", size=(600,300))
  panel = wx.Panel(frame)
  blist = ResolutionBinTable(panel, -1)
  blist.SetBins(table)
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__") :
  exercise()
