
from rstbx.viewer import results_base
import wx

class ListBase (wx.ListCtrl) :
  def Reset (self) :
    self.dataSource = results_base.EmptyData()
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
