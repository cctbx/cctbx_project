
from wx.lib.embeddedimage import PyEmbeddedImage
import wx
import sys

button_ok = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAACK0lEQVQ4jaWSXUhTYRjH/+/O"
    "3uPOOfs4W2oRRGI0jPYRStaIjMS6SyQoC+kD0rqJwAjqqo+bbpLoQorIhjaShIqRLUtzdaHD"
    "DywpmzZ0NIS0rDV1tg87e7uQRMps1R+em+fh91z8+RHGGP4nqr8FhFqi/ucHopNI8qip23RQ"
    "bplfMsbSGuEGpOVn5P72wVZ2wLmXqXbhIWMMJJ0OhGtEMoRk363DjbaGrjqMJcchJw1wP3jc"
    "PP8gq5FfBkZyJyoSvQthfT2RMob0PldVo83lv4nWTx5wMYpZLw0mZqIFKgDQ1JJM80d7YEPY"
    "0SHVk+0/4JX3Ja32ndHnOn7b1jTegF7aAS3VIeyJDUenI5Zoy2xEnX1RzC4WS/zVpdWmqBLF"
    "pCfiWXFXs1MStT10SNdZc+SyrfnzPQzS1xATOvibQiMklbImnrA4AHCKTK4XlTgKPwjviSwa"
    "UZCXT/u8A+UI0apzlWfNL+M9CHBvoEwR9DsDQS5LscbvzMEA5koUyjRP83av2movtPP5uo2g"
    "hIeBGjGovMII3mL6ywye1XQH+U1JS/ioElvYkRoAYu54CVemaoOMIrb5G2/m1yHMTWCCG4My"
    "BXRcfRHMqchcP1A8GsdPmbdKcad20D2qdllv2mLcZswQiIDUJPD8UlfQccJqebS28xd4UZEy"
    "9qnbyntLExfCp1nuqZzhk18rNUsJtuhSd0j0rjm/uq+OXRH+ZOhvD8fi+2k6iqel8lL5Dvo2"
    "ZEkd1D8WAAAAAElFTkSuQmCC")

blank = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAQAAAC1+jfqAAAAAXNSR0IArs4c6QAAAAJiS0dE"
    "AP+Hj8y/AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH2QsNFBwTBbkZpAAAABl0RVh0"
    "Q29tbWVudABDcmVhdGVkIHdpdGggR0lNUFeBDhcAAAAOSURBVCjPY2AYBaMAAQACEAABFMLA"
    "kgAAAABJRU5ErkJggg==")

class CheckListCtrl (wx.ListCtrl) :
  def __init__ (self, *args, **kwds) :
    wx.ListCtrl.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Bind(wx.EVT_CHAR, self.OnChar)
    self._checklist = []
    il = wx.ImageList(16, 16, True)
    il.Add(blank.GetBitmap())
    il.Add(button_ok.GetBitmap())
    self.AssignImageList(il, wx.IMAGE_LIST_SMALL)
    self._multiple_sele = not (self.GetWindowStyle() & wx.LC_SINGLE_SEL)
    tip = """\
You may check or uncheck items by double-clicking them.  You may also use \
the keyboard: arrow keys move up and down the list, Return checks items, \
backspace/delete unchecks them."""
    if self._multiple_sele :
      tip += """  You can change the status of multiple items at once by \
holding down the shift key in combination with mouse selection."""
    self.SetToolTip(wx.ToolTip(tip))

  def IsItemChecked (self, item) :
    return self._checklist[item]

  def GetCheckedItems (self) :
    i = 0
    items = []
    while (i < self.GetItemCount()) :
      if self._checklist[i] :
        items.append(i)
      i += 1
    return items

  def GetCheckedItemsText (self) :
    items = self.GetCheckedItems()
    strings = []
    for item in items :
      strings.append(self.GetItemText(item))
    return strings

  def InsertStringItem (self, *args, **kwds) :
    i = wx.ListCtrl.InsertStringItem(self, *args, **kwds)
    self._checklist.append(False)
    self.SetItemImage(i, 0)
    return i

  def InsertItem (self, *args, **kwds) :
    i = wx.ListCtrl.InsertItem(self, *args, **kwds)
    self._checklist.append(False)
    self.SetItemImage(i, 0)
    return i

  def CheckItem (self, item, check=True) :
    self._checklist[item] = check
    if check :
      self.SetItemImage(item, 1)
    else :
      self.SetItemImage(item, 0)

  def ToggleItems (self, items, check=None) :
    if (len(items) == 0) : return
    assert (check in [True, False, None])
    item_status = []
    for item in items :
      item_status.append(self.IsItemChecked(item))
    if (item_status.count(True) > 0) :
      for item in items :
        if (check is None) :
          self.CheckItem(item, False)
        else :
          self.CheckItem(item, check)
    else :
      for item in items :
        if (check is None) :
          self.CheckItem(item, True)
        else :
          self.CheckItem(item, check)

  def OnDoubleClick (self, event) :
    i, flags = self.HitTest(event.GetPosition())
    if (flags & wx.LIST_HITTEST_ONITEM) :
      if self._multiple_sele :
        item_list = self.GetSelectedItems()
        self.ToggleItems(item_list)
      else :
        is_checked = self.IsItemChecked(i)
        self.CheckItem(i, not is_checked)

  def GetSelectedItems (self) :
    item_list = []
    item = self.GetFirstSelected()
    while (item >= 0) :
      item_list.append(item)
      item = self.GetNextSelected(item)
    return item_list

  def OnChar (self, event) :
    key = event.GetKeyCode()
    if (key == wx.WXK_RETURN) :
      items = self.GetSelectedItems()
      self.ToggleItems(items, check=True)
    elif (key == wx.WXK_DELETE) or (key == wx.WXK_BACK) :
      items = self.GetSelectedItems()
      self.ToggleItems(items, check=False)
    elif (key == wx.WXK_UP) :
      item = self.GetFirstSelected()
      if (item > 0) :
        self.Select(item - 1)
        self.EnsureVisible(item - 1)
    elif (key == wx.WXK_DOWN) :
      items = self.GetSelectedItems()
      if (len(items) == 0) : return
      if (items[-1] < (self.GetItemCount() - 1)) :
        self.Select(items[-1] + 1)
        self.EnsureVisible(items[-1] + 1)

  def UpdateParentControls (self) :
    pass

  def OnSelect (self, event) :
    pass

  def OnDeSelect (self, event) :
    pass

def demo () :
  import random
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Test frame for CheckListCtrl")
  panel = wx.Panel(frame, -1)
  sizer = wx.BoxSizer(wx.VERTICAL)
  frame.SetSizer(sizer)
  sizer.Add(panel, 1, wx.EXPAND)
  sizer2 = wx.BoxSizer(wx.VERTICAL)
  panel.SetSizer(sizer2)
  txt = wx.StaticText(panel, -1, "Double-click an item to check/uncheck it.")
  sizer2.Add(txt, 0, wx.ALL, 5)
  ctrl1 = CheckListCtrl(panel, -1, size=(600,200), style=wx.LC_REPORT)
  ctrl1.InsertColumn(0, "Item name")
  ctrl1.SetColumnWidth(0, 580)
  sizer2.Add(ctrl1, 0, wx.ALL, 5)
  btn1 = wx.Button(panel, -1, "Print checked items")
  sizer2.Add(btn1, 0, wx.ALL, 5)
  def OnShow1 (event) :
    items = ctrl1.GetCheckedItemsText()
    print "Control 1 selection:"
    print "\n".join([ " " + item for item in items ])
  frame.Bind(wx.EVT_BUTTON, OnShow1, btn1)
  txt2 = wx.StaticText(panel, -1, "This control will only allow you to "+
    "check/uncheck a single item at a time.")
  sizer2.Add(txt2, 0, wx.ALL, 5)
  ctrl2 = CheckListCtrl(panel, -1, size=(600,200),
    style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
  ctrl2.InsertColumn(0, "Item name")
  ctrl2.InsertColumn(1, "Other data")
  ctrl2.SetColumnWidth(0, 380)
  ctrl2.SetColumnWidth(1, 200)
  sizer2.Add(ctrl2, 0, wx.ALL, 5)
  btn2 = wx.Button(panel, -1, "Print checked items")
  def OnShow2 (event) :
    items = ctrl2.GetCheckedItemsText()
    print "Control 2 selection:"
    print "\n".join([ " " + item for item in items ])
  sizer2.Add(btn2, 0, wx.ALL, 5)
  frame.Bind(wx.EVT_BUTTON, OnShow2, btn2)
  for i in range(10) :
    ctrl1.InsertStringItem(sys.maxsize, "Item %d" % i)
    item = ctrl2.InsertStringItem(sys.maxsize, "Item %d" % i)
    ctrl2.SetStringItem(item, 1, str(random.random() * 1000))
  sizer2.Layout()
  sizer.Layout()
  sizer.Fit(panel)
  frame.Fit()
  frame.Show()
  app.MainLoop()

if __name__ == "__main__" :
  demo()
