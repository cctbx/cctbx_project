
from wxtbx import metallicbutton
import wxtbx.bitmaps
import wx
import sys

class ListEditor (wx.Panel) :
  header_label = "Items"
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self._default_label = "---"
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.list = wx.ListCtrl(
      parent=self,
      id=-1,
      style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    self.list.InsertColumn(0, self.header_label, width=460)
    self.list.SetMinSize((480,160))
    self.list.SetItemSpacing(5)
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self.list)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self.list)
    szr.Add(self.list, 1, wx.EXPAND|wx.ALL, 5)
    btn_szr = self.CreateButtons()
    szr.Add(btn_szr, 0, wx.LEFT|wx.BOTTOM|wx.RIGHT, 5)
    edit_szr = wx.BoxSizer(wx.HORIZONTAL)
    edit_label = wx.StaticText(
      parent=self,
      label="Edit selected:")
    edit_szr.Add(edit_label, 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 5)
    self.edit = wx.TextCtrl(
      parent=self,
      size=(300,-1),
      style=wx.TE_PROCESS_ENTER)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, self.edit)
    edit_szr.Add(self.edit, 1, wx.EXPAND|wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Add(edit_szr, 0, wx.EXPAND|wx.TOP|wx.LEFT|wx.BOTTOM, 5)
    szr.Layout()

  def CreateButtons (self) :
    btn_szr = wx.BoxSizer(wx.HORIZONTAL)
    add_btn = metallicbutton.MetallicButton(
      parent=self,
      label="Add",
      bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "edit_add", 16),
      highlight_color=(200,220,240))
    self.Bind(wx.EVT_BUTTON, self.OnAdd, add_btn)
    del_btn = metallicbutton.MetallicButton(
      parent=self,
      label="Delete",
      bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "editdelete", 16),
      highlight_color=(200,220,240))
    self.Bind(wx.EVT_BUTTON, self.OnDelete, del_btn)
    update_btn = metallicbutton.MetallicButton(
      parent=self,
      label="Update item",
      bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "recur", 16),
      highlight_color=(200,220,240))
    self.Bind(wx.EVT_BUTTON, self.OnUpdate, update_btn)
    btn_szr.Add(add_btn, 0, wx.RIGHT, 5)
    btn_szr.Add(del_btn, 0, wx.RIGHT, 5)
    btn_szr.Add(update_btn, 0, wx.RIGHT, 5)
    return btn_szr

  def AddItem (self, item) :
    return self.list.InsertStringItem(sys.maxint, item)

  def OnAdd (self, event) :
    i = self.AddItem(self._default_label)
    self.list.Select(i, 1)
    self.edit.SetFocus()

  def OnDelete (self, event) :
    i = self.list.GetFirstSelected()
    self.list.DeleteItem(i)
    self.edit.SetValue("")

  def OnUpdate (self, event) :
    evt_type = event.GetEventType()
    txt = self.edit.GetValue()
    if (txt == "") or (txt is None) :
      txt = self._default_label
    i = self.list.GetFirstSelected()
    if (i == -1) :
      if (event.GetEventType() == wx.EVT_TEXT_ENTER.typeId) :
        i = self.AddItem(txt)
        self.list.Select(i, 1)
    else :
      self.list.SetItemText(i, txt)
    self.list.SetFocus()

  def OnSelect (self, event) :
    item = self.list.GetFirstSelected()
    txt = self.list.GetItemText(item)
    if (txt == self._default_label) :
      txt = ""
    self.edit.SetValue(txt)

  def OnDeSelect (self, event) :
    self.edit.SetValue("")

  def SetDefaultItemLabel (self, label) :
    self._default_label = label

  def GetValues (self) :
    items = []
    i = 0
    while (i < self.list.GetItemCount()) :
      txt = self.list.GetItemText(i)
      if (txt == self._default_label) :
        txt = None
      items.append(txt)
    return items

if __name__ == "__main__" :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Test frame")
  szr = wx.BoxSizer(wx.VERTICAL)
  frame.SetSizer(szr)
  panel = ListEditor(parent=frame)
  szr.Add(panel, 1, wx.EXPAND)
  szr.Layout()
  szr.Fit(panel)
  frame.Fit()
  frame.Show()
  app.MainLoop()
