
from wxtbx import metallicbutton
import wxtbx.bitmaps
import wx
import sys

class ListEditor (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self._default_label = "---"
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.sizer = szr
    self.CreateList()
    self.buttons = wx.BoxSizer(wx.HORIZONTAL)
    add_btn = self.AddControlButton(
      label="Add",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "edit_add", 16))
    self.Bind(wx.EVT_BUTTON, self.OnAdd, add_btn)
    del_btn = self.AddControlButton(
      label="Delete",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16))
    self.Bind(wx.EVT_BUTTON, self.OnDelete, del_btn)
    clear_btn = self.AddControlButton(
      label="Clear all",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "editdelete", 16))
    self.Bind(wx.EVT_BUTTON, self.OnDeleteAll, clear_btn)
    update_btn = self.AddControlButton(
      label="Update item",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "recur", 16))
    self.Bind(wx.EVT_BUTTON, self.OnUpdate, update_btn)
    szr.Add(self.buttons, 0, wx.LEFT|wx.BOTTOM|wx.RIGHT, 5)
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
    self.sizer = szr
    self._label = None
    self._callback = None

  def CreateList (self) :
    self.list = wx.ListCtrl(
      parent=self,
      id=-1,
      style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    self.list.InsertColumn(0, "Items", width=460)
    self.list.SetMinSize((480,160))
    if (hasattr(self.list, "SetItemSpacing")) :
      self.list.SetItemSpacing(5)
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self.list)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self.list)
    self.sizer.Add(self.list, 1, wx.EXPAND|wx.ALL, 5)

  def SetLabel (self, label, font_weight=wx.FONTWEIGHT_BOLD) :
    if (self._label is not None) :
      self._label.SetLabel(label)
    else :
      self._label = wx.StaticText(parent=self, label=label)
      font = self._label.GetFont()
      font.SetWeight(font_weight)
      self._label.SetFont(font)
      self.sizer.Insert(0, self._label, 0, wx.TOP|wx.LEFT, 5)
    self.sizer.Layout()

  def SetColumnHeader (self, header) :
    col = self.list.GetColumn(0)
    col.SetText(header)
    self.list.SetColumn(0, col)

  def SetToolTip (self, tool_tip) :
    if isinstance(tool_tip, str) :
      self.list.SetToolTip(wx.ToolTip(tool_tip))
    else :
      self.list.SetToolTip(tool_tip)

  def AddControlButton (self, label, bitmap) :
    btn = metallicbutton.MetallicButton(
      parent=self,
      label=label,
      bmp=bitmap,
      highlight_color=(200,220,240))
    self.buttons.Add(btn, 0, wx.RIGHT, 5)
    return btn

  def AddItem (self, item) :
    return self.list.InsertStringItem(sys.maxint, item)

  def OnAdd (self, event) :
    i = self.AddItem(self._default_label)
    self.list.Select(i, 1)
    self.edit.SetFocus()
    self.call_back()

  def OnDelete (self, event) :
    i = self.list.GetFirstSelected()
    if (i >= 0) :
      self.list.DeleteItem(i)
    self.edit.SetValue("")
    self.call_back()

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
    self.call_back()

  def OnSelect (self, event) :
    item = self.list.GetFirstSelected()
    txt = str(self.list.GetItemText(item))
    if (txt == self._default_label) :
      txt = ""
    self.edit.SetValue(txt)
    self.call_back()

  def OnDeSelect (self, event) :
    self.edit.SetValue("")
    self.call_back()

  def SetDefaultItemLabel (self, label) :
    self._default_label = label

  def GetValues (self) :
    items = []
    i = 0
    n = self.list.GetItemCount()
    while (i < n) :
      txt = str(self.list.GetItemText(i))
      if (txt == self._default_label) :
        txt = None
      items.append(txt)
      i += 1
    return items

  def DeleteAllItems (self) :
    self.list.DeleteAllItems()
    self.call_back()

  def OnDeleteAll (self, evt) :
    if (self.list.GetItemCount() == 0) :
      return False
    confirm = wx.MessageBox(caption="Confirm delete",
      message="Are you sure you want to delete all items in the list?")
    if (confirm == wx.OK) :
      self.DeleteAllItems()

  def SetSelectedValue (self, txt) :
    i = self.list.GetFirstSelected()
    if (i == -1) :
      return
    self.list.SetItemText(i, txt)
    self.edit.SetValue(txt)

  def GetSelectedValue (self) :
    i = self.list.GetFirstSelected()
    if (i == -1) :
      return None
    return str(self.list.GetItemText(i))

  def SetCallback (self, callback) :
    assert hasattr(callback, "__call__")
    self._callback = callback

  def call_back (self) :
    if (self._callback is not None) :
      self._callback()

if __name__ == "__main__" :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Test frame")
  szr = wx.BoxSizer(wx.VERTICAL)
  frame.SetSizer(szr)
  panel = ListEditor(parent=frame)
  panel.SetLabel("TLS groups:")
  panel.SetColumnHeader("Atom selection")
  szr.Add(panel, 1, wx.EXPAND)
  szr.Layout()
  szr.Fit(panel)
  frame.Fit()
  frame.Show()
  app.MainLoop()
