import wx
import wx.lib.stattext
import unicodedata


class Inspector(wx.Panel):

  def __init__(self, parent, id=wx.ID_ANY, label="", colour=wx.BLACK,
               pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
    assert parent.GetSizer() is not None
    super(Inspector, self).__init__(parent, id)
    self._normal_colour, self._pressed_colour = [
      wx.Colour(*c) for c in [(208,)*3, (180,)*3] ]
    if wx.Platform == "__WXGTK__":
      self._header = wx.lib.stattext.GenStaticText(self, -1, "",
                                                   style=wx.BORDER_NONE)
    elif wx.Platform == "__WXMAC__":
      self._header = wx.StaticText(self)
    self._header.SetBackgroundColour(self._normal_colour)
    self._header.Bind(wx.EVT_LEFT_DOWN, self._on_left_down)
    self._header.Bind(wx.EVT_LEFT_UP, self._on_left_up)
    self._header.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave_window)
    self._header.Bind(wx.EVT_ENTER_WINDOW, self._on_enter_window)
    self._pane = wx.Panel(self)
    self._right_triangle = unicodedata.lookup(
      'BLACK RIGHT-POINTING TRIANGLE')
    self._down_triangle = unicodedata.lookup(
      'BLACK DOWN-POINTING TRIANGLE')
    self._pressing = False
    self._expanded = True
    self._label = None
    self.on_click_in_header = None
    s = wx.BoxSizer(wx.VERTICAL)
    s.Add(self._header, 0, wx.EXPAND, 5)
    s.Add(self._pane, flag=wx.LEFT|wx.RIGHT|wx.TOP, border=5)
    self.SetSizer(s)
    self.SetLabel(label)

  def _on_left_down(self, event):
    self._header.SetBackgroundColour(self._pressed_colour)
    self._pressing = True
    event.Skip()

  def _on_left_up(self, event):
    self._header.SetBackgroundColour(self._normal_colour)
    self._pressing = False
    self._on_click_in_header()
    event.Skip()

  def _on_enter_window(self, event):
    if self._pressing:
      if event.LeftIsDown():
        self._header.SetBackgroundColour(self._pressed_colour)
      else:
        self._pressing = False
    event.Skip()

  def _on_leave_window(self, event):
    if self._pressing:
      self._header.SetBackgroundColour(self._normal_colour)
    event.Skip()

  def _on_click_in_header(self):
    self.Collapse(self._expanded)
    if self.on_click_in_header is not None:
      self.on_click_in_header()

  def IsExpanded(self):
    return self._expanded

  def Collapse(self, f):
    self._expanded = not f
    self._update()

  def SetLabel(self, lbl):
    if self._label == lbl: return
    self._label = lbl
    self._update()

  def _update(self):
    if self._expanded:
      arrow = self._down_triangle
      self._pane.Show()
    else:
      arrow = self._right_triangle
      self._pane.Hide()
    lbl = "%s\t%s" % (arrow, self._label)
    self._header.SetLabel(lbl)
    top = wx.GetTopLevelParent(self)
    top_sizer = top.GetSizer()
    top_sizer.SetSizeHints(top)
    if self._expanded:
      top.Fit()
    else:
      sz = top_sizer.CalcMin()
      top.SetClientSize(sz)

  def GetLabel(self):
    return self._label

  def GetPane(self):
    return self._pane

  def width_for_stacking(self):
    w_h = self._header.GetBestSize()[0]
    w_p = self._pane.GetBestSize()[0]
    w = max(w_h, w_p)
    return w


class InspectorToolFrame(wx.MiniFrame):

  def __init__(self, parent, id=wx.ID_ANY,
               pos=wx.DefaultPosition, size=wx.DefaultSize):
    super(InspectorToolFrame, self).__init__(parent, id, title=" ",
                                             pos=pos, size=size,
                                             style=wx.CAPTION|wx.CLOSE_BOX)
    s = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(s)
    if parent is not None:
      self.Bind(wx.EVT_CLOSE, lambda event: self.Hide())
      self.Bind(wx.EVT_MOVE, self._on_move)
      parent.Bind(wx.EVT_MOVE, self._on_move_or_resize_parent)
      parent.Bind(wx.EVT_SIZE, self._on_move_or_resize_parent)
      self._position_relative_to_parent = None
      self._parent_moved = False

  def _on_move(self, event):
    if self._parent_moved:
      self._parent_moved = False
      return
    x, y = event.GetPosition()
    x_p, y_p = self.Parent.GetPosition()
    self._position_relative_to_parent = (x - x_p, y - y_p)

  def _on_move_or_resize_parent(self, event):
    self.place()
    event.Skip()

  def place(self):
    gap = 5
    x, y, w, h = self.Parent.Rect
    w_t, h_t = self.Size
    w_s, h_s = wx.GetDisplaySize()
    w_s -= 50
    if self._position_relative_to_parent is None:
      x_t_rel, y_t_rel = w + gap, 0
    else:
      x_t_rel, y_t_rel = self._position_relative_to_parent
    y_t = y + y_t_rel
    x_t = min(x + x_t_rel, w_s - w_t)
    self._parent_moved = True
    self.MoveXY(x_t, y_t)

  def Layout(self):
    s = self.Sizer
    for i in self.Children:
      s.Add(i, 0, wx.EXPAND|wx.BOTTOM, 5)
    w = max([ i.width_for_stacking() for i in self.Children ])
    s.SetMinSize((w,-1))
    s.SetSizeHints(self)

if __name__ == '__main__':
  a = wx.App()
  w = InspectorToolFrame(None, pos=(500, 30))

  i1 = Inspector(w, label="Header")
  def f():
    print "Click in 'Header'"
  i1.on_click_in_header = f
  pane = i1.GetPane()
  s = wx.BoxSizer(wx.HORIZONTAL)
  s.Add(wx.Button(pane, label="Click me!"))
  s.Add(wx.CheckBox(pane, label="Check me!"))
  pane.SetSizer(s)

  i2 = Inspector(w, label="Longer header")
  def g():
    print "Click in 'Longer header'"
  i2.on_click_in_header = g
  pane = i2.GetPane()
  s = wx.BoxSizer(wx.HORIZONTAL)
  s.Add(wx.StaticText(pane, label="Static"))
  s.Add(wx.StaticText(pane, label="..."))
  pane.SetSizer(s)

  w.Layout()

  a.SetTopWindow(w)
  w.Show()
  a.MainLoop()
