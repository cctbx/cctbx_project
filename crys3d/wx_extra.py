import wx
import wx.lib.stattext

from libtbx.utils import group_args

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
    if pos is wx.DefaultPosition and wx.Platform == '__WXMAC__':
      pos = (0, 26)
    super(InspectorToolFrame, self).__init__(parent, id, title=" ",
                                             pos=pos, size=size,
                                             style=wx.CAPTION|wx.CLOSE_BOX)
    s = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(s)

  def move_parent_out_of_the_way(self):
    x,y,w,h = self.GetRect()
    xp,yp = self.Parent.GetPosition()
    if wx.Platform == '__WXGTK__': gap = 10
    else: gap = 5
    xp = max(xp, x + w + gap)
    self.Parent.Move((xp,yp))

  def Layout(self):
    s = self.Sizer
    for i in self.Children:
      s.Add(i, 0, wx.EXPAND|wx.BOTTOM, 5)
    w = max([ i.width_for_stacking() for i in self.Children ])
    s.SetMinSize((w,-1))
    s.SetSizeHints(self)


class Slider(wx.Slider):

  def __init__(self, parent, id=wx.ID_ANY,
               value=0, minValue=0, maxValue=100,
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style=wx.SL_HORIZONTAL,
               validator=wx.DefaultValidator,
               name=wx.SliderNameStr):
    kwds = group_args(**locals())
    w,h = kwds.size
    if w == -1: w = 150
    kwds.size = (w,h)
    wx.Slider.__init__(self, **kwds.__dict__)


if __name__ == '__main__':
  a = wx.App(redirect=False)

  w1 = wx.Frame(None, pos=(0, 0))
  s = wx.BoxSizer(wx.VERTICAL)
  w1.SetSizer(s)
  s.Add(Slider(w1, -1, 5, 0, 10, wx.DefaultPosition, (100,-1)))
  s.Add(Slider(w1, -1, 5, 0, 10, wx.DefaultPosition, (-1,-1)))
  s.Add(Slider(w1, size=(150,-1)))
  s.Add(Slider(w1, size=(-1,-1)))
  s.Add(Slider(w1))
  w1.Fit()
  a.SetTopWindow(w1)

  w = InspectorToolFrame(w1)

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

  w.move_parent_out_of_the_way()
  w.Show()
  w1.Show()

  a.MainLoop()
