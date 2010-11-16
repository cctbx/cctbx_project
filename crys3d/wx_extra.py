from __future__ import division
import wx

from libtbx import copy_init_args
import libtbx.object_oriented_patterns as oop

class _extended_wxDC(oop.injector, wx.DC):

  def FillWith3DGradient(self, rect, colour, direction, step=1):
    """ Based on Horst Puschmann's ImageTools.gradient_bgr in Olex """
    x,y,w,h = rect
    if direction == wx.SOUTH:
      d = h
      draw_line = lambda i: self.DrawLine(x, y+i, x+w, y+i)
    elif direction == wx.EAST:
      d = w
      draw_line = lambda i: self.DrawLine(x+i, y, x+i, y+h)
    elif direction == wx.NORTH:
      d = h
      draw_line = lambda i: self.DrawLine(x, y+d-1-i, x+w, y+d-1-i)
    elif direction == wx.WEST:
      d = w
      draw_line = lambda i: self.DrawLine(x+w-1-i, y, x+w-1-i, y+h)
    slope_breaks       = (0,     d//10+1,     d//5+1,       d)
    step_adjustements  = (  0.6,          1.2,        1.4)
    slopes = [ x*step/d for x in step_adjustements ]
    red_green_slopes = [ x*58 for x in slopes ]
    blue_slopes      = [ x*44 for x in slopes ]
    ranges = [ xrange(slope_breaks[i], slope_breaks[i+1])
               for i in xrange(len(slope_breaks)-1) ]
    for range, red_green_slope, blue_slope in zip(
      ranges, red_green_slopes, blue_slopes):
      for i in range:
        r = int(colour.Red()   - i*red_green_slope)
        g = int(colour.Green() - i*red_green_slope)
        b = int(colour.Blue()  - i*blue_slope)
        self.SetPen(wx.Pen(wx.Colour(r,g,b)))
        draw_line(i)


class MouseClickButtonMixin(object):

  def __init__(self):
    self._pressing = False
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.colour = self.normal_colour
    self.Bind(wx.EVT_LEAVE_WINDOW, self.OnLeaveWindow)
    self.Bind(wx.EVT_ENTER_WINDOW, self.OnEnterWindow)

  def OnLeftDown(self, event):
    if not self.IsEnabled(): return
    self._pressing = True
    self.CaptureMouse()
    self.colour = self.pressed_colour
    self.Refresh()

  def OnLeftUp(self, event):
    if not self.IsEnabled(): return
    if self.HasCapture():
      self.ReleaseMouse()
      if self._pressing and self.Rect.Contains(event.Position):
        self.OnClick(event)
      self._pressing = False
      self.colour = self.normal_colour
      self.Refresh()

  def OnLeaveWindow(self, event):
    if not self.IsEnabled(): return
    if self.HasCapture() and self._pressing:
      self.colour = self.normal_colour
      self.Refresh()

  def OnEnterWindow(self, event):
    if not self.IsEnabled(): return
    if self.HasCapture() and self._pressing:
      self.colour = self.pressed_colour
      self.Refresh()


class InspectorHeader(wx.PyControl, MouseClickButtonMixin):

  normal_colour = wx.Colour(237, 237, 235)
  pressed_colour = wx.Colour(200, 200, 198)
  horizontal_margin = 2
  vertical_margin = 2

  def __init__(self, parent, label=""):
    wx.PyControl.__init__(self, parent, style=wx.BORDER_NONE)
    MouseClickButtonMixin.__init__(self)

    self.Label = label
    self.SetFont(wx.NORMAL_FONT)

    self.Bind(wx.EVT_SIZING, self.OnResize)
    self.Bind(wx.EVT_SIZE, self.OnResize)

    # to reduce flicker: use BufferedPaintDC in handling EVT_PAINT
    # and ignore EVT_ERASE_BACKGROUND
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_ERASE_BACKGROUND, lambda event: None)

    self._opened = True

  def DoGetBestSize(self):
    try:
      return self._best_size
    except AttributeError:
      dc = wx.ClientDC(self)
      dc.SetFont(self.Font)
      w_txt, h_txt = dc.GetTextExtent(self.Label)
      self.triangle_base_length = 2*(h_txt//2) - 4
      h_txt += 2*self.vertical_margin
      w_txt += 4*self.vertical_margin + self.triangle_base_length
      self._best_size = (w_txt, h_txt)
      return self._best_size

  def OnPaint(self, event):
    dc = wx.BufferedPaintDC(self)
    dc.SetFont(self.Font)
    dc.FillWith3DGradient(self.Rect, self.colour,
                          direction=wx.SOUTH, step=0.6)
    w_t = self.triangle_base_length
    h_t = w_t//2
    h = self.GetSize()[1]
    dc.SetBrush(wx.BLACK_BRUSH)
    dc.SetPen(wx.BLACK_PEN)
    if self._opened:
      dc.DrawPolygon((wx.Point(0, -h_t), wx.Point(w_t, -h_t),
                              wx.Point(h_t, 0) ),
                     xoffset=self.horizontal_margin,
                     yoffset=(h-h_t)//2 + h_t)
    else:
      dc.DrawPolygon((wx.Point(0, -w_t), wx.Point(0, 0),
                              wx.Point(h_t, -h_t) ),
                     xoffset=self.horizontal_margin + h_t,
                     yoffset=(h-w_t)//2 + w_t)
    dc.DrawText(self.Label, w_t + 3*self.horizontal_margin,
                            self.vertical_margin)
    dc.Destroy()

  def OnResize(self, event):
    event.Skip()
    self.Refresh()

  def OnClick(self, event):
    self._opened = not self._opened
    self.on_click()


class Inspector(wx.PyPanel):

  def __init__(self, parent, id=wx.ID_ANY, label="", colour=wx.BLACK,
               pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
    assert parent.GetSizer() is not None
    super(Inspector, self).__init__(parent, id)
    self._header = InspectorHeader(self, label=label)
    self._header.on_click = self.Toggle
    self._pane = wx.Panel(self)
    self._expanded = True
    s = wx.BoxSizer(wx.VERTICAL)
    s.Add(self._header, 0, wx.EXPAND, 5)
    s.Add(self._pane, flag=wx.LEFT|wx.RIGHT|wx.TOP, border=5)
    self.SetSizer(s)

  def IsExpanded(self):
    return self._expanded

  def Collapse(self, f):
    self._expanded = not f
    self._update()

  def Toggle(self):
    self._expanded = not self._expanded
    self._update()

  def DoLayout(self):
    self.Sizer.Layout()

  def _update(self):
    if self._expanded:
      self._pane.Show()
    else:
      self._pane.Hide()
    top = wx.GetTopLevelParent(self)
    top_sizer = top.GetSizer()
    top_sizer.SetSizeHints(top)
    if self._expanded:
      top.Fit()
    else:
      sz = top_sizer.CalcMin()
      top.SetClientSize(sz)

  def GetLabel(self):
    return self._header.GetLabel()

  def GetPane(self):
    return self._pane

  def DoGetBestSize(self):
    return self.Sizer.CalcMin()


class InspectorToolFrame(wx.MiniFrame):

  def __init__(self, parent, id=wx.ID_ANY,
               pos=wx.DefaultPosition, size=wx.DefaultSize):
    if pos is wx.DefaultPosition and wx.Platform == '__WXMAC__':
      pos = (0, 26)
    super(InspectorToolFrame, self).__init__(
      parent, id, title=" ",
      pos=pos, size=size,
      style=wx.CAPTION|wx.CLOSE_BOX)
    s = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(s)
    self.Bind(wx.EVT_CLOSE, self.on_close)

  def on_close(self, event):
    self.Hide()

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
      s.Add(i, 0, wx.EXPAND)
    w = max([ i.GetBestSize()[0] for i in self.Children ])
    s.SetMinSize((w,-1))
    s.SetSizeHints(self)


class Slider(wx.Slider):

  def __init__(self, parent, id=wx.ID_ANY,
               value=0, minValue=0, maxValue=100,
               pos=wx.DefaultPosition, size=wx.DefaultSize,
               style=wx.SL_HORIZONTAL,
               validator=wx.DefaultValidator,
               name=wx.SliderNameStr):
    kwds = copy_init_args(locals())
    w,h = kwds.size
    if w == -1: w = 150
    kwds.size = (w,h)
    wx.Slider.__init__(self, **kwds.__dict__)


if __name__ == '__main__':
  a = wx.App(redirect=False)

  class GradientTestFrame(wx.Frame):

    def __init__(self):
      super(GradientTestFrame, self).__init__(None, pos=(0,-1))
      size = (120, 90)
      box = wx.BoxSizer(wx.VERTICAL)
      for direction in (wx.NORTH, wx.EAST, wx.SOUTH, wx.WEST):
        p = wx.Window(self, style=wx.BORDER_NONE)
        p.SetMinSize(size)
        box.Add(p)
        p.Bind(wx.EVT_PAINT, self.make_on_paint(p, direction))
      self.SetSizer(box)
      box.SetSizeHints(self)
      self.Show()

    def make_on_paint(window, direction):
      def on_paint(event):
        dc = wx.PaintDC(window)
        dc.FillWith3DGradient(window.ClientRect,
                              colour=wx.Colour(237, 237, 235),
                              direction=direction)
        del dc
      return on_paint
    make_on_paint = staticmethod(make_on_paint)

  gradient_test_frame = GradientTestFrame()


  class InspectorHeaderTestFrame(wx.Frame):

    def __init__(self):
      super(InspectorHeaderTestFrame, self).__init__(None, pos=(300,-1))
      box = wx.BoxSizer(wx.VERTICAL)
      inspector = InspectorHeader(self, label="Title")
      box.Add(inspector, flag=wx.EXPAND)
      self.log = wx.TextCtrl(self, style=wx.TE_MULTILINE|wx.BORDER_NONE,
                                   size=(80, 240))
      self.i = 0
      inspector.on_click = self.log_click
      box.Add(self.log, proportion=1, flag=wx.ALL|wx.EXPAND, border=5)
      self.SetSizer(box)
      box.SetSizeHints(self)
      self.Show()

    def log_click(self):
      self.i += 1
      self.log.AppendText("Click #%i\n" % self.i)

  inspector_header_test_frame = InspectorHeaderTestFrame()

  w1 = wx.Frame(None, pos=(0, 500))
  s = wx.BoxSizer(wx.VERTICAL)
  w1.SetSizer(s)
  s.Add(Slider(w1, -1, 5, 0, 10, wx.DefaultPosition, (100,-1)))
  s.Add(Slider(w1, -1, 5, 0, 10, wx.DefaultPosition, (-1,-1)))
  s.Add(Slider(w1, size=(150,-1)))
  s.Add(Slider(w1, size=(-1,-1)))
  s.Add(Slider(w1))
  w1.Fit()
  a.SetTopWindow(w1)

  w = InspectorToolFrame(w1, pos=(100, 500))

  i1 = Inspector(w, label="Header")
  def f():
    print "Click in 'Header'"
  i1.on_click_in_header = f
  pane = i1.GetPane()
  top = wx.BoxSizer(wx.VERTICAL)
  s = wx.BoxSizer(wx.HORIZONTAL)
  s.Add(wx.Button(pane, label="Click me!"))
  s.Add(wx.CheckBox(pane, label="Check me!"))
  top.Add(s, flag=wx.ALL, border=5)
  pane.SetSizer(top)

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
