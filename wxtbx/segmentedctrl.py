
# Copyright 2010 University of California
# All rights reserved

# This is an approximate mimic of the NSSegmentedControl in Cocoa.  The
# Apple control has many more styles than this, but the basic functionality
# is the same, with three selection modes (any, one, or none).  Either
# labels, bitmaps, or both are supported.

# TODO:
#   - Windows bugs
#   - additional styles?

import wx
import sys, os

SEGBTN_ROUNDED_CORNERS = 1
SEGBTN_VERTICAL = 2
SEGBTN_HORIZONTAL = 4
SEGBTN_INVERT_BITMAP = 8
SEGBTN_HIGHLIGHT = 16
SEGBTN_RADIO_ALLOW_NONE = 32

SEGBTN_DEFAULT_STYLE = SEGBTN_HORIZONTAL|SEGBTN_ROUNDED_CORNERS

if wx.Platform == '__WXMAC__' :
  _DEFAULT_PADDING = 8
elif wx.Platform == '__WXGTK__' :
  _DEFAULT_PADDING = 6
else :
  _DEFAULT_PADDING = 6

class SegmentedControl (wx.PyControl) :
  def __init__ (self,
                parent,
                id=wx.ID_ANY,
                pos=wx.DefaultPosition,
                size=wx.DefaultSize,
                style=SEGBTN_DEFAULT_STYLE,
                name=wx.ButtonNameStr,
                border=0,
                pad=_DEFAULT_PADDING) :
    wx.PyControl.__init__(self, parent, id, pos=pos, size=size,
      style=wx.NO_BORDER,
      name=name)
    self.SetBackgroundStyle(wx.BG_STYLE_CUSTOM)
    self.InheritAttributes()
    self.segments = []
    self.values = [] # for radio and toggle versions
    self._clicked_segment = None
    self._padding = _DEFAULT_PADDING
    self._border = border
    self._style = style
    self._border_color = (0, 0, 0)
    if wx.Platform == '__WXGTK__' :
      self.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))

    # Event Handlers
    self.Bind(wx.EVT_PAINT, lambda evt: self.__DrawButtons())
    self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnErase)
    #self.Bind(wx.EVT_SET_FOCUS, self.OnFocus)
    #self.Bind(wx.EVT_KILL_FOCUS, self.OnKillFocus)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    #self.Bind(wx.EVT_LEFT_DCLICK, lambda evt: self.ToggleState())
    self.Bind(wx.EVT_ENTER_WINDOW, self.OnEnter)
    self.Bind(wx.EVT_LEAVE_WINDOW, self.OnLeave)
    self.Bind(wx.EVT_MOTION, self.OnMotion)

  def AddSegment (self, label="", bitmap=None) :
    if label == "" and bitmap is None :
      raise RuntimeError("Must specify either a label or a bitmap or both "+
        "for SegmentControl.AddSegment().")
    self.segments.append((label, bitmap))
    self.values.append(False)

  def InsertSegment (self, index, label="", bitmap=None) :
    if label == "" and bitmap is None :
      raise RuntimeError("Must specify either a label or a bitmap or both "+
        "for SegmentControl.InsertSegment().")
    self.segments.insert(index, (label, bitmap))
    self.values.append(True)

  def Realize (self) :
    self.SetSize(self.DoGetBestSize())

  def SetBorderColor (self, color) :
    self._border_color = color

  def DoGetBestSize (self) :
    dc = wx.ClientDC(self)
    dc.SetFont(self.GetFont())
    n_seg = len(self.segments)
    buf = 1
    if n_seg == 0 :
      (w, h) = dc.GetTextExtent("XXX")
      w += self._padding * 2
      h += self._padding * 2
    else :
      txt_w, txt_h = self.__OverallSegmentContentSize(dc)
      seg_w = txt_w + (2 * self._padding)
      seg_h = txt_h + (2 * self._padding)
      if self._style & SEGBTN_VERTICAL :
        w = seg_w + 2
        h = 1 + ((seg_h + buf) * n_seg)
      else :
        w = 1 + ((seg_w + buf) * n_seg)
        h = seg_h + 2
    w += (self._border * 2)
    h += (self._border * 2)
    return (w, h)

  def __OverallSegmentContentSize (self, dc) :
    max_h = 0
    max_w = 0
    for (label, bitmap) in self.segments :
      (text_w, text_h) = (0, 0)
      (bmp_w, bmp_h) = (0, 0)
      if label != "" :
        (text_w, text_h) = dc.GetTextExtent(label)
      if bitmap is not None :
        (bmp_w, bmp_h) = bitmap.GetSize()
      if label != "" and bitmap is not None :
        seg_w = max(text_w, bmp_w)
        seg_h = text_h + bmp_h + 2
      elif label != "" :
        (seg_w, seg_h) = (text_w, text_h)
      else :
        (seg_w, seg_h) = (bmp_w, bmp_h)
      if seg_h > max_h :
        max_h = seg_h
      if seg_w > max_w :
        max_w = seg_w
    return (max_w, max_h)

  def __DrawButtons (self) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    if wx.Platform == '__WXMSW__' :
      dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
      dc.Clear()
    gc = wx.GCDC(dc)
    gc.SetBackgroundMode(wx.TRANSPARENT)
    gc.SetFont(self.GetFont())

    (w, h) = self.GetSize() #DoGetBestSize()
    w -= (self._border * 2)
    h -= (self._border * 2)
    n_seg = len(self.segments)
    rgc = gc.GetGraphicsContext()
    if self._style & SEGBTN_VERTICAL :
      # XXX a smaller gradient looks better for vertical controls
      gradient_brush = rgc.CreateLinearGradientBrush(0, 0, w, 0, (240,240,240),
        (200,200,200))
    else :
      gradient_brush = rgc.CreateLinearGradientBrush(0, 0, 0, h, (250,250,250),
        (185,185,185))
    rgc.SetBrush(gradient_brush)
    border_pen = wx.Pen(self._border_color, 1)
    # XXX wx.GCDC on GTK draws thin black lines as thicker dark grey lines
    # due to antialiasing.
    if wx.Platform == '__WXGTK__' :
      gc.SetPen(wx.TRANSPARENT_PEN)
      dc.SetPen(border_pen)
      dc.SetBrush(wx.TRANSPARENT_BRUSH)
    else :
      gc.SetPen(border_pen)
    if self._style & SEGBTN_ROUNDED_CORNERS :
      gc.DrawRoundedRectangle(self._border, self._border, w, h, 3)
      if wx.Platform == '__WXGTK__' :
        dc.DrawRoundedRectangle(self._border, self._border, w, h, 3)
    else :
      gc.DrawRectangle(self._border, self._border, w, h)
      if wx.Platform == '__WXGTK__' :
        dc.DrawRectangle(self._border, self._border, w, h)
    gc.SetPen(wx.TRANSPARENT_PEN)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    if n_seg == 0 :
      gc.SetTextForeground('red')
      gc.DrawText("XXX", self._padding, self._padding)
    else :
      gc.SetTextForeground('black')
      gc.SetPen(wx.Pen('black', 1)) #BLACK_PEN)
      buf = 1
      (txt_w, txt_h) = self.__OverallSegmentContentSize(gc)
      if self._style & SEGBTN_VERTICAL :
        seg_w = w - 2
        seg_h = ((h - 1) / n_seg) - buf
      else :
        seg_w = ((w - 1) / n_seg) - buf
        seg_h = h -  2
      current_x = 1 + self._border
      current_y = 1 + self._border
      for i, (label, bitmap) in enumerate(self.segments) :
        gc.SetPen(border_pen)
        dc.SetPen(border_pen)
        (txt_x, txt_y) = (self._padding + 1, self._padding + 1)
        if self._style & SEGBTN_VERTICAL :
          if i > 0 :
            dc.DrawLine(self._border, current_y, w + self._border - 1,
              current_y)
            current_y += 1
          txt_y += current_y
          txt_x += self._border
        else :
          if i > 0 :
            dc.DrawLine(current_x, self._border + 1, current_x,
              h + self._border - 1)
            current_x += 1
          txt_x += current_x
          txt_y += self._border
        if ((self.HasCapture() and i == self._clicked_segment) or
             self.values[i] == True) :
          xo = current_x + (seg_w / 2)
          yo = current_y + (seg_h / 2)
          radius  = max(seg_w, seg_h) - 1
          brush = rgc.CreateRadialGradientBrush(xo, yo, xo, yo, radius,
            (220,220,220), (120,120,120))
          rgc.SetBrush(brush)
          gc.SetPen(wx.TRANSPARENT_PEN)
          if self._style & SEGBTN_ROUNDED_CORNERS and (i == 0 or i == n_seg) :
            gc.DrawRoundedRectangle(current_x, current_y, seg_w, seg_h, 1)
          else :
            gc.DrawRectangle(current_x, current_y, seg_w, seg_h)
          #gc.DrawRectangle(current_x, current_y, total_w, total_h)
        rgc.SetBrush(gradient_brush)
        gc.SetBrush(wx.TRANSPARENT_BRUSH)
        gc.SetPen(wx.BLACK_PEN)
        label_w, label_h = gc.GetTextExtent(label)
        if bitmap is not None :
          (img_w, img_h) = bitmap.GetSize()
          img_x, img_y = txt_x, txt_y
          if label != "" :
            if wx.Platform == '__WXGTK__' :
              img_y -= 4
            extra_h_space = seg_w - img_w - (self._padding * 2)
            extra_v_space = seg_w - img_h - label_h - (self._padding * 2)
          else :
            extra_h_space = seg_w - img_w - (self._padding * 2)
            extra_v_space = seg_h - img_h - (self._padding * 2)
          if extra_h_space > 0 :
            img_x += extra_h_space / 2
          if extra_v_space > 0 :
            img_y += extra_v_space / 2
          gc.DrawBitmap(bitmap, img_x, img_y, bitmap.GetMask() != None)
          txt_y += img_h
        if label != "" :
          #print seg_w, label_w
          extra_w = seg_w - label_w - (self._padding*2)
          if extra_w > 0 :
            txt_x += extra_w / 2
          extra_h = seg_h - label_h - (self._padding*2)
          if extra_h > 0 :
            txt_y += extra_h / 2
          gc.DrawText(label, txt_x, txt_y)
        if self._style & SEGBTN_VERTICAL :
          current_y += seg_h
        else :
          current_x += seg_w

  def HandleClick (self) :
    pass

  def HitTest (self, event) :
    x, y = event.GetPositionTuple()
    dc = wx.ClientDC(self)
    n_seg = len(self.segments)
    buf = 1
    (w, h) = self.DoGetBestSize()
    if self._style & SEGBTN_VERTICAL :
      seg_w = w - 2
      seg_h = ((h - 1) / n_seg) - buf
    else :
      seg_w = ((w - 1) / n_seg) - buf
      seg_h = h -  2
    if n_seg == 0 :
      return None
    elif (x < 1 or x > (w - 1)) or (y < 1 or y > (h-1)) :
      return None
    elif self._style & SEGBTN_VERTICAL :
      current_y = 1
      for i in range(n_seg) :
        max_y = current_y + seg_h
        if y > current_y and y < max_y :
          return i
        current_y = max_y
    else :
      current_x = 1
      for i in range(n_seg) :
        max_x = current_x + seg_w -1
        if x >= current_x and x < max_x :
          return i
        current_x = max_x + 2
    return None

  def OnLeftDown (self, event) :
    segment_id = self.HitTest(event)
    self._have_mouse = True
    self._clicked_segment = segment_id
    self.CaptureMouse()
    self.SetFocus()
    self.Refresh()

  def OnLeftUp (self, event) :
    if self.HasCapture() :
      self.ReleaseMouse()
    segment_id = self.HitTest(event)
    if segment_id is not None and segment_id == self._clicked_segment :
      self.HandleClick()
    self._clicked_segment = None
    self.Refresh()

  def OnEnter (self, event) :
    if event.LeftDown() :
      self.CaptureMouse()
      self.Refresh()

  def OnLeave (self, event) :
    pass
    if self.HasCapture() :
      self.ReleaseMouse()
      self.Refresh()
#      self._clicked_segment = None

  def OnMotion (self, event) :
    return
    if not self.IsEnabled() :
      return
    x,y = event.GetPositionTuple()
    w,h = self.GetClientSizeTuple()
    if (x<0 or y<0 or x>=w or y>=h) and self.HasCapture() :
      #self._clicked_segment = None
      print -3
      self.ReleaseMouse()
      self.Refresh()
      return
    elif event.LeftDown() :
      print -2
      self.CaptureMouse()
      self.Refresh()
      return
    event.Skip()

  def OnErase (self, event) :
    pass

  def OnFocus (self, event) :
    pass

  def OnKillFocus (self, event) :
    pass

class SegmentedButtonControl (SegmentedControl) :
  def HandleClick (self) :
    bevt = wx.CommandEvent(wx.wxEVT_COMMAND_BUTTON_CLICKED, self.GetId())
    bevt.SetEventObject(self)
    bevt.SetString(self.GetLabel())
    bevt.SetClientData(self._clicked_segment)
    wx.PostEvent(self.GetParent(), bevt)
    self._clicked_segment = None

class SegmentedRadioControl (SegmentedControl) :
  def HandleClick (self) :
    index = self._clicked_segment
    assert index >= 0
    if self.values[index] == True and self._style & SEGBTN_RADIO_ALLOW_NONE :
      self.values[index] = False
    elif not self.values[index] :
      self.SetSelection(index)
    bevt = wx.CommandEvent(wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, self.GetId())
    bevt.SetEventObject(self)
    bevt.SetString(self.GetLabel())
    bevt.SetClientData(self.GetSelection())
    wx.PostEvent(self.GetParent(), bevt)
    self._clicked_segment = None

  def GetSelection (self) :
    assert self.values.count(True) <= 1
    if True in self.values :
      return self.values.index(True)
    return -1

  def SetSelection (self, index) :
    self.values = [False] * len(self.segments)
    self.values[index] = True

class SegmentedToggleControl (SegmentedControl) :
  def HandleClick (self) :
    index = self._clicked_segment
    assert index >= 0
    old_value = self.values[index]
    self.values[index] = not old_value
    bevt = wx.CommandEvent(wx.wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, self.GetId())
    bevt.SetEventObject(self)
    bevt.SetString(self.GetLabel())
    bevt.SetClientData(self._clicked_segment)
    wx.PostEvent(self.GetParent(), bevt)
    self._clicked_segment = None

  def GetValue (self, index) :
    return self.values[index]

  def SetValue (self, index, value) :
    self.values[index] = value

#-----------------------------------------------------------------------
if __name__ == "__main__" :
  def OnButton (evt) :
    print "Button clicked: %d" % evt.GetClientData()
  def OnRadio (evt) :
    print "Radio button selection: %d" % evt.GetClientData()
  def OnToggle (evt) :
    index = evt.GetClientData()
    print "Toggle button clicked: %d, %s" % (index,
      evt.GetEventObject().GetValue(index))
  try :
    from wxtbx import bitmaps
  except ImportError :
    bitmaps = None
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Test frame")
  toolbar = frame.CreateToolBar()
  panel = wx.Panel(frame, -1, size=(640,480))
  sizer = wx.BoxSizer(wx.HORIZONTAL)
  panel.SetSizer(sizer)
  v_sizer = wx.BoxSizer(wx.VERTICAL)
  btn = SegmentedButtonControl(panel,
    style=SEGBTN_HORIZONTAL) #|SEGBTN_AUTOSIZE)
  btn.AddSegment("Prev")
  btn.AddSegment("Home")
  btn.AddSegment("Search")
  btn.AddSegment("Help")
  btn.AddSegment("Next")
  frame.Bind(wx.EVT_BUTTON, OnButton, btn)
  v_sizer.Add(wx.StaticText(panel, -1, "Standard buttons"), 0, wx.ALL, 5)
  v_sizer.Add(btn, 0, wx.ALL, 5)
  sizer.Add(v_sizer)
  if True :
    tbtn = SegmentedButtonControl(toolbar,
      style=SEGBTN_HORIZONTAL|SEGBTN_ROUNDED_CORNERS,
      border=4)
    tbtn.AddSegment("Prev")
    tbtn.AddSegment("Home")
    tbtn.AddSegment("Search")
    tbtn.AddSegment("Help")
    tbtn.AddSegment("Next")
    tbtn.Realize()
    (w,h) = tbtn.GetSize()
    toolbar.SetToolBitmapSize((h,h))
    frame.Bind(wx.EVT_BUTTON, OnButton, tbtn)
    toolbar.AddControl(tbtn)
    #toolbar.AddControl(wx.Button(toolbar, -1, "Hello"))
    btn2 = SegmentedRadioControl(panel)
    btn2.AddSegment("Rotate")
    btn2.AddSegment("Pan")
    btn2.AddSegment("Zoom")
    btn2.SetSelection(0)
    frame.Bind(wx.EVT_RADIOBUTTON, OnRadio, btn2)
    v_sizer.Add(wx.StaticText(panel, -1, "Radio buttons"), 0, wx.ALL, 5)
    v_sizer.Add(btn2, 0, wx.ALL, 5)
    btn3 = SegmentedToggleControl(panel)
    btn3.AddSegment("Maps")
    btn3.AddSegment("Models")
    btn3.AddSegment("Selection")
    btn3.Realize()
    frame.Bind(wx.EVT_TOGGLEBUTTON, OnToggle, btn3)
    v_sizer.Add(wx.StaticText(panel, -1, "Toggle buttons"), 0, wx.ALL, 5)
    v_sizer.Add(btn3, 0, wx.ALL, 5)
    if (bitmaps is not None) and (bitmaps.icon_lib is not None) :
      bmp1 = bitmaps.fetch_icon_bitmap("actions", "1leftarrow", 16)
      bmp2 = bitmaps.fetch_icon_bitmap("actions", "gohome", 16)
      bmp3 = bitmaps.fetch_icon_bitmap("actions", "viewmag", 16)
      bmp4 = bitmaps.fetch_icon_bitmap("actions", "1rightarrow", 16)
      btn4 = SegmentedButtonControl(panel)
      btn4.AddSegment(bitmap=bmp1)
      btn4.AddSegment(bitmap=bmp2)
      btn4.AddSegment(bitmap=bmp3)
      btn4.AddSegment(bitmap=bmp4)
      v_sizer.Add(wx.StaticText(panel, -1, "Bitmap buttons"), 0, wx.ALL, 5)
      v_sizer.Add(btn4, 0, wx.ALL, 5)
      frame.Bind(wx.EVT_BUTTON, OnButton, btn4)
      btn5 = SegmentedButtonControl(panel,
        style=SEGBTN_ROUNDED_CORNERS)
      btn5.AddSegment(label="Back", bitmap=bmp1)
      btn5.AddSegment(label="Home", bitmap=bmp2)
      btn5.AddSegment(label="Search", bitmap=bmp3)
      btn5.AddSegment(label="Next", bitmap=bmp4)
      btn5.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))
      v_sizer.Add(wx.StaticText(panel, -1, "Bitmaps and text"), 0, wx.ALL, 5)
      v_sizer.Add(btn5, 0, wx.ALL, 5)
      frame.Bind(wx.EVT_BUTTON, OnButton, btn5)
    vbtn1 = SegmentedButtonControl(panel,
      style=SEGBTN_VERTICAL)
    vbtn1.AddSegment("Up")
    vbtn1.AddSegment("Home")
    vbtn1.AddSegment("Search")
    vbtn1.AddSegment("Down")
    frame.Bind(wx.EVT_BUTTON, OnButton, vbtn1)
    sizer.Add(vbtn1, 0, wx.ALL, 20)
    vbtn2 = SegmentedButtonControl(panel,
      style=SEGBTN_VERTICAL|SEGBTN_ROUNDED_CORNERS,
      border=8)
    vbtn2.AddSegment("Up")
    vbtn2.AddSegment("Home")
    vbtn2.AddSegment("Search")
    vbtn2.AddSegment("Down")
    frame.Bind(wx.EVT_BUTTON, OnButton, vbtn2)
    sizer.Add(vbtn2, 0, wx.ALL, 20)
  sizer.Layout()
  sizer.Fit(panel)
  frame.Fit()
  toolbar.Realize()
  frame.Show()
  app.MainLoop()
