
import wx
import mmtbx.polygon.output
from math import radians

class wx_renderer (mmtbx.polygon.output.renderer) :
  def draw_bin (self, out, start, end, angle, color) :
    gc = out
    path = gc.CreatePath()
    path.MoveToPoint(start[0], start[1])
    path.AddLineToPoint(end[0], end[1])
    path.CloseSubpath()
    gc.PushState()
    gc.SetPen(wx.Pen(color, 10))
    gc.StrokePath(path)
    gc.PopState()

  def draw_box (self, out, points, color) :
    gc = out
    path = gc.CreatePath()
    path.MoveToPoint(points[0][0], points[0][1])
    path.AddLineToPoint(points[1][0], points[1][1])
    path.AddLineToPoint(points[2][0], points[2][1])
    path.AddLineToPoint(points[3][0], points[3][1])
    path.AddLineToPoint(points[0][0], points[0][1])
    path.CloseSubpath()
    gc.PushState()
    gc.SetPen(wx.Pen(color, 1)) #TRANSPARENT_PEN)
    gc.SetBrush(wx.Brush(color))
    gc.FillPath(path)
    gc.PopState()

  def draw_solid_line (self, out, start, end, color) :
    gc = out
    line = gc.CreatePath()
    line.MoveToPoint(start[0], start[1])
    line.AddLineToPoint(end[0], end[1])
    line.CloseSubpath()
    gc.PushState()
    if self.color_model == "gray" :
      gc.SetPen(wx.Pen("red", 2))
    else :
      gc.SetPen(wx.Pen("black", 2))
    gc.StrokePath(line)
    gc.PopState()

  def draw_dashed_line (self, out, start, end, color) :
    pass

  def draw_labels (self, out, label, min, max, value, pos, angle) :
    gc = out
    label_font = wx.Font(14, wx.NORMAL, wx.NORMAL, wx.BOLD)
    stat_font = wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL)
    gc.PushState()
    gc.SetPen(wx.Pen("black", 1))
    gc.SetFont(label_font)
    (text_w, text_h) = gc.GetTextExtent(label)
    (anchor_x, anchor_y) = pos
    (text_x, text_y) = self.get_text_position(text_w, text_h, anchor_x,
      anchor_y, angle)
    gc.DrawText(label, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.RED))
    (text_w, text_h) = gc.GetTextExtent(min)
    text_y += text_h + 2
    gc.DrawText(min, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.BLACK))
    (text_w, text_h) = gc.GetTextExtent(value)
    text_y += text_h + 2
    gc.DrawText(value, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.RED))
    (text_w, text_h) = gc.GetTextExtent(max)
    text_y += text_h + 2
    gc.DrawText(max, text_x, text_y)
    gc.PopState()

  def get_text_position (self, w, h, x, y, angle) :
    if angle >= radians(60) and angle < radians(120) :
      text_x = x - (w/2) - 5
      text_y = y - h - 15
    elif angle >= radians(120) and angle < radians(240) :
      text_x = x - w - 15
      text_y = y - (h/2)
    elif angle >= radians(240) and angle < radians(300) :
      text_x = x - (w/2)
      text_y = y
    else : # 300 =< angle < 420
      text_x = x + 5
      text_y = y - (h/2)
    return (text_x, text_y)

class PolygonPanel (wx.Panel) :
  def __init__ (self, parent, renderer) :
    wx.Panel.__init__(self, parent, -1)
    self.renderer = renderer
    self.renderer.resize((640, 640))
    self.Bind(wx.EVT_PAINT, self.OnPaint)

  def OnPaint (self, event) :
    self.renderer.resize(self.GetSize())
    dc = wx.PaintDC(self)
    gc = wx.GraphicsContext.Create(dc)
    self.renderer.draw(gc)

  def draw_color_key (self, dc) :
    gc = wx.GraphicsContext.Create(dc)
    stat_font = wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL)
    x = 40
    y = self.h - 10
    i = 0
    for shade in bin_colors :
      gc.PushState()
      gc.SetPen(wx.Pen(shade, 10))
      path = gc.CreatePath()
      path.MoveToPoint(x, y)
      path.AddLineToPoint(x + 40, y)
      path.CloseSubpath()
      gc.StrokePath(path)
      gc.PopState()
      if i < len(self.cutoffs) :
        gc.PushState()
        gc.SetFont(gc.CreateFont(stat_font, wx.BLACK))
        gc.DrawText(str(self.cutoffs[i]), x + 50, y - 6)
        gc.PopState()
      x += 80
      i += 1

  def OnChar (self, event) :
    keycode = event.GetKeyCode()
    if keycode == 32 :
      self.OnSave()

  def OnSave (self, event=None) :
    rect = self.GetRect()
    bitmap = wx.EmptyBitmap(rect.width, rect.height)
    memory_dc = wx.MemoryDC()
    memory_dc.SelectObject(bitmap)
    memory_dc.SetBackgroundMode(wx.TRANSPARENT)
    gc = wx.GraphicsContext.Create(memory_dc)
    self.renderer.draw(gc)
    output_file = wx.FileSelector("Save image as:",
      default_filename="polygon.png",
      wildcard="PNG image (*.png)|*.png", flags=wx.SAVE)
    if output_file != "" :
      bitmap.SaveFile(output_file, wx.BITMAP_TYPE_PNG)
    if event is not None :
      event.Skip()

  def reset_layout (self) :
    pass
