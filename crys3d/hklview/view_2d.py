
# TODO:
#  - cached scenes

from crys3d import hklview
import cctbx.miller.display
import wx
from math import sqrt
import sys

class hklview_2d (wx.PyPanel, cctbx.miller.display.render_2d) :
  def __init__ (self, *args, **kwds) :
    wx.PyPanel.__init__(self, *args, **kwds)
    font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    self.SetFont(font)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftClick, self)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp, self)
    self.Bind(wx.EVT_MOTION, self.OnMouseMotion, self)
    self.scene = None
    self.miller_array = None
    self.parent = self.GetParent()
    self.settings = self.parent.settings
    self.was_dragged = False
    self.initLeft = None, None
    self._points_2d = []
    self._radii_2d = []
    self._clicked = None

  def GetSize (self) :
    return wx.PyPanel.GetSize(self)

  # XXX silent keyword 'zoom=False' is for compatibility with view_3d.py
  def set_miller_array (self, array, zoom=False) :
    self.miller_array = array
    if (array is not None) :
      self.construct_reciprocal_space()

  def construct_reciprocal_space (self) :
    self.scene = hklview.scene(miller_array=self.miller_array,
      settings=self.settings)
    self._clicked = None
    self.setup_colors()

  def update_settings (self) :
    self.construct_reciprocal_space()
    self.Refresh()

  def get_color (self, c) :
    return (int(c[0]*255), int(c[1]*255), int(c[2]*255))

  def draw_line (self, canvas, x1, y1, x2, y2) :
    gc = canvas
    x_axis = gc.CreatePath()
    x_axis.MoveToPoint(x1, y1)
    x_axis.AddLineToPoint(x2, y2)
    x_axis.CloseSubpath()
    gc.SetPen(wx.Pen(self.get_color(self._foreground)))
    gc.PushState()
    gc.StrokePath(x_axis)
    gc.PopState()

  def draw_text (self, canvas, text, x, y) :
    gc = canvas
    gc.SetPen(wx.Pen(self.get_color(self._foreground)))
    gc.DrawText(text, x, y)

  def draw_open_circle (self, canvas, x, y, radius, color=None) :
    gc = canvas
    path = gc.CreatePath()
    path.AddCircle(0, 0, radius)
    path.CloseSubpath()
    gc.PushState()
    gc.Translate(x,y)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    if (color is None) :
      color = self._foreground
    pen = wx.Pen(self.get_color(color))
    gc.SetPen(pen)
    gc.StrokePath(path)
    gc.PopState()

  def draw_filled_circle (self, canvas, x, y, radius, color) :
    gc = canvas
    path = gc.CreatePath()
    path.AddCircle(0, 0, radius)
    path.CloseSubpath()
    gc.PushState()
    gc.Translate(x,y)
    if (color is None) :
      color = self._foreground
    pen = wx.Pen(self.get_color(color))
    brush = wx.Brush(self.get_color(color))
    gc.SetPen(pen)
    gc.SetBrush(brush)
    gc.FillPath(path)
    gc.PopState()

  def paint (self, gc) :
    font = self.GetFont()
    font.SetFamily(wx.FONTFAMILY_MODERN)
    if (self.settings.black_background) :
      gc.SetFont(gc.CreateFont(font, (255,255,255)))
    else :
      gc.SetFont(gc.CreateFont(font, (0,0,0)))
    self.render(gc)

  def save_screen_shot (self, **kwds) :
    pass

  def process_pick_points (self, x, y) :
    context = wx.ClientDC( self )
    w, h = self.GetClientSize()
    bitmap = wx.EmptyBitmap( w, h, -1 )
    memory = wx.MemoryDC(bitmap)
    memory.SelectObject(bitmap)
    memory.Blit(0, 0, w, h, context, 0, 0)
    memory.SelectObject(wx.NullBitmap)
    if (wx.Platform == '__WXMAC__') :
      pixelData = wx.AlphaPixelData(bitmap)
      pixelAccessor = pixelData.GetPixels()
      pixelAccessor.MoveTo(pixelData, x, y)
      c = pixelAccessor.Get()
    else :
      c = memory.GetPixel(x, y)
    bg = self.GetBackgroundColour()
    self._clicked = None
    min_dist = sys.maxint
    closest_hkl = None
    for k, (x2,y2) in enumerate(self._points_2d) :
      dist = sqrt((x2-x)**2 + (y2-y)**2)
      if (dist <= (self._radii_2d[k] + 2)) :
        self._clicked = k
        break
    hkl = d_min = value = None
    if (self._clicked is not None) :
      hkl, d_min, value = self.scene.get_reflection_info(self._clicked)
    self.GetParent().update_clicked(hkl, d_min, value)
    self.Refresh()

  # placeholder - mimics gltbx.wx_viewer.wxGLWindow.fit_into_viewport
  def fit_into_viewport (self) :
    pass

  # mimics gltbx.wx_viewer.wxGLWindow.save_screen_shot
  def save_screen_shot (self, file_name, extensions=None) :
    rect = self.GetRect()
    bitmap = wx.EmptyBitmap(rect.width, rect.height)
    memory_dc = wx.MemoryDC()
    memory_dc.SelectObject(bitmap)
    #memory_dc.SetBackgroundMode(wx.TRANSPARENT)
    if (self.settings.black_background) :
      memory_dc.SetBackground(wx.BLACK_BRUSH)
    else :
      memory_dc.SetBackground(wx.WHITE_BRUSH)
    memory_dc.Clear()
    gc = wx.GraphicsContext.Create(memory_dc)
    self.paint(gc)
    bitmap.SaveFile(file_name, wx.BITMAP_TYPE_PNG)

  def OnPaint (self, event) :
    if (self.scene is None) :
      return
    if (self.settings.black_background) :
      self.SetBackgroundColour((0,0,0))
    else :
      self.SetBackgroundColour((255,255,255))
    dc = wx.AutoBufferedPaintDCFactory(self)
    if (self.settings.black_background) :
      dc.SetBackground(wx.BLACK_BRUSH)
    else :
      dc.SetBackground(wx.WHITE_BRUSH)
    dc.Clear()
    gc = wx.GraphicsContext.Create(dc)
    self.paint(gc)

  def OnLeftClick (self, evt) :
    self.initLeft = evt.GetX(), evt.GetY()

  def OnLeftUp (self, evt) :
    x = evt.GetX()
    y = evt.GetY()
    if (not self.was_dragged) :
      if (x == self.initLeft[0]) and (y == self.initLeft[1]) :
        self.process_pick_points(x,y)
    self.was_dragged = False

  def OnMouseMotion (self, evt) :
    if (not evt.Dragging()) :
      return
    elif (evt.LeftIsDown()) :
      self.was_dragged = True

  def OnChar (self, evt) :
    pass
