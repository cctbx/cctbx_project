
# TODO:
#  - cached scenes

from crys3d import hklview
import wx
from math import sqrt
import sys

class hklview_2d (wx.PyPanel) :
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

  def set_miller_array (self, array) :
    self.miller_array = array
    if (array is not None) :
      self.construct_reciprocal_space()

  def construct_reciprocal_space (self) :
    self.scene = hklview.scene(miller_array=self.miller_array,
      settings=self.settings)
    self._clicked = None

  def update_settings (self) :
    self.construct_reciprocal_space()
    self.Refresh()

  def get_center_and_radius (self) :
    w, h = self.GetSize()
    r = (min(w,h) // 2) - 20
    center_x = max(w // 2, r + 20)
    center_y = max(h // 2, r + 20)
    return center_x, center_y, r

  def paint (self, gc) :
    self._points_2d = []
    self._radii_2d = []
    assert (self.settings.slice_mode)
    if (self.settings.slice_axis == 0) :
      i_x, i_y = 1, 2
      axes = ("k", "l")
    elif (self.settings.slice_axis == 1) :
      i_x, i_y = 0, 2
      axes = ("h", "l")
    else :
      i_x, i_y = 0, 1
      axes = ("h", "k")
    center_x, center_y, r = self.get_center_and_radius()
    x_max = self.scene.axes[i_x][i_x] * 100.
    y_max = self.scene.axes[i_y][i_y] * 100.
    font = self.GetFont()
    font.SetFamily(wx.FONTFAMILY_MODERN)
    if (self.settings.black_background) :
      gc.SetFont(gc.CreateFont(font, (255,255,255)))
    else :
      gc.SetFont(gc.CreateFont(font, (0,0,0)))
    if (self.settings.show_axes) :
      # FIXME dimensions not right?
      x_end = self.scene.axes[i_x][i_x], self.scene.axes[i_x][i_y]
      y_end = self.scene.axes[i_y][i_x], self.scene.axes[i_y][i_y]
      x_len = sqrt(x_end[0]**2 + x_end[1]**2)
      y_len = sqrt(y_end[0]**2 + y_end[1]**2)
      x_scale = (r+10) / x_len
      y_scale = (r+10) / y_len
      x_end = (x_end[0] * x_scale, x_end[1] * x_scale)
      y_end = (y_end[0] * y_scale, y_end[1] * y_scale)
      x_axis = gc.CreatePath()
      x_axis.MoveToPoint(center_x, center_y)
      x_axis.AddLineToPoint(center_x + x_end[0], center_y - x_end[1])
      x_axis.CloseSubpath()
      #x_axis.MoveToPoint(center_x + r + 15, center_y - 5)
      #x_axis.AddLineToPoint(center_x + r + 20, center_y)
      #x_axis.CloseSubpath()
      #x_axis.MoveToPoint(center_x + r + 20, center_y)
      #x_axis.AddLineToPoint(center_x + r + 15, center_y + 5)
      #x_axis.CloseSubpath()
      if (self.settings.black_background) :
        gc.SetPen(wx.Pen('white'))
      else :
        gc.SetPen(wx.Pen('black'))
      gc.PushState()
      gc.StrokePath(x_axis)
      gc.PopState()
      y_axis = gc.CreatePath()
      y_axis.MoveToPoint(center_x, center_y)
      y_axis.AddLineToPoint(center_x + y_end[0], center_y - y_end[1])
      y_axis.CloseSubpath()
      #y_axis.MoveToPoint(center_x - 5, center_y - r - 15)
      #y_axis.AddLineToPoint(center_x, center_y - r - 20)
      #y_axis.CloseSubpath()
      #y_axis.MoveToPoint(center_x, center_y - r - 20)
      #y_axis.AddLineToPoint(center_x + 5, center_y - r - 15)
      #y_axis.CloseSubpath()
      gc.PushState()
      gc.StrokePath(y_axis)
      gc.PopState()
      gc.DrawText(axes[0], center_x + x_end[0] - 6, center_y - x_end[1] - 20)
      gc.DrawText(axes[1], center_x + y_end[0] + 6, center_y - y_end[1])
    gc.SetPen(wx.TRANSPARENT_PEN)
    if (self.settings.black_background) :
      main_pen = wx.Pen('white')
      main_brush = wx.Brush('white')
      missing_brush = wx.Brush((1,1,1))
      if (self.settings.color_scheme != 0) :
        missing_pen = wx.Pen('red')
      else :
        missing_pen = wx.Pen('white')
    else :
      main_pen = wx.Pen('black')
      main_brush = wx.Brush('black')
      missing_brush = wx.Brush((250,250,250))
      if (self.settings.color_scheme != 0) :
        missing_pen = wx.Pen('red')
      else :
        missing_pen = wx.Pen('black')
    max_radius = 40 / max(self.scene.unit_cell.parameters()[0:3])
    max_radius *= r / max(x_max, y_max)
    for k, hkl in enumerate(self.scene.points) :
      x_, y_ = hkl[i_x], hkl[i_y]
      x = center_x + r * x_ / x_max
      y = center_y - r * y_ / y_max
      r_point = self.scene.radii[k] * r / max(x_max, y_max)
      if (self.settings.uniform_size) :
        r_point = max_radius
      elif (self.settings.pad_radii) :
        r_point = min(1, r_point)
      self._points_2d.append((x,y))
      self._radii_2d.append(r_point)
      path = gc.CreatePath()
      path.AddCircle(0, 0, r_point)
      path.CloseSubpath()
      gc.PushState()
      gc.Translate(x,y)
      if (self.scene.missing[k]) :
        gc.SetBrush(wx.TRANSPARENT_BRUSH)
        gc.SetPen(missing_pen)
        gc.StrokePath(path)
      else :
        c = self.scene.colors[k]
        gc.SetBrush(wx.Brush((c[0]*255,c[1]*255,c[2]*255)))
        gc.FillPath(path)
      gc.PopState()
    if (self._clicked is not None) :
      gc.SetPen(main_pen)
      gc.DrawText("Clicked: ", 10, 10)
      w, h = gc.GetTextExtent("Clicked: ")
      gc.DrawText("d = %g" % self.scene.get_resolution_at_point(self._clicked),
        w+10, 30)
      if (self.settings.color_scheme == 0) :
        c = self.scene.colors[self._clicked]
        c = (c[0]*255, c[1]*255, c[2]*255)
        gc.SetPen(wx.Pen(c))
        gc.SetFont(gc.CreateFont(self.GetFont(),c))
      gc.DrawText("%d,%d,%d" % self.scene.indices[self._clicked], w+10, 10)

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
    hkl = resolution = None
    if (self._clicked is not None) :
      hkl = self.scene.indices[self._clicked]
      resolution = self.scene.get_resolution_at_point(self._clicked)
    self.GetParent().update_clicked(hkl, resolution)
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
    memory_dc.SetBackgroundMode(wx.TRANSPARENT)
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
