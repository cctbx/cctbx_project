
# TODO:
#  - handle 2theta properly
#  - resolution circles (2theta-dependent)
#  - measure reciprocal-space distance between spots

from rstbx.viewer import screen_params
import wx
user_callback = None

class XrayView (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    self._img = None
    super(XrayView, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.SetupEventHandlers()
    self.xmouse = None
    self.ymouse = None
    self.line_start = None
    self.line_end = None
    self.was_dragged = False
    self.shift_was_down = False
    self._last_zoom = 0
    self.zoom_level = None
    # miscellaneous non-user flags
    self.flag_always_show_predictions = False
    self.flag_spots_as_points = False

  def SetupEventHandlers (self) :
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.Bind(wx.EVT_MOTION, self.OnMotion)
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
    self.Bind(wx.EVT_MIDDLE_DOWN, self.OnMiddleDown)
    self.Bind(wx.EVT_MIDDLE_UP, self.OnMiddleUp)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
    self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
    self.Bind(wx.EVT_ENTER_WINDOW, self.OnEnter)
    self.Bind(wx.EVT_LEAVE_WINDOW, self.OnLeave)

  def set_image (self, image) :
    self._img = image
    self._img.set_screen_size(*(self.GetSize()))
    self.update_settings()

  def get_scale (self) :
    if (self.zoom_level is not None) :
      return self.zoom_level
    else :
      return self._img.get_scale()

  def update_settings (self, layout=True) :
    self.line = None
    scales = [0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
    zoom = scales[self.settings.zoom_level]
    self._img.set_zoom(zoom)
    self._img.update_settings(
      brightness=self.settings.brightness)
    if (layout) :
      self.line_start = None
      self.line_end = None
      self.OnSize(None)
    self.Refresh()
    if (self.GetParent().zoom_frame is not None) :
      self.GetParent().zoom_frame.Refresh()
    self.GetParent().settings_frame.refresh_thumbnail()

  # EVENTS
  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    w, h = self.GetSize()
    bitmap = self._img.get_bitmap()
    x, y = self._img.adjust_screen_coordinates(0, 0)
    dc.DrawBitmap(bitmap, x, y)
    if (self.settings.show_beam_center) :
      center_x, center_y = self._img.get_beam_center()
      xc, yc = self._img.image_coords_as_screen_coords(center_x, center_y)
      if (xc < w) and (yc < h) :
        dc.SetPen(wx.Pen('blue'))
        dc.DrawLine(xc - 10, yc, xc + 10, yc)
        dc.DrawLine(xc, yc - 10, xc, yc + 10)
    if (self.line_start is not None) and (self.line_end is not None) :
      dc.SetPen(wx.Pen('red', 2, wx.DOT))
      x1, y1 = self._img.image_coords_as_screen_coords(*(self.line_start))
      x2, y2 = self._img.image_coords_as_screen_coords(*(self.line_end))
      dc.DrawLine(x1, y1, x2, y2)
    if (self.settings.show_spotfinder_spots) :
      self.draw_spotfinder_spots(dc)
    if (self.settings.show_integration) :
      self.draw_integration_results(dc)
    if user_callback != None:
      user_callback(dc,self,wx)

  def draw_spotfinder_spots (self, dc) :
    spots = self._img.get_drawable_spots()
    dc.SetPen(wx.Pen('red'))
    if (self.flag_spots_as_points) :
      dc.DrawPoint(x, y)
    else :
      spot_scale = self._img.get_scale() * 5
      for x,y in spots :
        dc.DrawLine(x-spot_scale, y, x+spot_scale, y)
        dc.DrawLine(x, y-spot_scale, x, y+spot_scale)

  def draw_integration_results (self, dc) :
    scale = self.get_scale()
    if (scale >= 4) :
      bg_masks = self._img.get_drawable_background_mask()
      #dc.SetPen(wx.YELLOW_PEN)
      for (x, y) in bg_masks :
        dc.DrawCircle(x,y,1)
      int_masks = self._img.get_drawable_integration_mask()
      dc.SetPen(wx.CYAN_PEN)
      for (x, y) in int_masks :
        dc.DrawCircle(x,y,1)
    else :
      predictions = self._img.get_drawable_predictions()
      #dc.SetPen(wx.YELLOW_PEN)
      dc.SetBrush(wx.TRANSPARENT_BRUSH)
      for (x, y) in predictions :
        dc.DrawCircle(x, y, 5*scale)

  def OnSize (self, event) :
    if (self._img is not None) :
      w, h = self.GetSize()
      self._img.set_screen_size(w, h)

  def OnRecordMouse (self, event) :
    self.xmouse = event.GetX()
    self.ymouse = event.GetY()

  def OnMotion (self, event) :
    if (event.Dragging()) :
      self.was_dragged = True
      if (event.LeftIsDown()) :
        if (event.ShiftDown()) :
          self.OnMiddleDrag(event)
        else :
          self.OnLeftDrag(event)
      elif (event.MiddleIsDown()) :
        self.OnMiddleDrag(event)
      elif (event.RightIsDown()) :
        self.OnRightDrag(event)
    else :
      x, y = self._img.screen_coords_as_image_coords(event.GetX(),event.GetY())
      img_w, img_h = self._img.get_image_size()
      if (x < 0) or (x > img_w) or (y < 0) or (y > img_h) :
        self.GetParent().update_statusbar()
      else :
        info = self._img.get_point_info(x, y)
        self.GetParent().update_statusbar(info)

  def OnMiddleDown (self, event) :
    self.was_dragged = False
    self.OnRecordMouse(event)
    wx.SetCursor(wx.StockCursor(wx.CURSOR_HAND))

  def OnMiddleUp (self, event) :
    wx.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

  def OnLeftDown (self, event) :
    self.was_dragged = False
    if (event.ShiftDown()) :
      self.shift_was_down = True
      self.OnMiddleDown(event)
    else :
      self.line_end = None
      x, y = event.GetPositionTuple()
      self.line_start = self._img.screen_coords_as_image_coords(x, y)
      #self.OnRecordMouse(event)

  def OnDoubleClick (self, event) :
    pass

  def OnLeftDrag (self, event) :
    x, y = event.GetPositionTuple()
    self.line_end = self._img.screen_coords_as_image_coords(x, y)
    self.Refresh()

  def OnLeftUp (self, event) :
    if (self.shift_was_down) :
      pass
    elif (self.was_dragged) and (self.line_start is not None) :
      x, y = event.GetPositionTuple()
      self.line_end = self._img.screen_coords_as_image_coords(x, y)
      x1, y1 = self.line_start
      x2, y2 = self.line_end
      if (x1 <= x2) :
        line = self._img.line_between_points(x1, y1, x2, y2)
      else :
        line = self._img.line_between_points(x2, y2, x1, y1)
      self.GetParent().OnShowPlot(None)
      self.GetParent().plot_frame.show_plot(line)
    else :
      self.line = None
    self.Refresh()
    self.shift_was_down = False
    self.was_dragged = False

  def OnMiddleDrag (self, event) :
    self.OnTranslate(event)

  def OnRightDown (self, event) :
    self.was_dragged = False
    self.OnZoom(event)

  def OnRightDrag (self, event) :
    self.OnZoom(event)

  def OnZoom (self, event) :
    x, y = event.GetPositionTuple()
    img_x, img_y = self._img.screen_coords_as_image_coords(x, y)
    self.GetParent().OnShowZoom(None)
    self.GetParent().zoom_frame.set_zoom(img_x, img_y)

  def OnTranslate (self, event) :
    x, y = event.GetX(), event.GetY()
    delta_x = x - self.xmouse
    delta_y = y - self.ymouse
    self.OnRecordMouse(event)
    self.TranslateImage(delta_x, delta_y)

  def TranslateImage (self, delta_x, delta_y) :
    if (self.settings.zoom_level == 0) :
      return
    self._img.translate_image(delta_x, delta_y)
    self.Refresh()
    self.GetParent().settings_frame.refresh_thumbnail()

  def OnMouseWheel (self, event) :
    d_x = d_y = 0
    if (event.ShiftDown()) :
      d_x = - 10 * event.GetWheelRotation()
    else :
      d_y = - 10 * event.GetWheelRotation()
    self.TranslateImage(d_x, d_y)

  def OnEnter (self, event) :
    if (not event.MiddleIsDown()) and (not event.RightIsDown()) :
      wx.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

  def OnLeave (self, event) :
    self.was_dragged = False
    wx.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))

class ThumbnailView (XrayView) :
  def __init__ (self, *args, **kwds) :
    XrayView.__init__(self, *args, **kwds)
    self.flag_always_show_predictions = True

  def SetupEventHandlers (self) :
    self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)

  def set_image (self, image) :
    self._img = image
    self.SetSize(image.get_thumbnail_size())
    self.SetMinSize(image.get_thumbnail_size())
    self.GetParent().Layout()

  def OnPaint (self, event) :
    if (self._img is None) : return
    dc = wx.AutoBufferedPaintDCFactory(self)
    dc.SetBackground(wx.Brush((255,255,255)))
    dc.Clear()
    bitmap = self._img.get_thumbnail_bitmap()
    dc.SetBrush(wx.TRANSPARENT_BRUSH)
    dc.DrawBitmap(bitmap, 0, 0)
    x, y, w, h = self._img.get_thumbnail_box()
    dc.SetPen(wx.Pen('red', 2))
    dc.DrawRectangle(x, y, w-1, h-1)

  def OnLeftDown (self, event) :
    x, y = event.GetPositionTuple()
    self._img.center_view_from_thumbnail(x, y)
    self.Refresh()
    self.GetParent().refresh_main()

#class ZoomView (wx.Panel) :
class ZoomView (XrayView) :
  def __init__ (self, *args, **kwds) :
    super(ZoomView, self).__init__(*args, **kwds)
    self.SetSize((400,400))
    self.SetMinSize((400,400))
    self._img = None
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.x_center = None
    self.y_center = None
    self.zoom_level = 16
    self.screen = screen_params()
    self.screen.set_zoom(16)

  def SetupEventHandlers (self) :
    pass

  def set_zoom (self, x, y) :
    self.x_center = x
    self.y_center = y
    self.Refresh()

  def update_settings (self) :
    pass

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    if (not None in [self._img, self.x_center, self.y_center]) :
      wx_image = self._img.get_zoomed_bitmap(self.x_center, self.y_center)
      bitmap = wx_image.ConvertToBitmap()
      dc.DrawBitmap(bitmap, 0, 0)
    else :
      dc.SetPen(wx.Pen('red'))
      dc.DrawText("Right-click in the main image field to zoom.", 10, 10)
