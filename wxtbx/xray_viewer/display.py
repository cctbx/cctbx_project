
import wxtbx.xray_viewer
import wxtbx.plots
from wxtbx import icons
import wx
import os

class XrayView (wx.ScrolledWindow) :
  def __init__ (self, *args, **kwds) :
    self._img = None
    super(XrayView, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    self.Bind(wx.EVT_PAINT, self.OnPaint)
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
    self.SetMinSize((640,640))
    self.x_offset = None
    self.y_offset = None
    self.xmouse = None
    self.ymouse = None
    self.line_start = None
    self.line_end = None
    self.was_dragged = False
    self.was_double_clicked = False

  def set_image (self, image) :
    self._img = image
    self.update_settings()

  def update_settings (self, layout=True) :
    self.line = None
    self._img.update_settings(
      brightness=self.settings.brightness,
      zoom_level=self.settings.zoom_level,
      invert_beam_center=self.settings.invert_beam_center_axes,
      w=self.GetSize()[0] - 20,
      h=self.GetSize()[1] - 20)
    if (layout) :
      self.OnSize(None)
    self.Refresh()
    if (self.GetParent().zoom_frame is not None) :
      self.GetParent().zoom_frame.Refresh()

  def image_coords_as_screen_coords (self, x, y) :
    x_offset, y_offset = self.GetViewStart()
    return (x - x_offset + 10, y - y_offset + 10)

  def screen_coords_as_image_coords (self, x=None, y=None, event=None) :
    assert (event is not None) or (not None in [x,y])
    if (event is not None) :
      assert ([x,y] == [None, None])
      x, y = event.GetPositionTuple()
    x_offset, y_offset = self.GetViewStart()
    return (x + x_offset - 10, y + y_offset - 10)

  def DoGetBestSize (self) :
    if (self._img is None) :
      return (640, 640)
    else :
      x, y = self._img.get_size()
      return (x+20, y+20)

  # EVENTS
  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    x_offset, y_offset = self.GetViewStart()
    bitmap = self._img.get_bitmap()
    dc.DrawBitmap(bitmap, 10 - x_offset, 10 - y_offset)
    if (self.settings.show_beam_center) :
      center_x, center_y = self._img.get_beam_center()
      dc.SetPen(wx.Pen('red'))
      x0, y0 = self.image_coords_as_screen_coords(center_x, center_y)
    #  print center_x, center_y, x0, y0
      dc.DrawLine(x0 - 10, y0, x0 + 10, y0)
      dc.DrawLine(x0, y0 - 10, x0, y0 + 10)
    if (self.line_start is not None) and (self.line_end is not None) :
      dc.SetPen(wx.Pen('red'))
      x1, y1 = self.image_coords_as_screen_coords(*(self.line_start))
      x2, y2 = self.image_coords_as_screen_coords(*(self.line_end))
      dc.DrawLine(x1, y1, x2, y2)

  def OnSize (self, event) :
    w, h = self.DoGetBestSize()
    self.SetVirtualSize((w,h))
    self.SetScrollbars(1, 1, w, h)
    mouse = wx.GetMouseState()

  def OnRecordMouse (self, event) :
    self.x_offset, self.y_offset = self.GetViewStart()
    self.xmouse = event.GetX()
    self.ymouse = event.GetY()

  def OnMotion (self, event) :
    if (event.Dragging()) :
      self.was_dragged = True
      if (event.LeftIsDown()) :
        self.OnLeftDrag(event)
      elif (event.MiddleIsDown()) :
        self.OnMiddleDrag(event)
      elif (event.RightIsDown()) :
        self.OnRightDrag(event)
    elif (self.was_double_clicked) :
      x, y = event.GetPositionTuple()
      self.Refresh()
    else :
      x, y = self.screen_coords_as_image_coords(event=event)
      img_w, img_h = self._img.get_size()
      if (x < 0) or (x > img_w) or (y < 0) or (y > img_h) :
        self.GetParent().update_statusbar()
      else :
        info = self._img.get_point_info(x, y)
        self.GetParent().update_statusbar(info)

  def OnMiddleDown (self, event) :
    self.was_dragged = False
    self.was_double_clicked = False
    self.OnRecordMouse(event)
    wx.SetCursor(wx.StockCursor(wx.CURSOR_HAND))

  def OnMiddleUp (self, event) :
    wx.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

  def OnLeftDown (self, event) :
    self.was_dragged = False
    self.was_double_clicked = False
    self.line_end = None
    self.line_start = self.screen_coords_as_image_coords(event=event)
    #self.OnRecordMouse(event)

  def OnDoubleClick (self, event) :
    self.was_dragged = False
    self.was_double_clicked = True
    self.OnRecordMouse(event)

  def OnLeftDrag (self, event) :
    x, y = event.GetPositionTuple()
    self.line_end = self.screen_coords_as_image_coords(event=event)
    self.Refresh()

  def OnLeftUp (self, event) :
    if (self.was_dragged) and (self.line_start is not None) :
      self.line_end = self.screen_coords_as_image_coords(event=event)
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
    self.was_dragged = False
    #elif (self.was_double_clicked) and (self.line is not None) :
    #  x1, y1, x2, y2 = self.line
    #  print self._img.line_between_points(x1, y1, x2, y2,
    #    zoom=self.settings.zoom_level)

  def OnMiddleDrag (self, event) :
    self.OnTranslate(event)

  def OnRightDown (self, event) :
    self.was_dragged = False
    self.was_double_clicked = False
    self.OnZoom(event)

  def OnRightDrag (self, event) :
    self.OnZoom(event)

  def OnZoom (self, event) :
    x, y = event.GetPositionTuple()
    x_offset, y_offset = self.GetViewStart()
    img_x = x + x_offset - 10
    img_y = y + y_offset - 10
    self.GetParent().OnShowZoom(None)
    self.GetParent().zoom_frame.set_zoom(img_x, img_y)

  def OnTranslate (self, event) :
    x, y = event.GetX(), event.GetY()
    delta_x = x - self.xmouse
    delta_y = y - self.ymouse
    self.Scroll(self.x_offset - delta_x, self.y_offset - delta_y)

  def OnMouseWheel (self, event) :
    d_brightness = 10 * event.GetWheelRotation()
    self.GetParent().set_brightness(self.settings.brightness + d_brightness)

  def OnEnter (self, event) :
    if (not event.MiddleIsDown()) and (not event.RightIsDown()) :
      wx.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))

  def OnLeave (self, event) :
    self.was_dragged = False
    self.was_double_clicked = False
    wx.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))

########################################################################
class XrayFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    super(XrayFrame, self).__init__(*args, **kwds)
    self.settings = wxtbx.xray_viewer.settings()
    self._panel = XrayView(self, -1, size=(640,640))
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.sizer.Add(self._panel, 1, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.settings_frame = None
    self.zoom_frame = None
    self.plot_frame = None
    self._img = None
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Load file",
      bitmap=icons.hkl_file.GetBitmap(),
      shortHelp="Load file",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnLoadFile, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Settings",
      bitmap=icons.advancedsettings.GetBitmap(),
      shortHelp="Settings",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnShowSettings, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Zoom",
      bitmap=icons.search.GetBitmap(),
      shortHelp="Zoom",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnZoom, btn)
    self.toolbar.Realize()
    self.Fit()
    self.SetMinSize(self.GetSize())

  def load_image (self, file_name) :
    file_name = os.path.abspath(file_name)
    self._img = wxtbx.xray_viewer.image(file_name)
    self._panel.set_image(self._img)
    self.SetTitle(file_name)
    self.update_statusbar()
    self.Layout()

  def update_statusbar (self, info=None) :
    if (info is None) :
      self.statusbar.SetStatusText("Click and drag to plot intensity profile; "+
        "middle-click to pan, right-click to zoom")
    else :
      self.statusbar.SetStatusText(info.format())

  def update_settings (self, layout=True) :
    self._panel.update_settings(layout)

  def set_brightness (self, brightness) :
    if (brightness > 0) and (brightness <= 500) :
      self.settings.brightness = brightness
      if (self.settings_frame is not None) :
        self.settings_frame.update_controls()
      self._panel.update_settings(layout=False)

  def OnLoadFile (self, event) :
    file_name = wx.FileSelector("Reflections file",
      wildcard="Image files (*.img, *.ccd, *.mccd)|*.img;*.ccd;*.mccd",
      default_path="",
      flags=wx.OPEN)
    if (file_name != "") :
      self.load_image(file_name)

  def OnShowSettings (self, event) :
    if (self.settings_frame is None) :
      frame_rect = self.GetRect()
      display_rect = wx.GetClientDisplayRect()
      x_start = frame_rect[0] + frame_rect[2]
      if (x_start > (display_rect[2] - 400)) :
        x_start = display_rect[2] - 400
      y_start = frame_rect[1] + 200
      if (y_start > (display_rect[3] - 200)) :
        y_start = display_rect[3] - 200
      self.settings_frame = SettingsFrame(self, -1, "Settings",
        style=wx.CAPTION|wx.CLOSE_BOX, pos=(x_start, y_start))
    self.settings_frame.Show()

  def OnShowZoom (self, event) :
    if (self.zoom_frame is None) :
      self.zoom_frame = ZoomFrame(self, -1, "Zoom",
        style=wx.CAPTION|wx.CLOSE_BOX)
      self.zoom_frame.set_image(self._img)
      self.zoom_frame.Show()

  def OnShowPlot (self, event) :
    if (self.plot_frame is None) :
      self.plot_frame = PlotFrame(self, -1, "Intensity profile",
        style=wx.CAPTION|wx.CLOSE_BOX)
      self.plot_frame.Show()

  def OnZoom (self, event) :
    if (self.settings.zoom_level == 1) :
      self.settings.zoom_level = 0
    else :
      self.settings.zoom_level = 1
    self._panel.update_settings(layout=True)
    if (self.settings_frame is not None) :
      self.settings_frame.update_controls()

class SettingsFrame (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SettingsPanel(self, -1)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def OnDestroy (self, event) :
    self.GetParent().settings_frame = None

  def update_controls (self) :
    self.panel.zoom_ctrl.SetSelection(self.settings.zoom_level)
    self.panel.brightness_ctrl.SetValue(self.settings.brightness)

class SettingsPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.settings = self.GetParent().settings
    self._sizer = wx.BoxSizer(wx.VERTICAL)
    s = self._sizer
    self.SetSizer(self._sizer)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt1 = wx.StaticText(self, -1, "Zoom level:")
    box.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.zoom_ctrl = wx.Choice(self, -1,
      choices=["Auto", "100%", "50%", "25%"])
    self.zoom_ctrl.SetSelection(self.settings.zoom_level)
    box.Add(self.zoom_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._sizer.Fit(self)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt2 = wx.StaticText(self, -1, "Brightness")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl = wx.Slider(self, -1, size=(200,-1),
      style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    box.Add(self.brightness_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl.SetMin(10)
    self.brightness_ctrl.SetMax(500)
    self.brightness_ctrl.SetValue(self.settings.brightness)
    self.brightness_ctrl.SetTickFreq(25)
    self.center_ctrl = wx.CheckBox(self, -1, "Mark beam center")
    self.center_ctrl.SetValue(self.settings.show_beam_center)
    s.Add(self.center_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.invert_ctrl = wx.CheckBox(self, -1, "Invert beam center axes")
    self.invert_ctrl.SetValue(self.settings.invert_beam_center_axes)
    s.Add(self.invert_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.invert_ctrl)

  def collect_values (self) :
    self.settings.zoom_level = self.zoom_ctrl.GetSelection()
    self.settings.brightness = self.brightness_ctrl.GetValue()
    self.settings.show_beam_center = self.center_ctrl.GetValue()
    self.settings.invert_beam_center_axes = self.invert_ctrl.GetValue()

  def OnUpdate (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=True)

  def OnUpdateBrightness (self, event) :
    mouse = wx.GetMouseState()
    if (mouse.LeftDown()) : return
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)

  def OnUpdate2 (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)

class ZoomPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    super(ZoomPanel, self).__init__(*args, **kwds)
    self.SetSize((400,400))
    self._img = None
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.x_center = None
    self.y_center = None

  def set_image (self, image) :
    self._img = image

  def set_zoom (self, x, y) :
    self.x_center = x
    self.y_center = y
    self.Refresh()

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    if (not None in [self._img, self.x_center, self.y_center]) :
      wx_image = self._img.get_zoomed_region(self.x_center, self.y_center,
        zoom=self.GetParent().settings.zoom_level)
      bitmap = wx_image.ConvertToBitmap()
      dc.DrawBitmap(bitmap, 0, 0)
    else :
      dc.SetPen(wx.Pen('red'))
      dc.DrawText("Right-click in the main image field to zoom.", 10, 10)

class ZoomFrame (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(ZoomFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    self.panel = ZoomPanel(self, -1)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr.Add(self.panel, 1, wx.EXPAND)
    szr.Fit(self.panel)
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def __getattr__ (self, name) :
    return getattr(self.panel, name)

  def OnDestroy (self, event) :
    self.GetParent().zoom_frame = None

class PlotFrame (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(PlotFrame, self).__init__(*args, **kwds)
    self.plot = LinePlot(self, figure_size=(8,3))
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr.Add(self.plot, 1, wx.EXPAND)
    self.Fit()

  def __getattr__ (self, name) :
    return getattr(self.plot, name)

class LinePlot (wxtbx.plots.plot_container) :
  def show_plot (self, y_data) :
    self.figure.clear()
    ax = self.figure.add_subplot(111)
    x_data = range(len(y_data))
    ax.plot(x_data, y_data, 'b-', linewidth=1)
    ax.set_ylabel("Intensity")
    self.canvas.draw()
