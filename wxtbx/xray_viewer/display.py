
from wxtbx import icons
import wxtbx.xray_viewer
import wx
import os

class XrayView (wx.ScrolledWindow) :
  def __init__ (self, *args, **kwds) :
    self._img = None
    super(XrayView, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.Bind(wx.EVT_MOTION, self.OnMove)

  def DoGetBestSize (self) :
    if (self._img is None) :
      return (640, 640)
    else :
      x, y = self._img.get_size()
      return (x+20, y+20)

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    x_offset, y_offset = self.GetViewStart()
    bitmap = self._img.get_bitmap()
    dc.DrawBitmap(bitmap, 10 - x_offset, 10 - y_offset)
    center_x, center_y = self._img.get_beam_center()
    dc.SetPen(wx.Pen('red'))
    x0 = center_x - x_offset + 10
    y0 = center_y - y_offset + 10
    dc.DrawLine(x0 - 10, y0, x0 + 10, y0)
    dc.DrawLine(x0, y0 - 10, x0, y0 + 10)

  def OnSize (self, event) :
    w, h = self.DoGetBestSize()
    self.SetVirtualSize((w,h))
    self.SetScrollbars(1, 1, w, h)

  def OnMove (self, event) :
    x, y = event.GetPositionTuple()
    x_offset, y_offset = self.GetViewStart()
    img_x = x + x_offset - 10
    img_y = y + y_offset - 10
    img_w, img_h = self._img.get_size()
    if (img_x < 0) or (img_x > img_w) or (img_y < 0) or (img_y > img_h) :
      self.GetParent().display_resolution(None)
    else :
      center_x, center_y = self._img.get_beam_center()
      d_min = self._img.get_d_min_at_point(img_x, img_y)
      self.GetParent().display_resolution(d_min)

  def load_image (self, file_name) :
    wx.GetApp().Yield(True)
    self._img = wxtbx.xray_viewer.image(file_name, self.settings)
    self.OnSize(None)
    self.Refresh()

  def update_settings (self) :
    self._img.update_settings()
    self.OnSize(None)
    self.Refresh()

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
    self.toolbar.Realize()
    self.Fit()

  def load_image (self, file_name) :
    file_name = os.path.abspath(file_name)
    self._panel.load_image(file_name)
    self.SetTitle(file_name)
    self.display_resolution(None)
    self.Layout()

  def display_resolution (self, d_min) :
    if (d_min is None) :
      self.statusbar.SetStatusText("Move mouse over image to show resolution")
    else :
      self.statusbar.SetStatusText("d_min = %.2f A" % d_min)

  def update_settings (self) :
    self._panel.update_settings()

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

class SettingsFrame (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SettingsPanel(self, -1)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def OnDestroy (self, event) :
    self.GetParent().settings_frame = None

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
      choices=["100%", "50%", "25%"])
    self.zoom_ctrl.SetSelection(self.settings.zoom_level)
    box.Add(self.zoom_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._sizer.Fit(self)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt2 = wx.StaticText(self, -1, "Brightness")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl = wx.Slider(self, -1, size=(160,-1))
    box.Add(self.brightness_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl.SetMin(0)
    self.brightness_ctrl.SetMax(200)
    self.brightness_ctrl.SetValue(self.settings.brightness)
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdate, self.brightness_ctrl)

  def OnUpdate (self, event) :
    zooms = ["100%", "50%"] #, None, "25%"]
    self.settings.zoom_level = zooms.index(self.zoom_ctrl.GetStringSelection())
    self.settings.brightness = self.brightness_ctrl.GetValue()
    self.GetParent().GetParent().update_settings()
