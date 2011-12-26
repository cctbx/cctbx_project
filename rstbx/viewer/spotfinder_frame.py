import rstbx.viewer.display
import wx
from rstbx.viewer.frame import XrayFrame,SettingsFrame,SettingsPanel

def derived_draw_spotfinder_spots(self,dc):
  if self.settings.show_all_pix:
    for spot in self.spotfinder.images[self.frames[0]]["spots_total"]:
      for pxl in spot.bodypixels:
        x,y = self._img.image_coords_as_screen_coords(
          pxl.y,
          pxl.x)
        dc.SetPen(wx.Pen('green'))
        dc.SetBrush(wx.GREEN_BRUSH)
        dc.DrawCircle(x,y,1)

  if self.settings.show_max_pix:
    pink_brush = wx.Brush(colour=wx.Colour(255,192,203))
    for spot in self.spotfinder.images[self.frames[0]]["spots_total"]:
      x,y = self._img.image_coords_as_screen_coords(
        spot.max_pxl_y(),
        spot.max_pxl_x())
      dc.SetPen(wx.Pen('pink'))
      dc.SetBrush(pink_brush)
      dc.DrawCircle(x,y,1)

  if self.settings.show_ctr_mass:
    for spot in self.spotfinder.images[self.frames[0]]["spots_total"]:
      x,y = self._img.image_coords_as_screen_coords(
        spot.ctr_mass_y(),
        spot.ctr_mass_x())
      dc.SetPen(wx.Pen('red'))
      dc.SetBrush(wx.RED_BRUSH)
      dc.DrawCircle(x,y,1)
rstbx.viewer.display.XrayView.draw_spotfinder_spots = derived_draw_spotfinder_spots

class SpotFrame(XrayFrame) :
  def __init__ (self, *args, **kwds) :
    self.horizons_phil = kwds["horizons_phil"]
    self.spot_organizer = kwds["spot_organizer"]
    del kwds["horizons_phil"]; del kwds["spot_organizer"] #otherwise wx complains
    super(SpotFrame, self).__init__(*args, **kwds)
    self.viewer.horizons_phil = self.horizons_phil
    self.viewer.spotfinder = self.spot_organizer.S
    self.viewer.frames = self.spot_organizer.frames

  def OnShowSettings (self, event) :
    if (self.settings_frame is None) :
      frame_rect = self.GetRect()
      display_rect = wx.GetClientDisplayRect()
      x_start = frame_rect[0] + frame_rect[2]
      if (x_start > (display_rect[2] - 400)) :
        x_start = display_rect[2] - 400
      y_start = frame_rect[1]
      self.settings_frame = SpotSettingsFrame(self, -1, "Settings",
        style=wx.CAPTION|wx.MINIMIZE_BOX, pos=(x_start, y_start))
    self.settings_frame.Show()

class SpotSettingsFrame (SettingsFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SpotSettingsPanel(self, -1)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

class SpotSettingsPanel (SettingsPanel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.settings = self.GetParent().settings
    # CONTROLS 4: additional settings for derived class
    self.settings.show_ctr_mass = True
    self.settings.show_max_pix = True
    self.settings.show_all_pix = False
    self._sizer = wx.BoxSizer(wx.VERTICAL)
    s = self._sizer
    self.SetSizer(self._sizer)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt1 = wx.StaticText(self, -1, "Zoom level:")
    box.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # CONTROLS 1:  Lay them out on the display
    # Zoom control
    self.zoom_ctrl = wx.Choice(self, -1,
      choices=["Auto", "25%", "50%", "100%", "200%", "400%", "800%"])
    self.zoom_ctrl.SetSelection(self.settings.zoom_level)
    box.Add(self.zoom_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._sizer.Fit(self)
    box = wx.BoxSizer(wx.HORIZONTAL)
    s.Add(box)
    txt2 = wx.StaticText(self, -1, "Brightness")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Brightness control
    self.brightness_ctrl = wx.Slider(self, -1, size=(200,-1),
      style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    box.Add(self.brightness_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.brightness_ctrl.SetMin(10)
    self.brightness_ctrl.SetMax(500)
    self.brightness_ctrl.SetValue(self.settings.brightness)
    self.brightness_ctrl.SetTickFreq(25)

    # Center control
    self.center_ctrl = wx.CheckBox(self, -1, "Mark beam center")
    self.center_ctrl.SetValue(self.settings.show_beam_center)
    s.Add(self.center_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Center of mass control
    self.ctr_mass = wx.CheckBox(self, -1, "Spot centers of mass")
    self.ctr_mass.SetValue(self.settings.show_ctr_mass)
    s.Add(self.ctr_mass, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Max pixel control
    self.max_pix = wx.CheckBox(self, -1, "Spot max pixels")
    self.max_pix.SetValue(self.settings.show_max_pix)
    s.Add(self.max_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    # Spot pixels control
    self.all_pix = wx.CheckBox(self, -1, "Spot all pixels")
    self.all_pix.SetValue(self.settings.show_all_pix)
    s.Add(self.all_pix, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

    self.collect_values()

    # CONTROLS 3:  Bind events to actions
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.ctr_mass)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.max_pix)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdateCM, self.all_pix)

    txt3 = wx.StaticText(self, -1, "Thumbnail view:")
    s.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.thumb_panel = rstbx.viewer.display.ThumbnailView(
      parent=self,
      size=(256,256),
      style=wx.SUNKEN_BORDER)
    s.Add(self.thumb_panel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  # CONTROLS 2:  Fetch values from widgets
  def collect_values (self) :
    if self.settings.enable_collect_values:
      self.settings.zoom_level = self.zoom_ctrl.GetSelection()
      self.settings.brightness = self.brightness_ctrl.GetValue()
      self.settings.show_beam_center = self.center_ctrl.GetValue()
      self.settings.show_ctr_mass = self.ctr_mass.GetValue()
      self.settings.show_max_pix = self.max_pix.GetValue()
      self.settings.show_all_pix = self.all_pix.GetValue()

  def OnUpdateCM (self, event) :
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)

