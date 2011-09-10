import rstbx.viewer.display
import wx
import os
from rstbx.viewer.frame import XrayFrame,SettingsFrame,SettingsPanel

class SBFrame(XrayFrame) :
  def __init__ (self, *args, **kwds) :
    self.horizons_phil = kwds["horizons_phil"]
    del kwds["horizons_phil"]
    super(SBFrame, self).__init__(*args, **kwds)

  def OnShowSettings (self, event) :
    if (self.settings_frame is None) :
      frame_rect = self.GetRect()
      display_rect = wx.GetClientDisplayRect()
      x_start = frame_rect[0] + frame_rect[2]
      if (x_start > (display_rect[2] - 400)) :
        x_start = display_rect[2] - 400
      y_start = frame_rect[1]
      self.settings_frame = SBSettingsFrame(self, -1, "Settings",
        style=wx.CAPTION|wx.MINIMIZE_BOX, pos=(x_start, y_start))
    self.settings_frame.Show()

  def load_image_restricted (self) :
    #This method can probably be back-edited to rely more on the base class code...later

    file_name = os.path.abspath(self.path)
    x,y = self.viewer._img.last_thumb_x, self.viewer._img.last_thumb_y
    self._img = rstbx.viewer.image(self.path)
    #self.viewer.set_image(self._img)
    self.viewer._img = self._img
    self.viewer._img.set_screen_size(*(self.viewer.GetSize()))
    self.viewer.line = None
    scales = [0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
    zoom = scales[self.viewer.settings.zoom_level]
    self.viewer._img.set_zoom(zoom)
    self.viewer._img.update_settings(
      brightness=self.viewer.settings.brightness)
    self.viewer._img.center_view_from_thumbnail(x,y)
    self.viewer.Refresh()
    if (self.viewer.GetParent().zoom_frame is not None) :
      self.viewer.GetParent().zoom_frame.Refresh()
    self.viewer.GetParent().settings_frame.refresh_thumbnail()

    self.settings_frame.set_image(self._img)
    self.SetTitle(file_name)
    items = self.image_chooser.GetItems()
    if (not file_name in items) :
      items.append(file_name)
    self.image_chooser.SetItems(items)
    self.image_chooser.SetStringSelection(file_name)
    self.update_statusbar()
    self.Layout()

class SBSettingsFrame (SettingsFrame) :
  def __init__ (self, *args, **kwds) :
    super(SettingsFrame, self).__init__(*args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = SBSettingsPanel(self, -1)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

class SBSettingsPanel (SettingsPanel) :
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
      choices=["Auto", "25%", "50%", "100%", "200%", "400%", "800%"])
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

    box = wx.BoxSizer(wx.HORIZONTAL)
    from wxtbx.phil_controls.floatctrl import FloatCtrl
    from wxtbx.phil_controls import EVT_PHIL_CONTROL
    self.distance_ctrl = FloatCtrl(self, -1, pos=(300,180), size=(80,-1),
    value=80.00,
    name="Detector Distance")
    self.distance_ctrl.SetMax(1000)
    self.distance_ctrl.SetMin(5)
    self.distance_ctrl.SetOptional(False)
    box.Add(self.distance_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txtd = wx.StaticText(self, -1,  "Detector Distance",)
    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    s.Add(box)

    from wxtbx.phil_controls.intctrl import IntCtrl
    for s_quad in ["UL","UR","LL","LR"]:

      box = wx.BoxSizer(wx.HORIZONTAL)
      setattr(self,s_quad+"x_ctrl",IntCtrl(self, -1, pos=(300,180), size=(80,-1),value=0,name=s_quad+"x"))
      spinbtn = wx.SpinButton(self, -1)
      getattr(self,s_quad+"x_ctrl").AttachSpinner(spinbtn)
      box.Add(getattr(self,s_quad+"x_ctrl"), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box.Add(spinbtn, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
      txtd = wx.StaticText(self, -1,  s_quad+" x",)
      box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

      setattr(self,s_quad+"y_ctrl",IntCtrl(self, -1, pos=(300,180), size=(80,-1),value=0,name=s_quad+"y"))
      spinbtn = wx.SpinButton(self, -1)
      getattr(self,s_quad+"y_ctrl").AttachSpinner(spinbtn)
      box.Add(getattr(self,s_quad+"y_ctrl"), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box.Add(spinbtn, 0, wx.RIGHT|wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
      txtd = wx.StaticText(self, -1,  s_quad+" y",)
      box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      s.Add(box)

    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(EVT_PHIL_CONTROL, self.OnUpdateDist, self.distance_ctrl)
    for s_quad in ["UL","UR","LL","LR"]:
      self.Bind(EVT_PHIL_CONTROL, self.OnUpdateQuad, getattr(self,s_quad+"x_ctrl"))
      self.Bind(EVT_PHIL_CONTROL, self.OnUpdateQuad, getattr(self,s_quad+"y_ctrl"))
    txt3 = wx.StaticText(self, -1, "Thumbnail view:")
    s.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.thumb_panel = rstbx.viewer.display.ThumbnailView(
      parent=self,
      size=(256,256),
      style=wx.SUNKEN_BORDER)
    s.Add(self.thumb_panel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def collect_values (self) :
    if self.settings.enable_collect_values:
      self.settings.zoom_level = self.zoom_ctrl.GetSelection()
      self.settings.brightness = self.brightness_ctrl.GetValue()
      self.settings.show_beam_center = self.center_ctrl.GetValue()
      self.settings.distance = self.distance_ctrl.GetPhilValue()

  def OnUpdateQuad (self, event) :
    for iq,s_quad in enumerate(["UL","UR","LL","LR"]):
      setattr(self.settings,s_quad+"x",getattr(self,s_quad+"x_ctrl").GetPhilValue())
      setattr(self.settings,s_quad+"y",getattr(self,s_quad+"y_ctrl").GetPhilValue())
      self.GetParent().GetParent().horizons_phil.distl.quad_translations[2*iq+0]=getattr(self.settings,s_quad+"x")
      self.GetParent().GetParent().horizons_phil.distl.quad_translations[2*iq+1]=getattr(self.settings,s_quad+"y")
    self.GetParent().GetParent().load_image_restricted()
    #self.GetParent().GetParent().update_settings(layout=True)

  def OnUpdateDist (self, event) :
    print "update dist"
    self.collect_values()
    self.GetParent().GetParent().update_settings(layout=False)
