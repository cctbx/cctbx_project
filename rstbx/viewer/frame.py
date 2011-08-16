
import rstbx.viewer.display
import wxtbx.plots
from wxtbx import bitmaps
from wxtbx import icons
import wx
import os

class XrayFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    super(XrayFrame, self).__init__(*args, **kwds)
    self.settings = rstbx.viewer.settings()
    self.viewer = rstbx.viewer.display.XrayView(self, -1, size=(1024,640))
    self.viewer.SetMinSize((640,640))
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.sizer.Add(self.viewer, 1, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.settings_frame = None
    self.zoom_frame = None
    self.plot_frame = None
    self._img = None
    self._distl = None
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Load file",
      bitmap=icons.hkl_file.GetBitmap(),
      shortHelp="Load file",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnLoadFile, btn)
    #btn = self.toolbar.AddLabelTool(id=-1,
    #  label="Settings",
    #  bitmap=icons.advancedsettings.GetBitmap(),
    #  shortHelp="Settings",
    #  kind=wx.ITEM_NORMAL)
    #self.Bind(wx.EVT_MENU, self.OnShowSettings, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Zoom",
      bitmap=icons.search.GetBitmap(),
      shortHelp="Zoom",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnZoom, btn)
    txt = wx.StaticText(self.toolbar, -1, "Image:")
    self.toolbar.AddControl(txt)
    self.image_chooser = wx.Choice(self.toolbar, -1, size=(300,-1))
    self.toolbar.AddControl(self.image_chooser)
    self.Bind(wx.EVT_CHOICE, self.OnChooseImage, self.image_chooser)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Previous",
      bitmap=bitmaps.fetch_icon_bitmap("actions","1leftarrow"),
      shortHelp="Previous",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnPrevious, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Next",
      bitmap=bitmaps.fetch_icon_bitmap("actions","1rightarrow"),
      shortHelp="Next",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnNext, btn)
    self.toolbar.Realize()
    self.Fit()
    self.SetMinSize(self.GetSize())
    self.OnShowSettings(None)

  def load_image (self, file_name) :
    file_name = os.path.abspath(file_name)
    self._img = rstbx.viewer.image(file_name)
    self.viewer.set_image(self._img)
    self.settings_frame.set_image(self._img)
    self.SetTitle(file_name)
    items = self.image_chooser.GetItems()
    if (not file_name in items) :
      items.append(file_name)
    self.image_chooser.SetItems(items)
    self.image_chooser.SetStringSelection(file_name)
    self.update_statusbar()
    self.Layout()

  def load_distl_output (self, file_name) :
    from libtbx.easy_pickle import load
    distl = load(file_name)
    self._distl = distl
    img_files = []
    for img_id in sorted(distl.images.keys()) :
      img = distl.images[img_id]
      img_files.append(img['relpath'])
    if (len(img_files) == 0) :
      raise Sorry("No images in this result!")
    self.image_chooser.SetItems([ os.path.basename(f) for f in img_files ])
    self.image_chooser.SetSelection(0)
    self.load_image(img_files[0])
    self.annotate_image(img_files[0])

  def annotate_image (self, file_name) :
    assert (self._distl is not None)
    for img_id in sorted(self._distl.images.keys()) :
      img = self._distl.images[img_id]
      if (img['relpath'] == file_name) :
        spots = img['spotoutput']['inlier_spots']
        self._img.set_spots(spots)
        break

  def update_statusbar (self, info=None) :
    if (info is None) :
      self.statusbar.SetStatusText("Click and drag to plot intensity profile; "+
        "middle-click to pan, right-click to zoom")
    else :
      self.statusbar.SetStatusText(info.format())

  def update_settings (self, layout=True) :
    self.viewer.update_settings(layout)

  def set_brightness (self, brightness) :
    if (brightness > 0) and (brightness <= 500) :
      self.settings.brightness = brightness
      if (self.settings_frame is not None) :
        self.settings_frame.update_controls()
      self.viewer.update_settings(layout=False)

  def OnLoadFile (self, event) :
    file_name = wx.FileSelector("Reflections file",
      wildcard="Image files (*.img, *.ccd, *.mccd)|*.img;*.ccd;*.mccd",
      default_path="",
      flags=wx.OPEN)
    if (file_name != "") :
      self.load_image(file_name)

  def OnLoadLabelitResult (self, event) :
    file_name = wx.FileSelector("Labelit result",
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
      y_start = frame_rect[1]
      self.settings_frame = SettingsFrame(self, -1, "Settings",
        style=wx.CAPTION|wx.MINIMIZE_BOX, pos=(x_start, y_start))
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
    if (self.settings.zoom_level == 6) :
      self.settings.zoom_level = 0
    else :
      self.settings.zoom_level += 1
    self.viewer.update_settings(layout=True)
    if (self.settings_frame is not None) :
      self.settings_frame.update_controls()

  def OnChooseImage (self, event) :
    print "Not implemented"

  def OnPrevious (self, event) :
    print "Not implemented"

  def OnNext (self, event) :
    print "Not implemented"

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
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def OnDestroy (self, event) :
    self.GetParent().settings_frame = None

  def update_controls (self) :
    self.panel.zoom_ctrl.SetSelection(self.settings.zoom_level)
    self.panel.brightness_ctrl.SetValue(self.settings.brightness)

  def set_image (self, image) :
    self.panel.thumb_panel.set_image(image)
    self.panel.GetSizer().Layout()
    self.sizer.Fit(self.panel)
    self.Layout()
    self.Fit()

  def refresh_thumbnail (self) :
    self.panel.thumb_panel.Refresh()

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
    self.spots_ctrl = wx.CheckBox(self, -1, "Show predictions")
    self.spots_ctrl.SetValue(self.settings.show_predictions)
    s.Add(self.spots_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#    self.invert_ctrl = wx.CheckBox(self, -1, "Invert beam center axes")
#    self.invert_ctrl.SetValue(self.settings.invert_beam_center_axes)
#    s.Add(self.invert_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.zoom_ctrl)
    self.Bind(wx.EVT_SLIDER, self.OnUpdateBrightness, self.brightness_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.center_ctrl)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.spots_ctrl)
    txt3 = wx.StaticText(self, -1, "Thumbnail view:")
    s.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.thumb_panel = rstbx.viewer.display.ThumbnailView(
      parent=self,
      size=(256,256),
      style=wx.SUNKEN_BORDER)
    s.Add(self.thumb_panel, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#    self.Bind(wx.EVT_CHECKBOX, self.OnUpdate2, self.invert_ctrl)

  def collect_values (self) :
    self.settings.zoom_level = self.zoom_ctrl.GetSelection()
    self.settings.brightness = self.brightness_ctrl.GetValue()
    self.settings.show_beam_center = self.center_ctrl.GetValue()
    self.settings.show_predictions = self.spots_ctrl.GetValue()
#    self.settings.invert_beam_center_axes = self.invert_ctrl.GetValue()

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

  def refresh_main (self) :
    self.GetParent().GetParent().viewer.Refresh()

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
      wx_image = self._img.get_zoomed_bitmap(self.x_center, self.y_center)
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
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def __getattr__ (self, name) :
    return getattr(self.plot, name)

  def OnDestroy (self, event) :
    self.GetParent().plot_frame = None

class LinePlot (wxtbx.plots.plot_container) :
  def show_plot (self, y_data) :
    self.figure.clear()
    ax = self.figure.add_subplot(111)
    x_data = range(len(y_data))
    ax.plot(x_data, y_data, 'b-', linewidth=1)
    ax.set_ylabel("Intensity")
    self.canvas.draw()
