# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# Known issues: Recentering on resize and when switching between
# different image types.  Ring centre on image switch.
#
# $Id: frame.py 323 2012-05-16 21:40:27Z hattne $

import os
import wx

from rstbx.viewer.frame import EVT_EXTERNAL_UPDATE
from rstbx.viewer.frame import XrayFrame as XFBaseClass
from rstbx.viewer import settings as rv_settings, image as rv_image

from rstbx.slip_viewer.slip_display import AppFrame
class XrayFrame (AppFrame,XFBaseClass) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self,*args,**kwds)
    self.settings = rv_settings()

    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)

    self.init_pyslip_presizer()

    self.sizer.Add(self.viewer, 1, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.settings_frame = None
    self._calibration_frame = None
    self._ring_frame = None
    self.zoom_frame = None
    self.plot_frame = None

    self.metrology_matrices = None
    self.params = None

    # Currently displayed image.  XXX Can this be zapped?
    self._img = None

    self._distl = None
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.setup_toolbar()
    self.toolbar.Realize()
    self.mb = wx.MenuBar()
    self.setup_menus()
    self.SetMenuBar(self.mb)
    self.Fit()
    self.SetMinSize(self.GetSize())
    self.SetSize((720,720))
    self.OnShowSettings(None)
    self.Bind(EVT_EXTERNAL_UPDATE, self.OnExternalUpdate)

    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=wx.ID_BACKWARD)
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=wx.ID_FORWARD)
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=self._id_calibration)
    self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=self._id_ring)


  def init_pyslip_presizer(self):
    self.viewer  = wx.Panel(self, wx.ID_ANY, size=(1024,640))
    self.viewer.SetMinSize((640,640))
    self.viewer.SetBackgroundColour(wx.WHITE)
    self.viewer.ClearBackground()
    # create select event dispatch directory
    self.demo_select_dispatch = {}

    #self.tile_directory = None#"/Users/nksauter/rawdata/demo/insulin_1_001.img"

    # build the GUI
    self.make_gui(self.viewer)

    from rstbx.slip_viewer import pyslip
    # finally, bind events to handlers
    self.pyslip.Bind(pyslip.EVT_PYSLIP_SELECT, self.handle_select_event)
    self.pyslip.Bind(pyslip.EVT_PYSLIP_POSITION, self.handle_position_event)
    self.pyslip.Bind(pyslip.EVT_PYSLIP_LEVEL, self.handle_level_change)

  def init_pyslip_postsizer(self):
    self.pyslip.ZoomToLevel(-2)#tiles.zoom_level
    self.pyslip.GotoPosition(
      self.pyslip.tiles.get_initial_instrument_centering_within_picture_as_lon_lat()
    )

  def setup_menus (self) :
    file_menu = wx.Menu()
    self.mb.Append(file_menu, "File")
    item = file_menu.Append(-1, "Open integration results...")
    self.Bind(wx.EVT_MENU, self.OnLoadIntegration, item)
    item = file_menu.Append(-1, "Open image...")
    self.Bind(wx.EVT_MENU, self.OnLoadFile, item)
    actions_menu = wx.Menu()
    self.mb.Append(actions_menu, "Actions")
    item = actions_menu.Append(-1, "Change beam center...")
    self.Bind(wx.EVT_MENU, self.OnChangeBeamCenter, item)
    item = actions_menu.Append(-1, "Reset beam center to header value")
    self.Bind(wx.EVT_MENU, lambda evt: self.viewer.ResetBeamCenter(), item)
    item = actions_menu.Append(-1, "Save screenshot...")
    self.Bind(wx.EVT_MENU, self.OnScreenShot, item)

    # Known wxWidgets/wxPython issue
    # (http://trac.wxwidgets.org/ticket/12394): stock item ID is
    # expected for zero-length text.  Work around by making text
    # contain single space. XXX Placement
    self._id_calibration = wx.NewId()
    item = actions_menu.Append(self._id_calibration, " ")
    self.Bind(wx.EVT_MENU, self.OnCalibration, source=item)

    # XXX Placement
    self._id_ring = wx.NewId()
    item = actions_menu.Append(self._id_ring, " ")
    self.Bind(wx.EVT_MENU, self.OnRing, source=item)


  def load_image (self, file_name_or_data) :
    """The load_image() function displays the image from @p
    file_name_or_data.  The chooser is updated appropriately.
    """

    if (isinstance(file_name_or_data, dict)):
      img = rv_image(file_name_or_data)
      self.SetTitle(
        file_name_or_data.get("TIMESTAMP", "No timestamp available"))
    elif (isinstance(file_name_or_data, str)):
      img = rv_image(os.path.abspath(file_name_or_data))
      self.SetTitle(file_name_or_data)
    else:
      img = rv_image(os.path.abspath(file_name_or_data.encode("ascii")))
      self.SetTitle(file_name_or_data)

    # Update the selection in the chooser.
    i = self.add_file_name_or_data(file_name_or_data)
    self.image_chooser.SetSelection(i)

    self.pyslip.tiles.set_image(
      file_name_or_data=img, metrology_matrices=self.metrology_matrices)

    # Initialise position zoom level for first image.  XXX Why do we
    # have to coll ZoomToLevel to refresh subsequent images?
    if (self._img is None):
      self.init_pyslip_postsizer()
    else:
      self.pyslip.ZoomToLevel(self.pyslip.tiles.zoom_level)

    self._img = img # XXX

    self.settings_frame.set_image(self._img)
    self.update_statusbar() # XXX Not always working?
    self.Layout()

    beam_pixel_fast,beam_pixel_slow = self.pyslip.tiles.raw_image.get_beam_center_pixels_fast_slow()

    self.beam_center_cross_data = [
      ((self.pyslip.tiles.picture_fast_slow_to_map_relative(beam_pixel_fast+3.,beam_pixel_slow),
        self.pyslip.tiles.picture_fast_slow_to_map_relative(beam_pixel_fast-3.,beam_pixel_slow)),
        {'width': 2, 'color': '#0000FFA0', 'closed': False}),
      ((self.pyslip.tiles.picture_fast_slow_to_map_relative(beam_pixel_fast,beam_pixel_slow+3.),
        self.pyslip.tiles.picture_fast_slow_to_map_relative(beam_pixel_fast,beam_pixel_slow-3.)),
        {'width': 2, 'color': '#0000FFA0', 'closed': False})
                             ]
    # Unconditionally delete beam_layer and
    # spotfinder_layer--update_settings() will add them back if
    # appropriate.  This also creates the self.*_layer variables.
    if (hasattr(self, "beam_layer") and
        self.beam_layer is not None):
      self.pyslip.DeleteLayer(self.beam_layer)
    self.beam_layer = None

    if (hasattr(self, "spotfinder_layer") and
        self.spotfinder_layer is not None):
      self.pyslip.DeleteLayer(self.spotfinder_layer)
    self.spotfinder_layer = None

    if (hasattr(self, "tile_layer") and
        self.tile_layer is not None):
      self.pyslip.DeleteLayer(self.tile_layer)
    self.tile_layer = None

    if (hasattr(self, "tile_text_layer") and
        self.tile_text_layer is not None):
      self.pyslip.DeleteLayer(self.tile_text_layer)
    self.tile_text_layer = None

    self.update_settings()

    # Destroy the calibration frame if it present but unsupported for
    # this image.  XXX Need to do something about the ring tool too
    # when switching between different kinds of images.  XXX Centering
    # is broken when switching between different kinds of images.
    if (self._calibration_frame and
        not self.pyslip.tiles.raw_image.supports_quadrant_calibration()):
      self.OnCalibration(None)

  def update_settings (self, layout=True) :
    # XXX The zoom level from the settings panel are not taken into
    # account here.

    new_brightness = self.settings.brightness
    new_color_scheme = self.settings.color_scheme
    if new_brightness is not self.pyslip.tiles.current_brightness or \
       new_color_scheme is not self.pyslip.tiles.current_color_scheme:
      self.pyslip.tiles.update_brightness(new_brightness,new_color_scheme)
      self.pyslip.Update()#triggers redraw

    if (self.settings.show_beam_center):
      if (self.beam_layer is None):
        self.beam_layer = self.pyslip.AddPolygonLayer(
          self.beam_center_cross_data, name="<beam_layer>",
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    elif (self.beam_layer is not None):
      self.pyslip.DeleteLayer(self.beam_layer)
      self.beam_layer = None

    if (self.settings.show_spotfinder_spots):
      if (self.spotfinder_layer is None):
        tdata = self.pyslip.tiles.get_spotfinder_data(self.params)
        self.spotfinder_layer = self.pyslip.AddPointLayer(
          tdata, color="green", name="<spotfinder_layer>",
          radius=2,
          renderer = frame.pyslip.LightweightDrawPointLayer,
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
    elif (self.spotfinder_layer is not None):
      self.pyslip.DeleteLayer(self.spotfinder_layer)
      self.spotfinder_layer = None

    if (self.settings.show_effective_tiling):
      if (self.tile_layer is None):
        tdata, ttdata = self.pyslip.tiles.get_effective_tiling_data(self.params)
        self.tile_layer = self.pyslip.AddPolygonLayer(
          tdata, name="<tiling_layer>",
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5])
      if (self.tile_text_layer is None):
        self.tile_text_layer = self.pyslip.AddTextLayer(
          ttdata, name="<tiling_text_layer>",
          show_levels=[-2, -1, 0, 1, 2, 3, 4, 5],
          colour='#0000FFA0',
          textcolour='#0000FFA0',
          fontsize=30,
          placement='cc',
          radius=0)
    elif (self.tile_layer is not None) and (self.tile_text_layer is not None):
        self.pyslip.DeleteLayer(self.tile_layer)
        self.tile_layer = None
        self.pyslip.DeleteLayer(self.tile_text_layer)
        self.tile_text_layer = None

    if hasattr(self,"user_callback"):
      self.user_callback(self)

  def OnCalibration(self, event):
    from rstbx.slip_viewer.calibration_frame import SBSettingsFrame

    if (not self._calibration_frame):
      self._calibration_frame = SBSettingsFrame(
        self, -1, "Quadrant calibration", style=wx.CAPTION | wx.CLOSE_BOX)
      self._calibration_frame.Show()
      self._calibration_frame.Raise()
    else:
      self._calibration_frame.Destroy()


  def OnRing(self, event):
    from rstbx.slip_viewer.ring_frame import RingSettingsFrame

    if (not self._ring_frame):
      self._ring_frame = RingSettingsFrame(
        self, -1, "Ring tool", style=wx.CAPTION | wx.CLOSE_BOX)
      self._ring_frame.Show()
      self._ring_frame.Raise()
    else:
      self._ring_frame.Destroy()


  def OnUpdateUI(self, event):
    # Toggle the text of the menu items depending on the state of the
    # calibration and ring frames.  If quadrant calibration is not
    # supported for this image, disable the corresponding menu item.
    # Enable/disable previous and next buttons based on the image's
    # position in the list.

    eid = event.GetId()
    if (eid == self._id_calibration):

      if (self._calibration_frame):
        event.SetText("Hide quadrant calibration")
      else:
        event.SetText("Show quadrant calibration")

      if self.pyslip.tiles.raw_image.supports_quadrant_calibration():
        event.Enable(True)
      else:
        event.Enable(False)
      return

    elif (eid == self._id_ring):
      if (self._ring_frame):
        event.SetText("Hide ring tool")
      else:
        event.SetText("Show ring tool")
      return

    elif (eid == wx.ID_BACKWARD):
      if (self.image_chooser.GetSelection() - 1 >= 0):
        event.Enable(True)
      else:
        event.Enable(False)
      return

    elif (eid == wx.ID_FORWARD):
      if (self.image_chooser.GetSelection() + 1 <
          self.image_chooser.GetCount()):
        event.Enable(True)
      else:
        event.Enable(False)
      return


from rstbx.viewer.frame import SettingsFrame

def override_SF_set_image(self,image):
  self.Layout()
  self.Fit()
SettingsFrame.set_image = override_SF_set_image
