from __future__ import division
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id
from six.moves import range

"""
Plugin to calibrate kapton correction parameters for an XFEL image
"""

import wx
from dials.algorithms.integration.kapton_2019_correction import KaptonTape_2019
from scitbx.matrix import col

class KaptonSettingsFrame(wx.MiniFrame):
  def __init__ (self, *args, **kwds):
    super(KaptonSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.phil_params = args[0].params
    panel = KaptonSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class KaptonSettingsPanel(wx.Panel):
  def __init__ (self, *args, **kwds) :
    super(KaptonSettingsPanel, self).__init__(*args, **kwds)

    # self.phil_params = args[0].phil_params
    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    # Needed to draw and delete overlays
    self._pyslip = self.GetParent().GetParent().pyslip

    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Number of decimal digits for distances.
    self.digits = 4

    # Kapton tape controls and params from phil
    self._xtal_height_mm = 0.035 # LN84 shift 2
    self._tape_thickness_mm = 0.05 # LN84 kapton tape
    self._tape_half_width_mm = 3.175/4. # LN84 kapton tape
    self._tape_angle_deg = 178.8-180.0 # LN84 shift 2 tape drive setup
    self._center = [0, 0]
    self.frame          = self.GetParent().GetParent()
    self.detector       = self._pyslip.tiles.raw_image.get_detector()
    self.panels         = [p for p in self.detector]
    self.beam           = self._pyslip.tiles.raw_image.get_beam()
    self.wavelength_ang = self.beam.get_wavelength()
    self.pixels_size_mm = [p.get_pixel_size()[0] for p in self.panels]
    self.distance       = [p.get_directed_distance() for p in self.panels]




    # determine beam center and detector limits in pixels
    # data reads slow,fast but the dxtbx objects are organized fast,slow -- label explicitly
    self.panel_size = [p.get_image_size() for p in self.panels]
    #self.size_fast, self.size_slow = self.detector[0].get_image_size()
    origin = [p.get_parent_origin()[:2] for p in self.panels]
    #orig_fast, orig_slow = self.detector[0].get_parent_origin()[:2]
    #self.min_px=[map(int,self.panels[ii].millimeter_to_pixel(ori)) for ii,ori in enumerate(origin)]
    self.min_px = [
        [int(x) for x in self.panels[ii].millimeter_to_pixel(ori)]
        for ii, ori in enumerate(origin)
        ]
    #self.fast_min, self.slow_min = map(int, self.detector[0].millimeter_to_pixel((orig_fast, orig_slow)))
    self.max_px=[(self.panel_size[ii][0]+min_px[0]-1,self.panel_size[ii][1]+min_px[1]-1) for ii,min_px in enumerate(self.min_px)]
    #self.fast_max, self.slow_max = map(lambda i: i-1,
    #  (self.size_fast + self.fast_min, self.size_slow + self.slow_min))
    self.s0=[map(int, p.get_beam_centre_px(self.beam.get_s0())) for p in self.panels]
    #self.s0_fast, self.s0_slow = map(int, self.detector[0].get_beam_centre_px(self.beam.get_s0()))
    panel_id, beam_pixel_fast, beam_pixel_slow = self.frame.get_beam_center_px()

    if len(self.detector) > 1:
      beam_pixel_slow, beam_pixel_fast = self.frame.pyslip.tiles.flex_image.tile_readout_to_picture(
        panel_id, beam_pixel_slow - 0.5, beam_pixel_fast - 0.5)

    self.center_lon_lat = self._pyslip.tiles.picture_fast_slow_to_map_relative(
      beam_pixel_fast + self._center[0], beam_pixel_slow + self._center[1])

    #self.center_lon_lat = self._pyslip.tiles.picture_fast_slow_to_map_relative(self.s0_fast, self.s0_slow)

    self._absorption = None

    self._kapton_control_names = ["xtal_height_ctrl", "tape_thickness_ctrl", "tape_width_ctrl",
    "tape_angle_ctrl", "show_shadow_button", "show_contours_button", "show_edge_max_button",
    "save_params"]

    box = wx.BoxSizer(wx.HORIZONTAL)
    self.height = FloatSpin(
      self, digits=self.digits, name=self._kapton_control_names[0], value=self._xtal_height_mm)
    self.height.SetIncrement(0.005)
    self.height.SetRange(0.005, None)
    box.Add(self.height,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Crystal height (mm)"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinHeight, self.height)
    sizer.Add(box)


    box = wx.BoxSizer(wx.HORIZONTAL)
    self.thickness = FloatSpin(
      self, digits=self.digits, name=self._kapton_control_names[1], value=self._tape_thickness_mm)
    self.thickness.SetIncrement(0.001)
    self.thickness.SetRange(0.001, None)
    box.Add(self.thickness,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Kapton tape thickness (mm)"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinThickness, self.thickness)
    sizer.Add(box)


    box = wx.BoxSizer(wx.HORIZONTAL)
    self.width = FloatSpin(
      self, digits=self.digits, name=self._kapton_control_names[2], value=self._tape_half_width_mm)
    self.width.SetIncrement(0.005)
    self.width.SetRange(0.005, None)
    box.Add(self.width,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Kapton tape half width (mm)"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinWidth, self.width)
    sizer.Add(box)


    box = wx.BoxSizer(wx.HORIZONTAL)
    self.angle = FloatSpin(
      self, digits=self.digits, name=self._kapton_control_names[3], value=self._tape_angle_deg)
    self.angle.SetIncrement(0.05)
    self.angle.SetRange(-360, 360)
    box.Add(self.angle,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Kapton tape counterclockwise rotation from vertical (deg)"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinAngle, self.angle)
    sizer.Add(box)

    box = wx.BoxSizer(wx.HORIZONTAL)
    self.show_shadow = wx.Button(
      self, name=self._kapton_control_names[4], label="show kapton absorption as shadow")
    box.Add(self.show_shadow,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnShowShadow, self.show_shadow)
    sizer.Add(box)

    box = wx.BoxSizer(wx.HORIZONTAL)
    self.show_contours = wx.Button(
      self, name=self._kapton_control_names[5], label="show kapton absorption as contours")
    box.Add(self.show_contours,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnShowContours, self.show_contours)
    sizer.Add(box)

    box = wx.BoxSizer(wx.HORIZONTAL)
    self.show_edge_max = wx.Button(
      self, name=self._kapton_control_names[6], label="show kapton absorption edge and max")
    box.Add(self.show_edge_max,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnShowEdgeMax, self.show_edge_max)
    sizer.Add(box)

    box = wx.BoxSizer(wx.HORIZONTAL)
    self.save_phil = wx.Button(
      self, name=self._kapton_control_names[7], label="save params and go to next image")
    box.Add(self.save_phil,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnSavePhil, self.save_phil)
    sizer.Add(box)

  def update_absorption(self):
    self._absorption = KaptonTape_2019(self._xtal_height_mm,
                                        self._tape_thickness_mm,
                                        self._tape_half_width_mm,
                                        self._tape_angle_deg,
                                        wavelength_ang=self.wavelength_ang)

  def __del__(self):
    if (hasattr(self, "_shadow_layer") and self._shadow_layer is not None):
      self._pyslip.DeleteLayer(self._shadow_layer)
    if (hasattr(self, "_contours_layer") and self._contours_layer is not None):
      self._pyslip.DeleteLayer(self._contours_layer)
    if (hasattr(self, "_edge_max_layer") and self._edge_max_layer is not None):
      self._pyslip.DeleteLayer(self._edge_max_layer)

  def OnSpinHeight(self, event):
    obj = event.EventObject
    name = obj.GetName()

    self._xtal_height_mm = obj.GetValue()
    self.UpdateAbsorptionData()

  def OnSpinThickness(self, event):
    obj = event.EventObject
    name = obj.GetName()

    self._tape_thickness_mm = obj.GetValue()
    self.UpdateAbsorptionData()

  def OnSpinWidth(self, event):
    obj = event.EventObject
    name = obj.GetName()

    self._tape_half_width_mm = obj.GetValue()
    self.UpdateAbsorptionData()

  def OnSpinAngle(self, event):
    obj = event.EventObject
    name = obj.GetName()

    self._tape_angle_deg = obj.GetValue()
    self.UpdateAbsorptionData()

  def OnSpin(self, event):
    self.UpdateAbsorptionData()

  def OnShowShadow(self, event):
    if (hasattr(self, "_shadow_layer") and self._shadow_layer is not None):
      self._pyslip.DeleteLayer(self._shadow_layer)
      self._shadow_layer = None
    else:
      self.DrawShadow()

  def OnShowContours(self, event):
    if (hasattr(self, "_contours_layer") and self._contours_layer is not None):
      self._pyslip.DeleteLayer(self._contours_layer)
      self._contours_layer = None
    else:
      self.DrawContours()

  def OnShowEdgeMax(self, event):
    if (hasattr(self, "_edge_max_layer") and self._edge_max_layer is not None):
      self._pyslip.DeleteLayer(self._edge_max_layer)
      self._edge_max_layer = None
    else:
      self.DrawEdgeMax()

  def OnSavePhil(self, event):
    import os
    if not (hasattr(self, "_phil_destination_dir") and self._phil_destination_dir is not None):
      self._run = os.path.basename(self.frame._img._raw.path).split(".pickle")[0].split("-")[-1]
      name = "Kapton_correction_%s.phil" % self._run
      dialog = wx.FileDialog(
        self,
        defaultDir='',
        defaultFile=name,
        message="Choose Kapton parameter snippet destination",
        style=wx.FD_SAVE,
        wildcard="*")
      if dialog.ShowModal() == wx.ID_OK:
        dest_path = dialog.GetPath()
        if os.path.isdir(dest_path):
          self._phil_destination_dir = dest_path
          self._phil_destination_path = os.path.join(self._phil_destination_dir, name)
        else:
          self._phil_destination_dir = os.path.dirname(dest_path)
          self._phil_destination_path = dest_path
    else:
      prev_run = self._run
      self._run = os.path.basename(self.frame._img._raw.path).split(".pickle")[0].split("-")[-1]
      name_parts = os.path.basename(self._phil_destination_path).split(prev_run)
      name = name_parts[0] + self._run + name_parts[1]
      self._phil_destination_path = os.path.join(self._phil_destination_dir, name)
    wb = open(self._phil_destination_path, "wb")
    wb.write(b"integration {\n")
    wb.write(b"  absorption_correction {\n")
    wb.write(b"    apply = True\n")
    wb.write(b"    algorithm = fuller_kapton\n")
    wb.write(b"    fuller_kapton.xtal_height_above_kapton_mm {\n")
    wb.write(b"      value = %f\n" % self._xtal_height_mm)
    wb.write(b"      }\n")
    wb.write(b"    fuller_kapton.rotation_angle_deg {\n")
    wb.write(b"      value = %f\n" % self._tape_angle_deg)
    wb.write(b"      }\n")
    wb.write(b"    fuller_kapton.kapton_half_width_mm {\n")
    wb.write(b"       value = %f\n" % self._tape_half_width_mm)
    wb.write(b"      }\n")
    wb.write(b"    fuller_kapton.kapton_thickness_mm {\n")
    wb.write(b"       value = %f\n" % self._tape_thickness_mm)
    wb.write(b"      }\n")
    wb.write(b"    fuller_kapton.smart_sigmas = True\n")
    wb.write(b"  }\n")
    wb.write(b"}\n")
    wb.close()
    print("Wrote phil snippet to", self._phil_destination_path)
    if not hasattr(self, "_root_frame") or self._root_frame is None:
      self._root_frame = self._pyslip.GetParent().GetParent()
    self._root_frame.OnNext(event)

  def _draw_shadow_layer(self, dc, data, map_rel):
    """Draw a points layer.

    dc       the device context to draw on
    data     an iterable of point tuples:
             (x, y, place, radius, colour, x_off, y_off, pdata)
    map_rel  points relative to map if True, MUST BE TRUE for lightweight
    Assumes all points are the same colour, saving 100's of ms.
    """

    assert map_rel is True
    if len(data)==0:
      return
    (lon, lat, place, radius, colour, x_off, y_off, pdata) = data[0]

    scale = 2**self._pyslip.tiles.zoom_level

    # Draw points on map/view, using transparency if implemented.
    try:
      dc = wx.GCDC(dc)
    except NotImplementedError:
      pass
    dc.SetBrush(wx.Brush(colour, wx.TRANSPARENT))
    for (lon, lat, place, radius, colour, x_off, y_off, pdata) in data:
      dc.SetPen(wx.Pen(colour))
      (x, y) = self._pyslip.ConvertGeo2View((lon, lat))
      dc.DrawCircle(x, y, 1)
      if colour == "purple":
        dc.DrawCircle(x, y, 10)
      pass

  def _draw_contours_layer(self, dc, data, map_rel):
    """Draw a layer consisting of differently colored contours.

    dc       the device context to draw on
    data     a list of dictionaries:
             {"Paths":[Path, Path, ...], "level":level, "color":color}
             each Path is an iterable of line segments to be overlaid on the pyslip image
             level is the contour level of the Paths
             colour is the colour to be applied to the contour lines at this level
    map_rel  points relative to map if True, MUST BE TRUE for lightweight
    """

    assert map_rel is True
    if len(data) == 0:
      return
    scale = 2**self._pyslip.tiles.zoom_level

    # Draw line segments of the appropriate color on the pyslip map object
    try:
      dc = wx.GCDC(dc)
    except NotImplementedError:
      pass
    for contour_level in data:
      path_list = contour_level["Paths"]
      level = contour_level["level"]
      colour = wx.Colour(*[(1-decimal)*255 for decimal in contour_level["colour"]])
      dc.SetBrush(wx.Brush(colour))
      dc.SetPen(wx.Pen(colour))
      for path in path_list:
        segments = path.iter_segments()
        this_vertex = segments.next()[0]
        vertex_1 = self._pyslip.ConvertGeo2View(this_vertex) # x, y = Geo2View(lon, lat)
        while True:
          try:
            this_vertex = segments.next()[0]
            vertex_2 = self._pyslip.ConvertGeo2View(this_vertex) # x, y = Geo2View(lon, lat)
            dc.DrawLine(vertex_1[0], vertex_1[1], vertex_2[0], vertex_2[1])
            vertex_1 = vertex_2
          except StopIteration:
            break

  def _draw_edge_max_layer(self, dc, data, map_rel):
    """Draw a layer consisting of two line segments, identifying the Kapton absorption edge and max.

    dc       the device context to draw on
    data     a list of line segments (lon1, lat1, lon2, lat2)
    map_rel  points relative to map if True, MUST BE TRUE for lightweight
    """
    assert map_rel is True
    if len(data) == 0:
      return
    scale = 2**self._pyslip.tiles.zoom_level

    # Draw line segments on the pyslip map object
    try:
      dc = wx.GCDC(dc)
    except NotImplementedError:
      pass
    dc.SetBrush(wx.Brush("black"))
    dc.SetPen(wx.Pen("black"))
    for segment in data:
      lon1, lat1, lon2, lat2 = segment
      x1, y1 = self._pyslip.ConvertGeo2View((lon1, lat1))
      x2, y2 = self._pyslip.ConvertGeo2View((lon2, lat2))
      dc.DrawLine(x1, y1, x2, y2)

  def _map_coords(self, x, y, p):
      y, x = self._pyslip.tiles.flex_image.tile_readout_to_picture(
          p, y - 0.5, x - 0.5)
      return self._pyslip.tiles.picture_fast_slow_to_map_relative(
        x, y)

  def UpdateAbsorptionData(self, edge_max_mode=False):
    from scitbx.array_family import flex
    #from IPython import embed; embed(); exit()
    self.update_absorption() # basically resets the self._absorption to a new KaptonTape class
    self.absorption_data = []
    if True: #not hasattr(self, 'cache_mode'):
      self.cache_mode=True
      self.all_s1_flex=flex.vec3_double()
      #self.all_lon_lat=[]
    #FIXME
      detector=self._pyslip.tiles.raw_image.get_detector()
      for panel in detector:
        image_size = panel.get_image_size()
        for f_dir in range(0, image_size[0], 10):
          for s_dir in range(0, image_size[1], 10):
            s1=panel.get_pixel_lab_coord((f_dir, s_dir))
            #self.all_s1_flex.append(s1)
            lon, lat = self._map_coords(float(f_dir), float(s_dir),panel.index())
            #self.all_lon_lat.append((lon,lat))
      #self.cache_mode=True

            absorption_correction = self._absorption.abs_correction(s1)
            self.absorption_data.append((lon, lat, absorption_correction))
    #absorption_corrections = self._absorption.abs_correction_flex(self.all_s1_flex)
    #self.absorption_data=[(x[0],x[1],y) for x,y in zip(self.all_lon_lat, absorption_corrections)]

    # Redraw layers if present
    if (hasattr(self, "_shadow_layer") and self._shadow_layer is not None):
      self.DrawShadow()
    if (hasattr(self, "_contours_layer") and self._contours_layer is not None):
      self.DrawContours()
    if (hasattr(self, "_edge_max_layer") and self._edge_max_layer is not None):
      self.DrawEdgeMax()

  def DrawShadow(self):
    if not hasattr(self, 'absorption_data'):
      self.UpdateAbsorptionData()
    shadow_data = []
    shadow_data.append((self.center_lon_lat[0], self.center_lon_lat[1], {"colour":"purple"}))
    for lon, lat, absorption_correction in self.absorption_data:
      if absorption_correction == 1:
        colour = "white"
      elif absorption_correction < 1.2:
        colour = "blue"
      elif absorption_correction < 1.4:
        colour = "green"
      elif absorption_correction < 1.6:
        colour = "yellow"
      elif absorption_correction < 1.8:
        colour = "orange"
      else:
        colour = "red"
      shadow_data.append((lon, lat, {"colour":colour}))

    # Remove the old shadow layer, and draw a new one.
    if (hasattr(self, "_shadow_layer") and self._shadow_layer is not None):
      self._pyslip.DeleteLayer(self._shadow_layer)
      self._shadow_layer = None
    self._shadow_layer = self._pyslip.AddPointLayer(
      shadow_data,
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      renderer=self._draw_shadow_layer,
      name="<shadow_layer>")

  def DrawContours(self):
    if not hasattr(self, 'absorption_data'):
      self.UpdateAbsorptionData()
    import matplotlib.pyplot as plt
    print ('DrawContours does not work currently. Please try the other options')
    return
    def convert_to_tile(slow_fast_ab_tup):
      slow, fast, ab = slow_fast_ab_tup
      local_lon, local_lat = self._pyslip.tiles.picture_fast_slow_to_map_relative(fast, slow)
      return (local_lon, local_lat, ab)
    lon, lat, absorptions = zip(*map(convert_to_tile, self.absorption_data))
    n_lon = len(set(lon))
    n_lat = len(set(lat))
    indices = []
    if lon[0] == lon[1]:
      indices = [[i_lon*n_lat + i_lat
        for i_lat in range(0, n_lat, 1)]
        for i_lon in range(0, n_lon, 1)]
    elif lat[0] == lat[1]:
      indices = [[i_lat*n_lon + i_lon
        for i_lon in range(0, n_lon, 1)]
        for i_lat in range(0, n_lat, 1)]
    else: # appears to be a non-rectangular array
      indices = [range(n_lon*n_lat)]
    lon_array = [[lon[i] for i in l] for l in indices]
    lat_array = [[lat[i] for i in l] for l in indices]
    abs_array = [[absorptions[i] for i in l] for l in indices]
    contours = plt.contour(lon_array, lat_array, abs_array)
    levels = list(contours.levels)
    def convert_level_to_rgb_colour(level):
      r, g, b, a = contours.cmap(level - 1)
      return (r, g, b)
    colours = map(convert_level_to_rgb_colour, levels)
    contour_paths = [{"Paths":contours.collections[i].get_paths(),
                     "level":levels[i],
                     "colour":colours[i]} for i in range(len(contours.collections))]

    # Remove the old contours layer, and draw a new one.
    if (hasattr(self, "_contours_layer") and self._contours_layer is not None):
      self._pyslip.DeleteLayer(self._contours_layer)
      self._contours_layer = None
    self._contours_layer = self._pyslip.AddLayer(
      self._draw_contours_layer,
      contour_paths,
      True, #map_rel
      True, #visible
      [-3, -2, -1, 0, 1, 2, 3, 4, 5], #show_levels
      False, #selectable
      "<contours_layer>", #name
      "contours")

  def DrawEdgeMax(self):
    if not hasattr(self, 'absorption_data'):
      self.UpdateAbsorptionData(edge_max_mode=True)
    detector=self._pyslip.tiles.raw_image.get_detector()
    beam=self._pyslip.tiles.raw_image.get_beam()
    xrayframe=self.GetParent().GetParent()
    int_with_det = self._absorption.abs_bounding_lines_on_image(detector)
    s0=beam.get_s0()
    panel_id, beam_pixel_fast, beam_pixel_slow=xrayframe.get_beam_center_px()
    if len(detector) > 1:
      beam_pixel_slow, beam_pixel_fast = xrayframe.pyslip.tiles.flex_image.tile_readout_to_picture(
        panel_id, beam_pixel_slow - 0.5, beam_pixel_fast - 0.5)
    int_with_det_px = []
    edge_max_data_slow_fast = []

    def panel_of_intersection(elem, detector):
      """Takes in  x1,y1,x2,y2 and returns the panel where these points lie. If outside any panel then returns None"""
      x1,y1,x2,y2=elem
      # FIXME tolerance levels, a bit hacky
      d1=3.09
      d2=-0.09
      int_panels=[None, None]
      for panel in detector:
        # Convert x,y to f,s for that panel
        z_panel=panel.get_distance()
        ori=col(panel.get_origin())
        r1 = col((x1,y1,-z_panel))-ori
        r2 = col((x2,y2,-z_panel))-ori
        f_u = col(panel.get_fast_axis())
        s_u = col(panel.get_slow_axis())
        f1 = r1.dot(f_u)
        s1 = r1.dot(s_u)
        f2 = r2.dot(f_u)
        s2 = r2.dot(s_u)
        detz=panel.get_distance()
        px1=panel.get_ray_intersection_px((f1,s1, -detz))
        px2=panel.get_ray_intersection_px((f2,s2, -detz))
        if px1[0]/panel.get_image_size()[0]<=d1 and px1[0]/panel.get_image_size()[0] >=d2 and \
           px1[1]/panel.get_image_size()[1]<=d1 and px1[1]/panel.get_image_size()[1] >=d2:
           int_panels[0]=panel
        if px2[0]/panel.get_image_size()[0]<=d1 and px2[0]/panel.get_image_size()[0] >=d2 and \
           px2[1]/panel.get_image_size()[1]<=d1 and px2[1]/panel.get_image_size()[1] >=d2:
           int_panels[1]=panel
      #from IPython import embed; embed(); exit()
      # FIXME
      if int_panels[0] is None and int_panels[1] is None:
        from libtbx.utils import Sorry
        raise Sorry("No intersecting panels ?? this is awful .. get out of here")
      if int_panels[0] is None: int_panels[0]=int_panels[1]
      if int_panels[1] is None: int_panels[1]=int_panels[0]
      return int_panels

    for elem in int_with_det:
      # Determine panel number for pts 1 and 2
      p1,p2=panel_of_intersection(elem, detector)
      pixel_size_1=p1.get_pixel_size()[0]
      pixel_size_2=p2.get_pixel_size()[0]
      if True:
        # for panel 1
        z_panel=p1.get_distance()
        ori=col(p1.get_origin())
        r1 = col((elem[0], elem[1], -z_panel))-ori
        f1_u = col(p1.get_fast_axis())
        s1_u = col(p1.get_slow_axis())
        f1=r1.dot(f1_u)
        s1=r1.dot(s1_u)
        # For panel 2
        z_panel=p2.get_distance()
        ori=col(p2.get_origin())
        r2 = col((elem[2], elem[3], -z_panel))-ori
        f2_u = col(p2.get_fast_axis())
        s2_u = col(p2.get_slow_axis())
        f2=r2.dot(f2_u)
        s2=r2.dot(s2_u)

        if True:
          # FIXME is this correct to comment out beam_pixel_fast / beam_pixel_slow
          px_f1=f1/pixel_size_1#+beam_pixel_fast
          px_s1=s1/pixel_size_1#+beam_pixel_slow
          px_f2=f2/pixel_size_2#+beam_pixel_fast
          px_s2=s2/pixel_size_2#+beam_pixel_slow
      # If debugging
      if False:
        px_f1=beam_pixel_fast
        px_s1=beam_pixel_slow
        px_f2=beam_pixel_fast+100
        px_s2=beam_pixel_slow
      edge_max_data_slow_fast.append((px_s1, px_f1, p1.index(), px_s2, px_f2, p2.index()))

    edge_max_data_lonlat = []
    kapton_data=[]

    for segment in edge_max_data_slow_fast:
      slow1, fast1, p1, slow2, fast2, p2= segment
      lon1, lat1 = self._map_coords(fast1,slow1,p1)#self._pyslip.tiles.picture_fast_slow_to_map_relative(fast1, slow1)
      lon2, lat2 = self._map_coords(fast2, slow2, p2)#self._pyslip.tiles.picture_fast_slow_to_map_relative(fast2, slow2)
      edge_max_data_lonlat.append((lon1, lat1, lon2, lat2))
      shoebox_dict = {'width': 2, 'color': '#0000FFA0', 'closed': False}

      if False:
        x0_, y0_=lon1,lat1
        x1_, y1_=lon2,lat2
        my_attrs = dict(shoebox_dict)
        lines = [(((x0_, y0_), (x0_, y1_)), my_attrs),
                 (((x0_, y1_), (x1_, y1_)), my_attrs),
                 (((x1_, y1_), (x1_, y0_)), my_attrs),
                 (((x1_, y0_), (x0_, y0_)), my_attrs)]
        kapton_data.extend(lines)

    if False:
      self.edge_max_layer = self._pyslip.AddPolygonLayer(
          kapton_data, map_rel=True, visible=True,
          show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
          selectable=False,
          name='<kapton_layer>', update=False)

    # Remove the old edge_max layer, and draw a new one.
    if (hasattr(self, "_edge_max_layer") and self._edge_max_layer is not None):
      self._pyslip.DeleteLayer(self._edge_max_layer)
      self._edge_max_layer = None
    if True:
      self._edge_max_layer = self._pyslip.AddLayer(
        self._draw_edge_max_layer,
        edge_max_data_lonlat,
        map_rel=True, #map_rel
        visible=True, #visible
        show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5], #show_levels
        selectable=False, #selectable
        name="<edge_max_layer>",
        type="edge_max")


class PluginHelper(object):
  _plugin_layer = "_kapton_layer"
  _plugin_title = "Kapton 2019 correction tool"
  _plugin_hide_text = "Hide 2019 kapton tool"
  _plugin_show_text = "Show 2019 kapton tool"
  _plugin_settings_frame = KaptonSettingsFrame
  _plugin_settings_panel = KaptonSettingsPanel
