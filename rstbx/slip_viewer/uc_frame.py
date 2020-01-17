from __future__ import absolute_import, division, print_function
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id

import wx, math

class UCSettingsFrame(wx.MiniFrame):
  def __init__(self, *args, **kwds):
    super(UCSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.phil_params = args[0].params
    panel = UCSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class UCSettingsPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    super(UCSettingsPanel, self).__init__(*args, **kwds)

    self.phil_params = args[0].phil_params
    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    # Needed to draw and delete the rings.  XXX Applies to
    # calibration_frame as well?
    self._pyslip = self.GetParent().GetParent().pyslip

    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Number of decimal digits for distances.
    self.digits = 2

    # Wavelength control.
    beam = self._pyslip.tiles.raw_image.get_beam()
    self._wavelength = beam.get_wavelength()

    # Unit cell controls.
    if self.phil_params.calibrate_unitcell.unitcell is not None:
      self._cell = list(self.phil_params.calibrate_unitcell.unitcell.parameters())
    else:
      self._cell = [4.18,4.72,58.38,89.44,89.63,75.85]

    if self.phil_params.calibrate_unitcell.spacegroup is not None:
      self._spacegroup = self.phil_params.calibrate_unitcell.spacegroup
    else:
      self._spacegroup = "P1"

    self._cell_control_names = ["uc_a_ctrl","uc_b_ctrl","uc_c_ctrl",
                                "uc_alpha_ctrl","uc_beta_ctrl","uc_gamma_ctrl"]

    box = wx.BoxSizer(wx.HORIZONTAL)

    self.uc_a = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[0], value=self._cell[0])
    box.Add(self.uc_a,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="a"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_a)

    self.uc_alpha = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[3], value=self._cell[3])
    box.Add(self.uc_alpha,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="alpha"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_alpha)

    sizer.Add(box)


    box = wx.BoxSizer(wx.HORIZONTAL)

    self.uc_b = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[1], value=self._cell[1])
    box.Add(self.uc_b,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="b"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_b)

    self.uc_beta = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[4], value=self._cell[4])
    box.Add(self.uc_beta,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="beta"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_beta)

    sizer.Add(box)


    box = wx.BoxSizer(wx.HORIZONTAL)

    self.uc_c = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[2], value=self._cell[2])
    box.Add(self.uc_c,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="c"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_c)

    self.uc_gamma = FloatSpin(
      self, digits=self.digits, name=self._cell_control_names[5], value=self._cell[5])
    box.Add(self.uc_gamma,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="gamma"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCell, self.uc_gamma)

    sizer.Add(box)

    # Space group control
    box = wx.BoxSizer(wx.HORIZONTAL)

    self.space_group_ctrl =  wx.TextCtrl(
      self, name="space group", value=self._spacegroup)
    box.Add(self.space_group_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Space group"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_TEXT, self.OnSpaceGroup, self.space_group_ctrl)

    sizer.Add(box)


    # Distance control
    img = self.GetParent().GetParent()._img
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.distance_ctrl = FloatSpin(
          self, digits=self.digits, name="Detector Distance", value=img.get_detector_distance())
    self.distance_ctrl.SetIncrement(0.5)
    box.Add(self.distance_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)

    txtd = wx.StaticText(self, label="Detector Distance")
    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpin, self.distance_ctrl)
    sizer.Add(box)

    # Wavelength control
    img = self.GetParent().GetParent()._img
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.wavelength_ctrl = FloatSpin(
          self, digits=4, name="Wavelength", value=img.get_wavelength())
    self.wavelength_ctrl.SetIncrement(0.05)
    box.Add(self.wavelength_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)

    txtw = wx.StaticText(self, label="Wavelength")
    box.Add(txtw, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpin, self.wavelength_ctrl)
    sizer.Add(box)


    # d_min control
    if self.phil_params.calibrate_unitcell.d_min is not None:
      self.d_min = self.phil_params.calibrate_unitcell.d_min
    else:
      self.d_min = 3.5
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.d_min_ctrl = FloatSpin(
          self, digits=self.digits, name="d_min", value=self.d_min)
    box.Add(self.d_min_ctrl,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)

    txtd = wx.StaticText(self, label="Highest resolution for ring display")
    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpin, self.d_min_ctrl)
    sizer.Add(box)


    # Centering controls.
    self._center = [0, 0]
    box = wx.BoxSizer(wx.HORIZONTAL)

    self.spinner_fast = FloatSpin(
      self, digits=self.digits, name="fast_ctrl", value=self._center[0])
    box.Add(self.spinner_fast,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Center fast"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCenter, self.spinner_fast)

    self.spinner_slow = FloatSpin(
      self, digits=self.digits, name="slow_ctrl", value=self._center[1])
    box.Add(self.spinner_slow,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(wx.StaticText(self, label="Center slow"),
            0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpinCenter, self.spinner_slow)

    sizer.Add(box)

    self.DrawRings()


  def __del__(self):
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)


  def OnSpinCenter(self, event):
    obj = event.EventObject
    name = obj.GetName()

    if (name == "fast_ctrl"):
      self._center[0] = obj.GetValue()
    elif (name == "slow_ctrl"):
      self._center[1] = obj.GetValue()

    self.DrawRings()


  def OnSpinCell(self, event):
    obj = event.EventObject
    name = obj.GetName()

    self._cell[self._cell_control_names.index(name)] = obj.GetValue()

    self.DrawRings()


  def OnSpin(self, event):
    self.DrawRings()


  def OnSpaceGroup(self, event):
    obj = event.EventObject
    self._spacegroup = obj.GetValue()

    self.DrawRings()


  def _draw_rings_layer(self, dc, data, map_rel):
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
    dc.SetPen(wx.Pen(colour))
    dc.SetBrush(wx.Brush(colour, wx.TRANSPARENT))
    for (lon, lat, place, radius, colour, x_off, y_off, pdata) in data:
      (x, y) = self._pyslip.ConvertGeo2View((lon, lat))
      dc.DrawCircle(x, y, radius * scale)


  def DrawRings(self):
    from cctbx.crystal import symmetry
    import cctbx.miller

    frame = self.GetParent().GetParent()

    try:
      uc = symmetry(unit_cell=self._cell, space_group_symbol=str(self._spacegroup))
      hkl_list = cctbx.miller.build_set(uc, False, d_min=self.d_min_ctrl.GetValue())
    except Exception as e:
      frame.update_statusbar(str(e))
      return

    frame.update_statusbar("%d %d %d %d %d %d, "%tuple(self._cell) + "number of indices: %d"%len(hkl_list.indices()))

    spacings = list(hkl_list.d_spacings())
    print("Printing spacings, len: %s"%len(spacings))

    def cmp(a,b):
      if a[1] > b[1]: return 1
      elif a[1] < b[1]: return -1
      return 0

    spacings = sorted(spacings, cmp=cmp, reverse=True)

    for d in spacings:
      print(d)

    detector = self._pyslip.tiles.raw_image.get_detector()
    beam     = self._pyslip.tiles.raw_image.get_beam()

    wavelength = float(self.wavelength_ctrl.GetValue())
    distance = float(self.distance_ctrl.GetValue())
    pixel_size = detector[0].get_pixel_size()[0] # FIXME assumes square pixels, and that all panels use same pixel size

    twotheta = hkl_list.two_theta(wavelength = wavelength)
    L_mm = []
    L_pixels = []
    for tt in twotheta: L_mm.append(distance * math.tan(tt[1]))
    for lmm in L_mm: L_pixels.append(lmm/pixel_size)

    xrayframe = self.GetParent().GetParent()
    panel_id, beam_pixel_fast, beam_pixel_slow = xrayframe.get_beam_center_px()

    if len(detector) > 1:
      beam_pixel_slow, beam_pixel_fast = xrayframe.pyslip.tiles.flex_image.tile_readout_to_picture(
        panel_id, beam_pixel_slow - 0.5, beam_pixel_fast - 0.5)

    center = self._pyslip.tiles.picture_fast_slow_to_map_relative(
      beam_pixel_fast + self._center[0], beam_pixel_slow + self._center[1])

    # XXX Transparency?
    ring_data = [(center[0], center[1], {"colour": "red", "radius": pxl}) for pxl in L_pixels]

    # Remove the old ring layer, and draw a new one.
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)
      self._ring_layer = None
    self._ring_layer = self._pyslip.AddPointLayer(
      ring_data,
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      renderer=self._draw_rings_layer,
      name="<ring_layer>")
