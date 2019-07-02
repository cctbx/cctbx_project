from __future__ import absolute_import, division, print_function
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

import wx
from six.moves import range


class RingSettingsFrame(wx.MiniFrame):
  def __init__(self, *args, **kwds):
    super(RingSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = RingSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class RingSettingsPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    # XXX Support several rings.  Plot radial distribution somewhere
    # (not here), but maybe distribution along ring.  Drop-down menu
    # for ring center, and button to reset to beam center.

    super(RingSettingsPanel, self).__init__(*args, **kwds)

    # Needed to draw and delete the rings.  XXX Applies to
    # calibration_frame as well?
    self._pyslip = self.GetParent().GetParent().pyslip

    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Number of decimal digits for distances.
    self.digits = 2

    # Distance control XXX probably does not belong here
#    print "DISTANCE",self.GetParent().GetParent().viewer._img
#    box = wx.BoxSizer(wx.HORIZONTAL)
#    from wxtbx.phil_controls.floatctrl import FloatCtrl
#    from wxtbx.phil_controls import EVT_PHIL_CONTROL
#    self.distance_ctrl = FloatCtrl(self, pos=(300,180), size=(80,-1),
#    value=80.00,
#    name="Detector Distance")
#    self.distance_ctrl.SetMax(1000)
#    self.distance_ctrl.SetMin(5)
#    self.distance_ctrl.SetOptional(False)
#    box.Add(self.distance_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#    txtd = wx.StaticText(self, label="Detector Distance")
#    box.Add(txtd, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
#    s.Add(box)

    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    # XXX Should really make value be in Aangstroem resolution, and
    # have a non-linear slider.
    self._radius = 100
    self._center = [0, 0]
    radius_max = 2000
    radius_min = 10

    # Radius controls.
    box = wx.BoxSizer(wx.HORIZONTAL)

    self.slider = wx.Slider(self, maxValue=radius_max,
                            minValue=radius_min, size=(250, -1),
                            style=wx.SL_AUTOTICKS | wx.SL_HORIZONTAL,
                            value=self._radius)
    box.Add(self.slider,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_SLIDER, self.OnSlide, self.slider)

    self.spinner = FloatSpin(self, digits=self.digits, max_val=radius_max,
                             min_val=radius_min, value=self._radius)
    box.Add(self.spinner,
            0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(EVT_FLOATSPIN, self.OnSpin, self.spinner)

    self.auto = wx.Button(self, label="Auto fit")
    self.Bind(wx.EVT_BUTTON, self.OnAutoFit, self.auto)
    box.Add(self.auto, 0, wx.RIGHT | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)

    sizer.Add(box)

    # Centering controls.
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

    self.DrawRing()


  def __del__(self):
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)


  def OnSlide(self, event):
    # Keep slider and spinner synchronized.
    obj = event.EventObject # XXX Is this construct sane?  See below
                            # and in calibration_frame, too!
    self._radius = obj.GetValue()
    self.spinner.SetValue(self._radius)

    self.DrawRing()

  def OnAutoFit(self, event):
    import math
    jitter = 6

    detector = self._pyslip.tiles.raw_image.get_detector()
    beam     = self._pyslip.tiles.raw_image.get_beam()
    # FIXME assumes all detector elements use the same millimeter-to-pixel convention
    if detector[0].get_distance() > 0:
      if len(detector) > 1:
        h = detector.hierarchy()
        if len(h) > 0:
          beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(
            detector.hierarchy().get_beam_centre(beam.get_s0()))
        else:
          beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(
            detector[0].get_beam_centre(beam.get_s0()))
      else:
        beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(
          detector[0].get_beam_centre(beam.get_s0()))

    avg_distance = -sum([p.get_distance() for p in detector])/len(detector)

    beam_pixel_fast += self._center[0]
    beam_pixel_slow += self._center[1]

    def PointsOnCircle(center, radius, count):
      for r in range(count):
        t = (r/count)*2*math.pi
        yield (center[0] + (radius*math.cos(t)),
               center[1] + (radius*math.sin(t)))

    best = float("-inf")
    bestc = [self._center[0],self._center[1]]
    bestr = self._radius

    raw_data = self._pyslip.tiles.raw_image.get_raw_data()
    if not isinstance(raw_data, tuple):
      raw_data = (raw_data,)

    for j in range(-jitter, jitter, 1):
      j /= 2
      for i in range(-jitter, jitter, 1):
        i /= 2
        for r in range(-jitter, jitter, 1):
          r /= 2
          total = 0.0
          for point in PointsOnCircle((beam_pixel_fast+i,beam_pixel_slow+j),self._radius+r,360):
            mm = detector[0].pixel_to_millimeter(point)
            mm = (mm[0],mm[1],avg_distance)
            pid = detector.get_panel_intersection(mm)
            if pid >= 0:
              px = detector[pid].get_ray_intersection_px(mm)
              px = [int(round(px[0])),int(round(px[1]))]
              data = raw_data[pid]
              if px[0] >= 0 and px[0] < data.focus()[1] and px[1] >= 0 and px[1] < data.focus()[0]:
                total += data[px[1],px[0]]
          if total > best:
            best = total
            bestc = [self._center[0]+i,self._center[1]+j]
            bestr = self._radius+r
          print("r: % 3.1f, i: % 3.1f, j: % 3.1f, best: %f"%(r, i, j, best))
    print("DONE", bestc, bestr)
    self._radius = bestr
    self._center = bestc

    self.spinner.SetValue(bestr)
    self.spinner_fast.SetValue(bestc[0])
    self.spinner_slow.SetValue(bestc[1])

    self.DrawRing()

  def OnSpin(self, event):
    # Keep slider and spinner synchronized.  XXX OnSpinRadius()?
    obj = event.EventObject
    self._radius = obj.GetValue()
    self.slider.SetValue(self._radius)

    self.DrawRing()


  def OnSpinCenter(self, event):
    obj = event.EventObject
    name = obj.GetName()

    if (name == "fast_ctrl"):
      self._center[0] = obj.GetValue()
    elif (name == "slow_ctrl"):
      self._center[1] = obj.GetValue()

    self.DrawRing()


  def _draw_ring_layer(self, dc, data, map_rel):
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


  def DrawRing(self):
    detector = self._pyslip.tiles.raw_image.get_detector()
    beam     = self._pyslip.tiles.raw_image.get_beam()

    xrayframe = self.GetParent().GetParent()
    panel_id, beam_pixel_fast, beam_pixel_slow = xrayframe.get_beam_center_px()

    if len(detector) > 1:
      beam_pixel_slow, beam_pixel_fast = xrayframe.pyslip.tiles.flex_image.tile_readout_to_picture(
        panel_id, beam_pixel_slow - 0.5, beam_pixel_fast - 0.5)

    center = self._pyslip.tiles.picture_fast_slow_to_map_relative(
      beam_pixel_fast + self._center[0], beam_pixel_slow + self._center[1])

    # XXX Transparency?
    ring_data = [(center[0], center[1],
                  {"colour": "red", "radius": self._radius})]

    # Remove the old ring layer, and draw a new one.  XXX Why
    # disappears at highest levels?
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)
      self._ring_layer = None
    self._ring_layer = self._pyslip.AddPointLayer(
      ring_data,
      map_rel=True,
      visible=True,
      show_levels=[-3, -2, -1, 0, 1, 2, 3, 4, 5],
      renderer=self._draw_ring_layer,
      name="<ring_layer>")
