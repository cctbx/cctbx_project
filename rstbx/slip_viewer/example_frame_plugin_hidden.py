from __future__ import absolute_import, division, print_function
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: ring_frame.py 18950 2013-12-20 20:23:08Z phyy-nx $

import wx

### Enable the plugin by renaming to end in "_frame_plugin.py"

class ExampleSettingsFrame(wx.MiniFrame):
  def __init__(self, *args, **kwds):
    super(ExampleSettingsFrame, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    panel = ExampleSettingsPanel(self)
    self.SetSizer(szr)
    szr.Add(panel, 1, wx.EXPAND)
    szr.Fit(panel)
    self.panel = panel
    self.sizer = szr
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt : self.Destroy(), self)


class ExampleSettingsPanel(wx.Panel):
  def __init__(self, *args, **kwds):

    super(ExampleSettingsPanel, self).__init__(*args, **kwds)

    # Needed to draw and delete the overlay
    self._pyslip = self.GetParent().GetParent().pyslip

    sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(sizer)

    # Number of decimal digits
    self.digits = 2

    from wx.lib.agw.floatspin import EVT_FLOATSPIN, FloatSpin

    # Set initial values
    self._radius = 100
    self._center = [0, 0]
    radius_max = 2000
    radius_min = 10

    # Bind data to controls -- duplicate this section for each control
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

    sizer.Add(box)
    # end duplicate section

    # Update
    self.DrawRing()


  def __del__(self):
    # Delete layer method
    if (hasattr(self, "_ring_layer") and self._ring_layer is not None):
      self._pyslip.DeleteLayer(self._ring_layer)


  def OnSlide(self, event):
    # Keep slider and spinner synchronized
    obj = event.EventObject
    self._radius = obj.GetValue()
    self.spinner.SetValue(self._radius)

    self.DrawRing()

  def OnSpin(self, event):
    # Keep slider and spinner synchronized
    obj = event.EventObject
    self._radius = obj.GetValue()
    self.slider.SetValue(self._radius)

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

    if len(detector) > 1:
      beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(  # FIXME assumes all detector elements use the same
        detector.hierarchy().get_beam_centre(beam.get_s0()))               # millimeter-to-pixel convention
    else:
      beam_pixel_fast, beam_pixel_slow = detector[0].millimeter_to_pixel(
        detector[0].get_beam_centre(beam.get_s0()))

    center = self._pyslip.tiles.picture_fast_slow_to_map_relative(
      beam_pixel_fast + self._center[0], beam_pixel_slow + self._center[1])

    ring_data = [(center[0], center[1],
                  {"colour": "red", "radius": self._radius})]

    # Remove the old ring layer, and draw a new one.
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

class PluginHelper(object):
  _plugin_layer = "_example_layer"
  _plugin_title = "Example ring-drawing tool"
  _plugin_hide_text = "Hide example tool"
  _plugin_show_text = "Show example tool"
  _plugin_settings_frame = ExampleSettingsFrame
  _plugin_settings_panel = ExampleSettingsPanel
