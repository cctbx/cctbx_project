
from libtbx.utils import Sorry
from libtbx.math_utils import ifloor, iceil
from libtbx.str_utils import format_value
import math

class image (object) :
  def __init__ (self, file_name) :
    self.file_name = file_name
    from iotbx.detectors import ImageFactory
    img = ImageFactory(file_name)
    img.read()
    self._raw = img
    print img.show_header()
    self._invert_beam_center = False
    #self.update_image()
    #self.convert_to_bitmap()

  def create_flex_image (self,
                         brightness=100) :
    # FIXME
    try :
      from labelit.detectors import FlexImage
    except ImportError, e :
      raise Sorry("Labelit not installed or not configured.")
    fi = FlexImage(
      rawdata=self._raw.linearintdata,
      binning=1,
      vendortype=self._raw.vendortype,
      brightness=brightness / 100.,
      saturation=self._raw.saturation)
    #from scitbx.array_family import flex
    #print flex.max(self._raw.linearintdata), flex.min(self._raw.linearintdata)
    fi.setWindow(0.0, 0.0, 1)
    fi.adjust()
    fi.prep_string()
    return fi

  def update_settings (self, **kwds) :
    self._bmp = None
    self._wx_img = None
    self.update_image(**kwds)

  def update_image (self, brightness=100, zoom_level=1, w=None, h=None,
      invert_beam_center=False) :
    self._invert_beam_center = invert_beam_center
    self._img = self.create_flex_image(brightness)
    x = self._img.ex_size2()
    y = self._img.ex_size1()
    self._size = (x, y)
    import wx
    wx_image = wx.EmptyImage(*(self._size))
    wx_image.SetData(self._img.export_string)
    self._full_mag = wx_image
    if (zoom_level == 0) and (not None in [w, h]) :
      scale = min(w / x, h / y)
      wx_image = wx_image.Scale(int(x*scale), int(y*scale),
        wx.IMAGE_QUALITY_NORMAL)
    elif (zoom_level == 1) :
      pass
    elif (zoom_level == 2) :
      x /= 2
      y /= 2
      wx_image = wx_image.Scale(x, y, wx.IMAGE_QUALITY_NORMAL)
    elif (zoom_level == 3) :
      x /= 4
      y /= 4
      wx_image = wx_image.Scale(x, y, wx.IMAGE_QUALITY_NORMAL)
    bmp = wx_image.ConvertToBitmap()
    self._bmp = bmp
    self._wx_img = wx_image

  def get_bitmap (self) :
    if (self._bmp is None) :
      self.convert_to_bitmap()
    return self._bmp

  def get_zoomed_region (self, x, y, zoom=1, mag=10, n_pixels=40) :
    import wx
    x, y = self.image_coords_as_pixel_coords(x, y)
    x0 = x - (n_pixels / 2)
    y0 = y - (n_pixels / 2)
    img = self._full_mag.GetSubImage((x0, y0, n_pixels, n_pixels))
    return img.Scale(n_pixels * mag, n_pixels * mag, wx.IMAGE_QUALITY_NORMAL)

  def get_size (self) :
    return self._bmp.GetSize()

  def get_original_size (self) :
    return self._size

  def get_zoom_level (self) :
    w, h = self.get_size()
    w0,h0 = self.get_original_size()
    return w / w0

  def get_detector_dimensions (self) :
    pixel_size = self._raw.parameters['PIXEL_SIZE']
    w = pixel_size * self._raw.parameters['SIZE2']
    h = pixel_size * self._raw.parameters['SIZE1']
    return (w, h)

  def image_coords_as_detector_coords (self, x, y) :
    dw, dh = self.get_detector_dimensions()
    w, h = self.get_size()
    x_frac = x / w
    y_frac = y / h
    x_detector = x_frac * dw
    y_detector = (1.0 - y_frac) * dh
    return x_detector, y_detector

  def detector_coords_as_image_coords (self, x, y) :
    dw, dh = self.get_detector_dimensions()
    w, h = self.get_size()
    x_frac = x / dw
    y_frac = - ((y / dh) - 1.0)
    x_point = x_frac * w
    y_point = y_frac * h
    return (x_point, y_point)

  def image_coords_as_pixel_coords (self, x, y) :
    zoom = self.get_zoom_level()
    x_ = ifloor(x / zoom)
    y_ = ifloor(y / zoom)
    return x_, y_

  def image_coords_as_array_coords (self, x, y) :
    x_, y_ = self.image_coords_as_pixel_coords(x, y)
    return y_, x_

  def get_beam_center (self) :
    # FIXME Pilatus and ADSC images appear to have different conventions???
    if (self._invert_beam_center) :
      center_x = self._raw.parameters['BEAM_CENTER_Y']
      center_y = self._raw.parameters['BEAM_CENTER_X']
    else :
      center_x = self._raw.parameters['BEAM_CENTER_X']
      center_y = self._raw.parameters['BEAM_CENTER_Y']
    x_, y_ = self.image_coords_as_detector_coords(center_x, center_y)
    print "beam center:", x_, y_, center_x, center_y
    return self.detector_coords_as_image_coords(center_x, center_y)

  def get_point_info (self, x, y) :
    x_point, y_point = self.image_coords_as_detector_coords(x, y)
    x0, y0 = self.detector_coords_as_image_coords(x_point, y_point)
    if (self._invert_beam_center) :
      center_x = self._raw.parameters['BEAM_CENTER_Y']
      center_y = self._raw.parameters['BEAM_CENTER_X']
    else :
      center_x = self._raw.parameters['BEAM_CENTER_X']
      center_y = self._raw.parameters['BEAM_CENTER_Y']
    dist = self._raw.parameters['DISTANCE']
    wavelength = self._raw.parameters['WAVELENGTH']
    r = math.sqrt((center_x - x_point)**2 + (center_y - y_point)**2)
    two_theta = math.atan(r / dist)
    if (two_theta == 0.0) :
      d_min = None
    else :
      d_min = wavelength / (2 * math.sin(two_theta / 2))
    slow, fast = self.image_coords_as_array_coords(x, y)
    try :
      intensity = self._raw.linearintdata[slow, fast]
    except IndexError :
      return None
    else :
      return point_info(slow, fast, intensity, d_min)

  def line_between_points (self, x1, y1, x2, y2, n_values=100) :
    x1_, y1_ = self.image_coords_as_array_coords(x1, y1)
    x2_, y2_ = self.image_coords_as_array_coords(x2, y2)
    delta_x = (x2_ - x1_) / (n_values - 1)
    delta_y = (y2_ - y1_) / (n_values - 1)
    vals = []
    d = self._raw.linearintdata
    for n in range(n_values) :
      x = x1_ + (n * delta_x)
      y = y1_ + (n * delta_y)
      x_1 = ifloor(x)
      x_2 = iceil(x)
      y_1 = ifloor(y)
      y_2 = iceil(y)
      v11 = d[(x_1, y_1)]
      v12 = d[(x_1, y_2)]
      v21 = d[(x_2, y_1)]
      v22 = d[(x_2, y_2)]
      if (x_2 == x_1) :
        if (y_2 == y_1) :
          vxy = v11
        else :
          vxy = ((v12 * (y - y_1)) + (v11 * (y_2 - y))) / (y_2 - y_1)
      elif (y_2 == y_1) :
        vxy =  ((v21 * (x - x_1)) + (v11 * (x_2 - x))) / (x_2 - x_1)
      else :
        dxdy = (y_2 - y_1) * (x_2 - x_1)
        vxy = ((v11 / dxdy) * (x_2 - x) * (y_2 - y)) + \
              ((v21 / dxdy) * (x - x_1) * (y_2 - y)) + \
              ((v12 / dxdy) * (x_2 - x) * (y - y_1)) + \
              ((v22 / dxdy) * (x - x_1) * (y - y_1))
      vals.append(vxy)
    return vals

class point_info (object) :
  def __init__ (self, slow, fast, intensity, d_min) :
    self.slow = slow
    self.fast = fast
    self.intensity = intensity
    self.d_min = d_min

  def format (self) :
    return "resolution = %s  intensity = %.2f  slow=%d  fast=%d" % (
      format_value("%.2f A", self.d_min), self.intensity, self.slow, self.fast)

class settings (object) :
  def __init__ (self) :
    self.zoom_level = 0
    self.brightness = 100
    self.show_beam_center = True
    self.invert_beam_center_axes = False
