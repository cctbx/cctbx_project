
from libtbx.utils import Sorry
from libtbx.math_utils import ifloor, iceil
import math

class image (object) :
  def __init__ (self, file_name) :
    self.file_name = file_name
    from iotbx.detectors import ImageFactory
    img = ImageFactory(file_name)
    img.read()
    self._raw = img
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
    fi.setWindow(0.0, 0.0, 1)
    fi.adjust()
    fi.prep_string()
    return fi

  def update_settings (self, **kwds) :
    self._bmp = None
    self._wx_img = None
    self.update_image(**kwds)

  def update_image (self, brightness=100, zoom_level=1, w=None, h=None) :
    self._img = self.create_flex_image(brightness)
    x = self._img.ex_size2()
    y = self._img.ex_size1()
    self._size = (x, y)
    import wx
    wx_image = wx.EmptyImage(*(self._size))
    wx_image.SetData(self._img.export_string)
    self._full_mag = wx_image
    if (zoom_level == 0) and (not None in [w, h]) :
      wx_image = wx_image.Scale(min(w,h), min(w,h), wx.IMAGE_QUALITY_NORMAL)
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
    x, y = self.image_coords_as_array_coords(x, y, zoom)
    x0 = x - (n_pixels / 2)
    y0 = y - (n_pixels / 2)
    img = self._full_mag.GetSubImage((x0, y0, n_pixels, n_pixels))
    return img.Scale(n_pixels * mag, n_pixels * mag, wx.IMAGE_QUALITY_NORMAL)

  def get_size (self) :
    return self._bmp.GetSize()

  def image_coords_as_detector_coords (self, x, y) :
    pixel_size = self._raw.parameters['PIXEL_SIZE']
    detector_width = self._raw.parameters['SIZE2']
    detector_height = self._raw.parameters['SIZE1']
    screen_width = pixel_size * detector_width
    screen_height = pixel_size * detector_height
    w, h = self.get_size()
    x_frac = x / w
    y_frac = y / h
    x_detector = x_frac * screen_width
    y_detector = (1 - y_frac) * screen_height
    return x_detector, y_detector

  def image_coords_as_array_coords (self, x, y, zoom=1) :
    if (zoom == 0) :
      w, h = self._bmp.GetSize()
      w0, h0 = self._size
      x_ = x / (w / w0)
      y_ = y / (h / h0)
    else :
      x_ = x / zoom
      y_ = y / zoom
    return x_, y_

  def detector_coords_as_image_coords (self, x, y) :
    pixel_size = self._raw.parameters['PIXEL_SIZE']
    detector_width = self._raw.parameters['SIZE2']
    detector_height = self._raw.parameters['SIZE1']
    screen_width = pixel_size * detector_width
    screen_height = pixel_size * detector_height
    x_frac = x / screen_width
    y_frac = 1 - (y / screen_height)
    w, h = self.get_size()
    return (x_frac * w), (y_frac * h)

  def get_beam_center (self) :
    center_x = self._raw.parameters['BEAM_CENTER_X']
    center_y = self._raw.parameters['BEAM_CENTER_Y']
    return self.detector_coords_as_image_coords(center_x, center_y)

  def get_d_min_at_point (self, x, y) :
    x_point, y_point = self.image_coords_as_detector_coords(x, y)
    center_x = self._raw.parameters['BEAM_CENTER_X']
    center_y = self._raw.parameters['BEAM_CENTER_Y']
    dist = self._raw.parameters['DISTANCE']
    wavelength = self._raw.parameters['WAVELENGTH']
    r = math.sqrt((center_x - x_point)**2 + (center_y - y_point)**2)
    two_theta = math.atan(r / dist)
    d_min = wavelength / (2 * math.sin(two_theta / 2))
    return d_min

  def line_between_points (self, x1, y1, x2, y2, n_values=100, zoom=1) :
    x1_, y1_ = self.image_coords_as_array_coords(x1, y1, zoom)
    x2_, y2_ = self.image_coords_as_array_coords(x2, y2, zoom)
    print x1_, y1_, x2_, y2_
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
      dxdy = (y_2 - y_1) * (x_2 - x_1)
      vxy = ((v11 / dxdy) * (x_2 - x) * (y_2 - y)) + \
            ((v21 / dxdy) * (x - x_1) * (y_2 - y)) + \
            ((v12 / dxdy) * (x_2 - x) * (y - y_1)) + \
            ((v22 / dxdy) * (x - x_1) * (y - y_1))
      vals.append(vxy)
    return vals

class settings (object) :
  def __init__ (self) :
    self.zoom_level = 0
    self.brightness = 100
    self.show_beam_center = True
