
from libtbx.utils import Sorry
from libtbx.math_utils import ifloor, iceil
from libtbx.str_utils import format_value
import math

class screen_params (object) :
  def __init__ (self, img_w=None, img_h=None) :
    self.img_w = img_w
    self.img_h = img_h
    self.screen_w = None
    self.screen_h = None
    self.thumb_w = None
    self.thumb_h = None
    self.thumb_ratio = None
    self.img_x_offset = 0
    self.img_y_offset = 0
    self.screen_x_start = 0
    self.screen_y_start = 0
    self.detector_pixel_size = 0
    self.zoom = 0
    self._invert_x = False
    self._invert_y = False

  def set_zoom (self, zoom) :
    assert (zoom >= 0)
    self.zoom = zoom

  def set_screen_size (self, w, h) :
    self.screen_w = w
    self.screen_h = h
    scale = self.get_scale()
    w_img_pixels = int(w / scale)
    h_img_pixels = int(h / scale)
    if (w_img_pixels > (self.img_w - self.img_x_offset)) :
      self.img_x_offset = max(0, self.img_w - w_img_pixels)
    if (h_img_pixels > (self.img_h - self.img_y_offset)) :
      self.img_y_offset = max(0, self.img_h - w_img_pixels)

  def set_image_size (self, w, h) :
    self.img_w = w
    self.img_h = h

  def set_thumbnail_size (self, w, h, ratio) :
    self.thumb_w = w
    self.thumb_h = h
    self.thumb_ratio = ratio

  def get_image_size (self) :
    assert (not None in [self.img_w, self.img_h])
    return self.img_w, self.img_h

  def get_thumbnail_size (self) :
    assert (not None in [self.thumb_w, self.thumb_h])
    return (self.thumb_w, self.thumb_h)

  def set_detector_resolution (self, p) :
    self.detector_pixel_size = p

  def get_scale (self) :
    assert (not None in [self.screen_w, self.screen_h, self.img_w, self.img_h])
    if (self.zoom == 0) :
      return min(self.screen_w / self.img_w, self.screen_h / self.img_h)
    else :
      return self.zoom

  def get_detector_dimensions (self) :
    assert (not None in [self.img_w, self.img_h, self.detector_pixel_size])
    return (self.img_w * self.detector_pixel_size,
            self.img_h * self.detector_pixel_size)

  def adjust_screen_coordinates (self, x, y) :
    xi, yi, w, h = self.get_bitmap_params()
    scale = self.get_scale()
    x_ = x + max(0, (self.screen_w - (w*scale)) / 2)
    y_ = y + max(0, (self.screen_h - (h*scale)) / 2)
    return (int(x_), int(y_))

  def get_bitmap_params (self) :
    scale = self.get_scale()
    x0 = max(self.img_x_offset, 0)
    y0 = max(self.img_y_offset, 0)
    w_scaled = min(self.img_w, (self.screen_w / scale)) #- x0
    h_scaled = min(self.img_h, (self.screen_h / scale))
    return (x0, y0, w_scaled, h_scaled)

  def get_thumbnail_box (self) :
    x, y, w, h = self.get_bitmap_params()
    tr = self.thumb_ratio
    return int(x / tr), int(y / tr), int(w / tr), int(h / tr)

  def get_zoom_box (self, x, y, boxsize=400, mag=16) :
    assert ((boxsize % mag) == 0)
    n_pixels = boxsize / mag
    x0 = min(self.img_w - n_pixels, x - (n_pixels / 2))
    y0 = min(self.img_h - n_pixels, y - (n_pixels / 2))
    return (x0, y0, n_pixels, n_pixels)

  def translate_image (self, delta_x, delta_y) :
    scale = self.get_scale()
    x_new = max(0, ifloor(self.img_x_offset - (delta_x / scale)))
    y_new = max(0, ifloor(self.img_y_offset - (delta_y / scale)))
    max_x = ifloor(self.img_w - (self.screen_w / scale))
    max_y = ifloor(self.img_h - (self.screen_h / scale))
    self.img_x_offset = min(x_new, max_x)
    self.img_y_offset = min(y_new, max_y)

  def center_view_from_thumbnail (self, x, y) :
    if (self.zoom == 0) : return
    x0, y0, w, h = self.get_bitmap_params()
    img_x = max(0, ifloor((x * self.thumb_ratio) - (w / 2)))
    img_y = max(0, ifloor((y * self.thumb_ratio) - (h / 2)))
    scale = self.get_scale()
    max_x = ifloor(self.img_w - (self.screen_w / scale))
    max_y = ifloor(self.img_h - (self.screen_h / scale))
    self.img_x_offset = min(img_x, max_x)
    self.img_y_offset = min(img_y, max_y)

  def image_coords_as_screen_coords (self, x, y) :
    scale = self.get_scale()
    x1 = self.screen_x_start + ((x+0.5) - self.img_x_offset) * scale
    y1 = self.screen_y_start + ((y+0.5) - self.img_y_offset) * scale
    xi, yi, w, h = self.get_bitmap_params()
    x2 = x1 + max(0, (self.screen_w - (w*scale)) / 2)
    y2 = y1 + max(0, (self.screen_h - (h*scale)) / 2)
    return ((x2), (y2))

  def detector_coords_as_image_coords (self, x, y) :
    dw = self.img_w * self.detector_pixel_size
    dh = self.img_h * self.detector_pixel_size
    x_frac = x / dw
    if (self._invert_y) :
      y_frac = - ((y / dh) - 1.0)
    else :
      y_frac = y / dh
    x_point = x_frac * self.img_w
    y_point = y_frac * self.img_h
    return (int(x_point), int(y_point))

  def image_coords_as_detector_coords (self, x, y) :
    dw, dh = self.get_detector_dimensions()
    w, h = self.get_image_size()
    x_frac = x / w
    y_frac = y / h
    x_detector = x_frac * dw
    if (self._invert_y) :
      y_detector = (1.0 - y_frac) * dh
    else :
      y_detector = y_frac * dh
    return x_detector, y_detector

  def screen_coords_as_image_coords (self, x, y) :
    scale = self.get_scale()
    xi, yi, w, h = self.get_bitmap_params()
    x1 = x - max(0, (self.screen_w - (w*scale)) / 2)
    y1 = y - max(0, (self.screen_h - (h*scale)) / 2)
    x2 = self.img_x_offset + (x1 / scale)
    y2 = self.img_y_offset + (y1 / scale)
    return (ifloor(x2) - 1, ifloor(y2) - 1)

  def image_coords_as_array_coords (self, x, y) :
    return y-1, x-1

class image (screen_params) :
  def __init__ (self, file_name) :
    screen_params.__init__(self)
    self.file_name = file_name
    from iotbx.detectors import ImageFactory
    img = ImageFactory(file_name)
    img.read()
    self._raw = img
    print img.show_header()
    self._invert_beam_center = False
    self.set_image_size(
      w=self._raw.parameters['SIZE2'],
      h=self._raw.parameters['SIZE1'])
    self.set_detector_resolution(self._raw.parameters['PIXEL_SIZE'])
    from spotfinder.command_line.signal_strength import master_params
    from iotbx.detectors.context.config_detector import \
      beam_center_convention_from_image_object
    params = master_params.extract()
    bc = beam_center_convention_from_image_object(img,params)
    print "beam center convention: %d" % bc
    # FIXME what about 2-4 & 6-7?
    if (bc == 0) :
      self._invert_beam_center = True
      self._invert_y = True
    elif (bc == 1) :
      self._invert_y = False
    elif (bc == 5) :
      self._invert_y = True
    #self.update_image()
    #self.convert_to_bitmap()

  def create_flex_image (self,
                         brightness=100,
                         binning=1) :
    # FIXME
    try :
      from labelit.detectors import FlexImage
    except ImportError, e :
      raise Sorry("Labelit not installed or not configured.")
    saturation = getattr(self._raw, "saturation", 65535)
    fi = FlexImage(
      rawdata=self._raw.linearintdata,
      binning=binning,
      vendortype=self._raw.vendortype,
      brightness=brightness / 100.,
      saturation=int(saturation))
    #from scitbx.array_family import flex
    #print flex.max(self._raw.linearintdata), flex.min(self._raw.linearintdata)
    fi.setWindow(0.0, 0.0, 1)
    fi.adjust()
    fi.prep_string()
    return fi

  def update_settings (self, **kwds) :
    self._wx_img = None
    self.update_image(**kwds)

  def update_image (self, brightness=100) :
    import wx
    self._img = self.create_flex_image(brightness)
    w = self._img.ex_size2()
    h = self._img.ex_size1()
    self.set_image_size(w, h)
    wx_image = wx.EmptyImage(w, h)
    wx_image.SetData(self._img.export_string)
    self._wx_img = wx_image
    fi_thumb = self.create_flex_image(brightness=brightness,
      binning=8)
    w = fi_thumb.ex_size2()
    h = fi_thumb.ex_size1()
    wx_thumb = wx.EmptyImage(w, h)
    wx_thumb.SetData(fi_thumb.export_string)
    ratio = min(w / 256, h / 256)
    wt, ht = ifloor(w/ratio), ifloor(h / ratio)
    wx_thumb = wx_thumb.Rescale(wt, ht)
    self.set_thumbnail_size(wt, ht, ratio * 8)
    self._wx_thumb = wx_thumb
    self._wx_thumb_bmp = wx_thumb.ConvertToBitmap()

  def get_bitmap (self) :
    import wx
    x, y, w, h = self.get_bitmap_params()
    scale = self.get_scale()
    img = self._wx_img.GetSubImage((x, y, w, h))
    img = img.Scale(w * scale, h * scale, wx.IMAGE_QUALITY_NORMAL)
    return img.ConvertToBitmap()

  def get_thumbnail_bitmap (self) :
    return self._wx_thumb_bmp #.ConvertToBitmap()

  def get_zoomed_bitmap (self, x, y, boxsize=400, mag=16) :
    import wx
    x0, y0, w, h = self.get_zoom_box(x, y, boxsize, mag)
    assert (w == h)
    img = self._wx_img.GetSubImage((x0, y0, w, h))
    return img.Scale(boxsize, boxsize, wx.IMAGE_QUALITY_NORMAL)

  def get_beam_center (self) :
    # FIXME Pilatus and ADSC images appear to have different conventions???
    if (self._invert_beam_center) :
      center_x = self._raw.parameters['BEAM_CENTER_Y']
      center_y = self._raw.parameters['BEAM_CENTER_X']
    else :
      center_x = self._raw.parameters['BEAM_CENTER_X']
      center_y = self._raw.parameters['BEAM_CENTER_Y']
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
    n_values = ifloor(math.sqrt((x2_-x1_)**2 + (y2_-y1_)**2))
    delta_x = (x2_ - x1_) / (n_values - 1)
    delta_y = (y2_ - y1_) / (n_values - 1)
    vals = []
    d = self._raw.linearintdata
    # TODO remarkably, this is reasonably fast in Python, but it would
    # probably be more at home in scitbx.math
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
