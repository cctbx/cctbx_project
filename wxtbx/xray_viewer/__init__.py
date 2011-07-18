
from libtbx.utils import Sorry
import math

class image (object) :
  def __init__ (self, file_name) :
    from iotbx.detectors import ImageFactory
    try :
      from labelit.detectors import FlexImage # FIXME
    except ImportError, e :
      raise Sorry("Labelit not installed or not configured.")
    try :
      import Image
    except ImportError, e :
      raise Sorry("Python Imaging Library not installed - you can get it "+
        "at http://www.pythonware.com/products/pil/.")
    img = ImageFactory(file_name)
    img.read()
    print img.show_header()
    fi = FlexImage(
      rawdata=img.linearintdata,
      binning=2,
      vendortype=img.vendortype,
      brightness=1.0)
    fi.setWindow(0.0, 0.0, 1)
    fi.adjust()
    fi.prep_string()
    data_string = fi.export_string
    self._img = fi
    x = fi.ex_size2()
    y = fi.ex_size1()
    pil_image = Image.fromstring("RGB", (x,y),  data_string)
    self._fi = fi
    self._raw = img
    # XXX this doesn't work - why?  would be nice to remove the Labelit
    # dependency...
    #a = img.linearintdata.as_numpy_array()
    #pil_image = Image.fromstring("L", (a.shape[1], a.shape[0]), a.tostring())
    self._img = pil_image
    self._bmp = None

  def convert_to_bitmap (self) :
    import wx
    wx_image = wx.EmptyImage(*(self._img.size))
    wx_image.SetData(self._img.convert("RGB").tostring())
    bmp = wx_image.ConvertToBitmap() # wx.BitmapFromImage(image)
    self._bmp = bmp

  def get_bitmap (self) :
    if (self._bmp is None) :
      self.convert_to_bitmap()
    return self._bmp

  def get_size (self) :
    if (self._bmp is not None) :
      return self._bmp.GetSize()
    else :
      return self._img.size

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
