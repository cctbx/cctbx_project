from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from iotbx.detectors import generic_flex_image

GenericFlexImage = generic_flex_image

class GenericDetector(object):
  def __init__(self,filename):
    self.filename = filename
    self.size2 = 200
    self.size1 = 250
    self.pixel_size = 0.1
    self.vendortype = "npy_raw"
    self.beamx = 0.
    self.beamy = 0.

  def readHeader(self):
    return
    from six.moves import cPickle as pickle
    G = open(self.filename,"rb")
    phil_stuff=pickle.load(G)
    data_stuff=pickle.load(G)

  def show_header(self):
    return "Generic detector with nothing in it"

  def read(self): pass

  def get_flex_image(self,brightness,**kwargs):
    # no kwargs supported at present
    rawdata = flex.random_double(200*250)
    rawdata.reshape(flex.grid(250,200))
    self.data = rawdata
    return GenericFlexImage(
      rawdata=rawdata,
      size1_readout=250,
      size2_readout=200,
      brightness=brightness,
      saturation=256.0)

  def initialize_viewer_properties(self,dummy_params):
    self.image_size_fast = self.size2
    self.image_size_slow = self.size1
    self.pixel_resolution = self.pixel_size

  def detector_coords_as_image_coords_float (self, x, y) :
    """
    Convert absolute detector position (in mm) to floating-value image pixel coordinates.
    """
    return x * self.pixel_resolution, \
           y * self.pixel_resolution

  def detector_coords_as_image_coords (self, x, y) :
    """
    Convert absolute detector position (in mm) to integer-value image pixel coordinates.
    """
    x_point,y_point = self.detector_coords_as_image_coords_float(x,y)
    return (int(x_point), int(y_point))

  def image_coords_as_detector_coords (self, x, y, readout=None) :
    """
    Convert image pixel coordinates to absolute position on the detector
    (in mm).
    """
    x_detector = x / self.pixel_resolution
    y_detector = y / self.pixel_resolution
    return x_detector, y_detector

  def get_beam_center_mm (self) :
    center_x = self.beamx
    center_y = self.beamy
    return center_x, center_y

  def get_beam_center_pixels_fast_slow(self):
    center_x, center_y = self.get_beam_center_mm()
    return self.detector_coords_as_image_coords_float(center_x, center_y)

  def get_pixel_intensity(self,coords):
    try:
      return self.data[(int(round(coords[0],0)), (int(round(coords[1],0))))]
    except IndexError:
      return None
