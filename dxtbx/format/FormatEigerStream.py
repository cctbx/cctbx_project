from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
from dxtbx.model import Beam # import dependency
from dxtbx.model import Detector # import dependency
from dxtbx.model import Goniometer # import dependency
from dxtbx.model import Scan # import dependency

injected_data = {}

class FormatEigerStream(FormatMultiImage, Format):
  '''
  A format class to understand an EIGER stream

  '''
  @staticmethod
  def understand(image_file):
    return bool(injected_data)

  def get_num_images(*args):
    return 1

  def __init__(self, image_file, **kwargs):
    from dxtbx import IncorrectFormatError
    if not injected_data:
      raise IncorrectFormatError(self, image_file)

    import json
    self.header = { 'configuration': json.loads(injected_data.get('header2', '')),
                    'info': json.loads(injected_data.get('streamfile_2', '')),
    }

    self._goniometer_instance = None
    self._detector_instance = None
    self._beam_instance = None
    self._scan_instance = None

    FormatMultiImage.__init__(self, **kwargs)
    Format.__init__(self, image_file, **kwargs)

    self.setup()

  def _detector(self):
    '''
    Create the detector model
    '''
    from scitbx import matrix
    from dxtbx.model.detector import DetectorFactory
    configuration = self.header['configuration']

    # Set the trusted range
#   trusted_range = 0, configuration['countrate_correction_count_cutoff']
    trusted_range = 0, 2 ** configuration['bit_depth_readout']

    # Get the sensor material and thickness
    sensor_material = str(configuration['sensor_material'])
    sensor_thickness = configuration['sensor_thickness']

    beam_center = configuration['beam_center_x'], configuration['beam_center_y']
    distance = configuration['detector_distance']

    # Get the pixel and image sizes
    pixel_size = (
      configuration['x_pixel_size'],
      configuration['y_pixel_size'])
    image_size = (
      configuration['x_pixels_in_detector'],
      configuration['y_pixels_in_detector'])

    # Get the detector axes
    # TODO Spec doesn't have detector_orientation
    # TODO THIS NEEDS FIXING
    #fast_axis = (1, 0, 0)
    #slow_axis = (0, 1, 0)
    #origin = (0, 0, -1)
    fast_axis = configuration['detector_orientation'][0:3]
    slow_axis = configuration['detector_orientation'][3:6]
    origin = matrix.col(configuration['detector_translation']) + matrix.col((0,0,-distance))

    # Create the detector model
    return DetectorFactory.make_detector(
      'SENSOR_PAD',
      fast_axis,
      slow_axis,
      origin,
      pixel_size,
      image_size,
      trusted_range,
      px_mm         = None,
      name          = 'Panel',
      thickness     = sensor_thickness,
      material      = sensor_material,
      mu            = 0.0)

  def _beam(self):
    '''
    Create the beam model

    '''
    from dxtbx.model.beam import BeamFactory
    configuration = self.header['configuration']
    return BeamFactory.simple(configuration['wavelength'])

  def _goniometer(self):
    '''
    Create the goniometer model

    '''
    from dxtbx.model.goniometer import GoniometerFactory
    return GoniometerFactory.single_axis()

  def _scan(self):
    '''
    Create the scan object
    '''
    from dxtbx.model.scan import ScanFactory
    configuration = self.header['configuration']
#   kappa_start = configuration['kappa_start']
#   kappa_increment = configuration['kappa_increment']
#   phi_start = configuration['phi_start']
#   phi_increment = configuration['phi_increment']
#   omega_start = configuration['omega_start']
#   omega_increment = configuration['omega_increment']
#   two_theta_start = configuration['two_theta_start']
#   two_theta_increment = configuration['two_theta_increment']
    kappa_start = 0
    kappa_increment = 0
    phi_start = 0
    phi_increment = 0
    omega_start = 0
    omega_increment = 0
    two_theta_start = 0
    two_theta_increment = 0
#   nimages = configuration['nimages']
    nimages = 1
    return ScanFactory.make_scan(
      image_range    = (1, nimages),
      exposure_times = [0] * nimages,
      oscillation    = (phi_start, phi_increment),
      epochs         = [0] * nimages)

  def get_raw_data(self, index):
    '''
    Get the raw data from the image
    '''
    import numpy as np
    from scitbx.array_family import flex

    info = self.header['info']
    data = injected_data['streamfile_3']
    if info["encoding"] == "lz4<":
      data = self.readLZ4(data, info["shape"], info["type"], info['size'])
    elif info["encoding"] == "bs16-lz4<":
      data = self.readBS16LZ4(data, info["shape"], info["type"], info['size'])
    elif info["encoding"] == "bs32-lz4<":
      data = self.readBSLZ4(data, info["shape"], info["type"], info['size'])
    else:
      raise IOError("encoding %s is not implemented" % info["encoding"])

    data = np.array(data,ndmin=3) # handle data, must be 3 dim
    data = data.reshape(data.shape[1:3]).astype("int32")
    return flex.int(data)

  def readBSLZ4(self, data, shape, dtype, size):
    """
    Unpack bitshuffle-lz4 compressed frame and return np array image data
    """
    import numpy as np
    import bitshuffle
    blob = np.fromstring(data[12:], dtype=np.uint8)
    # blocksize is big endian uint32 starting at byte 8, divided by element size
    blocksize = np.ndarray(shape=(), dtype=">u4", buffer=data[8:12])/4
    imgData = bitshuffle.decompress_lz4(blob, shape[::-1], np.dtype(dtype), blocksize)
    return imgData

  def readBS16LZ4(self, data, shape, dtype, size):
    """
    Unpack bitshuffle-lz4 compressed 16 bit frame and return np array image data
    """
    import bitshuffle
    import numpy
    blob = numpy.fromstring(data[12:], dtype=numpy.uint8)
    return bitshuffle.decompress_lz4(blob, shape[::-1], numpy.dtype(dtype))

  def readLZ4(self, data, shape, dtype, size):
    """
    Unpack lz4 compressed frame and return np array image data
    """
    import numpy as np
    import lz4

    dtype = np.dtype(dtype)
    data = lz4.loads(data)

    return np.reshape(np.fromstring(data, dtype=dtype), shape[::-1])

  def get_vendortype(self):
    from dxtbx.format.FormatPilatusHelpers import get_vendortype_eiger
    return get_vendortype_eiger(self.get_detector())
