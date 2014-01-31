from __future__ import division
from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5

class FormatHDF5Nexus(FormatHDF5):
  
  @staticmethod
  def understand(image_file):
    try:
      tag = FormatHDF5.open_file(image_file, 'rb').read(8)
    except IOError, e:
      return False

    return tag == "\211HDF\r\n\032\n"

  def __init__(self, image_file):
    assert(self.understand(image_file))
    FormatHDF5.__init__(self, image_file)

  def _start(self):
    import h5py 
    self._h5_handle = h5py.File(self.get_image_file())

  def _goniometer(self):
    ''' Get the rotation axis. '''
    entry = self._h5_handle['entry']
    sample = entry['sample']
    pose = sample['pose']
    axis = pose['CBF_axis_omega'].attrs['vector']
    return self._goniometer_factory.known_axis(axis) 

  def _detector(self):
    from scitbx import matrix

    # Get the detector geometry stuff
    entry = self._h5_handle['entry']
    instrument = entry['instrument']
    detector = instrument['detector']
    pose = detector['pose']
    translation = pose['translation']
    rotation = pose['rotation']
   
    # Get the translation
    offset = translation.attrs['offset']
    trans = translation[0]
    vector = matrix.col(translation.attrs['vector']).normalize()

    # Initialise detector frame
    fast = matrix.col((1.0, 0.0, 0.0))
    slow = matrix.col((0.0, 1.0, 0.0))
    orig = matrix.col((offset[0] + trans * vector[0],
                       offset[1] + trans * vector[1],
                       offset[2] + trans * vector[2]))

    # Next comes a rotation about an axis
    vector = matrix.col(rotation.attrs['vector']).normalize()
    angle = rotation[0]
    m_rot = vector.axis_and_angle_as_r3_rotation_matrix(angle, deg=True)

    # Transform detector frame
    fast = (m_rot * fast).normalize() 
    slow = (m_rot * slow).normalize()
    orig = m_rot * orig

    # Get the pixel and image size
    pixel_size = detector['x_pixel_size'].value, detector['y_pixel_size'].value
    image_size = len(detector['x_pixel_offset']), len(detector['y_pixel_offset'])
    trusted_range = (0, detector['saturation_value'][0])
   
    # Make the detector
    return self._detector_factory.make_detector(
      "", fast, slow, orig,
      pixel_size, image_size, trusted_range)


  def _beam(self):
    ''' Nexus defines beam along z axis. '''    
    entry = self._h5_handle['entry']
    sample = entry['sample']
    beam = sample['beam']
    wavelength = beam['wavelength']
    return self._beam_factory.simple(wavelength[0])

  def _scan(self):
    ''' Get the scan. '''
    import time
    entry = self._h5_handle['entry']
    sample = entry['sample']
    pose = sample['pose']
    angles = pose['CBF_axis_omega']
    oscillation = (angles[0], angles[1] - angles[0])
    image_range = (1, len(angles))
    instrument = entry['instrument']
    detector = instrument['detector']
    exposure_times = detector['count_time']
    
    # Create the epochs
    frame_time = detector['frame_time']
    start_time = entry['start_time']
    time_ssec = start_time.value.split('.')
    time_struct = time.strptime(time_ssec[0], "%Y-%m-%dT%H:%M:%S")
    start_time = time.mktime(time_struct) + float('0.%s' % time_ssec[1])
    epochs = {0 : start_time}
    for i, t in enumerate(frame_time[:-1]):
      epochs[i] = epochs[i-1] + f
    
    # Create the scan  
    return self._scan_factory.make_scan(
      image_range, 
      list(exposure_times), 
      oscillation, 
      list(epochs), 
      deg=True)

  def get_num_images(self):
    entry = self._h5_handle['entry']
    data = entry['data']['data']
    return data.shape[0]

  def get_goniometer(self, index=None):
    return Format.get_goniometer(self)

  def get_detector(self, index=None):
    return Format.get_detector(self)

  def get_beam(self, index=None):
    return Format.get_beam(self)

  def get_scan(self, index=None):
    return Format.get_scan(self)

  def get_raw_data(self, index):
    from scitbx.array_family import flex
    entry = self._h5_handle['entry']
    data = entry['data']['data']
    return flex.int(data[index,:,:])

  def get_image_file(self, index=None):
    return Format.get_image_file(self)

