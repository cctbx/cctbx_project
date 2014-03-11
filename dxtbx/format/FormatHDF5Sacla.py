from __future__ import division
from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatStill import FormatStill

class FormatHDF5Sacla(FormatHDF5, FormatStill):

  @staticmethod
  def understand(image_file):
    try:
      tag = FormatHDF5.open_file(image_file, 'rb').read(8)
    except IOError, e:
      return False

    # sorry all for the moment this is causing me problems trying to work out
    # the NeXus HDF5 format...

    if True:
      return False

    return tag == "\211HDF\r\n\032\n"

  def __init__(self, image_file):
    assert(self.understand(image_file))
    FormatHDF5.__init__(self, image_file)

  def _start(self):
    import h5py
    self._h5_handle = h5py.File(self.get_image_file())

    file_info = self._h5_handle['file_info']
    run_number_list = file_info['run_number_list']

    run_str = "run_%d"%run_number_list[0]

    self._run = self._h5_handle[run_str]

    event_info = self._run['event_info']
    tag_number_list = event_info['tag_number_list']
    self._images = ["tag_%d"%tag for tag in tag_number_list]

  def _detector(self, index=None):
    from scitbx import matrix

    # Get the pixel and image size
    detector_2d_assembled_1 = self._run['detector_2d_assembled_1']
    detector_info = detector_2d_assembled_1['detector_info']
    pixel_size = detector_info['pixel_size_in_micro_meter'][0]/1000, \
                 detector_info['pixel_size_in_micro_meter'][1]/1000
    tag = detector_2d_assembled_1[self._images[0]]
    data = tag['detector_data']
    image_size = data.shape
    trusted_range = (0, 200000)

    # Initialise detector frame
    fast = matrix.col((1.0, 0.0, 0.0))
    slow = matrix.col((0.0, -1.0, 0.0))
    orig = matrix.col((-image_size[0]*pixel_size[0]/2,
                        image_size[1]*pixel_size[1]/2, -100.0))


    # Make the detector
    return self._detector_factory.make_detector(
      "", fast, slow, orig,
      pixel_size, image_size, trusted_range)


  def _beam(self, index=None):
    run_info = self._run['run_info']
    sacla_config = run_info['sacla_config']
    eV = sacla_config['photon_energy_in_eV'].value

    return self._beam_factory.simple(12398.4/eV)

  def get_num_images(self):
    return len(self._images)

  def get_raw_data(self, index=0):
    from scitbx.array_family import flex
    detector_2d_assembled_1 = self._run['detector_2d_assembled_1']
    tag = detector_2d_assembled_1[self._images[index]]
    data = tag['detector_data']
    return flex.float(data[:,:]).as_double()

  def get_detectorbase(self, index=None):
    #self.detectorbase.img_number = index
    return Format.get_detectorbase(self)

  def get_image_file(self, index=None):
    return Format.get_image_file(self)

  def get_detector(self, index=None):
    return self._detector_instance

  def get_beam(self, index=None):
    return self._beam_instance
