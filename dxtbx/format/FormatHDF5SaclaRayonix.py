from __future__ import absolute_import, division, print_function
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatStill import FormatStill

class FormatHDF5SaclaRayonix(FormatHDF5, FormatStill):
  '''
  Class to handle multi-event HDF5 files from Rayonix detector (MX300-HS)
  preprocessed by Cheetah SFX pipeline at SACLA.
  '''

  @staticmethod
  def understand(image_file):
    import h5py

    h5_handle = h5py.File(image_file, 'r')
    if h5_handle['metadata/detector'].value == 'Rayonix MX300HS':
      return True
    return False

  def __init__(self, image_file, index = 0, reconst_mode = False, **kwargs):
    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)
    self._raw_data = None
    self.index = index
    self.image_filename = image_file
    FormatHDF5.__init__(self, image_file, **kwargs)

    self.PIXEL_SIZE = 78.2 / 1000 # um
    self.RECONST_SIZE = 3840
    # This hard-coded value can be overwritten
    # by RAYONIX_DISTANCE

    self.distance = 240.0 # mm
    self.mask = None

    # Read metadata if possible
    self.read_metadata()
    # Override by environmental variables
    import os

    if 'RAYONIX_DISTANCE' in os.environ:
      self.distance = float(os.environ['RAYONIX_DISTANCE'])

  def _start(self):
    import h5py
    h5_handle = h5py.File(self.image_filename, 'r')

    self._images = sorted([tag for tag in h5_handle if tag.startswith("tag-")])
    self.tag = self._images[self.index]
    h5_handle.close()

  def read_metadata(self):
    pass

  def get_image_file(self, index=None):
    return self.image_filename

  def set_index(self, index):
    assert(index < len(self._images))

    self.index = index
    self.tag = self._images[self.index]
    self._raw_data = None

  def _detector(self, index=None):
    wavelength = self.get_beam(index).get_wavelength()

    return self._detector_factory.simple(
      sensor = 'CCD',
      distance = self.distance,
      beam_centre = (self.RECONST_SIZE / 2 * self.PIXEL_SIZE,
                     self.RECONST_SIZE / 2 * self.PIXEL_SIZE),
      fast_direction = '-x',
      slow_direction = '-y',
      pixel_size = (self.PIXEL_SIZE,
                    self.PIXEL_SIZE),
      image_size = (self.RECONST_SIZE,
                    self.RECONST_SIZE),
      trusted_range = (-1, 65535),
      #px_mm = px_mm,
      mask = [])

  def _beam(self):
    import h5py
    h5_handle = h5py.File(self.image_filename, 'r')
    eV = h5_handle[self.tag]['photon_energy_ev'].value
    h5_handle.close()

    return self._beam_factory.simple(12398.4 / eV)

  def get_num_images(self):
    return len(self._images)

  def get_raw_data(self, index=None):
    if index is not None and self.index != index:
      self.set_index(index)

    if self._raw_data is None:
      from scitbx.array_family import flex
      import numpy
      import h5py

      h5_handle = h5py.File(self.image_filename, 'r')
      data = h5_handle[self.tag]["data"][()].astype(numpy.int32)
      h5_handle.close()
      self._raw_data = flex.int(data)

    return self._raw_data

  def get_detector(self, index=None):
    if self._detector_instance is None:
      self._detector_instance = self._detector()

    return self._detector_instance

  def get_mask(self, index=None, goniometer=None):
    # This means when the pixel mask is present, trusted region is ignored.
    # The used provided masks (if any) will be automatically merged.
    # see https://github.com/dials/dials/issues/236
    return self.mask

  def get_beam(self, index=None):
    if index is not None and self.index != index:
      self.set_index(index)
      self._beam_instance = None

    if self._beam_instance is None:
      self._beam_instance = self._beam()

    return self._beam_instance

if __name__ == '__main__':
  import sys
  print(FormatHDF5SaclaRayonix.understand(sys.argv[1]))
  FormatHDF5SaclaRayonix(sys.argv[1])
