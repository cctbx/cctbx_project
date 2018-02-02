from __future__ import absolute_import, division
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatStill import FormatStill

# 151028: deepcopying this class causes crash in h5py
#         temporary fix by closing the file in every methods(!)
# 161003: updated to follow dxtbx changes
#         removed iotbx support, which was incomplete anyway
# 161005: get wavelength from the file

class FormatHDF5SaclaMPCCD(FormatHDF5, FormatStill):
  '''
  Class to handle multi-event HDF5 files from MPCCD
  preprocessed by Cheetah SFX pipeline at SACLA.
  To handle reassembled images from "DataConvert3 -reconst"
  (old pipeline), use FormatHDF5Sacla.
  To override metrology, use the following environmental variables.
     MPCCD_GEOMETRY, MPCCD_DISTANCE
     MPCCD_RECONST_MODE

  You can also specify reference_geometry in dials.stills_process.
  '''

  @staticmethod
  def understand(image_file):
    import h5py

    h5_handle = h5py.File(image_file, 'r')

    for elem in h5_handle:
      if elem.startswith("tag-"):
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

    self.PIXEL_SIZE = 50 / 1000 # 50 um
    self.RECONST_SIZE = 2398 # compatible with DataConvert3 -reconst mode

    # These hard-coded values can be overwritten
    # by MPCCD_GEOMETRY and MPCCD_DISTANCE
    #
    # These values can be retrieved from SACLA API.
    # Alternatively, you can get it from a CrystFEL geometry file by
    # awk '/corner_x/{x=50*$3} /corner_y/{y=50*$3; printf x","y","rot","}
    #      /\/ss/{rot=-atan2($3, $4)/3.141592*180}' input.geom

    self.distance = 50.0 # mm
    self.panel_origins = [(-1755.000000, 51711.000000, 0.000000),
                          (-1711.000000, 24944.000000, 0.000000),
                          (817.000000, -1808.000000, 0.000000),
                          (812.000000, -28466.000000, 0.000000),
                          (-792.000000, 28544.000000, 0.000000),
                          (-781.000000, 1840.000000, 0.000000),
                          (1650.000000, -24900.000000, 0.000000),
                          (1655.000000, -51626.000000, 0.000000)] # um
    self.panel_rotations = [-89.906197, -89.915802, -89.980003, -89.929298,
                            89.963097, 89.880798, 90.000000, 90.029503]

    import os

    if 'MPCCD_RECONST_MODE' in os.environ:
      reconst_mode = bool(os.environ['MPCCD_RECONST_MODE'])
    self.RECONST_MODE = reconst_mode
    self.RECONST_64 = True # Set False if you want to keep panels completely horizontal
                           # But this makes errors bigger.

    if 'MPCCD_GEOMETRY' in os.environ:
      try:
        tmp = map(float, os.environ['MPCCD_GEOMETRY'].split(","))
        if len(tmp) != 24:
          raise
        for i in range(8):
          self.panel_origins[i] = (-tmp[i * 3], tmp[i * 3 + 1], 0)
          self.panel_rotations[i] = tmp[i * 3 + 2]
      except Exception:
        raise "Invalid MPCCD Geomtry"
    if 'MPCCD_DISTANCE' in os.environ:
      self.distance = float(os.environ['MPCCD_DISTANCE'])

  def _start(self):
    import h5py
    h5_handle = h5py.File(self.image_filename, 'r')

    self._images = sorted([tag for tag in h5_handle if tag.startswith("tag-")])
    self.tag = self._images[self.index]
    h5_handle.close()

  def get_image_file(self, index=None):
    return self.image_filename

  def set_index(self, index):
    assert(index < len(self._images))

    self.index = index
    self.tag = self._images[self.index]
    self._raw_data = None

  def _detector(self, index=None):
    from dxtbx.model.detector import Detector
    from scitbx import matrix
    import math

    wavelength = self.get_beam(index).get_wavelength()

    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    t0 = 0.050  # sensor thickness in um
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

    if self.RECONST_MODE:
      return self._detector_factory.simple(
        sensor = 'PAD',
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
        px_mm = px_mm,
        mask = [])  # TODO: add gaps

    detector = Detector()
    root = detector.hierarchy()
    root.set_frame(
      (-1, 0, 0),
      ( 0, 1, 0),
      ( 0, 0, - self.distance))

    for i in range(8):
      angle = math.pi * self.panel_rotations[i] / 180.0
      fast = matrix.col((math.cos(angle), math.sin(angle), 0))
      slow = matrix.col((-math.sin(angle), math.cos(angle), 0))
      normal = fast.cross(slow)

      origin = matrix.col((-self.panel_origins[i][0],
                            self.panel_origins[i][1],
                            self.panel_origins[i][2])) / 1000.0
      p = root.add_panel()
      p.set_type("SENSOR_PAD")
      p.set_name('Panel%d' % i)
      p.set_image_size((512, 1024))
      p.set_trusted_range((-1, 65535))
      p.set_pixel_size((self.PIXEL_SIZE, self.PIXEL_SIZE))
      p.set_thickness(t0)
      p.set_local_frame(
        fast.elems,
        slow.elems,
        origin.elems)
      p.set_px_mm_strategy(px_mm)
      p.set_gain(10)

    return detector

  def _beam(self):
    import h5py
    h5_handle = h5py.File(self.image_filename, 'r')
    eV = h5_handle[self.tag]['photon_energy_ev'].value
    h5_handle.close()

    return self._beam_factory.simple(12398.4 / eV)

  def get_num_images(self):
    return len(self._images)

  def get_raw_data(self, index=None):
    import numpy

    if index is not None and self.index != index:
      self.set_index(index)

    if self._raw_data is None:
      from scitbx.array_family import flex

      if self.RECONST_MODE:
        self._raw_data = flex.int(self.reconst_image())

      else:
        import h5py
        h5_handle = h5py.File(self.image_filename, 'r')

        data = h5_handle[self.tag]["data"][()].astype(numpy.int32)
        # [()] forces conversion to ndarray
        # this is 8192x512 (slow/fast) tiled image
        tmp = []

        h5_handle.close()

        for i in range(8):
          xmin, ymin, xmax, ymax = 0, i * 1024, 512, (i + 1) * 1024
          # To avoid "numpy.ndarray instance is not contiguous"
          # TODO: Is this the right way?
          source = numpy.ascontiguousarray(data[ymin:ymax,xmin:xmax])
          tmp.append(flex.int(source))

        self._raw_data = tuple(tmp)

    return self._raw_data

  def get_active_areas(self):
    assert(self.RECONST_MODE)

    return self.active_areas

  def reconst_image(self):
    import numpy
    from scitbx import matrix
    import math
    import h5py

    det = numpy.empty((self.RECONST_SIZE, self.RECONST_SIZE), dtype="int32")
    det.fill(-1)

    h5_handle = h5py.File(self.image_filename, 'r')
    data = h5_handle[self.tag]["data"][()].astype(numpy.int32)
    h5_handle.close()

    self.active_areas = []

    for i in range(8):
      angle = math.pi * self.panel_rotations[i] / 180.0
      fast = matrix.col((math.cos(angle), math.sin(angle)))
      slow = matrix.col((-math.sin(angle), math.cos(angle)))
      origin = matrix.col((-self.panel_origins[i][0] / self.PIXEL_SIZE / 1000 + self.RECONST_SIZE / 2,
                           -self.panel_origins[i][1] / self.PIXEL_SIZE / 1000 + self.RECONST_SIZE / 2))

      if self.RECONST_64:
        size_fast = 256
        size_slow = 256

        for j in range(2):
          for k in range(4):
            xmin, ymin, xmax, ymax = j * size_slow, (i * 4 + k) * size_fast, (j + 1) * size_slow, (i * 4 + k + 1) * size_fast
            source = data[ymin:ymax,xmin:xmax].transpose()

            subpanel_origin = origin - j * size_fast * fast + k * size_slow * slow
            if abs(round(self.panel_rotations[i]) + 90) < 1:
              det[round(subpanel_origin[1]):round(subpanel_origin[1] + size_slow),
                  round(subpanel_origin[0]):round(subpanel_origin[0] + size_fast)] = source
              # TODO: Is the border inclusive?
              self.active_areas.extend([round(subpanel_origin[1]), round(subpanel_origin[0]),
                                        round(subpanel_origin[1] + size_slow), round(subpanel_origin[0] + size_fast)])
            elif abs(round(self.panel_rotations[i]) - 90) < 1:
              det[round(subpanel_origin[1]):round(subpanel_origin[1] - size_slow):-1,
                  round(subpanel_origin[0]):round(subpanel_origin[0] - size_fast):-1] = source
              self.active_areas.extend([round(subpanel_origin[1] - size_slow), round(subpanel_origin[0] - size_fast),
                                        round(subpanel_origin[1]),round(subpanel_origin[0])])
            else:
              raise "Panel angle deviation is too large! Do not use reconst mode!"

      else:
        size_fast = 1024
        size_slow = 512

        xmin, ymin, xmax, ymax = 0, i * size_fast, size_slow, (i + 1) * size_fast
        source = data[ymin:ymax,xmin:xmax].transpose()

        if abs(round(self.panel_rotations[i]) + 90) < 1:
          det[round(origin[1]):round(origin[1] + size_slow),
              round(origin[0]):round(origin[0] + size_fast)] = source
        elif abs(round(self.panel_rotations[i]) - 90) < 1:
          det[round(origin[1]):round(origin[1] - size_slow):-1,
              round(origin[0]):round(origin[0] - size_fast):-1] = source
        else:
          raise "Panel angle deviation is too large! Do not use reconst mode!"

    self.active_areas = map(int, self.active_areas)
    return det

  def get_detector(self, index=None):
    if self._detector_instance is None:
      self._detector_instance = self._detector()

    return self._detector_instance

  def get_beam(self, index=None):
    if index is not None and self.index != index:
      self.set_index(index)
      self._beam_instance = None

    if self._beam_instance is None:
      self._beam_instance = self._beam()

    return self._beam_instance

if __name__ == '__main__':
  import sys
  print FormatHDF5SaclaMPCCD.understand(sys.argv[1])
  FormatHDF5SaclaMPCCD(sys.argv[1])
