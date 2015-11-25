from __future__ import division
from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatStill import FormatStill

class FormatHDF5SaclaMPCCD(FormatHDF5, FormatStill):
  '''
  Class to handle multi-event HDF5 files from MPCCD
  preprocessed by Cheetah SFX pipeline at SACLA.

  To handle reassembled images from "DataConvert3 -reconst"
  (old pipeline), use FormatHDF5Sacla.
  '''

  @staticmethod
  def understand(image_file):
    import h5py
    h5_handle = h5py.File(image_file, 'r')

    for elem in h5_handle:
      if elem.startswith("tag-"):
        return True

    return False

  def __init__(self, image_file):
    assert(self.understand(image_file))
    FormatHDF5.__init__(self, image_file)
    self._raw_data = None

    self.pixel_size = 50 / 1000 # 50 um
    self.RECONST_SIZE = 2398 # compatible with DataConvert3 -reconst mode

    self.RECONST_MODE = True
    self.RECONST_64 = True

    self.distance = 52.0 # mm

    # TODO: These should be read from geometry file.
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

  def _start(self):
    import h5py
    self._h5_handle = h5py.File(self.get_image_file(), 'r')

    self._images = [tag for tag in self._h5_handle if tag.startswith("tag-")]
    self.tag = self._images[0]

  def _detector(self, index=None):
    from dxtbx.model.detector import HierarchicalDetector
    from scitbx import matrix
    import math

    # Ignore this for now.
    # TODO: must consider when a new MPCCD with 300um thickness
    #       is released.

    self.wavelength = 1.77 # A, this must be got from _beam
    #    thickness = 60 # um
    #
    #    from cctbx.eltbx import attenuation_coefficient
    #    table = attenuation_coefficient.get_table("Si")
    #    mu = table.mu_at_angstrom(wavelength) / 10.0
    #    t0 = thickness
    #    px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

    detector = HierarchicalDetector()
    root = detector.hierarchy()
    root.set_frame(
      (1, 0, 0),
      (0, 1, 0),
      (0, 0, - self.distance))

    if self.RECONST_MODE:
      return self._detector_factory.simple(
        sensor = 'PAD',
        distance = self.distance,
        beam_centre = (self.RECONST_SIZE * self.pixel_size / 2,
                       self.RECONST_SIZE * self.pixel_size / 2),
        fast_direction = '+x',
        slow_direction = '-y',
        pixel_size = (self.pixel_size,
                      self.pixel_size),
        image_size = (self.RECONST_SIZE,
                      self.RECONST_SIZE),
        trusted_range = (-1, 1000000),
        mask = [])  # a list of dead rectangles

    for i in range(8):
      angle = math.pi * self.panel_rotations[i] / 180.0
      fast = matrix.col((math.cos(angle), math.sin(angle), 0))
      slow = matrix.col((-math.sin(angle), math.cos(angle), 0))
      normal = fast.cross(slow)

      origin = matrix.col((-self.panel_origins[i][0],
                            self.panel_origins[i][1],
                            self.panel_origins[i][2])) / 1000.0

      p = detector.add_panel()

      # OBS! you need to set the panel to a root before set local frame...
      root.add_panel(p)
      p.set_name('panel-%01d' % i)
      p.set_image_size((512, 1024))
      p.set_trusted_range((-1, 1000000))
      p.set_pixel_size((self.pixel_size, self.pixel_size))
      p.set_local_frame(
        fast.elems,
        slow.elems,
        origin.elems)
      # p.set_px_mm_strategy(px_mm)

    return detector

  def _beam(self, index=None):
    eV = self._h5_handle[self.tag]['photon_energy_ev'].value

    return self._beam_factory.simple(12398.4/eV)

  def get_num_images(self):
    return len(self._images)

  def get_raw_data(self, index=0):
    import numpy

    self.tag = self._images[index]

    if self._raw_data is None:
      from scitbx.array_family import flex

      if self.RECONST_MODE:
        return flex.int(self.reconst_image())

      else:
        print "get_raw_data(%d) for %s" % (index, self.tag)
        data = self._h5_handle[self.tag]["data"][()].astype(numpy.int32)
        # [()] forces conversion to ndarray
        # this is 8192x512 (slow/fast) tiled image
        self._raw_data = []

        for i in range(8):
          xmin, ymin, xmax, ymax = 0, i * 1024, 512, (i + 1) * 1024
          # To avoid "numpy.ndarray instance is not contiguous"
          # TODO: Is this the right way?
          source = numpy.ascontiguousarray(data[ymin:ymax,xmin:xmax])
          self._raw_data.append(flex.int(source))

    if index is not None:
      return self._raw_data[index]

    return self._raw_data[0]

  def reconst_image(self):
    # TODO:
    # subdivide a panel into eight, otherwise coordinate error can be as big as 2px!

    import numpy
    from scitbx import matrix
    import math

    det = numpy.empty((self.RECONST_SIZE, self.RECONST_SIZE), dtype="int32") # was uint16
    det.fill(-1)

    data = self._h5_handle[self.tag]["data"][()].astype(numpy.int32)

    for i in range(8):
      angle = math.pi * self.panel_rotations[i] / 180.0
      fast = matrix.col((math.cos(angle), math.sin(angle)))
      slow = matrix.col((-math.sin(angle), math.cos(angle)))
      origin = matrix.col((-self.panel_origins[i][0] / self.pixel_size / 1000 + self.RECONST_SIZE / 2,
                           -self.panel_origins[i][1] / self.pixel_size / 1000 + self.RECONST_SIZE / 2))

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
            elif abs(round(self.panel_rotations[i]) - 90) < 1:
              det[round(subpanel_origin[1]):round(subpanel_origin[1] - size_slow):-1,
                  round(subpanel_origin[0]):round(subpanel_origin[0] - size_fast):-1] = source
            else:
              raise "Panel angle deviation is too large!"

      else:
        size_fast = 1024
        size_slow = 512

        xmin, ymin, xmax, ymax = 0, i * size_fast, size_slow, (i + 1) * size_fast
        source = data[ymin:ymax,xmin:xmax].transpose()

        print origin
        if abs(round(self.panel_rotations[i]) + 90) < 1:
          det[round(origin[1]):round(origin[1] + size_slow),
              round(origin[0]):round(origin[0] + size_fast)] = source
        elif abs(round(self.panel_rotations[i]) - 90) < 1:
          det[round(origin[1]):round(origin[1] - size_slow):-1,
              round(origin[0]):round(origin[0] - size_fast):-1] = source
        else:
          raise "Panel angle deviation is too large!"

    return det

  def get_image_file(self, index=None):
    return Format.get_image_file(self)

  def get_detector(self, index=None):
    if self._detector_instance is None:
      self._detector_instance = self._detector()

    return self._detector_instance

  def get_beam(self, index=None):
    if self._beam_instance is None:
      self._beam_instance = self._beam()

    return self._beam_instance

  def get_detectorbase(self, index=None):
    return None
