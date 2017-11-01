#!/usr/bin/env python
# FormatCBFMiniPilatusDLS12M.py
#
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, for the P12M-DLS

from __future__ import absolute_import, division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model import ParallaxCorrectedPxMmStrategy

class FormatCBFMiniPilatusDLS12M(FormatCBFMiniPilatus):

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 120-0100' in header:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    import libtbx
    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    # if multi_panel == False, then interpret data as 24 panels, where each row
    # of 5 panels is grouped as one "panel"
    # elif multi_panel == True, then interpret data as 120 panels,
    # 24 rows * 5 columns
    self._dynamic_shadowing = kwargs.get('dynamic_shadowing', False)
    if self._dynamic_shadowing in (libtbx.Auto, 'Auto'):
      self._dynamic_shadowing = True
    self._multi_panel = kwargs.get('multi_panel', False)
    FormatCBFMiniPilatus.__init__(self, image_file, **kwargs)

    self._raw_data = None


    return

  def _detector(self):

    # module positions from detector blueprints - modelling at the moment as
    # 24 modules, each consisting of 5 sensors (the latter is ignored)

    from dxtbx.model import Detector
    from scitbx import matrix
    import math

    x = matrix.col((-1, 0, 0))
    y = matrix.col((0, 1, 0))
    z = matrix.col((0, 0, 1))

    beam_xy = self._cif_header_dictionary['Beam_xy']
    beam_xy = beam_xy.replace('(', '').replace(')', '').replace(',', '').split()[:2]
    obs_beam_x, obs_beam_y = [float(f) for f in beam_xy]

    ideal_beam_x = 1075
    ideal_beam_y = 2594

    beam_shift_x = 0.172 * (ideal_beam_x - obs_beam_x)
    beam_shift_y = 0.172 * (ideal_beam_y - obs_beam_y)

    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0]) * 1000.0

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[2]) * 1000.0

    off_x = 184.9

    detector = Detector()
    root = detector.hierarchy()
    root.set_frame(
      x.elems,
      y.elems,
      (-distance * z + (beam_shift_x * x) + (beam_shift_y * y)).elems)

    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness
    px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

    self.coords = {}

    for j in range(24):
      shift_y = 195 + 17
      ymin, ymax = j * shift_y, j * shift_y + 195

      angle = math.pi * (-12.2 + 0.5 * 7.903 + j * (7.903 + 0.441)) / 180.0
      fast = matrix.col((1, 0, 0))
      slow = matrix.col((0, math.sin(angle), math.cos(angle)))
      normal = fast.cross(slow)

      row_origin = 250.0 * normal - off_x * fast - 16.8 * slow

      if not self._multi_panel:
        xmin, xmax = 0, 2463

        # OK two calls to add_panel here for detector like things => two
        # copies of the panel then? https://github.com/dials/dials/issues/189
        # ... this is also not the source of the leak

        # OBS! you need to set the panel to a root before set local frame...
        p = root.add_panel()
        p.set_type('SENSOR_PAD')
        p.set_name('row-%02d' % j)
        p.set_raw_image_offset((xmin, ymin))
        p.set_image_size((2463, 195))
        p.set_trusted_range((-1, 1000000))
        p.set_pixel_size((0.172, 0.172))
        p.set_local_frame(
          fast.elems,
          slow.elems,
          row_origin.elems)
        p.set_thickness(thickness)
        p.set_material('Si')
        p.set_mu(mu)
        p.set_px_mm_strategy(px_mm)
        p.set_raw_image_offset((xmin,ymin))
        self.coords[p.get_name()] = (xmin,ymin,xmax,ymax)

      else:
        shift_x = 487 + 7

        for i in range(5):
          xmin, xmax = i * shift_x, i * shift_x + 487
          origin = row_origin + i * (487+7) * 0.172 * fast

          # OBS! you need to set the panel to a root before set local frame...
          p = root.add_panel()
          p.set_type('SENSOR_PAD')
          p.set_name('row-%02d-col-%02d' % (j, i))
          p.set_raw_image_offset((xmin, ymin))
          p.set_image_size((487, 195))
          p.set_trusted_range((-1, 1000000))
          p.set_pixel_size((0.172, 0.172))
          p.set_local_frame(
            fast.elems,
            slow.elems,
            origin.elems)
          p.set_thickness(thickness)
          p.set_material('Si')
          p.set_mu(mu)
          p.set_px_mm_strategy(px_mm)
          p.set_raw_image_offset((xmin,ymin))
          self.coords[p.get_name()] = (xmin,ymin,xmax,ymax)

    return detector

  def read_cbf_image(self, cbf_image):
    from cbflib_adaptbx import uncompress
    import binascii

    start_tag = binascii.unhexlify('0c1a04d5')

    data = self.open_file(cbf_image, 'rb').read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[:data_offset - 4]

    fast = 0
    slow = 0
    length = 0

    for record in cbf_header.split('\n'):
      if 'X-Binary-Size-Fastest-Dimension' in record:
        fast = int(record.split()[-1])
      elif 'X-Binary-Size-Second-Dimension' in record:
        slow = int(record.split()[-1])
      elif 'X-Binary-Number-of-Elements' in record:
        length = int(record.split()[-1])
      elif 'X-Binary-Size:' in record:
        size = int(record.split()[-1])

    assert(length == fast * slow)

    pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                              fast = fast, slow = slow)

    return pixel_values

  def get_raw_data(self):
    if self._raw_data is None:
      raw_data = self.read_cbf_image(self._image_file)
      self._raw_data = []

      for panel in self.get_detector():
        xmin, ymin = panel.get_raw_image_offset()
        xmax = xmin + panel.get_image_size()[0]
        ymax = ymin + panel.get_image_size()[1]
        self._raw_data.append(raw_data[ymin:ymax,xmin:xmax])

    return tuple(self._raw_data)

  def get_goniometer_shadow_masker(self, goniometer=None):
    from dials.util.masking import GoniometerShadowMaskGenerator
    from scitbx.array_family import flex
    import math

    coords = flex.vec3_double((
      (0,0,0),
    ))

    alpha = flex.double_range(0, 190, step=10) * math.pi / 180
    r = flex.double(alpha.size(), 40)
    x = flex.double(r.size(), 107.61)
    y = -r*flex.sin(alpha)
    z = -r*flex.cos(alpha)
    coords.extend(flex.vec3_double(x, y, z))

    coords.extend(flex.vec3_double((
      # fixed
      (107.49, 7.84, 39.49),
      (107.39, 15.69, 38.97),
      (107.27, 23.53, 38.46),
      (107.16, 31.37, 37.94),
      (101.76, 33.99, 36.25),
      (96.37, 36.63, 34.56),
      (90.98, 39.25, 33.00),
      (85.58, 41.88, 31.18),
      (80.89, 47.06, 31.00),
      (76.55, 51.51, 31.03),
      (72.90, 55.04, 31.18),
      (66.86, 60.46, 31.67),
      (62.10, 64.41, 32.25),
    )))

    alpha = flex.double_range(180, 370, step=10) * math.pi / 180
    r = flex.double(alpha.size(), 33)
    x = (flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02**2) * flex.cos((50 * math.pi/180) - flex.atan(r/89.02 * flex.sin(alpha))))
    y = (flex.sqrt(flex.pow2(r * flex.sin(alpha)) + 89.02**2) * flex.sin((50 * math.pi/180) - flex.atan(r/89.02 * flex.sin(alpha))))
    z = -r*flex.cos(alpha)
    coords.extend(flex.vec3_double(x, y, z))

    coords.extend(flex.vec3_double((
      # fixed
      (62.10, 64.41, -32.25),
      (66.86, 60.46, -31.67),
      (72.90, 55.04, -31.18),
      (76.55, 51.51, -31.03),
      (80.89, 47.06, -31.00),
      (85.58, 41.88, -31.18),
      (90.98, 39.25, -33.00),
      (96.37, 36.63, -34.56),
      (101.76, 33.99, -36.25),
      (107.16, 31.37, -37.94),
      (107.27, 23.53, -38.46),
      (107.39, 15.69, -38.97),
      (107.49, 7.84, -39.49),
      (107.61, 0.00, -40.00)
    )))

    # I23 end station coordinate system:
    #   X-axis: positive direction is facing away from the storage ring (from
    #           sample towards goniometer)
    #   Y-axis: positive direction is vertically up
    #   Z-axis: positive direction is in the direction of the beam (from
    #           sample towards detector)
    #   K-axis (kappa): at an angle of +50 degrees from the X-axis
    #   K & phi rotation axes: clockwise rotation is positive (right hand
    #           thumb rule)
    #   Omega-axis: along the X-axis; clockwise rotation is positive

    # End station x-axis is parallel to ImgCIF x-axis
    # End station z-axis points in opposite direction to ImgCIF definition
    # (ImgCIF: The Z-axis is derived from the source axis which goes from
    # the sample to the source)
    # Consequently end station y-axis (to complete set following right hand
    # rule) points in opposite direction to ImgCIF y-axis.
    # Kappa arm aligned with -y in ImgCIF convention

    from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
    from scitbx import matrix
    R = align_reference_frame(matrix.col((1,0,0)), matrix.col((1,0,0)),
                              matrix.col((0,1,0)), matrix.col((0,-1,0)))
    coords = R.elems * coords

    if goniometer is None:
      goniometer = self.get_goniometer()
    return GoniometerShadowMaskGenerator(
      goniometer, coords, flex.size_t(len(coords), 1))

  def get_mask(self, goniometer=None):
    from dxtbx.model import MultiAxisGoniometer
    mask = super(FormatCBFMiniPilatusDLS12M, self).get_mask()
    if (isinstance(self.get_goniometer(), MultiAxisGoniometer) and
        self._dynamic_shadowing):
      gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
      scan = self.get_scan()
      detector = self.get_detector()
      shadow_mask = gonio_masker.get_mask(detector, scan.get_oscillation()[0])
      assert len(mask) == len(shadow_mask)
      for m, sm in zip(mask, shadow_mask):
        if sm is not None:
          m &= sm
    return mask

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header.'''

    from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
    timestamp = get_pilatus_timestamp(
        self._cif_header_dictionary['timestamp'])
    # Goniometer changed from reverse phi to conventional rotation direction
    # on this date:
    # calendar.timegm(time.strptime('2016-04-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    if timestamp < 1459468800:
      return self._goniometer_factory.single_axis_reverse()

    alpha = 50.0
    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])
    else:
      phi_value = 0.0

    if 'Kappa' in self._cif_header_dictionary:
      kappa_value = float(self._cif_header_dictionary['Kappa'].split()[0])
    else:
      kappa_value = 0.0

    if 'Omega' in self._cif_header_dictionary:
      omega_value = float(self._cif_header_dictionary['Omega'].split()[0])
    else:
      omega_value = 0.0

    return self._goniometer_factory.make_kappa_goniometer(
      alpha, omega_value, kappa_value, phi_value, '-y', 'omega')


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusDLS12M.understand(arg)
