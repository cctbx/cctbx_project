#!/usr/bin/env cctbx.python
#
# Biostruct-X Data Reduction Use Case 1.2:
#
# Given UB matrix, centring operation, generate a list of predictions as
# H K L x y phi. Also requires (clearly) a model for the detector positions
# and the crystal lattice type. This is aimed to help with identifying
# locations on the images.
#
# Requires:
#
# Determine maximum resolution limit.
# Generate full list of reflections to given resolution limit.
# Compute intersection angles for all reflections given UB matrix etc.
# Determine which of those will be recorded on the detector.

import sys
import math

from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
from rstbx.diffraction import rotation_angles, reflection_prediction
from rstbx.diffraction import full_sphere_indices
from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.uctbx import unit_cell
from rstbx.bpcx.detector_model.instrument_specifics import detector_factory_from_cfc

def Py_generate_indices(unit_cell_constants, resolution_limit):
    '''Generate all possible reflection indices out to a given resolution
    limit, ignoring symmetry and centring.'''

    uc = unit_cell(unit_cell_constants)

    maxh, maxk, maxl = uc.max_miller_indices(resolution_limit)

    indices = []

    for h in range(-maxh, maxh + 1):
        for k in range(-maxk, maxk + 1):
            for l in range(-maxl, maxl + 1):

                if h == 0 and k == 0 and l == 0:
                    continue

                if uc.d((h, k, l)) < resolution_limit:
                    continue

                indices.append((h, k, l))

    return indices

def Py_remove_absent_indices(indices, space_group_number):
    '''From the given list of indices, remove those reflections which should
    be systematic absences according to the given space group.'''

    sg = space_group(space_group_symbols(space_group_number).hall())

    present = []

    for hkl in indices:
        if not sg.is_sys_absent(hkl):
            present.append(hkl)

    return present

def parse_xds_xparm_scan_info(xparm_file):
    '''Read an XDS XPARM file, get the scan information.'''

    values = map(float, open(xparm_file).read().split())

    assert(len(values) == 42)

    img_start = values[0]
    osc_start = values[1]
    osc_range = values[2]

    return img_start, osc_start, osc_range

class python_reflection_prediction:
    def __init__(self, axis, s0, ub, detector_origin,
                 detector_fast, detector_slow,
                 f_min, f_max, s_min, s_max):
        self._axis = axis
        self._s0 = s0
        self._ub = ub
        self._detector_origin = detector_origin
        self._detector_fast = detector_fast
        self._detector_slow = detector_slow
        self._limits = f_min, f_max, s_min, s_max

        return

    def predict(self, indices, angles):

        detector_normal = self._detector_fast.cross(self._detector_slow)
        distance = self._detector_origin.dot(detector_normal)

        observed_reflection_positions = []

        for hkl, angle in zip(indices, angles):
            s = (self._ub * hkl).rotate(self._axis, angle)
            q = (s + self._s0).normalize()

            # check if diffracted ray parallel to detector face

            q_dot_n = q.dot(detector_normal)

            if q_dot_n == 0:
                continue

            r = (q * distance / q_dot_n) - self._detector_origin

            x = r.dot(self._detector_fast)
            y = r.dot(self._detector_slow)

            if x < self._limits[0] or y < self._limits[2]:
                continue
            if x > self._limits[1] or y > self._limits[3]:
                continue

            observed_reflection_positions.append((hkl, x, y, angle))

        return observed_reflection_positions

class make_prediction_list:

  def __init__(self, configuration_file, img_range, dmin = None,
               rocking_curve = "none", mosaicity_deg = 0.0):
      self._configuration_file = configuration_file
      self._img_range = img_range
      self._dmin = dmin
      self._rocking_curve = rocking_curve
      self._mosaicity_deg = mosaicity_deg
      return

  def predict_observations(self):
    '''Actually perform the prediction calculations.'''

    d2r = math.pi / 180.0
    cfc = coordinate_frame_converter(self._configuration_file)

    self.img_start, self.osc_start, self.osc_range = parse_xds_xparm_scan_info(
        self._configuration_file)

    if self._dmin is None:
        self._dmin = cfc.derive_detector_highest_resolution()

    phi_start = ((self._img_range[0] - self.img_start) * self.osc_range + \
                 self.osc_start) * d2r
    phi_end = ((self._img_range[1] - self.img_start + 1) * self.osc_range + \
               self.osc_start) * d2r
    self.phi_start_rad = phi_start
    self.phi_end_rad = phi_end
    # in principle this should come from the crystal model - should that
    # crystal model record the cell parameters or derive them from the
    # axis directions?

    A = cfc.get_c('real_space_a')
    B = cfc.get_c('real_space_b')
    C = cfc.get_c('real_space_c')

    cell = (A.length(), B.length(), C.length(), B.angle(C, deg = True),
            C.angle(A, deg = True), A.angle(B, deg = True))
    self.uc = unit_cell(cell)

    # generate all of the possible indices, then pull out those which should
    # be systematically absent

    sg = cfc.get('space_group_number')

    indices = full_sphere_indices(
      unit_cell = self.uc,
      resolution_limit = self._dmin,
      space_group = space_group(space_group_symbols(sg).hall()))

    # then get the UB matrix according to the Rossmann convention which
    # is used within the Labelit code.

    u, b = cfc.get_u_b(convention = cfc.ROSSMANN)
    axis = cfc.get('rotation_axis', convention = cfc.ROSSMANN)
    ub = u * b

    wavelength = cfc.get('wavelength')
    self.wavelength = wavelength

    # work out which reflections should be observed (i.e. pass through the
    # Ewald sphere)
    ra = rotation_angles(self._dmin, ub, wavelength, axis)

    obs_indices, obs_angles = ra.observed_indices_and_angles_from_angle_range(
        phi_start_rad = phi_start,
        phi_end_rad = phi_end,
        indices = indices)

    # convert all of these to full scattering vectors in a laboratory frame
    # (for which I will use the CBF coordinate frame) and calculate which
    # will intersect with the detector

    u, b = cfc.get_u_b()
    axis = cfc.get_c('rotation_axis')
    # must guarantee that sample_to_source vector is normalized so that
    #  s0 has length of 1/wavelength.
    sample_to_source_vec = cfc.get_c('sample_to_source').normalize()
    s0 = (- 1.0 / wavelength) * sample_to_source_vec
    ub = u * b

    # need some detector properties for this as well... starting to
    # abstract these to a detector model.
    df = detector_factory_from_cfc(cfc)
    d = df.build()

    # the Use Case assumes the detector consists of a single sensor
    sensor = d.sensors()[0]

    self.pixel_size_fast, self.pixel_size_slow = d.px_size_fast(), \
        d.px_size_slow()

    # used for polarization correction
    self.distance = sensor.distance

    rp = reflection_prediction(axis, s0, ub, sensor.origin,
                                   sensor.dir1,
                                   sensor.dir2,
                                   sensor.lim1[0], sensor.lim1[1],
                                   sensor.lim2[0], sensor.lim2[1])
    if self._rocking_curve is not None:
      assert self._rocking_curve != "none"
      rp.set_rocking_curve(self._rocking_curve)
      rp.set_mosaicity(self._mosaicity_deg, degrees = True)
    return rp.predict(obs_indices, obs_angles)

def test(configuration_file, img_range, dmin = None):
    '''Perform the calculations needed for use case 1.1.'''

    mp = make_prediction_list(configuration_file, img_range, dmin, None)
    obs_hkl,obs_fast,obs_slow,obs_angle = mp.predict_observations()

    r2d = 180.0 / math.pi

    for iobs in xrange(len(obs_hkl)):
      hkl = obs_hkl[iobs]
      f = obs_fast[iobs]
      s = obs_slow[iobs]
      angle = obs_angle[iobs]
#FIXME this is not a class, should not see 'self'
      print '%5d %5d %5d' % hkl, '%11.4f %11.4f %9.2f' % (
            f / self.pixel_size_fast, s / self.pixel_size_slow,
            (self.img_start - 1) + ((angle * r2d) - self.osc_start) / \
          self.osc_range)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        msg = "Requires 3 arguments: path/to/xparm.xds start_image_no end_image_no"
        sys.exit(msg)
    if len(sys.argv) == 4:
        test(sys.argv[1], (int(sys.argv[2]), int(sys.argv[3])))
    else:
        test(sys.argv[1], (int(sys.argv[2]), int(sys.argv[3])),
            float(sys.argv[4]))
