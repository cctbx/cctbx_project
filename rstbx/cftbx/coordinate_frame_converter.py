import math
import sys

from scitbx import matrix
from cctbx import uctbx

from coordinate_frame_helpers import is_xds_xparm, import_xds_xparm
from coordinate_frame_helpers import is_xds_integrate_hkl, \
    import_xds_integrate_hkl

class coordinate_frame_converter:
    '''A class which is instantiated from a supported file (initially an
    imgCIF image or an XDS XPARM / INTEGRATE.HKL / XDS_ASCII.HKL file) and
    will make available the rotation axis, beam vector, detector position
    and attitude, detector origin, fast and slow directions and so on in
    a range of different program specific coordinate frames.'''

    CBF = 'CBF'
    ROSSMANN = 'Rossmann'
    MOSFLM = 'Mosflm'

    def __init__(self, configuration_file):
        '''Construct a coordinate frame converter from a configuration file.'''

        if is_xds_xparm(configuration_file):
            self._coordinate_frame_information = import_xds_xparm(
                configuration_file)

        elif is_xds_integrate_hkl(configuration_file):
            self._coordinate_frame_information = import_xds_integrate_hkl(
                configuration_file)

        else:
            raise RuntimeError, 'unknown configuration file %s' % \
                  configuration_file

        return

    def get(self, parameter, convention = CBF):
        '''Get a parameter, in a given reference frame if a vector quantity,
        as a Python basic type.'''

        parameter_value = self._coordinate_frame_information.get(parameter)

        if not hasattr(parameter_value, 'elems'):
            return parameter_value

        if convention == coordinate_frame_converter.CBF:
            R = self._coordinate_frame_information.R_to_CBF()
            return (R * parameter_value).elems
        elif convention == coordinate_frame_converter.ROSSMANN:
            R = self._coordinate_frame_information.R_to_Rossmann()
            return (R * parameter_value).elems
        elif convention == coordinate_frame_converter.MOSFLM:
            R = self._coordinate_frame_information.R_to_Mosflm()
            return (R * parameter_value).elems
        else:
            raise RuntimeError, 'convention %s not currently supported' % \
                  convention

        return

    def get_c(self, parameter, convention = CBF):
        '''Get the parameter, in the correct coordinate convention if a
        vector, as a cctbx matrix.col or a floating point value.'''

        parameter_value = self._coordinate_frame_information.get(parameter)

        if not hasattr(parameter_value, 'elems'):
            return parameter_value

        if convention == coordinate_frame_converter.CBF:
            R = self._coordinate_frame_information.R_to_CBF()
            return R * parameter_value
        elif convention == coordinate_frame_converter.ROSSMANN:
            R = self._coordinate_frame_information.R_to_Rossmann()
            return R * parameter_value
        elif convention == coordinate_frame_converter.MOSFLM:
            R = self._coordinate_frame_information.R_to_Mosflm()
            return R * parameter_value
        else:
            raise RuntimeError, 'convention %s not currently supported' % \
                  convention

        return

    def get_u_b(self, convention = CBF):
        '''Get the [U] and [B] matrices in the requested coordinate system.'''

        cfi = self._coordinate_frame_information

        if not cfi.get_real_space_a() or not cfi.get_real_space_b() or \
           not cfi.get_real_space_c():
            raise RuntimeError, 'orientation matrix information missing'

        axis_a = cfi.get_real_space_a()
        axis_b = cfi.get_real_space_b()
        axis_c = cfi.get_real_space_c()

        A = matrix.sqr(axis_a.elems +  axis_b.elems + axis_c.elems).inverse()

        a = axis_a.length()
        b = axis_b.length()
        c = axis_c.length()

        alpha = axis_b.angle(axis_c, deg = True)
        beta = axis_c.angle(axis_a, deg = True)
        gamma = axis_a.angle(axis_b, deg = True)

        uc = uctbx.unit_cell((a, b, c, alpha, beta, gamma))

        B = matrix.sqr(uc.fractionalization_matrix()).transpose()

        U = A * B.inverse()

        if convention == coordinate_frame_converter.CBF:
            R = cfi.R_to_CBF()
        elif convention == coordinate_frame_converter.ROSSMANN:
            R = cfi.R_to_Rossmann()
        elif convention == coordinate_frame_converter.MOSFLM:
            R = cfi.R_to_Mosflm()
        else:
            raise RuntimeError, 'convention %s not currently supported' % \
                  convention

        return R * U, B

    def get_unit_cell(self):
        '''Get the unit cell.'''

        cfi = self._coordinate_frame_information

        if not cfi.get_real_space_a() or not cfi.get_real_space_b() or \
           not cfi.get_real_space_c():
            raise RuntimeError, 'orientation matrix information missing'

        axis_a = cfi.get_real_space_a()
        axis_b = cfi.get_real_space_b()
        axis_c = cfi.get_real_space_c()

        A = matrix.sqr(axis_a.elems +  axis_b.elems + axis_c.elems).inverse()

        a = axis_a.length()
        b = axis_b.length()
        c = axis_c.length()

        alpha = axis_b.angle(axis_c, deg = True)
        beta = axis_c.angle(axis_a, deg = True)
        gamma = axis_a.angle(axis_b, deg = True)

        uc = uctbx.unit_cell((a, b, c, alpha, beta, gamma))

        return uc

    def derive_beam_centre_pixels_fast_slow(self):
        '''Derive the pixel position at which the direct beam would intersect
        with the detector plane, and return this in terms of fast and slow.'''

        cfi = self._coordinate_frame_information

        detector_origin = cfi.get('detector_origin')
        detector_fast = cfi.get('detector_fast')
        detector_slow = cfi.get('detector_slow')
        sample_to_source = cfi.get('sample_to_source')
        pixel_size_fast, pixel_size_slow = cfi.get(
            'detector_pixel_size_fast_slow')

        detector_normal = detector_fast.cross(detector_slow)

        if not sample_to_source.dot(detector_normal):
            raise RuntimeError, 'beam parallel to detector'

        distance = detector_origin.dot(detector_normal)

        sample_to_detector = sample_to_source * distance / \
                             sample_to_source.dot(detector_normal)

        beam_centre = sample_to_detector - detector_origin

        beam_centre_fast_mm = beam_centre.dot(detector_fast)
        beam_centre_slow_mm = beam_centre.dot(detector_slow)

        return beam_centre_fast_mm / pixel_size_fast, \
               beam_centre_slow_mm / pixel_size_slow

    def derive_detector_highest_resolution(self):
        '''Determine the highest resolution recorded on the detector, which
        should be at one of the corners.'''

        cfi = self._coordinate_frame_information

        detector_origin = cfi.get('detector_origin')
        detector_fast = cfi.get('detector_fast')
        detector_slow = cfi.get('detector_slow')
        sample_to_source = cfi.get('sample_to_source')
        pixel_size_fast, pixel_size_slow = cfi.get(
            'detector_pixel_size_fast_slow')
        size_fast, size_slow = cfi.get(
            'detector_size_fast_slow')

        F = detector_origin + detector_fast * pixel_size_fast * size_fast
        S = detector_origin + detector_slow * pixel_size_slow * size_slow
        FS = F + S - detector_origin

        detector_normal = detector_fast.cross(detector_slow)
        distance = detector_origin.dot(detector_normal)

        sample_to_detector = sample_to_source * distance / \
                             sample_to_source.dot(detector_normal)

        theta = 0.5 * max([sample_to_detector.angle(detector_origin),
                           sample_to_detector.angle(F),
                           sample_to_detector.angle(S),
                           sample_to_detector.angle(FS)])

        return cfi.get_wavelength() / (2.0 * math.sin(theta))

    def __str__(self):

        return '\n'.join([
            'detector origin: %.3f %.3f %.3f' % self.get('detector_origin'),
            'detector fast: %6.3f %6.3f %6.3f' % self.get('detector_fast'),
            'detector slow: %6.3f %6.3f %6.3f' % self.get('detector_slow'),
            'rotation axis: %6.3f %6.3f %6.3f' % self.get('rotation_axis'),
            '- s0 vector:   %6.3f %6.3f %6.3f' % self.get('sample_to_source')
            ])

if __name__ == '__main__':

    if len(sys.argv) < 2:
        raise RuntimeError, '%s configuration-file mosflm-matrix' % sys.argv[0]

    configuration_file = sys.argv[1]

    cfc = coordinate_frame_converter(configuration_file)

    print 'Maximum resolution: %.2f' % cfc.derive_detector_highest_resolution()

    mosflm_matrix = matrix.sqr(
        map(float, open(sys.argv[2]).read().split()[:9]))

    u, b = cfc.get_u_b(convention = cfc.MOSFLM)

    wavelength = cfc.get('wavelength')

    mosflm_matrix = (1.0 / wavelength) * mosflm_matrix

    matrix_format = '%8.5f %8.5f %8.5f\n%8.5f %8.5f %8.5f\n%8.5f %8.5f %8.5f'

    print matrix_format % mosflm_matrix.elems
    print matrix_format % (u * b).elems
