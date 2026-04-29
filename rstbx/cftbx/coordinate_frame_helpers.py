from __future__ import absolute_import, division, print_function

import io
import math
import random

from scitbx import matrix
from cctbx import sgtbx
from six.moves import range, zip

class coordinate_frame_information:
    '''A bucket class to store coordinate frame information.'''

    def __init__(self, detector_origin, detector_fast, detector_slow,
                 detector_size_fast_slow, detector_pixel_size_fast_slow,
                 rotation_axis, sample_to_source, wavelength,
                 real_space_a=None, real_space_b=None,
                 real_space_c=None, space_group_number=None,
                 sigma_divergence=None,
                 mosaicity=None,
                 starting_angle=None, oscillation_range=None,
                 starting_frame=None,
                 original_rotation=None,
                 data_range=None,
                 panel_offset=None,
                 panel_size=None,
                 panel_origin=None,
                 panel_fast=None,
                 panel_slow=None):

        self._detector_origin = detector_origin
        self._detector_fast = detector_fast
        self._detector_slow = detector_slow
        self._detector_size_fast_slow = detector_size_fast_slow
        self._detector_pixel_size_fast_slow = detector_pixel_size_fast_slow
        self._rotation_axis = rotation_axis
        self._sample_to_source = sample_to_source
        self._wavelength = wavelength
        self._real_space_a = real_space_a
        self._real_space_b = real_space_b
        self._real_space_c = real_space_c
        self._space_group_number = space_group_number
        self._sigma_divergence = sigma_divergence
        self._mosaicity = mosaicity
        self._starting_angle = starting_angle
        self._oscillation_range = oscillation_range
        self._starting_frame = starting_frame
        self._data_range = data_range
        self._original_rotation = original_rotation
        self._panel_offset = panel_offset
        self._panel_size = panel_size
        self._panel_origin = panel_origin
        self._panel_fast = panel_fast
        self._panel_slow = panel_slow

        self._R_to_CBF = None
        self._R_to_Rossmann = None
        self._R_to_Mosflm = None

        return

    def get_detector_origin(self):
        return self._detector_origin

    def get_detector_fast(self):
        return self._detector_fast

    def get_detector_slow(self):
        return self._detector_slow

    def get_rotation_axis(self):
        return self._rotation_axis

    def get_sample_to_source(self):
        return self._sample_to_source

    def get_wavelength(self):
        return self._wavelength

    def get_real_space_a(self):
        return self._real_space_a

    def get_real_space_b(self):
        return self._real_space_b

    def get_real_space_c(self):
        return self._real_space_c

    def get_space_group_number(self):
        return self._space_group_number

    def get_original_rotation(self):
        return self._original_rotation

    def get(self, parameter_name):
        if not hasattr(self, '_%s' % parameter_name):
            raise RuntimeError('no parameter %s' % parameter_name)
        return getattr(self, '_%s' % parameter_name)

    def R_to_CBF(self):

        if not self._R_to_CBF:
            self._R_to_CBF = align_reference_frame(
                self._rotation_axis, (1.0, 0.0, 0.0),
                self._sample_to_source, (0.0, 0.0, 1.0))

        return self._R_to_CBF

    def R_to_Rossmann(self):

        if not self._R_to_Rossmann:
            self._R_to_Rossmann = align_reference_frame(
                self._sample_to_source, (0.0, 0.0, - 1.0),
                self._rotation_axis, (0.0, 1.0, 0.0))

        return self._R_to_Rossmann

    def R_to_Mosflm(self):

        if not self._R_to_Mosflm:
            self._R_to_Mosflm = align_reference_frame(
                self._sample_to_source, (- 1.0, 0.0, 0.0),
                self._rotation_axis, (0.0, 0.0, 1.0))

        return self._R_to_Mosflm

def orthogonal_component(reference, changing):
    '''Return unit vector corresponding to component of changing orthogonal to
    reference.'''

    r = reference.normalize()
    c = changing.normalize()

    return (c - c.dot(r) * r).normalize()

def align_reference_frame(primary_axis, primary_target,
                          secondary_axis, secondary_target):
    '''Compute a rotation matrix R: R x primary_axis = primary_target and
    R x secondary_axis places the secondary_axis in the plane perpendicular
    to the primary_target, as close as possible to the secondary_target.
    Require: primary_target orthogonal to secondary_target, primary axis
    not colinear with secondary axis.'''

    if type(primary_axis) == type(()) or type(primary_axis) == type([]):
        primary_axis = matrix.col(primary_axis).normalize()
    else:
        primary_axis = primary_axis.normalize()

    if type(primary_target) == type(()) or type(primary_target) == type([]):
        primary_target = matrix.col(primary_target).normalize()
    else:
        primary_target = primary_target.normalize()

    if type(secondary_axis) == type(()) or type(secondary_axis) == type([]):
        secondary_axis = matrix.col(secondary_axis).normalize()
    else:
        secondary_axis = secondary_axis.normalize()

    if type(secondary_target) == type(()) or \
           type(secondary_target) == type([]):
        secondary_target = matrix.col(secondary_target).normalize()
    else:
        secondary_target = secondary_target.normalize()

    # check properties of input axes

    assert math.fabs(primary_axis.angle(secondary_axis) % math.pi) > 0.001
    p_dot_s = primary_target.dot(secondary_target)
    assert p_dot_s < 0.001, p_dot_s

    if primary_target.angle(primary_axis) % math.pi:
      axis_p = primary_target.cross(primary_axis)
      angle_p = - primary_target.angle(primary_axis)
      Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    elif primary_target.dot(primary_axis) < 0:
      axis_p = primary_axis.ortho().normalize()
      angle_p = math.pi
      Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    else:
      Rprimary = matrix.identity(3)

    if math.fabs(secondary_target.angle(Rprimary * secondary_axis)) < 1.0e-6:
      Rsecondary = matrix.identity(3)
    else:
      axis_r = secondary_target.cross(Rprimary * secondary_axis)
      axis_s = primary_target

      if (axis_r.angle(primary_target, value_if_undefined=0) > 0.5 * math.pi):
        angle_s = orthogonal_component(axis_s, secondary_target).angle(
          orthogonal_component(axis_s, Rprimary * secondary_axis))
      else:
        angle_s = - orthogonal_component(axis_s, secondary_target).angle(
          orthogonal_component(axis_s, Rprimary * secondary_axis))

      Rsecondary = axis_s.axis_and_angle_as_r3_rotation_matrix(angle_s)

    return Rsecondary * Rprimary

def is_xds_inp(putative_xds_inp_file):
    '''See if this file looks like an XDS.INP file.'''
    from iotbx.xds import xds_inp
    return xds_inp.reader.is_xds_inp_file(putative_xds_inp_file)

def is_xds_xparm(putative_xds_xparm_file):
    '''See if this file looks like an XDS XPARM file i.e. it consists of 42
    floating point values and nothing else.'''
    from iotbx.xds import xparm
    return xparm.reader.is_xparm_file(putative_xds_xparm_file)

def is_xds_integrate_hkl(putative_integrate_hkl_file):
    '''See if this looks like an XDS INTEGRATE.HKL file.'''

    with io.open(putative_integrate_hkl_file, encoding="ascii") as fh:
      try:
        first_record = fh.readline()
        return '!OUTPUT_FILE=INTEGRATE.HKL' in first_record
      except UnicodeDecodeError:
        return False

def is_xds_ascii_hkl(putative_xds_ascii_hkl_file):
    '''See if this looks like an XDS INTEGRATE.HKL file.'''

    with io.open(putative_xds_ascii_hkl_file, encoding="ascii") as f:
      try:
        f.readline()
        return "!OUTPUT_FILE=XDS_ASCII.HKL" in f.readline()
      except UnicodeDecodeError:
        return False

def is_recognized_file(filename):
    ''' Check if the file is recognized.'''
    if is_xds_xparm(filename):
        return True
    elif is_xds_integrate_hkl(filename):
        return True
    elif is_xds_ascii_hkl(filename):
        return True
    elif is_xds_inp(filename):
        return True

    # Not recognized
    return False

def import_xds_integrate_hkl(integrate_hkl_file):
    '''Read an XDS INTEGRATE.HKL file, transform the parameters contained
    therein into the standard coordinate frame, record this as a dictionary.'''

    assert(is_xds_integrate_hkl(integrate_hkl_file))

    header = []

    with open(integrate_hkl_file) as fh:
        for record in fh:
            if not record.startswith('!'):
                break

            header.append(record)

    # now need to dig out the values I want, convert and return

    distance = None

    for record in header:
        if record.startswith('!ROTATION_AXIS='):
            axis = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!INCIDENT_BEAM_DIRECTION='):
            beam = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!DIRECTION_OF_DETECTOR_X-AXIS='):
            x = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!DIRECTION_OF_DETECTOR_Y-AXIS='):
            y = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_A-AXIS='):
            a = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_B-AXIS='):
            b = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_C-AXIS='):
            c = [float(r) for r in  record.split()[-3:]]
            continue
        if record.startswith('!X-RAY_WAVELENGTH='):
            wavelength = float(record.split()[-1])
            continue
        if record.startswith('!DETECTOR_DISTANCE='):
            distance = float(record.split()[-1])
            continue
        if record.startswith('!SPACE_GROUP_NUMBER='):
            space_group_number = int(record.split()[-1])
            continue
        if record.startswith('!BEAM_DIVERGENCE_E.S.D.'):
            sigma_divergence = float(record.split()[-1])
            continue
        if record.startswith('!REFLECTING_RANGE_E.S.D.'):
            mosaicity = float(record.split()[-1])
            continue
        if record.startswith('!NX='):
            nx = int(record.split()[1])
            ny = int(record.split()[3])
            px = float(record.split()[5])
            py = float(record.split()[7])
            continue
        if record.startswith('!ORGX='):
            ox = float(record.split()[1])
            oy = float(record.split()[3])
            try:
                distance = float(record.split()[5])
            except IndexError: # Older versions of this file do not contain distance on this line.
                pass
            continue
        if record.startswith('!STARTING_FRAME'):
            starting_frame = int(record.split()[-1])
            continue
        if record.startswith('!STARTING_ANGLE'):
            starting_angle = float(record.split()[-1])
            continue
        if record.startswith('!OSCILLATION_RANGE'):
            oscillation_range = float(record.split()[-1])
            continue

    # check parameters set
    assert not distance is None

    # XDS defines the beam vector as s0 rather than from sample -> source.
    # Keep in mind that any inversion of a vector needs to be made with great
    # care!

    B = - matrix.col(beam).normalize()
    A = matrix.col(axis).normalize()

    X = matrix.col(x).normalize()
    Y = matrix.col(y).normalize()
    N = X.cross(Y)

    _X = matrix.col([1, 0, 0])
    _Y = matrix.col([0, 1, 0])
    _Z = matrix.col([0, 0, 1])

    R = align_reference_frame(A, _X, B, _Z)

    # Need to subtract 0.5 because XDS seems to do centroids in fortran coords
    ox = ox - 0.5
    oy = oy - 0.5

    detector_origin = R * (distance * N - ox * px * X - oy * py * Y)
    detector_fast = R * X
    detector_slow = R * Y
    rotation_axis = R * A
    sample_to_source = R * B
    real_space_a = R * matrix.col(a)
    real_space_b = R * matrix.col(b)
    real_space_c = R * matrix.col(c)

    return coordinate_frame_information(
        detector_origin, detector_fast, detector_slow, (nx, ny), (px, py),
        rotation_axis, sample_to_source, wavelength,
        real_space_a, real_space_b, real_space_c, space_group_number,
        sigma_divergence, mosaicity,
        starting_angle, oscillation_range, starting_frame, original_rotation = R)

def import_xds_ascii_hkl(xds_ascii_hkl_file):
    '''Read an XDS INTEGRATE.HKL file, transform the parameters contained therein
    into the standard coordinate frame, record this as a dictionary.'''

    assert(is_xds_ascii_hkl(xds_ascii_hkl_file))

    header = []

    for record in open(xds_ascii_hkl_file):
        if not record.startswith('!'):
            break

        header.append(record)

    # now need to dig out the values I want, convert and return

    for record in header:
        if record.startswith('!ROTATION_AXIS='):
            axis = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!INCIDENT_BEAM_DIRECTION='):
            beam = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!DIRECTION_OF_DETECTOR_X-AXIS='):
            x = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!DIRECTION_OF_DETECTOR_Y-AXIS='):
            y = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_A-AXIS='):
            a = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_B-AXIS='):
            b = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!UNIT_CELL_C-AXIS='):
            c = [float(r) for r in record.split()[-3:]]
            continue
        if record.startswith('!X-RAY_WAVELENGTH='):
            wavelength = float(record.split()[-1])
            continue
        if record.startswith('!DETECTOR_DISTANCE='):
            distance = float(record.split()[-1])
            continue
        if record.startswith('!SPACE_GROUP_NUMBER='):
            space_group_number = int(record.split()[-1])
            continue
        if record.startswith('!BEAM_DIVERGENCE_E.S.D.'):
            sigma_divergence = float(record.split()[-1])
            continue
        if record.startswith('!REFLECTING_RANGE_E.S.D.'):
            mosaicity = float(record.split()[-1])
            continue
        if record.startswith('!NX='):
            nx = int(record.split()[1])
            ny = int(record.split()[3])
            px = float(record.split()[5])
            py = float(record.split()[7])
            continue
        if record.startswith('!ORGX='):
            ox = float(record.split()[1])
            oy = float(record.split()[3])
            continue
        if record.startswith('!STARTING_FRAME'):
            starting_frame = int(record.split()[-1])
            continue
        if record.startswith('!STARTING_ANGLE'):
            starting_angle = float(record.split()[-1])
            continue
        if record.startswith('!OSCILLATION_RANGE'):
            oscillation_range = float(record.split()[-1])
            continue

    # XDS defines the beam vector as s0 rather than from sample -> source.
    # Keep in mind that any inversion of a vector needs to be made with great
    # care!

    B = - matrix.col(beam).normalize()
    A = matrix.col(axis).normalize()

    X = matrix.col(x).normalize()
    Y = matrix.col(y).normalize()
    N = X.cross(Y)

    _X = matrix.col([1, 0, 0])
    _Y = matrix.col([0, 1, 0])
    _Z = matrix.col([0, 0, 1])

    R = align_reference_frame(A, _X, B, _Z)

    # Need to subtract 0.5 because XDS seems to do centroids in fortran coords
    ox = ox - 0.5
    oy = oy - 0.5


    detector_origin = R * (distance * N - ox * px * X - oy * py * Y)
    detector_fast = R * X
    detector_slow = R * Y
    rotation_axis = R * A
    sample_to_source = R * B
    real_space_a = R * matrix.col(a)
    real_space_b = R * matrix.col(b)
    real_space_c = R * matrix.col(c)

    return coordinate_frame_information(
        detector_origin, detector_fast, detector_slow, (nx, ny), (px, py),
        rotation_axis, sample_to_source, wavelength,
        real_space_a, real_space_b, real_space_c, space_group_number,
        sigma_divergence, mosaicity,
        starting_angle, oscillation_range, starting_frame, original_rotation = R)

def import_xds_inp(xds_inp_file):
    '''Read an XDS XPARM file, transform the parameters contained therein
    into the standard coordinate frame, record this as a dictionary.'''
    from iotbx.xds import xds_inp

    handle = xds_inp.reader()
    handle.read_file(xds_inp_file)

    # first determine the rotation R from the XDS coordinate frame used in
    # the processing to the central (i.e. imgCIF) coordinate frame. N.B.
    # if the scan was e.g. a PHI scan the resulting frame could well come out
    # a little odd...

    axis = handle.rotation_axis
    beam = handle.incident_beam_direction
    x, y = handle.direction_of_detector_x_axis, handle.direction_of_detector_y_axis

    # XDS defines the beam vector as s0 rather than from sample -> source.

    B = - matrix.col(beam).normalize()
    A = matrix.col(axis).normalize()

    X = matrix.col(x).normalize()
    Y = matrix.col(y).normalize()
    N = X.cross(Y)

    _X = matrix.col([1, 0, 0])
    _Y = matrix.col([0, 1, 0])
    _Z = matrix.col([0, 0, 1])

    R = align_reference_frame(A, _X, B, _Z)

    # now transform contents of the XPARM file to the form which we want to
    # return...

    nx, ny = handle.nx, handle.ny
    px, py = handle.px, handle.py

    distance = handle.detector_distance
    ox, oy = handle.orgx, handle.orgy

    # Need to subtract 0.5 because XDS seems to do centroids in fortran coords
    ox = ox - 0.5
    oy = oy - 0.5

    detector_origin = R * (distance * N - ox * px * X - oy * py * Y)
    detector_fast = R * X
    detector_slow = R * Y
    rotation_axis = R * A
    sample_to_source = R * B
    wavelength = handle.xray_wavelength
    real_space_a, real_space_b, real_space_c = None, None, None
    space_group_number = handle.space_group_number
    starting_angle = handle.starting_angle
    oscillation_range = handle.oscillation_range
    starting_frame = handle.starting_frame
    data_range = handle.data_range

    panel_offset = None
    panel_size = None
    panel_origins = None
    panel_fast_axes = None
    panel_slow_axes = None

    if handle.num_segments > 1:
        # Now, for each detector segment the following two lines of information
        # are provided.

        # The 5 numbers of this line, iseg x1 x2 y1 y2, define the pixel numbers
        # IX,IY belonging to segment #iseg as x1<=IX<=x2, y1<=IY<=y2.
        # The 9 numbers of this line, ORGXS ORGYS FS EDS(:,1) EDS(:,2), describe
        # origin and orientation of segment #iseg with respect to the detector
        # coordinate system.

        panel_offset = []
        panel_size = []
        panel_origins = []
        panel_fast_axes = []
        panel_slow_axes = []
        for i in range(handle.num_segments):
            x1, x2, y1, y2 = handle.segment[i]
            # XDS panel limits inclusive range
            panel_offset.append((x1-1, y1-1))
            panel_size.append((x2-x1+1, y2-y1+1))
            panel_fast = matrix.col(handle.direction_of_segment_x_axis[i])
            panel_slow = matrix.col(handle.direction_of_segment_y_axis[i])
            # local basis vectors
            fl = matrix.col(panel_fast)
            sl = matrix.col(panel_slow)
            nl = fl.cross(sl)

            orgxs = handle.segment_orgx[i]
            orgys = handle.segment_orgy[i]
            fs = handle.segment_distance[i]
            panel_origin = R * (- (orgxs - x1 + 1) * px * fl \
                                - (orgys - y1 + 1) * py * sl \
                                + fs * nl ) \
                + detector_origin

            # detector to laboratory transformation
            ED = matrix.sqr(list(X) + list(Y) + list(N))

            panel_normal = (R * panel_fast).cross(R * panel_slow)
            panel_origins.append(panel_origin)
            panel_fast_axes.append(R * ED * panel_fast)
            panel_slow_axes.append(R * ED * panel_slow)

    return coordinate_frame_information(
        detector_origin, detector_fast, detector_slow, (nx, ny), (px, py),
        rotation_axis, sample_to_source, wavelength,
        real_space_a, real_space_b, real_space_c, space_group_number,
        None, None, starting_angle, oscillation_range, starting_frame,
        original_rotation=R,
        data_range=data_range,
        panel_offset=panel_offset,
        panel_size=panel_size,
        panel_origin=panel_origins,
        panel_fast=panel_fast_axes,
        panel_slow=panel_slow_axes)


def import_xds_xparm(xparm_file):
    '''Read an XDS XPARM file, transform the parameters contained therein
    into the standard coordinate frame, record this as a dictionary.'''
    from iotbx.xds import xparm

    handle = xparm.reader()
    handle.read_file(xparm_file)

    # first determine the rotation R from the XDS coordinate frame used in
    # the processing to the central (i.e. imgCIF) coordinate frame. N.B.
    # if the scan was e.g. a PHI scan the resulting frame could well come out
    # a little odd...

    axis = handle.rotation_axis
    beam = handle.beam_vector
    x, y = handle.detector_x_axis, handle.detector_y_axis

    # XDS defines the beam vector as s0 rather than from sample -> source.

    B = - matrix.col(beam).normalize()
    A = matrix.col(axis).normalize()

    X = matrix.col(x).normalize()
    Y = matrix.col(y).normalize()
    N = X.cross(Y)

    _X = matrix.col([1, 0, 0])
    _Y = matrix.col([0, 1, 0])
    _Z = matrix.col([0, 0, 1])

    R = align_reference_frame(A, _X, B, _Z)

    # now transform contents of the XPARM file to the form which we want to
    # return...

    nx, ny = handle.detector_size
    px, py = handle.pixel_size

    distance = handle.detector_distance
    ox, oy = handle.detector_origin

    a = handle.unit_cell_a_axis
    b = handle.unit_cell_b_axis
    c = handle.unit_cell_c_axis

    # Need to subtract 0.5 because XDS seems to do centroids in fortran coords
    ox = ox - 0.5
    oy = oy - 0.5

    detector_origin = R * (distance * N - ox * px * X - oy * py * Y)
    detector_fast = R * X
    detector_slow = R * Y
    rotation_axis = R * A
    sample_to_source = R * B
    wavelength = handle.wavelength
    real_space_a = R * matrix.col(a)
    real_space_b = R * matrix.col(b)
    real_space_c = R * matrix.col(c)
    space_group_number = handle.space_group
    starting_angle = handle.starting_angle
    oscillation_range = handle.oscillation_range
    starting_frame = handle.starting_frame

    panel_offset = None
    panel_size = None
    panel_origins = None
    panel_fast_axes = None
    panel_slow_axes = None

    if handle.num_segments > 1:
        # Now, for each detector segment the following two lines of information
        # are provided.

        # The 5 numbers of this line, iseg x1 x2 y1 y2, define the pixel numbers
        # IX,IY belonging to segment #iseg as x1<=IX<=x2, y1<=IY<=y2.
        # The 9 numbers of this line, ORGXS ORGYS FS EDS(:,1) EDS(:,2), describe
        # origin and orientation of segment #iseg with respect to the detector
        # coordinate system.

        panel_offset = []
        panel_size = []
        panel_origins = []
        panel_fast_axes = []
        panel_slow_axes = []
        for segment, orientation,  in zip(handle.segments, handle.orientation):
            iseg, x1, x2, y1, y2 = segment
            # XDS panel limits inclusive range
            panel_offset.append((x1-1, y1-1))
            panel_size.append((x2-x1+1, y2-y1+1))
            panel_fast = matrix.col(orientation[3:6])
            panel_slow = matrix.col(orientation[6:9])
            # local basis vectors
            fl = matrix.col(panel_fast)
            sl = matrix.col(panel_slow)
            nl = fl.cross(sl)

            orgxs, orgys, fs = orientation[:3]
            panel_origin = R * (- (orgxs - x1 + 1) * px * fl \
                                - (orgys - y1 + 1) * py * sl \
                                + fs * nl ) \
                + detector_origin

            # detector to laboratory transformation
            ED = matrix.sqr(list(X) + list(Y) + list(N))

            panel_normal = (R * panel_fast).cross(R * panel_slow)
            panel_origins.append(panel_origin)
            panel_fast_axes.append(R * ED * panel_fast)
            panel_slow_axes.append(R * ED * panel_slow)

    return coordinate_frame_information(
        detector_origin, detector_fast, detector_slow, (nx, ny), (px, py),
        rotation_axis, sample_to_source, wavelength,
        real_space_a, real_space_b, real_space_c, space_group_number,
        None, None, starting_angle, oscillation_range, starting_frame,
        original_rotation=R,
        panel_offset=panel_offset,
        panel_size=panel_size,
        panel_origin=panel_origins,
        panel_fast=panel_fast_axes,
        panel_slow=panel_slow_axes)

def test_align_reference_frame():

    _i = (1, 0, 0)
    _j = (0, 1, 0)
    _k = (0, 0, 1)

    primary_axis = _i
    primary_target = _i
    secondary_axis = _k
    secondary_target = _k

    m = align_reference_frame(primary_axis, primary_target,
                              secondary_axis, secondary_target)

    i = matrix.identity(3)

    for j in range(9):
        assert(math.fabs(m.elems[j] - i.elems[j]) < 0.001)

    primary_axis = _j
    primary_target = _i
    secondary_axis = _k
    secondary_target = _k

    m = align_reference_frame(primary_axis, primary_target,
                              secondary_axis, secondary_target)

    for j in range(3):
        assert(math.fabs((m * primary_axis).elems[j] -
                         matrix.col(primary_target).elems[j]) < 0.001)

def test_align_reference_frame_dw():

    s = math.sqrt(0.5)

    pa = [s, s, 0]
    pt = [1, 0, 0]
    sa = [- s, s, 0]
    st = [0, 1, 0]

    R = align_reference_frame(pa,pt,sa,st)

    print(R * pa)
    print(pt)
    print(R * sa)
    print(st)

def random_orthogonal_vectors():
    v1 = matrix.col((random.random(), random.random(),
                     random.random())).normalize()
    v2 = v1.ortho().normalize()

    return v1, v2

def test_align_reference_frame_brute():

    for j in range(10000):
        m = random_orthogonal_vectors()
        t = random_orthogonal_vectors()

        assert(math.fabs(m[0].dot(m[1])) < 0.001)
        assert(math.fabs(t[0].dot(t[1])) < 0.001)

        R = align_reference_frame(m[0], t[0],
                                  m[1], t[1])

        r = (R * m[0], R * m[1])

        assert(math.fabs(r[0].dot(t[0])) > 0.999)
        assert(math.fabs(r[1].dot(t[1])) > 0.999)

    return

def find_closest_matrix(moving, target):
    '''Work through lattice permutations to try to align moving with target,
    with the metric of trace(inverse(moving) * target).'''

    trace = 0.0
    reindex = matrix.identity(3)

    for op in sgtbx.space_group_info('P422').type().group().all_ops():
        moved = matrix.sqr(op.r().as_double()) * moving
        if (moved.inverse() * target).trace() > trace:
            trace = (moved.inverse() * target).trace()
            reindex = matrix.sqr(op.r().as_double())

    return reindex

def work():
    import sys
    import_xds_integrate_hkl(sys.argv[1])
    print('OK')

if __name__ == '__main__':
    work()
