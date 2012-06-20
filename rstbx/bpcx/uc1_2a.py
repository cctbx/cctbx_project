#!/usr/bin/env cctbx.python
#
# Use case 1.2A - going back to the use case for predictions and making sure that
# the predicted indices are actually identical to those from XDS, from both a
# brute force calculation and the Reeke method. These will start from an
# integrated data set from XDS.

# want methods to: read indices from XDS INTEGRATE.HKL, transform to our "standard
# frame" which is -
#
# h k l xcen ycen phicen
#
# where xcen, ycen are on the camera in mm in direction fast, slow respectively
# with origin outer corner, so range goes from [0, width], [0, height]. N.B.
# that this is identical to XDS frame modulo pixel size. Also phicen in degrees.

# 1 read INTEGRATE.HKL - this frame, create map (x, y, phi) -> HKL. Can then use
# ANN to 'find' a prediction and verify HKL later on.

def import_xds_integrate_hkl(xds_integrate_hkl_file):
    from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
    cfc = coordinate_frame_converter(xds_integrate_hkl_file)

    px, py = cfc.get('detector_pixel_size_fast_slow')

    # read header, get out phi0, frame0, dphi, so can transform frame# to
    # phi value

    phi0 = None
    dphi = None
    frame0 = None

    for record in open(xds_integrate_hkl_file):
        if not record.startswith('!'):
            break
        if record.startswith('!STARTING_FRAME'):
            frame0 = int(record.split()[-1])
            continue
        if record.startswith('!STARTING_ANGLE'):
            phi0 = float(record.split()[-1])
            continue
        if record.startswith('!OSCILLATION_RANGE'):
            dphi = float(record.split()[-1])
            continue

    assert(not dphi is None)
    assert(not phi0 is None)
    assert(not frame0 is None)

    xyz_to_hkl = { }

    for record in open(xds_integrate_hkl_file):
        if record.startswith('!'):
            continue

        values = record.split()

        hkl = map(int, values[:3])
        xyz = map(float, values[5:8])

        xyz_mod = (xyz[0] * px, xyz[1] * py,
                   (xyz[2] - frame0) * dphi + phi0)

        xyz_to_hkl[xyz_mod] = hkl

    return xyz_to_hkl

# 2 read coordinate frame from same INTEGRATE.HKL, provide number of frames as
# input parameter, predict all reflections given camera position in same frame as
# above, store in same structure. Verify that all reflections in INTEGRATE.HKL
# appear in prediction list. N.B. lists not identical, as some reflections on dead
# regions of camera. Compare HKL, (x, y, phi)

def verify_predictions_against_integrate_hkl(xds_integrate_hkl_file,
                                             phi_range):
    from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
    from rstbx.diffraction import rotation_angles, reflection_prediction
    from rstbx.diffraction import full_sphere_indices
    from cctbx.sgtbx import space_group, space_group_symbols
    from cctbx.uctbx import unit_cell
    import math

    cfc = coordinate_frame_converter(xds_integrate_hkl_file)

    d2r = math.pi / 180.0

    dmin = cfc.derive_detector_highest_resolution()

    A = cfc.get_c('real_space_a')
    B = cfc.get_c('real_space_b')
    C = cfc.get_c('real_space_c')

    cell = (A.length(), B.length(), C.length(), B.angle(C, deg = True),
            C.angle(A, deg = True), A.angle(B, deg = True))

    uc = unit_cell(cell)
    sg = cfc.get('space_group_number')

    indices = full_sphere_indices(
        unit_cell = uc,
        resolution_limit = dmin,
        space_group = space_group(space_group_symbols(sg).hall()))

    u, b = cfc.get_u_b(convention = cfc.ROSSMANN)
    axis = cfc.get('rotation_axis', convention = cfc.ROSSMANN)
    ub = u * b

    wavelength = cfc.get('wavelength')

    ra = rotation_angles(dmin, ub, wavelength, axis)

    obs_indices, obs_angles = ra.observed_indices_and_angles_from_angle_range(
        phi_start_rad = phi_range[0] * d2r,
        phi_end_rad = phi_range[1] * d2r,
        indices = indices)

    # in here work in internal (i.e. not Rossmann) coordinate frame

    u, b = cfc.get_u_b()
    axis = cfc.get_c('rotation_axis')
    sample_to_source_vec = cfc.get_c('sample_to_source').normalize()
    s0 = (- 1.0 / wavelength) * sample_to_source_vec
    ub = u * b

    detector_origin = cfc.get_c('detector_origin')
    detector_fast = cfc.get_c('detector_fast')
    detector_slow = cfc.get_c('detector_slow')
    detector_normal = detector_fast.cross(detector_slow)
    distance = detector_origin.dot(detector_normal.normalize())
    nx, ny = cfc.get('detector_size_fast_slow')
    px, py = cfc.get('detector_pixel_size_fast_slow')

    limits = [0, nx * px, 0, ny * py]

    xyz_to_hkl = { }

    for hkl, angle in zip(obs_indices, obs_angles):
        s = (ub * hkl).rotate(axis, angle)
        q = (s + s0).normalize()

        # check if diffracted ray parallel to detector face

        q_dot_n = q.dot(detector_normal)

        if q_dot_n == 0:
            continue

        r = (q * distance / q_dot_n) - detector_origin

        x = r.dot(detector_fast)
        y = r.dot(detector_slow)

        if x < limits[0] or y < limits[2]:
            continue
        if x > limits[1] or y > limits[3]:
            continue

        xyz = (x, y, angle / d2r)
        xyz_to_hkl[xyz] = map(int, hkl)

    xyz_to_hkl_xds = import_xds_integrate_hkl(xds_integrate_hkl_file)

    # construct ann to perform search...

    from cctbx.array_family import flex
    from annlib_ext import AnnAdaptor as ann_adaptor

    reference = flex.double()

    xyzs = [xyz for xyz in xyz_to_hkl]

    for xyz in xyzs:
        reference.append(xyz[0])
        reference.append(xyz[1])
        reference.append(xyz[2])

    ann = ann_adaptor(data = reference, dim = 3, k = 1)

    n_correct = 0
    n_wrong = 0

    for xyz in xyz_to_hkl_xds:
        query = flex.double(xyz)
        ann.query(query)
        nnxyz = xyzs[ann.nn[0]]
        if xyz_to_hkl_xds[xyz] == xyz_to_hkl[nnxyz]:
            n_correct += 1
        else:
            n_wrong += 1

    return n_correct, n_wrong


# 3 same as #2 verify that the predictions are all correct.

# 4 define 'matrix' and 'hkl' file formats - or follow d*TREK model:
#
# {
# HEADER_SIZE=NNNN;
# PARAMETER=VALUE;
# ...
# }
# H K L x y phi
# ...
#
# + read / write class.

def work():
    import sys
    n_correct, n_wrong = verify_predictions_against_integrate_hkl(
        sys.argv[1], (0, 90))
    assert(float(n_correct) / float(n_correct + n_wrong) > 0.999)

    print 'OK'

if __name__ == '__main__':
    work()
