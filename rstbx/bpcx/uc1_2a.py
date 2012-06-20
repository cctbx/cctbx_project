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
    xyz_to_hkl = import_xds_integrate_hkl(sys.argv[1])
    print 'OK'

if __name__ == '__main__':
    work()
