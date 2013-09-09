from __future__ import division
from dxtbx.model import Panel, Detector

def tst_get_pixel_lab_coord(detector):
    from scitbx import matrix
    eps = 1e-7

    # Check lab coordinates at the origin
    orig = detector[0].get_pixel_lab_coord((0, 0))
    dorig = abs(matrix.col(orig) - matrix.col(detector[0].get_origin()))
    assert(dorig < eps)

    # Check lab coordinate at the opposite corner
    corner = detector[0].get_pixel_lab_coord((512, 512))
    corner2 = (512 * 0.172, 512 * 0.172, 200)
    dcorner = abs(matrix.col(corner) - matrix.col(corner))
    assert(dcorner < eps)
    print "OK"

def tst_get_image_size_mm(detector):
    from scitbx import matrix
    eps = 1e-7
    size = detector[0].get_image_size_mm()
    size2 = (512 * 0.172, 512 * 0.172)
    dsize = abs(matrix.col(size) - matrix.col(size2))
    assert(dsize < eps)
    print "OK"

def tst_is_value_in_trusted_range(detector):
    """Check values are either inside or outside trusted range."""
    assert(detector[0].is_value_in_trusted_range(-1) == False)
    assert(detector[0].is_value_in_trusted_range(0) == True)
    assert(detector[0].is_value_in_trusted_range(999) == True)
    assert(detector[0].is_value_in_trusted_range(1000) == False)
    print "OK"

def tst_is_coord_valid(detector):
    """Check points are either inside or outside detector range."""
    assert(detector[0].is_coord_valid((-1, 256)) == False)
    assert(detector[0].is_coord_valid((256, 256)) == True)
    assert(detector[0].is_coord_valid((512, 256)) == False)
    assert(detector[0].is_coord_valid((256, -1)) == False)
    assert(detector[0].is_coord_valid((256, 256)) == True)
    assert(detector[0].is_coord_valid((256, 513)) == False)
    print "OK"

def tst_pixel_to_millimeter_to_pixel(detector):

    from scitbx import matrix
    from random import random
    eps = 1e-7

    # Pick some random pixels and check that px -> mm -> px give px == px
    w, h = detector[0].get_image_size()
    random_pixel = lambda: (random() * w, random() * h)
    for i in range(100):
        xy = random_pixel()
        xy_mm = detector[0].pixel_to_millimeter(xy)
        xy_px = detector[0].millimeter_to_pixel(xy_mm)
        assert(abs(matrix.col(xy_px) - matrix.col(xy)) < eps)

    # Test Passed
    print "OK"

def tst_parallax_correction(detector):
    from random import uniform
    from scitbx import matrix
    random_coord = lambda: (
        uniform(-1000, 1000),
        uniform(-1000, 1000))
    for i in range(10000):
        mm = random_coord()
        px = detector.millimeter_to_pixel(mm)
        mm2 = detector.pixel_to_millimeter(px)
        assert(abs(matrix.col(mm) - matrix.col(mm2)) < 1e-3)

    print 'OK'

def tst_get_names(detector):
    names = detector.get_names()
    assert(len(names) == 1)
    assert(names[0] == 'Panel')
    print 'OK'

def tst_detector():
    from dxtbx.model import ParallaxCorrectedPxMmStrategy

    # Create the detector
    detector = Detector(Panel(
        "",                 # Type
        "Panel",            # Name
        (10, 0, 0),         # Fast axis
        (0, 10, 0),         # Slow axis
        (0, 0, 200),        # Origin
        (0.172, 0.172),     # Pixel size
        (512, 512),         # Image size
        (0, 1000)))         # Trusted range

    # Perform some tests
    tst_get_pixel_lab_coord(detector)
    tst_get_image_size_mm(detector)
    tst_is_value_in_trusted_range(detector)
    tst_is_coord_valid(detector)
    tst_pixel_to_millimeter_to_pixel(detector)
    tst_get_names(detector)

    # Attenuation length
    la = 0.252500934883

    # Create the detector
    detector = Detector(Panel(
        "",                 # Type
        "",                 # Name
        (10, 0, 0),         # Fast axis
        (0, 10, 0),         # Slow axis
        (0, 0, 200),        # Origin
        (0.172, 0.172),     # Pixel size
        (512, 512),         # Image size
        (0, 1000),          # Trusted range
        ParallaxCorrectedPxMmStrategy(la)))

    tst_parallax_correction(detector)

def run():
    tst_detector()

if __name__ == '__main__':
    run()
