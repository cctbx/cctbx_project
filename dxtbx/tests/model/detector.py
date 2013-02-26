from __future__ import division
from dxtbx.model import FlatPanelDetector, MultiFlatPanelDetector

def tst_get_pixel_lab_coord(detector):
    from scitbx import matrix
    eps = 1e-7

    # Check lab coordinates at the origin
    orig = detector.get_pixel_lab_coord((0, 0))
    dorig = abs(matrix.col(orig) - matrix.col(detector.origin))
    assert(dorig < eps)

    # Check lab coordinate at the opposite corner
    corner = detector.get_pixel_lab_coord((512, 512))
    corner2 = (512 * 0.172, 512 * 0.172, 200)
    dcorner = abs(matrix.col(corner) - matrix.col(corner))
    assert(dcorner < eps)
    print "OK"

def tst_get_image_rectangle(detector):
    from scitbx import matrix
    eps = 1e-7
    orig = detector.origin
    corner = (512 * 0.172, 512 * 0.172, 200)
    rect1 = (orig[0], orig[1], orig[2], corner[0], corner[1], corner[2])
    rect2 = detector.get_image_rectangle()
    drect = abs(matrix.col(rect1) - matrix.col(rect2))
    assert(drect < eps)
    print "OK"

def tst_get_image_size_mm(detector):
    from scitbx import matrix
    eps = 1e-7
    size = detector.get_image_size_mm()
    size2 = (512 * 0.172, 512 * 0.172)
    dsize = abs(matrix.col(size) - matrix.col(size2))
    assert(dsize < eps)
    print "OK"

def tst_is_value_in_trusted_range(detector):
    """Check values are either inside or outside trusted range."""
    assert(detector.is_value_in_trusted_range(-1) == False)
    assert(detector.is_value_in_trusted_range(0) == True)
    assert(detector.is_value_in_trusted_range(999) == True)
    assert(detector.is_value_in_trusted_range(1000) == False)
    print "OK"

def tst_is_coord_valid(detector):
    """Check points are either inside or outside detector range."""
    assert(detector.is_coord_valid((-1, 256)) == False)
    assert(detector.is_coord_valid((256, 256)) == True)
    assert(detector.is_coord_valid((512, 256)) == False)
    assert(detector.is_coord_valid((256, -1)) == False)
    assert(detector.is_coord_valid((256, 256)) == True)
    assert(detector.is_coord_valid((256, 513)) == False)
    print "OK"

def tst_pixel_to_millimeter_to_pixel(detector):

    from scitbx import matrix
    from random import random
    eps = 1e-7

    # Pick some random pixels and check that px -> mm -> px give px == px
    w, h = detector.image_size
    random_pixel = lambda: (random() * w, random() * h)
    for i in range(100):
        xy = random_pixel()
        xy_mm = detector.pixel_to_millimeter(xy)
        xy_px = detector.millimeter_to_pixel(xy_mm)
        assert(abs(matrix.col(xy_px) - matrix.col(xy)) < eps)

    # Test Passed
    print "OK"

def tst_flat_panel_detector():

    # Create the detector
    detector = FlatPanelDetector(
        "",                 # Type
        (10, 0, 0),         # Fast axis
        (0, 10, 0),         # Slow axis
        (0, 0, 200),        # Origin
        (0.172, 0.172),     # Pixel size
        (512, 512),         # Image size
        (0, 1000))          # Trusted range

    # Perform some tests
    tst_get_pixel_lab_coord(detector)
    tst_get_image_rectangle(detector)
    tst_get_image_size_mm(detector)
    tst_is_value_in_trusted_range(detector)
    tst_is_coord_valid(detector)
    tst_pixel_to_millimeter_to_pixel(detector)

def tst_multi_flat_panel_detector():
    print "OK"

def run():
    tst_flat_panel_detector()
    tst_multi_flat_panel_detector()

if __name__ == '__main__':
    run()
