from __future__ import division
from dxtbx.model import Scan

def tst_is_angle_valid(scan):
    """Check that the is_angle_valid function behaves properly."""
    oscillation_range = scan.get_oscillation_range()
    os1 = int(oscillation_range[0])
    os2 = int(oscillation_range[1])
    for i in range(0, os1):
        assert(scan.is_angle_valid(i) == False)
    for i in range(os1, os2):
        assert(scan.is_angle_valid(i) == True)
    for i in range(os2+1, 360):
        assert(scan.is_angle_valid(i) == False)
    print "OK"

def tst_is_frame_valid(scan):
    """Check that the is_frame_valid function behaves properly."""
    image_range = scan.get_image_range()
    for i  in range(image_range[0] - 100, image_range[0]):
        assert(scan.is_image_index_valid(i) == False)
    for i in range(image_range[0], image_range[1]+1):
        assert(scan.is_image_index_valid(i) == True)
    for i in range(image_range[1]+1, image_range[1] + 100):
        assert(scan.is_image_index_valid(i) == False)
    print "OK"

def tst_get_angle_from_frame(scan):
    pass

def tst_get_frame_from_angle(scan):
    pass

def tst_get_frames_with_angle(scan):
    pass

def run():
    image_range = (0, 1000)
    oscillation = (0, 0.1)
    scan = Scan(image_range, oscillation, 0.0)
    tst_is_angle_valid(scan)
    tst_is_frame_valid(scan)
    tst_get_angle_from_frame(scan)
    tst_get_frame_from_angle(scan)
    tst_get_frames_with_angle(scan)

if __name__ == '__main__':
    run()
