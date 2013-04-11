from __future__ import division
def print_sweep(list_of_images):

    from dxtbx.sweep import SweepFactory

    s = SweepFactory.sweep(list_of_images)

    print s.get_detector()
    print s.get_beam()
    print s.get_goniometer()
    print s.get_scan()

if __name__ == '__main__':
    import sys

    if len(sys.argv) == 2:
        print_sweep(sys.argv[1])
    else:
        print_sweep(sys.argv[1:])
