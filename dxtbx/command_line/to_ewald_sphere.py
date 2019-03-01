from __future__ import absolute_import, division, print_function


def to_ewald_sphere(list_of_images):

    from dxtbx.sweep import SweepFactory
    from dxtbx import ImageToEwaldSphere

    sweep = SweepFactory.sweep(list_of_images)
    beam = sweep.get_beam()
    detector = sweep.get_detector()
    gonio = sweep.get_goniometer()
    scan = sweep.get_scan()
    start, end = sweep.get_array_range()
    image_to_ewald_sphere = ImageToEwaldSphere(beam, detector, gonio, scan)

    for frame in xrange(start, end):
        intensity = sweep[frame]
        x_list = image_to_ewald_sphere(frame)
        for k, x in enumerate(x_list):
            i = intensity[k]
            print("%f %f %f %d" % (x[0], x[1], x[2], i))


if __name__ == "__main__":
    import sys

    to_ewald_sphere(sys.argv[1:])
