from __future__ import division
def to_ewald_sphere(list_of_images):

    from dxtbx.sweep import SweepFactory
    from scitbx import matrix

    s = SweepFactory.sweep(list_of_images)

    start, end = s.get_array_range()
    size = s.get_detector().get_image_size()

    F = matrix.col(s.get_detector().get_fast_axis())
    S = matrix.col(s.get_detector().get_slow_axis())
    O = matrix.col(s.get_detector().get_origin())

    pixel_size = s.get_detector().get_pixel_size()

    F *= pixel_size[0]
    S *= pixel_size[1]

    wavelength = s.get_beam().get_wavelength()

    phi_start = s.get_scan().get_oscillation()[0]
    phi_width = s.get_scan().get_oscillation()[1]
    axis = matrix.col(s.get_goniometer().get_rotation_axis())

    for frame in range(start, end):
        phi = (frame - 0.5) * phi_width + phi_start
        R = axis.axis_and_angle_as_r3_rotation_matrix(phi, deg = True)
        for j in range(size[1]):
            for i in range(size[0]):
                intensity = s[frame][(j, i)]
                x = R * (O + F * i + S * j).normalize() / wavelength
                print '%f %f %f %d' % (
                    x.elems[0], x.elems[1], x.elems[2], intensity)


if __name__ == '__main__':
    import sys

    to_ewald_sphere(sys.argv[1:])
