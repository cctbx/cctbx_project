from __future__ import division
def tst_dxtbx():
    import libtbx.load_env
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    import os
    #from boost.python import streambuf
    #from dxtbx import read_uint16
    from dxtbx.format.Registry import Registry

    from dials_regression.image_examples.get_all_working_images import \
        get_all_working_images

    for directory, image in get_all_working_images():
        file_path = os.path.join(dials_regression, 'image_examples',
                                 directory, image)
        format = Registry.find(file_path)
        i = format(file_path)
        size = i.get_detector().get_image_size()
        b = i.get_beam()
        g = i.get_goniometer()
        s = i.get_scan()

    print 'OK'

def tst_dxtbx_models():

    from dxtbx.tests.test_beam import test_beam
    test_beam()

    from dxtbx.tests.test_detector import test_detector
    test_detector()

    from dxtbx.tests.test_goniometer import test_goniometer
    test_goniometer()

    from dxtbx.tests.test_scan import test_scan
    test_scan()

    return

if __name__ == '__main__':
    tst_dxtbx()
    tst_dxtbx_models()
