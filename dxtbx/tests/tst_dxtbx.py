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
        det = i.get_detector()
        if det is not None:
            size = det[0].get_image_size()
        b = i.get_beam()
        g = i.get_goniometer()
        s = i.get_scan()
        try:
            d = i.get_raw_data()
        except IOError:
            pass

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

def tst_sweep():

    from dxtbx.sweep_filenames import template_regex

    questions_answers = {
        'foo_bar_001.img':'foo_bar_###.img',
        'foo_bar001.img':'foo_bar###.img',
        'foo_bar_1.8A_001.img':'foo_bar_1.8A_###.img',
        'foo_bar.001':'foo_bar.###',
        'foo_bar_001.img1000':'foo_bar_###.img1000',
        'foo_bar_00001.img':'foo_bar_#####.img'
        }

    for filename in questions_answers:
        answer = template_regex(filename)
        assert answer[0] == questions_answers[filename]

    print 'OK'

if __name__ == '__main__':
    tst_dxtbx()
    tst_dxtbx_models()  # these tests are failing in the misc_build.  Disable them until they can be debugged.
