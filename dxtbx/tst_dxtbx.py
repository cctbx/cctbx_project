def tst_dxtbx():
    import libtbx.load_env
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    import os
    from boost.python import streambuf
    from dxtbx import read_uint16
    from dxtbx.format.Registry import Registry

    from dials_regression.image_examples.get_all_working_images import \
        get_all_working_images

    for directory, image in get_all_working_images():
        file_path = os.path.join(dials_regression, 'image_examples',
                                 directory, image)
        format = Registry.find(file_path)
        i = format(file_path)
        size = i.get_detector().get_image_size()

    print 'OK'

tst_dxtbx()
