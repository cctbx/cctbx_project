from __future__ import division

def tst_sweep():
    import os
    from dxtbx.sweep import SweepFactory
    import libtbx.load_env

    # Try to find the dials regression directory
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    # Get the filenames
    filenames = [os.path.join(dials_regression, 'centroid_test_data',
                           'centroid_%04d.cbf' % j) for j in range(1, 10)]

    # Create the sweep
    sweep = SweepFactory.sweep(filenames)
    assert(len(sweep) == 9)
    print "OK"

    # Get the models from the sweep
    sweep.get_beam()
    sweep.get_detector()
    sweep.get_goniometer()
    sweep.get_scan()

    # Get the image size
    image_size = sweep.get_detector().get_image_size()

    # Get a couple of sub-sweeps
    sub_sweep_1 = sweep[0:7]
    sub_sweep_2 = sweep[3:9]

    # Check sweep length
    assert(len(sub_sweep_1) == 7)
    assert(len(sub_sweep_2) == 6)
    print "OK"

    # Get arrays from sub sweeps
    array_0 = sweep.to_array()
    array_1 = sub_sweep_1.to_array()
    array_2 = sub_sweep_2.to_array()

    # Check 3d array sizes
    assert(array_0.all() == (9, image_size[1], image_size[0]))
    assert(array_1.all() == (7, image_size[1], image_size[0]))
    assert(array_2.all() == (6, image_size[1], image_size[0]))
    print "OK"

    # Check sweep is valid
    assert(sweep.is_valid() == True)
    print "OK"

    # Get sub-sections of sweep
    array_1 = sweep.to_array((2, 5))
    array_2 = sweep.to_array((2, 5, 100, 105, 200, 240))

    assert(array_1.all() == (3, image_size[1], image_size[0]))
    assert(array_2.all() == (3, 5, 40))
    print "OK"

    # Loop through all the images in a sweep
    count = 0
    for image in sweep:
        count += 1
    assert(count == 9)
    print "OK"

    # Loop through all the images in a sub sweep
    count = 0
    for image in sweep[3:7]:
        count += 1
    assert(count == 4)
    print "OK"

if __name__ == '__main__':
    tst_sweep()
