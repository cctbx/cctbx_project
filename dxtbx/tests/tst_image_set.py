from __future__ import division

def run():
    import libtbx.load_env
    import os

    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    path = os.path.join(dials_regression, 'centroid_test_data')

    # Non-sequential Filenames and image indices
    filenames = []
    image_indices = [7, 4, 6, 2, 1]
    for i in image_indices:
        filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

    # Read the image set
    from dxtbx.image_set import ImageSetFactory
    imageset = ImageSetFactory.image_set(filenames)

    # Check imageset is the right length
    assert(len(imageset) == len(filenames))

    # Check images are in given order
    assert(imageset.indices() == image_indices)

    # Read through all the images
    for image in imageset:
        pass

    # Try by index
    for i in range(len(imageset)):
        image = imageset[i]

    # Try to get a subset
    sub_imageset = imageset[1:4]
    assert(len(sub_imageset) == 3)
    assert(sub_imageset.indices() == [4, 6, 2])
    for image in sub_imageset:
        pass
    for i in range(len(sub_imageset)):
        image = sub_imageset[i]

    # Check image set is valid
    assert(imageset.is_valid() == True)

if __name__ == '__main__':
    run()
