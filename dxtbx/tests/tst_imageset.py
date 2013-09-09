from __future__ import division

class TestMultiFileState(object):

    def __init__(self):
        pass

    def get_filenames(self):
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
        image_indices = range(1, 10)
        for i in image_indices:
            filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

        return filenames

    def run(self):
        from dxtbx.imageset import MultiFileState
        from dxtbx.format.Registry import Registry

        # Get the filenames
        filenames = self.get_filenames()

        # Get the parameters we need
        format_class = Registry.find(filenames[0])

        # Create the state object
        state = MultiFileState(format_class)

        # Run a load of tests
        self.tst_format_class(state, format_class)
        self.tst_load_file(state, filenames)

    def tst_format_class(self, state, format_class):
        assert(state.format_class() == format_class)
        print 'OK'

    def tst_load_file(self, state, filenames):
        for f in filenames:
            state.load_file(f)
        print 'OK'


class TestSweepFileList(object):
    def __init__(self):
        pass

    def run(self):
        from dxtbx.imageset import SweepFileList

        # Create the template and array range
        template = template = '%s%%0%dd%s' % ('filename', 5, '.cbf')
        array_range = (20, 50)

        # Create the sweep file list class
        filelist = SweepFileList(template, array_range)

        # Call a load of tests
        self.tst_len(filelist, 30)
        self.tst_template(filelist, template)
        self.tst_array_range(filelist, array_range)
        self.tst_indices(filelist, range(*array_range))
        self.tst_is_index_in_range(filelist)
        self.tst_get_filename(filelist)
        self.tst_get_item(filelist)
        self.tst_iter(filelist)
        self.tst_cmp(filelist)

    def tst_len(self, filelist, length):
        assert(len(filelist) == length)
        print 'OK'

    def tst_template(self, filelist, template):
        assert(filelist.template() == template)
        print 'OK'

    def tst_array_range(self, filelist, array_range):
        assert(filelist.array_range() == array_range)
        print 'OK'

    def tst_indices(self, filelist, indices):
        assert(filelist.indices() == indices)
        print 'OK'

    def tst_is_index_in_range(self, filelist):
        array_range = filelist.array_range()
        for i in range(array_range[0] - 10, array_range[0]):
            assert(filelist.is_index_in_range(i) == False)
        for i in range(array_range[0], array_range[1]):
            assert(filelist.is_index_in_range(i) == True)
        for i in range(array_range[1], array_range[1] + 10):
            assert(filelist.is_index_in_range(i) == False)
        print 'OK'

    def tst_get_filename(self, filelist):
        template = filelist.template()
        array_range = filelist.array_range()
        too_low = array_range[0] - 1
        too_high = array_range[1]
        just_right = (array_range[0] + array_range[1]) // 2
        try:
            filename = filelist.get_filename(too_low)
            assert(False)
        except AssertionError:
            raise
        except Exception:
            pass
        try:
            filename = filelist.get_filename(too_high)
            assert(False)
        except AssertionError:
            raise
        except Exception:
            pass
        filename = filelist.get_filename(just_right)
        assert(filename == template % (just_right + 1))
        print 'OK'

    def tst_get_item(self, filelist):
        template = filelist.template()
        array_range = filelist.array_range()
        too_low = -1
        too_high = len(filelist)
        just_right = len(filelist) / 2
        try:
            filename = filelist[too_low]
            assert(False)
        except AssertionError:
            raise
        except Exception:
            pass
        try:
            filename = filelist[too_high]
            assert(False)
        except AssertionError:
            raise
        except Exception:
            pass
        filename = filelist[just_right]
        assert(filename == template % (array_range[0] + just_right + 1))
        print 'OK'

    def tst_iter(self, filelist):
        template = filelist.template()
        array_range = filelist.array_range()
        for index, filename in enumerate(filelist):
            assert(filename == template % (array_range[0] + index + 1))
        print 'OK'

    def tst_cmp(self, filelist):
        from dxtbx.imageset import SweepFileList

        assert(filelist == filelist)

        template = template = '%s%%0%dd%s' % ('filename', 5, '.cbf')
        array_range = (0, 10)
        filelist2 = SweepFileList(template, array_range)
        assert(filelist != filelist2)

        template = template = '%s%%0%dd%s' % ('filename2', 5, '.cbf')
        array_range = (20, 50)
        filelist3 = SweepFileList(template, array_range)
        assert(filelist != filelist3)

        template = template = '%s%%0%dd%s' % ('filename', 5, '.cbf')
        array_range = (20, 50)
        filelist4 = SweepFileList(template, array_range)
        assert(filelist == filelist4)

        print 'OK'


class TestMultiFileReader(object):

    def __init__(self):
        pass

    def normal_file_list(self):
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
        image_indices = range(1, 10)
        for i in image_indices:
            filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

        return filenames

    def sweep_file_list(self):
        import libtbx.load_env
        import os
        from dxtbx.imageset import SweepFileList

        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        path = os.path.join(dials_regression, 'centroid_test_data')

        # Non-sequential Filenames and image indices
        template = os.path.join(path, 'centroid_%04d.cbf')
        array_range = (0, 9)

        filenames = SweepFileList(template, array_range)

        return filenames

    def run(self):
        self.run_tests(self.normal_file_list())
        self.run_tests(self.sweep_file_list())

    def run_tests(self, filenames):

        from dxtbx.imageset import MultiFileReader
        from dxtbx.format.Registry import Registry

        # Get the parameters we need
        format_class = Registry.find(filenames[0])

        # Create the reader
        reader = MultiFileReader(format_class, filenames)

        # Run a load of tests
        self.tst_get_image_paths(reader, filenames)
        self.tst_get_format_class(reader, format_class)
        self.tst_get_image_size(reader)
        self.tst_get_path(reader, filenames)
        self.tst_get_detectorbase(reader)
        self.tst_get_models(reader)
        self.tst_read(reader)
        self.tst_get_format(reader)
        self.tst_is_valid(reader)

    def tst_get_image_paths(self, reader, filenames):
        filenames2 = reader.get_image_paths()
        assert(len(filenames2) == len(filenames))
        for f1, f2 in zip(filenames, filenames2):
            assert(f1 == f2)
        print 'OK'

    def tst_get_format_class(self, reader, format_class):
        assert(reader.get_format_class() == format_class)
        print 'OK'

    def tst_get_image_size(self, reader):
        image_size = reader.get_image_size()
        for index in range(0, len(reader.get_image_paths())):
            reader.read(index)
            assert(image_size == reader.get_image_size())

    def tst_get_path(self, reader, filenames):
        for index in range(len(filenames)):
            assert(reader.get_path(index) == filenames[index])
        print 'OK'

    def tst_get_detectorbase(self, reader):
        reader.get_detectorbase()
        reader.get_detectorbase(0)
        reader.get_detectorbase(4)
        try:
            reader.get_detectorbase(9)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models(self, reader):
        self.tst_get_models_index(reader, None)
        self.tst_get_models_index(reader, 0)
        self.tst_get_models_index(reader, 4)
        try:
            self.tst_get_models_index(reader, 9)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models_index(self, reader, index):
        reader.get_detector(index)
        reader.get_goniometer(index)
        reader.get_beam(index)
        reader.get_scan(index)

    def tst_read(self, reader):
        for index in range(0, len(reader.get_image_paths())):
            reader.read(index)
        try:
            reader.read(9)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_format(self, reader):
        reader.get_format()
        reader.get_format(0)
        reader.get_format(4)
        try:
            reader.get_format(9)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_is_valid(self, reader):
        assert(reader.is_valid() == True)
        assert(reader.is_valid([0, 1, 2, 3]) == True)
        assert(reader.is_valid([5, 6, 7, 8]) == True)
        assert(reader.is_valid([9, 10]) == False)


class TestImageSet(object):

    def __init__(self):
        pass

    def get_file_list(self):
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
        image_indices = range(1, 10)
        for i in image_indices:
            filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

        return filenames

    def run(self):
        from dxtbx.imageset import MultiFileReader, ImageSet
        from dxtbx.format.Registry import Registry

        # Get the filenames
        filenames = self.get_file_list()

        # Create the format class
        format_class = Registry.find(filenames[0])

        # Create the reader
        reader = MultiFileReader(format_class, filenames)

        # Create the imageset
        imageset = ImageSet(reader)

        # Run a load of tests
        self.tst_get_item(imageset)
        self.tst_len(imageset, len(filenames))
        self.tst_iter(imageset)
        self.tst_indices(imageset, range(0, 9))
        self.tst_paths(imageset, filenames)
        self.tst_is_valid(imageset)
        self.tst_get_detectorbase(imageset, range(len(filenames)), 9)
        self.tst_get_models(imageset, range(len(filenames)), 9)

    def tst_get_item(self, imageset):
        image = imageset[0]
        try:
            image = imageset[9]
            assert(False)
        except IndexError:
            pass

        imageset2 = imageset[3:7]
        image = imageset2[0]
        try:
            image = imageset2[5]
            assert(False)
        except IndexError:
            pass

        self.tst_len(imageset2, 4)
        self.tst_is_valid(imageset2)
        self.tst_get_detectorbase(imageset2, range(0, 4), 5)
        self.tst_get_models(imageset2, range(0, 4), 5)
        self.tst_paths(imageset2, imageset.paths()[3:7])
        self.tst_indices(imageset2, range(3, 7))
        self.tst_iter(imageset2)

        imageset2 = imageset[3:7:2]
        image = imageset2[0]
        try:
            image = imageset2[2]
            assert(False)
        except IndexError:
            pass

        self.tst_len(imageset2, 2)
        self.tst_is_valid(imageset2)
        self.tst_get_detectorbase(imageset2, range(0, 2), 2)
        self.tst_get_models(imageset2, range(0, 2), 2)
        self.tst_paths(imageset2, imageset.paths()[3:7:2])
        self.tst_indices(imageset2, range(3, 7, 2))
        self.tst_iter(imageset2)

        print 'OK'

    def tst_len(self, imageset, length):
        assert(len(imageset) == length)
        print 'OK'

    def tst_iter(self, imageset):
        for image in imageset:
            pass
        print 'OK'

    def tst_indices(self, imageset, indices2):
        indices1 = imageset.indices()
        for i1, i2 in zip(indices1, indices2):
            assert(i1 == i2)
        print 'OK'

    def tst_paths(self, imageset, filenames1):
        filenames2 = imageset.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert(f1 == f2)
        print 'OK'

    def tst_is_valid(self, imageset):
        assert(imageset.is_valid() == True)
        print 'OK'

    def tst_get_detectorbase(self, imageset, indices, outside_index):
        imageset.get_detectorbase()
        for i in indices:
            imageset.get_detectorbase(i)

        try:
            imageset.get_detectorbase(outside_index)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models(self, imageset, indices, outside_index):
        self.tst_get_models_index(imageset)
        for i in indices:
            self.tst_get_models_index(imageset, i)

        try:
            self.tst_get_models_index(imageset, outside_index)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models_index(self, imageset, index=None):
        imageset.get_detector(index)
        imageset.get_beam(index)

    def tst_set_models(self, imageset):
        from dxtbx.model import Beam, Detector, Panel

        # Create some other models
        beam = Beam((1, 0, 0), 0.5)
        detector = Detector(Panel("UNKNOWN", "Panel",
                                  (1, 0, 0), (0, 1, 0), (0, 0, 1),
                                  (0.1, 0.1), (1000, 1000), (0, 1)))

        # Override sweep models
        imageset.set_beam(beam)
        imageset.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in imageset:
            pass

        # Get the models back and check they're ok
        beam2 = imageset.get_beam()
        detector2 = imageset.get_detector()
        assert(beam2 == beam)
        assert(detector2 == detector)

        # Get the models from an index back and check they're not the same
        beam2 = imageset.get_beam(0)
        detector2 = imageset.get_detector(0)
        assert(beam2 != beam)
        assert(detector2 != detector)


class TestImageSweep(object):

    def __init__(self):
        pass

    def get_file_list(self):
        import libtbx.load_env
        import os
        from dxtbx.imageset import SweepFileList

        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        path = os.path.join(dials_regression, 'centroid_test_data')

        # Non-sequential Filenames and image indices
        template = os.path.join(path, 'centroid_%04d.cbf')
        array_range = (0, 9)

        filenames = SweepFileList(template, array_range)

        return filenames

    def run(self):
        from dxtbx.imageset import MultiFileReader, ImageSweep
        from dxtbx.format.Registry import Registry

        # Get the filenames
        filenames = self.get_file_list()

        # Create the format class
        format_class = Registry.find(filenames[0])

        # Create the reader
        reader = MultiFileReader(format_class, filenames)

        # Create the sweep
        sweep = ImageSweep(reader)

        # Run a load of tests
        self.tst_get_item(sweep)
        self.tst_len(sweep, len(filenames))
        self.tst_iter(sweep)
        self.tst_indices(sweep, range(0, 9))
        self.tst_paths(sweep, filenames)
        self.tst_is_valid(sweep)
        self.tst_get_detectorbase(sweep, range(len(filenames)), 9)
        self.tst_get_models(sweep, range(len(filenames)), 9)
        self.tst_get_array_range(sweep, (0, 9))
        self.tst_to_array(sweep, (3, 7), (3, 7, 50, 100, 100, 200))
        self.tst_set_models(sweep)

    def tst_get_item(self, sweep):
        image = sweep[0]
        try:
            image = sweep[9]
            assert(False)
        except IndexError:
            pass

        sweep2 = sweep[3:7]
        image = sweep2[0]
        try:
            image = sweep2[5]
            assert(False)
        except IndexError:
            pass

        self.tst_len(sweep2, 4)
        self.tst_is_valid(sweep2)
        self.tst_get_detectorbase(sweep2, range(0, 4), 5)
        self.tst_get_models(sweep2, range(0, 4), 5)
        self.tst_paths(sweep2, sweep.paths()[3:7])
        self.tst_indices(sweep2, range(3, 7))
        self.tst_iter(sweep2)
        self.tst_get_array_range(sweep2, (3, 7))
        self.tst_to_array(sweep2, (1, 3), (1, 3, 50, 100, 100, 200))

        try:
            sweep2 = sweep[3:7:2]
            assert(False)
        except IndexError:
            pass

        print 'OK'

    def tst_len(self, sweep, length):
        assert(len(sweep) == length)
        print 'OK'

    def tst_iter(self, sweep):
        for image in sweep:
            pass
        print 'OK'

    def tst_indices(self, sweep, indices2):
        indices1 = sweep.indices()
        for i1, i2 in zip(indices1, indices2):
            assert(i1 == i2)
        print 'OK'

    def tst_paths(self, sweep, filenames1):
        filenames2 = sweep.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert(f1 == f2)
        print 'OK'

    def tst_is_valid(self, sweep):
        assert(sweep.is_valid() == True)
        print 'OK'

    def tst_get_detectorbase(self, sweep, indices, outside_index):
        sweep.get_detectorbase()
        for i in indices:
            sweep.get_detectorbase(i)

        try:
            sweep.get_detectorbase(outside_index)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models(self, sweep, indices, outside_index):
        self.tst_get_models_index(sweep)
        for i in indices:
            self.tst_get_models_index(sweep, i)

        try:
            self.tst_get_models_index(sweep, outside_index)
            assert(False)
        except IndexError:
            pass
        print 'OK'

    def tst_get_models_index(self, sweep, index=None):
        sweep.get_detector(index)
        sweep.get_beam(index)
        sweep.get_goniometer(index)
        sweep.get_scan(index)

        # Ensure state at zero
        sweep[0]
        scan1 = sweep.get_scan()
        # Put sweep to end
        sweep[len(sweep)-1]
        scan2 = sweep.get_scan()
        assert(scan1 == scan2)
        print 'OK'

    def tst_get_array_range(self, sweep, array_range):
        assert(sweep.get_array_range() == array_range)
        print 'OK'

    def tst_to_array(self, sweep, array_range, sub_section):
        self.tst_to_array_all(sweep)
        self.tst_to_array_num_frames(sweep, array_range)
        self.tst_to_array_sub_section(sweep, sub_section)
        print 'OK'

    def tst_to_array_all(self, sweep):
        volume = sweep.to_array()
        image_size = sweep.get_image_size()
        num_frames = len(sweep)
        assert(volume.all() == (num_frames, image_size[1], image_size[0]))
        print 'OK'

    def tst_to_array_num_frames(self, sweep, array_range):
        volume = sweep.to_array(array_range)
        image_size = sweep.get_image_size()
        assert(volume.all() == (array_range[1] - array_range[0],
            image_size[1], image_size[0]))

    def tst_to_array_sub_section(self, sweep, sub_section):
        volume = sweep.to_array(sub_section)
        size = (sub_section[1]-sub_section[0],
                sub_section[3]-sub_section[2],
                sub_section[5]-sub_section[4])
        assert(volume.all() == size)

    def tst_set_models(self, sweep):
        from dxtbx.model import Beam, Detector, Panel

        # Get some models
        beam = sweep.get_beam()
        gonio = sweep.get_goniometer()
        detector = sweep.get_detector()

        # Modify the geometry
        beam.set_direction((1, 0, 0))
        gonio.set_rotation_axis((0, 1, 0))
        detector.set_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))

        # Override sweep models
        sweep.set_beam(beam)
        sweep.set_goniometer(gonio)
        sweep.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in sweep:
            pass

        # Get the models back and check they're ok
        beam2 = sweep.get_beam()
        gonio2 = sweep.get_goniometer()
        detector2 = sweep.get_detector()
        assert(beam2 == beam)
        assert(gonio2 == gonio)
        assert(detector2 == detector)

        # Get the models from an index back and check they're not the same
        beam2 = sweep.get_beam(0)
        gonio2 = sweep.get_goniometer(0)
        detector2 = sweep.get_detector(0)
        assert(beam2 != beam)
        assert(gonio2 != gonio)
        assert(detector2 != detector)

        # Get a sub sweep
        sub_sweep = sweep[3:7]

        # Get the models back and check they're ok
        beam2 = sub_sweep.get_beam()
        gonio2 = sub_sweep.get_goniometer()
        detector2 = sub_sweep.get_detector()
        assert(beam2 == beam)
        assert(gonio2 == gonio)
        assert(detector2 == detector)

        # Get the models from an index back and check they're not the same
        beam2 = sub_sweep.get_beam(0)
        gonio2 = sub_sweep.get_goniometer(0)
        detector2 = sub_sweep.get_detector(0)
        assert(beam2 != beam)
        assert(gonio2 != gonio)
        assert(detector2 != detector)


class TestImageSetFactory(object):

    def __init__(self):
        pass

    def get_file_list(self):
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
        image_indices = range(1, 10)
        for i in image_indices:
            filenames.append(os.path.join(path, 'centroid_000{0}.cbf'.format(i)))

        return filenames

    def run(self):
        from dxtbx.imageset import ImageSetFactory, ImageSweep

        filenames = self.get_file_list()

        sweep = ImageSetFactory.new(filenames)

        assert(isinstance(sweep[0], ImageSweep) == True)

        print 'OK'


class TestRunner(object):

    def __init__(self):
        pass

    def run(self):

        # Test the multi file state object
        test = TestMultiFileState()
        test.run()

        # Test the sweep file list
        test = TestSweepFileList()
        test.run()

        # Test the multi file reader class
        test = TestMultiFileReader()
        test.run()

        # Test the image set class
        test = TestImageSet()
        test.run()

        # The the sweep class
        test = TestImageSweep()
        test.run()

        # Test the ImageSetFactory class
        test = TestImageSetFactory()
        test.run()


if __name__ == '__main__':
    runner = TestRunner()
    runner.run()
