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
        from dxtbx.imageset2 import MultiFileState
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
        from dxtbx.imageset2 import SweepFileList

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
        except:
            pass
        try:
            filename = filelist.get_filename(too_high)
            assert(False)
        except:
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
        except:
            pass
        try:
            filename = filelist[too_high]
            assert(False)
        except:
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
        from dxtbx.imageset2 import SweepFileList

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
        from dxtbx.imageset2 import SweepFileList

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

        from dxtbx.imageset2 import MultiFileReader
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

if __name__ == '__main__':
    runner = TestRunner()
    runner.run()
