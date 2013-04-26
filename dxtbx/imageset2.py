
class ReaderBase(object):
    '''The imaegset reader base class.'''

    def __init__(self):
        pass

    def __cmp__(self, other):
        pass

    def get_image_paths(self, indices=None):
        pass

    def get_image_size(self):
        pass

    def get_format(self, index=None):
        pass

    def get_path(self, index=None):
        pass

    def is_valid(self, indices=None):
        pass

    def read_image(self, index=None):
        pass

    def get_detectorbase(self, index=None):
        pass

    def get_detector(self, index=None):
        pass

    def get_beam(self, index=None):
        pass

    def get_goniometer(self, index=None):
        pass

    def get_scan(self, index=None):
        pass


class MultiFileState(object):

    def __init__(self, format_class):
        self._format_class = format_class
        self._current_format_instance = None

    def get_format_class(self):
        return self._format_class

    def load_file(self, filename):
        format_instance = self._format_class(filename)

        if not self._is_format_valid(format_instance):
            RuntimeError("Format is invalid.")

        self._current_format_instance = format_instance

    def get_format(self):
        return self._current_format_instance

    def _is_format_valid(self, format_instance):
        """ Check if the format object is valid. """
        return format_instance.understand(format_instance.get_image_file())


class MultiFileReader(ReaderBase):

    def __init__(self, format_class, filenames, formatchecker=None):
        import os

        # Ensure we have enough images and format has been specified
        assert(format_class != None)
        assert(len(filenames) > 0)

        # Save the image indices
        self._filenames = filenames

        # Handle the state of the MultiFileReader class
        self._state = MultiFileState(format_class)

        # A function object to check formats are valid
        if formatchecker != None:
            self._is_format_valid = formatchecker
        else:
            self._is_format_valid = lambda fmt: True

    def __cmp__(self, other):
        return (self.get_format_class() == other.get_format_class() and
                self.get_image_paths() == other.get_image_paths())

    def get_image_paths(self, indices=None):
        if indices == None:
            return self._filenames
        else:
            return [self._filenames[i] for i in indices]

    def get_format_class(self):
        return self._state.get_format_class()

    def get_image_size(self):
        return self.get_format().get_detector().get_image_size()

    def get_path(self, index=None):
        if index == None:
            return self.get_format().get_image_file()
        else:
            return self._filenames[index]

    def get_detectorbase(self, index=None):
        return self.get_format(index).get_detectorbase()

    def get_detector(self, index=None):
        return self.get_format(index).get_detector()

    def get_beam(self, index=None):
        return self.get_format(index).get_beam()

    def get_goniometer(self, index=None):
        return self.get_format(index).get_goniometer()

    def get_scan(self, index=None):
        return self.get_format(index).get_scan()

    def read(self, index=None):
        return self.get_format(index).get_raw_data()

    def get_format(self, index=None):
        return self._update_state(index).get_format()

    def _update_state(self, index=None):
        if index:
            self._state.load_file(self.get_path(index))
        elif self._state.get_format() == None:
            self._state.load_file(self.get_path(0))
        return self._state

    def is_valid(self, indices=None):
        import os

        # If no indices, get indices of all filenames
        if indices == None:
            indices = range(len(self._filenames))

        # Loop through all the images
        for index in indices:

            # Read and try to cache the format, if this fails, the
            # format is invalid, so return false.
            try:
                format_instance = self.get_format(index)
            except RuntimeError:
                return False

            # Check the format experimental models
            if not self._is_format_valid(format_instance):
                return False

        # All images valid
        return True


class SweepFileList(object):
    def __init__(self, template, array_range):
        self._template = template
        self._array_range = array_range

    def __getitem__(self, index):
        return self._get_filename(index)

    def __iter__(self):
        for i in self.indices():
            yield self._get_filename(i)

    def __str__(self):
        return str([filename for filename in self])

    def __len__(self):
        return self._array_range[1] - self._array_range[0]

    def __cmp__(self, other):
        return (self.get_template() == other.get_template() and
                self.get_array_range() == other.get_array_range())

    def get_template(self):
        return self._template

    def array_range(self):
        return self._array_range

    def indices(self):
        return range(*self._array_range)

    def get_filename(self, index):
        if not self._is_index_in_range(index):
            raise IndexError('Image file index out of range')

        return self._template % (index + 1)

    def is_index_in_range(self, index):
        return self._array_range[0] <= index < self._array_range[1]
