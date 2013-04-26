
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
    '''A class to keep track of multi file reader state.'''

    def __init__(self, format_class):
        '''Initialise with format class.'''
        self._format_class = format_class
        self._current_format_instance = None

    def format_class(self):
        '''Get the format class.'''
        return self._format_class

    def load_file(self, filename):
        '''Load the file with the given filename.'''

        # Check the current format is the one we need
        if (self.get_format() == None or
            filename != self.get_format().get_image_file()):

            # Read the format instance
            format_instance = self._format_class(filename)

            # Check the format instance is valid
            if not self._is_format_valid(format_instance):
                RuntimeError("Format is invalid.")

            # Set the current format instance
            self._current_format_instance = format_instance

    def get_format(self):
        '''Get the current format instance.'''
        return self._current_format_instance

    def _is_format_valid(self, format_instance):
        '''Check if the format object is valid.'''
        return format_instance.understand(format_instance.get_image_file())


class MultiFileReader(ReaderBase):
    '''A multi file reader class implementing the ReaderBase interface.'''

    def __init__(self, format_class, filenames, formatchecker=None):
        '''Initialise the reader with the format and list of filenames.'''
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
        '''Compare the reader by format class and filename list.'''
        return (self.get_format_class() == other.get_format_class() and
                self.get_image_paths() == other.get_image_paths())

    def get_image_paths(self, indices=None):
        '''Get the list of image paths.'''
        if indices == None:
            return self._filenames
        else:
            return [self._filenames[i] for i in indices]

    def get_format_class(self):
        '''Get the format class.'''
        return self._state.format_class()

    def get_image_size(self):
        '''Get the image size.'''
        return self.get_format().get_detector().get_image_size()

    def get_path(self, index=None):
        '''Get the path the given index.'''
        if index == None:
            return self.get_format().get_image_file()
        else:
            return self._filenames[index]

    def get_detectorbase(self, index=None):
        '''Get the detector base instance at given index.'''
        return self.get_format(index).get_detectorbase()

    def get_detector(self, index=None):
        '''Get the detector instance at given index.'''
        return self.get_format(index).get_detector()

    def get_beam(self, index=None):
        '''Get the beam instance at given index.'''
        return self.get_format(index).get_beam()

    def get_goniometer(self, index=None):
        '''Get the goniometer instance at given index.'''
        return self.get_format(index).get_goniometer()

    def get_scan(self, index=None):
        '''Get the scan instance at given index.'''
        return self.get_format(index).get_scan()

    def read(self, index=None):
        '''Read the image frame at the given index.'''
        return self.get_format(index).get_raw_data()

    def get_format(self, index=None):
        '''Get the format at the given index.'''
        return self._update_state(index).get_format()

    def _update_state(self, index=None):
        '''Update the state and load file at given index.'''
        if index:
            self._state.load_file(self.get_path(index))
        elif self._state.get_format() == None:
            self._state.load_file(self.get_path(0))

        return self._state

    def is_valid(self, indices=None):
        '''Ensure imageset is valid.'''
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
            except IndexError, RuntimeError:
                return False

            # Check the format experimental models
            if not self._is_format_valid(format_instance):
                return False

        # All images valid
        return True


class SweepFileList(object):
    '''Class implementing a file list interface for sweep templates.'''

    def __init__(self, template, array_range):
        '''Initialise the class with the template and array range.'''

        assert(array_range[0] >= 0)
        assert(array_range[0] <= array_range[1])
        self._template = template
        self._array_range = array_range

    def __getitem__(self, index):
        '''Get the filename at that array index.'''
        return self.get_filename(self._array_range[0] + index)

    def __iter__(self):
        '''Iterate through the filenames.'''
        for i in range(len(self)):
            yield self.__getitem__(i)

    def __str__(self):
        '''Get the string representation of the file list.'''
        return str([filename for filename in self])

    def __len__(self):
        '''Get the length of the file list.'''
        return self._array_range[1] - self._array_range[0]

    def __eq__(self, other):
        '''Compare filelist by template and array range.'''
        return (self.template() == other.template() and
                self.array_range() == other.array_range())

    def template(self):
        '''Get the template.'''
        return self._template

    def array_range(self):
        '''Get the array range.'''
        return self._array_range

    def indices(self):
        '''Get the image indices.'''
        return range(*self._array_range)

    def get_filename(self, index):
        '''Get the filename at the given index.'''
        if not self.is_index_in_range(index):
            raise IndexError('Image file index out of range')

        return self._template % (index + 1)

    def is_index_in_range(self, index):
        '''Ensure that the index is within the array range.'''
        return self._array_range[0] <= index < self._array_range[1]
