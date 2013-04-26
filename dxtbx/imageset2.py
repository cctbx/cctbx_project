#!/usr/bin/env python
#
# imageset.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


class ReaderBase(object):
    '''The imageset reader base class.'''

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


class SingleFileReader(ReaderBase):
    '''The single file reader class.'''

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
            return list(self._filenames)
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


class ImageSet(object):
    ''' A class exposing the external image set interface. '''

    def __init__(self, reader, indices=None):
        ''' Initialise the ImageSet object.

        Params:
            reader The reader object
            array_range The image range (first, last)

        '''
        # If no reader is set then throw an exception
        if not reader:
            raise ValueError("ImageSet needs a reader!")

        # Set the reader
        self._reader = reader

        # Set the array range or get the range from the reader
        if indices:
            self._indices = indices
        else:
            self._indices = range(len(self.reader().get_image_paths()))

    def __getitem__(self, item):
        ''' Get an item from the image set stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new ImageSet object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new ImageSet object

        '''
        if isinstance(item, slice):
            indices = self._indices[item]
            return ImageSet(self.reader(), indices)
        else:
            return self.reader().read(self._indices[item])

    def __len__(self):
        ''' Return the number of images in this image set. '''
        return len(self._indices)

    def __str__(self):
        ''' Return the array indices of the image set as a string. '''
        return str(self.paths())

    def __iter__(self):
        ''' Iterate over the array indices and read each image in turn. '''
        for f in self._indices:
            yield self.reader().read(f)

    def __cmp__(self, other):
        ''' Compare this image set to another. '''
        return self.reader() == other.reader()

    def indices(self):
        ''' Return the indices in the image set. '''
        return list(self._indices)

    def paths(self):
        ''' Return a list of filenames referenced by this set. '''
        filenames = self.reader().get_image_paths()
        return [filenames[i] for i in self._indices]

    def is_valid(self):
        ''' Validate all the images in the image set. Can take a long time. '''
        return self.reader().is_valid(self._indices)

    def get_detector(self, index=None):
        ''' Get the detector. '''
        return self.reader().get_detector(self._image_index(index))

    def get_beam(self, index=None):
        ''' Get the beam. '''
        return self.reader().get_beam(self._image_index(index))

    def get_image_size(self):
        ''' Get the image size. '''
        return self.reader().get_image_size()

    def get_detectorbase(self, index=None):
        ''' Get the detector base instance for the given index. '''
        return self.reader().get_detectorbase(self._image_index(index))

    def reader(self):
        ''' Return the image set reader. '''
        return self._reader

    def _image_index(self, index=None):
        ''' Convert image set index to image index.'''
        if index == None:
            return None
        elif index < 0 or index >= len(self._indices):
            raise IndexError('Index out of range')
        return self._indices[index]


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


class Sweep(ImageSet):
    ''' A class exposing the external sweep interface. '''

    def __init__(self, reader, indices=None):
        ImageSet.__init__(self, reader, indices)

    def __getitem__(self, item):
        ''' Get an item from the sweep stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new Sweep object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new Sweep object

        '''
        if isinstance(item, slice):
            if item.step != None:
                raise IndexError('Sweeps must be sequential')
            indices = self._indices[item]
            return Sweep(self.reader(), indices)
        else:
            return self.reader().read(self._indices[item])

    def get_array_range(self):
        ''' Get the array range. '''
        return (self._indices[0], self._indices[-1] + 1)

    def get_goniometer(self, index=None):
        ''' Get goniometer, '''
        return self.reader().get_goniometer(self._image_index(index))

    def get_scan(self, index=None):
        ''' Get the scan. '''
        return self.reader().get_scan(self._image_index(index))

    def to_array(self, item=None):
        ''' Read all the files in the sweep and convert them into an array
        of the appropriate dimensions.

        The required array is allocated first, this has he useful property
        that if you've been lazy and just got a sweep and extracted the
        array without consideration for the amount of memory available on
        your machine, you'll get an exception straight away.

        TODO:
            Currently uses numpy for fast copying of arrays, try to do
            this using flex arrays.

        Params:
            item The index item (frame 0, frame n), (z0, z1, y0, y1, x0, x1)

        Returns:
            The sweep image data as an array.

        '''
        if not item:
            return self._to_array_all()
        else:
            return self._to_array_w_range(item)

    def _to_array_all(self):
        ''' Get the array from all the sweep elements. '''

        # FIXME this should be using flex arrays not numpy ones as we want
        # to be able to pass the data to cctbx C++ code...

        from scitbx.array_family import flex
        from numpy import zeros, int32

        # Get the image dimensions
        size_z = len(self)
        size_y = self.reader().get_image_size()[1]
        size_x = self.reader().get_image_size()[0]

        # Check sizes are valid
        if size_z <= 0 or size_y <= 0 or size_x <= 0:
            raise RuntimeError("Invalid dimensions")

        # Allocate the array
        array = zeros(shape=(size_z, size_y, size_x), dtype=int32)

        # Loop through all the images and set the image data
        for k, image in enumerate(self):
            array[k,:,:] = image.as_numpy_array()

        # Return the array
        return flex.int(array)

    def _to_array_w_range(self, item):
        ''' Get the array from the user specified range. '''
        from scitbx.array_family import flex
        from numpy import zeros, int32

        # Get the range from the given index item
        z0, z1, y0, y1, x0, x1 = self._get_data_range(item)

        # Get the image dimensions
        size_z = z1 - z0
        size_y = y1 - y0
        size_x = x1 - x0

        # Check sizes are valid
        if size_z <= 0 or size_y <= 0 or size_x <= 0:
            raise RuntimeError("Invalid dimensions")

        # Allocate the array
        array = zeros(shape=(size_z, size_y, size_x), dtype=int32)

        # Loop through all the images and set the image data
        for k, index in enumerate(self.indices()[z0:z1]):
            image = self.reader().read(index)
            array[k,:,:] = image.as_numpy_array()[y0:y1, x0:x1]

        # Return the array
        return flex.int(array)

    def _get_data_range(self, item):
        ''' Get the range from the user specified index item. '''

        # Ensure item is a tuple
        if isinstance(item, tuple):

            # Just the range of images given
            if len(item) == 2:
                z0, z1 = item
                y0, y1 = (0, self.reader().get_image_size()[1])
                x0, x1 = (0, self.reader().get_image_size()[0])
                return self._truncate_range((z0, z1, y0, y1, x0, x1))

            # The range in each direction given
            elif len(item) == 6:
                return self._truncate_range(item)

        # Raise index error
        raise IndexError("bad index")

    def _truncate_range(self, data_range):
        ''' Truncate the range to the size of available data. '''

        # Get items from range
        z0, z1, y0, y1, x0, x1 = data_range

        # Get the number of frames and image size
        size_z = len(self)
        size_x, size_y = self.reader().get_image_size()

        # Ensure range is valid
        if z0 < 0: z0 = 0
        if y0 < 0: y0 = 0
        if x0 < 0: x0 = 0
        if z1 > size_z: z1 = size_z
        if y1 > size_y: y1 = size_y
        if x1 > size_x: x1 = size_x

        # Return truncated range
        return (z0, z1, y0, y1, x0, x1)


class FilenameAnalyser(object):
    '''Group images by filename into image sets.'''

    def __init__(self):
        '''Initialise the class.'''
        pass

    def __call__(self, filenames):
        '''Group the filenames by imageset.

        Params:
            filenames The list of filenames

        Returns:
            A list of (template, [indices], is_sweep)

        '''
        from sweep_filenames import group_files_by_imageset

        # Analyse filenames to figure out how many imagesets we have
        filelist_per_imageset = group_files_by_imageset(filenames)

        # Label each group as either an imageset or a sweep.
        file_groups = []
        for template, indices in filelist_per_imageset.iteritems():

            # Check if this imageset is a sweep
            is_sweep = self._is_imageset_a_sweep(template, indices)

            # Append the items to the group list
            file_groups.append((template, indices, is_sweep))

        # Return the groups of files
        return file_groups

    def _is_imageset_a_sweep(self, template, indices):
        ''' Return True/False if the imageset is a sweep or not.

        Where more than 1 image that follow sequential numbers are given
        the images are catagorised as belonging to a sweep, otherwise they
        belong to an image set.

        '''
        if len(indices) <= 1:
            return False
        else:
            indices = sorted(indices)
            if self._indices_sequential(indices):
                return True
            else:
                return False

    def _indices_sequential(self, indices):
        ''' Determine if indices are sequential.'''
        prev = indices[0]
        for curr in indices[1:]:
            if curr != prev + 1:
                return False
            prev = curr

        return True


class ImageSetFactory(object):
    ''' Factory to create imagesets and sweeps. '''

    @staticmethod
    def new(filenames, check_headers=False):
        ''' Create an imageset or sweep

        Params:
            filenames A list of filenames
            check_headers Check the headers to ensure all images are valid

        Returns:
            A list of imagesets

        '''
        from sweep_filenames import find_matching_images

        # Ensure we have enough images
        if isinstance(filenames, list):
            assert(len(filenames) > 0)
        elif isinstance(filenames, str):
            filenames = [filenames]
        else:
            raise RuntimeError, 'unknown argument passed to ImageSetFactory'

        # Analyse the filenames and group the images into imagesets.
        analyse_files = FilenameAnalyser()
        filelist_per_imageset = analyse_files(filenames)

        # For each file list denoting an image set, create the imageset
        # and return as a list of imagesets. N.B sweeps and image sets are
        # returned in the same list.
        return [ImageSetFactory._create_imageset_or_sweep(
            filelist, check_headers) for filelist in filelist_per_imageset]

    @staticmethod
    def _create_imageset_or_sweep(filelist, check_headers):
        '''Create either an imageset of sweep.'''
        if filelist[2] == True:
            return ImageSetFactory._create_sweep(filelist, check_headers)
        else:
            return ImageSetFactory._create_imageset(filelist, check_headers)

    @staticmethod
    def _create_imageset(filelist, check_headers):
        '''Create an image set'''
        from dxtbx.format.Registry import Registry

        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count('#')
        if count > 0:
            pfx = template.split('#')[0]
            sfx = template.split('#')[-1]
            template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = Registry.find(filenames[0])

        # Create the image set object
        image_set = ImageSet(MultiFileReader(format_class, filenames))

        # Check the image set is valid
        if check_headers and not image_set.is_valid():
            raise RuntimeError('Invalid ImageSet')

        # Return the image set
        return image_set

    @staticmethod
    def _create_sweep(filelist, check_headers):
        '''Create a sweep'''
        import os
        from dxtbx.format.Registry import Registry

        # Extract info from filelist
        template, indices, is_sweep = filelist

        # Get the template format
        count = template.count('#')
        if count > 0:
            pfx = template.split('#')[0]
            sfx = template.split('#')[-1]
            template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
            filenames = [template_format % index for index in indices]
        else:
            filenames = [template]

        # Sort the filenames
        filenames = sorted(filenames)

        # Get the format object
        format_class = Registry.find(filenames[0])

        # Get the first image and our understanding
        first_image = filenames[0]

        # Get the directory and first filename and set the template format
        directory, first_image_name = os.path.split(first_image)
        first_image_number = indices[0]

        # Get the template format
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
        template = os.path.join(directory, template_format)

        # Set the image range
        array_range = (min(indices) - 1, max(indices))

        # Create the sweep file list
        filenames = SweepFileList(template, array_range)

        # Create the sweep object
        sweep = Sweep(MultiFileReader(format_class, filenames))

        # Check the sweep is valid
        if check_headers and not sweep.is_valid():
            raise RuntimeError('Invalid sweep of images')

        # Return the sweep
        return sweep
