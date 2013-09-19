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

        # Save some overidden models. If these models are not set, then
        # getting the models defers to the reader. Otherwise, if they are
        # set and the get functions are called without an index then
        # these models will be returned
        self._beam = None
        self._detector = None
        self._goniometer = None

    def __cmp__(self, other):
        pass

    def get_image_paths(self, indices=None):
        pass

    def get_image_size(self, panel=0):
        pass

    def get_format(self, index=None):
        pass

    def get_format_class(self, index=None):
        pass

    def get_path(self, index=None):
        pass

    def is_valid(self, indices=None):
        pass

    def read(self, index=None):
        pass

    def get_detectorbase(self, index=None):
        pass

    def get_detector(self, index=None):
        if index is None and self._detector is not None:
            return self._detector

        return self.read_detector(index)

    def set_detector(self, detector):

        # Get the current detector
        curr_detector = self.get_detector()

        # Ensure the same number of panels
        assert(len(curr_detector) == len(detector))

        # Loop through all the panels
        for i in range(len(curr_detector)):

            # Override the geometry in the detector
            curr_detector[i].set_frame(
                detector[i].get_fast_axis(),
                detector[i].get_slow_axis(),
                detector[i].get_origin())

        # Set the new detector
        self._detector = curr_detector

    def get_beam(self, index=None):
        if index is None and self._beam is not None:
            return self._beam

        return self.read_beam(index)

    def set_beam(self, beam):

        # Get the current beam
        curr_beam = self.get_beam()

        # Set beam geometry
        curr_beam.set_direction(beam.get_direction())
        curr_beam.set_divergence(beam.get_divergence())
        curr_beam.set_sigma_divergence(beam.get_sigma_divergence())

        # Set the new beam
        self._beam = curr_beam

    def get_goniometer(self, index=None):
        if index is None and self._goniometer is not None:
            return self._goniometer

        return self.read_goniometer(index)

    def set_goniometer(self, goniometer):

        # Get the current goniometer
        curr_goniometer = self.get_goniometer()

        # Set goniometer geometry
        curr_goniometer.set_rotation_axis(goniometer.get_rotation_axis())
        curr_goniometer.set_fixed_rotation(goniometer.get_fixed_rotation())

        # Set the new goniometer
        self._goniometer = curr_goniometer

    def get_scan(self, index=None):
        pass

    def read_detector(self, index=None):
        pass

    def read_goniometer(self, index=None):
        pass

    def read_beam(self, index=None):
        pass

class SingleFileReader(ReaderBase):
    '''The single file reader class.'''

    def __init__(self, format_instance):
        '''Initialise the reader class.'''
        ReaderBase.__init__(self)

        # Set the format instance
        self._format = format_instance

    def __cmp__(self, other):
        '''Compare the reader to another reader.'''
        return self._format == other._format

    def get_image_paths(self, indices=None):
        '''Get the image paths within the file.'''

        # Get paths for each file
        filenames = [self._format.get_image_file(i)
            for i in range(self._format.get_num_images())]

        # Return within the given range
        if indices == None:
            return filenames
        else:
            return [filenames[i] for i in indices]

    def get_image_size(self, panel=0):
        '''Get the image size from the detector.'''
        return self._format.get_detector()[panel].get_image_size()

    def get_format(self, index=None):
        '''Get the format instance'''
        return self._format

    def get_format_class(self, index=None):
        '''Get the format class'''
        return self._format.__class__

    def get_path(self, index=None):
        '''Get the image file for the given index.'''
        return self._format.get_image_file(index)

    def is_valid(self, indices=None):
        '''Ensure the reader is valid.'''
        return True

    def read(self, index=None):
        '''Get the image data.'''
        return self._format.get_raw_data(index)

    def get_detectorbase(self, index=None):
        '''Get the detector base instance.'''
        return self._format.get_detectorbase(index)

    def read_detector(self, index=None):
        '''Get the detector instance.'''
        return self._format.get_detector(index)

    def read_beam(self, index=None):
        '''Get the beam instance.'''
        return self._format.get_beam(index)

    def read_goniometer(self, index=None):
        '''Get the goniometer instance.'''
        return self._format.get_goniometer(index)

    def get_scan(self, index=None):
        '''Get the scan instance.'''
        if isinstance(index, tuple):
            scan = self._format.get_scan(index[0])
            scan.set_image_range((index[0] + 1, index[1]))
        else:
            scan = self._format.get_scan(index)

        return scan


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
        ReaderBase.__init__(self)

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

    def get_image_size(self, panel=0):
        '''Get the image size.'''
        return self.get_format().get_detector()[panel].get_image_size()

    def get_path(self, index=None):
        '''Get the path the given index.'''
        if index == None:
            return self.get_format().get_image_file()
        else:
            return self._filenames[index]

    def get_detectorbase(self, index=None):
        '''Get the detector base instance at given index.'''
        return self.get_format(index).get_detectorbase()

    def read_detector(self, index=None):
        '''Get the detector instance at given index.'''
        return self.get_format(index).get_detector()

    def read_beam(self, index=None):
        '''Get the beam instance at given index.'''
        return self.get_format(index).get_beam()

    def read_goniometer(self, index=None):
        '''Get the goniometer instance at given index.'''
        return self.get_format(index).get_goniometer()

    def get_scan(self, index=None):
        '''Get the scan instance at given index.'''
        if isinstance(index, tuple):
            scan = self.get_format(index[0]).get_scan()
            scan.set_image_range((index[0] + 1, index[1]))
        else:
            scan = self.get_format(index).get_scan()

        return scan

    def read(self, index=None):
        '''Read the image frame at the given index.'''

        # Get the format instance and the number of panels
        format_instance = self.get_format(index)
        npanels = len(format_instance.get_detector())

        # Return a flex array for single panels and a tuple of flex arrays
        # for multiple panels
        assert(npanels > 0)
        if npanels == 1:
            return format_instance.get_raw_data()
        else:
            return tuple([format_instance.get_raw_data(i) for i in range(npanels)])

    def get_format(self, index=None):
        '''Get the format at the given index.'''
        return self._update_state(index).get_format()

    def _update_state(self, index=None):
        '''Update the state and load file at given index.'''
        if index is not None:
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

    def set_detector(self, detector):
        ''' Set the detector model.'''
        self.reader().set_detector(detector)

    def get_beam(self, index=None):
        ''' Get the beam. '''
        return self.reader().get_beam(self._image_index(index))

    def set_beam(self, beam):
        ''' Set the beam model.'''
        self.reader().set_beam(beam)

    def get_image_size(self, index=0):
        ''' Get the image size. '''
        return self.reader().get_image_size(index)

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

        #assert(array_range[0] >= 0)
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


class ImageSweep(ImageSet):
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
            return ImageSweep(self.reader(), indices)
        else:
            return self.reader().read(self._indices[item])

    def get_template(self):
        ''' Get the template. '''
        from dxtbx.sweep_filenames import template_format_to_string
        return template_format_to_string(self.reader()._filenames.template())

    def get_array_range(self):
        ''' Get the array range. '''
        return (self._indices[0], self._indices[-1] + 1)

    def get_goniometer(self, index=None):
        ''' Get goniometer, '''
        return self.reader().get_goniometer(self._image_index(index))

    def set_goniometer(self, goniometer):
        ''' Set the goniometer model '''
        self.reader().set_goniometer(goniometer)

    def get_scan(self, index=None):
        ''' Get the scan. '''
        if index == None:
            index = self.get_array_range()
        else:
            index = self._image_index(index)

        return self.reader().get_scan(index)

    def to_array(self, item=None, panel=0):
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
        if item is None:
            return self._to_array_all(panel)
        else:
            return self._to_array_w_range(item, panel)

    def _to_array_all(self, panel=0):
        ''' Get the array from all the sweep elements. '''

        # FIXME this should be using flex arrays not numpy ones as we want
        # to be able to pass the data to cctbx C++ code...

        from scitbx.array_family import flex

        # Get the image dimensions
        size_z = len(self)
        size_y = self.reader().get_image_size(panel)[1]
        size_x = self.reader().get_image_size(panel)[0]

        # Check sizes are valid
        if size_z <= 0 or size_y <= 0 or size_x <= 0:
            raise RuntimeError("Invalid dimensions")

        # Allocate the array
        array = flex.int(flex.grid(size_z, size_y, size_x))

        # Loop through all the images and set the image data
        for k, image in enumerate(self):
            if not isinstance(image, tuple):
                image = (image,)
            im = image[panel]
            im.reshape(flex.grid(1, *im.all()))
            array[k:k+1,:,:] = im

        # Return the array
        return array

    def _to_array_w_range(self, item, panel=0):
        ''' Get the array from the user specified range. '''
        from scitbx.array_family import flex

        # Get the range from the given index item
        z0, z1, y0, y1, x0, x1 = self._get_data_range(item, panel)

        # Get the image dimensions
        size_z = z1 - z0
        size_y = y1 - y0
        size_x = x1 - x0

        # Check sizes are valid
        if size_z <= 0 or size_y <= 0 or size_x <= 0:
            raise RuntimeError("Invalid dimensions")

        # Allocate the array
        array = flex.int(flex.grid(size_z, size_y, size_x))

        # Loop through all the images and set the image data
        for k, index in enumerate(self.indices()[z0:z1]):
            image = self.reader().read(index)
            if not isinstance(image, tuple):
                image = (image,)
            im = image[panel]
            im.reshape(flex.grid(1, *im.all()))
            array[k:k+1,:,:] = im[0:1, y0:y1, x0:x1]

        # Return the array
        return array

    def _get_data_range(self, item, panel=0):
        ''' Get the range from the user specified index item. '''

        # Ensure item is a tuple
        if isinstance(item, tuple):

            # Just the range of images given
            if len(item) == 2:
                z0, z1 = item
                y0, y1 = (0, self.reader().get_image_size(panel)[1])
                x0, x1 = (0, self.reader().get_image_size(panel)[0])
                return self._truncate_range((z0, z1, y0, y1, x0, x1))

            # The range in each direction given
            elif len(item) == 6:
                return self._truncate_range(item, panel)

        # Raise index error
        raise IndexError("bad index")

    def _truncate_range(self, data_range, panel=0):
        ''' Truncate the range to the size of available data. '''

        # Get items from range
        z0, z1, y0, y1, x0, x1 = data_range

        # Get the number of frames and image size
        size_z = len(self)
        size_x, size_y = self.reader().get_image_size(panel)

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
            if self._indices_sequential_ge_zero(indices):
                return True
            else:
                return False

    def _indices_sequential_ge_zero(self, indices):
        ''' Determine if indices are sequential.'''
        prev = indices[0]
        if prev < 0:
            return False
        for curr in indices[1:]:
            if curr != prev + 1:
                return False
            prev = curr

        return True


class ImageSetFactory(object):
    ''' Factory to create imagesets and sweeps. '''

    @staticmethod
    def new(filenames, check_headers=False, ignore_unknown=False):
        ''' Create an imageset or sweep

        Params:
            filenames A list of filenames
            check_headers Check the headers to ensure all images are valid
            ignore_unknown Ignore unknown formats

        Returns:
            A list of imagesets

        '''
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
        imagesetlist = []
        for filelist in filelist_per_imageset:
            try:
                imagesetlist.append(ImageSetFactory._create_imageset_or_sweep(
                    filelist, check_headers))
            except Exception, e:
                if not ignore_unknown:
                    raise e

        # Return the imageset list
        return imagesetlist

    @staticmethod
    def from_template(template, image_range=None, check_headers=False):
        '''Create a new sweep from a template.

        Params:
            template The template argument
            image_range The image range
            check_headers Check the headers to ensure all images are valid

        Returns:
            A list of sweeps

        '''
        import os
        from dxtbx.format.Registry import Registry
        from dxtbx.sweep_filenames import template_image_range

        # Check the template is valid
        if template.count('#') < 1:
            raise ValueError("Invalid template")

        # Get the template format
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)

        # Get the template image range
        if image_range is None:
            image_range = template_image_range(template)


        # Set the image range
        array_range = (image_range[0] - 1, image_range[1])

        # Create the sweep file list
        filenames = SweepFileList(template_format, array_range)

        # Get the format class
        format_class = Registry.find(filenames[0])

        # Create the sweep object
        sweep = ImageSweep(MultiFileReader(format_class, filenames))

        # Check the sweep is valid
        if check_headers and not sweep.is_valid():
            raise RuntimeError('Invalid sweep of images')

        # Return the sweep
        return [sweep]

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
        from format.FormatMultiImage import FormatMultiImage
        if issubclass(format_class, FormatMultiImage):
            assert len(filenames) == 1
            format_instance = format_class(filenames[0])
            image_set = ImageSet(SingleFileReader(format_instance))
        else:
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

        # Set the image range
        array_range = (min(indices) - 1, max(indices))

        # Create the sweep file list
        filenames = SweepFileList(template_format, array_range)

        # Create the sweep object
        sweep = ImageSweep(MultiFileReader(format_class, filenames))

        # Check the sweep is valid
        if check_headers and not sweep.is_valid():
            raise RuntimeError('Invalid sweep of images')

        # Return the sweep
        return sweep
