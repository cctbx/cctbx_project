#!/usr/bin/env python
# dxtbx/sweep.py
#
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Top level interface to the dxtbx: treatment of a sweep of images as a single
# unit, incorporation of the registry and so on... N.B. factory functions will
# be defined which can look up information from the file system etc.
#
# The classes in here should not be instantiated directly as they will depend
# on the form of the data (a single file, a sequence of images) instead
# accessed *only* from the factories.

from __future__ import division

class SweepReader(object):
    '''Definition for a sweep of images, defined to be a set of diffraction
    images in files with matching templates, or a volume formatted file.'''

    def __init__(self, format_class, template_format, image_range):
        """ Construct a sweep object from the given format.

        Params:
            format_class The format class object
            template_format The image template format
            image_range The range of the images in the full sweep

        """
        # Ensure we have enough images and format has been specified
        assert(format_class != None)

        # Set the internal format class
        self._format = format_class

        # Set the image range
        self._template_format = template_format

        # Get the first format instance and get the models
        first_filename = self.get_filename(image_range[0])
        format_instance = self._format(first_filename)
        self._understanding = self._format.understand(first_filename)
        self._detector = format_instance.get_detector()
        self._goniometer = format_instance.get_goniometer()
        self._beam = format_instance.get_beam()
        self._scan = format_instance.get_scan()
        self._scan.set_image_range(image_range)

        # Create the format cache and cache the first format
        self._format_cache = {}
        self._format_cache[0] = format_instance

    def __cmp__(self, other):
        """ Compare this reader to another. """
        return (self.get_template() == other.get_template() and
                self.get_format() == other.get_format())

    def get_format(self):
        """ Get the format instance. """
        return self._format

    def get_template(self):
        """ Get the template. """
        return self._template_format

    def get_filename(self, image):
        """ Get the filename for the given image. """
        return self._template_format % image

    def is_valid(self, array_range=None):
        """ Check that all the images in the given range are valid.

        Params:
            image_range The range of images to check (default all)

        Returns:
            True/False the images are valid

        """
        import os

        # If no range set then use all
        if array_range:
            image_range = map(lambda x: x + 1, array_range)
        else:
            image_range = self._scan.get_image_range()

        # Loop through all the images
        for image in range(*image_range):

            # If format is in cache, it's valid
            try:
                self._format_cache[image]
            except KeyError:

                # Read and try to cache the format, if this fails, the
                # format is invalid, so return false.
                try:
                    self._read_and_cache_format(image)
                except RuntimeError:
                    return False

        # All images valid
        return True

    def read(self, array_index):
        """ Read an image at the given array index. """

        # Convert from array to image numbering and ensure is valid
        image_index = array_index + 1
        if not self._scan.is_image_index_valid(image_index):
            raise IndexError, 'array index out of range'

        # Try to read format from cache, otherwise read from file
        try:
            format_instance = self._format_cache[image_index]
        except KeyError:
            format_instance = self._read_and_cache_format(image_index)

        # Return the raw data
        return format_instance.get_raw_data()

    def get_array_range(self):
        '''Get the useful array range for the [low, high) limits for what
        can be called as sweep[j].'''
        return self._scan.get_array_range()

    def get_image_size(self):
        """ Get the image size. """
        return self._detector.get_image_size()

    def get_detector(self):
        """ Get the detector model. """
        return self._detector

    def get_goniometer(self):
        """ Get the goniometer model. """
        return self._goniometer

    def get_beam(self):
        """ Get the beam model. """
        return self._beam

    def get_scan(self):
        """ Get the scan model. """
        return self._scan

    def _read_and_cache_format(self, image_index):
        """ Read and cache the format object at the given image index. """

        # Get the format instance
        format_instance = self._format(self.get_filename(image_index))

        # If the format instance is valid then cache
        if self._is_format_valid(format_instance):
            pass
            self._format_cache[image_index] = format_instance
        else:
            RuntimeError("Format is invalid.")

        # Return the format instance
        return format_instance

    def _is_format_valid(self, format_instance):
        """ Check if the format object is valid. """

        # Get the image filename
        filename = format_instance.get_image_file()

        # Check we understand in the same way
        if self._format.understand(filename) != self._understanding:
            return False

        # Check that the models are all the same
        if format_instance.get_beam() != self._beam:
            return False
        if format_instance.get_goniometer() != self._goniometer:
            return False
        if format_instance.get_detector() != self._detector:
            return False

        # Format is valid
        return True


class BufferedSweepReader(SweepReader):
    """ Class to provide cached reading of images. """

    def __init__(self, format_class, template_format, image_range, max_cache=None):
        """ Initialise the SweepBuffer class.

        Params:
            reader The reader class that does the actual reading
            max_cache The maximum images to cache

        """
        from collections import OrderedDict

        # Initialise the base class
        SweepReader.__init__(self, format_class, template_format, image_range)

        # Set the maximum images to cache
        if max_cache:
            self._max_cache = max_cache
        else:
            self._max_cache = 20

        # Create the image cache
        self._image_cache = OrderedDict()

    def set_max_cache(self, max_cache):
        """ Set the maximum cache size. """
        self._max_cache = max_cache

    def get_max_cache(self):
        """ Get the maximum cache size. """
        return self._max_cache

    def read(self, index):
        """ Read an image. Try to find image in cache, otherwise read a
        new image and cache.

        Params:
            index The array index

        Returns:
            The raw image data

        """
        try:
            return self._image_cache[index]
        except KeyError:
            return self._read_image_and_cache(index)

    def cached(self):
        """ Return the array indices currently in the cache. """
        return self._image_cache.keys()

    def _read_image_and_cache(self, index):
        """ Read an image an insert into the cache.

        Params:
            index The array index

        Returns:
            The raw image data.

        """
        # Read an image from file
        image = SweepReader.read(self, index)

        # If the size of the cache is > max cache then remove an image
        # from the cache in a FIFO manner.
        if len(self._image_cache) >= self._max_cache:
            self._image_cache.pop(self._image_cache.keys()[0])

        # Add the image to the cache
        self._image_cache[index] = image

        # Return the image image
        return image


class Sweep(object):
    """ A class exposing the external sweep interface. """

    def __init__(self, reader, array_range=None):
        """ Initialise the Sweep object.

        Params:
            reader The reader object
            array_range The image range (first, last)

        """
        # If no reader is set then throw an exception
        if not reader:
            raise ValueError("Sweep needs a reader!")

        # Set the reader
        self._reader = reader

        # Set the array range or get the range from the reader
        if array_range:
            if len(array_range) == 2:
                self._array_range = array_range
            else:
                raise ValueError("array_range should have two elements")
        else:
            self._array_range = self._reader.get_array_range()

    def __getitem__(self, item):
        """ Get an item from the sweep stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new Sweep object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new Sweep object

        """
        if isinstance(item, slice):
            indices = self.indices()[item]
            if len(indices) == 0:
                raise IndexError("Bad slice")
            return Sweep(self._reader, (indices[0], indices[-1] + 1))
        else:
            return self._reader.read(self.indices()[item])

    def __len__(self):
        """ Return the number of images in this sweep. """
        return self._array_range[1] - self._array_range[0]

    def __str__(self):
        """ Return the array indices of the sweep as a string. """
        return str(range(*self._array_range))

    def __iter__(self):
        """ Iterate over the array indices and read each image in turn. """
        for f in self.indices():
            yield self._reader.read(f)

    def __cmp__(self, other):
        """ Compare this sweep to another. """
        return self.reader() == other.reader()

    def get_array_range(self):
        """ Get the array range. """
        return self._array_range

    def indices(self):
        """ Return the array indices in the sweep. """
        return range(*self._array_range)

    def to_array(self, item=None):
        """ Read all the files in the sweep and convert them into an array
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

        """
        if not item:
            return self._to_array_all()
        else:
            return self._to_array_w_range(item)

    def is_valid(self):
        """ Validate all the images in the sweep. Can take a long time. """
        return self._reader.is_valid(self._array_range)

    def get_detector(self):
        """ Get the detector. """
        return self._reader.get_detector()

    def get_goniometer(self):
        """ Get the goniometer. """
        return self._reader.get_goniometer()

    def get_beam(self):
        """ Get the beam. """
        return self._reader.get_beam()

    def get_scan(self):
        """ Get the scan. """
        return self._reader.get_scan()

    def reader(self):
        """ Return the sweep reader. """
        return self._reader

    def get_image_size(self):
        """ Get the image size. """
        return self._reader.get_image_size()

    def _to_array_all(self):
        """ Get the array from all the sweep elements. """

        # FIXME this should be using flex arrays not numpy ones as we want
        # to be able to pass the data to cctbx C++ code...

        from scitbx.array_family import flex
        from numpy import zeros, int32

        # Get the image dimensions
        size_z = len(self)
        size_y = self._reader.get_image_size()[1]
        size_x = self._reader.get_image_size()[0]

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
        """ Get the array from the user specified range. """
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
            image = self._reader.read(index)
            array[k,:,:] = image.as_numpy_array()[y0:y1, x0:x1]

        # Return the array
        return flex.int(array)

    def _get_data_range(self, item):
        """ Get the range from the user specified index item. """

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
        """ Truncate the range to the size of available data. """

        # Get items from range
        z0, z1, y0, y1, x0, x1 = data_range

        # Get the number of frames and image size
        size_z = len(self)
        size_x, size_y = self._reader.get_image_size()

        # Ensure range is valid
        if z0 < 0: z0 = 0
        if y0 < 0: y0 = 0
        if x0 < 0: x0 = 0
        if z1 > size_z: z1 = size_z
        if y1 > size_y: y1 = size_y
        if x1 > size_x: x1 = size_x

        # Return truncated range
        return (z0, z1, y0, y1, x0, x1)

class SweepFactory:
    """ The factory class to create a sweep object. """

    @staticmethod
    def sweep(argument, check_headers = False):
        import os
        from dxtbx.format.Registry import Registry
        from dxtbx.sweep_filenames import template_regex

        if type(argument) == type([]):
            filenames = argument

        elif type(argument) == type('str'):
            from dxtbx.sweep_filenames import find_matching_images
            filenames = find_matching_images(argument)
        else:
            raise RuntimeError, 'unknown argument passed to sweep factory'

        # Ensure we have enough images and format has been specified
        assert(len(filenames) > 0)

        # Get the format object
        filenames = sorted(filenames)
        format_class = Registry.find(filenames[0])

        # Get the first image and our understanding
        first_image = filenames[0]

        # Get the directory and first filename and set the template format
        directory, first_image_name = os.path.split(first_image)
        template, first_image_number = template_regex(first_image_name)
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
        template_format = os.path.join(directory, template_format)

        # Check the input filenames are valid
        image_numbers = [first_image_number]
        for j, image in enumerate(filenames[1:]):
            path, image_name = os.path.split(image)
            assert(path == directory)
            image_number = int(image_name.replace(pfx, '').replace(sfx, ''))
            assert(image_number == first_image_number + 1 + j)
            image_numbers.append(image_number)

        # Set the image range
        image_range = (min(image_numbers), max(image_numbers))

        # Create the sweep object
        sweep = Sweep(BufferedSweepReader(format_class,
            template_format, image_range))

        # Check the sweep is valid
        if check_headers and not sweep.is_valid():
            raise RuntimeError('Invalid sweep of images')

        # Return the sweep
        return sweep
