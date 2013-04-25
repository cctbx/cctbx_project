#!/usr/bin/env python
# dxtbx/image_set.py
#
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Top level interface to the dxtbx: treatment of a image set of images as a single
# unit, incorporation of the registry and so on... N.B. factory functions will
# be defined which can look up information from the file system etc.
#
# The classes in here should not be instantiated directly as they will depend
# on the form of the data (a single file, a sequence of images) instead
# accessed *only* from the factories.

from __future__ import division

class ImageSetReader(object):
    '''Definition for a image set of images, defined to be a set of diffraction
    images in files with matching templates, or a volume formatted file.'''

    def __init__(self, format_class, template, image_indices):
        """ Construct a image set object from the given format.

        Params:
            format_class The format class object
            template_format The image template format
            image_range The range of the images in the full image set

        """
        import os

        # Ensure we have enough images and format has been specified
        assert(format_class != None)

        # Set the internal format class
        self._format = format_class

        # Set the template
        self._template = template

        # Get the template format
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        self._template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)

        # Get the first format instance and get the models
        first_filename = self.get_filename(image_indices[0])
        format_instance = self._format(first_filename)
        self._understanding = self._format.understand(first_filename)
        self._detector = format_instance.get_detector()
        self._beam = format_instance.get_beam()

        # Create the format cache and cache the first format
        self._model_cache = {}
        self._cache_models(image_indices[0], format_instance)
        self._current_format_instance = None

        # Save the image indices
        self._indices = image_indices

    def __cmp__(self, other):
        """ Compare this reader to another. """
        return (self.get_template() == other.get_template() and
                self.get_format() == other.get_format())

    def get_format(self):
        """ Get the format instance. """
        return self._format

    def get_template(self):
        """ Get the template. """
        return self._template

    def get_filename(self, image):
        """ Get the filename for the given image. """
        return self._template_format % image

    def get_image_indices(self):
        """ Get the image indices."""
        return self._indices

    def is_valid(self, image_indices=None):
        """ Check that all the images in the given range are valid.

        Params:
            image_range The range of images to check (default all)

        Returns:
            True/False the images are valid

        """
        import os

        # If no range set then use all
        if not image_indices:
            image_indices = self.get_image_indices()

        # Loop through all the images
        for image in image_indices:

            # If format is in cache, it's valid
            try:
                self._model_cache[image]
            except KeyError:

                # Read and try to cache the format, if this fails, the
                # format is invalid, so return false.
                try:
                    self._read_format_and_cache_models(image)
                except RuntimeError:
                    return False

        # All images valid
        return True

    def read(self, image_index):
        """ Read an image at the given array index. """

        # Convert from array to image numbering and ensure is valid
        if image_index not in self.get_image_indices():
            raise IndexError, 'Invalid image index'

        # Try to read format from cache, otherwise read from file
        format_instance = self._read_format_and_cache_models(image_index)

        # Set the current format instance
        self._current_format_instance = format_instance

        # Return the raw data
        return format_instance.get_raw_data()

    def get_image_size(self):
        """ Get the image size. """
        return self._detector.get_image_size()

    def get_detectorbase(self, array_index=None):
        ''' Get the instance of the detector base at this index.'''

        # Get the format instance
        if array_index == None and self._current_format_instance != None:
            format_instance = self._current_format_instance
        else:

            # Get the image index
            if array_index == None:
                image_index = 1
            else:
                image_index = array_index + 1

            # Read the format instance
            format_instance = self._read_format_and_cache_models(image_index)

        # Return the instance of detector base
        return format_instance.get_detectorbase()

    def get_detector(self):
        """ Get the detector model. """
        return self._detector

    def get_beam(self):
        """ Get the beam model. """
        return self._beam

    def _read_format_and_cache_models(self, image_index):
        """ Read and cache the format object at the given image index. """

        # Get the format instance
        format_instance = self._format(self.get_filename(image_index))

        # If the format instance is valid then cache
        if self._is_format_valid(format_instance):
            self._cache_models(image_index, format_instance)
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
        if format_instance.get_detector() != self._detector:
            return False

        # Format is valid
        return True

    def _cache_models(self, image_index, format_instance):
        """ Cache the models. """
        b = format_instance.get_beam()
        d = format_instance.get_detector()
        g = format_instance.get_goniometer()
        s = format_instance.get_scan()
        self._model_cache[image_index] = (b, d, g, s)


class BufferedImageSetReader(ImageSetReader):
    """ Class to provide cached reading of images. """

    def __init__(self, format_class, template, image_indices, max_cache=None):
        """ Initialise the BufferedImageSetReader class.

        Params:
            reader The reader class that does the actual reading
            max_cache The maximum images to cache

        """
        from collections import OrderedDict

        # Initialise the base class
        ImageSetReader.__init__(self, format_class, template, image_indices)

        # Set the maximum images to cache
        if max_cache:
            self._max_cache = max_cache
        else:
            self._max_cache = 1

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
        image = ImageSetReader.read(self, index)

        # If the size of the cache is > max cache then remove an image
        # from the cache in a FIFO manner.
        if len(self._image_cache) >= self._max_cache:
            self._image_cache.pop(self._image_cache.keys()[0])

        # Add the image to the cache
        self._image_cache[index] = image

        # Return the image image
        return image


class ImageSet(object):
    """ A class exposing the external image set interface. """

    def __init__(self, reader, indices=None):
        """ Initialise the ImageSet object.

        Params:
            reader The reader object
            array_range The image range (first, last)

        """
        # If no reader is set then throw an exception
        if not reader:
            raise ValueError("ImageSet needs a reader!")

        # Set the reader
        self._reader = reader

        # Set the array range or get the range from the reader
        if indices:
            if isinstance(indices, list):
                self._indices = indices
            else:
                raise ValueError("indices should be a list")
        else:
            self._indices = self._reader.get_image_indices()

    def __getitem__(self, item):
        """ Get an item from the image set stream.

        If the item is an index, read and return the image at the given index.
        Otherwise, if the item is a slice, then create a new ImageSet object
        with the given number of array indices from the slice.

        Params:
            item The index or slice

        Returns:
            An image or new ImageSet object

        """
        if isinstance(item, slice):
            indices = self.indices()[item]
            if len(indices) == 0:
                raise IndexError("Bad slice")
            return ImageSet(self._reader, indices)
        else:
            return self._reader.read(self.indices()[item])

    def __len__(self):
        """ Return the number of images in this image set. """
        return len(self.indices())

    def __str__(self):
        """ Return the array indices of the image set as a string. """
        return str(self.indices())

    def __iter__(self):
        """ Iterate over the array indices and read each image in turn. """
        for f in self.indices():
            yield self._reader.read(f)

    def __cmp__(self, other):
        """ Compare this image set to another. """
        return self.reader() == other.reader()

    def indices(self):
        """ Return the indices in the image set. """
        return self._indices

    def is_valid(self):
        """ Validate all the images in the image set. Can take a long time. """
        return self._reader.is_valid(self._indices)

    def get_detector(self):
        """ Get the detector. """
        return self._reader.get_detector()

    def get_beam(self):
        """ Get the beam. """
        return self._reader.get_beam()

    def reader(self):
        """ Return the image set reader. """
        return self._reader

    def get_image_size(self):
        """ Get the image size. """
        return self._reader.get_image_size()

    def get_detectorbase(self, index=None):
        """ Get the detector base instance for the given index. """
        return self._reader.get_detectorbase(index)


class ImageSetFactory:
    """ The factory class to create a image set object. """

    @staticmethod
    def image_set(argument, check_headers = False):
        '''Get an image set.'''
        import os
        from dxtbx.format.Registry import Registry
        from dxtbx.sweep_filenames import template_regex

        # If a list if given just use given items, otherwise if a string
        # has been given then find matching filenames
        if isinstance(argument, list):
            filenames = argument

        elif isinstance(argument, str):
            from dxtbx.sweep_filenames import find_matching_images
            filenames = find_matching_images(argument)

        else:
            raise RuntimeError, 'unknown argument passed to image set factory'

        # Ensure we have enough images and format has been specified
        assert(len(filenames) > 0)

        # Get the format object
        format_class = Registry.find(filenames[0])

        # Get the first image and our understanding
        first_image = filenames[0]

        # Get the directory and first filename and set the template format
        directory, first_image_name = os.path.split(first_image)
        template, first_image_index = template_regex(first_image_name)

        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
        template = os.path.join(directory, template)

        # Check the input filenames are valid
        image_indices = [first_image_index]
        assert(first_image_index > 0)
        for j, image in enumerate(filenames[1:]):
            path, image_name = os.path.split(image)
            assert(path == directory)
            image_index = int(image_name.replace(pfx, '').replace(sfx, ''))
            assert(image_index > 0)
            image_indices.append(image_index)

        # Create the image set object
        image_set = ImageSet(BufferedImageSetReader(format_class,
            template, image_indices))

        # Check the image set is valid
        if check_headers and not image_set.is_valid():
            raise RuntimeError('Invalid image set of images')

        # Return the image set
        return image_set
