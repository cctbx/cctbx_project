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
from __future__ import absolute_import, division


# class ReaderBase(object):
#   '''The imageset reader base class.'''

#   def __init__(self):
#     pass

#   def __eq__(self, other):
#     pass

#   def get_image_paths(self, indices=None):
#     pass

#   def get_image_size(self, panel=0):
#     pass

#   def get_format(self, index=None):
#     pass

#   def get_format_class(self, index=None):
#     pass

#   def get_path(self, index=None):
#     pass

#   def is_valid(self, indices=None):
#     pass

#   def read(self, index=None):
#     pass

#   def get_detectorbase(self, index=None):
#     pass

#   def get_vendortype(self, index=None):
#     pass

#   def get_detector(self, index=None):
#     pass

#   def get_goniometer(self, index=None):
#     pass

#   def get_beam(self, index=None):
#     pass

#   def get_scan(self, index=None):
#     pass


# class NullReader(ReaderBase):
#   ''' A placeholder reader. '''

#   def __init__(self, filenames, single_file=False):
#     ReaderBase.__init__(self)
#     self._filenames = filenames
#     self._single_file = single_file

#   def is_single_file_reader(self):
#     ''' Return if single file reader '''
#     return self._single_file

#   def __eq__(self, other):
#     ''' Compare with another reader. '''
#     return isinstance(other, NullReader)

#   def get_image_paths(self, indices=None):
#     ''' Get the image paths. '''
#     if indices == None:
#       return list(self._filenames)
#     return self._filenames(indices)

#   def get_format(self, index=None):
#     ''' Get the format. '''
#     return None

#   def get_format_class(self, index=None):
#     ''' Get the format class. '''
#     return None

#   def get_path(self, index=None):
#     ''' Get an image path. '''
#     if index == None or self._single_file:
#       return self._filenames[0]
#     return self._filenames[index]

#   def is_valid(self, indices=None):
#     ''' Return whether the reader is valid. '''
#     return True

#   def read(self, index):
#     raise RuntimeError('NullReader doesn\'t have image data')

#   def get_detector(self, index=None):
#     '''Get the detector instance.'''
#     raise RuntimeError('NullReader doesn\'t have detector data')

#   def get_beam(self, index=None):
#     '''Get the beam instance.'''
#     raise RuntimeError('NullReader doesn\'t have beam data')

#   def get_goniometer(self, index=None):
#     '''Get the goniometer instance.'''
#     raise RuntimeError('NullReader doesn\'t have goniometer data')

#   def get_scan(self, index=None):
#     '''Get the scan instance.'''
#     raise RuntimeError('NullReader doesn\'t have scan data')


# class SingleFileReader(ReaderBase):
#   '''The single file reader class.'''

#   def __init__(self, format_instance = None):
#     '''Initialise the reader class.'''
#     ReaderBase.__init__(self)

#     # Set the format instance
#     self._format = format_instance
#     self._recover = self.__getstate__()

#   def __eq__(self, other):
#     '''Compare the reader to another reader.'''
#     return self.get_format() == other.get_format()

#   def __getstate__(self):
#     return self.get_format().__class__, self.get_format().get_image_file()

#   def __setstate__(self, state):
#     self._format = state[0](state[1])

#   def nullify_format_instance(self):
#     self._recover = self.__getstate__()
#     self._format = None

#   def get_image_paths(self, indices=None):
#     '''Get the image paths within the file.'''

#     # Get paths for each file
#     filenames = [self.get_format().get_image_file(i)
#                  for i in range(self.get_format().get_num_images())]

#     # Return within the given range
#     if indices == None:
#       return filenames
#     else:
#       return [filenames[i] for i in indices]

#   def get_format(self, index=None):
#     '''Get the format instance'''
#     if self._format is None and self._recover is not None:
#       self.__setstate__(self._recover)
#     return self._format

#   def get_format_class(self, index=None):
#     '''Get the format class'''
#     return self.get_format().__class__

#   def get_path(self, index=None):
#     '''Get the image file for the given index.'''
#     return self.get_format().get_image_file(index)

#   def is_valid(self, indices=None):
#     '''Ensure the reader is valid.'''
#     return True

#   def read(self, index=None):
#     '''Get the image data.'''
#     return self.get_format().get_raw_data(index)

#   def get_detectorbase(self, index=None):
#     '''Get the detector base instance.'''
#     return self.get_format().get_detectorbase(index)

#   def get_vendortype(self, index=None):
#     return self.get_format().get_vendortype()

#   def get_detector(self, index=None):
#     '''Get the detector instance.'''
#     return self.get_format().get_detector(index)

#   def get_beam(self, index=None):
#     '''Get the beam instance.'''
#     return self.get_format().get_beam(index)

#   def get_goniometer(self, index=None):
#     '''Get the goniometer instance.'''
#     return self.get_format().get_goniometer(index)

#   def get_scan(self, index=None):
#     '''Get the scan instance.'''
#     return self.get_format().get_scan(index)

#   def is_single_file_reader(self):
#     ''' Return if single file reader '''
#     return True

#   def __deepcopy__(self, memo):
#     '''
#     Override deep copy behaviour to use same format instance

#     '''
#     return SingleFileReader(self.get_format())


# class MultiFileState(object):
#   '''A class to keep track of multi file reader state.'''

#   def __init__(self, format_class, format_kwargs=None):
#     '''Initialise with format class.'''
#     self._format_class = format_class
#     self._current_format_instance = None

#     if format_kwargs is None:
#       self._format_kwargs = {}
#     else:
#       self._format_kwargs = format_kwargs

#   def format_class(self):
#     '''Get the format class.'''
#     return self._format_class

#   def load_file(self, filename):
#     '''Load the file with the given filename.'''

#     # Check the current format is the one we need
#     if (self.get_format() == None or
#         filename != self.get_format().get_image_file()):

#       # Read the format instance
#       format_instance = self._format_class(filename, **self._format_kwargs)

#       # Check the format instance is valid
#       if not self._is_format_valid(format_instance):
#         RuntimeError("Format is invalid.")

#       # Set the current format instance
#       self._current_format_instance = format_instance

#   def get_format(self):
#     '''Get the current format instance.'''
#     return self._current_format_instance

#   def _is_format_valid(self, format_instance):
#     '''Check if the format object is valid.'''
#     return format_instance.understand(format_instance.get_image_file())

#   def __getstate__(self):
#     ''' Save the current image and format class for pickling. '''
#     if self._current_format_instance is not None:
#       current_filename = self._current_format_instance.get_image_file()
#     else:
#       current_filename = None
#     return { 'format_class' : self._format_class,
#              'current_filename' : current_filename,
#              'format_kwargs' : self._format_kwargs }

#   def __setstate__(self, state):
#     ''' Set the format class and load the image. '''
#     self._format_class = state['format_class']
#     self._format_kwargs = state['format_kwargs']
#     self._current_format_instance = None
#     if state['current_filename'] is not None:
#       self.load_file(state['current_filename'])


# class NullFormatChecker(object):
#   def __call__(self, fmt):
#     return True


# class MultiFileReader(ReaderBase):
#   '''A multi file reader class implementing the ReaderBase interface.'''

#   def __init__(self, format_class, filenames, formatchecker=None,
#                format_kwargs=None):
#     '''Initialise the reader with the format and list of filenames.'''
#     ReaderBase.__init__(self)

#     import os

#     # Ensure we have enough images and format has been specified
#     assert(format_class != None)
#     assert(len(filenames) > 0)

#     # Save the image indices
#     self._filenames = filenames

#     # Handle the state of the MultiFileReader class
#     self._state = MultiFileState(format_class, format_kwargs=format_kwargs)

#     # A function object to check formats are valid
#     if formatchecker != None:
#       self._is_format_valid = formatchecker
#     else:
#       self._is_format_valid = NullFormatChecker()

#   def __eq__(self, other):
#     '''Compare the reader by format class and filename list.'''
#     return (self.get_format_class() == other.get_format_class() and
#             self.get_image_paths() == other.get_image_paths())

#   def get_image_paths(self, indices=None):
#     '''Get the list of image paths.'''
#     if indices == None:
#       return list(self._filenames)
#     else:
#       return [self._filenames[i] for i in indices]

#   def get_format_class(self):
#     '''Get the format class.'''
#     return self._state.format_class()

#   def get_image_size(self, panel=0):
#     '''Get the image size.'''
#     return self.get_format().get_detector()[panel].get_image_size()

#   def get_path(self, index=None):
#     '''Get the path the given index.'''
#     if index == None:
#       return self.get_format().get_image_file()
#     else:
#       return self._filenames[index]

#   def get_detectorbase(self, index=None):
#     '''Get the detector base instance at given index.'''
#     return self.get_format(index).get_detectorbase()

#   def get_vendortype(self, index=None):
#     return self.get_format(index).get_vendortype()

#   def get_detector(self, index=None):
#     '''Get the detector instance at given index.'''
#     return self.get_format(index).get_detector()

#   def get_beam(self, index=None):
#     '''Get the beam instance at given index.'''
#     return self.get_format(index).get_beam()

#   def get_goniometer(self, index=None):
#     '''Get the goniometer instance at given index.'''
#     return self.get_format(index).get_goniometer()

#   def get_scan(self, index=None):
#     '''Get the scan instance at given index.'''
#     return self.get_format(index).get_scan()

#   def read(self, index=None):
#     '''Read the image frame at the given index.'''

#     # Get the format instance
#     format_instance = self.get_format(index)

#     return format_instance.get_raw_data()

#   def get_format(self, index=None):
#     '''Get the format at the given index.'''
#     return self._update_state(index).get_format()

#   def _update_state(self, index=None):
#     '''Update the state and load file at given index.'''
#     if index is not None:
#       self._state.load_file(self.get_path(index))
#     elif self._state.get_format() == None:
#       self._state.load_file(self.get_path(0))

#     return self._state

#   def is_valid(self, indices=None):
#     '''Ensure imageset is valid.'''
#     import os

#     # If no indices, get indices of all filenames
#     if indices == None:
#       indices = range(len(self._filenames))

#     # Loop through all the images
#     for index in indices:

#       # Read and try to cache the format, if this fails, the
#       # format is invalid, so return false.
#       try:
#         format_instance = self.get_format(index)
#       except IndexError, RuntimeError:
#         return False

#       # Check the format experimental models
#       if not self._is_format_valid(format_instance):
#         return False

#     # All images valid
#     return True

#   def is_single_file_reader(self):
#     ''' Return if single file reader '''
#     return False

# class MemReader(ReaderBase):
#   '''A reader for data already loaded in memory'''

#   def __init__(self, images):
#     self._images = images

#   def __eq__(self, other):
#     return self.get_format_class() == other.get_format_class()

#   def get_image_paths(self, indices=None):
#     raise NotImplementedError("MemReader has no image paths")

#   def get_image_size(self, panel=0):
#     return self.get_format().get_detector()[panel].get_image_size()

#   def get_format(self, index=None):
#     return self._images[index]

#   def get_format_class(self, index=None):
#     return self._images[index].__class__

#   def get_path(self, index=None):
#     raise NotImplementedError("MemReader has no image paths")

#   def is_valid(self, indices=None):
#     return True

#   def read(self, index=None):
#     # Get the format instance
#     format_instance = self.get_format(index)

#     return format_instance.get_raw_data()

#   def get_detectorbase(self, index=None):
#     return self._images[index].get_detectorbase()

#   def get_vendortype(self, index=None):
#     return self._images[index].get_vendortype()

#   def get_detector(self, index=None):
#     return self._images[index].get_detector()

#   def get_goniometer(self, index=None):
#     return self._images[index].get_goniometer()

#   def get_beam(self, index=None):
#     return self._images[index].get_beam()

#   def get_scan(self, index=None):
#     return self._images[index].get_scan()

class ExternalLookupItem(object):
  '''
  Save an external lookup item.

  '''
  def __init__(self, data=None, filename=None):
    self.data = data
    self.filename = filename

class ExternalLookup(object):
  '''
  A class to hold some external lookup data

  '''
  def __init__(self):
    self.mask = ExternalLookupItem()
    self.gain = ExternalLookupItem()
    self.pedestal = ExternalLookupItem()


class ImageSet(object):
  ''' A class exposing the external image set interface. '''

  def __init__(self, 
               reader,
               masker=None,
               properties={},
               detectorbase_factory=None):
    ''' Initialise the ImageSet object.

    Params:
        reader The reader object

    '''


    # If no reader is set then throw an exception
    if not reader:
      raise ValueError("ImageSet needs a reader!")

    # Check input
    if masker:
      assert len(masker) == len(reader)
    if detectorbase_factory:
      assert len(detectorbase_factory) == len(reader)

    # Set the reader
    self._reader = reader
    self._masker = masker
    self._beam_list = dict()
    self._detector_list = dict()
    self._goniometer_list = dict()
    self._scan_list = dict()
    self._properties = properties
    self._detectorbase_factory = detectorbase_factory

    # Image cache
    self.image_cache = None

    # Some static stuff
    self.external_lookup = ExternalLookup()

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

      # Get the filenames
      filenames = self.paths()[item]
      
      # Get reader
      reader = self._reader.copy(filenames) 

      # Get masker
      if self._masker:
        masker = self._masker.copy(filenames)
      else:
        masker = None

      # Get detector base factory
      if self._detectorbase_factory:
        dbfact = self._detectorbase_factory.copy(filenames)
      else:
        dbfact = None

      # Create new imageset
      subset = ImageSet(
        reader,
        masker = masker,
        properties = self._properties,
        detectorbase_factory = dbfact)

      # Set the models
      for i in range(len(self))[item]:
        subset.set_beam(index, self.get_beam(index))
        subset.set_detector(index, self.get_detector(index))
        subset.set_goniometer(index, self.get_goniometer(index))
        subset.set_scan(index, self.get_scan(index))
      
      # Set external lookup maps
      subset.external_lookup = self.external_lookup
      
      return subset
    else:
      return self.get_corrected_data(item)

  def get_raw_data(self, index):
    '''
    Get the image at the given index

    '''
    if self.image_cache is not None and self.image_cache[0] == index:
      image = self.image_cache[1]
    else:
      image = self.reader().read(index)
      if not isinstance(image, tuple):
        image = (image,)
      self.image_cache = (index, image)
    return image

  def get_corrected_data(self, index):
    '''
    Get the corrected data: (raw_data - pedestal) * gain

    '''
    data = self.get_raw_data(index)
    gain = self.get_gain(index)
    pedestal = self.get_pedestal(index)
    if gain is None:
      gain = [None] * len(data)
    if pedestal is None:
      pedestal = [None] * len(data)
    result = []
    for d, p, g in zip(data, pedestal, gain):
      r = d.as_double()
      if p is not None:
        r = r - p
      if g is not None:
        r = r / g
      result.append(r)
    return tuple(result)

  def get_gain(self, index):
    '''
    Get the gain map

    '''
    from scitbx.array_family import flex
    gain = [p.get_gain() for p in self.get_detector(0)]
    if all([g > 0 for g in gain]):
      return gain
    return self.external_lookup.gain.data

  def get_pedestal(self, index):
    '''
    Get the pedestal

    '''
    from scitbx.array_family import flex
    return self.external_lookup.pedestal.data

  def get_mask(self, index, goniometer=None):
    '''
    Get the mask at the given index.
    Queries a format object for a dynamic mask if it exists.
    Otherwise uses image and trusted range.

    '''

    # Compute the trusted range mask
    image = self.get_raw_data(index)
    detector = self.get_detector(index)
    assert(len(image) == len(detector))
    mask = []
    for im, panel in zip(image, detector):
      mask.append(panel.get_trusted_range_mask(im))
    mask = tuple(mask)

    # Check for a dynamic mask
    if goniometer is None:
      goniometer = self.get_goniometer(index)
    dyn_mask = self._masker.get(index, goniometer=gonioneter)
    if dyn_mask is not None:
      mask = tuple([m1 & m2 for m1, m2 in zip(dyn_mask, mask)])

    # Get the external mask
    ext_mask = self.external_lookup.mask.data

    # Return a combination mask
    if ext_mask is not None:
      mask = tuple([m1 & m2 for m1, m2 in zip(mask, ext_mask)])
    return mask

  def __len__(self):
    ''' Return the number of images in this image set. '''
    return len(self._reader)

  def __str__(self):
    ''' Return the array indices of the image set as a string. '''
    return str(self.paths())

  def __iter__(self):
    ''' Iterate over the array indices and read each image in turn. '''
    for i in len(self):
      yield self._reader.read(i)

  def __eq__(self, other):
    ''' Compare this image set to another. '''
    if other is None:
      return False
    if other is self:
      return True
    return self.paths() == other.paths()

  def paths(self):
    ''' Return a list of filenames referenced by this set. '''
    return self._reader.paths()

  def get_detector(self, index=None):
    ''' Get the detector. '''
    return self._detector_list[index]

  def set_detector(self, detector, index=None):
    ''' Set the detector model.'''
    self._detector_list[index] = detector

  def get_beam(self, index=None):
    ''' Get the beam. '''
    return self._beam_list[index]

  def set_beam(self, beam, index=None):
    ''' Set the beam model.'''
    self._beam_list[index] = beam

  def get_goniometer(self, index=None):
    ''' Get the goniometer model. '''
    return self._goniometer_list[index]

  def set_goniometer(self, goniometer, index=None):
    ''' Set the goniometer model. '''
    self._goniometer_list[index] = goniometer

  def get_scan(self, index=None):
    ''' Get the scan model. '''
    return self._scan[index]

  def set_scan(self, scan, index=None):
    ''' Set the scan model. '''
    self._scan_list[index] = scan

  def get_detectorbase(self, index):
    ''' Get the detector base instance for the given index. '''
    return self._detectorbase_factory(index)

  def get_vendortype(self, index):
    ''' Get the vendor information. '''
    return self._properties['vendor']

  def get_image_identifier(self, index):
    ''' Get the path for the index '''
    return self._reader.identifiers()[index]


class ImageGrid(ImageSet):
  '''
  A class implementing an interface useful for processing grid scans

  '''
  def __init__(self, 
               reader, 
               masker=None,
               properties={},
               detectorbase_factory=None,
               grid_size=None):
    ''' Initialise the ImageSet object.

    Params:
        reader The reader object
        array_range The image range (first, last)

    '''
    super(ImageGrid, self).__init__(
      reader, 
      masker,
      properties=properties,
      detectorbase_factory=detectorbase_factory)

    # Set the grid size
    num = grid_size[0] * grid_size[1]
    assert num == len(indices)
    self._grid_size = grid_size

  def get_grid_size(self):
    '''
    Return the grid size

    '''
    return self._grid_size

  @classmethod
  def from_imageset(cls, imageset, grid_size):
    '''
    Convert an imageset into an image grid

    '''
    return cls(
      imageset._reader,
      imageset._masker,
      imageset._properties,
      imageset._detectorbase_factory,
      grid_size)


# class MemImageSet(ImageSet):
#   ''' A class exposing the external image set interface, but instead of a file list, uses
#   an already instantiated list of Format objects. Derives from ImageSet for clarity and for
#   the dials importer, but overrides all of ImageSet's methods '''

#   def __init__(self, images, indices=None, format_kwargs=None):
#     ''' Initialise the MemImageSet object.

#     Params:
#         images: list of Format objects
#         indices: list of indices into the images list

#     '''
#     # If no list of images is set then throw an exception
#     if images is None:
#       raise ValueError("MemImageSet needs a list of images!")

#     # Save a reader
#     self._reader = MemReader(images)

#     # Save the images
#     self._images = images

#     # Set the array range or get the range from the list of images
#     if indices is not None:
#       self._indices = indices
#     else:
#       self._indices = range(len(images))

#     # Image cache
#     self.image_cache = None

#     # Some static stuff
#     self.external_lookup = ExternalLookup()

#     # The format kwargs
#     if format_kwargs is None:
#       self._format_kwargs = {}
#     else:
#       self._format_kwargs = format_kwargs

#   def __getitem__(self, item):
#     ''' Get an item from the image set.

#     If the item is an index, read and return the image at the given index.
#     Otherwise, if the item is a slice, then create a new MemImageSet object
#     with the given number of array indices from the slice.

#     Params:
#         item The index or slice

#     Returns:
#         An image or new ImageSet object

#     '''
#     if isinstance(item, slice):
#       indices = self._indices[item]
#       subset = MemImageSet(
#         self._images,
#         indices,
#         format_kwargs=self.format_kwargs())
#       subset.external_lookup = self.external_lookup
#       return subset
#     else:
#       return self.get_corrected_data(item)

#   def __len__(self):
#     ''' Return the number of images in this image set. '''
#     return len(self._indices)

#   def __str__(self):
#     ''' Return the array indices of the image set as a string. '''
#     return str(self._indices)

#   def __iter__(self):
#     ''' Iterate over the array indices and read each image in turn. '''
#     for j in self._indices:
#       img = self._images[j]

#       # Yield a tuple of flex arrays
#       yield img.get_raw_data()

#   def __eq__(self, other):
#     ''' Compare this image set to another. '''
#     if other is None:
#       return False
#     if other is self:
#       return True
#     return self._images == other._images

#   def indices(self):
#     ''' Return the indices in the image set. '''
#     return list(self._indices)

#   def paths(self):
#     ''' Return a list of filenames referenced by this set. '''
#     raise NotImplementedError("No path list for an in-memory image set")

#   def is_valid(self):
#     ''' Validate all the images in the image set. Can take a long time. '''
#     return self.reader().is_valid(self._indices)

#   def get_detector(self, index=None):
#     ''' Get the detector. '''
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_detector()

#   def set_detector(self, detector):
#     ''' Set the detector model.'''
#     for img in self._images:
#       img._detector = detector

#   def get_beam(self, index=None):
#     ''' Get the beam. '''
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_beam()

#   def set_beam(self, beam):
#     ''' Set the beam model.'''
#     for img in self._images:
#       img._beam = beam

#   def get_goniometer(self, index=None):
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_goniometer()

#   def get_scan(self, index=None):
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_scan()

#   def get_image_size(self, index=0):
#     ''' Get the image size. '''
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_image_size()

#   def get_detectorbase(self, index=None):
#     ''' Get the detector base instance for the given index. '''
#     if index is None:
#       index = self._indices[0]
#     return self._images[index].get_detectorbase()

#   def reader(self):
#     ''' Return the image set reader. '''
#     return self._reader

#   def _image_index(self, index=None):
#     ''' Convert image set index to image index.'''
#     if index == None:
#       return None
#     elif index < 0 or index >= len(self._indices):
#       raise IndexError('Index out of range')
#     return self._indices[index]

#   def get_image_models(self, index=None, no_read=False):
#     ''' Get the models for the image.'''
#     image_index = self._image_index(index)
#     models = {}
#     if not no_read:
#       try:
#         models['detector'] = self.get_detector(image_index)
#       except Exception:
#         models['detector'] = None
#       try:
#         models['goniometer'] = self.get_goniometer(image_index)
#       except Exception:
#         models['goniometer'] = None
#       try:
#         models['beam'] = self.get_beam(image_index)
#       except Exception:
#         models['beam'] = None
#       try:
#         models['scan'] = self.get_scan(image_index)
#       except Exception:
#         models['scan'] = None
#     else:
#       models['detector'] = None
#       models['beam'] = None
#       models['goniometer'] = None
#       models['scan'] = None
#     return models

#   def complete_set(self):
#     ''' Return the set of all images (i.e. not just the subset). '''
#     return MemImageSet(self._images)

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
    for i in xrange(len(self)):
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

  def __init__(self, 
               reader, 
               masker = None,
               properties = {},
               detectorbase_factory=None,
               beam=None, 
               goniometer=None,
               detector=None, 
               scan=None):
    ''' Create the sweep.

    If the models are given here. They are used, otherwise the models
    are read from the files themselves, with the beam, detector and
    goniometer taken from the first image and the scan read from the
    whole range of images. The scan must be consistent with the indices
    given if both are specified.

    Params:
        reader The reader class
        indices The list of image indices
        beam The beam model
        goniometer The goniometer model
        detector The detector model
        scan The scan model

    '''
    ImageSet.__init__(self, 
                      reader, 
                      masker,
                      properties,
                      detectorbase_factory)
    self._beam = beam
    self._goniometer = goniometer
    self._detector = detector
    self._scan = scan

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

      # Get the filenames
      filenames = self.paths()[item]
      
      # Get reader
      reader = self._reader.copy(filenames) 

      # Get masker
      if self._masker:
        masker = self._masker.copy(filenames)
      else:
        masker = None

      # Get detector base factory
      if self._detectorbase_factory:
        dbfact = self._detectorbase_factory.copy(filenames)
      else:
        dbfact = None
      
      if self._scan is None:
        scan = None
      else:
        scan = self._scan[item]

      # Create new imageset
      subset = ImageSweep(
        reader,
        masker = masker,
        properties = self._properties,
        detectorbase_factory = dbfact,
        beam = self._beam,
        detector = self._detector,
        goniometer = self._goniometer,
        scan = scan)

      # Set external lookup maps
      subset.external_lookup = self.external_lookup
      
      return subset
    else:
      return self.get_corrected_data(item)

  def get_array_range(self):
    ''' Get the array range. '''
    return self.get_scan().get_array_range()

  def get_beam(self, index=None):
    ''' Get the beam. '''
    return self._beam

  def get_detector(self, index=None):
    ''' Get the detector. '''
    return self._detector

  def get_goniometer(self, index=None):
    ''' Get goniometer, '''
    return self._goniometer

  def get_scan(self, index=None):
    ''' Get the scan.'''
    return self._scan

  def set_beam(self, beam):
    ''' Set the beam. '''
    self._beam = beam

  def set_goniometer(self, goniometer):
    ''' Set the goniometer model '''
    self._goniometer = goniometer

  def set_detector(self, detector):
    ''' Set the detector model. '''
    self._detector = detector

  def set_scan(self, scan):
    ''' Set the scan model. '''
    self._scan = scan


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
    from dxtbx.sweep_filenames import group_files_by_imageset

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
  def from_template(template, image_range=None, check_headers=False,
                    check_format=True):
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

    if not check_format: assert not check_headers

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
    if check_format:
      format_class = Registry.find(filenames[0])
      from dxtbx.format.FormatMultiImage import FormatMultiImage
      if issubclass(format_class, FormatMultiImage):
        assert len(filenames) == 1
        format_instance = format_class(filenames[0])
        reader = SingleFileReader(format_instance)
      else:
        reader = MultiFileReader(format_class, filenames)
    else:
      reader = NullReader(filenames)

    # Create the sweep object
    sweep = ImageSweep(reader)

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
    from dxtbx.format.FormatMultiImage import FormatMultiImage
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


  @staticmethod
  def make_imageset(filenames, format_class=None, check_format=True,
                    single_file_indices=None, format_kwargs=None):
    '''Create an image set'''
    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage

    # Get the format object
    if format_class == None and check_format:
      format_class = Registry.find(filenames[0])
    if format_class is None:
      reader = NullReader(filenames, single_file_indices is not None)
    else:
      if issubclass(format_class, FormatMultiImage):
        assert len(set(filenames)) == 1
        if format_kwargs is None:
          format_kwargs = {}
        format_instance = format_class(filenames[0], **format_kwargs)
        reader = SingleFileReader(format_instance)
      else:
        reader = MultiFileReader(format_class, filenames,
                                 format_kwargs=format_kwargs)

    # Return the imageset
    return ImageSet(reader, indices=single_file_indices,
                    format_kwargs=format_kwargs)

  @staticmethod
  def make_sweep(template, indices, format_class=None, beam=None,
                 detector=None, goniometer=None, scan=None,
                 check_format=True, format_kwargs=None):
    '''Create a sweep'''
    import os
    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage

    indices = sorted(indices)

    # Get the template format
    count = template.count('#')
    if count > 0:
      pfx = template.split('#')[0]
      sfx = template.split('#')[-1]
      template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
      filenames = [template_format % index for index in indices]
    else:
      template_format = None
      filenames = [template]

    # Sort the filenames
    filenames = sorted(filenames)

    # Get the first image and our understanding
    first_image = filenames[0]

    # Get the directory and first filename and set the template format
    directory, first_image_name = os.path.split(first_image)
    first_image_number = indices[0]

    # Set the image range
    array_range = (min(indices) - 1, max(indices))
    if scan is not None:
      assert(array_range == scan.get_array_range())

    # Get the format object and reader
    if format_class is None and check_format:
      format_class = Registry.find(filenames[0])

    # Create the reader
    indices = None
    if format_class is None:
      if template_format is not None:
        filenames = SweepFileList(template_format, array_range)
      reader = NullReader(filenames)
    else:
      if issubclass(format_class, FormatMultiImage):
        if format_kwargs is None:
          format_kwargs = {}
        assert len(filenames) == 1
        format_instance = format_class(filenames[0], **format_kwargs)
        reader = SingleFileReader(format_instance)
        indices = list(range(*array_range))
      else:
        assert(template_format is not None)
        filenames = SweepFileList(template_format, array_range)
        reader = MultiFileReader(format_class, filenames,
                                 format_kwargs=format_kwargs)

    # Create the sweep object
    sweep = ImageSweep(
      reader,
      indices=indices,
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
      format_kwargs=format_kwargs)

    # Return the sweep
    return sweep
