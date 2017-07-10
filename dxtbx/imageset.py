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


class MemReader(object):
  '''A reader for data already loaded in memory'''

  def __init__(self, images):
    self._images = images

  def paths(self):
    raise NotImplementedError("MemReader has no image paths")

  def identifiers(self):
    raise NotImplementedError("MemReader has no image paths")

  def __len__(self):
    return len(self._images)

  def read(self, index):
    format_instance = self._images[index]
    return format_instance.get_raw_data()

  def is_single_file_reader(self):
    return False

  def master_path(self):
    return ''


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



class ImageSetData(object):

  def __init__(self,
               reader,
               masker=None,
               beam = None,
               detector = None,
               goniometer = None,
               scan = None,
               properties={}):

    # If no reader is set then throw an exception
    if reader is None:
      raise ValueError("ImageSet needs a reader!")

    # Check input
    if masker is not None:
      assert len(masker) == len(reader)

    # Set the data
    self.reader = reader
    self.masker = masker
    if beam is None:
      self.beam = dict()
    else:
      self.beam = beam
    if detector is None:
      self.detector = dict()
    else:
      self.detector = detector
    if goniometer is None:
      self.goniometer = dict()
    else:
      self.goniometer = goniometer
    if scan is None:
      self.scan = dict()
    else:
      self.scan = scan
    self.properties = properties

  def data(self, index):
    return self.reader.read(index)

  def mask(self, index, goniometer=None):
    return self.masker.get(index, goniometer=goniometer)

  def path(self, index):
    print self.paths(), index
    return self.paths()[index]

  def identifier(self, index):
    return self.reader.identifier()[index]

  def paths(self):
    return self.reader.paths()

  def identifiers(self):
    return self.reader.identifiers()


class ImageSet(object):

  def __init__(self,
               data,
               indices = None):

    # Set the imageset data
    self._data = data

    # Check if the indices have been set
    if indices:
      assert min(indices) >= 0
      assert max(indices) < len(data.reader)
      self._indices = indices
    else:
      self._indices = list(range(len(data.reader)))

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
      subset = ImageSet(self._data, self._indices[item])
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
      image = self._data.data(self._indices[index])
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
    gain = [p.get_gain() for p in self.get_detector(index)]
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
    dyn_mask = self._data.mask(self._indices[index], goniometer=goniometer)
    if dyn_mask is not None:
      mask = tuple([m1 & m2 for m1, m2 in zip(dyn_mask, mask)])

    # Get the external mask
    ext_mask = self.external_lookup.mask.data

    # Return a combination mask
    if ext_mask is not None:
      mask = tuple([m1 & m2 for m1, m2 in zip(mask, ext_mask)])
    return mask

  def data(self):
    return self._data

  def indices(self):
    ''' Return the indices '''
    return self._indices

  def __len__(self):
    ''' Return the number of images in this image set. '''
    return len(self._indices)

  def __str__(self):
    ''' Return the array indices of the image set as a string. '''
    return str(self.paths())

  def __iter__(self):
    ''' Iterate over the array indices and read each image in turn. '''
    for i in range(len(self)):
      yield self.get_raw_data(i)

  def __eq__(self, other):
    ''' Compare this image set to another. '''
    if other is None:
      return False
    if other is self:
      return True
    return self.paths() == other.paths()

  def get_detector(self, index=0):
    ''' Get the detector. '''
    return self._data.detector[self._indices[index]]

  def set_detector(self, detector, index=0):
    ''' Set the detector model.'''
    self._data.detector[self._indices[index]] = detector

  def get_beam(self, index=0):
    ''' Get the beam. '''
    return self._data.beam[self._indices[index]]

  def set_beam(self, beam, index=0):
    ''' Set the beam model.'''
    self._data.beam[self._indices[index]] = beam

  def get_goniometer(self, index=0):
    ''' Get the goniometer model. '''
    return self._data.goniometer[self._indices[index]]

  def set_goniometer(self, goniometer, index=0):
    ''' Set the goniometer model. '''
    self._data.goniometer[self._indices[index]] = goniometer

  def get_scan(self, index=0):
    ''' Get the scan model. '''
    return self._data.scan[self._indices[index]]

  def set_scan(self, scan, index=0):
    ''' Set the scan model. '''
    self._data.scan[self._indices[index]] = scan

  def get_vendortype(self, index):
    ''' Get the vendor information. '''
    return self._data.properties['vendor']

  def paths(self):
    ''' Return a list of filenames referenced by this set. '''
    paths = self._data.paths()
    if self._data.reader.is_single_file_reader():
      assert len(paths) == 1
      return [paths[0] for i in self._indices]
    return [paths[i] for i in self._indices]

  def get_path(self, index):
    return self.paths()[index]

  def get_image_identifier(self, index):
    ''' Get the path for the index '''
    return self._data.identifiers(self._indices[index])

  def get_format_class(self):
    ''' Get format class name '''
    return self._data.properties['format']

  def reader(self):
    return self._data.reader

  def masker(self):
    return self._data.masker

  def params(self):
    return self._data.properties['params']

  def complete_set(self):
    '''
    Return an imageset with all images

    '''
    return ImageSet(self._data)



def get_detectorbase(self, index):
  '''
  A function to be injected into the imageset to get the detectorbase instance

  '''
  kwargs = self.params()
  if self.reader().is_single_file_reader():
    format_instance = self.get_format_class().get_instance(self.reader().master_path(), **kwargs)
    return format_instance.get_detectorbase(self.indices()[index])
  else:
    format_instance = self.get_format_class().get_instance(self.paths()[index], **kwargs)
    return format_instance.get_detectorbase()

# Inject the function
ImageSet.get_detectorbase = get_detectorbase



class ImageGrid(ImageSet):
  '''
  A class implementing an interface useful for processing grid scans

  '''
  def __init__(self,
               reader,
               masker=None,
               properties={},
               detectorbase_reader=None,
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
      detectorbase_reader=detectorbase_reader)

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
      imageset._detectorbase_reader,
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


class ImageSweep(ImageSet):
  ''' A class exposing the external sweep interface. '''

  def __init__(self,
               data,
               indices=None,
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

    if indices:
      assert min(indices) >= 0
      assert max(indices) < len(data.reader)
      assert all(i1+1 == i2 for i1, i2 in zip(indices[:-1], indices[1:]))

    ImageSet.__init__(self, data, indices)
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

      if self._scan is None:
        scan = None
      else:
        scan = self._scan[item]

      # Create new imageset
      subset = ImageSweep(
        self._data,
        self._indices[item],
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

  def get_template(self):
    ''' Return the template '''
    return self._data.properties['template']


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



# FIXME Lots of duplication in this class, need to tidy up
class ImageSetFactory(object):
  ''' Factory to create imagesets and sweeps. '''

  @staticmethod
  def new(filenames,
          check_headers=False,
          ignore_unknown=False):
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
        if filelist[2] == True:
          iset = ImageSetFactory._create_sweep(filelist, check_headers)
        else:
          iset = ImageSetFactory._create_imageset(filelist, check_headers)
        imagesetlist.append(iset)
      except Exception, e:
        if not ignore_unknown:
          raise

    # Return the imageset list
    return imagesetlist

  @staticmethod
  def from_template(template,
                    image_range=None,
                    check_headers=False,
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
    from dxtbx.format.Format import Format

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
    filenames = [template_format % (i+1) for i in range(*array_range)]

    # Get the format class
    if check_format:
      format_class = Registry.find(filenames[0])
    else:
      format_class = Format

    # Create the sweep object
    sweep = format_class.get_imageset(
      filenames,
      template=template,
      as_sweep=True,
      check_format=check_format)

    # Return the sweep
    return [sweep]

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

    # Create the imageset
    imageset = format_class.get_imageset(filenames, as_imageset=True)

    # Return the image set
    return imageset

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
    filenames = [template_format % (i+1) for i in range(*array_range)]

    sweep = format_class.get_imageset(filenames, template=template,
                                      as_sweep=True)

    # Return the sweep
    return sweep


  @staticmethod
  def make_imageset(filenames,
                    format_class=None,
                    check_format=True,
                    single_file_indices=None,
                    format_kwargs=None):
    '''Create an image set'''
    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    from dxtbx.format.Format import Format

    # Get the format object
    if format_class == None:
      if check_format:
        format_class = Registry.find(filenames[0])
      else:
        format_class = Format
    else:
      format_class = format_class

    imageset = format_class.get_imageset(
      filenames,
      single_file_indices = single_file_indices,
      as_imageset   = True,
      format_kwargs = format_kwargs,
      check_format  = check_format)

    # Return the imageset
    return imageset

  @staticmethod
  def make_sweep(template,
                 indices,
                 format_class=None,
                 beam=None,
                 detector=None,
                 goniometer=None,
                 scan=None,
                 check_format=True,
                 format_kwargs=None):
    '''Create a sweep'''
    import os
    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    from dxtbx.format.Format import Format

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
    if format_class == None:
      if check_format:
        format_class = Registry.find(filenames[0])
      else:
        format_class = Format
    else:
      format_class = format_class

    # Done require template to be vaid if not checking format
    try:
      filenames = [template_format % (i+1) for i in range(*array_range)]
    except Exception:
      if check_format:
        raise
      else:
        filenames = []

    sweep = format_class.get_imageset(
      filenames,
      beam          = beam,
      detector      = detector,
      goniometer    = goniometer,
      scan          = scan,
      format_kwargs = format_kwargs,
      template      = template,
      as_sweep      = True,
      check_format  = check_format)

    # Return the sweep
    return sweep
