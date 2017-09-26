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
import boost.python
import dxtbx.format.image
ext = boost.python.import_ext("dxtbx_ext")
from dxtbx_imageset_ext import *


class MemReader(object):
  '''A reader for data already loaded in memory'''

  def __init__(self, images):
    self._images = images

  def paths(self):
    return ["" for im in self._images]

  def identifiers(self):
    return self.paths()

  def __len__(self):
    return len(self._images)

  def read(self, index):
    format_instance = self._images[index]
    return format_instance.get_raw_data()

  def is_single_file_reader(self):
    return False

  def master_path(self):
    return ''

class MemMasker(object):

  def __init__(self, images):
    self._images = images

  def get(self, index, goniometer=None):
    format_instance = self._images[index]
    return format_instance.get_mask(goniometer=goniometer)

  def paths(self):
    return ["" for im in self._images]

  def identifiers(self):
    return self.paths()

  def __len__(self):
    return len(self._images)


class ImageSetAux(boost.python.injector, ImageSet):
  '''
  A class to inject additional methods into the imageset class

  '''

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
      if item.start is None:
        start = 0
      else:
        start = item.start
      if item.stop is None:
        stop = len(self)
      else:
        stop = item.stop
      if item.step is not None and item.step != 1:
        raise IndexError("Step must be 1")
      return self.partial_set(start, stop)
    else:
      return self.get_corrected_data(item)

  def __iter__(self):
    ''' Iterate over the array indices and read each image in turn. '''
    for i in range(len(self)):
      yield self[i]

  def get_vendortype(self, index):
    ''' Get the vendor information. '''
    return self.data().get_vendor()

  def get_format_class(self):
    ''' Get format class name '''
    return self.data().get_format_class()

  def params(self):
    ''' Get the parameters '''
    return self.data().get_params()

  def get_detectorbase(self, index):
    '''
    A function to be injected into the imageset to get the detectorbase instance

    '''
    kwargs = self.params()
    if self.data().has_single_file_reader():
      format_instance = self.get_format_class().get_instance(
        self.data().get_master_path(), **kwargs)
      return format_instance.get_detectorbase(self.indices()[index])
    else:
      format_instance = self.get_format_class().get_instance(
        self.get_path(index), **kwargs)
      return format_instance.get_detectorbase()

  def reader(self):
    '''
    Return the reader

    '''
    return self.data().reader()

  def masker(self):
    '''
    Return the masker

    '''
    return self.data().masker()

  def paths(self):
    '''
    Return the list of paths

    '''
    return [self.get_path(i) for i in range(len(self))]

#class ImageSet(object):

  # def __init__(self,
  #              data,
  #              indices = None):

  #   # Set the imageset data
  #   self._data = data

  #   # Check if the indices have been set
  #   if indices:
  #     assert min(indices) >= 0
  #     assert max(indices) < len(data.reader)
  #     self._indices = indices
  #   else:
  #     self._indices = list(range(len(data.reader)))

  #   # Image cache
  #   self.image_cache = None

  #   # Some static stuff
  #   self.external_lookup = ExternalLookup()

  # def __getitem__(self, item):
  #   ''' Get an item from the image set stream.

  #   If the item is an index, read and return the image at the given index.
  #   Otherwise, if the item is a slice, then create a new ImageSet object
  #   with the given number of array indices from the slice.

  #   Params:
  #       item The index or slice

  #   Returns:
  #       An image or new ImageSet object

  #   '''
  #   if isinstance(item, slice):
  #     subset = ImageSet(self._data, self._indices[item])
  #     subset.external_lookup = self.external_lookup
  #     return subset
  #   else:
  #     return self.get_corrected_data(item)

  # def get_raw_data(self, index):
  #   '''
  #   Get the image at the given index

  #   '''
  #   if self.image_cache is not None and self.image_cache[0] == index:
  #     image = self.image_cache[1]
  #   else:
  #     image = self._data.data(self._indices[index])
  #     if not isinstance(image, tuple):
  #       image = (image,)
  #     self.image_cache = (index, image)
  #   return image

  # def get_corrected_data(self, index):
  #   '''
  #   Get the corrected data: (raw_data - pedestal) * gain

  #   '''
  #   data = self.get_raw_data(index)
  #   gain = self.get_gain(index)
  #   pedestal = self.get_pedestal(index)
  #   if gain is None:
  #     gain = [None] * len(data)
  #   if pedestal is None:
  #     pedestal = [None] * len(data)
  #   result = []
  #   for d, p, g in zip(data, pedestal, gain):
  #     r = d.as_double()
  #     if p is not None:
  #       r = r - p
  #     if g is not None:
  #       r = r / g
  #     result.append(r)
  #   return tuple(result)

  # def get_gain(self, index):
  #   '''
  #   Get the gain map

  #   '''
  #   from scitbx.array_family import flex
  #   gain = [p.get_gain() for p in self.get_detector(index)]
  #   if all([g > 0 for g in gain]):
  #     return gain
  #   return self.external_lookup.gain.data

  # def get_pedestal(self, index):
  #   '''
  #   Get the pedestal

  #   '''
  #   from scitbx.array_family import flex
  #   return self.external_lookup.pedestal.data

  # def get_mask(self, index, goniometer=None):
  #   '''
  #   Get the mask at the given index.
  #   Queries a format object for a dynamic mask if it exists.
  #   Otherwise uses image and trusted range.

  #   '''

  #   # Compute the trusted range mask
  #   image = self.get_raw_data(index)
  #   detector = self.get_detector(index)
  #   assert(len(image) == len(detector))
  #   mask = []
  #   for im, panel in zip(image, detector):
  #     mask.append(panel.get_trusted_range_mask(im))
  #   mask = tuple(mask)

  #   # Check for a dynamic mask
  #   if goniometer is None:
  #     goniometer = self.get_goniometer(index)
  #   dyn_mask = self._data.mask(self._indices[index], goniometer=goniometer)
  #   if dyn_mask is not None:
  #     mask = tuple([m1 & m2 for m1, m2 in zip(dyn_mask, mask)])

  #   # Get the external mask
  #   ext_mask = self.external_lookup.mask.data

  #   # Return a combination mask
  #   if ext_mask is not None:
  #     mask = tuple([m1 & m2 for m1, m2 in zip(mask, ext_mask)])
  #   return mask

  # def data(self):
  #   return self._data

  # def indices(self):
  #   ''' Return the indices '''
  #   return self._indices

  # def __len__(self):
  #   ''' Return the number of images in this image set. '''
  #   return len(self._indices)

  # def __iter__(self):
  #   ''' Iterate over the array indices and read each image in turn. '''
  #   for i in range(len(self)):
  #     yield self.get_raw_data(i)

  # def __eq__(self, other):
  #   ''' Compare this image set to another. '''
  #   if other is None:
  #     return False
  #   if other is self:
  #     return True
  #   return self.paths() == other.paths()

  # def get_detector(self, index=0):
  #   ''' Get the detector. '''
  #   return self._data.detector[self._indices[index]]

  # def set_detector(self, detector, index=0):
  #   ''' Set the detector model.'''
  #   self._data.detector[self._indices[index]] = detector

  # def get_beam(self, index=0):
  #   ''' Get the beam. '''
  #   return self._data.beam[self._indices[index]]

  # def set_beam(self, beam, index=0):
  #   ''' Set the beam model.'''
  #   self._data.beam[self._indices[index]] = beam

  # def get_goniometer(self, index=0):
  #   ''' Get the goniometer model. '''
  #   if self._data.goniometer is None or len(self._data.goniometer) == 0:
  #     return None
  #   return self._data.goniometer[self._indices[index]]

  # def set_goniometer(self, goniometer, index=0):
  #   ''' Set the goniometer model. '''
  #   self._data.goniometer[self._indices[index]] = goniometer

  # def get_scan(self, index=0):
  #   ''' Get the scan model. '''
  #   if self._data.goniometer is None or len(self._data.goniometer) == 0:
  #     return None
  #   return self._data.scan[self._indices[index]]

  # def set_scan(self, scan, index=0):
  #   ''' Set the scan model. '''
  #   if scan is not None:
  #     assert len(scan) == 1
  #   self._data.scan[self._indices[index]] = scan

  # def get_vendortype(self, index):
  #   ''' Get the vendor information. '''
  #   return self._data.properties['vendor']

#   def paths(self):
#     ''' Return a list of filenames referenced by this set. '''
#     paths = self._data.paths()
#     if self._data.reader.is_single_file_reader():
#       assert len(paths) == 1
#       return [paths[0] for i in self._indices]
#     return [paths[i] for i in self._indices]

  # def get_path(self, index):
  #   return self.paths()[index]

  # def get_image_identifier(self, index):
  #   ''' Get the path for the index '''
  #   return self._data.identifiers(self._indices[index])

#   def get_format_class(self):
#     ''' Get format class name '''
#     return self._data.properties['format']

#   def reader(self):
#     return self._data.reader

#   def masker(self):
#     return self._data.masker

#   def params(self):
#     return self._data.properties.get('params', {})

  # def complete_set(self):
  #   '''
  #   Return an imageset with all images

  #   '''
  #   return ImageSet(self._data)

  # def as_imageset(self):
  #   '''
  #   Return itself

  #   '''
  #   return self


class ImageSweepAux(boost.python.injector, ImageSweep):

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
      if item.start is None:
        start = 0
      else:
        start = item.start
      if item.stop is None:
        stop = len(self)
      else:
        stop = item.stop
      if item.step != None:
        raise IndexError('Sweeps must be sequential')
      return self.partial_set(start, stop)
    else:
      return self.get_corrected_data(item)

  def get_template(self):
    ''' Return the template '''
    return self.data().get_template()



# class ImageSweep(ImageSet):
#   ''' A class exposing the external sweep interface. '''

#   def __init__(self,
#                data,
#                indices=None,
#                beam=None,
#                goniometer=None,
#                detector=None,
#                scan=None):
#     ''' Create the sweep.

#     If the models are given here. They are used, otherwise the models
#     are read from the files themselves, with the beam, detector and
#     goniometer taken from the first image and the scan read from the
#     whole range of images. The scan must be consistent with the indices
#     given if both are specified.

#     Params:
#         reader The reader class
#         indices The list of image indices
#         beam The beam model
#         goniometer The goniometer model
#         detector The detector model
#         scan The scan model

#     '''

#     if indices:
#       assert min(indices) >= 0
#       assert max(indices) < len(data.reader)
#       assert all(i1+1 == i2 for i1, i2 in zip(indices[:-1], indices[1:]))
#       assert scan is None or indices is None or len(indices) == len(scan)

#     ImageSet.__init__(self, data, indices)
#     self._beam = beam
#     self._goniometer = goniometer
#     self._detector = detector
#     self._scan = scan
#     for i in range(len(self)):
#       ImageSet.set_beam(self, self._beam, i)
#       ImageSet.set_detector(self, self._detector, i)
#       ImageSet.set_goniometer(self, self._goniometer, i)
#       if self._scan is not None:
#         s = self._scan[i]
#       else:
#         s = None
#       ImageSet.set_scan(self, s, i)

#   def __getitem__(self, item):
#     ''' Get an item from the sweep stream.

#     If the item is an index, read and return the image at the given index.
#     Otherwise, if the item is a slice, then create a new Sweep object
#     with the given number of array indices from the slice.

#     Params:
#         item The index or slice

#     Returns:
#         An image or new Sweep object

#     '''
#     if isinstance(item, slice):
#       if item.step != None:
#         raise IndexError('Sweeps must be sequential')

#       if self._scan is None:
#         scan = None
#       else:
#         scan = self._scan[item]

#       # Create new imageset
#       subset = ImageSweep(
#         self._data,
#         self._indices[item],
#         beam = self._beam,
#         detector = self._detector,
#         goniometer = self._goniometer,
#         scan = scan)

#       # Set external lookup maps
#       subset.external_lookup = self.external_lookup

#       return subset
#     else:
#       return self.get_corrected_data(item)

#   def get_array_range(self):
#     ''' Get the array range. '''
#     return self.get_scan().get_array_range()

#   def get_beam(self, index=None):
#     ''' Get the beam. '''
#     return self._beam

#   def get_detector(self, index=None):
#     ''' Get the detector. '''
#     return self._detector

#   def get_goniometer(self, index=None):
#     ''' Get goniometer, '''
#     return self._goniometer

#   def get_scan(self, index=None):
#     ''' Get the scan.'''
#     if index is not None:
#       return self._scan[index]
#     return self._scan

#   def set_beam(self, beam):
#     ''' Set the beam. '''
#     self._beam = beam
#     for i in range(len(self)):
#       ImageSet.set_beam(self, self._beam, i)

#   def set_goniometer(self, goniometer):
#     ''' Set the goniometer model '''
#     self._goniometer = goniometer
#     for i in range(len(self)):
#       ImageSet.set_goniometer(self, self._goniometer, i)

#   def set_detector(self, detector):
#     ''' Set the detector model. '''
#     self._detector = detector
#     for i in range(len(self)):
#       ImageSet.set_detector(self, self._detector, i)

  # def set_scan(self, scan):
  #   ''' Set the scan model. '''
  #   if len(scan) != len(self):
  #     i0, i1 = scan.get_array_range()
  #     j0, j1 = self.get_scan().get_array_range()
  #     assert i0 >= j0
  #     assert i1 > i0
  #     k0 = i0 - j0
  #     k1 = i1 - j0
  #     index0 = self._indices[k0]
  #     index1 = index0 + (i1 - i0)
  #     self._indices = list(range(index0, index1))
  #   self._scan = scan
  #   for i in range(len(self)):
  #     ImageSet.set_scan(self, self._scan[i], i)

#   def get_template(self):
#     ''' Return the template '''
#     return self._data.properties['template']

#   def as_imageset(self):
#     '''
#     Return as an imageset

#     '''
#     return ImageSet(self._data, self._indices)


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
                    check_format=True,
                    beam=None,
                    detector=None,
                    goniometer=None,
                    scan=None):
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
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
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
    # try:
    #   filenames = [template_format % (i+1) for i in range(*array_range)]
    # except Exception:
    #   if check_format:
    #     raise
    #   else:
    #     filenames = []
    sweep = format_class.get_imageset(
      filenames,
      beam          = beam,
      detector      = detector,
      goniometer    = goniometer,
      scan          = scan,
      format_kwargs = format_kwargs,
      template      = template,
      as_sweep      = True,
      check_format  = check_format,
      single_file_indices = range(*array_range))

    # Return the sweep
    return sweep
