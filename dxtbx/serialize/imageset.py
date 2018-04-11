#!/usr/bin/env python
#
# dxtbx.serialize.imageset.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division

def filename_to_absolute(filename):
  ''' Convert filenames to absolute form. '''
  from os.path import abspath
  if isinstance(filename, list):
    return [abspath(f) for f in filename]

  return abspath(filename)

def filename_or_none(filename):
  if filename is None or filename == "":
    return None
  return filename_to_absolute(filename)

def basic_imageset_to_dict(imageset):
  ''' Convert an imageset to a dictionary

  Params:
      imageset The imageset

  Returns:
      A dictionary of the parameters

  '''
  from libtbx.containers import OrderedDict

  # Return the dictionary representation
  return OrderedDict([
      ("__id__", "imageset"),
      ("filenames", filename_to_absolute(imageset.paths())),
      ("mask", filename_or_none(imageset.external_lookup.mask.filename)),
      ("gain", filename_or_none(imageset.external_lookup.gain.filename)),
      ("pedestal", filename_or_none(imageset.external_lookup.pedestal.filename)),
      ("beam", imageset.get_beam(0).to_dict()),
      ("detector", imageset.get_detector(0).to_dict())])

def imagesweep_to_dict(sweep):
  ''' Convert a sweep to a dictionary

  Params:
      sweep The sweep

  Returns:
      A dictionary of the parameters

  '''
  from libtbx.containers import OrderedDict

  # Return the dictionary representation
  return OrderedDict([
      ("__id__", "imageset"),
      ("template", filename_to_absolute(sweep.get_template())),
      ("mask", filename_or_none(sweep.external_lookup.mask.filename)),
      ("gain", filename_or_none(sweep.external_lookup.gain.filename)),
      ("pedestal", filename_or_none(sweep.external_lookup.pedestal.filename)),
      ("beam", sweep.get_beam().to_dict()),
      ("detector", sweep.get_detector().to_dict()),
      ("goniometer", sweep.get_goniometer().to_dict()),
      ("scan", sweep.get_scan().to_dict())])

def imageset_to_dict(imageset):
  ''' Convert the imageset to a dictionary

  Params:
      imageset The imageset

  Returns:
      A dictionary of the parameters

  '''
  from dxtbx.imageset import ImageSet, ImageSweep
  from dxtbx.format.image import ImageBool, ImageDouble

  # If this is an imageset then return a list of filenames
  if isinstance(imageset, ImageSweep):
    return imagesweep_to_dict(imageset)
  elif isinstance(imageset, ImageSet):
    return basic_imageset_to_dict(imageset)
  else:
    raise TypeError("Unknown ImageSet Type")

def basic_imageset_from_dict(d, directory=None):
  ''' Construct an ImageSet class from the dictionary.'''
  from dxtbx.model import BeamFactory, DetectorFactory
  from dxtbx.imageset import ImageSetFactory
  from dxtbx.serialize.filename import load_path

  # Get the filename list and create the imageset
  filenames = map(lambda p: load_path(p, directory=directory), map(str, d['filenames']))
  imageset = ImageSetFactory.new(filenames)[0]

  # Set some external lookups
  if 'mask' in d and d['mask'] is not None and d['mask'] is not "":
    path = load_path(d['mask'], directory=directory)
    with open(path) as infile:
      imageset.external_lookup.mask.filename = path
      imageset.external_lookup.mask.data = ImageBool(pickle.load(infile))
  if 'gain' in d and d['gain'] is not None and d['gain'] is not "":
    path = load_path(d['gain'], directory=directory)
    with open(path) as infile:
      imageset.external_lookup.gain.filename = path
      imageset.external_lookup.gain.data = ImageDouble(pickle.load(infile))
  if 'pedestal' in d and d['pedestal'] is not None and d['pedestal'] is not "":
    path = load_path(d['pedestal'], directory=directory)
    with open(path) as infile:
      imageset.external_lookup.pedestal.filename = path
      imageset.external_lookup.pedestal.data = ImageDouble(pickle.load(infile))

  # Get the existing models as dictionaries
  beam_dict = imageset.get_beam(0).to_dict()
  detector_dict = imageset.get_detector(0).to_dict()

  # Set models
  imageset.set_beam(BeamFactory.from_dict(d.get('beam'), beam_dict))
  imageset.set_detector(DetectorFactory.from_dict(d.get('detector'), detector_dict))

  # Return the imageset
  return imageset

def imagesweep_from_dict(d, check_format=True, directory=None):
  '''Construct and image sweep from the dictionary.'''
  from dxtbx.imageset import ImageSetFactory
  from dxtbx.model import BeamFactory, DetectorFactory, GoniometerFactory, ScanFactory
  from dxtbx.serialize.filename import load_path

  # Get the template (required)
  template = load_path(str(d['template']), directory=directory)

  # If the scan isn't set, find all available files
  scan_dict = d.get('scan')
  if scan_dict is None:
    image_range = None
  else:
    image_range = scan_dict.get('image_range')

  # Set the models with the exisiting models as templates
  beam = BeamFactory.from_dict(d.get('beam'))
  goniometer = GoniometerFactory.from_dict(d.get('goniometer'))
  detector = DetectorFactory.from_dict(d.get('detector'))
  scan = ScanFactory.from_dict(d.get('scan'))

  # Construct the sweep
  try:
    sweep = ImageSetFactory.from_template(
      template,
      image_range,
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
      check_format=check_format)[0]
  except Exception:
    indices = range(image_range[0], image_range[1] + 1)
    sweep = ImageSetFactory.make_sweep(
      template,
      indices,
      beam=beam,
      detector=detector,
      goniometer=goniometer,
      scan=scan,
      check_format=check_format)

  # Set some external lookups
  if 'mask' in d and d['mask'] is not None and d['mask'] is not "":
    path = load_path(d['mask'], directory=directory)
    with open(path) as infile:
      sweep.external_lookup.mask.filename = path
      sweep.external_lookup.mask.data = ImageBool(pickle.load(infile))
  if 'gain' in d and d['gain'] is not None and d['gain'] is not "":
    path = load_path(d['gain'], directory=directory)
    with open(path) as infile:
      sweep.external_lookup.gain.filename = path
      sweep.external_lookup.gain.data = ImageDouble(pickle.load(infile))
  if 'pedestal' in d and d['pedestal'] is not None and d['pedestal'] is not "":
    path = load_path(d['pedestal'], directory=directory)
    with open(path) as infile:
      sweep.external_lookup.pedestal.filename = path
      sweep.external_lookup.pedestal.data = ImageDouble(pickle.load(infile))

  # Return the sweep
  return sweep

def imageset_from_dict(d, check_format=True, directory=None):
  ''' Convert the dictionary to a sweep

  Params:
      d The dictionary of parameters

  Returns:
      The sweep

  '''
  # Check the input
  if d == None:
    return None

  # Check the version and id
  if str(d['__id__']) != "imageset":
    raise ValueError("\"__id__\" does not equal \"imageset\"")

  if "filenames" in d:
    return basic_imageset_from_dict(d, directory=directory)
  elif "template" in d:
    return imagesweep_from_dict(d, check_format=check_format, directory=directory)
  else:
    raise TypeError("Unable to deserialize given imageset")
