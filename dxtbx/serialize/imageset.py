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
  if filename is None:
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

  # If this is an imageset then return a list of filenames
  if isinstance(imageset, ImageSweep):
    return imagesweep_to_dict(imageset)
  elif isinstance(imageset, ImageSet):
    return basic_imageset_to_dict(imageset)
  else:
    raise TypeError("Unknown ImageSet Type")

def basic_imageset_from_dict(d):
  ''' Construct an ImageSet class from the dictionary.'''
  from dxtbx.model import BeamFactory, DetectorFactory
  from dxtbx.imageset import ImageSetFactory
  from dxtbx.serialize.filename import load_path

  # Get the filename list and create the imageset
  filenames = map(lambda p: load_path(p), map(str, d['filenames']))
  imageset = ImageSetFactory.new(filenames)[0]

  # Set some external lookups
  if 'mask' in d and d['mask'] is not None:
    with open(d['mask']) as infile:
      imageset.external_lookup.mask.filename = d['mask']
      imageset.external_lookup.mask.data = pickle.load(infile)
  if 'gain' in d and d['gain'] is not None:
    with open(d['gain']) as infile:
      imageset.external_lookup.gain.filename = d['gain']
      imageset.external_lookup.gain.data = pickle.load(infile)
  if 'pedestal' in d and d['pedestal'] is not None:
    with open(d['pedestal']) as infile:
      imageset.external_lookup.pedestal.filename = d['pedestal']
      imageset.external_lookup.pedestal.data = pickle.load(infile)

  # Get the existing models as dictionaries
  beam_dict = imageset.get_beam(0).to_dict()
  detector_dict = imageset.get_detector(0).to_dict()

  # Set models
  imageset.set_beam(BeamFactory.from_dict(d.get('beam'), beam_dict))
  imageset.set_detector(DetectorFactory.from_dict(d.get('detector'), detector_dict))

  # Return the imageset
  return imageset

def imagesweep_from_dict(d, check_format=True):
  '''Construct and image sweep from the dictionary.'''
  from dxtbx.imageset import ImageSetFactory
  from dxtbx.model import BeamFactory, DetectorFactory, GoniometerFactory, ScanFactory
  from dxtbx.serialize.filename import load_path

  # Get the template (required)
  template = load_path(str(d['template']))

  # If the scan isn't set, find all available files
  scan_dict = d.get('scan')
  if scan_dict is None:
    image_range = None
  else:
    image_range = scan_dict.get('image_range')

  # Construct the sweep
  try:
    sweep = ImageSetFactory.from_template(
      template, image_range, check_format=check_format)[0]

    # Get the existing models as dictionaries
    beam_dict = sweep.get_beam().to_dict()
    gonio_dict = sweep.get_goniometer().to_dict()
    detector_dict = sweep.get_detector().to_dict()
    scan_dict = sweep.get_scan().to_dict()
  except Exception:
    indices = range(image_range[0], image_range[1] + 1)
    sweep = ImageSetFactory.make_sweep(
      template, indices, check_format=False)
    beam_dict = None
    gonio_dict = None
    detector_dict = None
    scan_dict = None

  # Set some external lookups
  if 'mask' in d and d['mask'] is not None:
    with open(d['mask']) as infile:
      sweep.external_lookup.mask.filename = d['mask']
      sweep.external_lookup.mask.data = pickle.load(infile)
  if 'gain' in d and d['gain'] is not None:
    with open(d['gain']) as infile:
      sweep.external_lookup.gain.filename = d['gain']
      sweep.external_lookup.gain.data = pickle.load(infile)
  if 'pedestal' in d and d['pedestal'] is not None:
    with open(d['pedestal']) as infile:
      sweep.external_lookup.pedestal.filename = d['pedestal']
      sweep.external_lookup.pedestal.data = pickle.load(infile)

  # Set the models with the exisiting models as templates
  sweep.set_beam(BeamFactory.from_dict(d.get('beam'), beam_dict))
  sweep.set_goniometer(GoniometerFactory.from_dict(d.get('goniometer'), gonio_dict))
  sweep.set_detector(DetectorFactory.from_dict(d.get('detector'), detector_dict))
  sweep.set_scan(ScanFactory.from_dict(d.get('scan'), scan_dict))

  # Return the sweep
  return sweep

def imageset_from_dict(d, check_format=True):
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
    return basic_imageset_from_dict(d)
  elif "template" in d:
    return imagesweep_from_dict(d, check_format=check_format)
  else:
    raise TypeError("Unable to deserialize given imageset")

