#!/usr/bin/env python
#
# dxtbx.serialize.load.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def _decode_list(data):
  '''Decode a list to str from unicode. '''
  rv = []
  for item in data:
    if isinstance(item, unicode):
      item = item.encode('utf-8')
    elif isinstance(item, list):
      item = _decode_list(item)
    elif isinstance(item, dict):
      item = _decode_dict(item)
    rv.append(item)
  return rv

def _decode_dict(data):
  ''' Decode a dict to str from unicode. '''
  rv = {}
  for key, value in data.iteritems():
    if isinstance(key, unicode):
      key = key.encode('utf-8')
    if isinstance(value, unicode):
      value = value.encode('utf-8')
    elif isinstance(value, list):
      value = _decode_list(value)
    elif isinstance(value, dict):
      value = _decode_dict(value)
    rv[key] = value
  return rv

def imageset_from_string(string):
  ''' Load the string and return the models.

  Params:
      string The JSON string to load

  Returns:
      The models

  '''
  import json
  from dxtbx.serialize.imageset import imageset_from_dict
  return imageset_from_dict(json.loads(string, object_hook=_decode_dict))

def imageset(filename):
  ''' Load the given JSON file.

  Params:
      infile The input filename

  Returns:
      The models

  '''
  from dxtbx.serialize.filename import temp_chdir
  from os.path import abspath, dirname
  # If the input is a string then open and read from that file
  filename = abspath(filename)
  with temp_chdir(dirname(filename)):
    with open(filename, 'r') as infile:
      return imageset_from_string(infile.read())

def datablock(filename, check_format=True):
  ''' Load a given JSON or pickle file.

  Params:
    filename The input filename

  Returns:
    The datablock

  '''
  from dxtbx.datablock import DataBlockFactory
  return DataBlockFactory.from_serialized_format(
    filename, check_format=check_format)

def crystal(infile):
  ''' Load the given JSON file.

  Params:
      infile The input filename or file object

  Returns:
      The models

  '''
  from dxtbx.serialize.crystal import from_string
  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'r') as infile:
      return from_string(infile.read())

  # Otherwise assume the input is a file and read from it
  else:
    return from_string(infile.read())

def experiment_list(infile, check_format=True):
  ''' Load an experiment list from a serialzied format. '''
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_serialized_format(
    infile, check_format=check_format)
