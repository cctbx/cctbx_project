#!/usr/bin/env python
#
# cctbx.crystal.crystal_model.serialize.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


def load_crystal(infile):
  ''' Load the given JSON file.

  Params:
      infile The input filename or file object

  Returns:
      The models

  '''
  # If the input is a string then open and read from that file
  if isinstance(infile, str):
    with open(infile, 'r') as infile:
      return crystal_from_string(infile.read())

  # Otherwise assume the input is a file and read from it
  else:
    return crystal_from_string(infile.read())


def dump_crystal(obj, outfile, compact=False):
  ''' Dump the given object to file.

  Params:
      obj The crystal to dump
      outfile The output file name or file object
      compact Write in compact representation

  '''
  # If the input is a string then open and write to that file
  if isinstance(outfile, str):
    with open(outfile, 'w') as outfile:
      outfile.write(crystal_to_string(obj, compact))

  # Otherwise assume the input is a file and write to it
  else:
    outfile.write(crystal_to_string(obj, compact))


def crystal_from_string(string):
  ''' Load the string and return the models.

  Params:
      string The JSON string to load

  Returns:
      The models

  '''
  import json
  return crystal_from_dict(json.loads(string))


def crystal_to_string(obj, compact=False):
  ''' Dump the given object to string.

  Params:
      obj The crystal model
      compact Write in compact representation

  Returns:
      The JSON string

  '''
  import json
  from dials.model.serialize.crystal import crystal_to_dict
  from dxtbx.serialize.dump import compact_simple_lists

  # Return as a JSON string
  if compact == False:
    string = json.dumps(crystal_to_dict(obj), indent=2)

    # Hack to make more readable
    string = compact_simple_lists(string)

  else:
    string = json.dumps(crystal_to_dict(obj), separators=(',',':'))

  # Return the string
  return string


def crystal_to_dict(crystal):
  ''' Convert the crystal model to a dictionary

  Params:
      crystal The crystal model

  Returns:
      A dictionary of the parameters

  '''
  from collections import OrderedDict

  # Get the real space vectors
  A = crystal.get_A().inverse()
  real_space_a = (A[0], A[1], A[2])
  real_space_b = (A[3], A[4], A[5])
  real_space_c = (A[6], A[7], A[8])

  # Get the space group Hall symbol
  hall = crystal.get_space_group().info().type().hall_symbol()

  # Get the mosaicity
  mosaicity = crystal.get_mosaicity()

  # Collect the information as a python dictionary
  xl_dict = OrderedDict([
    ('__id__', 'crystal'),
    ('real_space_a', real_space_a),
    ('real_space_b', real_space_b),
    ('real_space_c', real_space_c),
    ('space_group_hall_symbol', hall),
    ('mosaicity', mosaicity)])

  # Add in scan points if present
  if crystal.num_scan_points > 0:
    A_at_scan_points = tuple([crystal.get_A_at_scan_point(i).elems \
                              for i in range(crystal.num_scan_points)])
    xl_dict['A_at_scan_points'] = A_at_scan_points

  return xl_dict


def crystal_from_dict(d):
  ''' Convert the dictionary to a crystal model

  Params:
      d The dictionary of parameters

  Returns:
      The crystal model

  '''
  from cctbx.crystal.crystal_model import crystal_model

  # If None, return None
  if d is None:
    return None

  # Check the version and id
  if str(d['__id__']) != "crystal":
    raise ValueError("\"__id__\" does not equal \"crystal\"")

  # Extract from the dictionary
  real_space_a = d['real_space_a']
  real_space_b = d['real_space_b']
  real_space_c = d['real_space_c']
  # str required to force unicode to ascii conversion
  space_group  = str("Hall:" + d['space_group_hall_symbol'])
  mosaicity    = d['mosaicity']
  xl = crystal_model(real_space_a, real_space_b, real_space_c,
                     space_group_symbol=space_group,
                     mosaicity=mosaicity)

  # Extract scan point setting matrices, if present
  try:
    A_at_scan_points = d['A_at_scan_points']
    xl.set_A_at_scan_points(A_at_scan_points)
  except KeyError:
    pass

  return xl
