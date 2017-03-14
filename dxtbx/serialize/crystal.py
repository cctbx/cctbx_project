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

from __future__ import absolute_import, division


# def from_string(string):
#   ''' Load the string and return the models.

#   Params:
#       string The JSON string to load

#   Returns:
#       The models

#   '''
#   import json
#   return from_dict(json.loads(string))


# def to_string(obj, compact=False):
#   ''' Dump the given object to string.

#   Params:
#       obj The crystal model
#       compact Write in compact representation

#   Returns:
#       The JSON string

#   '''
#   import json
#   from dxtbx.serialize.dump import compact_simple_lists

#   # Return as a JSON string
#   if compact == False:
#     string = json.dumps(to_dict(obj), indent=2)

#     # Hack to make more readable
#     string = compact_simple_lists(string)

#   else:
#     string = json.dumps(to_dict(obj), separators=(',',':'))

#   # Return the string
#   return string



# def from_dict(d):
#   ''' Convert the dictionary to a crystal model

#   Params:
#       d The dictionary of parameters

#   Returns:
#       The crystal model

#   '''
#   from dxtbx.model import Crystal

#   # If None, return None
#   if d is None:
#     return None

#   # Check the version and id
#   if str(d['__id__']) != "crystal":
#     raise ValueError("\"__id__\" does not equal \"crystal\"")

#   # Extract from the dictionary
#   real_space_a = d['real_space_a']
#   real_space_b = d['real_space_b']
#   real_space_c = d['real_space_c']
#   # str required to force unicode to ascii conversion
#   space_group  = str("Hall:" + d['space_group_hall_symbol'])
#   xl = Crystal(real_space_a, real_space_b, real_space_c,
#                      space_group_symbol=space_group)
#   # New parameters for maximum likelihood values
#   try:
#     xl._ML_half_mosaicity_deg = d['ML_half_mosaicity_deg']
#   except KeyError:
#     pass
#   try:
#     xl._ML_domain_size_ang = d['ML_domain_size_ang']
#   except KeyError:
#     pass

#   # Isoforms used for stills
#   try:
#     xl.identified_isoform = d['identified_isoform']
#   except KeyError:
#     pass

#   # Extract scan point setting matrices, if present
#   try:
#     A_at_scan_points = d['A_at_scan_points']
#     xl.set_A_at_scan_points(A_at_scan_points)
#   except KeyError:
#     pass

#   # Extract covariance of B, if present
#   try:
#     cov_B = d['B_covariance']
#     xl.set_B_covariance(cov_B)
#   except KeyError:
#     pass

#   # Extract mosaicity, if present
#   try:
#     mosaicity = d['mosaicity']
#     xl.set_mosaicity(mosaicity)
#   except KeyError:
#     pass

#   return xl
