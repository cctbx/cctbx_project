from __future__ import division
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.print_pickle
#

"""
Simple utility for printing the contents of a cctbx.xfel pickle file
"""

from libtbx import easy_pickle
import sys, os
from xfel.detector_formats import detector_format_version as detector_format_function
from xfel.detector_formats import reverse_timestamp
from cctbx import sgtbx # import dependency

args = sys.argv[1:]
if "--break" in args:
  args.remove("--break")
  dobreak = True
else:
  dobreak = False

for path in sys.argv[1:]:
  if not os.path.isfile(path):
    print "Not a file:", path
    continue

  data = easy_pickle.load(path)
  if not isinstance(data, dict):
    print "Not a dictionary pickle"
    continue

  if data.has_key('TIMESTAMP'):
    # this is how FormatPYunspecified guesses the address
    if not "DETECTOR_ADDRESS" in data:
      # legacy format; try to guess the address
      LCLS_detector_address = 'CxiDs1-0|Cspad-0'
      if "DISTANCE" in data and data["DISTANCE"] > 1000:
        # downstream CS-PAD detector station of CXI instrument
        LCLS_detector_address = 'CxiDsd-0|Cspad-0'
    else:
      LCLS_detector_address = data["DETECTOR_ADDRESS"]

    detector_format_version = detector_format_function(
      LCLS_detector_address, reverse_timestamp(data['TIMESTAMP'])[0])
    print "Detector format version:", detector_format_version
  else:
    print "Not an image pickle"

  for key in data:
    if key == 'ACTIVE_AREAS':
      print int(len(data[key])/4), "active areas, first one: ", list(data[key][0:4])
    elif key == 'observations':
      print key, data[key], "Showing unit cell/spacegroup:"
      obs = data[key][0]
      obs.unit_cell().show_parameters()
      obs.space_group().info().show_summary()
    elif key == 'mapped_predictions':
      print key, data[key][0][0], "(only first shown)"
    elif key == 'correction_vectors':
      print key, data[key][0][0], "(only first shown)"
    elif key == "DATA":
      from cctbx.array_family import flex
      print key,"len=%d max=%f min=%f"%(data[key].size(),flex.max(data[key]),flex.min(data[key]))
    else:
      print key, data[key]

  if dobreak:
    print "Entering break. The pickle is loaded in the variable 'data'"
    import pdb; pdb.set_trace()
