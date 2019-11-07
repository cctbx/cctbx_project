from __future__ import absolute_import, division, print_function
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.detector_format_versions
#
# Utility for printing information about cspad detector format versions
#

import sys
from iotbx.detectors.cspad_detector_formats import _detector_format_version_dict

def print_version(version):
  """
  Given a detector format version, print its vital stats
  """
  from xfel.cxi.cspad_ana.cspad_tbx import evt_timestamp
  print("%14s "%version, "%17s"%_detector_format_version_dict[version]['address'], end=' ')
  if _detector_format_version_dict[version]['start_time'] is None:
    print("None                   ", end=' ')
  else:
    print(evt_timestamp((_detector_format_version_dict[version]['start_time'], 0)), end=' ')

  if _detector_format_version_dict[version]['start_time'] is None:
    print("None")
  else:
    print(evt_timestamp((_detector_format_version_dict[version]['end_time'], 0)))

if __name__=='__main__':

  if len(sys.argv) <= 1:
    print("Listing of all known detector formats. Use cxi.detector_format_versions <version> to show the quadrant and unit pixel translations for a given format version")
    print()
    print("Format version   Det. address     Start time              End time")

    for key in sorted(_detector_format_version_dict):
      print_version(key)
  else:
    from spotfinder.applications.xfel.cxi_phil import cxi_versioned_extract
    for version in sys.argv[1:]:
      if version not in _detector_format_version_dict:
        print("Version %s not found.  Did you use quotes?"%version)
        continue
      print("Showing info for %s detector format version"%version)
      print()
      print("Format version   Det. address     Start time              End time")
      print_version(version)
      print()

      phil = cxi_versioned_extract(["distl.detector_format_version=%s"%version])

      print("Quad translations:", phil.distl.quad_translations)
      print("Tile translations:", phil.distl.tile_translations)
