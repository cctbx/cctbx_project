from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.cspad_pinwheel
#
# Removes all but the central sensors from a CSPAD CBF
#

import dxtbx, sys, os
import libtbx.option_parser
from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict, write_cspad_cbf
from libtbx.utils import Usage
from six.moves import range

def run(argv=None):
  if argv is None:
    argv = sys.argv[1:]

  command_line = libtbx.option_parser.option_parser(
    usage="%s files" % libtbx.env.dispatcher_name).process(args=argv)

  paths = command_line.args
  if len(paths) <= 0:
    raise Usage("No files specified")

  for path in paths:
    # Load the metrology dictionary, containting basis shifts for each item in the hierarchy
    metro = cbf_file_to_basis_dict(path)

    # Remove from the hiearchy all but the central sensors (sensor 1 of each quadrant).
    # Need to remove the sesnor basis shifts and the corresponding asic shifts
    for quad in range(4):
      for sensor in [0,2,3,4,5,6,7]:
        metro.pop((0,quad,sensor))
        for asic in range(2):
          metro.pop((0,quad,sensor,asic))

    # Renumber the sensors to 0 instead of 1
    for key in metro:
      if len(key) == 3:
        detector, quad, sensor = key
        metro[(detector,quad,0)] = metro.pop(key)
      elif len(key) == 4:
        detector, quad, sensor, asic = key
        metro[(detector,quad,0,asic)] = metro.pop(key)

    # Build the tiles dictionary for only sensor 1 of each quadrant.  Rename that sensor to zero.
    img = dxtbx.load(path)
    tiles = {}
    for quad in range(4):
      src_sensor = 1
      dest_sensor = 0
      for asic in range(2):
        tiles[(0,quad,dest_sensor,asic)] = img.get_raw_data()[(quad*16)+(src_sensor*2)+asic] # FIXME get the panel ID from dxtbx

    destpath = os.path.splitext(path)[0] + "_pinwheel.cbf"

    hierarchy = img.get_detector().hierarchy()
    beam = img.get_beam()

    # write the result.  Have to call abs on the root distance because of a bug with the hierarchy matrix.
    write_cspad_cbf(tiles, metro, 'cbf', None, destpath, beam.get_wavelength(), hierarchy.get_distance())

if (__name__ == "__main__") :
  run(sys.argv[1:])
