from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.gain_map_cbf
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT
#
# $Id
#
# Convert an LCLS gain mask to a CBF for use with converting XTC streams
# to gain-corrected CSPAD CBF images
#

import sys, os, numpy

from libtbx.option_parser import option_parser
from scitbx.array_family import flex
from libtbx.utils import Usage
from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict, write_cspad_cbf

def run(args):
  command_line = (option_parser()
                  .option("-o", "--output_filename",
                          action="store",
                          type="string",
                          help="Filename for the output cbf file",
                          default="gain_map.cbf")
                  .option("-m", "--metrology",
                          action="store",
                          type="string",
                          help="CBF or DEF file",
                          default=None)
                  .option("-d", "--distance",
                          action="store",
                          type="int",
                          help="Detector distance put into the gain cbf file. Not needed for processing.",
                          default="0")
                  .option("-w", "--wavelength",
                          action="store",
                          type="float",
                          help="Incident beam wavelength put into the gain cbf file. Not needed for processing.",
                          default="0")
                     ).process(args=args)

  output_filename = command_line.options.output_filename
  metrology = command_line.options.metrology
  assert metrology is not None and os.path.isfile(metrology)

  args = command_line.args

  assert len(args) == 1
  if args[0].endswith('.txt') or args[0].endswith('.gain'):
    raw_data = numpy.loadtxt(args[0])
    assert raw_data.shape in [(5920, 388), (11840, 194)]
    tiles = convert_detector(raw_data)
  else:
    raise Usage("Gain input file should be a text file with extension .txt or .gain")

  metro = cbf_file_to_basis_dict(metrology)
  write_cspad_cbf(tiles, metro, 'cbf', None, output_filename,
                  command_line.options.wavelength, command_line.options.distance)



def convert_detector(raw_data):
  # https://confluence.slac.stanford.edu/display/PCDS/CSPad+metrology+and+calibration+files%2C+links
  data3d = []
  if raw_data.shape == (5920,388):
    raise NotImplementedError("Please contact the authors.  This is an older gain map format.")

  asic_start = 0

  # input dictionary for writing the CBF
  tiles = {}

  # iterate through and extract the data, converting it to asic-sized chunks
  for i_quad in range(4):
    asic_size = 185 * 194
    section_size = asic_size * 4
    quad_start = i_quad * section_size * 4
    quad_asics = []

    # the data is laid out by quadrants, and then by 4 sections per quadrant, each with 2 sensors,
    # each with two asics
    for i_2x2 in range(4):
      for i_sensor in range(2):
        asic_end = asic_start + 185
        a = raw_data[asic_start:asic_end, :]
        asic_start = asic_end
        tiles[(0,i_quad,(i_2x2*2)+i_sensor,0)] = flex.int(a.astype(numpy.int32))

        asic_end = asic_start + 185
        b = raw_data[asic_start:asic_end, :]
        asic_start = asic_end
        tiles[(0,i_quad,(i_2x2*2)+i_sensor,1)] = flex.int(b.astype(numpy.int32))

  return tiles

if __name__ == '__main__':
  run(sys.argv[1:])
