"input pdb -> structure factors -> x-plor map files"

import iotbx.pdb
import iotbx.xplor.map
from libtbx.math_utils import ifloor, iceil
import libtbx.option_parser
import sys

def run(args):
  if (len(args) == 0): args = ["--help"]
  command_line = (libtbx.option_parser.option_parser(
    usage="iotbx.python pdb_to_map_simple.py [options] pdb_file...")
    .option(None, "--d_min",
      type="float",
      default=3,
      help="high-resolution limit for structure-factor calculation",
      metavar="FLOAT")
  ).process(args=args)
  d_min = command_line.options.d_min
  assert d_min > 0
  for file_name in command_line.args:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    xray_structure = pdb_inp.xray_structure_simple()
    xray_structure.show_summary()
    print
    print "d_min:", d_min
    f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
    f_calc.show_summary()
    print
    fft_map = f_calc.fft_map()
    n = fft_map.n_real()
    print "unit cell gridding:", n
    fft_map.as_xplor_map(file_name="unit_cell.map")
    print
    block_first = tuple([ifloor(i*0.2) for i in n])
    block_last = tuple([max(f+10, iceil(i*0.7)) for f,i in zip(block_first, n)])
    print "block first:", block_first
    print "block last: ", block_last
    fft_map.as_xplor_map(
      file_name="block.map",
      gridding_first=block_first,
      gridding_last=block_last)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
