from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.find_reticular_twin_laws

import sys
from cctbx.sgtbx import reticular_twin_laws
from iotbx.option_parser import option_parser

def run(args, command_name="phenix.find_reticular_twin_laws"):
  command_line = (option_parser(
    usage=command_name+" [options] ",
    description="Example: %s data1.mtz" % command_name)
    .enable_show_defaults()
    .enable_symmetry_comprehensive()
    .option(None, "--max_index",
            action="store", dest="max_index",default=3,type="int",metavar="INT")
    .option(None,"--max_delta", dest="max_delta", default=5.0,type="float",metavar="FLOAT")
  ).process(args=args)
  co = command_line.options
  if (len(args) == 0):
    command_line.parser.show_help()
    return


  max_delta = co.max_delta
  max_index = co.max_index
  xs = command_line.symmetry

  print("----- Input symmetry -----")
  xs.show_summary()
  print()
  print(" Settings:")
  print("   max_delta: ", max_delta)
  print("   max_index: ", max_index)
  print()
  print("----- Finding reticular twin laws -----")
  print()
  print("""



   Reticular twin laws are grouped by their associated sublattice.
   The matrix M acting on the input (base) lattice is listed, as well as the metric R value (%)
   between the symmetrized sublattice and the unsymmetrized sublattice.
   reticular twin laws of twin index 1, are 'normal' twin laws


  """)
  tl = reticular_twin_laws.reticular_twin_laws(xs,max_delta,max_index)
  print("--------------- Twin law listing ---------------")
  print()
  tl.show()





if __name__ == "__main__":
  run(sys.argv[1:])
