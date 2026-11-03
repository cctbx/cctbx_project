from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.development.the_run_down

# import sys

from mmtbx.programs import next_gen_residue_typing
from iotbx.cli_parser import run_program

# =============================================================================
if __name__ == '__main__':
  #run(sys.argv[1:])
  run_program(program_class=next_gen_residue_typing.Program, hide_parsing_output=True)
