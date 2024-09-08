# LIBTBX_SET_DISPATCHER_NAME phenix.start_molstar
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_molstar
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program, CCTBXParser
from qttbx.programs import start_molstar

class MolstarParser(CCTBXParser):
  """
  This seems to be required to allow running programs without args..
  """
  def __init__(self,*args,**kwargs):
    super().__init__(*args,**kwargs)

  def parse_args(self, args, skip_help = True):
    # Allow running without args
    return super().parse_args(args,skip_help=skip_help)


if __name__ == '__main__':
  run_program(program_class=start_molstar.Program,parser_class=MolstarParser)
