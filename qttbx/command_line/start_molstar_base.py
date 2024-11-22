# LIBTBX_SET_DISPATCHER_NAME phenix.start_molstar_base
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_molstar_base
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program, CCTBXParser
from qttbx.programs.start_molstar_base import MolstarBaseProgram

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
  run_program(program_class=MolstarBaseProgram,parser_class=MolstarParser)
