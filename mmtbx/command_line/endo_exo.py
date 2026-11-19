# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.endo_exo
from iotbx.cli_parser import run_program
from mmtbx.programs import endo_exo

if __name__ == '__main__':
  run_program(program_class=endo_exo.Program)
