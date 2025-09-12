from __future__ import absolute_import, division, print_function
import sys
from iotbx.cli_parser import run_program
from mmtbx.programs.elbow_services import Program_where_is_that_cif_file

def custom_args_proc(cli_parser):
  wf = cli_parser.working_phil.extract()
  if len(cli_parser.namespace.unknown) > 0:
    print("What is unknown: %s" % cli_parser.namespace.unknown)
    print("Curr selection: '%s'" % wf.where.code)
    res = []
    for unk in cli_parser.namespace.unknown:
      res += unk.replace(',',' ').split()
    wf.where.code = cli_parser.namespace.unknown[0]
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)

if (__name__ == '__main__'):
  print(sys.argv)
  results = run_program(program_class=Program_where_is_that_cif_file,
                        custom_process_arguments=custom_args_proc)
