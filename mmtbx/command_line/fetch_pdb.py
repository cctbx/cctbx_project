"""Download structure and optionally data from PDB"""
from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.fetch_pdb

from iotbx.cli_parser import run_program
from mmtbx.programs import fetch

def custom_args_proc(cli_parser):
  wf = cli_parser.working_phil.extract()
  if len(cli_parser.namespace.unknown) > 0:
    # print("What is unknown: %s" % cli_parser.namespace.unknown)
    # print("Curr selection: '%s'" % wf.fetch.pdb_ids)
    res = []
    for unk in cli_parser.namespace.unknown:
      res += unk.replace(',',' ').split()
    wf.fetch.pdb_ids = res
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)

if (__name__ == "__main__"):
  run_program(program_class=fetch.Program, custom_process_arguments=custom_args_proc)

