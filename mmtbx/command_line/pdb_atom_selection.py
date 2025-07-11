"""Extract selected atoms from PDB file (useful for experimenting with atom selections)"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

from iotbx.cli_parser import CCTBXParser, run_program
from mmtbx.programs import atom_selection

def custom_args_proc(cli_parser):
  """
  This is going to put inselection positional argument into phil scope of the
  program.
  Useful things are:
  cli_parser.namespace
  cli_parser.data_manager
  cli_parser.namespace.unknown - where these argument is going to be
  cli_parser.working_phil
  """
  wf = cli_parser.working_phil.extract()

  # Since we have phil parameters for these as well, we want to make
  # sure that we don't overwrite them if nothing was provided via
  # command-line.
  if cli_parser.namespace.cryst1_replacement_buffer_layer is not None:
    wf.atom_selection_program.cryst1_replacement_buffer_layer = \
        cli_parser.namespace.cryst1_replacement_buffer_layer
  if cli_parser.namespace.write_pdb_file is not None:
    wf.atom_selection_program.write_pdb_file = \
        cli_parser.namespace.write_pdb_file

  if len(cli_parser.namespace.unknown) > 0:
    # print("What is unknown: %s" % cli_parser.namespace.unknown)
    # print("Curr selection: '%s'" % wf.atom_selection_program.inselection)
    wf.atom_selection_program.inselection = cli_parser.namespace.unknown
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)

class SelectionParser(CCTBXParser):

  def add_default_options(self):
    super(SelectionParser, self).add_default_options()

    # Stick in additional keyword args
    # They going to end up in namespace.cryst1_replacement_buffer_layer etc
    # file_name is obsolet, parser takes care of it putting it in data manager
    # inselections is going to be handled by custom_args_proc function
    # because it is intended to be positional argument
    #
    # !!! This is done here for legacy support and illustrative purposes.
    # Don't do it anywhere else, since we strive to have the same
    # command-line flags across all programs, like --overwrite etc.

    self.add_argument(
      "--cryst1-replacement-buffer-layer",
      action="store",
      type=float,
      help="replace CRYST1 with pseudo unit cell covering the selected"
        " atoms plus a surrounding buffer layer",
      default=None)
    self.add_argument(
      "--write-pdb-file", "--write_pdb_file", "--write_pdb-file", "--write-pdb_file",
      action="store",
      help="write selected atoms to new PDB file",
      default=None)

if (__name__ == "__main__"):
  run_program(parser_class=SelectionParser,
              program_class=atom_selection.Program,
              custom_process_arguments=custom_args_proc)

