from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.development.valence_simple

from iotbx.cli_parser import run_program
from mmtbx.programs import atom_valence_simple

run_program(atom_valence_simple.Program)
