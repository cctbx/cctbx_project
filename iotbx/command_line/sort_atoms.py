"""Sort atoms in residues so they will be in the same order in all residues"""
# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.sort_atoms

from __future__ import absolute_import, division, print_function
from libtbx.utils import Usage
import sys
import iotbx.pdb
import mmtbx.model

master_phil_str = """
file_name = None
  .type = path
  .multiple = False
  .optional = False
  .style = hidden
"""

def show_usage():
  help_msg = """\
iotbx.pdb.sort_atoms model.pdb

Sort atoms in residues so they will be in the same order in all residues.
Also renumbers atoms (atom serial number field 7-11 columns)."""

  raise Usage(help_msg)

def run(args):
  if len(args) == 0:
    show_usage()
    return
  inp_fn = args[0]
  pdb_input = iotbx.pdb.input(
      file_name=inp_fn,
      source_info=None,
      raise_sorry_if_format_error=True)
  model = mmtbx.model.manager(
      model_input = pdb_input)

  out_fn_prefix = inp_fn
  if inp_fn.endswith(".pdb") or inp_fn.endswith(".cif"):
    out_fn_prefix = inp_fn[:-4]
  out_fn = out_fn_prefix + "_sorted"

  txt = ""
  if model.input_model_format_cif():
    out_fn += ".cif"
    txt = model.model_as_mmcif()
  else:
    out_fn += ".pdb"
    txt = model.model_as_pdb()
  with open(out_fn, 'w') as f:
    f.write(txt)

if (__name__ == "__main__"):
  run(sys.argv[1:])
