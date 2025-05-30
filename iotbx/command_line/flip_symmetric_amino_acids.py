"""Flip symmetric amino acids in a model"""
# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.flip_symmetric_amino_acids

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
iotbx.pdb.flip_symmetric_amino_acids model.pdb

"""

  raise Usage(help_msg)

def run(args):
  if len(args) == 0:
    show_usage()
    return
  inp_fn = args[0]
  import time
  t0=time.time()
  pdb_input = iotbx.pdb.input(
      file_name=inp_fn,
      source_info=None,
      raise_sorry_if_format_error=True)
  t0=time.time()
  model = mmtbx.model.manager(
      model_input = pdb_input)
  pdb_h = model.get_hierarchy()
  info = pdb_h.flip_symmetric_amino_acids()
  print(info)
  model.set_sites_cart_from_hierarchy()

  out_fn_prefix = inp_fn
  if inp_fn.endswith(".pdb") or inp_fn.endswith(".cif"):
    out_fn_prefix = inp_fn[:-4]

  if model.input_model_format_cif():
    out_fn = out_fn_prefix + "_iupac.cif"
    txt = model.model_as_mmcif()
  else:
    out_fn = out_fn_prefix + "_iupac.pdb"
    txt = model.model_as_pdb()
  with open(out_fn, 'w') as f:
    f.write(txt)

if (__name__ == "__main__"):
  run(sys.argv[1:])
