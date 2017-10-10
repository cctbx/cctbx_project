# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.flip_symmetric_amino_acids

from __future__ import division
from __future__ import print_function
from libtbx.utils import Usage
import sys
import iotbx.pdb
from iotbx.pdb import write_whole_pdb_file

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
  pdb_h = pdb_input.construct_hierarchy()
  info = pdb_h.flip_symmetric_amino_acids()
  print(info)

  out_fn_prefix = inp_fn
  if inp_fn.endswith(".pdb") or inp_fn.endswith(".cif"):
    out_fn_prefix = inp_fn[:-4]
  out_fn = out_fn_prefix + "_iupac.pdb"

  if hasattr(pdb_input, "extract_secondary_structure"):
    ss_annotation = pdb_input.extract_secondary_structure()
    write_whole_pdb_file(
        file_name=out_fn,
        output_file=None,
        processed_pdb_file=None,
        pdb_hierarchy=pdb_h,
        crystal_symmetry=pdb_input.crystal_symmetry(),
        ss_annotation=ss_annotation,
        atoms_reset_serial_first_value=None,
        link_records=None)
  else:
    # This was a mmcif file, so outputting mmcif
    pdb_h.write_mmcif_file(
        file_name = out_fn,
        crystal_symmetry=pdb_input.crystal_symmetry(),
    )

if (__name__ == "__main__") :
  run(sys.argv[1:])
