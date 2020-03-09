from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.helix_sheet_recs_as_pdb_files

import sys
import iotbx.pdb
from libtbx.utils import Sorry

legend = """phenix.helix_sheet_recs_as_pdb_files:
  Given PDB file with HELIX/SHEET records output PDB files corresponding to
  each individual HELIX/SHEET record.

How to run:
  phenix.helix_sheet_recs_as_pdb_files model.pdb

Feedback:
  PAfonine@lbl.gov
  phenixbb@phenix-online.org"""

def run(args):
  if(len(args)!=1): raise Sorry("PDB file must be provided.")
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  h = pdb_inp.construct_hierarchy()
  asc = h.atom_selection_cache()
  sso = pdb_inp.extract_secondary_structure()
  for rec in sso.sheets+sso.helices:
    file_name = "_".join(rec.as_pdb_str().split())
    file_name = file_name[:min(36, len(file_name))]
    file_name += ".pdb"
    sel_list = rec.as_atom_selections()
    assert type(sel_list) == list
    if(len(sel_list) == 1):
      sel_str=sel_list[0]
    else:
      sel_str=" or ".join( ["(%s)"%s for s in rec.as_atom_selections()] )
    sel = asc.selection(string=sel_str)
    h_selected = h.select(sel)
    h_selected.write_pdb_file(file_name=file_name)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
