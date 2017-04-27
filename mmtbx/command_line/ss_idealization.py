from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME mmtbx.ss_idealization

from mmtbx.secondary_structure import build as ssb
import iotbx.pdb
import os, sys



def run(args):
  pdb_inp = iotbx.pdb.input(source_info=None,
    file_name=args[0])
  ss_ann = pdb_inp.secondary_structure()
  pdb_h = pdb_inp.construct_hierarchy()
  rm = ssb.substitute_ss(
    real_h=pdb_h,
    xray_structure=pdb_inp.xray_structure_simple(),
    use_plane_peptide_bond_restr=True,
    ss_annotation=ss_ann)
  iotbx.pdb.write_whole_pdb_file(
    file_name="%s_ss_ideal.pdb" % args[0],
    output_file=None,
    processed_pdb_file=None,
    pdb_hierarchy=pdb_h,
    crystal_symmetry=None,
    ss_annotation=ss_ann,
    append_end=True,
    atoms_reset_serial_first_value=None,
    link_records=None)

if __name__ == "__main__" :
  run(sys.argv[1:])