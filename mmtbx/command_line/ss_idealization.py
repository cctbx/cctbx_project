from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.ss_idealization

from mmtbx.secondary_structure import build as ssb
import iotbx.pdb
import os, sys



def run(args):
  pdb_inp = iotbx.pdb.input(source_info=None,
    file_name=args[0])
  ss_ann = pdb_inp.extract_secondary_structure()
  pdb_h = pdb_inp.construct_hierarchy()
  params = ssb.master_phil.extract()
  params.ss_idealization.file_name_before_regularization="before_reg.pdb"
  params.ss_idealization.enabled=True
  rm = ssb.substitute_ss(
    real_h=pdb_h,
    xray_structure=pdb_inp.xray_structure_simple(),
    use_plane_peptide_bond_restr=True,
    log=sys.stdout,
    ss_annotation=ss_ann,
    params=params)
  out_fname = "%s_ss_ideal.pdb" % os.path.basename(args[0])
  iotbx.pdb.write_whole_pdb_file(
    file_name=out_fname,
    output_file=None,
    processed_pdb_file=None,
    pdb_hierarchy=pdb_h,
    crystal_symmetry=None,
    ss_annotation=ss_ann,
    append_end=True,
    atoms_reset_serial_first_value=None,
    link_records=None)
  print("File saved: %s" % out_fname)
  print("All done.")

if __name__ == "__main__" :
  run(sys.argv[1:])