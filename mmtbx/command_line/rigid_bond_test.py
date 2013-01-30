
from __future__ import division
from libtbx.str_utils import make_header
from libtbx.utils import Usage, Sorry
import libtbx.phil
import sys

master_phil = libtbx.phil.parse("""
model = None
  .type = path
restraints = None
  .type = path
  .multiple = True
""")

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("mmtbx.rigid_bond_test model.pdb")
  from mmtbx.monomer_library import pdb_interpretation
  import mmtbx.restraints
  import mmtbx.model
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="model",
    cif_file_def="restraints")
  params = cmdline.work.extract()
  validate_params(params)
  processed_pdb_file = pdb_interpretation.run(
    args=[params.model] + params.restraints)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = True)
  restraints_manager = mmtbx.restraints.manager(
    geometry = geometry,
    normalization = True)
  model = mmtbx.model.manager(
    xray_structure = processed_pdb_file.xray_structure(),
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    restraints_manager = restraints_manager,
    log                = out)
  make_header("Rigid-bond test", out=out)
  model.show_rigid_bond_test(out=out,
    use_id_str=True,
    prefix="  ")

def validate_params (params) :
  if (params.model is None) :
    raise Sorry("Please specify a PDB file.")
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
