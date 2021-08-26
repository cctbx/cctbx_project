
from __future__ import absolute_import, division, print_function
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

def run(args, out=sys.stdout):
  if (len(args) == 0) or ("--help" in args):
    raise Usage("mmtbx.rigid_bond_test model.pdb")
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
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name = params.model))
  model.process(make_restraints=True)
  model.get_xray_structure()
  model.show_rigid_bond_test(
    out=out,
    use_id_str=True,
    prefix="  ")

def validate_params(params):
  if (params.model is None):
    raise Sorry("Please specify a PDB file.")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
