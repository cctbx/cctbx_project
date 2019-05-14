
from __future__ import absolute_import, division, print_function
import sys

def get_master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string="""
clash_min = 0.2
  .type = float
""")

def run(args, out=sys.stdout):
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=get_master_phil(),
    require_data=False,
    process_pdb_file=True,
    out=out)
  from mmtbx.building.alternate_conformations import altloc_labeling
  altloc_labeling.show_altloc_clashes(
    pdb_hierarchy=cmdline.pdb_hierarchy,
    xray_structure=cmdline.xray_structure,
    geometry_restraints_manager=cmdline.geometry,
    clash_min=cmdline.params.clash_min,
    out=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
