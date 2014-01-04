
from __future__ import division
import sys

def master_phil () :
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string="""
selection = None
  .type = atom_selection
nproc = Auto
  .type = int
output {
  file_name = refined.pdb
    .type = path
  verbose = False
    .type = bool
}
include scope mmtbx.building.disorder.single_residue.master_phil_str
""")

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  from mmtbx.building.disorder import single_residue
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=True,
    create_fmodel=True,
    out=out)
  params = cmdline.params
  pdb_hierarchy = single_residue.build_cycle(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    params = params,
    cif_objects=cmdline.cif_objects,
    selection=params.selection,
    nproc=params.nproc,
    verbose=params.output.verbose,
    out=out)
  # TODO real-space refinement of multi-conformer model
  f = open(params.output.file_name, "w")
  f.write(pdb_hierarchy.as_pdb_string(
    crystal_symmetry=cmdline.fmodel.xray_structure))
  f.close()
  print >> out, "Wrote %s" % params.output.file_name

if (__name__ == "__main__") :
  run(sys.argv[1:])
