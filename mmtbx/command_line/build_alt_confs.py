
from __future__ import division
from libtbx.str_utils import make_header, make_sub_header
import sys

def master_phil () :
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string="""
selection = None
  .type = atom_selection
nproc = Auto
  .type = int
building {
  include scope mmtbx.building.disorder.single_residue.params_str
}
output {
  file_name = refined.pdb
    .type = path
}
prefilter {
  include scope mmtbx.building.disorder.filter_params_str
}
""")

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  from mmtbx.building.disorder import single_residue
  from mmtbx.building import disorder
  from mmtbx import restraints
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=True,
    create_fmodel=True,
    out=out)
  params = cmdline.params
  fmodel = cmdline.fmodel
  hierarchy = cmdline.pdb_hierarchy
  geometry_restraints_manager = cmdline.geometry
  hd_sel = fmodel.xray_structure.hd_selection()
  n_hydrogen = hd_sel.count(True)
  if (n_hydrogen > 0) :
    print >> out, "WARNING: %d hydrogen atoms will be removed!" % n_hydrogen
    non_hd_sel = ~hd_sel
    hierarchy = hierarchy.select(non_hd_sel)
    xray_structure = fmodel.xray_structure.select(non_hd_sel)
    fmodel.update_xray_structure(xray_structure)
    geometry_restraints_manager = geometry_restraints_manager.select(non_hd_sel)
  sele_cache = hierarchy.atom_selection_cache()
  selection = None
  if (params.selection is not None) :
    selection = sele_cache.selection(params.selection)
  make_header("Build cycle", out=out)
  candidate_residues = disorder.filter_before_build(
    pdb_hierarchy=hierarchy,
    fmodel=fmodel,
    geometry_restraints_manager=geometry_restraints_manager,
    selection=selection,
    params=params.prefilter,
    verbose=True,
    log=out)
  restraints_manager = restraints.manager(
    geometry=geometry_restraints_manager,
    normalization=True)
  make_sub_header("Finding alternate conformations", out=out)
  fmodel.info().show_rfactors_targets_scales_overall(out=out)
  new_residues = single_residue.find_all_alternates(
    residues=candidate_residues,
    pdb_hierarchy=hierarchy,
    restraints_manager=restraints_manager,
    fmodel=fmodel,
    params=params.building,
    nproc=params.nproc,
    log=out).results
  single_residue.process_results(
    pdb_hierarchy=hierarchy,
    fmodel=fmodel,
    residues_in=candidate_residues,
    new_residues=new_residues,
    params=params.building,
    log=out)
  # TODO real-space refinement of multi-conformer model
  f = open(params.output.file_name, "w")
  f.write(hierarchy.as_pdb_string(crystal_symmetry=fmodel.xray_structure))
  f.close()
  print >> out, "Wrote %s" % params.output.file_name

if (__name__ == "__main__") :
  run(sys.argv[1:])
