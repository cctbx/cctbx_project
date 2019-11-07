
from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import os.path
import sys

def get_master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_twin_law=True,
    enable_experimental_phases=True,
    enable_pdb_interpretation_params=True,
    enable_stop_for_unknowns=False,
    enable_full_geometry_params=True,
    phil_string="""
selection = None
  .type = atom_selection
target_map = 2mFo-DFc
  .type = str
rsr_after_anneal = False
  .type = bool
resolution_factor = 0.25
  .type = float
reference_sigma = 0.5
  .type = float
  .help = Sigma for harmonic restraints for neighboring atoms not included \
    in the selection of interest.
output {
  file_name = annealed.pdb
    .type = path
  save_map_coeffs = False
    .type = path
  verbose = False
    .type = bool
  debug = False
    .type = bool
}
simulated_annealing {
  include scope mmtbx.dynamics.simulated_annealing.master_params
}
""")

def run(args, out=sys.stdout):
  import mmtbx.command_line
  import mmtbx.building
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=get_master_phil(),
    process_pdb_file=True,
    create_fmodel=True,
    out=out,
    usage_string="""\
mmtbx.anneal_real_space model.pdb data.mtz selection="resname ATP"

Perform simulatead annealing against a real-space (map) target for a selection
of atoms.  For development purposes and experimentation only.
""")
  params = cmdline.params
  fmodel = cmdline.fmodel
  validate_params(params)
  pdb_hierarchy = cmdline.pdb_hierarchy
  a_c_p = cmdline.processed_pdb_file.all_chain_proxies
  selection = a_c_p.selection(params.selection)
  assert selection.count(True) > 0
  map_coeffs = fmodel.map_coefficients(
    map_type=params.target_map,
    exclude_free_r_reflections=True)
  if (params.output.save_map_coeffs):
    file_base = os.path.basename(os.path.splitext(params.output.file_name)[0])
    map_file = file_base + "_%s.mtz" % params.target_map
    mtz = map_coeffs.as_mtz_dataset(column_root_label=params.target_map)
    mtz.mtz_object().write(map_file)
    print("Wrote map coefficients to %s" % map_file, file=out)
  map = map_coeffs.fft_map(resolution_factor=params.resolution_factor
    ).apply_sigma_scaling().real_map_unpadded()
  sites_ref = pdb_hierarchy.atoms().extract_xyz().deep_copy()
  sites_updated = mmtbx.building.run_real_space_annealing(
    xray_structure=fmodel.xray_structure,
    pdb_hierarchy=pdb_hierarchy,
    selection=selection,
    target_map=map,
    d_min=fmodel.f_obs().d_min(),
    processed_pdb_file=cmdline.processed_pdb_file,
    cif_objects=(),
    resolution_factor=params.resolution_factor,
    params=params.simulated_annealing,
    rsr_after_anneal=params.rsr_after_anneal,
    out=out,
    debug=params.output.debug)
  pdb_hierarchy.atoms().set_xyz(sites_updated)
  rmsd = sites_updated.select(selection).rms_difference(
    sites_ref.select(selection))
  print("RMSD (selected atoms) = %.3f" % rmsd, file=out)
  f = open(params.output.file_name, "w")
  f.write(pdb_hierarchy.as_pdb_string(
    crystal_symmetry=fmodel.xray_structure))
  f.close()
  print("Wrote annealed model to %s" % params.output.file_name, file=out)
  return rmsd

def validate_params(params):
  if (params.selection is None):
    raise Sorry("You must specificy an atom selection to anneal.")

if (__name__ == "__main__"):
  run(sys.argv[1:])
