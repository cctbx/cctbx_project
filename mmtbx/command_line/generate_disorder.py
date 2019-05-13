
from __future__ import division
from __future__ import print_function
from libtbx.utils import Sorry
from libtbx.str_utils import make_header
from libtbx import Auto, adopt_init_args
from libtbx import easy_mp
import os.path
import random
import time
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
whole_residues = False
  .type = bool
  .help = If the initial selection includes partial residues, expand it to \
    include each residue in its entirety.
selection_delete = None
  .type = atom_selection
  .help = Designates atoms to be removed from the structure before \
    calculating the target map.
target_map = *mFo-DFc 2mFo-DFc
  .type = choice
occ = 0.5
  .type = float
  .help = Partial occupancy for selected atoms for map calculation.
rsr_after_anneal = False
  .type = bool
resolution_factor = 0.25
  .type = float
negate_surrounding_sites = False
  .type = bool
  .help = Set map values to negative around atoms outside of the target \
    selection.
reference_sigma = 0.5
  .type = float
  .help = Sigma for harmonic restraints for neighboring atoms not included \
    in the selection of interest.
wc = 1
  .type = float
  .help = Geometry restraints weight
n_confs = 1
  .type = int
nproc = Auto
  .type = int
random_seed = None
  .type = int
output {
  file_name = disordered.pdb
    .type = path
  include_starting_model = True
    .type = bool
  map_file_name = None
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
  from mmtbx.building import alternate_conformations
  import mmtbx.command_line
  import mmtbx.building
  import iotbx.pdb.hierarchy
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=get_master_phil(),
    process_pdb_file=True,
    create_fmodel=True,
    out=out,
    usage_string="""\
mmtbx.generate_disorder model.pdb data.mtz selection="resname ATP" [occ=0.6]

Perform simulatead annealing against an mFo-DFc map to generate possible
alternate conformations for a selection of atoms.  For development purposes
and experimentation only.
""")
  params = cmdline.params
  fmodel = cmdline.fmodel
  validate_params(params)
  pdb_hierarchy = cmdline.pdb_hierarchy
  make_header("Generating disorder", out=out)
  a_c_p = cmdline.processed_pdb_file.all_chain_proxies
  selection = a_c_p.selection(params.selection)
  if (params.whole_residues):
    selection = iotbx.pdb.atom_selection.expand_selection_to_entire_atom_groups(
      selection=selection,
      pdb_atoms=pdb_hierarchy.atoms())
  n_sel = selection.count(True)
  assert (n_sel > 0)
  print("%d atoms selected" % n_sel, file=out)
  selection_delete = None
  if (params.selection_delete is not None):
    selection_delete = a_c_p.selection(params.selection_delete)
  two_fofc_map, fofc_map = alternate_conformations.get_partial_omit_map(
    fmodel=fmodel.deep_copy(),
    selection=selection,
    selection_delete=selection_delete,
    negate_surrounding=params.negate_surrounding_sites,
    map_file_name=params.output.map_file_name,
    partial_occupancy=params.occ,
    resolution_factor=params.resolution_factor)
  target_map = fofc_map
  if (params.target_map == "2mFo-DFc"):
    target_map = two_fofc_map
  annealer = annealing_manager(
    xray_structure=fmodel.xray_structure,
    pdb_hierarchy=pdb_hierarchy,
    processed_pdb_file=cmdline.processed_pdb_file,
    target_map=target_map,
    two_fofc_map=two_fofc_map,
    d_min=fmodel.f_obs().d_min(),
    params=params,
    selection=selection,
    resolution_factor=params.resolution_factor,
    out=out,
    debug=params.output.debug)
  sites_ref = pdb_hierarchy.atoms().extract_xyz().deep_copy()
  sites_all = easy_mp.pool_map(
    fixed_func=annealer,
    iterable=range(params.n_confs),
    processes=params.nproc)
  ensemble = iotbx.pdb.hierarchy.root()
  if (params.output.include_starting_model):
    sites_all.insert(0, sites_ref)
  rmsds = []
  for i_conf, sites_new in enumerate(sites_all):
    assert (sites_new is not None)
    model = pdb_hierarchy.only_model().detached_copy()
    model.atoms().set_xyz(sites_new)
    model.id = str(i_conf+1)
    rmsd = sites_new.select(selection).rms_difference(
      sites_ref.select(selection))
    print("Model %d: rmsd=%.3f" % (i_conf+1, rmsd), file=out)
    rmsds.append(rmsd)
    ensemble.append_model(model)
  f = open(params.output.file_name, "w")
  f.write(ensemble.as_pdb_string(
    crystal_symmetry=fmodel.xray_structure))
  f.close()
  print("Wrote ensemble model to %s" % params.output.file_name, file=out)
  return rmsds

class annealing_manager(object):
  def __init__(self,
      xray_structure,
      pdb_hierarchy,
      processed_pdb_file,
      target_map,
      two_fofc_map,
      d_min,
      resolution_factor,
      params,
      selection,
      debug=False,
      out=None):
    adopt_init_args(self, locals())

  def __call__(self, i_proc):
    import mmtbx.building
    from scitbx.array_family import flex
    seed = self.params.random_seed
    if (seed is None):
      seed = int(time.time() / os.getpid())
    random.seed(seed)
    flex.set_random_seed(seed)
    sites_new = mmtbx.building.run_real_space_annealing(
      xray_structure=self.xray_structure.deep_copy_scatterers(),
      pdb_hierarchy=self.pdb_hierarchy.deep_copy(),
      processed_pdb_file=self.processed_pdb_file,
      selection=self.selection,
      target_map=self.target_map,
      d_min=self.d_min,
      resolution_factor=self.resolution_factor,
      params=self.params.simulated_annealing,
      target_map_rsr=self.two_fofc_map,
      rsr_after_anneal=self.params.rsr_after_anneal,
      reference_sigma=self.params.reference_sigma,
      wc=self.params.wc,
      out=self.out,
      debug=self.debug)
    return sites_new

def validate_params(params):
  if (params.selection is None):
    raise Sorry("You must specificy an atom selection to anneal.")

if (__name__ == "__main__"):
  run(sys.argv[1:])
