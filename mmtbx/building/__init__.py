
from __future__ import division
import sys

# XXX this can be moved elsewhere if more appropriate
def run_real_space_annealing (
    xray_structure,
    pdb_hierarchy,
    selection,
    target_map,
    d_min,
    processed_pdb_file=None,
    cif_objects=(),
    resolution_factor=0.25,
    params=None,
    wc=1,
    target_map_rsr=None,
    rsr_after_anneal=False,
    out=None,
    debug=False) :
  """
  Run simulated annealing with a real-space target.  For maximum flexibility,
  a separate map may be used for the initial real-space refinement step.
  """
  import mmtbx.refinement.real_space.individual_sites
  from mmtbx.command_line import real_space_refine
  from mmtbx.dynamics import simulated_annealing
  import mmtbx.restraints
  import mmtbx.maps.utils
  import mmtbx.utils
  import iotbx.phil
  from scitbx.array_family import flex
  if (out is None) :
    out = sys.stdout
  if (params is None) :
    params = iotbx.phil.parse(simulated_annealing.master_params_str).extract()
  if (processed_pdb_file is None) :
    processed_pdb_file = reprocess_pdb(
      pdb_hierarchy=pdb_hierarchy,
      cif_objects=cif_objects,
      crystal_symmetry=xray_structure,
      out=out)
  iselection = selection.iselection()
  xrs = xray_structure
  if (target_map_rsr is None) :
    target_map_rsr = target_map
  grm = real_space_refine.get_geometry_restraints_manager(
    processed_pdb_file = processed_pdb_file,
    xray_structure     = xrs)
  selection_within = xrs.selection_within(
      radius    = 5,#selection_buffer_radius,
      selection = selection)
  selection_in_box = selection.select(selection_within)
  assert (len(iselection) == selection_in_box.count(True))
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure   = xrs,
    pdb_hierarchy    = pdb_hierarchy,
    map_data         = target_map,
    map_data_2       = target_map_rsr,
    selection        = selection_within,
    box_cushion      = 2)
  new_unit_cell = box.xray_structure_box.unit_cell()
  geo_box = grm.geometry.select(box.selection_within).discard_symmetry(
    new_unit_cell=new_unit_cell)
  target_map_box = box.map_box
  target_map_rsr_box = box.map_box_2
  sites_cart_box = box.xray_structure_box.sites_cart()
  selection_all_box = flex.bool(sites_cart_box.size(), True)
  # we want to keep the rest of the atoms more or less fixed, otherwise
  # they simply fly apart during annealing.  this will have obvious
  # implications for the fragment being refined, as it will remain tightly
  # tethered to the atoms on either side (if any).
  # XXX may need a more flexible approach
  reference_selection = (~selection_in_box).iselection()
  reference_sites = sites_cart_box.select(reference_selection)
  geo_box.generic_restraints_manager.reference_manager.\
    add_coordinate_restraints(
      sites_cart = reference_sites,
      selection  = reference_selection,
      sigma      = 0.05)
  rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
    target_map                  = target_map_box,
    selection                   = selection_in_box,
    real_space_gradients_delta  = d_min*resolution_factor,
    max_iterations              = 150,
    geometry_restraints_manager = geo_box)
  refined = mmtbx.refinement.real_space.individual_sites.refinery(
    refiner                  = rsr_simple_refiner,
    xray_structure           = box.xray_structure_box,
    start_trial_weight_value = 1.0,
    rms_bonds_limit          = 0.01,
    rms_angles_limit         = 1.5)
  print >> out, ""
  print >> out, \
    "  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
      (refined.rms_bonds_final, refined.rms_angles_final)
  sites_cart_box_refined = refined.sites_cart_result
  box.xray_structure_box.set_sites_cart(sites_cart_box_refined)
  box_restraints_manager = mmtbx.restraints.manager(
    geometry      = geo_box,
    normalization = True)
  simulated_annealing.run(
    params             = params,
    fmodel             = None,
    xray_structure     = box.xray_structure_box,
    real_space         = True,
    target_map         = target_map_box,
    restraints_manager = box_restraints_manager,
    wx                 = refined.weight_final,# * 0.1,
    wc                 = wc,
    log                = out,
    verbose            = True)
  sites_cart_box_refined = box.xray_structure_box.sites_cart()
  # XXX this seems to work poorly
  if (rsr_after_anneal) :
    rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = target_map_rsr_box, # XXX which map?
      selection                   = selection_all_box,
      real_space_gradients_delta  = d_min*resolution_factor,
      max_iterations              = 150,
      geometry_restraints_manager = geo_box)
    refined = mmtbx.refinement.real_space.individual_sites.refinery(
      refiner                  = rsr_simple_refiner,
      xray_structure           = box.xray_structure_box,
      start_trial_weight_value = 1,
      rms_bonds_limit          = 0.015,
      rms_angles_limit         = 1.5)
    sites_cart_box_refined = refined.sites_cart_result
    print >> out, ""
    print >> out, \
      "  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
        (refined.rms_bonds_final, refined.rms_angles_final)
  sites_cart_box_refined_shifted_back = \
      sites_cart_box_refined + box.shift_to_map_boxed_sites_back
  sites_cart_refined = sites_cart_box_refined_shifted_back.select(
    selection_in_box)
  assert (len(iselection) == len(sites_cart_refined))
  sites_cart_moving = xrs.sites_cart().set_selected(
     iselection, sites_cart_refined)
  if (debug) :
    box.write_ccp4_map()
    box.pdb_hierarchy_box.atoms().set_xyz(box.xray_structure_box.sites_cart())
    box.write_pdb_file("box.pdb")
  return sites_cart_moving

def reprocess_pdb (pdb_hierarchy, crystal_symmetry, cif_objects, out) :
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx.monomer_library import server
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  for cif_object in cif_objects :
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(cif_object=cif_object)
  return pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry),
    crystal_symmetry=crystal_symmetry,
    log=out)
