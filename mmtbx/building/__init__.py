
from __future__ import division
from libtbx import group_args
import sys

#-----------------------------------------------------------------------
# MODEL UTILITIES
def residue_id_str (residue) :
  chain = residue.parent().parent()
  resid = residue.parent().resid()
  id_str="%2s %3s%5s"%(chain.id,residue.resname,resid)
  return id_str

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

def get_nearby_water_selection (
    pdb_hierarchy,
    xray_structure,
    selection,
    selection_buffer_radius=5) :
  sel_cache = pdb_hierarchy.atom_selection_cache()
  water_selection = sel_cache.selection("resname HOH")
  selection_within = xray_structure.selection_within(
    radius=selection_buffer_radius,
    selection=selection)
  return (selection_within & water_selection)

def iter_residue_groups (hierarchy) :
  for chain in hierarchy.only_model().chains() :
    for residue_group in chain.residue_groups() :
      yield residue_group

#-----------------------------------------------------------------------
# MAP-RELATED
class local_density_quality (object) :
  def __init__ (O,
      fofc_map,
      two_fofc_map,
      atom_selection,
      xray_structure,
      radius=2.5) :
    from cctbx import maptbx
    from scitbx.array_family import flex
    assert ((len(fofc_map) == len(two_fofc_map)) and
            (fofc_map.focus() == two_fofc_map.focus()))
    O.atom_selection = atom_selection # assumes no HD already
    assert (len(O.atom_selection) > 0)
    O.sites_cart = xray_structure.sites_cart().select(O.atom_selection)
    O.map_sel = maptbx.grid_indices_around_sites(
      unit_cell=xray_structure.unit_cell(),
      fft_n_real=fofc_map.focus(),
      fft_m_real=fofc_map.all(),
      sites_cart=O.sites_cart,
      site_radii=flex.double(O.sites_cart.size(), radius))
    assert (len(O.map_sel) > 0)
    O.fofc_map_sel = fofc_map.select(O.map_sel)
    O.two_fofc_map_sel = two_fofc_map.select(O.map_sel)
    O.difference_map_levels = flex.double()
    for site_cart in O.sites_cart :
      site_frac = xray_structure.unit_cell().fractionalize(site_cart=site_cart)
      fofc_map_value = fofc_map.tricubic_interpolation(site_frac)
      O.difference_map_levels.append(fofc_map_value)

  def number_of_atoms_in_difference_holes (O, sigma_cutoff=-2.5) :
    return (O.difference_map_levels <= sigma_cutoff).count(True)

  def fraction_of_nearby_grid_points_above_cutoff (O, sigma_cutoff=2.5) :
    return (O.fofc_map_sel >= sigma_cutoff).count(True) / len(O.fofc_map_sel)

def get_difference_maps (fmodel, resolution_factor=0.25) :
  fofc_coeffs = fmodel.map_coefficients(map_type="mFo-DFc",
    exclude_free_r_reflections=True)
  fofc_fft = fofc_coeffs.fft_map(
    resolution_factor=resolution_factor)
  fofc_map = fofc_fft.apply_sigma_scaling().real_map_unpadded()
  two_fofc_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True)
  two_fofc_fft = two_fofc_coeffs.fft_map(resolution_factor=resolution_factor)
  two_fofc_map = two_fofc_fft.apply_sigma_scaling().real_map_unpadded()
  return two_fofc_map, fofc_map

# XXX this should probably be merged with something else, e.g. the
# local_density_quality class
def get_model_map_stats (
    selection,
    target_map,
    model_map,
    unit_cell,
    sites_cart,
    pdb_atoms,
    local_sampling=False) :
  """
  Collect basic statistics for a model map and some target map (usually an
  mFo-DFc map), including CC, mean, and minimum density at the atomic
  positions.
  """
  assert (len(target_map) == len(model_map))
  iselection = selection
  if (type(selection).__name__ == 'bool') :
    iselection = selection.iselection()
  from scitbx.array_family import flex
  sites_cart_refined = sites_cart.select(selection)
  sites_selected = flex.vec3_double()
  map1 = flex.double()
  map2 = flex.double()
  min_density = sys.maxint
  sum_density = n_sites = 0
  worst_atom = None
  # XXX I'm not sure the strict density cutoff is a good idea here
  for i_seq, xyz in zip(iselection, sites_cart_refined) :
    if (pdb_atoms[i_seq].element.strip() != "H") :
      sites_selected.append(xyz)
      site_frac = unit_cell.fractionalize(site_cart=xyz)
      target_value = target_map.tricubic_interpolation(site_frac)
      if (target_value < min_density) :
        min_density = target_value
        worst_atom = pdb_atoms[i_seq]
      sum_density += target_value
      n_sites += 1
      if (not local_sampling) :
        map1.append(target_value)
        map2.append(model_map.tricubic_interpolation(site_frac))
  assert (n_sites > 0)
  if (local_sampling) :
    from cctbx import maptbx
    map_sel = maptbx.grid_indices_around_sites(
      unit_cell=unit_cell,
      fft_n_real=target_map.focus(),
      fft_m_real=target_map.all(),
      sites_cart=sites_selected,
      site_radii=flex.double(sites_selected.size(), 1.0))
    map1 = target_map.select(map_sel)
    map2 = model_map.select(map_sel)
  assert (len(map1) > 0) and (len(map1) == len(map2))
  cc = flex.linear_correlation(x=map1, y=map2).coefficient()
  return group_args(
    cc=cc,
    min=min_density,
    mean=sum_density/n_sites)

#-----------------------------------------------------------------------
# OPTIMIZATION
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
    reference_sigma=0.5, # XXX if this is too tight, nothing moves
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
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure   = xrs,
    pdb_hierarchy    = pdb_hierarchy,
    map_data         = target_map,
    map_data_2       = target_map_rsr,
    selection        = selection_within,
    box_cushion      = 2)
  selection_in_box = selection.select(box.selection_within)
  assert (len(iselection) == selection_in_box.count(True))
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
  if (reference_sigma is not None) :
    assert (reference_sigma > 0)
    reference_selection = (~selection_in_box).iselection()
    reference_sites = sites_cart_box.select(reference_selection)
    geo_box.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
        sites_cart = reference_sites,
        selection  = reference_selection,
        sigma      = reference_sigma)
  rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
    target_map                  = target_map_rsr_box,
    selection                   = selection_all_box,
    real_space_gradients_delta  = d_min*resolution_factor,
    max_iterations              = 150,
    geometry_restraints_manager = geo_box)
  refined = mmtbx.refinement.real_space.individual_sites.refinery(
    refiner                  = rsr_simple_refiner,
    xray_structure           = box.xray_structure_box,
    start_trial_weight_value = 1.0,
    rms_bonds_limit          = 0.01,
    rms_angles_limit         = 1.0)
  print >> out, ""
  print >> out, \
    "  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
      (refined.rms_bonds_final, refined.rms_angles_final)
  sites_cart_box_refined = refined.sites_cart_result
  box.xray_structure_box.set_sites_cart(sites_cart_box_refined)
  if (debug) :
    box.write_pdb_file("box_start.pdb")
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
      selection                   = selection_in_box,
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
