
from __future__ import division
from libtbx import group_args, adopt_init_args
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

class box_build_refine_base (object) :
  """
  Base class for handling functions associated with rebuilding and refinement
  of atoms in a box, with the building methods implemented elsewhere.
  """
  def __init__ (O,
      pdb_hierarchy,
      xray_structure,
      processed_pdb_file,
      target_map,
      selection,
      d_min,
      out,
      cif_objects=(),
      resolution_factor=0.25,
      target_map_rsr=None,
      debug=False) :
    adopt_init_args(O, locals())
    from mmtbx.command_line import real_space_refine
    import mmtbx.restraints
    import mmtbx.utils
    from scitbx.array_family import flex
    if (O.processed_pdb_file is None) :
      O.processed_pdb_file = reprocess_pdb(
        pdb_hierarchy=pdb_hierarchy,
        cif_objects=cif_objects,
        crystal_symmetry=xray_structure,
        out=out)
    O.iselection = selection.iselection()
    grm = real_space_refine.get_geometry_restraints_manager(
      processed_pdb_file = O.processed_pdb_file,
      xray_structure     = O.xray_structure)
    O.selection_within = xray_structure.selection_within(
      radius    = 5,#selection_buffer_radius,
      selection = selection)
    O.box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure   = xray_structure,
      pdb_hierarchy    = pdb_hierarchy,
      map_data         = target_map,
      map_data_2       = target_map_rsr,
      selection        = O.selection_within,
      box_cushion      = 2)
    O.hd_selection_box = O.box.xray_structure_box.hd_selection()
    O.n_sites_box = O.box.xray_structure_box.sites_cart().size()
    O.target_map_box = O.box.map_box
    O.target_map_rsr_box = O.box.map_box_2
    if (O.target_map_rsr_box is None) :
      O.target_map_rsr_box = O.target_map_box
    new_unit_cell = O.box.xray_structure_box.unit_cell()
    O.geo_box = grm.geometry.select(O.box.selection_within).discard_symmetry(
      new_unit_cell=new_unit_cell)
    O.reference_manager_box = \
      O.geo_box.generic_restraints_manager.reference_manager
    O.box_restraints_manager = mmtbx.restraints.manager(
      geometry      = O.geo_box,
      normalization = True)
    O.selection_all_box = flex.bool(O.n_sites_box, True)
    O.selection_in_box = selection.select(O.selection_within).iselection()
    O.others_in_box = flex.bool(O.n_sites_box, True).set_selected(
      O.selection_in_box, False).iselection()
    assert (len(O.selection_in_box.intersection(O.others_in_box)) == 0)
    O.box_selected_hierarchy = O.box.pdb_hierarchy_box.select(
      O.selection_in_box)
    O.box_selected_hierarchy.atoms().reset_i_seq()
    sites_cart_box = O.box.xray_structure_box.sites_cart()
    O.atom_id_mapping = {}
    for atom in O.box_selected_hierarchy.atoms() :
      O.atom_id_mapping[atom.id_str()] = atom.i_seq

  def get_selected_sites (O, selection=None, hydrogens=True) :
    from scitbx.array_family import flex
    if (selection is None) :
      selection = O.selection_in_box
    if (type(selection).__name__ == 'size_t') :
      selection = flex.bool(O.n_sites_box, False).set_selected(selection, True)
    if (not hydrogens) :
      selection &= ~O.hd_selection_box
    return O.box.xray_structure_box.sites_cart().select(selection)

  def mean_density_at_sites (O, selection=None) :
    """
    Compute the mean density of the target map at coordinates for non-hydrogen
    atoms.
    """
    if (selection is None) :
      selection = O.selection_in_box
    unit_cell = O.box.xray_structure_box.unit_cell()
    sum_density = n_heavy = 0
    for site in O.get_selected_sites(selection=selection, hydrogens=False) :
      site_frac = unit_cell.fractionalize(site)
      sum_density += O.target_map_box.tricubic_interpolation(site_frac)
      n_heavy += 1
    assert (n_heavy > 0)
    return sum_density / n_heavy

  def restrain_atoms (O, selection, reference_sigma) :
    """
    Apply harmonic reference restraints to the selected atoms, wiping out
    any previously existing harmonic restraints.
    """
    assert (reference_sigma > 0)
    sites_cart_box = O.box.xray_structure_box.sites_cart()
    reference_sites = sites_cart_box.select(selection)
    O.reference_manager_box.remove_coordinate_restraints()
    O.reference_manager_box.add_coordinate_restraints(
        sites_cart = reference_sites,
        selection  = selection,
        sigma      = reference_sigma)

  def real_space_refine (O, selection=None) :
    """
    Run simple real-space refinement, returning the optimized weight for the
    map target (wx).
    """
    import mmtbx.refinement.real_space.individual_sites
    if (selection is None) :
      selection = O.selection_all_box
    rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = O.target_map_rsr_box,
      selection                   = selection,
      real_space_gradients_delta  = O.d_min*O.resolution_factor,
      max_iterations              = 150,
      geometry_restraints_manager = O.geo_box)
    refined = mmtbx.refinement.real_space.individual_sites.refinery(
      refiner                  = rsr_simple_refiner,
      xray_structure           = O.box.xray_structure_box,
      start_trial_weight_value = 1.0,
      rms_bonds_limit          = 0.01,
      rms_angles_limit         = 1.0)
    O.update_coordinates(refined.sites_cart_result)
    print >> O.out, ""
    print >> O.out, \
      "  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
        (refined.rms_bonds_final, refined.rms_angles_final)
    return refined.weight_final

  def update_coordinates (O, sites_cart) :
    """
    Set coordinates for the boxed xray_structure and pdb_hierarchy objects.
    """
    O.box.xray_structure_box.set_sites_cart(sites_cart)
    O.box.pdb_hierarchy_box.atoms().set_xyz(sites_cart)

  def update_original_coordinates (O) :
    """
    Copies the new coordinates from the box to the selected atoms in the
    original xray_structure object, and returns the (complete) list of
    coordinates.
    """
    sites_cart_box_refined_shifted_back = \
      O.box.xray_structure_box.sites_cart() + \
      O.box.shift_to_map_boxed_sites_back
    sites_cart_refined = sites_cart_box_refined_shifted_back.select(
      O.selection_in_box)
    assert (len(O.iselection) == len(sites_cart_refined))
    sites_cart_moving = O.xray_structure.sites_cart().set_selected(
       O.iselection, sites_cart_refined)
    return sites_cart_moving

  def geometry_minimization (O,
      selection=None,
      bond=True,
      nonbonded=True,
      angle=True,
      dihedral=True,
      chirality=True,
      planarity=True,
      generic_restraints=True,
      max_number_of_iterations=500,
      number_of_macro_cycles=5,
      rmsd_bonds_termination_cutoff  = 0,
      rmsd_angles_termination_cutoff = 0) :
    """
    Perform geometry minimization on a selection of boxed atoms.
    """
    from mmtbx.command_line import geometry_minimization
    from cctbx import geometry_restraints
    from scitbx.array_family import flex
    import scitbx.lbfgs
    if (selection is None) :
      selection = O.selection_in_box
    if (type(selection).__name__ == 'size_t') :
      selection = flex.bool(O.n_sites_box, False).set_selected(selection, True)
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_number_of_iterations)
    geometry_restraints_flags = geometry_restraints.flags.flags(
      bond               = bond,
      nonbonded          = nonbonded,
      angle              = angle,
      dihedral           = dihedral,
      chirality          = chirality,
      planarity          = planarity,
      generic_restraints = generic_restraints)
    site_labels = O.box.xray_structure_box.scatterers().extract_labels()
    for i in xrange(number_of_macro_cycles):
      sites_cart = O.box.xray_structure_box.sites_cart()
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True)
      minimized = geometry_minimization.lbfgs(
        sites_cart                  = sites_cart,
        geometry_restraints_manager = O.geo_box,
        geometry_restraints_flags   = geometry_restraints_flags,
        lbfgs_termination_params    = lbfgs_termination_params,
        lbfgs_exception_handling_params = exception_handling_params,
        sites_cart_selection        = selection,
        rmsd_bonds_termination_cutoff = rmsd_bonds_termination_cutoff,
        rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff,
        site_labels=site_labels)
      O.update_coordinates(sites_cart)

class anneal_box  (box_build_refine_base) :
  """
  Class for running simulated annealing against the real-space target, with
  optional minimization.
  """
  def __init__ (O, params, *args, **kwds) :
    O.params = params
    box_build_refine_base.__init__(O, *args, **kwds)

  def run (O, rsr_after_anneal=False) :
    from mmtbx.dynamics import simulated_annealing
    import mmtbx.utils
    wx = O.real_space_refine(selection=O.selection_all_box)
    if (O.debug) :
      O.box.write_pdb_file("box_start.pdb")
    states_collector = None
    if (O.debug) :
      states_collector = mmtbx.utils.states(
        xray_structure=O.box.xray_structure_box,
        pdb_hierarchy=O.box.pdb_hierarchy_box)
    simulated_annealing.run(
      params             = O.params,
      fmodel             = None,
      xray_structure     = O.box.xray_structure_box,
      real_space         = True,
      target_map         = O.target_map_box,
      restraints_manager = O.box_restraints_manager,
      wx                 = wx,
      wc                 = 1.0,
      log                = O.out,
      verbose            = True,
      states_collector   = states_collector)
    if (states_collector is not None) :
      states_collector.write("box_traj.pdb")
    if (rsr_after_anneal) :
      O.real_space_refine(selection=O.selection_in_box)
    return O.update_original_coordinates()

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
  from mmtbx.dynamics import simulated_annealing
  import iotbx.phil
  if (params is None) :
    params = iotbx.phil.parse(simulated_annealing.master_params_str).extract()
  O = anneal_box(
    params=params,
    xray_structure=xray_structure,
    pdb_hierarchy=pdb_hierarchy,
    selection=selection,
    target_map=target_map,
    d_min=d_min,
    processed_pdb_file=processed_pdb_file,
    cif_objects=cif_objects,
    resolution_factor=resolution_factor,
    target_map_rsr=target_map_rsr,
    out=out,
    debug=debug)
  if (reference_sigma is not None) :
    O.restrain_atoms(
      selection=O.others_in_box,
      reference_sigma=reference_sigma)
  return O.run(rsr_after_anneal=rsr_after_anneal)
