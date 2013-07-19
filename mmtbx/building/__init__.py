
from __future__ import division
from libtbx import group_args, adopt_init_args
import sys

#-----------------------------------------------------------------------
# MODEL UTILITIES
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
  def __init__ (self,
      fofc_map,
      two_fofc_map,
      atom_selection,
      xray_structure,
      radius=2.5) :
    from cctbx import maptbx
    from scitbx.array_family import flex
    assert ((len(fofc_map) == len(two_fofc_map)) and
            (fofc_map.focus() == two_fofc_map.focus()))
    self.atom_selection = atom_selection # assumes no HD already
    assert (len(self.atom_selection) > 0)
    self.sites_cart = xray_structure.sites_cart().select(self.atom_selection)
    self.map_sel = maptbx.grid_indices_around_sites(
      unit_cell=xray_structure.unit_cell(),
      fft_n_real=fofc_map.focus(),
      fft_m_real=fofc_map.all(),
      sites_cart=self.sites_cart,
      site_radii=flex.double(self.sites_cart.size(), radius))
    assert (len(self.map_sel) > 0)
    self.fofc_map_sel = fofc_map.select(self.map_sel)
    self.two_fofc_map_sel = two_fofc_map.select(self.map_sel)
    self.difference_map_levels = flex.double()
    for site_cart in self.sites_cart :
      site_frac = xray_structure.unit_cell().fractionalize(site_cart=site_cart)
      fofc_map_value = fofc_map.tricubic_interpolation(site_frac)
      self.difference_map_levels.append(fofc_map_value)

  def number_of_atoms_in_difference_holes (self, sigma_cutoff=-2.5) :
    return (self.difference_map_levels <= sigma_cutoff).count(True)

  def fraction_of_nearby_grid_points_above_cutoff (self, sigma_cutoff=2.5) :
    return (self.fofc_map_sel >= sigma_cutoff).count(True) / len(self.fofc_map_sel)

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
  def __init__ (self,
      pdb_hierarchy,
      xray_structure,
      processed_pdb_file,
      target_map,
      selection,
      d_min,
      out,
      cif_objects=(),
      resolution_factor=0.25,
      selection_buffer_radius=5,
      box_cushion=2,
      target_map_rsr=None,
      debug=False) :
    adopt_init_args(self, locals())
    import mmtbx.restraints
    import mmtbx.utils
    from scitbx.array_family import flex
    if (self.processed_pdb_file is None) :
      self.processed_pdb_file = reprocess_pdb(
        pdb_hierarchy=pdb_hierarchy,
        cif_objects=cif_objects,
        crystal_symmetry=xray_structure,
        out=out)
    self.iselection = selection.iselection()
    grm = mmtbx.restraints.manager(
      geometry=self.processed_pdb_file.geometry_restraints_manager(show_energies=False),
      normalization = True)
    self.selection_within = xray_structure.selection_within(
      radius    = selection_buffer_radius,
      selection = selection)
    self.box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure   = xray_structure,
      pdb_hierarchy    = pdb_hierarchy,
      map_data         = target_map,
      map_data_2       = target_map_rsr,
      selection        = self.selection_within,
      box_cushion      = box_cushion)
    self.hd_selection_box = self.box.xray_structure_box.hd_selection()
    self.n_sites_box = self.box.xray_structure_box.sites_cart().size()
    self.target_map_box = self.box.map_box
    self.target_map_rsr_box = self.box.map_box_2
    if (self.target_map_rsr_box is None) :
      self.target_map_rsr_box = self.target_map_box
    new_unit_cell = self.box.xray_structure_box.unit_cell()
    self.geo_box = grm.geometry.select(self.box.selection_within
      ).discard_symmetry(new_unit_cell=new_unit_cell)
    self.reference_manager_box = \
      self.geo_box.generic_restraints_manager.reference_manager
    self.box_restraints_manager = mmtbx.restraints.manager(
      geometry      = self.geo_box,
      normalization = True)
    self.selection_all_box = flex.bool(self.n_sites_box, True)
    self.selection_in_box = selection.select(self.selection_within).iselection()
    self.others_in_box = flex.bool(self.n_sites_box, True).set_selected(
      self.selection_in_box, False).iselection()
    assert (len(self.selection_in_box.intersection(self.others_in_box)) == 0)
    self.box_selected_hierarchy = self.box.pdb_hierarchy_box.select(
      self.selection_in_box)
    self.box_selected_hierarchy.atoms().reset_i_seq()
    sites_cart_box = self.box.xray_structure_box.sites_cart()
    self.atom_id_mapping = {}
    for atom in self.box_selected_hierarchy.atoms() :
      self.atom_id_mapping[atom.id_str()] = atom.i_seq

  def get_selected_sites (self, selection=None, hydrogens=True) :
    from scitbx.array_family import flex
    if (selection is None) :
      selection = self.selection_in_box
    if (type(selection).__name__ == 'size_t') :
      selection = flex.bool(self.n_sites_box, False).set_selected(selection, True)
    if (not hydrogens) :
      selection &= ~self.hd_selection_box
    return self.box.xray_structure_box.sites_cart().select(selection)

  def only_residue (self) :
    """
    Returns the atom_group object corresponding to the original selection.
    Will raise an exception if the selection covers more than one atom_group.
    """
    residue = None
    box_atoms = self.box.pdb_hierarchy_box.atoms()
    sel_atoms = box_atoms.select(self.selection_in_box)
    for atom in sel_atoms :
      parent = atom.parent()
      if (residue is None) :
        residue = parent
      else :
        assert (parent == residue)
    return residue

  def mean_density_at_sites (self, selection=None) :
    """
    Compute the mean density of the target map at coordinates for non-hydrogen
    atoms.
    """
    if (selection is None) :
      selection = self.selection_in_box
    unit_cell = self.box.xray_structure_box.unit_cell()
    sum_density = n_heavy = 0
    for site in self.get_selected_sites(selection=selection, hydrogens=False) :
      site_frac = unit_cell.fractionalize(site)
      sum_density += self.target_map_box.tricubic_interpolation(site_frac)
      n_heavy += 1
    assert (n_heavy > 0)
    return sum_density / n_heavy

  def cc_model_map (self, selection=None, radius=1.5) :
    """
    Calculate the correlation coefficient for the current model (in terms of
    F(calc) from the xray structure) and the target map, calculated at atomic
    positions rather than grid points.  This will be much
    less accurate than the CC calculated in the original crystal environment,
    with full F(model) including bulk solvent correction.
    """
    from scitbx.array_family import flex
    if (selection is None) :
      selection = self.selection_in_box
    fcalc = self.box.xray_structure_box.structure_factors(d_min=self.d_min).f_calc()
    fc_fft_map = fcalc.fft_map(resolution_factor=self.resolution_factor)
    fc_map = fc_fft_map.apply_sigma_scaling().real_map_unpadded()
    sites_selected = self.get_selected_sites(selection, hydrogens=False)
    assert (len(sites_selected) > 0)
    fc_values = flex.double()
    map_values = flex.double()
    unit_cell = self.box.xray_structure_box.unit_cell()
    for site in sites_selected :
      site_frac = unit_cell.fractionalize(site)
      fc_values.append(fc_map.tricubic_interpolation(site_frac))
      map_values.append(self.target_map_box.tricubic_interpolation(site_frac))
    return flex.linear_correlation(
      x=map_values,
      y=fc_values).coefficient()

  def restrain_atoms (self,
        selection,
        reference_sigma,
        limit=1.0,
        top_out_potential=False) :
    """
    Apply harmonic reference restraints to the selected atoms, wiping out
    any previously existing harmonic restraints.
    """
    assert (reference_sigma > 0)
    sites_cart_box = self.box.xray_structure_box.sites_cart()
    reference_sites = sites_cart_box.select(selection)
    self.reference_manager_box.remove_coordinate_restraints()
    self.reference_manager_box.add_coordinate_restraints(
        sites_cart = reference_sites,
        selection  = selection,
        sigma      = reference_sigma,
        limit      = limit,
        top_out_potential=top_out_potential)

  def real_space_refine (self, selection=None) :
    """
    Run simple real-space refinement, returning the optimized weight for the
    map target (wx).
    """
    import mmtbx.refinement.real_space.individual_sites
    if (selection is None) :
      selection = self.selection_all_box
    rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
      target_map                  = self.target_map_rsr_box,
      selection                   = selection,
      real_space_gradients_delta  = self.d_min*self.resolution_factor,
      max_iterations              = 150,
      geometry_restraints_manager = self.geo_box)
    refined = mmtbx.refinement.real_space.individual_sites.refinery(
      refiner                  = rsr_simple_refiner,
      xray_structure           = self.box.xray_structure_box,
      start_trial_weight_value = 1.0,
      rms_bonds_limit          = 0.01,
      rms_angles_limit         = 1.0)
    self.update_coordinates(refined.sites_cart_result)
    print >> self.out, \
      "  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
        (refined.rms_bonds_final, refined.rms_angles_final)
    return refined.weight_final

  def fit_residue_in_box (self, mon_lib_srv=None, rotamer_manager=None) :
    import mmtbx.refinement.real_space.fit_residue
    if (mon_lib_srv is None) :
      import mmtbx.monomer_library.server
      mon_lib_srv = mmtbx.monomer_library.server.server()
    if (rotamer_manager is None) :
      from mmtbx.rotamer import rotamer_eval
      rotamer_manager = rotamer_eval.RotamerEval()
    residue = self.only_residue()
    self.restrain_atoms(
      selection=self.others_in_box,
      reference_sigma=0.1)
    self.real_space_refine()
    mmtbx.refinement.real_space.fit_residue.manager(
      target_map           = self.target_map_box,
      mon_lib_srv          = mon_lib_srv,
      special_position_settings = \
        self.box.xray_structure_box.special_position_settings(),
      residue              = residue,
      sites_cart_all       = None,
      rotamer_manager      = rotamer_manager,
      debug                = False,# XXX
      torsion_search_all_start = 0,
      torsion_search_all_stop  = 360,
      torsion_search_all_step  = 5)

  def update_coordinates (self, sites_cart) :
    """
    Set coordinates for the boxed xray_structure and pdb_hierarchy objects.
    """
    self.box.xray_structure_box.set_sites_cart(sites_cart)
    self.box.pdb_hierarchy_box.atoms().set_xyz(sites_cart)

  def update_original_coordinates (self) :
    """
    Copies the new coordinates from the box to the selected atoms in the
    original xray_structure object, and returns the (complete) list of
    coordinates.
    """
    sites_cart_box_refined_shifted_back = \
      self.box.xray_structure_box.sites_cart() + \
      self.box.shift_to_map_boxed_sites_back
    sites_cart_refined = sites_cart_box_refined_shifted_back.select(
      self.selection_in_box)
    assert (len(self.iselection) == len(sites_cart_refined))
    sites_cart_moving = self.xray_structure.sites_cart().set_selected(
       self.iselection, sites_cart_refined)
    return sites_cart_moving

  def geometry_minimization (self,
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
    from mmtbx.refinement import geometry_minimization
    from cctbx import geometry_restraints
    from scitbx.array_family import flex
    import scitbx.lbfgs
    if (selection is None) :
      selection = self.selection_in_box
    if (type(selection).__name__ == 'size_t') :
      selection = flex.bool(self.n_sites_box, False).set_selected(selection, True)
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
    site_labels = self.box.xray_structure_box.scatterers().extract_labels()
    for i in xrange(number_of_macro_cycles):
      sites_cart = self.box.xray_structure_box.sites_cart()
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True)
      minimized = geometry_minimization.lbfgs(
        sites_cart                  = sites_cart,
        geometry_restraints_manager = self.geo_box,
        geometry_restraints_flags   = geometry_restraints_flags,
        lbfgs_termination_params    = lbfgs_termination_params,
        lbfgs_exception_handling_params = exception_handling_params,
        sites_cart_selection        = selection,
        rmsd_bonds_termination_cutoff = rmsd_bonds_termination_cutoff,
        rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff,
        site_labels=site_labels)
      self.update_coordinates(sites_cart)

  def anneal (self,
              simulated_annealing_params=None,
              start_temperature=None,
              cool_rate=None,
              number_of_steps=50) :
    """
    Run real-space simulated annealing using the target map (not the RSR map,
    if this is different).  In practice, the non-selection atoms in the box
    should almost always be restrained to their current positions, but the
    setup is left to the calling code.
    """
    from mmtbx.dynamics import simulated_annealing
    import mmtbx.utils
    wx = self.real_space_refine(selection=self.selection_all_box)
    if (self.debug) :
      self.box.write_pdb_file("box_start.pdb")
    states_collector = None
    if (self.debug) :
      states_collector = mmtbx.utils.states(
        xray_structure=self.box.xray_structure_box,
        pdb_hierarchy=self.box.pdb_hierarchy_box)
    if (simulated_annealing_params is None) :
      simulated_annealing_params = simulated_annealing.master_params().extract()
    if (start_temperature is not None) :
      simulated_annealing_params.start_temperature = start_temperature
    if (cool_rate is not None) :
      simulated_annealing_params.cool_rate = cool_rate
    if (number_of_steps is not None) :
      simulated_annealing_params.number_of_steps = number_of_steps
    simulated_annealing.run(
      params             = simulated_annealing_params,
      fmodel             = None,
      xray_structure     = self.box.xray_structure_box,
      real_space         = True,
      target_map         = self.target_map_box,
      restraints_manager = self.box_restraints_manager,
      wx                 = wx,
      wc                 = 1.0,
      log                = self.out,
      verbose            = True,
      states_collector   = states_collector)
    if (states_collector is not None) :
      states_collector.write("box_traj.pdb")
    self.update_coordinates(self.box.xray_structure_box.sites_cart())

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
  box = box_build_refine_base(
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
    box.restrain_atoms(
      selection=box.others_in_box,
      reference_sigma=reference_sigma)
  box.anneal(simulated_annealing_params=params)
  if (rsr_after_anneal) :
    box.real_space_refine(selection=box.selection_in_box)
  return box.update_original_coordinates()

#-----------------------------------------------------------------------
# SAMPLING UTILITIES
def generate_sidechain_clusters (residue, mon_lib_srv) :
  """
  Extract Chi angle indices (including rotation axis) from the atom_group
  """
  from mmtbx.refinement.real_space import fit_residue
  from mmtbx.utils import rotatable_bonds
  from scitbx.array_family import flex
  atoms = residue.atoms()
  axes_and_atoms_aa_specific = \
      rotatable_bonds.axes_and_atoms_aa_specific(residue = residue,
        mon_lib_srv = mon_lib_srv)
  result = []
  if(axes_and_atoms_aa_specific is not None):
    for i_aa, aa in enumerate(axes_and_atoms_aa_specific):
      n_heavy = 0
      for i_seq in aa[1] :
        if (atoms[i_seq].element.strip() != "H") :
          n_heavy += 1
      if (n_heavy == 0) : continue
      if(i_aa == len(axes_and_atoms_aa_specific)-1):
        result.append(fit_residue.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t(aa[1]),
          vector=None)) # XXX
      else:
        result.append(fit_residue.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t([aa[1][0]]),
          vector=None)) # XXX
  return result
