
from __future__ import absolute_import, division, print_function
from libtbx import group_args, adopt_init_args
import sys
import mmtbx.refinement.real_space
from six.moves import range

#-----------------------------------------------------------------------
# MODEL UTILITIES
def show_chain_resseq_ranges(residues, out=sys.stdout, prefix=""):
  """
  Given a list of residues (either residue_group or atom_group objects),
  print a summary of the chains and residue ranges present.
  """
  ranges = {}
  prev_chain = prev_resseq = None
  for i_res, residue in enumerate(residues):
    rg = residue
    if (type(residue).__name__ == 'atom_group'):
      rg = residue.parent()
    chain = rg.parent()
    if (not chain.id in ranges):
      ranges[chain.id] = []
      prev_resseq = None
    resseq = rg.resseq_as_int()
    if (prev_resseq is None) or (resseq > prev_resseq + 1):
      ranges[chain.id].append([ (resseq, rg.icode), ])
    else :
      ranges[chain.id][-1].append( (resseq, rg.icode) )
    prev_resseq = resseq
  def format_resid(resseq, icode):
    return ("%s%s" % (resseq, icode)).strip()
  for chain_id in sorted(ranges.keys()):
    clauses = []
    for range in ranges[chain_id] :
      if (len(range) >= 2):
        clauses.append("%s-%s" % (format_resid(*(range[0])),
          format_resid(*(range[-1]))))
      else :
        clauses.append(format_resid(*(range[0])))
    print(prefix+"chain '%s': %s" % (chain_id, ",".join(clauses)), file=out)

def reprocess_pdb(pdb_hierarchy, crystal_symmetry, cif_objects, out):
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

def get_nearby_water_selection(
    pdb_hierarchy,
    xray_structure,
    selection,
    selection_buffer_radius=5):
  sel_cache = pdb_hierarchy.atom_selection_cache()
  water_selection = sel_cache.selection("resname HOH")
  selection_within = xray_structure.selection_within(
    radius=selection_buffer_radius,
    selection=selection)
  return (selection_within & water_selection)

def iter_residue_groups(hierarchy):
  for chain in hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      yield residue_group

def extract_iselection(pdb_objects):
  from scitbx.array_family import flex
  isel = flex.size_t()
  for pdb_obj in pdb_objects :
    isel.extend(pdb_obj.atoms().extract_i_seq())
  assert not isel.all_eq(0)
  return isel

def remove_sidechain_atoms(pdb_objects):
  for pdb_obj in pdb_objects :
    atoms = pdb_obj.atoms()
    for atom in atoms :
      if (not atom.name.strip() in ["CA","C","N","O","CB"]):
        atom.parent().remove_atom(atom)

def is_stub_residue(atom_group):
  if (atom_group.resname in ["GLY","ALA"]):
    return True
  has_sidechain_atoms = False
  for atom in atom_group.atoms():
    if (not atom.name.strip() in ["H","HA","C","CA","CB","N","O",]):
      has_sidechain_atoms = True
  return not has_sidechain_atoms

def get_non_hydrogen_atom_indices(pdb_object):
  from scitbx.array_family import flex
  i_seqs_no_hd = flex.size_t()
  for atom in pdb_object.atoms():
    if (atom.element.strip() not in ["H","D"]):
      i_seqs_no_hd.append(atom.i_seq)
  return i_seqs_no_hd

def get_window_around_residue(residue, window_size=2):
  residue_group = residue.parent()
  chain = residue_group.parent()
  all_residues = chain.residue_groups()
  sel_residues = []
  for i_res, other_rg in enumerate(all_residues):
    if (other_rg == residue_group):
      for j_res in range(-window_size, window_size+1):
        k_res = i_res + j_res
        if (k_res >= 0) and (k_res < len(all_residues)):
          sel_residues.append(all_residues[k_res])
  assert (len(sel_residues) > 0)
  return sel_residues

def atom_group_as_hierarchy(atom_group):
  """
  Convert an iotbx.pdb.hierarchy.atom_group object to a separate hierarchy
  (leaving the original hierarchy untouched), which allows us to use the
  pickling capabilities of iotbx.pdb.hierarchy.root.
  """
  import iotbx.pdb.hierarchy
  old_rg = atom_group.parent()
  assert (old_rg is not None)
  old_chain = old_rg.parent()
  root = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  chain = iotbx.pdb.hierarchy.chain(id=old_chain.id)
  residue_group = iotbx.pdb.hierarchy.residue_group(
    resseq=old_rg.resseq,
    icode=old_rg.icode)
  residue_group.append_atom_group(atom_group.detached_copy())
  chain.append_residue_group(residue_group)
  model.append_chain(chain)
  root.append_model(model)
  atoms = root.atoms()
  atoms.reset_i_seq()
  atoms.reset_serial()
  return root

def residues_are_adjacent(residue1, residue2, max_sep=2.5):
  if (type(residue1).__name__ == "residue_group"):
    residue1 = residue1.atom_groups()[0]
  if (type(residue2).__name__ == "residue_group"):
    residue2 = residue2.atom_groups()[0]
  c_pep = n_pep = None
  for atom in residue1.atoms():
    if (atom.name.strip() == "C"):
      c_pep = atom
      break
  for atom in residue2.atoms():
    if (atom.name.strip() == "N"):
      n_pep = atom
      break
  if (None in [c_pep, n_pep]):
    return False
  return c_pep.distance(n_pep) < max_sep

#-----------------------------------------------------------------------
# MAP-RELATED
class local_density_quality(object):
  def __init__(self,
      fofc_map,
      two_fofc_map,
      sites_cart=None,
      atom_names=None,
      unit_cell=None,
      atom_selection=None,
      xray_structure=None,
      radius=2.5):
    from cctbx import maptbx
    from scitbx.array_family import flex
    assert ((len(fofc_map) == len(two_fofc_map)) and
            (fofc_map.focus() == two_fofc_map.focus()))
    self.atom_selection = atom_selection # assumes no HD already
    self.sites_cart = sites_cart
    if (atom_names is not None) and (type(atom_names).__name__!='std_string'):
      atom_names = flex.std_string(atom_names)
    self.atom_names = atom_names
    self.unit_cell = unit_cell
    if (atom_selection is not None):
      assert (xray_structure is not None) and (len(self.atom_selection) > 0)
      self.sites_cart = xray_structure.sites_cart().select(self.atom_selection)
      self.unit_cell = xray_structure.unit_cell()
      labels = xray_structure.scatterers().extract_labels()
      self.atom_names = labels.select(self.atom_selection)
    self.map_sel = maptbx.grid_indices_around_sites(
      unit_cell=self.unit_cell,
      fft_n_real=fofc_map.focus(),
      fft_m_real=fofc_map.all(),
      sites_cart=self.sites_cart,
      site_radii=flex.double(self.sites_cart.size(), radius))
    assert (len(self.map_sel) > 0)
    self.fofc_map_sel = fofc_map.select(self.map_sel)
    self.two_fofc_map_sel = two_fofc_map.select(self.map_sel)
    self.fofc_map_levels = flex.double()
    self.two_fofc_map_levels = flex.double()
    for site_cart in self.sites_cart :
      site_frac = self.unit_cell.fractionalize(site_cart=site_cart)
      fofc_map_value = fofc_map.tricubic_interpolation(site_frac)
      two_fofc_map_value = two_fofc_map.tricubic_interpolation(site_frac)
      self.fofc_map_levels.append(fofc_map_value)
      self.two_fofc_map_levels.append(two_fofc_map_value)
    self._atom_lookup = { n.strip():i for i,n in enumerate(self.atom_names) }

  def density_at_atom(self, atom_name):
    """Return a simple object containing density levels for the named atom,
    or None if no match was found."""
    i_atom = self._atom_lookup.get(atom_name.strip())
    if (i_atom is not None):
      return group_args(
        two_fofc=self.two_fofc_map_levels[i_atom],
        fofc=self.fofc_map_levels[i_atom])
    return None

  def atoms_below_fofc_map_level(self, sigma_cutoff=-2.5):
    selection = self.fofc_map_levels <= sigma_cutoff
    return self.atom_names.select(selection)

  def atoms_below_two_fofc_map_level(self, sigma_cutoff=-2.5):
    selection = self.two_fofc_map_levels <= sigma_cutoff
    return self.atom_names.select(selection)

  def _atoms_outside_density_selection(self, fofc_cutoff=3.0,
      two_fofc_cutoff=1.0,
      require_difference_map_peaks=False):
    sel_outside_fofc = (self.fofc_map_levels <= fofc_cutoff)
    if (not require_difference_map_peaks):
      sel_outside_two_fofc = (self.two_fofc_map_levels < two_fofc_cutoff)
    else :
      sel_outside_two_fofc = self.two_fofc_map_levels < float(sys.maxsize)
    return (sel_outside_fofc & sel_outside_two_fofc)

  def atoms_outside_density(self, **kwds):
    selection = self._atoms_outside_density_selection(**kwds)
    return self.atom_names.select(selection)

  def fraction_of_nearby_grid_points_above_cutoff(self, sigma_cutoff=2.5):
    return (self.fofc_map_sel >= sigma_cutoff).count(True) / len(self.fofc_map_sel)

  def number_of_atoms_below_fofc_map_level(self, *args, **kwds):
    return len(self.atoms_below_fofc_map_level(*args, **kwds))

  def number_of_atoms_below_two_fofc_map_level(self, *args, **kwds):
    return len(self.atoms_below_two_fofc_map_level(*args, **kwds))

  def number_of_atoms_outside_density(self, *args, **kwds):
    return len(self.atoms_outside_density(*args, **kwds))

  def show_atoms_outside_density(self, **kwds):
    out = kwds.pop("out", sys.stdout)
    prefix = kwds.pop("prefix", "")
    isel = self._atoms_outside_density_selection(**kwds).iselection()
    for i_seq in isel :
      print(prefix+"%s: 2mFo-DFc=%6.2f  mFo-DFc=%6.2f" % \
        (self.atom_names[i_seq], self.two_fofc_map_levels[i_seq],
         self.fofc_map_levels[i_seq]), file=out)
    return len(isel)

  def max_fofc_value(self):
    from scitbx.array_family import flex
    return flex.max(self.fofc_map_sel)

def residue_density_quality(
    atom_group,
    unit_cell,
    two_fofc_map,
    fofc_map):
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double()
  atom_names = flex.std_string()
  for atom in atom_group.atoms():
    if (not atom.element.strip() in ["H","D"]):
      sites_cart.append(atom.xyz)
      atom_names.append(atom.name)
  if (len(sites_cart) == 0):
    raise RuntimeError("No non-hydrogen atoms in %s" % atom_group.id_str())
  density = local_density_quality(
    fofc_map=fofc_map,
    two_fofc_map=two_fofc_map,
    sites_cart=sites_cart,
    atom_names=atom_names,
    unit_cell=unit_cell)
  return density

def get_difference_maps(fmodel, resolution_factor=0.25):
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
def get_model_map_stats(
    selection,
    target_map,
    model_map,
    unit_cell,
    sites_cart,
    pdb_atoms,
    local_sampling=False):
  """
  Collect basic statistics for a model map and some target map (usually an
  mFo-DFc map), including CC, mean, and minimum density at the atomic
  positions.
  """
  assert (len(target_map) == len(model_map))
  iselection = selection
  if (type(selection).__name__ == 'bool'):
    iselection = selection.iselection()
  from scitbx.array_family import flex
  sites_cart_refined = sites_cart.select(selection)
  sites_selected = flex.vec3_double()
  map1 = flex.double()
  map2 = flex.double()
  min_density = sys.maxsize
  sum_density = n_sites = 0
  worst_atom = None
  # XXX I'm not sure the strict density cutoff is a good idea here
  for i_seq, xyz in zip(iselection, sites_cart_refined):
    if (pdb_atoms[i_seq].element.strip() != "H"):
      sites_selected.append(xyz)
      site_frac = unit_cell.fractionalize(site_cart=xyz)
      target_value = target_map.tricubic_interpolation(site_frac)
      if (target_value < min_density):
        min_density = target_value
        worst_atom = pdb_atoms[i_seq]
      sum_density += target_value
      n_sites += 1
      if (not local_sampling):
        map1.append(target_value)
        map2.append(model_map.tricubic_interpolation(site_frac))
  assert (n_sites > 0)
  if (local_sampling):
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

class box_build_refine_base(object):
  """
  Base class for handling functions associated with rebuilding and refinement
  of atoms in a box, with the building methods implemented elsewhere.
  """
  def __init__(self,
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
      geometry_restraints_manager=None,
      debug=False):
    adopt_init_args(self, locals())
    import mmtbx.restraints
    import mmtbx.utils
    from scitbx.array_family import flex
    self.iselection = selection.iselection()
    if (geometry_restraints_manager is None):
      if (self.processed_pdb_file is None):
        self.processed_pdb_file = reprocess_pdb(
          pdb_hierarchy=pdb_hierarchy,
          cif_objects=cif_objects,
          crystal_symmetry=xray_structure,
          out=out)
      self.geometry_restraints_manager = \
        self.processed_pdb_file.geometry_restraints_manager(show_energies=False)
    grm = mmtbx.restraints.manager(
      geometry=self.geometry_restraints_manager,
      normalization = True)
    self.selection_within = xray_structure.selection_within(
      radius    = selection_buffer_radius,
      selection = selection)
    self.box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure   = xray_structure,
      map_data         = target_map,
      selection        = self.selection_within,
      box_cushion      = box_cushion)
    self.hd_selection_box = self.box.xray_structure_box.hd_selection()
    self.n_sites_box = self.box.xray_structure_box.sites_cart().size()
    self.target_map_box = self.box.map_box
    self.target_map_rsr_box = None
    if (self.target_map_rsr_box is None):
      self.target_map_rsr_box = self.target_map_box
    self.unit_cell_box = self.box.xray_structure_box.unit_cell()
    self.geo_box = grm.geometry.select(self.box.selection
      ).discard_symmetry(new_unit_cell=self.unit_cell_box)
    self.box_restraints_manager = mmtbx.restraints.manager(
      geometry      = self.geo_box,
      normalization = True)
    self.selection_all_box = flex.bool(self.n_sites_box, True)
    self.selection_in_box = selection.select(self.selection_within).iselection()
    self.others_in_box = flex.bool(self.n_sites_box, True).set_selected(
      self.selection_in_box, False).iselection()
    assert (len(self.selection_in_box.intersection(self.others_in_box)) == 0)

    self.box_hierarchy = pdb_hierarchy.select(self.selection_within)
    self.box_hierarchy.adopt_xray_structure(self.box.xray_structure_box)

    self.box_selected_hierarchy = self.box_hierarchy.select(
      self.selection_in_box)
    self.box_selected_hierarchy.atoms().reset_i_seq()
    sites_cart_box = self.box.xray_structure_box.sites_cart()
    self.atom_id_mapping = {}
    for atom in self.box_selected_hierarchy.atoms():
      self.atom_id_mapping[atom.id_str()] = atom.i_seq

  def get_selected_sites(self, selection=None, hydrogens=True):
    from scitbx.array_family import flex
    if (selection is None):
      selection = self.selection_in_box
    if (type(selection).__name__ == 'size_t'):
      selection = flex.bool(self.n_sites_box, False).set_selected(selection, True)
    if (not hydrogens):
      selection &= ~self.hd_selection_box
    return self.box.xray_structure_box.sites_cart().select(selection)

  def only_residue(self):
    """
    Returns the atom_group object corresponding to the original selection.
    Will raise an exception if the selection covers more than one atom_group.
    """
    residue = None
    box_atoms = self.box_hierarchy.atoms()
    sel_atoms = box_atoms.select(self.selection_in_box)
    for atom in sel_atoms :
      parent = atom.parent()
      if (residue is None):
        residue = parent
      else :
        assert (parent == residue)
    return residue

  def update_sites_from_pdb_atoms(self):
    """
    Update the xray_structure coordinates after manipulating the PDB hierarchy.
    """
    box_atoms = self.box.pdb_hierarchy_box.atoms()
    self.box.xray_structure_box.set_sites_cart(box_atoms.extract_xyz())

  def mean_density_at_sites(self, selection=None):
    """
    Compute the mean density of the target map at coordinates for non-hydrogen
    atoms.
    """
    if (selection is None):
      selection = self.selection_in_box
    unit_cell = self.box.xray_structure_box.unit_cell()
    sum_density = n_heavy = 0
    for site in self.get_selected_sites(selection=selection, hydrogens=False):
      site_frac = unit_cell.fractionalize(site)
      sum_density += self.target_map_box.tricubic_interpolation(site_frac)
      n_heavy += 1
    assert (n_heavy > 0)
    return sum_density / n_heavy

  def cc_model_map(self, selection=None, radius=1.5):
    """
    Calculate the correlation coefficient for the current model (in terms of
    F(calc) from the xray structure) and the target map, calculated at atomic
    positions rather than grid points.  This will be much
    less accurate than the CC calculated in the original crystal environment,
    with full F(model) including bulk solvent correction.
    """
    from scitbx.array_family import flex
    if (selection is None):
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

  def restrain_atoms(self,
        selection,
        reference_sigma,
        limit=1.0,
        reset_first=True,
        top_out_potential=False):
    """
    Apply harmonic reference restraints to the selected atoms, wiping out
    any previously existing harmonic restraints.
    """
    assert (reference_sigma > 0)
    if isinstance(selection, str):
      selection = self.box_hierarchy.atom_selection_cache().selection(
        selection).iselection()
      assert (len(selection) > 0)
    sites_cart_box = self.box.xray_structure_box.sites_cart()
    reference_sites = sites_cart_box.select(selection)
    if (reset_first):
      self.geo_box.remove_reference_coordinate_restraints_in_place()
    self.geo_box.add_reference_coordinate_restraints_in_place(
        pdb_hierarchy=self.box_hierarchy,
        selection=selection,
        sigma=reference_sigma,
        limit=limit,
        top_out=top_out_potential)

  def real_space_refine(self, selection=None):
    """
    Run simple real-space refinement, returning the optimized weight for the
    map target (wx).
    """
    import mmtbx.refinement.real_space.individual_sites
    if (selection is None):
      selection = self.selection_all_box
    elif isinstance(selection, str):
      selection = self.box_hierarchy.atom_selection_cache().selection(
        selection).iselection()
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
    print("  after real-space refinement: rms_bonds=%.3f  rms_angles=%.3f" % \
        (refined.rms_bonds_final, refined.rms_angles_final), file=self.out)
    return refined.weight_final

  def fit_residue_in_box(self,
      mon_lib_srv=None,
      rotamer_manager=None,
      backbone_sample_angle=20):
    import mmtbx.refinement.real_space.fit_residue
    if (mon_lib_srv is None):
      import mmtbx.monomer_library.server
      mon_lib_srv = mmtbx.monomer_library.server.server()
    if (rotamer_manager is None):
      from mmtbx.rotamer import rotamer_eval
      rotamer_manager = rotamer_eval.RotamerEval(mon_lib_srv=mon_lib_srv)
    residue = self.only_residue()
    self.restrain_atoms(
      selection=self.others_in_box,
      reference_sigma=0.1)
    self.real_space_refine()
    mmtbx.refinement.real_space.fit_residue.run_with_minimization(
      target_map=self.target_map_box,
      residue=residue,
      xray_structure=self.box.xray_structure_box,
      mon_lib_srv=mon_lib_srv,
      rotamer_manager=rotamer_manager,
      geometry_restraints_manager=self.box_restraints_manager.geometry,
      real_space_gradients_delta=self.d_min*self.resolution_factor,
      rms_bonds_limit=0.01,
      rms_angles_limit=1.0,
      backbone_sample_angle=backbone_sample_angle)

  def update_coordinates(self, sites_cart):
    """
    Set coordinates for the boxed xray_structure and pdb_hierarchy objects.
    """
    self.box.xray_structure_box.set_sites_cart(sites_cart)
    self.box_hierarchy.atoms().set_xyz(sites_cart)

  def update_original_coordinates(self):
    """
    Copies the new coordinates from the box to the selected atoms in the
    original xray_structure object, and returns the (complete) list of
    coordinates.
    """
    sites_cart_box_refined_shifted_back = \
      self.box.xray_structure_box.sites_cart() + \
      self.box.shift_cart
    sites_cart_refined = sites_cart_box_refined_shifted_back.select(
      self.selection_in_box)
    assert (len(self.iselection) == len(sites_cart_refined))
    sites_cart_moving = self.xray_structure.sites_cart().set_selected(
       self.iselection, sites_cart_refined)
    return sites_cart_moving

  def geometry_minimization(self,
      selection=None,
      bond=True,
      nonbonded=True,
      angle=True,
      dihedral=True,
      chirality=True,
      planarity=True,
      max_number_of_iterations=500,
      number_of_macro_cycles=5,
      rmsd_bonds_termination_cutoff  = 0,
      rmsd_angles_termination_cutoff = 0):
    """
    Perform geometry minimization on a selection of boxed atoms.
    """
    from mmtbx.refinement import geometry_minimization
    from cctbx import geometry_restraints
    from scitbx.array_family import flex
    import scitbx.lbfgs
    if (selection is None):
      selection = self.selection_in_box
    if (type(selection).__name__ == 'size_t'):
      selection = flex.bool(self.n_sites_box, False).set_selected(selection, True)
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_number_of_iterations)
    geometry_restraints_flags = geometry_restraints.flags.flags(
      bond               = bond,
      nonbonded          = nonbonded,
      angle              = angle,
      dihedral           = dihedral,
      chirality          = chirality,
      planarity          = planarity)
    site_labels = self.box.xray_structure_box.scatterers().extract_labels()
    for i in range(number_of_macro_cycles):
      sites_cart = self.box.xray_structure_box.sites_cart()
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True)
      minimized = geometry_minimization.lbfgs(
        sites_cart                  = sites_cart,
        correct_special_position_tolerance = 1.0,
        geometry_restraints_manager = self.geo_box,
        geometry_restraints_flags   = geometry_restraints_flags,
        lbfgs_termination_params    = lbfgs_termination_params,
        lbfgs_exception_handling_params = exception_handling_params,
        sites_cart_selection        = selection,
        rmsd_bonds_termination_cutoff = rmsd_bonds_termination_cutoff,
        rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff,
        site_labels=site_labels)
      self.update_coordinates(sites_cart)

  def anneal(self,
              simulated_annealing_params=None,
              start_temperature=None,
              cool_rate=None,
              number_of_steps=50):
    """
    Run real-space simulated annealing using the target map (not the RSR map,
    if this is different).  In practice, the non-selection atoms in the box
    should almost always be restrained to their current positions, but the
    setup is left to the calling code.
    """
    from mmtbx.dynamics import simulated_annealing
    import mmtbx.utils
    wx = self.real_space_refine(selection=self.selection_all_box)
    if (self.debug):
      self.box.write_pdb_file("box_start.pdb")
    states_collector = None
    if (self.debug):
      states_collector = mmtbx.utils.states(
        xray_structure=self.box.xray_structure_box,
        pdb_hierarchy=self.box.pdb_hierarchy_box)
    if (simulated_annealing_params is None):
      simulated_annealing_params = simulated_annealing.master_params().extract()
    if (start_temperature is not None):
      simulated_annealing_params.start_temperature = start_temperature
    if (cool_rate is not None):
      simulated_annealing_params.cool_rate = cool_rate
    if (number_of_steps is not None):
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
    if (states_collector is not None):
      states_collector.write("box_traj.pdb")
    self.update_coordinates(self.box.xray_structure_box.sites_cart())

def run_real_space_annealing(
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
    debug=False):
  """
  Run simulated annealing with a real-space target.  For maximum flexibility,
  a separate map may be used for the initial real-space refinement step.
  """
  from mmtbx.dynamics import simulated_annealing
  import iotbx.phil
  if (params is None):
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
  if (reference_sigma is not None):
    box.restrain_atoms(
      selection=box.others_in_box,
      reference_sigma=reference_sigma)
  box.anneal(simulated_annealing_params=params)
  if (rsr_after_anneal):
    box.real_space_refine(selection=box.selection_in_box)
  return box.update_original_coordinates()

#-----------------------------------------------------------------------
# SAMPLING UTILITIES
def generate_sidechain_clusters(residue, mon_lib_srv):
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
      #for i_seq in aa[1] :
      #  if (atoms[i_seq].element.strip() != "H"):
      #    n_heavy += 1
      #if (n_heavy == 0) : continue
      if(i_aa == len(axes_and_atoms_aa_specific)-1):
        result.append(mmtbx.refinement.real_space.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t(aa[1]),
          vector=None)) # XXX
      else:
        result.append(mmtbx.refinement.real_space.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t([aa[1][0]]),
          vector=None)) # XXX
  return result
