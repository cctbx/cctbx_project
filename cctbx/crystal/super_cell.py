from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
import cctbx.crystal
import iotbx.pdb.utils
import boost_adaptbx.boost.python as bp
import libtbx
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

# A site whose symmetry mate is within this distance is taken to lie on a
# symmetry element and is moved onto it.  Coordinates written to three
# decimals leave such a site within 0.002 A of its mate.
min_distance_sym_equiv_coincident = 1.e-2

class sym_equiv_tables(libtbx.slots_getstate_setstate):
  """The tables get_siiu reads.

  Both must be built from the same sites; neither records which, so a
  mismatched pair gives wrong equivalents.
  """

  __slots__ = ["pair_sym_table", "site_symmetry_table"]

  def __init__(self, pair_sym_table, site_symmetry_table):
    self.pair_sym_table = pair_sym_table
    self.site_symmetry_table = site_symmetry_table

def get_sym_equiv_tables(sites_cart, crystal_symmetry, radius,
                         min_distance_sym_equiv=None):
  """Pair and site symmetry tables covering every contact within radius.

  Neither depends on the selection, so one build serves many of them.

  Parameters
  ----------
  sites_cart : flex.vec3_double
  crystal_symmetry : cctbx.crystal.symmetry
  radius : float
      In Angstrom.
  min_distance_sym_equiv : float, optional
      min_distance_sym_equiv_coincident by default.

  Returns
  -------
  sym_equiv_tables
  """
  # An empty crystal_symmetry is truthy, so test the members.
  assert crystal_symmetry is not None, "no crystal_symmetry"
  assert crystal_symmetry.unit_cell() is not None, "no unit cell"
  assert crystal_symmetry.space_group() is not None, "no space group"
  assert 0 <= radius < float("inf"), "radius is %s" % radius
  if (min_distance_sym_equiv is None):
    min_distance_sym_equiv = min_distance_sym_equiv_coincident
  sps = crystal_symmetry.special_position_settings(
    min_distance_sym_equiv=min_distance_sym_equiv)
  site_symmetry_table = sps.site_symmetry_table(sites_cart=sites_cart)
  # add_all_pairs searches to distance_cutoff*(1+is_inside_epsilon()), which
  # defaults to 1e-6; the buffer matches.
  pair_asu_table = sps.pair_asu_table(
    distance_cutoff               = radius,
    sites_cart                    = sites_cart,
    site_symmetry_table           = site_symmetry_table,
    asu_mappings_buffer_thickness = radius * (1 + 1.e-6))
  # all_interactions_from_inside_asu keeps every operator of a group of
  # symmetry-equivalent interactions.
  return sym_equiv_tables(
    pair_sym_table = pair_asu_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq       = False,
      all_interactions_from_inside_asu = True),
    site_symmetry_table = site_symmetry_table)

def get_siiu(pdb_hierarchy=None, crystal_symmetry=None,
             select_within_radius=None, sites_cart=None, selection=None,
             buffer=1, min_distance_sym_equiv=0.5, symmetry_tables=None):
  """Symmetry equivalents within select_within_radius+buffer of the selection.

  Parameters
  ----------
  pdb_hierarchy : iotbx.pdb.hierarchy.root, optional
      Supplies sites_cart when that is not given.
  crystal_symmetry : cctbx.crystal.symmetry
  select_within_radius : float
      In Angstrom.
  sites_cart : flex.vec3_double, optional
  selection : iterable of int, optional
      The sites to search around; all of them by default.
  buffer : float
  min_distance_sym_equiv : float
      Pass min_distance_sym_equiv_coincident to keep a site written on a
      symmetry element where it is.
  symmetry_tables : sym_equiv_tables, optional
      From get_sym_equiv_tables, to search several selections without
      rebuilding them.

  Returns
  -------
  tuple
      ``({j_seq: [rt_mx_ji]}, [rt_mx_ji])``.  Applying rt_mx_ji to the
      fractional coordinates of site j_seq places that equivalent near a
      selected site, plus the distinct operators among them.  Operators
      carry denominators (1, 12), so equal operators hash equally.  One
      operator per equivalent, and an equivalent coincident with its site
      is dropped.
  """
  if (symmetry_tables is None):
    if (sites_cart is None):
      assert pdb_hierarchy is not None, "no sites_cart and no pdb_hierarchy"
      sites_cart = pdb_hierarchy.atoms().extract_xyz()
    # +1 is nonbonded buffer, to match nonbonded_distance_cutoff
    symmetry_tables = get_sym_equiv_tables(
      sites_cart             = sites_cart,
      crystal_symmetry       = crystal_symmetry,
      radius                 = select_within_radius + buffer,
      min_distance_sym_equiv = min_distance_sym_equiv)
  if (selection is None):
    selection = range(symmetry_tables.pair_sym_table.size())
  matrices_by_j_seq = {}
  normalised = {}
  by_j_seq = {}
  pair_sym_table = symmetry_tables.pair_sym_table
  site_symmetry_table = symmetry_tables.site_symmetry_table
  assert site_symmetry_table.indices().size() == pair_sym_table.size(), (
    "the two tables cover different numbers of sites")
  for i_seq in sorted(set(selection)):
    assert 0 <= i_seq < pair_sym_table.size(), "selection index is %s" % i_seq
    for j_seq, rt_mx_ji_list in pair_sym_table[i_seq].items():
      by_key = None
      # stl_vector_rt_mx has no __iter__, so a for-in loop falls back to the
      # __getitem__ protocol and costs one C++ exception per vector.
      for i_op in range(len(rt_mx_ji_list)):
        rt_mx_ji = rt_mx_ji_list[i_op]
        # Site j's own symmetry operators are recorded as the identity, so
        # this drops them as well.
        if (rt_mx_ji.is_unit_mx()): continue
        if (by_key is None):
          # matrices() builds a fresh tuple on each call, so cache it; None
          # records a general position.
          if (j_seq in matrices_by_j_seq):
            site_symmetry_matrices = matrices_by_j_seq[j_seq]
          elif (site_symmetry_table.is_special_position(j_seq)):
            site_symmetry_matrices = site_symmetry_table.get(j_seq).matrices()
            matrices_by_j_seq[j_seq] = site_symmetry_matrices
          else:
            site_symmetry_matrices = matrices_by_j_seq[j_seq] = None
          by_key = by_j_seq.get(j_seq)
          if (by_key is None):
            by_key = by_j_seq[j_seq] = {}
        if (site_symmetry_matrices is not None):
          # The canonical coset member, so a union of selections gives the
          # union of the results.  min() orders by rt_mx.__lt__, the
          # ordering of space_group.make_tidy().
          rt_mx_ji = min([rt_mx_ji.multiply(matrix)
                          for matrix in site_symmetry_matrices])
        # Equal operators hash equally only on equal denominators, so the
        # normalised operator is the key.  A miss on an equal operator
        # normalises it twice, which is harmless.
        key = normalised.get(rt_mx_ji)
        if (key is None):
          key = normalised[rt_mx_ji] = rt_mx_ji.new_denominators(1, 12)
        by_key[key] = None
  siiu = {}
  for j_seq, by_key in by_j_seq.items():
    siiu[j_seq] = list(by_key)
  # Every operator that reached a by_key is in normalised, on denominators
  # (1, 12) and so hashing consistently.
  return siiu, list(dict.fromkeys(normalised.values()))

def sym_equiv_sites_cart(sites_cart, unit_cell, rt_mx, selection=None):
  """The symmetry-equivalent sites rt_mx generates.

  Parameters
  ----------
  sites_cart : flex.vec3_double
  unit_cell : cctbx.uctbx.unit_cell
  rt_mx : cctbx.sgtbx.rt_mx
  selection : flex.size_t, optional
      The sites to transform; all of them by default.

  Returns
  -------
  flex.vec3_double
  """
  if (selection is not None):
    sites_cart = sites_cart.select(selection)
  return ((unit_cell.matrix_cart(rt_mx.r()) * sites_cart)
          + unit_cell.orthogonalize(rt_mx.t().as_double()))

def apply_symop_inplace_chain(chain, op, unit_cell):
  chain.atoms().set_xyz(sym_equiv_sites_cart(
    sites_cart=chain.atoms().extract_xyz(), unit_cell=unit_cell, rt_mx=op))

class manager(object):
  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               select_within_radius,
               box_buffer_layer=3,
               siiu=None,
               debug_files=False,
               box_super_sphere = False):
    # Results to be available outside
    self.super_cell_hierarchy = None # XXX sites are NOT in box center!!!
    self.super_sphere_hierarchy = None # XXX sites are NOT in box center!!!
    self.cs_super_sphere      = None
    self.siiu                 = siiu
    self.pdb_hierarchy        = pdb_hierarchy
    self.chain_op_dict        = {}
    self.unit_cell = crystal_symmetry.unit_cell()
    #
    assert pdb_hierarchy.models_size() == 1, "one model is expected"
    # Get operators
    self.siiu, all_ops = get_siiu(
      pdb_hierarchy        = pdb_hierarchy,
      crystal_symmetry     = crystal_symmetry,
      select_within_radius = select_within_radius)
    # Get unique chain IDs for new symmetry copies
    all_chains = iotbx.pdb.utils.all_chain_ids()
    focus_chain_ids = pdb_hierarchy.chain_ids()
    for cid in focus_chain_ids:
      while cid in all_chains:
        all_chains.remove(cid)
    focus_ss = " or ".join(["chain %s"%it for it in focus_chain_ids])
    # Create super_cell: apply symmetry to chains and add new copies to master
    # hierarchy
    self.super_cell_hierarchy = pdb_hierarchy.deep_copy()
    cntr=0
    new_chains = []
    for op in all_ops:
      model_dc = self.super_cell_hierarchy.models()[0].detached_copy()
      new_ids = []
      for chain in model_dc.chains():
        chain_dc = chain.detached_copy()
        new_id = all_chains[cntr]
        cntr+=1
        chain_dc.id = new_id
        new_ids.append(new_id)
        apply_symop_inplace_chain(
          chain=chain_dc, op=op, unit_cell=self.unit_cell)
        new_chains.append(chain_dc)
      self.chain_op_dict[op] = new_ids
    for c in new_chains:
      self.super_cell_hierarchy.models()[0].append_chain(c)
    # Initiate selection to remove (outside the sphere)
    self.super_cell_hierarchy.atoms().reset_i_seq()
    self.keep_selection = flex.bool(self.super_cell_hierarchy.atoms().size(), True)
    # Move everything to the center of a big box
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart   = self.super_cell_hierarchy.atoms().extract_xyz(),
      buffer_layer = box_buffer_layer)
    self.super_cell_hierarchy.atoms().set_xyz(box.sites_cart)
    if(debug_files):
      self.super_cell_hierarchy.write_pdb_file(file_name="ss.pdb",
        crystal_symmetry = box.crystal_symmetry())
    # Create mask around focus atoms
    radii = flex.double(self.super_cell_hierarchy.atoms().size(), select_within_radius)
    uc = box.crystal_symmetry().unit_cell()
    fm_ss = box.crystal_symmetry().unit_cell().fractionalization_matrix()
    a,b,c = uc.parameters()[:3]
    step = 1.
    n_real = [int(a/step), int(b/step), int(c/step)]

    self.scasc = self.super_cell_hierarchy.atom_selection_cache()

    sel_focus = self.super_cell_hierarchy.atom_selection_cache().\
      selection(focus_ss)
    sites_frac_focus = fm_ss * self.super_cell_hierarchy.select(
      sel_focus).atoms().extract_xyz()

    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = sites_frac_focus,
      unit_cell                   = uc,
      n_real                      = n_real,
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = radii,
      wrapping                    = False)
    # Select the inside of the mask (discard the rest). This is super-sphere.
    # The mask is read for every atom at once, rather than one call per atom.
    nx, ny, nz = n_real
    fx, fy, fz = (fm_ss*self.super_cell_hierarchy.atoms().extract_xyz()).parts()
    grid_index = ((((fx*nx).iround() % nx) * ny
                   + ((fy*ny).iround() % ny)) * nz
                  + ((fz*nz).iround() % nz))
    inside_atom = mask.as_1d().select(grid_index.as_size_t()) == 1
    for chain in self.super_cell_hierarchy.chains():
      for rg in chain.residue_groups():
        i_seqs = rg.atoms().extract_i_seq()
        if(not inside_atom.select(i_seqs).count(True)):
          #chain.remove_residue_group(rg)
          self.keep_selection.set_selected(i_seqs, False)
    if(debug_files):
      self.super_cell_hierarchy.write_pdb_file(file_name="ss_cut.pdb",
        crystal_symmetry = box.crystal_symmetry())
    # Shift back
    self.super_cell_hierarchy.atoms().set_xyz(
      self.super_cell_hierarchy.atoms().extract_xyz()-box.shift_vector)
    # Extract super-sphere out of super-cell
    self.super_sphere_hierarchy = self.super_cell_hierarchy.select(
      self.keep_selection, copy_atoms=True)
    self.super_sphere_hierarchy.atoms().reset_i_seq()
    #
    if(debug_files):
      self.super_sphere_hierarchy.write_pdb_file(file_name="ss_cut_shifted.pdb")
    # Box final model
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart   = self.super_sphere_hierarchy.atoms().extract_xyz(),
      buffer_layer = box_buffer_layer)
    self.cs_super_sphere = box.crystal_symmetry()
    self.super_sphere_hierarchy.atoms().reset_i_seq()
    if box_super_sphere:
      self.super_sphere_hierarchy.atoms().set_xyz(box.sites_cart)

    #
    # DEBUG ONLY
    # To trigger the bug: change copy_atoms to False above, uncomment the code
    # below and run modules/cctbx_project/cctbx/crystal/tst_super_cell.py
    #
    #import mmtbx.model
    #from libtbx.utils import null_out
    #model = mmtbx.model.manager(
    #  model_input       = None,
    #  pdb_hierarchy     = self.super_cell_hierarchy,
    #  crystal_symmetry  = self.cs_super_sphere,
    #  log               = null_out())
    #model.process(make_restraints=True, grm_normalization=True)
    #STOP()

  def update(self, sites_cart, debug=False):
    """
    Change coordinates in the master (focus) hierarchy and propagate the change
    to the super-cell/sphere (self.super_cell_hierarchy,
    self.super_sphere_hierarchy). This changes self.super_cell_hierarchy and
    self.super_sphere_hierarchy inplace
    """
    assert self.pdb_hierarchy.atoms().size() == sites_cart.size()
    if(debug):
      dist1 = flex.sqrt((
         self.pdb_hierarchy.atoms().extract_xyz() -  sites_cart).dot())
      BEFORE = self.super_cell_hierarchy.atoms().extract_xyz()
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    all_xyz = self.super_cell_hierarchy.atoms().extract_xyz()
    all_xyz.set_selected(flex.size_t(range(sites_cart.size())), sites_cart)
    for op, cids in zip(self.chain_op_dict.keys(), self.chain_op_dict.values()):
      sel_str = " or ".join(["chain %s"%it for it in cids])
      sel = self.scasc.selection(sel_str)
      sites_cart_copy = sym_equiv_sites_cart(
        sites_cart=sites_cart, unit_cell=self.unit_cell, rt_mx=op)
      all_xyz.set_selected(sel, sites_cart_copy)
    self.super_cell_hierarchy.atoms().set_xyz(all_xyz)
    self.super_sphere_hierarchy = self.super_cell_hierarchy.select(
      self.keep_selection, copy_atoms=True)
    self.super_sphere_hierarchy.atoms().reset_i_seq()
    if(debug):
      AFTER = self.super_cell_hierarchy.atoms().extract_xyz()
      dist2 = flex.sqrt((BEFORE - AFTER).dot())
      assert abs(flex.max(dist1)-flex.max(dist2))<1.e-4
    return self.super_sphere_hierarchy
