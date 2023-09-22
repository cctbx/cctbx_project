from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
import cctbx.crystal
import iotbx.pdb.utils
import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

def get_siiu(pdb_hierarchy, crystal_symmetry, select_within_radius):
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  sst = crystal_symmetry.special_position_settings().site_symmetry_table(
    sites_cart = sites_cart)
  siiu = {}
  # +1 is nonbonded buffer, to match nonbonded_distance_cutoff
  cutoff = select_within_radius+1
  conn_asu_mappings = crystal_symmetry.special_position_settings().\
    asu_mappings(buffer_thickness=cutoff)
  conn_asu_mappings.process_sites_cart(
    original_sites      = sites_cart,
    site_symmetry_table = sst)
  conn_pair_asu_table = cctbx.crystal.pair_asu_table(
    asu_mappings=conn_asu_mappings)
  conn_pair_asu_table.add_all_pairs(
    distance_cutoff=cutoff)
  pair_generator = cctbx.crystal.neighbors_fast_pair_generator(
    conn_asu_mappings,
    distance_cutoff=cutoff)
  all_ops_mat = []
  all_ops = []
  for pair in pair_generator:
    rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
    rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
    rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
    #print(rt_mx_ji, str(rt_mx_ji))
    if str(rt_mx_ji)=="x,y,z": continue
    siiu.setdefault(pair.j_seq, []).append(rt_mx_ji)
    as_xyz = rt_mx_ji.as_xyz()
    if(not as_xyz in all_ops):
      all_ops.append(as_xyz)
      all_ops_mat.append(rt_mx_ji)
  for k,v in zip(siiu.keys(), siiu.values()): # remove duplicates!
    siiu[k] = list(set(v))
  return siiu, all_ops_mat

def apply_symop_sites_cart(sites_cart, op, fm, om):
  sites_frac = fm * sites_cart
  sites_frac_copy = flex.vec3_double(
    [op*site_frac for site_frac in sites_frac])
  return om*sites_frac_copy

def apply_symop_inplace_chain(chain, op, fm, om):
  sites_cart = chain.atoms().extract_xyz()
  sites_cart_copy = apply_symop_sites_cart(sites_cart, op, fm, om)
  chain.atoms().set_xyz(sites_cart_copy)

class manager(object):
  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               select_within_radius,
               box_buffer_layer=3,
               siiu=None,
               debug_files=False):
    # Results to be available outside
    self.super_cell_hierarchy = None # XXX sites are NOT in box center!!!
    self.super_sphere_hierarchy = None # XXX sites are NOT in box center!!!
    self.cs_super_sphere      = None
    self.siiu                 = siiu
    self.pdb_hierarchy        = pdb_hierarchy
    self.chain_op_dict        = {}
    self.fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    self.om = crystal_symmetry.unit_cell().orthogonalization_matrix()
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
        apply_symop_inplace_chain(chain=chain_dc, op=op, fm=self.fm, om=self.om)
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
    # Select the inside of the mask (discard the rest). This is super-sphere
    for chain in self.super_cell_hierarchy.chains():
      for rg in chain.residue_groups():
        sites_frac = fm_ss * rg.atoms().extract_xyz()
        inside = False
        for site_frac in sites_frac:
          if(mask.value_at_closest_grid_point(site_frac)==1):
            inside = True
            break
        if(not inside):
          #chain.remove_residue_group(rg)
          self.keep_selection.set_selected(rg.atoms().extract_i_seq(), False)
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
      sites_cart_copy = apply_symop_sites_cart(sites_cart, op, self.fm, self.om)
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
