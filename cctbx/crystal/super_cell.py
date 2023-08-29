from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
import cctbx.crystal
from libtbx import group_args
import math
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
    #if str(rt_mx_ji)=="x,y,z": continue REF1
    siiu.setdefault(pair.j_seq, []).append(rt_mx_ji)
    as_xyz = rt_mx_ji.as_xyz()
    if(not as_xyz in all_ops):
      all_ops.append(as_xyz)
      all_ops_mat.append(rt_mx_ji)
  for k,v in zip(siiu.keys(), siiu.values()): # remove duplicates!
    siiu[k] = list(set(v))
  return siiu, all_ops_mat

def run(pdb_hierarchy,
        crystal_symmetry,
        select_within_radius,
        box_buffer_layer=3,
        siiu=None,
        debug_files=False):
  assert pdb_hierarchy.models_size() == 1, "one model is expected"
  # Get operators
  siiu, all_ops = get_siiu(
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
  fm = crystal_symmetry.unit_cell().fractionalization_matrix()
  om = crystal_symmetry.unit_cell().orthogonalization_matrix()
  super_cell_hierarchy = pdb_hierarchy.deep_copy()
  cntr=0
  for chain in super_cell_hierarchy.chains():
    for op in all_ops:
      chain_dc = chain.detached_copy()
      chain_dc.id = all_chains[cntr]
      cntr+=1
      # Apply symmetry operator
      sites_cart = chain_dc.atoms().extract_xyz()
      sites_frac = fm * sites_cart
      sites_frac_copy = flex.vec3_double(
        [op*site_frac for site_frac in sites_frac])
      sites_cart_copy = om*sites_frac_copy
      #
      chain_dc.atoms().set_xyz(sites_cart_copy)
      super_cell_hierarchy.models()[0].append_chain(chain_dc)
  # Move everything to the center of a big box
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = super_cell_hierarchy.atoms().extract_xyz(),
    buffer_layer = box_buffer_layer)
  super_cell_hierarchy.atoms().set_xyz(box.sites_cart)
  if(debug_files):
    super_cell_hierarchy.write_pdb_file(file_name="ss.pdb",
      crystal_symmetry = box.crystal_symmetry())
  # Create mask around focus atoms
  radii = flex.double(super_cell_hierarchy.atoms().size(), select_within_radius)
  uc = box.crystal_symmetry().unit_cell()
  fm_ss = box.crystal_symmetry().unit_cell().fractionalization_matrix()
  a,b,c = uc.parameters()[:3]
  step = 1.
  n_real = [int(a/step), int(b/step), int(c/step)]

  sel_focus = super_cell_hierarchy.atom_selection_cache().\
    selection(focus_ss)
  sites_frac_focus = fm_ss * super_cell_hierarchy.select(
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
  for chain in super_cell_hierarchy.chains():
    for rg in chain.residue_groups():
      sites_frac = fm_ss * rg.atoms().extract_xyz()
      inside = False
      for site_frac in sites_frac:
        if(mask.value_at_closest_grid_point(site_frac)==1):
          inside = True
          break
      if(not inside):
        chain.remove_residue_group(rg)
  if(debug_files):
    super_cell_hierarchy.write_pdb_file(file_name="ss_cut.pdb",
      crystal_symmetry = box.crystal_symmetry())
  # Shift back
  super_cell_hierarchy.atoms().set_xyz(
    super_cell_hierarchy.atoms().extract_xyz()-box.shift_vector)
  # Filter out overlaps (very inefficient).
  # There must be a better way of doing it!
  for i, ci in enumerate(super_cell_hierarchy.chains()):
    for j, cj in enumerate(super_cell_hierarchy.chains()):
      if j<=i or cj.id in focus_chain_ids: continue
      if not ci.id in focus_chain_ids: continue
      overlap = False
      for ai in ci.atoms():
        if overlap: break
        for aj in cj.atoms():
          if ai.distance(aj)<1.e-2:
            overlap=True
            break
      if overlap:
        super_cell_hierarchy.models()[0].remove_chain(chain=cj)
  if(debug_files):
    super_cell_hierarchy.write_pdb_file(file_name="ss_cut_shifted.pdb")
  # Box final model
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = super_cell_hierarchy.atoms().extract_xyz(),
    buffer_layer = box_buffer_layer)
  cs_super_sphere = box.crystal_symmetry()
  return group_args(
    hierarchy        = super_cell_hierarchy,# XXX sites are NOT in box center!!!
    crystal_symmetry = cs_super_sphere,
    siiu             = siiu)
