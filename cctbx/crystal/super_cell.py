from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
import cctbx.crystal
from libtbx import group_args
import math
import iotbx.pdb.utils

def run(pdb_hierarchy,
        crystal_symmetry,
        select_within_radius,
        box_buffer_layer=10,
        link_min=1.0,
        link_max=2.0,
        siiu=None):
  """
  Returns symmetry-expanded hierarchy within select_within_radius. The whole
  residue is included if at least one of its atoms falls within
  select_within_radius. New crystal_symmetry P1 box with dimensions bases on
  box_buffer_layer is also returned for the expanded hierarchy.
  """
  #
  def dist(r1,r2):
    return math.sqrt((r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2)
  #
  if(siiu is None):
    # This is equivalent to (but much faster):
    #
    #siiu = grm.geometry.pair_proxies().nonbonded_proxies.\
    #  get_symmetry_interacting_indices_unique(
    #    sites_cart = pdb_hierarchy.atoms().extract_xyz())
    #
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
    for pair in pair_generator:
      rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      #if str(rt_mx_ji)=="x,y,z": continue REF1
      siiu.setdefault(pair.j_seq, []).append(rt_mx_ji)
    for k,v in zip(siiu.keys(), siiu.values()): # remove duplicates!
      siiu[k] = list(set(v))
  c = iotbx.pdb.hierarchy.chain(id = "SS") # all symmetry related full residues
  fm = crystal_symmetry.unit_cell().fractionalization_matrix()
  om = crystal_symmetry.unit_cell().orthogonalization_matrix()
  selection = siiu.keys()
  for rg in pdb_hierarchy.residue_groups():
    keep=False
    for i in rg.atoms().extract_i_seq():
      if(i in selection):
        keep=True
        break
    if(keep):
      ops = siiu[i]
      for op in ops:
        rg_ = rg.detached_copy()
        xyz = rg_.atoms().extract_xyz()
        new_xyz = flex.vec3_double()
        for xyz_ in xyz:
          t1 = fm*flex.vec3_double([xyz_])
          t2 = op*t1[0]
          t3 = om*flex.vec3_double([t2])
          new_xyz.append(t3[0])
        rg_.atoms().set_xyz(new_xyz)
        rg_.link_to_previous=True
        # REF1: have to do this if comment that line
        overlap=False
        for a1 in rg.atoms():
          for a2 in rg_.atoms():
            d = dist(a1.xyz, a2.xyz)
            if(d<0.1): overlap=True
        # END REF1
        if(not overlap):
          c.append_residue_group(rg_)
  resseq = 0
  for rg in c.residue_groups():
    rg.resseq = "%4d" % resseq
    resseq+=1
  # Now we need to re-order residues and group by chains so that they can be
  # linked by pdb interpretation.
  super_sphere_hierarchy = pdb_hierarchy.deep_copy()
  #
  all_chains = iotbx.pdb.utils.all_chain_ids()
  all_chains = [chid.strip() for chid in all_chains]
  for cid in list(set([i.id for i in pdb_hierarchy.chains()])):
    assert len(cid.strip())>0, "Chain ID cannot be blank."
    all_chains.remove(cid)
  all_chains = iter(all_chains)
  #
  def get_atom(atoms, name):
    for a in atoms:
      if(a.name.strip().lower()==name.strip().lower()):
        return a.xyz
  residue_groups_i = list(c.residue_groups())
  residue_groups_j = list(c.residue_groups())
  residue_groups_k = list() # standard aa residues only
  #
  def grow_chain(residue_groups_j, chunk, ci):
    n_linked = 0
    for rgj in residue_groups_j:
      nj = get_atom(atoms=rgj.atoms(), name="N")
      if(nj is None): continue
      d_ci_nj = dist(ci,nj)
      if(d_ci_nj<link_max and d_ci_nj>link_min):
        n_linked += 1
        chunk.append(rgj)
        residue_groups_j.remove(rgj)
        break
    return n_linked, rgj
  # Find all isolated single residues, and also save starters
  starters = []
  for rgi in residue_groups_i:
    ci = get_atom(atoms=rgi.atoms(), name="C")
    ni = get_atom(atoms=rgi.atoms(), name="N")
    if(ci is None or ni is None): # collect non-proteins
      c = iotbx.pdb.hierarchy.chain(id = next(all_chains))
      c.append_residue_group(rgi)
      super_sphere_hierarchy.models()[0].append_chain(c)
      continue
    residue_groups_k.append(rgi)
    c_linked = 0
    n_linked = 0
    for rgj in residue_groups_j:
      cj = get_atom(atoms=rgj.atoms(), name="C")
      nj = get_atom(atoms=rgj.atoms(), name="N")
      if(cj is None or nj is None): continue
      d_ci_nj = dist(ci,nj)
      d_ni_cj = dist(ni,cj)
      if(d_ci_nj<link_max and d_ci_nj>link_min):
        n_linked += 1
      if(d_ni_cj<link_max and d_ni_cj>link_min):
        c_linked += 1
    assert c_linked in [1, 0], c_linked # Either linked or not!
    assert n_linked in [1, 0], n_linked # Either linked or not!
    if(c_linked==0 and n_linked==0): # isolated single residue
      c = iotbx.pdb.hierarchy.chain(id = next(all_chains))
      rgi.link_to_previous=True
      c.append_residue_group(rgi)
      super_sphere_hierarchy.models()[0].append_chain(c)
    elif(c_linked==0): # collect c-terminal residues
      starters.append(rgi)
  # Grow continuous chains from remaining residues starting from c-terminals
  for rgi in starters:
    ci = get_atom(atoms=rgi.atoms(), name="C")
    chunk = [rgi]
    n_linked=1
    while n_linked==1:
      n_linked, rgj = grow_chain(residue_groups_k, chunk, ci)
      if(n_linked==1):
        ci = get_atom(atoms=rgj.atoms(), name="C")
    c = iotbx.pdb.hierarchy.chain(id = next(all_chains))
    for i, rg in enumerate(chunk):
      rg.resseq = "%4d" % i
      c.append_residue_group(rg)
    super_sphere_hierarchy.models()[0].append_chain(c)
  #
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = super_sphere_hierarchy.atoms().extract_xyz(),
    buffer_layer = box_buffer_layer)
  cs_super_sphere = box.crystal_symmetry()
  return group_args(
    hierarchy        = super_sphere_hierarchy,
    crystal_symmetry = cs_super_sphere,
    siiu             = siiu)
