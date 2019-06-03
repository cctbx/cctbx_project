from __future__ import division
import iotbx.pdb
import mmtbx.model
import math, sys, os
from libtbx import group_args
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def get_pair_generator(crystal_symmetry, buffer_thickness, sites_cart):
  sst = crystal_symmetry.special_position_settings().site_symmetry_table(
    sites_cart = sites_cart)
  from cctbx import crystal
  conn_asu_mappings = crystal_symmetry.special_position_settings().\
    asu_mappings(buffer_thickness = buffer_thickness)
  conn_asu_mappings.process_sites_cart(
    original_sites      = sites_cart,
    site_symmetry_table = sst)
  conn_pair_asu_table = crystal.pair_asu_table(asu_mappings = conn_asu_mappings)
  conn_pair_asu_table.add_all_pairs(distance_cutoff = buffer_thickness)
  pair_generator = crystal.neighbors_fast_pair_generator(
    conn_asu_mappings, distance_cutoff = buffer_thickness)
  return group_args(
    pair_generator    = pair_generator,
    conn_asu_mappings = conn_asu_mappings)

def apply_symop_to_copy(atom, rt_mx_ji, fm, om):
  atom = atom.detached_copy()
  t1 = fm*flex.vec3_double([atom.xyz])
  t2 = rt_mx_ji*t1[0]
  t3 = om*flex.vec3_double([t2])
  atom.set_xyz(t3[0])
  return atom

def make_atom_id(atom, index):
  return group_args(
    id_str = atom.id_str().replace("pdb=",""),
    index  = index,
    name   = atom.name,
    chain  = atom.parent().parent().parent().id,
    resseq = atom.parent().parent().resseq,
    altloc = atom.parent().altloc)

class find(object):
  """
     Y
      \
       A
        .
         .
         H
         |
         D
        / \

    A = O, N, S
    D = O, N, S
    90 <= a_YAH <= 180
    a_DHA >= 120
    1.4 <= d_HA <= 3.0
    2.5 <= d_DA <= 3.5
  """
  def __init__(self,
        model,
        Hs           = ["H", "D"],
        As           = ["O","N","S","F","CL"],
        Ds           = ["O","N","S"],
        d_HA_cutoff  = [1.4, 3.0], # original: [1.4, 2.4],
        d_DA_cutoff  = [2.5, 3.5], # not used
        a_DHA_cutoff = 120,        # should be greater than this
        a_YAH_cutoff = [90, 180],  # should be within this interval
        protein_only = False,
        pair_proxies = None
        ):
    self.result = []
    self.model = model
    self.pair_proxies = pair_proxies
    self.external_proxies = False
    if(self.pair_proxies is not None):
      self.external_proxies = True
    atoms = self.model.get_hierarchy().atoms()
    geometry = self.model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
      sites_cart = self.model.get_sites_cart())
    h_bonded_to = {}
    a_bonded_to = {}
    for p in bond_proxies_simple:
      i, j = p.i_seqs
      ei, ej = atoms[p.i_seqs[0]].element, atoms[p.i_seqs[1]].element
      if(ei in Hs): h_bonded_to[i] = atoms[j]
      if(ej in Hs): h_bonded_to[j] = atoms[i]
      # collect all bonded to A to make sure A has only one bonded
      if(ei in As): a_bonded_to.setdefault(i, []).append(atoms[j])
      if(ej in As): a_bonded_to.setdefault(j, []).append(atoms[i])
    sites_cart = self.model.get_sites_cart()
    crystal_symmetry = self.model.crystal_symmetry()
    fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    om = crystal_symmetry.unit_cell().orthogonalization_matrix()
    pg = get_pair_generator(
      crystal_symmetry = crystal_symmetry,
      buffer_thickness = d_HA_cutoff[1],
      sites_cart       = sites_cart)
    get_class = iotbx.pdb.common_residue_names_get_class
    # Find proxies if not provided!
    if(self.pair_proxies is None):
      pp = []
      self.pair_proxies = []
      pp = [p for p in pg.pair_generator]
    else:
      pp = self.pair_proxies

    for p in pp:
      i, j = p.i_seq, p.j_seq
      ei, ej = atoms[i].element, atoms[j].element
      altloc_i = atoms[i].parent().altloc
      altloc_j = atoms[j].parent().altloc
      resseq_i = atoms[i].parent().parent().resseq
      resseq_j = atoms[j].parent().parent().resseq
      # pre-screen candidates begin
      one_is_Hs = ei in Hs or ej in Hs
      other_is_acceptor = ei in As or ej in As
      is_candidate = one_is_Hs and other_is_acceptor and \
        altloc_i == altloc_j and resseq_i != resseq_j
      if(protein_only):
        for it in [i,j]:
          resname = atoms[it].parent().resname
          is_candidate &= get_class(name=resname) == "common_amino_acid"
      if(not is_candidate): continue
      if(ei in Hs and not h_bonded_to[i].element in As): continue
      if(ej in Hs and not h_bonded_to[j].element in As): continue
      # pre-screen candidates end
      # symop tp map onto symmetry related
      rt_mx_i = pg.conn_asu_mappings.get_rt_mx_i(p)
      rt_mx_j = pg.conn_asu_mappings.get_rt_mx_j(p)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      #
      Y = None
      if(ei in Hs):
        H = atoms[i]
        D = atoms[h_bonded_to[H.i_seq].i_seq]
        A = atoms[j]
        if(len(a_bonded_to[j])==1):
          Y = atoms[a_bonded_to[j][0].i_seq]
        atom_i = make_atom_id(atom = H, index = i)
        atom_j = make_atom_id(atom = A, index = j)
        if(str(rt_mx_ji) != "x,y,z"):
          A = apply_symop_to_copy(A, rt_mx_ji, fm, om)
          if(Y is not None):
            Y = apply_symop_to_copy(Y, rt_mx_ji, fm, om)
      if(ej in Hs):
        H = atoms[j]
        D = atoms[h_bonded_to[H.i_seq].i_seq]
        A = atoms[i]
        if(len(a_bonded_to[i])==1):
          Y = atoms[a_bonded_to[i][0].i_seq]
        atom_i = make_atom_id(atom = A, index = i)
        atom_j = make_atom_id(atom = H, index = j)
        if(str(rt_mx_ji) != "x,y,z"):
          H = apply_symop_to_copy(H, rt_mx_ji, fm, om)
          D = apply_symop_to_copy(D, rt_mx_ji, fm, om)
      d_HA = A.distance(H)
      if(not self.external_proxies):
        assert d_HA <= d_HA_cutoff[1]
        assert approx_equal(math.sqrt(p.dist_sq), d_HA, 1.e-3)
#      assert H.distance(D) < 1.15, [H.distance(D), H.name, D.name]
      # filter by a_DHA
      a_DHA = H.angle(A, D, deg=True)
      if(not self.external_proxies):
        if(a_DHA < a_DHA_cutoff): continue
      # filter by a_YAH
      a_YAH = None
      if(Y is not None):
        a_YAH = A.angle(Y, H, deg=True)
        if(not self.external_proxies):
          if(not (a_YAH >= a_YAH_cutoff[0] and a_YAH <= a_YAH_cutoff[1])):
            continue
      #
      assert approx_equal(d_HA, H.distance(A), 1.e-3)
      self.result.append(group_args(
        i       = i,
        j       = j,
        atom_i  = atom_i,
        atom_j  = atom_j,
        symop   = rt_mx_ji,
        d_HA    = d_HA,
        a_DHA   = a_DHA,
        a_YAH   = a_YAH,
        d_AD    = A.distance(D)
      ))
      if(not self.external_proxies):
        self.pair_proxies.append(p)

  def show(self, log = sys.stdout):
    for r in self.result:
      ids_i = r.atom_i.id_str
      ids_j = r.atom_j.id_str
      print >> log, "%4d %4d"%(r.i,r.j), "%s<>%s"%(ids_i, ids_j), \
        "d_HA=%5.3f"%r.d_HA, "d_AD=%5.3f"%r.d_AD, "a_DHA=%7.3f"%r.a_DHA, \
        "symop: %s"%str(r.symop)

  def as_pymol(self, prefix="hbonds_pymol"):
    pdb_file_name = "%s.pdb"%prefix
    of = open(pdb_file_name, "w")
    print >> of, self.model.model_as_pdb()
    of.close()
    of = open("%s.pml"%prefix, "w")
    print >> of, "load", "/".join([os.getcwd(), pdb_file_name])
    for r in self.result:
      if(str(r.symop) != "x,y,z"): continue
      ai = r.atom_i
      aj = r.atom_j
      one = "chain %s and resi %s and name %s and alt '%s'"%(
        ai.chain, ai.resseq, ai.name, ai.altloc)
      two = "chain %s and resi %s and name %s and alt '%s'"%(
        aj.chain, aj.resseq, aj.name, aj.altloc)
      print >> of, "dist %s, %s"%(one, two)
    of.close()

  def as_restraints(self, file_name, distance_ideal=None, sigma_dist=0.1,
       angle_ideal = None, sigma_angle=2):
    base = """bond{
      atom_selection_1 = %s
      atom_selection_2 = %s
      symmetry_operation = None
      distance_ideal = %f
      sigma = %f

    }
    angle {
      atom_selection_1 = %s
      atom_selection_2 = %s
      atom_selection_3 = %s
      angle_ideal = %f
      sigma = %f
    }
    """
    top = """refinement{
  geometry_restraints.edits{
    %s
  }
}
    """
    raise Sorry("Not Implemented.")
