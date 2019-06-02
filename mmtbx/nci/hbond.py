from __future__ import division
import iotbx.pdb
import mmtbx.model
import math, sys
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
        a_YAH_cutoff = [90, 180],  # should be within this interval, not used
        protein_only = False):
    self.result = []
    self.atoms = model.get_hierarchy().atoms()
    geometry = model.get_restraints_manager()
    bond_proxies_simple, asu = geometry.geometry.get_all_bond_proxies(
      sites_cart = model.get_sites_cart())
    h_bonded_to = {}
    for p in bond_proxies_simple:
      i, j = p.i_seqs
      ei, ej = self.atoms[p.i_seqs[0]].element, self.atoms[p.i_seqs[1]].element
      if(ei in Hs): h_bonded_to[i] = self.atoms[j]
      if(ej in Hs): h_bonded_to[j] = self.atoms[i]
    #
    sites_cart = model.get_sites_cart()
    crystal_symmetry = model.crystal_symmetry()
    fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    om = crystal_symmetry.unit_cell().orthogonalization_matrix()
    pg = get_pair_generator(
      crystal_symmetry = crystal_symmetry,
      buffer_thickness = d_HA_cutoff[1],
      sites_cart       = sites_cart)
    get_class = iotbx.pdb.common_residue_names_get_class
    for p in pg.pair_generator:
      i, j = p.i_seq, p.j_seq
      ei, ej = self.atoms[i].element, self.atoms[j].element
      altloc_i = self.atoms[i].parent().altloc
      altloc_j = self.atoms[j].parent().altloc
      resseq_i = self.atoms[i].parent().parent().resseq
      resseq_j = self.atoms[j].parent().parent().resseq
      # pre-screen candidates begin
      one_is_Hs = ei in Hs or ej in Hs
      other_is_acceptor = ei in As or ej in As
      d_HA = math.sqrt(p.dist_sq)
      assert d_HA <= d_HA_cutoff[1]
      is_candidate = one_is_Hs and other_is_acceptor and \
        d_HA >= d_HA_cutoff[0] and \
        altloc_i == altloc_j and resseq_i != resseq_j
      if(protein_only):
        for it in [i,j]:
          resname = self.atoms[it].parent().resname
          is_candidate &= get_class(name=resname) == "common_amino_acid"
      if(not is_candidate): continue
      if(ei in Hs and not h_bonded_to[i].element in As): continue
      if(ej in Hs and not h_bonded_to[j].element in As): continue
      # pre-screen candidates end
      rt_mx_i = pg.conn_asu_mappings.get_rt_mx_i(p)
      rt_mx_j = pg.conn_asu_mappings.get_rt_mx_j(p)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      #
      if(ei in Hs):
        H = self.atoms[i]
        D = self.atoms[h_bonded_to[H.i_seq].i_seq]
        A = self.atoms[j]
        if(str(rt_mx_ji) != "x,y,z"):
          A = apply_symop_to_copy(A, rt_mx_ji, fm, om)
      if(ej in Hs):
        H = self.atoms[j]
        D = self.atoms[h_bonded_to[H.i_seq].i_seq]
        A = self.atoms[i]
        if(str(rt_mx_ji) != "x,y,z"):
          H = apply_symop_to_copy(H, rt_mx_ji, fm, om)
          D = apply_symop_to_copy(D, rt_mx_ji, fm, om)
      assert H.distance(D) < 1.15
      # filter by a_DHA
      a_DHA = H.angle(A, D, deg=True)
      if(a_DHA < a_DHA_cutoff): continue
      #
      assert approx_equal(d_HA, H.distance(A), 1.e-3)
      self.result.append(group_args(
        i       = i,
        j       = j,
        symop   = rt_mx_ji,
        d_HA    = d_HA,
        a_DHA   = a_DHA,
        d_AD    = A.distance(D)
      ))

  def show(self, log = sys.stdout):
    for r in self.result:
      ids_i = self.atoms[r.i].id_str().replace("pdb=","")
      ids_j = self.atoms[r.j].id_str().replace("pdb=","")
      print >> log, "%4d %4d"%(r.i,r.j), ids_i, ids_j, "d_HA=%5.3f"%r.d_HA,\
        "d_AD=%5.3f"%r.d_AD, "a_DHA=%7.3f"%r.a_DHA, "symop: %s"%str(r.symop)

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
