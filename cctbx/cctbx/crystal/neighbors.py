from iotbx.kriber import strudat
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from libtbx.itertbx import count
from libtbx.test_utils import eps_eq
import math
import sys

class contact:

  def __init__(self, i, label, site, rt_mx, dist, i_sym_equiv):
    adopt_init_args(self, locals())

  def __str__(self):
    result = "  %2d: %-6s (%7.4f %7.4f %7.4f) %-25s %8.3f" % (
        (self.i, self.label)
      + self.site
      + (self.rt_mx, self.dist))
    if (self.i_sym_equiv is not None):
      result += " eq %2d" % self.i_sym_equiv
    return result

class contacts:

  def __init__(self):
    self.table = []
    self.unique = {}
    self._sort_keys = flex.double()

  def append(self, c):
    self.table.append(c)
    if (c.i_sym_equiv is None):
      i = len(self.unique)
      self.unique[c.i] = i
      key = 1000 * i
    else:
      assert c.i_sym_equiv < 1000
      key = self.unique[c.i_sym_equiv]*1000 + c.i_sym_equiv + 1
    self._sort_keys.append(key)

  def sort_in_place(self):
    perm = flex.sort_permutation(self._sort_keys)
    self.table = flex.select(self.table, perm)
    self._sort_keys = self._sort_keys.select(perm)
    rev_perm = flex.sort_permutation(perm.as_double())
    for i,c in zip(count(), self.table):
      c.i = i
      if (c.i_sym_equiv is not None):
        c.i_sym_equiv = rev_perm[c.i_sym_equiv]
    new_unique = {}
    for k,v in self.unique.items():
      new_unique[rev_perm[k]] = v
    self.unique = new_unique

def show_distances(structure, distance_cutoff=5, max_show_neighbors=15):
  asu_mappings = structure.asu_mappings(
    buffer_thickness=distance_cutoff)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    full_matrix=0001)
  distances_list = [flex.double()
    for i in xrange(structure.scatterers().size())]
  pairs_list = [[]
    for i in xrange(structure.scatterers().size())]
  for pair in pair_generator:
    distances_list[pair.i_seq].append(pair.dist_sq)
    pairs_list[pair.i_seq].append(pair)
  for i,distances,pairs in zip(count(),distances_list,pairs_list):
    perm = flex.sort_permutation(distances)
    pairs_list[i] = flex.select(pairs, perm)
  site_symmetries = []
  for scatterer in structure.scatterers():
    site_symmetries.append(structure.site_symmetry(scatterer.site))
  scatterers = structure.scatterers()
  for pairs in pairs_list:
    if (len(pairs) > 0):
      site_i_orig = scatterers[pairs[0].i_seq].site
      site_symmetry_i = site_symmetries[pairs[0].i_seq]
      rt_mx_i_inverse = asu_mappings.get_rt_mx(
        i_seq=pairs[0].i_seq, i_sym=0).inverse()
      print "%s: site symmetry: %s, multiplicity: %d" % (
        scatterers[pairs[0].i_seq].label,
        site_symmetry_i.point_group_type(),
        site_symmetry_i.multiplicity())
    contacts_ = contacts()
    all_neighbors = []
    for i_pair,pair in zip(count(), pairs[:max_show_neighbors]):
      rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
      rt_mx_j_contact_i_orig = rt_mx_i_inverse.multiply(rt_mx_j)
      site_j_contact_i_orig = rt_mx_j_contact_i_orig \
                            * scatterers[pair.j_seq].site
      special_op_j = rt_mx_j_contact_i_orig.multiply(
        site_symmetries[pair.j_seq].special_op())
      mismatch = abs(
        structure.unit_cell().distance(site_i_orig, site_j_contact_i_orig)**2
        -pair.dist_sq)
      assert mismatch < asu_mappings.sym_equiv_tolerance()**2, mismatch**.5
      i_sym_equiv = None
      prev_len_all_neighbors = len(all_neighbors)
      for m in site_symmetry_i.matrices():
        ms = m * site_j_contact_i_orig
        for i_neighbor,neighbor in zip(count(), all_neighbors):
          i,j,r,s = neighbor
          if (j == pair.j_seq):
            d = structure.unit_cell().distance(s, ms)
            if (d < asu_mappings.sym_equiv_minimum_distance()):
              x1 = r * scatterers[pair.j_seq].site
              x2 = m.multiply(special_op_j) * scatterers[pair.j_seq].site
              assert eps_eq(x1, x2)
              assert r == m.multiply(special_op_j)
              if (i_neighbor < prev_len_all_neighbors):
                i_sym_equiv = i
                break
            else:
              assert r != m.multiply(special_op_j)
        if (i_sym_equiv is not None):
          break
        else:
          all_neighbors.append(
            (i_pair,pair.j_seq,m.multiply(special_op_j),ms))
      contacts_.append(contact(
        i=i_pair,
        label=scatterers[pair.j_seq].label,
        site=site_j_contact_i_orig,
        rt_mx=rt_mx_j_contact_i_orig,
        dist=pair.dist_sq**.5,
        i_sym_equiv=i_sym_equiv))
    contacts_.sort_in_place()
    for c in contacts_.table:
      print c

def test_hcp():
  a = 2
  c = a*math.sqrt(8/3.)
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(a,a,c,90,90,120),
        space_group_symbol="P63/mmc"),
      min_distance_sym_equiv=0.1))
  for i,site in zip(count(1),[(1/3.,2/3.,0.25)]):
    structure.add_scatterer(
      xray.scatterer(label="S%d" % i, site=site))
  structure.show_summary().show_scatterers()
  show_distances(structure)

def run():
  for file_name in sys.argv[1:]:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      print "strudat tag:", entry.tag
      structure = entry.as_xray_structure()
      structure.show_summary().show_scatterers()
      show_distances(structure=structure)
      print

if (__name__ == "__main__"):
  run()
