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

def show_distances(structure, distance_cutoff=5):
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
      i_seq = pairs[0].i_seq
      scatterer_i = scatterers[i_seq]
      site_i = scatterer_i.site
      site_symmetry_i = site_symmetries[i_seq]
      rt_mx_i_inverse = asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0).inverse()
      print "%s: site symmetry: %s, multiplicity: %d" % (
        scatterer_i.label,
        site_symmetry_i.point_group_type(),
        site_symmetry_i.multiplicity())
    contacts = []
    covered_already = {}
    for pair in pairs:
      assert pair.i_seq == i_seq
      scatterer_j = scatterers[pair.j_seq]
      if (pair.j_seq not in covered_already):
        covered_already[pair.j_seq] = {}
      rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
      rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
      rt_mx_jis = rt_mx_ji.multiply(site_symmetries[pair.j_seq].special_op())
      if (str(rt_mx_jis) not in covered_already[pair.j_seq]):
        ms_dict = {}
        i_sym_equiv = None
        for m in site_symmetry_i.matrices():
          ms = m.multiply(rt_mx_jis)
          ms_str = str(ms)
          if (not ms_str in ms_dict):
            ms_dict[ms_str] = 0
            covered_already[pair.j_seq][ms_str] = 0
            contact_site = ms * scatterer_j.site
            contacts.append(contact(
              i=len(contacts),
              label=scatterer_j.label,
              site=contact_site,
              rt_mx=m.multiply(rt_mx_ji),
              dist=pair.dist_sq**.5,
              i_sym_equiv=i_sym_equiv))
            if (i_sym_equiv is None):
              i_sym_equiv = contacts[-1].i
            mismatch = abs(
              structure.unit_cell().distance(site_i, contact_site)**2
              -pair.dist_sq)
            assert mismatch < asu_mappings.sym_equiv_tolerance()**2
    for c in contacts:
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
  show_distances(structure, distance_cutoff=2)

def run():
  if (0):
    test_hcp()
    return
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
