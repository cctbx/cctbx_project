from iotbx.kriber import strudat
import iotbx.pdb
from iotbx.option_parser import iotbx_option_parser
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from scitbx import matrix
from libtbx.itertbx import count
from libtbx.test_utils import approx_equal
import math
import sys

def is_sym_equiv_interaction_simple(unit_cell,
                                    i_seq,
                                    site_i,
                                    j_seq,
                                    site_j,
                                    special_op_j,
                                    rt_mx_ji_1,
                                    rt_mx_ji_2):
  f = unit_cell.shortest_vector_sq()**.5*.1
  trial_shifts = [f*x for x in [math.sqrt(2),math.sqrt(3),math.sqrt(5)]]
  frac = unit_cell.fractionalize
  orth = unit_cell.orthogonalize
  dist = unit_cell.distance
  for shifts in [[0,0,0], trial_shifts]:
    site_j_mod = special_op_j * frac([x+s for x,s in zip(orth(site_j),shifts)])
    if (shifts == [0,0,0] or j_seq != i_seq):
      site_i_mod = site_i
    else:
      site_i_mod = site_j_mod
    d1 = dist(rt_mx_ji_1 * site_j_mod, site_i_mod)
    d2 = dist(rt_mx_ji_2 * site_j_mod, site_i_mod)
    if (shifts == [0,0,0]):
      assert abs(d1-d2) < 4.e-4
  return abs(d1-d2) < 1.e-3

class ata_homogeneous:

  def __init__(self, mm, rt, special_op):
    unit_mx = matrix.rt(([1,0,0,0,1,0,0,0,1], [0,0,0]))
    self.a = float((unit_mx - rt.as_rational()) * special_op.as_rational())
    self.rtgr = (self.a.r.transpose() * mm * self.a.r)
    self.rtgt = (self.a.r.transpose() * mm * self.a.t)
    self.ttgt = (self.a.t.transpose() * mm * self.a.t)(0,0)+1

  def as_tuple(self):
    return self.rtgr.elems + self.rtgt.elems + (self.ttgt,)

def is_sym_equiv_interaction_homogeneous(unit_cell,
                                         special_op_j,
                                         rt_mx_ji_1,
                                         rt_mx_ji_2):
  """
     <rt_mx_ji_1 * special_op_j * (site_j+x) - (site_j+x)>
  == <rt_mx_ji_2 * special_op_j * (site_j+x) - (site_j+x)>
     for all x element of R3
  """
  mm = matrix.sym(unit_cell.metrical_matrix())
  atas = [ata_homogeneous(mm, rt, special_op_j)
    for rt in [rt_mx_ji_1, rt_mx_ji_2]]
  return approx_equal(atas[0].as_tuple(), atas[1].as_tuple())

class ata_heterogeneous:

  def __init__(self, site_i, mm, rt, special_op):
    rtso = float((rt.multiply(special_op).as_rational()))
    site_i = matrix.col(site_i)
    # change of basis to move the origin to site_i:
    #   (I,-site_i)*(rtso.r,rtso.t)*(I,site_i)
    #   = (rtso.r, rtso.t-site_i)*(I,site_i)
    #   = (rtso.r, rtso.r*site_i+rtso.t-site_i)
    self.a = matrix.rt((rtso.r, rtso.r*site_i+rtso.t-site_i))
    self.rtgr = (self.a.r.transpose() * mm * self.a.r)
    self.rtgt = (self.a.r.transpose() * mm * self.a.t)
    self.ttgt = (self.a.t.transpose() * mm * self.a.t)(0,0)+1

  def as_tuple(self):
    return self.rtgr.elems + self.rtgt.elems + (self.ttgt,)

def is_sym_equiv_interaction_heterogeneous(unit_cell,
                                           site_i,
                                           special_op_j,
                                           rt_mx_ji_1,
                                           rt_mx_ji_2):
  """
     <rt_mx_ji_1 * special_op_j * (site_j+x) - site_i>
  == <rt_mx_ji_2 * special_op_j * (site_j+x) - site_i>
  for all x element of R3
  """
  mm = matrix.sym(unit_cell.metrical_matrix())
  atas = [ata_heterogeneous(site_i, mm, rt, special_op_j)
    for rt in [rt_mx_ji_1, rt_mx_ji_2]]
  return approx_equal(atas[0].as_tuple(), atas[1].as_tuple())

def is_sym_equiv_interaction(unit_cell,
                             i_seq,
                             site_i,
                             j_seq,
                             site_j,
                             special_op_j,
                             rt_mx_ji_1,
                             rt_mx_ji_2):
  if (i_seq == j_seq):
    return is_sym_equiv_interaction_homogeneous(
      unit_cell, special_op_j, rt_mx_ji_1, rt_mx_ji_2)
  else:
    return is_sym_equiv_interaction_heterogeneous(
      unit_cell, site_i, special_op_j, rt_mx_ji_1, rt_mx_ji_2)

class contact:

  def __init__(self, i, j_seq, label, site, rt_mx, dist):
    adopt_init_args(self, locals())

  def __str__(self):
    result = " %3d: %-6s (%7.4f %7.4f %7.4f) %-25s %8.3f" % (
        (self.i, self.label)
      + self.site
      + (self.rt_mx, self.dist))
    return result

def show_distances(structure, distance_cutoff=3.5):
  distance_cutoff_plus = distance_cutoff * (1+1.e-4)
  distance_cutoff_minus = distance_cutoff * (1-1.e-4)
  asu_mappings = structure.asu_mappings(
    buffer_thickness=distance_cutoff_plus)
  cb_op_to_niggli_cell = structure.change_of_basis_op_to_niggli_cell()
  pair_generator_simple = crystal.neighbors_simple_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff_plus)
  pair_generator_fast = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff_plus)
  n_simple = pair_generator_simple.count_pairs()
  n_fast = pair_generator_fast.count_pairs()
  assert n_fast == n_simple, (n_fast, n_simple)
  pair_generator = pair_generator_fast
  pair_generator.restart()
  distances_list = [flex.double()
    for i in xrange(structure.scatterers().size())]
  pairs_list = [[]
    for i in xrange(structure.scatterers().size())]
  for pair in pair_generator:
    diff_vec_frac_niggli = matrix.col(
        cb_op_to_niggli_cell.c()
      * structure.unit_cell().fractionalize(pair.diff_vec))
    if (diff_vec_frac_niggli.each_abs().max() > 1.+1.e-6):
      continue
    distances_list[pair.i_seq].append(pair.dist_sq)
    pairs_list[pair.i_seq].append(pair)
  for i,distances,pairs in zip(count(),distances_list,pairs_list):
    perm = flex.sort_permutation(distances)
    pairs_list[i] = flex.select(pairs, perm)
  site_symmetries = []
  scatterers = structure.scatterers()
  for scatterer in scatterers:
    site_symmetries.append(structure.site_symmetry(scatterer.site))
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
    for i_pair,pair in zip(count(), pairs):
      assert pair.i_seq == i_seq
      if (pair.dist_sq**.5 > distance_cutoff):
        break
      scatterer_j = scatterers[pair.j_seq]
      rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
      rt_mx_ji = rt_mx_i_inverse.multiply(rt_mx_j)
      contacts.append(contact(
        i=len(contacts),
        j_seq=pair.j_seq,
        label=scatterer_j.label,
        site=rt_mx_ji*scatterer_j.site,
        rt_mx=rt_mx_ji,
        dist=pair.dist_sq**.5))
      if (pair.dist_sq**.5 > distance_cutoff_minus):
        for pair_fwd in pairs[i_pair+1:]:
          if (    pair_fwd.j_seq == pair.j_seq
              and abs(pair_fwd.dist_sq**.5-pair.dist_sq**.5) < 1.e-4):
            rt_mx_j_fwd = asu_mappings.get_rt_mx(
              i_seq=pair_fwd.j_seq, i_sym=pair_fwd.j_sym)
            rt_mx_ji_fwd = rt_mx_i_inverse.multiply(rt_mx_j_fwd)
            if (is_sym_equiv_interaction_simple(
                  unit_cell=structure.unit_cell(),
                  i_seq=i_seq,
                  site_i=site_i,
                  j_seq=pair.j_seq,
                  site_j=scatterer_j.site,
                  special_op_j=site_symmetries[pair.j_seq].special_op(),
                  rt_mx_ji_1=rt_mx_ji,
                  rt_mx_ji_2=rt_mx_ji_fwd)):
              contacts.append(contact(
                i=len(contacts),
                j_seq=pair_fwd.j_seq,
                label=scatterer_j.label,
                site=rt_mx_ji_fwd*scatterer_j.site,
                rt_mx=rt_mx_ji_fwd,
                dist=pair.dist_sq**.5))
    prev_c = None
    for curr_c in contacts:
      print curr_c
      if (1 and prev_c is not None
            and prev_c.j_seq == curr_c.j_seq
            and abs(prev_c.dist-curr_c.dist)<1.e-4):
        is_sym_equiv_results = []
        for func in [is_sym_equiv_interaction_simple,is_sym_equiv_interaction]:
          is_sym_equiv_results.append(func(
                unit_cell=structure.unit_cell(),
                i_seq=i_seq,
                site_i=site_i,
                j_seq=curr_c.j_seq,
                site_j=scatterers[curr_c.j_seq].site,
                special_op_j=site_symmetries[curr_c.j_seq].special_op(),
                rt_mx_ji_1=prev_c.rt_mx,
                rt_mx_ji_2=curr_c.rt_mx))
        if (is_sym_equiv_results.count(00000) != 2):
          print "is_sym_equiv_interaction"
        if (is_sym_equiv_results.count(00000) % 2 != 0):
          print "is_sym_equiv mismatch:", is_sym_equiv_results
      prev_c = curr_c

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
  show_distances(structure, distance_cutoff=1.9999999999999992)

def run():
  command_line = (iotbx_option_parser(
    usage="python neighbors.py [options] studat_file [...]",
    description="Example: python neighbors.py strudat --tag=SOD")
    .enable_symmetry_comprehensive()
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
  ).process()
  if (len(sys.argv) == 1):
    test_hcp()
    return
  for file_name in sys.argv[1:]:
    try:
      strudat_entries = strudat.read_all_entries(open(file_name))
    except:
      strudat_entries = None
    if (strudat_entries is not None and len(strudat_entries.entries) > 0):
      for entry in strudat_entries.entries:
        print "strudat tag:", entry.tag
        structure = entry.as_xray_structure()
        structure.show_summary().show_scatterers()
        show_distances(structure=structure)
        print
    else:
      try:
        structure = iotbx.pdb.as_xray_structure(
          file_name=file_name,
          crystal_symmetry=command_line.symmetry)
      except:
        raise RuntimeError("Coordinate file %s: unknown format." % file_name)
      else:
        structure.show_summary().show_scatterers()
        show_distances(structure=structure)

if (__name__ == "__main__"):
  run()
