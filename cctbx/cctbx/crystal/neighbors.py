from iotbx.kriber import strudat
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from libtbx.itertbx import count
from libtbx.test_utils import eps_eq, approx_equal
import math
import sys

# XXX
from mmtbx.stereochemistry.tst_interaction import pml_stick, write_pml
from cctbx import sgtbx
from scitbx import matrix

def unit_cell_sticks(unit_cell, colors=[[1,1,1]]*2, width=0.05):
  ortho = unit_cell.orthogonalize
  points_frac = []
  points_cart = []
  for i in (0,1):
    for j in (0,1):
      for k in (0,1):
        v = (i,j,k)
        points_frac.append(v)
        points_cart.append(ortho(v))
  sticks = []
  for i in xrange(0,7):
    for j in xrange(i+1,8):
      v1 = points_frac[i]
      v2 = points_frac[j]
      n = 0
      for k in xrange(3):
        if (v1[k] != v2[k]): n += 1
      if (n == 1):
        sticks.append(
          pml_stick(
            begin=points_cart[i],
            end=points_cart[j],
            colors=colors,
            width=width))
  return sticks

class contact:

  def __init__(self, i, j_seq, label, site, rt_mx, dist, i_sym_equiv):
    adopt_init_args(self, locals())

  def __str__(self):
    result = " %3d: %-6s (%7.4f %7.4f %7.4f) %-25s %8.3f" % (
        (self.i, self.label)
      + self.site
      + (self.rt_mx, self.dist))
    if (self.i_sym_equiv is not None):
      result += " eq %2d" % self.i_sym_equiv
    return result

class ata_group:

  def __init__(self, mm, rt, special_op):
    unit_mx = matrix.rt(([1,0,0,0,1,0,0,0,1], [0,0,0]))
    self.a = float((unit_mx - rt.as_rational()) * special_op)
    self.rtgr = (self.a.r.transpose() * mm * self.a.r)
    self.rtgt = (self.a.r.transpose() * mm * self.a.t)
    self.ttgt = (self.a.t.transpose() * mm * self.a.t)(0,0)+1

def show_distances(structure, distance_cutoff=8):
  asu_mappings = structure.asu_mappings(
    buffer_thickness=distance_cutoff)
  cb_op_to_niggli_cell = structure.change_of_basis_op_to_niggli_cell()
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    full_matrix=0001)
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
  unit_mx = matrix.rt(([1,0,0,0,1,0,0,0,1], [0,0,0]))
  orth_rt = matrix.rt(
    (structure.unit_cell().orthogonalization_matrix(), [0,0,0]))
  mm = matrix.sym(structure.unit_cell().metrical_matrix())
  dist = structure.unit_cell().distance
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
      special_op_i = site_symmetries[i_seq].special_op().as_rational()
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
              j_seq=pair.j_seq,
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
        if (pair.j_seq == i_seq):
          s = rt_mx_ji
          while 1:
            print "sym eq inter:", s
            if (s.r().is_unit_mx()):
              break
            s = rt_mx_ji.multiply(s)
          print
        if (0 and pair.j_seq == i_seq):
          ag_ji = ata_group(mm, rt_mx_ji, special_op_i)
          site_ji = rt_mx_ji * site_i
          for rt_mx_k in structure.space_group():
            ag_k = ata_group(mm, rt_mx_k, special_op_i)
            if (approx_equal(ag_ji.rtgr, ag_k.rtgr)):
              print "CANDIDATE", i_seq, rt_mx_ji, pair.dist_sq**.5
              site_k = rt_mx_k * site_i
              print "  ji:", (matrix.col(site_ji) - matrix.col(site_i)).elems
              print "  ki:", (matrix.col(site_k) - matrix.col(site_i)).elems
              t_den = rt_mx_k.t().den()
              for u in flex.nested_loop(begin=[-3]*3,
                                        end=[3]*3, open_range=00000):
                tu = [t+x*t_den for t,x in zip(rt_mx_k.t().num(), u)]
                rt_mx_ku = sgtbx.rt_mx(rt_mx_k.r(), sgtbx.tr_vec(tu, t_den))
                if (str(rt_mx_ku) != str(rt_mx_ji)):
                  ag_ku = ata_group(mm, rt_mx_ku, special_op_i)
                  assert approx_equal(ag_ji.rtgr, ag_ku.rtgr)
                  if (    approx_equal(ag_ji.rtgt, ag_ku.rtgt)
                      and approx_equal(ag_ji.ttgt, ag_ku.ttgt)):
                    d = dist(site_i, rt_mx_ku*site_i)
                    print "HIT", rt_mx_ku, d
                    assert approx_equal(d, pair.dist_sq**.5)
    orth = structure.unit_cell().orthogonalize
    frac = structure.unit_cell().fractionalize
    mm_aug = orth_rt.as_augmented_matrix()
    mm_aug = mm_aug.transpose() * mm_aug
    color_list = []
    prev_c = None
    for c in contacts:
      print c
      if (0): continue
      colors=[[0,1,0]]*2
      if (    prev_c != None
          and c.j_seq == prev_c.j_seq
          and abs(prev_c.dist - c.dist) < 1.e-5
          and c.i_sym_equiv is None):
        if (c.j_seq == i_seq):
          atas = []
          fatas = []
          for rt in (prev_c.rt_mx, c.rt_mx):
            fdm = float((unit_mx - rt.as_rational()) * special_op_i)
            fa = fdm.as_augmented_matrix()
            fata = fa.transpose() * mm_aug * fa
            fatas.append(fata)
            ag = ata_group(mm, rt, special_op_i)
            if (0):
              print "RtGR:", ag.rtgr.elems
              print "xR:", fata.extract_block(stop=(3,3)).elems
              print "RtGT:", ag.rtgt.elems
              print "xT:", fata.extract_block(start=(0,3),stop=(3,4)).elems
              print "TtGT:", ag.ttgt
              print "x33:", fata(3,3)
            assert approx_equal(ag.rtgr, fata.extract_block(stop=(3,3)))
            assert approx_equal(ag.rtgt,
              fata.extract_block(start=(0,3),stop=(3,4)).elems)
            assert approx_equal(ag.ttgt, fata(3,3))
            dm = orth_rt * fdm
            a = dm.as_augmented_matrix()
            ata = a.transpose() * a
            atas.append(ata)
          if (0):
            for fata in fatas: print "fata:", list(fata)
            for ata in atas: print "ata:", list(ata)
          equiv_fata = approx_equal(fatas[0], fatas[1])
          equiv_ata = approx_equal(atas[0], atas[1])
          if (equiv_fata != equiv_ata):
            print "ata mismatch", structure.space_group_info()
        for shift in [0, distance_cutoff*0.1]:
          site_j_mod = site_symmetries[c.j_seq].special_op() \
            * frac([x+shift for x in orth(scatterers[c.j_seq].site)])
          if (shift == 0 or c.j_seq != i_seq):
            site_i_mod = site_i
          else:
            site_i_mod = site_j_mod
          dv1 = matrix.col(c.rt_mx * site_j_mod) - matrix.col(site_i_mod)
          dv2 = matrix.col(prev_c.rt_mx * site_j_mod) - matrix.col(site_i_mod)
          d1 = dist(c.rt_mx * site_j_mod, site_i_mod)
          d2 = dist(prev_c.rt_mx * site_j_mod, site_i_mod)
          assert abs(d1-structure.unit_cell().length(dv1)) < 1.e-6
          assert abs(d2-structure.unit_cell().length(dv2)) < 1.e-6
          if (shift == 0):
            assert abs(d1-d2) < 1.e-4
        if (abs(d1-d2) < 1.e-3):
          if (0):
            if (c.j_seq != i_seq):
              print "HARD",
          if (1):
            print "LOOK", scatterer_i.label,
          if (0):
            print site_symmetries[i_seq].point_group_type(),
            print structure.space_group_info()
            print "d1:", d1
            print "d2:", d2
          colors=[[1,0,0]]*2
          color_list[-1] = colors
          if (c.j_seq == i_seq):
            if (not equiv_ata):
              print "equiv_ata false mismatch"
            elif (1):
              print "equiv_ata true confirmation"
        else:
          if (0 or c.j_seq == i_seq):
            print "FALSE ALARM"
          if (c.j_seq == i_seq):
            if (equiv_ata):
              print "equiv_ata true mismatch"
            elif (0):
              print "equiv_ata false confirmation"
      color_list.append(colors)
      prev_c = c
    sticks = []
    for c,colors in zip(contacts, color_list):
      sticks.append(
        pml_stick(
          begin=structure.unit_cell().orthogonalize(scatterer_i.site),
          end=structure.unit_cell().orthogonalize(c.site),
          colors=colors,
          width=0.1))
    write_pml(
      f=open("contacts_%s.pml" % scatterer_i.label, "w"),
      label=scatterer_i.label,
      sticks=sticks)
  sticks = unit_cell_sticks(structure.unit_cell())
  write_pml(
    f=open("cell.pml", "w"),
    label="cell",
    sticks=sticks)

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
