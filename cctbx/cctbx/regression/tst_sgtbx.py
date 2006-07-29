from cctbx import sgtbx
import cctbx.sgtbx.bravais_types
from cctbx.array_family import flex
import scitbx.math
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random
import pickle
import sys

def exercise_space_group_info():
  i = sgtbx.space_group_info("P 1")
  assert i.type().number() == 1
  i = sgtbx.space_group_info("P -1")
  assert i.type().number() == 2
  i = sgtbx.space_group_info("P 2", "I")
  assert str(i) == "P 1 1 2"
  i = sgtbx.space_group_info("P 2", "a")
  assert str(i) == "P 1 2 1"
  assert i.group() == i.type().group()
  assert i.reciprocal_space_asu().reference_as_string() \
      == "k>=0 and (l>0 or (l==0 and h>=0))"
  assert str(i.brick()) == "0<=x<=1/2; 0<=y<1; 0<=z<1"
  assert i.wyckoff_table().space_group_type().group() == i.type().group()
  assert len(i.structure_seminvariants().vectors_and_moduli()) == 3
  for sg_number in (1,3,15,75,143,195):
    assert approx_equal(
      sgtbx.space_group_info(sg_number).any_compatible_unit_cell(100).volume(),
      100)
  s = pickle.dumps(i)
  j = pickle.loads(s)
  assert str(i) == str(j)
  i = sgtbx.space_group_info("B 2", "i")
  assert not i.is_reference_setting()
  assert str(i.change_of_basis_op_to_reference_setting().c()) == "-x,z,y"
  assert str(i.reference_setting()) == "C 1 2 1"
  assert str(i.as_reference_setting()) == "C 1 2 1"
  assert str(i.primitive_setting()) == "C 1 2 1 (-x+y,z,x+y)"
  asu = i.direct_space_asu()
  assert len(asu.facets) == 6
  assert sgtbx.space_group(asu.hall_symbol) == i.group()
  j = i.primitive_setting()
  asu = j.direct_space_asu()
  assert len(asu.facets) == 6
  assert sgtbx.space_group(asu.hall_symbol) == j.group()
  i = sgtbx.space_group_info(number=19)
  assert [str(sgtbx.space_group_info(group=group))
    for group in i.reflection_intensity_equivalent_groups()] == [
      "P 2 2 2",
      "P 2 2 21",
      "P 21 2 2",
      "P 2 21 2",
      "P 21 21 2",
      "P 2 21 21",
      "P 21 2 21",
      "P 21 21 21"]
  assert len(i.reflection_intensity_equivalent_groups(anomalous_flag=False)) \
      == 127

def test_enantiomorphic_pairs():
  pairs = []
  done = [False for i in xrange(231)]
  for i in xrange(1, 231):
    a = sgtbx.space_group_info(i)
    b = a.change_hand()
    assert a.type().is_enantiomorphic() == b.type().is_enantiomorphic()
    assert (a.group() == b.group()) == (not a.type().is_enantiomorphic())
    done[i] = True
    if (a.type().is_enantiomorphic()):
      j = b.type().number()
      if (not done[j]):
        assert j > i
        pairs.append((i,j))
        done[i] = True
      else:
        assert j < i
        assert (j,i) in pairs
  assert pairs == [(76, 78), (91, 95), (92, 96),
                   (144, 145), (151, 153), (152, 154),
                   (169, 170), (171, 172), (178, 179), (180, 181),
                   (212, 213)]

def exercise_monoclinic_cell_choices_core(space_group_number, verbose):
  # transformation matrices for cell choices
  # columns are basis vectors "new in terms of old"
  # see Int. Tab. Vol. A, p. 22, Fig. 2.2.6.4.
  b1 = (1, 0, 0,
        0, 1, 0,
        0, 0, 1)
  b2 = (-1, 0, 1,
         0, 1, 0,
        -1, 0, 0)
  b3 = (0, 0, -1,
        0, 1,  0,
        1, 0, -1)
  flip = (0,  0, 1,
          0, -1, 0,
          1,  0, 0)
  p3s = sgtbx.space_group("P 3*")
  done = {}
  ref = sgtbx.space_group_info(number=space_group_number)
  ref_uhm = ref.type().universal_hermann_mauguin_symbol()
  for i_fl,fl in enumerate([b1, flip]):
    rfl = sgtbx.rot_mx(fl)
    cfl = sgtbx.change_of_basis_op(sgtbx.rt_mx(rfl))
    for i_rt,rt in enumerate(p3s):
      rp3 = rt.r()
      cp3 = sgtbx.change_of_basis_op(sgtbx.rt_mx(rp3))
      for i_cs,cs in enumerate([b1,b2,b3]):
        rcs = sgtbx.rot_mx(cs).inverse()
        ccs = sgtbx.change_of_basis_op(sgtbx.rt_mx(rcs))
        cb_all = cp3 * cfl * ccs
        refcb = ref.change_basis(cb_all)
        refcb2 = sgtbx.space_group_info(symbol=ref_uhm+"("+str(cb_all.c())+")")
        assert refcb2.group() == refcb.group()
        s = sgtbx.space_group_symbols(str(refcb))
        q = s.qualifier()
        hm = str(refcb)
        if (0 or verbose): print hm, q, cb_all.c()
        if (i_fl == 0):
          assert q[0] == "bca"[i_rt]
          if (len(q) == 2): assert q[1] == "123"[i_cs]
        elif (q[0] == "-"):
          assert q[1] == "bca"[i_rt]
          if (len(q) == 3): assert q[2] == "123"[i_cs]
        else:
          assert q[0] == "bca"[i_rt]
          if (len(q) == 2 and q[1] != "123"[i_cs]):
            assert done[hm] == 1
        done.setdefault(hm, 0)
        done[hm] += 1
  assert len(done) in [3, 9, 18]
  assert done.values() == [18/len(done)]*len(done)
  if (0 or verbose): print
  return done

def exercise_monoclinic_cell_choices(verbose=0):
  done = {}
  for space_group_number in xrange(3,16):
    done.update(exercise_monoclinic_cell_choices_core(
      space_group_number=space_group_number, verbose=verbose))
  assert len(done) == 105
  assert done.values().count(0) == 0
  n = 0
  for s in sgtbx.space_group_symbol_iterator():
    if (s.number() < 3): continue
    if (s.number() > 15): break
    done[s.universal_hermann_mauguin()] = 0
    n += 1
  assert n == 105
  assert done.values().count(0) == 105

def exercise_orthorhombic_hm_qualifier_as_cb_symbol():
  cb_symbols = {
    "cab": ["c,a,b", "z,x,y"],
    "a-cb": ["a,-c,b", "x,-z,y"],
    "-cba": ["-c,b,a", "-z,y,x"],
    "bca": ["b,c,a", "y,z,x"],
    "ba-c": ["b,a,-c", "y,x,-z"]}
  for sgsyms1 in sgtbx.space_group_symbol_iterator():
    n = sgsyms1.number()
    if (n < 16 or n > 74): continue
    q = sgsyms1.qualifier()
    if (len(q) == 0): continue
    e = sgsyms1.extension()
    if (e == "\0"): e = ""
    ehm = sgtbx.space_group_symbols(
      space_group_number=n, extension=e).universal_hermann_mauguin()
    cabc, cxyz = cb_symbols[q]
    assert sgtbx.change_of_basis_op(cxyz).as_abc() == cabc
    assert sgtbx.change_of_basis_op(cabc).as_xyz() == cxyz
    uhm_xyz = ehm + " ("+cxyz+")"
    sgsyms2 = sgtbx.space_group_symbols(symbol=uhm_xyz)
    assert sgsyms2.change_of_basis_symbol() == cxyz
    assert sgsyms2.extension() == sgsyms1.extension()
    assert sgsyms2.universal_hermann_mauguin() == uhm_xyz
    g1 = sgtbx.space_group(space_group_symbols=sgsyms1)
    g2 = sgtbx.space_group(space_group_symbols=sgsyms2)
    assert g2 == g1
    g2 = sgtbx.space_group(
      sgtbx.space_group_symbols(symbol=ehm)).change_basis(
        sgtbx.change_of_basis_op(sgtbx.rt_mx(cxyz)))
    assert g2 == g1
    for c in [cxyz, cabc]:
      g2 = sgtbx.space_group_info(
        group=sgtbx.space_group(
          sgtbx.space_group_symbols(symbol=ehm))).change_basis(c).group()
      assert g2 == g1
    cit = sgtbx.rt_mx(cxyz).r().inverse().transpose()
    cit_xyz = cit.as_xyz()
    g2 = sgtbx.space_group_info(
      group=sgtbx.space_group(
        sgtbx.space_group_symbols(symbol=ehm))).change_basis(cit_xyz).group()
    assert g2 == g1
    assert cit.as_xyz(False, "abc") == cabc
    uhm_abc = ehm + " ("+cabc+")"
    sgsyms2 = sgtbx.space_group_symbols(symbol=uhm_abc)
    assert sgsyms2.change_of_basis_symbol() == cxyz
    assert sgsyms2.extension() == sgsyms1.extension()
    assert sgsyms2.universal_hermann_mauguin() == uhm_xyz
    g2 = sgtbx.space_group(space_group_symbols=sgsyms2)
    assert g2 == g1

def python_tensor_constraints(self, reciprocal_space):
  """row-reduced echelon form of coefficients
       r.transpose() * t * r - t = 0
     Mathematica code:
       r={{r0,r1,r2},{r3,r4,r5},{r6,r7,r8}}
       t={{t0,t3,t4},{t3,t1,t5},{t4,t5,t2}}
       FortranForm[Expand[Transpose[r].t.r - t]]
  """
  result = flex.int()
  for s in self.smx():
    r = s.r()
    if (reciprocal_space):
      r = r.transpose()
    r0,r1,r2,r3,r4,r5,r6,r7,r8 = r.num()
    result.extend(flex.int((
      r0*r0-1, r3*r3,   r6*r6,   2*r0*r3, 2*r0*r6, 2*r3*r6,
      r1*r1,   r4*r4-1, r7*r7,   2*r1*r4, 2*r1*r7, 2*r4*r7,
      r2*r2,   r5*r5,   r8*r8-1, 2*r2*r5, 2*r2*r8, 2*r5*r8,
      r0*r1, r3*r4, r6*r7, r1*r3+r0*r4-1, r1*r6+r0*r7,   r4*r6+r3*r7,
      r0*r2, r3*r5, r6*r8, r2*r3+r0*r5,   r2*r6+r0*r8-1, r5*r6+r3*r8,
      r1*r2, r4*r5, r7*r8, r2*r4+r1*r5,   r2*r7+r1*r8,   r5*r7+r4*r8-1)))
  result.resize(flex.grid(result.size()//6,6))
  scitbx.math.row_echelon_form(result)
  return result

def exercise_tensor_constraints_core(crystal_symmetry):
  from cctbx import crystal
  from cctbx import adptbx
  from scitbx import matrix
  site_symmetry = crystal.special_position_settings(
    crystal_symmetry).site_symmetry(site=(0,0,0))
  unit_cell = crystal_symmetry.unit_cell()
  group = crystal_symmetry.space_group()
  assert site_symmetry.n_matrices() == group.order_p()
  for reciprocal_space in [False, True]:
    c_tensor_constraints = sgtbx.tensor_rank_2_constraints(
      space_group=group,
      reciprocal_space=reciprocal_space).row_echelon_form()
    p_tensor_constraints = python_tensor_constraints(
      self=group, reciprocal_space=reciprocal_space)
    assert c_tensor_constraints.all_eq(p_tensor_constraints)
  adp_constraints = group.adp_constraints()
  u_cart_p1 = adptbx.random_u_cart()
  u_star_p1 = adptbx.u_cart_as_u_star(unit_cell, u_cart_p1)
  u_star = site_symmetry.average_u_star(u_star_p1)
  f = unit_cell.volume()**(2/3.)
  assert approx_equal(
    list(matrix.col(group.average_u_star(u_star=u_star_p1))*f),
    list(matrix.col(u_star)*f))
  independent_params = adp_constraints.independent_params(u_star)
  assert adp_constraints.n_independent_params() == len(independent_params)
  assert adp_constraints.n_independent_params() \
       + adp_constraints.n_dependent_params() == 6
  u_star_vfy = adp_constraints.all_params(independent_params)
  u_cart = adptbx.u_star_as_u_cart(unit_cell, u_star)
  u_cart_vfy = adptbx.u_star_as_u_cart(unit_cell, list(u_star_vfy))
  assert approx_equal(u_cart_vfy, u_cart)

def exercise_tensor_constraints():
  from cctbx import crystal
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    crystal_symmetry = crystal.symmetry(
      unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
      space_group_info=space_group_info)
    exercise_tensor_constraints_core(crystal_symmetry)
    exercise_tensor_constraints_core(crystal_symmetry.minimum_cell())

def run(args):
  exercise_space_group_info()
  test_enantiomorphic_pairs()
  exercise_monoclinic_cell_choices(verbose="--verbose" in args)
  exercise_orthorhombic_hm_qualifier_as_cb_symbol()
  exercise_tensor_constraints()
  print format_cpu_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
