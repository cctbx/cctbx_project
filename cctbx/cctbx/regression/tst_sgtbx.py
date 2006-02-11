from cctbx import sgtbx
import cctbx.sgtbx.bravais_types
from cctbx.array_family import flex
import scitbx.math
from libtbx.test_utils import approx_equal
import random
import pickle

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
  assert str(i.primitive_setting()) == "Hall:  C 2y (-x+y,z,x+y)"
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

def run():
  exercise_space_group_info()
  test_enantiomorphic_pairs()
  exercise_tensor_constraints()
  print "OK"

if (__name__ == "__main__"):
  run()
