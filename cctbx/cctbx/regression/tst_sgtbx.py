from cctbx import sgtbx
import cctbx.sgtbx.bravais_types
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

def exercise_tensor_constraints_core(crystal_symmetry):
  from cctbx import crystal
  from cctbx import adptbx
  site_symmetry = crystal.special_position_settings(
    crystal_symmetry).site_symmetry(site=(0,0,0))
  unit_cell = crystal_symmetry.unit_cell()
  group = crystal_symmetry.space_group()
  assert site_symmetry.n_matrices() == group.order_p()
  adp_constraints = group.adp_constraints()
  u_cart = adptbx.random_rotate_ellipsoid(
    u_cart=[random.random() for i in xrange(3)] + [0,0,0])
  u_star = adptbx.u_cart_as_u_star(unit_cell, u_cart)
  u_star = site_symmetry.average_u_star(u_star)
  independent_params = adp_constraints.independent_params(u_star)
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
