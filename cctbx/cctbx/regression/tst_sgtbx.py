from cctbx import sgtbx
from scitbx.test_utils import approx_equal
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
      == "k>=0 and (l>0 or (l=0 and h>=0))"
  assert str(i.brick()) == "0<=x<=1/2; 0<=y<1; 0<=z<1"
  assert i.wyckoff_table().space_group_type().group() == i.type().group()
  assert len(i.structure_seminvariant().vectors_and_moduli()) == 3
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
  assert str(i.primitive_setting()) == "Hall:  C 2y (-x+y,z,x+y)"

def run():
  exercise_space_group_info()
  print "OK"

if (__name__ == "__main__"):
  run()
