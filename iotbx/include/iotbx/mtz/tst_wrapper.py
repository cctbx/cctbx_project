from iotbx import mtz
import iotbx.mtz.wrapper
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys, os

def exercise_basic():
  mtz_object = mtz.wrapper.object()
  assert mtz_object.n_batches() == 0
  assert mtz_object.n_reflections() == 0
  assert mtz_object.space_group_number() == 0
  assert mtz_object.max_min_resolution() == (-1, -1)
  assert mtz_object.n_crystals() == 0
  assert mtz_object.n_active_crystals() == 0
  mtz_object = mtz.wrapper.object(n_datasets_for_each_crystal=flex.int([3,2,3]))
  assert mtz_object.n_crystals() == 3
  file_name = os.path.expandvars(
    "$LIBTBX_DIST_ROOT/regression/reflection_files/dano.mtz")
  if (os.path.isfile(file_name)):
    mtz_object = mtz.wrapper.object(file_name=file_name)
    assert mtz_object.n_batches() == 0
    assert mtz_object.n_reflections() == 165
    assert mtz_object.space_group_number() == 19
    assert approx_equal(mtz_object.max_min_resolution(),
      (0.0044435299932956696, 0.00253282580524683))
    assert mtz_object.n_crystals() == 4
    assert mtz_object.n_active_crystals() == 3
    crystal = mtz.wrapper.crystal(mtz_object, 0)
    assert crystal.object().n_reflections() == 165
    assert crystal.i_crystal() == 0
    assert mtz_object.crystals().size() == mtz_object.n_crystals()
    for i_crystal,crystal in enumerate(mtz_object.crystals()):
      assert crystal.object().n_reflections() == 165
      assert crystal.i_crystal() == i_crystal

def exercise():
  forever = "--forever" in sys.argv[1:]
  while 1:
    exercise_basic()
    if (not forever): break
  print "OK"

if (__name__ == "__main__"):
  exercise()
