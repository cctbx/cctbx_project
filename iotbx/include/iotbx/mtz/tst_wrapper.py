from iotbx import mtz
import iotbx.mtz.wrapper
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys, os

def exercise_basic():
  mtz_object = mtz.wrapper.object()
  assert mtz_object.title() == ""
  assert mtz_object.history().size() == 0
  assert mtz_object.space_group_name() == ""
  assert mtz_object.point_group_name() == ""
  assert mtz_object.space_group().order_z() == 1
  assert mtz_object.n_batches() == 0
  assert mtz_object.n_reflections() == 0
  assert mtz_object.space_group_number() == 0
  assert mtz_object.max_min_resolution() == (-1, -1)
  assert mtz_object.n_crystals() == 0
  assert mtz_object.n_active_crystals() == 0
  mtz_object = mtz.wrapper.object(
    n_datasets_for_each_crystal=flex.int([3,2,3]))
  assert mtz_object.n_crystals() == 3
  file_name = os.path.expandvars(
    "$LIBTBX_DIST_ROOT/regression/reflection_files/dano.mtz")
  if (os.path.isfile(file_name)):
    mtz_object = mtz.wrapper.object(file_name=file_name)
    assert mtz_object.title() == "......"
    assert mtz_object.history().size() == 17
    assert mtz_object.space_group_name() == "P212121"
    assert mtz_object.point_group_name() == "PG222"
    assert mtz_object.space_group().type().lookup_symbol() == "P 21 21 21"
    assert mtz_object.n_batches() == 0
    assert mtz_object.n_reflections() == 165
    assert mtz_object.space_group_number() == 19
    assert approx_equal(mtz_object.max_min_resolution(),
      (0.0044435299932956696, 0.00253282580524683))
    assert mtz_object.n_crystals() == 4
    assert mtz_object.n_active_crystals() == 3
    crystal = mtz.wrapper.crystal(mtz_object=mtz_object, i_crystal=1)
    assert crystal.mtz_object().n_reflections() == 165
    assert crystal.i_crystal() == 1
    assert mtz_object.crystals().size() == mtz_object.n_crystals()
    assert crystal.id() == 2
    assert crystal.name() == "unknown"
    assert crystal.project_name() == "unknown"
    assert approx_equal(crystal.unit_cell_parameters(),
      (84.511, 104.308, 174.103, 90, 90, 90))
    assert approx_equal(crystal.unit_cell().parameters(),
      (84.511, 104.308, 174.103, 90, 90, 90))
    assert crystal.n_datasets() == 1
    dataset = mtz.wrapper.dataset(mtz_crystal=crystal, i_dataset=0)
    assert dataset.mtz_crystal().i_crystal() == 1
    assert dataset.i_dataset() == 0
    assert dataset.id() == 1
    assert dataset.name() == "unknown230103:23:14:49"
    assert dataset.wavelength() == 0
    column = mtz.wrapper.column(mtz_dataset=dataset, i_column=0)
    assert column.mtz_dataset().mtz_crystal().i_crystal() == 1
    assert column.i_column() == 0
    assert column.mtz_crystal().i_crystal() == 1
    assert column.mtz_object().n_reflections() == 165
    assert column.label() == "H"
    assert column.type() == "H"
    assert column.is_active()
    assert column.path() == "/unknown/unknown230103:23:14:49/H"
    assert column.lookup_other("H").i_column() == 0
    assert column.lookup_other("K").i_column() == 1
    assert column.lookup_other("L").i_column() == 2
    assert column.valid_indices().size() == 165
    column = mtz_object.lookup_column("F*")
    assert column.label() == "Frem"
    expected_dataset_ids = iter(range(4))
    expected_dataset_names = iter([
      "HKL_base",
      "unknown230103:23:14:49",
      "unknown230103:23:14:21",
      "unknown230103:23:13:49"])
    expected_n_columns = iter([0,8,5,5])
    expected_column_labels = iter([
      "H", "K", "L",
      "Frem", "SIGFrem", "DANOrem", "SIGDANOrem", "ISYMrem",
      "Finf", "SIGFinf", "DANOinf", "SIGDANOinf", "ISYMinf",
      "Fabs", "SIGFabs", "DANOabs", "SIGDANOabs", "ISYMabs"])
    expected_column_types = iter("HHHFQDQYFQDQYFQDQY")
    expected_valid_indices_size = iter([
      165, 165, 165, 163, 163, 163, 163, 163, 165,
      165, 164, 164, 165, 165, 165, 164, 164, 165])
    for i_crystal,crystal in enumerate(mtz_object.crystals()):
      assert crystal.mtz_object().n_reflections() == 165
      assert crystal.i_crystal() == i_crystal
      assert crystal.n_datasets() == 1
      for i_dataset,dataset in enumerate(crystal.datasets()):
        assert dataset.mtz_crystal().i_crystal() == i_crystal
        assert dataset.i_dataset() == i_dataset
        assert dataset.id() == expected_dataset_ids.next()
        assert dataset.name() == expected_dataset_names.next()
        assert dataset.wavelength() == 0
        assert dataset.n_columns() == expected_n_columns.next()
        for i_column,column in enumerate(dataset.columns()):
          assert column.mtz_dataset().i_dataset() == i_dataset
          assert column.i_column() == i_column
          assert column.label() == expected_column_labels.next()
          assert column.type() == expected_column_types.next()
          assert column.is_active()
          assert column.path().endswith(column.label())
          assert column.valid_indices().size() \
              == expected_valid_indices_size.next()
          assert column.valid_values().size() == column.valid_indices().size()
          if (column.type() in ["H", "B", "Y", "I"]):
            assert column.valid_integers().size() \
                == column.valid_indices().size()
          lookup_column = mtz_object.lookup_column(column.label())
          assert lookup_column.label() == column.label()

def exercise():
  forever = "--forever" in sys.argv[1:]
  while 1:
    exercise_basic()
    if (not forever): break
  print "OK"

if (__name__ == "__main__"):
  exercise()
