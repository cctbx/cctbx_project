from iotbx import mtz
import iotbx.mtz.wrapper
from iotbx.option_parser import iotbx_option_parser
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from cStringIO import StringIO
import sys, os

def exercise_basic():
  assert mtz.wrapper.ccp4_liberr_verbosity(-1) == 0
  assert mtz.wrapper.ccp4_liberr_verbosity(1) == 1
  assert mtz.wrapper.ccp4_liberr_verbosity(-1) == 1
  assert mtz.wrapper.ccp4_liberr_verbosity(0) == 0
  assert mtz.wrapper.ccp4_liberr_verbosity(-1) == 0
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
    assert mtz_object.valid_indices(
      column_label="Frem").size() == 163
    assert mtz_object.valid_values(
      column_label="Frem").size() == 163
    assert mtz_object.valid_integers(
      column_label="ISYMabs").size() == 165
    assert mtz_object.valid_indices_anomalous(
      column_label_plus="Frem",
      column_label_minus="DANOrem").size() == 326
    assert mtz_object.valid_values_anomalous(
      column_label_plus="Frem",
      column_label_minus="DANOrem").size() == 326
    assert mtz_object.valid_complex(
      column_label_ampl="Frem",
      column_label_phi="DANOrem").size() == 163
    assert mtz_object.valid_complex_anomalous(
      column_label_ampl_plus="Frem",
      column_label_phi_plus="DANOrem",
      column_label_ampl_minus="Frem",
      column_label_phi_minus="DANOrem").size() == 326
    array_group = mtz_object.valid_delta_anomalous(
      column_label_f_data="Frem",
      column_label_f_sigmas="SIGFrem",
      column_label_d_data="DANOrem",
      column_label_d_sigmas="SIGDANOrem")
    assert array_group.indices.size() == 326
    assert array_group.data.size() == 326
    assert array_group.sigmas.size() == 326
    assert mtz_object.valid_hl(
      column_label_a="Frem",
      column_label_b="DANOrem",
      column_label_c="Frem",
      column_label_d="DANOrem").size() == 163

def exercise_valid_any(file_name, out):
  mtz_object = mtz.wrapper.object(file_name=file_name)
  mtz_object.show_summary(out=out)
  types = "".join(mtz_object.column_types())
  i = types.find("FP")
  if (i >= 0):
    labels = mtz_object.column_labels()[i:i+2]
    mtz_object.valid_complex(
      column_label_ampl=labels[0],
      column_label_phi=labels[1])
  for type_group in ["GLGL", "KMKM"]:
    i = types.find(type_group)
    if (i >= 0):
      labels = mtz_object.column_labels()[i:i+4]
      indices = mtz_object.valid_indices_anomalous(
        column_label_plus=labels[0],
        column_label_minus=labels[2])
      data = mtz_object.valid_values_anomalous(
        column_label_plus=labels[0],
        column_label_minus=labels[2])
      assert indices.size() == data.size()
      data = mtz_object.valid_complex_anomalous(
        column_label_ampl_plus=labels[0],
        column_label_phi_plus=labels[1],
        column_label_ampl_minus=labels[2],
        column_label_phi_minus=labels[3])
      assert indices.size() == data.size()
  i = types.find("FQDQ")
  if (i >= 0):
    labels = mtz_object.column_labels()[i:i+4]
    mtz_object.valid_delta_anomalous(
      column_label_f_data=labels[0],
      column_label_f_sigmas=labels[1],
      column_label_d_data=labels[2],
      column_label_d_sigmas=labels[3])
  i = types.find("AAAA")
  if (i >= 0):
    labels = mtz_object.column_labels()[i:i+4]
    mtz_object.valid_hl(
      column_label_a=labels[0],
      column_label_b=labels[1],
      column_label_c=labels[2],
      column_label_d=labels[3])

def walk_callback(arg, top, names):
  out = arg
  for name in names:
    if (not name.endswith(".mtz")): continue
    file_name = os.path.normpath(os.path.join(top, name))
    print >> out, "Processing:", file_name
    exercise_valid_any(file_name=file_name, out=out)

def exercise_walk(root_dir, verbose=False):
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  os.path.walk(top=root_dir, func=walk_callback, arg=out)

def exercise():
  command_line = (iotbx_option_parser()
    .option(None, "--verbose",
      action="store_true",
      dest="verbose")
    .option(None, "--forever",
      action="store_true",
      dest="forever",
      help="Infinite loop, for detection of memory leaks")
    .option(None, "--walk",
      action="store",
      type="string",
      dest="walk",
      metavar="ROOT_DIR",
      help="Find and process all MTZ files under ROOT_DIR")
  ).process(args=sys.argv[1:])
  exercise_basic()
  for file_name in command_line.args:
    exercise_valid_any(file_name=file_name, out=sys.stdout)
  if (command_line.options.walk is not None):
    exercise_walk(
      root_dir=command_line.options.walk,
      verbose=command_line.options.verbose)
  while (command_line.options.forever):
    exercise_basic()
  print "OK"

if (__name__ == "__main__"):
  exercise()
