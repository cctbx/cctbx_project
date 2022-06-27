from __future__ import absolute_import, division, print_function
import libtbx.load_env
from six.moves import range
from six.moves import zip
if (libtbx.env.has_module("ccp4io")):
  from iotbx import mtz
else:
  mtz = None
from iotbx.option_parser import option_parser
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from itertools import count
from six.moves import cStringIO as StringIO
import sys, os

def exercise_read_corrupt():
  for i_trial in range(5):
    f = open("tmp_iotbx_mtz_ext.mtz", "w")
    if (i_trial > 0):
      f.write("\0"*(40*i_trial))
    f.close()
    try: mtz.object(file_name="tmp_iotbx_mtz_ext.mtz")
    except RuntimeError as e:
      assert str(e) == "cctbx Error: MTZ file read error: tmp_iotbx_mtz_ext.mtz"
    else: raise Exception_expected

def exercise_setting_nref_etc():
  m = mtz.object()
  assert m.n_reflections() == 0
  m.adjust_column_array_sizes(10)
  m.set_n_reflections(10)
  assert m.n_reflections() == 10

def exercise_basic():
  assert mtz.ccp4_liberr_verbosity(-1) == 0
  assert mtz.ccp4_liberr_verbosity(1) == 1
  assert mtz.ccp4_liberr_verbosity(-1) == 1
  assert mtz.ccp4_liberr_verbosity(0) == 0
  assert mtz.ccp4_liberr_verbosity(-1) == 0
  mtz_object = mtz.object()
  assert mtz_object.title() == ""
  assert mtz_object.history().size() == 0
  assert mtz_object.space_group_name() == ""
  assert mtz_object.space_group_number() == 0
  assert mtz_object.n_symmetry_matrices() == 0
  assert mtz_object.space_group_confidence() == "\x00"
  assert mtz_object.space_group().order_z() == 1
  assert mtz_object.point_group_name() == ""
  assert mtz_object.lattice_centring_type() == "\0"
  assert mtz_object.n_batches() == 0
  assert mtz_object.batches().size() == 0
  assert mtz_object.n_reflections() == 0
  assert mtz_object.max_min_resolution() == (-1, -1)
  assert mtz_object.n_crystals() == 0
  assert mtz_object.n_active_crystals() == 0
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/dano.mtz",
    test=os.path.isfile)
  if (file_name is None):
    print("Skipping dano.mtz test: input file not available")
  else:
    mtz_object = mtz.object(file_name=file_name)
    assert mtz_object.title() == "......"
    assert mtz_object.history().size() == 17
    assert mtz_object.space_group_name() == "P212121"
    assert mtz_object.space_group_number() == 19
    assert mtz_object.point_group_name() == "PG222"
    assert mtz_object.lattice_centring_type() == "P"
    assert mtz_object.n_symmetry_matrices() == 4
    assert mtz_object.space_group().type().lookup_symbol() == "P 21 21 21"
    assert mtz_object.n_batches() == 0
    assert mtz_object.batches().size() == 0
    assert mtz_object.n_reflections() == 165
    assert approx_equal(mtz_object.max_min_resolution(),
      (19.869975507347792, 15.001543055390009))
    assert mtz_object.n_crystals() == 4
    assert mtz_object.n_active_crystals() == 3
    assert mtz_object.has_crystal("unknown")
    assert not mtz_object.has_crystal("abc")
    assert mtz_object.has_column("H")
    assert not mtz_object.has_column("abc")
    crystal = mtz.crystal(mtz_object=mtz_object, i_crystal=1)
    assert crystal.mtz_object().n_reflections() == 165
    assert crystal.i_crystal() == 1
    assert mtz_object.crystals().size() == mtz_object.n_crystals()
    assert crystal.id() == 2
    assert crystal.name() == "unknown"
    assert crystal.set_name(new_name="abc") is crystal
    assert crystal.name() == "abc"
    assert crystal.set_name("abc") is crystal
    assert crystal.name() == "abc"
    try: crystal.set_name("unknown3")
    except RuntimeError as e:
      assert str(e) == 'mtz::crystal::set_name(new_name="unknown3"):' \
        ' new_name is used already for another crystal.'
    else: raise Exception_expected
    assert crystal.name() == "abc"
    assert crystal.set_name("unknown") is crystal
    assert crystal.name() == "unknown"
    assert crystal.project_name() == "unknown"
    assert crystal.set_project_name(new_project_name="abc") is crystal
    assert crystal.project_name() == "abc"
    assert crystal.set_project_name(new_project_name="unknown") is crystal
    assert crystal.project_name() == "unknown"
    assert approx_equal(crystal.unit_cell_parameters(),
      (84.511, 104.308, 174.103, 90, 90, 90))
    assert approx_equal(crystal.unit_cell().parameters(),
      (84.511, 104.308, 174.103, 90, 90, 90))
    assert crystal.n_datasets() == 1
    dataset = mtz.dataset(mtz_crystal=crystal, i_dataset=0)
    assert dataset.mtz_crystal().i_crystal() == 1
    assert dataset.i_dataset() == 0
    assert dataset.mtz_object().n_crystals() == mtz_object.n_crystals()
    assert dataset.id() == 1
    assert dataset.name() == "unknown230103:23:14:49"
    assert dataset.set_name(new_name="abc") is dataset
    assert dataset.name() == "abc"
    assert dataset.set_name(new_name="abc") is dataset
    assert dataset.name() == "abc"
    assert dataset.set_name("unknown230103:23:14:49") is dataset
    assert dataset.name() == "unknown230103:23:14:49"
    assert dataset.wavelength() == 0
    assert dataset.set_wavelength(new_wavelength=0.12) is dataset
    assert approx_equal(dataset.wavelength(), 0.12)
    assert dataset.set_wavelength(0) is dataset
    assert dataset.wavelength() == 0
    column = mtz.column(mtz_dataset=dataset, i_column=0)
    assert column.mtz_dataset().mtz_crystal().i_crystal() == 1
    assert column.i_column() == 0
    assert column.mtz_crystal().i_crystal() == 1
    assert column.mtz_object().n_reflections() == 165
    assert column.label() == "H"
    assert column.set_label(new_label="New") is column
    assert column.label() == "New"
    try: column.set_label("a,b,c")
    except RuntimeError as e:
      assert str(e) == 'mtz::column::set_label(new_label="a,b,c"):' \
        ' new_label must not include commas.'
    else: raise Exception_expected
    assert column.label() == "New"
    assert column.set_label(new_label="New") is column
    assert column.label() == "New"
    try: column.set_label(new_label="K")
    except RuntimeError as e:
      assert str(e) == 'mtz::column::set_label(new_label="K"):' \
        ' new_label is used already for another column.'
    else: raise Exception_expected
    assert column.set_label("H") is column
    assert column.label() == "H"
    assert column.type() == "H"
    assert column.set_type(new_type="Nw") is column
    assert column.type() == "Nw"
    assert column.set_type("H") is column
    assert column.type() == "H"
    assert column.is_active()
    if (column.source() is None):
      assert column.group_name() is None
      assert column.group_type() is None
      assert column.group_position() == -1
    else:
      assert column.source() == ""
      assert column.set_source(new_source="NsRc") is column
      assert column.source() == "NsRc"
      assert column.set_source(new_source="") is column
      assert column.source() == ""
      assert column.group_name() == ""
      assert column.set_group_name(new_group_name="NgN") is column
      assert column.group_name() == "NgN"
      assert column.set_group_name(new_group_name="") is column
      assert column.group_name() == ""
      assert column.group_type() == ""
      assert column.set_group_type(new_group_type="NgT") is column
      assert column.group_type() == "NgT"
      assert column.set_group_type(new_group_type="") is column
      assert column.group_type() == ""
      assert column.group_position() == -1
      assert column.set_group_position(new_group_position=23) is column
      assert column.group_position() == 23
      assert column.set_group_position(new_group_position=-1) is column
      assert column.group_position() == -1
    assert column.array_size() == 165
    assert column.array_capacity() == 200
    assert column.path() == "/unknown/unknown230103:23:14:49/H"
    assert column.get_other("H").i_column() == 0
    assert column.get_other("K").i_column() == 1
    assert column.get_other("L").i_column() == 2
    assert column.n_valid_values() == 165
    valid_values = column.extract_valid_values()
    assert valid_values.size() == 165
    assert approx_equal(flex.min(valid_values), 0)
    assert approx_equal(flex.max(valid_values), 5)
    assert approx_equal(flex.mean(valid_values), 2.41818189621)
    assert column.selection_valid().count(True) == 165
    assert approx_equal(flex.mean(column.extract_values()), 2.41818189621)
    assert approx_equal(flex.mean(column.extract_values(
      not_a_number_substitute=10)), 2.41818189621)
    column = mtz_object.get_column("F*")
    assert column.label() == "Frem"
    assert column.n_valid_values() == 163
    valid_values = column.extract_valid_values()
    assert valid_values.size() == 163
    assert approx_equal(flex.min(valid_values), 32.5101776123)
    assert approx_equal(flex.max(valid_values), 2711.84350586)
    assert approx_equal(flex.mean(valid_values), 615.060852051)
    assert column.selection_valid().count(True) == 163
    assert approx_equal(flex.mean(column.extract_values()),
      615.060852051*163/165)
    assert approx_equal(flex.mean(column.extract_values(13)),
      (615.060852051*163+2*13)/165)
    v,s = column.extract_values_and_selection_valid(
      not_a_number_substitute=-97).as_tuple()
    assert v.size() == 165
    assert s.count(True) == 163
    assert approx_equal(v.select(~s), [-97]*2)
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
    for i_crystal,crystal in enumerate(mtz_object.crystals()):
      assert crystal.mtz_object().n_reflections() == 165
      assert crystal.i_crystal() == i_crystal
      assert crystal.n_datasets() == 1
      for i_dataset,dataset in enumerate(crystal.datasets()):
        assert dataset.mtz_crystal().i_crystal() == i_crystal
        assert dataset.i_dataset() == i_dataset
        assert dataset.id() == next(expected_dataset_ids)
        assert dataset.name() == next(expected_dataset_names)
        assert dataset.wavelength() == 0
        assert dataset.n_columns() == next(expected_n_columns)
        for i_column,column in enumerate(dataset.columns()):
          assert column.mtz_dataset().i_dataset() == i_dataset
          assert column.i_column() == i_column
          assert column.label() == next(expected_column_labels)
          assert column.type() == next(expected_column_types)
          assert column.is_active()
          assert column.array_size() == 165
          assert column.array_capacity() == 200
          assert column.path().endswith(column.label())
          get_column = mtz_object.get_column(column.label())
          assert get_column.label() == column.label()
    group = mtz_object.extract_integers(
      column_label="ISYMabs")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 165
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    data = mtz_object.extract_integers(
      mtz_reflection_indices=group.mtz_reflection_indices,
      column_label="ISYMabs")
    assert data.all_eq(group.data)
    group = mtz_object.extract_integers_anomalous(
      column_label_plus="ISYMabs",
      column_label_minus="ISYMabs")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 330
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_reals(
      column_label="Frem")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 163
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    data = mtz_object.extract_reals(
      mtz_reflection_indices=group.mtz_reflection_indices,
      column_label="Frem")
    assert data.all_eq(group.data)
    group = mtz_object.extract_reals_anomalous(
      column_label_plus="Frem",
      column_label_minus="DANOrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 326
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_hendrickson_lattman(
      column_label_a="Frem",
      column_label_b="DANOrem",
      column_label_c="Frem",
      column_label_d="DANOrem")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 163
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_hendrickson_lattman_ab_only(
      column_label_a="Frem",
      column_label_b="DANOrem")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 163
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_hendrickson_lattman_anomalous(
      column_label_a_plus="Frem",
      column_label_b_plus="DANOrem",
      column_label_c_plus="Frem",
      column_label_d_plus="DANOrem",
      column_label_a_minus="Frem",
      column_label_b_minus="DANOrem",
      column_label_c_minus="Frem",
      column_label_d_minus="DANOrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 326
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_hendrickson_lattman_anomalous_ab_only(
      column_label_a_plus="Frem",
      column_label_b_plus="DANOrem",
      column_label_a_minus="Frem",
      column_label_b_minus="DANOrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 326
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_observations(
      column_label_data="Frem",
      column_label_sigmas="SIGFrem")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 163
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    assert group.sigmas.size() == group.indices.size()
    group = mtz_object.extract_observations_anomalous(
      column_label_data_plus="Frem",
      column_label_sigmas_plus="SIGFrem",
      column_label_data_minus="DANOrem",
      column_label_sigmas_minus="SIGDANOrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 272
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    assert group.sigmas.size() == group.indices.size()
    group = mtz_object.extract_delta_anomalous(
      column_label_f_data="Frem",
      column_label_f_sigmas="SIGFrem",
      column_label_d_data="DANOrem",
      column_label_d_sigmas="SIGDANOrem",
      column_label_isym="ISYMrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 272
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    assert group.sigmas.size() == group.indices.size()
    group = mtz_object.extract_complex(
      column_label_ampl="Frem",
      column_label_phi="SIGFrem")
    assert not group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 163
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()
    group = mtz_object.extract_complex_anomalous(
      column_label_ampl_plus="Frem",
      column_label_phi_plus="SIGFrem",
      column_label_ampl_minus="DANOrem",
      column_label_phi_minus="SIGDANOrem")
    assert group.anomalous_flag
    assert group.mtz_reflection_indices.size() == 326
    assert group.indices.size() == group.mtz_reflection_indices.size()
    assert group.data.size() == group.mtz_reflection_indices.size()

class QuickStop(Exception): pass

class exercise_extract_any(object):

  def __init__(self, full=True):
    self.full = full
    self.counters = {
      "extract_integers": 0,
      "extract_reals": 0,
      "extract_reals_anomalous": 0,
      "extract_hendrickson_lattman": 0,
      "extract_observations": 0,
      "extract_observations_anomalous": 0,
      "extract_delta_anomalous": 0,
      "extract_complex": 0,
      "extract_complex_anomalous": 0}

  def all_tests_ran_at_least_once(self):
    return min(self.counters.values()) > 0

  def raise_if_all_tests_ran_at_least_once(self):
    if (self.full): return
    if (self.all_tests_ran_at_least_once()): raise QuickStop

  def __call__(self, file_name, out):
    mtz_object = mtz.object(file_name=file_name)
    mtz_object.show_summary(out=out)
    types = "".join(mtz_object.column_types())
    for type_group in ["B", "H", "I", "Y"]:
      i = types.find(type_group)
      if (i >= 0):
        label = mtz_object.column_labels()[i]
        group = mtz_object.extract_integers(
          column_label=label)
        assert group.data.size() == group.indices.size()
        self.counters["extract_integers"] += 1
    for type_group in ["FQ", "JQ"]:
      i = types.find(type_group)
      if (i >= 0):
        labels = mtz_object.column_labels()[i:i+2]
        group = mtz_object.extract_reals(
          column_label=labels[0])
        assert group.data.size() == group.indices.size()
        self.counters["extract_reals"] += 1
        group = mtz_object.extract_observations(
          column_label_data=labels[0],
          column_label_sigmas=labels[1])
        assert group.data.size() == group.indices.size()
        assert group.sigmas.size() == group.indices.size()
        self.counters["extract_observations"] += 1
    for type_group in ["FQFQ", "JQJQ"]:
      i = types.find(type_group)
      if (i >= 0):
        labels = mtz_object.column_labels()[i:i+4]
        group = mtz_object.extract_observations_anomalous(
          column_label_data_plus=labels[0],
          column_label_sigmas_plus=labels[1],
          column_label_data_minus=labels[2],
          column_label_sigmas_minus=labels[3])
        assert group.data.size() == group.indices.size()
        assert group.sigmas.size() == group.indices.size()
        self.counters["extract_observations_anomalous"] += 1
    i = types.find("FQDQ")
    if (i >= 0):
      labels = mtz_object.column_labels()[i:i+4]
      group = mtz_object.extract_delta_anomalous(
        column_label_f_data=labels[0],
        column_label_f_sigmas=labels[1],
        column_label_d_data=labels[2],
        column_label_d_sigmas=labels[3],
        column_label_isym=mtz_object.next_isym_column_starting_at(
          i_column=i+4, return_label=True))
      assert group.data.size() == group.indices.size()
      self.counters["extract_delta_anomalous"] += 1
    i = types.find("FP")
    if (i >= 0):
      labels = mtz_object.column_labels()[i:i+2]
      group = mtz_object.extract_complex(
        column_label_ampl=labels[0],
        column_label_phi=labels[1])
      assert group.data.size() == group.indices.size()
      self.counters["extract_complex"] += 1
    for type_group in ["GLGL", "KMKM"]:
      i = types.find(type_group)
      if (i >= 0):
        labels = mtz_object.column_labels()[i:i+4]
        group = mtz_object.extract_reals_anomalous(
          column_label_plus=labels[0],
          column_label_minus=labels[2])
        assert group.data.size() == group.indices.size()
        self.counters["extract_reals_anomalous"] += 1
        group = mtz_object.extract_observations_anomalous(
          column_label_data_plus=labels[0],
          column_label_sigmas_plus=labels[1],
          column_label_data_minus=labels[2],
          column_label_sigmas_minus=labels[3])
        assert group.data.size() == group.indices.size()
        assert group.sigmas.size() == group.indices.size()
        self.counters["extract_observations_anomalous"] += 1
        # work around lack of FPFP in MTZ files available
        group = mtz_object.extract_complex_anomalous(
          column_label_ampl_plus=labels[0],
          column_label_phi_plus=labels[1],
          column_label_ampl_minus=labels[2],
          column_label_phi_minus=labels[3])
        assert group.data.size() == group.indices.size()
        self.counters["extract_complex_anomalous"] += 1
    i = types.find("AAAA")
    if (i >= 0):
      labels = mtz_object.column_labels()[i:i+4]
      group = mtz_object.extract_hendrickson_lattman(
        column_label_a=labels[0],
        column_label_b=labels[1],
        column_label_c=labels[2],
        column_label_d=labels[3])
      assert group.data.size() == group.indices.size()
      self.counters["extract_hendrickson_lattman"] += 1

def walk_callback(arg, top, names):
  exercise_function, out = arg
  for name in names:
    if (not name.lower().endswith(".mtz")): continue
    file_name = os.path.normpath(os.path.join(top, name))
    print("Processing:", file_name, file=out)
    exercise_function(file_name=file_name, out=out)
    exercise_function.raise_if_all_tests_ran_at_least_once()

def exercise_walk(root_dir, full, verbose=False):
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  exercise_function = exercise_extract_any(full=full)
  try:
    if sys.version_info.major == 3:
      for root, dirs, files in os.walk(root_dir):
        walk_callback(
            arg=(exercise_function, out),
            top=root,
            names=files)
    else:
      os.path.walk(
        top=root_dir, func=walk_callback, arg=(exercise_function, out))
  except QuickStop:
    pass
  if (verbose):
    print(exercise_function.counters)

def exercise_modifiers(verbose=0):
  if (verbose):
    out = sys.stdout
  mtz_object = mtz.object()
  mtz_object.set_title(title="012345678")
  assert mtz_object.title() == "012345678"
  mtz_object.set_title(title="012345678", append=True)
  assert mtz_object.title() == "012345678 012345678"
  mtz_object.set_title(title="0123456789"*10+"012345678", append=True)
  assert mtz_object.title() == "012345678 "*2 + "0123456789"*5
  mtz_object.set_title("0123456789"*100)
  assert mtz_object.title() == "0123456789"*7
  mtz_object.set_title(title="")
  assert mtz_object.title() == ""
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "abc"
  mtz_object.set_title(title="def", append=True)
  assert mtz_object.title() == "abc def"
  mtz_object.set_title(title="a")
  assert mtz_object.title() == "a"
  mtz_object.set_title(title="bc", append=True)
  assert mtz_object.title() == "a bc"
  mtz_object.set_title(title=" "*70)
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "abc"
  mtz_object.set_title(title="z"*70)
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "z"*70
  mtz_object.set_title(title="z"*69)
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "z"*69
  mtz_object.set_title(title="z"*68)
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "z"*68 + " a"
  mtz_object.set_title(title="z"*67)
  mtz_object.set_title(title="abc", append=True)
  assert mtz_object.title() == "z"*67 + " ab"
  mtz_object.add_history(lines=flex.std_string(["a1", "a2"]))
  assert list(mtz_object.history()) == ["a1", "a2"]
  mtz_object.add_history(lines=flex.std_string(["b1", "b2"]))
  assert list(mtz_object.history()) == ["b1", "b2", "a1", "a2"]
  mtz_object.add_history(line="c1")
  assert list(mtz_object.history()) == ["c1", "b1", "b2", "a1", "a2"]
  mtz_object.set_space_group_name(name="sg"*100)
  assert mtz_object.space_group_name() == "sgsgsgsgsgsgsgsgsgsg"
  mtz_object.set_space_group_number(number=12)
  assert mtz_object.space_group_number() == 12
  mtz_object.set_point_group_name(name="pg"*100)
  assert mtz_object.point_group_name() == "pgpgpgpgpg"
  mtz_object.set_lattice_centring_type(symbol="C")
  assert mtz_object.lattice_centring_type() == "C"
  for space_group_symbols in sgtbx.space_group_symbol_iterator():
    space_group = sgtbx.space_group(space_group_symbols)
    mtz_object.set_space_group(space_group)
    assert mtz_object.space_group() == space_group
    assert mtz_object.n_symmetry_matrices() == space_group.order_z()
  assert mtz_object.xml() is None
  assert mtz_object.unknown_headers() is None
  assert mtz_object.number_of_unknown_headers() == 0
  mtz_object = mtz.object() \
    .set_title(title="exercise") \
    .add_history(lines=flex.std_string(["h2"])) \
    .add_history(line="h1") \
    .set_space_group_name("sg") \
    .set_space_group_number(123) \
    .set_point_group_name("pg") \
    .set_space_group(sgtbx.space_group_info(number=123).group())
  assert mtz_object.title() == "exercise"
  assert list(mtz_object.history()) == ["h1", "h2"]
  for stage in [0,1]:
    for i_crystal in range(3):
      if (stage == 0):
        if (i_crystal % 2 == 0):
          crystal = mtz_object.add_crystal(
            name="crystal_%d"%i_crystal,
            project_name="project_%d"%i_crystal,
            unit_cell_parameters=(10+i_crystal,20,20,90,90,120))
        else:
          crystal = mtz_object.add_crystal(
            name="crystal_%d"%i_crystal,
            project_name="project_%d"%i_crystal,
            unit_cell=uctbx.unit_cell((10+i_crystal,20,20,90,90,120)))
      else:
        crystal = mtz_object.crystals()[i_crystal]
      assert crystal.i_crystal() == i_crystal
      assert crystal.name() == "crystal_%d"%i_crystal
      assert crystal.project_name() == "project_%d"%i_crystal
      assert approx_equal(crystal.unit_cell_parameters(),
        (10+i_crystal,20,20,90,90,120))
  if (not verbose): out = StringIO()
  assert mtz_object.show_summary(out=out) is mtz_object
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Title: exercise
Space group symbol from file: sg
Space group number from file: 123
Space group from matrices: P 4/m m m (No. 123)
Point group symbol from file: pg
Number of crystals: 3
Number of Miller indices: 0
History:
  h1
  h2
Crystal 1:
  Name: crystal_0
  Project: project_0
  Id: 1
  Unit cell: (10, 20, 20, 90, 90, 120)
  Number of datasets: 0
Crystal 2:
  Name: crystal_1
  Project: project_1
  Id: 2
  Unit cell: (11, 20, 20, 90, 90, 120)
  Number of datasets: 0
Crystal 3:
  Name: crystal_2
  Project: project_2
  Id: 3
  Unit cell: (12, 20, 20, 90, 90, 120)
  Number of datasets: 0
""")
  mtz_object.crystals()[1].set_unit_cell_parameters([13,21,23,81,82,83])
  assert approx_equal(mtz_object.crystals()[1].unit_cell_parameters(),
    [13,21,23,81,82,83])
  mtz_object.crystals()[1].set_unit_cell_parameters([11,20,20,90,90,120])
  assert approx_equal(mtz_object.crystals()[1].unit_cell_parameters(),
    [11,20,20,90,90,120])
  for stage in [0,1]:
    for i_crystal,crystal in enumerate(mtz_object.crystals()):
      for i_dataset in range(5-i_crystal):
        if (stage == 0):
          new_name = "dataset_%d" % i_dataset
          assert not crystal.has_dataset(name=new_name)
          dataset = crystal.add_dataset(name=new_name, wavelength=10-i_dataset)
          assert crystal.has_dataset(name=new_name)
        else:
          dataset = crystal.datasets()[i_dataset]
        assert dataset.name() == "dataset_%d"%i_dataset
        assert approx_equal(dataset.wavelength(), 10-i_dataset)
  if (not verbose): out = StringIO()
  mtz_object.show_summary(out=out)
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Title: exercise
Space group symbol from file: sg
Space group number from file: 123
Space group from matrices: P 4/m m m (No. 123)
Point group symbol from file: pg
Number of crystals: 3
Number of Miller indices: 0
History:
  h1
  h2
Crystal 1:
  Name: crystal_0
  Project: project_0
  Id: 1
  Unit cell: (10, 20, 20, 90, 90, 120)
  Number of datasets: 5
  Dataset 1:
    Name: dataset_0
    Id: 1
    Wavelength: 10
    Number of columns: 0
  Dataset 2:
    Name: dataset_1
    Id: 2
    Wavelength: 9
    Number of columns: 0
  Dataset 3:
    Name: dataset_2
    Id: 3
    Wavelength: 8
    Number of columns: 0
  Dataset 4:
    Name: dataset_3
    Id: 4
    Wavelength: 7
    Number of columns: 0
  Dataset 5:
    Name: dataset_4
    Id: 5
    Wavelength: 6
    Number of columns: 0
Crystal 2:
  Name: crystal_1
  Project: project_1
  Id: 2
  Unit cell: (11, 20, 20, 90, 90, 120)
  Number of datasets: 4
  Dataset 1:
    Name: dataset_0
    Id: 6
    Wavelength: 10
    Number of columns: 0
  Dataset 2:
    Name: dataset_1
    Id: 7
    Wavelength: 9
    Number of columns: 0
  Dataset 3:
    Name: dataset_2
    Id: 8
    Wavelength: 8
    Number of columns: 0
  Dataset 4:
    Name: dataset_3
    Id: 9
    Wavelength: 7
    Number of columns: 0
Crystal 3:
  Name: crystal_2
  Project: project_2
  Id: 3
  Unit cell: (12, 20, 20, 90, 90, 120)
  Number of datasets: 3
  Dataset 1:
    Name: dataset_0
    Id: 10
    Wavelength: 10
    Number of columns: 0
  Dataset 2:
    Name: dataset_1
    Id: 11
    Wavelength: 9
    Number of columns: 0
  Dataset 3:
    Name: dataset_2
    Id: 12
    Wavelength: 8
    Number of columns: 0
""")
  #
  dataset_0_0 = mtz_object.crystals()[0].datasets()[0]
  assert dataset_0_0.name() == "dataset_0"
  assert dataset_0_0.set_name(new_name="dataset_x") is dataset_0_0
  assert dataset_0_0.name() == "dataset_x"
  assert dataset_0_0.set_name(new_name="dataset_0") is dataset_0_0
  assert dataset_0_0.name() == "dataset_0"
  try: dataset_0_0.set_name(new_name="dataset_1")
  except RuntimeError as e:
    assert str(e) == 'mtz::dataset::set_name(new_name="dataset_1"):' \
      ' new_name is used already for another dataset.'
  else: raise Exception_expected
  assert dataset_0_0.name() == "dataset_0"
  #
  for stage in [0,1]:
    i_seq_iter = count()
    for i_crystal,crystal in enumerate(mtz_object.crystals()):
      for i_dataset,dataset in enumerate(crystal.datasets()):
        for i_column in range((i_crystal+i_dataset) % 3 + 1):
          i_seq = next(i_seq_iter)
          col_label = "column_%d"%i_seq
          col_type = "FB?"[(i_crystal-i_dataset+i_column) % 3]
          if (stage == 0):
            column = dataset.add_column(label=col_label, type=col_type)
          else:
            column = dataset.columns()[i_column]
          assert column.label() == col_label
          assert column.type() == col_type
  if (not verbose): out = StringIO()
  mtz_object.show_summary(out=out)
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Title: exercise
Space group symbol from file: sg
Space group number from file: 123
Space group from matrices: P 4/m m m (No. 123)
Point group symbol from file: pg
Number of crystals: 3
Number of Miller indices: 0
History:
  h1
  h2
Crystal 1:
  Name: crystal_0
  Project: project_0
  Id: 1
  Unit cell: (10, 20, 20, 90, 90, 120)
  Number of datasets: 5
  Dataset 1:
    Name: dataset_0
    Id: 1
    Wavelength: 10
    Number of columns: 1
    label    #valid %valid  min  max type
    column_0      0  0.00% None None F: amplitude
  Dataset 2:
    Name: dataset_1
    Id: 2
    Wavelength: 9
    Number of columns: 2
    label    #valid %valid  min  max type
    column_1      0  0.00% None None ?: *** UNDEFINED column type ***
    column_2      0  0.00% None None F: amplitude
  Dataset 3:
    Name: dataset_2
    Id: 3
    Wavelength: 8
    Number of columns: 3
    label    #valid %valid  min  max type
    column_3      0  0.00% None None B: BATCH number
    column_4      0  0.00% None None ?: *** UNDEFINED column type ***
    column_5      0  0.00% None None F: amplitude
  Dataset 4:
    Name: dataset_3
    Id: 4
    Wavelength: 7
    Number of columns: 1
    label    #valid %valid  min  max type
    column_6      0  0.00% None None F: amplitude
  Dataset 5:
    Name: dataset_4
    Id: 5
    Wavelength: 6
    Number of columns: 2
    label    #valid %valid  min  max type
    column_7      0  0.00% None None ?: *** UNDEFINED column type ***
    column_8      0  0.00% None None F: amplitude
Crystal 2:
  Name: crystal_1
  Project: project_1
  Id: 2
  Unit cell: (11, 20, 20, 90, 90, 120)
  Number of datasets: 4
  Dataset 1:
    Name: dataset_0
    Id: 6
    Wavelength: 10
    Number of columns: 2
    label     #valid %valid  min  max type
    column_9       0  0.00% None None B: BATCH number
    column_10      0  0.00% None None ?: *** UNDEFINED column type ***
  Dataset 2:
    Name: dataset_1
    Id: 7
    Wavelength: 9
    Number of columns: 3
    label     #valid %valid  min  max type
    column_11      0  0.00% None None F: amplitude
    column_12      0  0.00% None None B: BATCH number
    column_13      0  0.00% None None ?: *** UNDEFINED column type ***
  Dataset 3:
    Name: dataset_2
    Id: 8
    Wavelength: 8
    Number of columns: 1
    label     #valid %valid  min  max type
    column_14      0  0.00% None None ?: *** UNDEFINED column type ***
  Dataset 4:
    Name: dataset_3
    Id: 9
    Wavelength: 7
    Number of columns: 2
    label     #valid %valid  min  max type
    column_15      0  0.00% None None B: BATCH number
    column_16      0  0.00% None None ?: *** UNDEFINED column type ***
Crystal 3:
  Name: crystal_2
  Project: project_2
  Id: 3
  Unit cell: (12, 20, 20, 90, 90, 120)
  Number of datasets: 3
  Dataset 1:
    Name: dataset_0
    Id: 10
    Wavelength: 10
    Number of columns: 3
    label     #valid %valid  min  max type
    column_17      0  0.00% None None ?: *** UNDEFINED column type ***
    column_18      0  0.00% None None F: amplitude
    column_19      0  0.00% None None B: BATCH number
  Dataset 2:
    Name: dataset_1
    Id: 11
    Wavelength: 9
    Number of columns: 1
    label     #valid %valid  min  max type
    column_20      0  0.00% None None B: BATCH number
  Dataset 3:
    Name: dataset_2
    Id: 12
    Wavelength: 8
    Number of columns: 2
    label     #valid %valid  min  max type
    column_21      0  0.00% None None F: amplitude
    column_22      0  0.00% None None B: BATCH number
""")
  for column in mtz_object.columns():
    assert column.array_size() == 2000
    assert column.array_capacity() == 2402
  mtz_object.reserve(5000)
  for column in mtz_object.columns():
    assert column.array_size() == 2000
    assert column.array_capacity() == 5000
  mtz_object.reserve(100)
  for column in mtz_object.columns():
    assert column.array_size() == 2000
    assert column.array_capacity() == 5000
  #
  mtz_object = mtz.object() \
    .set_title(title="exercise") \
    .set_space_group_name("sg") \
    .set_space_group_number(123) \
    .set_point_group_name("pg") \
    .set_lattice_centring_type("pg") \
    .set_space_group(sgtbx.space_group_info(number=123).group())
  unit_cell = uctbx.unit_cell((10,10,10,90,90,90))
  mtz_object.set_hkl_base(unit_cell=unit_cell)
  dataset = mtz_object.add_crystal(
    name="crystal_1",
    project_name="crystal_1",
    unit_cell=unit_cell).add_dataset(
      name="crystal_1",
      wavelength=0)
  try: dataset.add_column(label="a,b,c", type="H")
  except RuntimeError as e:
    assert str(e) == 'mtz::dataset::add_column(label="a,b,c", ...):' \
      ' label must not include commas.'
  else: raise Exception_expected
  for label in "HKL":
    dataset.add_column(label=label, type="H")
  column = dataset.add_column(label="F", type="F")
  mtz_reflection_indices = column.set_reals(
    miller_indices=flex.miller_index([(1,2,3),(2,3,4),(3,4,5)]),
    data=flex.double([10,20,30]))
  assert list(mtz_reflection_indices) == [0,1,2]
  column = dataset.add_column(label="SigF", type="Q")
  column.set_reals(
    mtz_reflection_indices=mtz_reflection_indices,
    data=flex.double([1,2,3]))
  group = mtz_object.extract_observations(
    column_label_data="F",
    column_label_sigmas="SigF")
  assert list(group.indices) == [(1, 2, 3), (2, 3, 4), (3, 4, 5)]
  assert approx_equal(group.data, [10, 20, 30])
  assert approx_equal(group.sigmas, [1, 2, 3])
  column = dataset.add_column(label="I", type="F")
  mtz_reflection_indices = column.set_reals(
    miller_indices=flex.miller_index([(2,3,5),(1,2,3),(3,4,5)]),
    data=flex.double([11,21,31]))
  assert list(mtz_reflection_indices) == [3, 0, 2]
  column = dataset.add_column(label="SigI", type="Q")
  column.set_reals(
    mtz_reflection_indices=mtz_reflection_indices,
    data=flex.double([4,5,6]))
  group = mtz_object.extract_observations(
    column_label_data="I",
    column_label_sigmas="SigI")
  assert list(group.indices) == [(1, 2, 3), (3, 4, 5), (2, 3, 5)]
  assert approx_equal(group.data, [21, 31, 11])
  assert approx_equal(group.sigmas, [5, 6, 4])
  if (not verbose): out = StringIO()
  mtz_object.show_summary(out=out)
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Title: exercise
Space group symbol from file: sg
Space group number from file: 123
Space group from matrices: P 4/m m m (No. 123)
Point group symbol from file: pg
Number of crystals: 2
Number of Miller indices: 4
Resolution range: 2.67261 1.41421
History:
Crystal 1:
  Name: HKL_base
  Project: HKL_base
  Id: 0
  Unit cell: (10, 10, 10, 90, 90, 90)
  Number of datasets: 1
  Dataset 1:
    Name: HKL_base
    Id: 0
    Wavelength: 0
    Number of columns: 0
Crystal 2:
  Name: crystal_1
  Project: crystal_1
  Id: 1
  Unit cell: (10, 10, 10, 90, 90, 90)
  Number of datasets: 1
  Dataset 1:
    Name: crystal_1
    Id: 1
    Wavelength: 0
    Number of columns: 7
    label #valid  %valid   min   max type
    H          4 100.00%  1.00  3.00 H: index h,k,l
    K          4 100.00%  2.00  4.00 H: index h,k,l
    L          4 100.00%  3.00  5.00 H: index h,k,l
    F          3  75.00% 10.00 30.00 F: amplitude
    SigF       3  75.00%  1.00  3.00 Q: standard deviation
    I          3  75.00% 11.00 31.00 F: amplitude
    SigI       3  75.00%  4.00  6.00 Q: standard deviation
""")
  if (not verbose): out = StringIO()
  assert mtz_object.show_column_data(out=out) is mtz_object
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Column data:
-------------------------------------------------------------------------------
                       F            SigF               I            SigI

 1  2  3              10               1              21               5
 2  3  4              20               2            None            None
 3  4  5              30               3              31               6
 2  3  5            None            None              11               4
-------------------------------------------------------------------------------
""")
  mtz_object.write(file_name="tmp_iotbx_mtz_ext.mtz")
  if (not verbose): out = StringIO()
  mtz.object(file_name="tmp_iotbx_mtz_ext.mtz").show_summary(out=out)
  if (not verbose):
    assert not show_diff(out.getvalue(), """\
Title: exercise
Space group symbol from file: sg
Space group number from file: 123
Space group from matrices: P 4/m m m (No. 123)
Point group symbol from file: pg
Number of crystals: 2
Number of Miller indices: 4
Resolution range: 2.67261 1.41421
History:
Crystal 1:
  Name: HKL_base
  Project: HKL_base
  Id: 0
  Unit cell: (10, 10, 10, 90, 90, 90)
  Number of datasets: 1
  Dataset 1:
    Name: HKL_base
    Id: 0
    Wavelength: 0
    Number of columns: 0
Crystal 2:
  Name: crystal_1
  Project: crystal_1
  Id: 2
  Unit cell: (10, 10, 10, 90, 90, 90)
  Number of datasets: 1
  Dataset 1:
    Name: crystal_1
    Id: 1
    Wavelength: 0
    Number of columns: 7
    label #valid  %valid   min   max type
    H          4 100.00%  1.00  3.00 H: index h,k,l
    K          4 100.00%  2.00  4.00 H: index h,k,l
    L          4 100.00%  3.00  5.00 H: index h,k,l
    F          3  75.00% 10.00 30.00 F: amplitude
    SigF       3  75.00%  1.00  3.00 Q: standard deviation
    I          3  75.00% 11.00 31.00 F: amplitude
    SigI       3  75.00%  4.00  6.00 Q: standard deviation
""")
  #
  original_miller_indices = mtz_object.extract_miller_indices()
  assert list(original_miller_indices) \
      == [(1, 2, 3), (2, 3, 4), (3, 4, 5), (2, 3, 5)]
  new_miller_indices = flex.miller_index(
    [(3, -1, 2), (-4, 2, -3), (5, -3, 4), (-5, 2, -3)])
  assert not mtz_object.extract_miller_indices().all_eq(new_miller_indices)
  mtz_object.replace_miller_indices(miller_indices=new_miller_indices)
  assert mtz_object.extract_miller_indices().all_eq(new_miller_indices)
  mtz_object.replace_miller_indices(miller_indices=original_miller_indices)
  assert not mtz_object.extract_miller_indices().all_eq(new_miller_indices)
  assert mtz_object.extract_miller_indices().all_eq(original_miller_indices)
  #
  c = mtz_object.get_column(label="F")
  s = c.selection_valid()
  v = c.extract_values(not_a_number_substitute=-1)
  assert list(s) == [True, True, True, False]
  assert approx_equal(v, [10.0, 20.0, 30.0, -1.0])
  c.set_values(values=flex.float([5,9,3,7]), selection_valid=None)
  v = c.extract_values(not_a_number_substitute=-1)
  assert approx_equal(v, [5.0, 9.0, 3.0, 7.0])
  c.set_values(values=flex.float([7,8,2,0]))
  v = c.extract_values(not_a_number_substitute=-1)
  assert approx_equal(v, [7.0, 8.0, 2.0, 0.0])
  c.set_values(
    values=flex.float([5,9,3,7]),
    selection_valid=flex.bool([False]*4))
  v = c.extract_values(not_a_number_substitute=-1)
  assert approx_equal(v, [-1]*4)
  for i_trial in range(10):
    s = flex.random_bool(size=4, threshold=0.5)
    v = flex.float(list(flex.random_double(size=4)*10-5))
    c.set_values(values=v, selection_valid=s)
    sx = c.selection_valid()
    vx = c.extract_values(not_a_number_substitute=99)
    assert list(s) == list(sx)
    assert vx.select(s).all_eq(v.select(s))
    assert vx.select(~s).all_ne(v.select(~s))
    assert vx.select(~s).all_eq(99)
  #
  values_in = count()
  values_out = count()
  for i_batch in range(10):
    batch = mtz_object.add_batch()
    assert batch.num() == i_batch+1
    assert batch.set_num(value=next(values_in)) is batch
    assert batch.num() == next(values_out)
    assert batch.set_num(value=i_batch+1) is batch
    assert batch.title() == " "
    assert batch.set_title("Hello MTZ") is batch
    assert batch.title() == "Hello MTZ"
    assert batch.set_title("Hello MTZ"*10) is batch
    assert len(batch.title()) == 70
    assert list(batch.gonlab()) == ["", "", ""]
    assert batch.set_gonlab(
      flex.std_string(["what", "ever", "this_is....."])) is batch
    assert list(batch.gonlab()) == ["what", "ever", "this_is"]
    assert batch.iortyp() == 0
    assert batch.set_iortyp(value=next(values_in)) is batch
    assert batch.iortyp() == next(values_out)
    assert list(batch.lbcell()) == [0, 0, 0, 0, 0, 0]
    assert batch.set_lbcell(flex.int(range(3,9))) is batch
    assert list(batch.lbcell()) == list(range(3,9))
    assert batch.misflg() == 0
    assert batch.set_misflg(value=next(values_in)) is batch
    assert batch.misflg() == next(values_out)
    assert batch.jumpax() == 0
    assert batch.set_jumpax(value=next(values_in)) is batch
    assert batch.jumpax() == next(values_out)
    assert batch.ncryst() == 0
    assert batch.set_ncryst(value=next(values_in)) is batch
    assert batch.ncryst() == next(values_out)
    assert batch.lcrflg() == 0
    assert batch.set_lcrflg(value=next(values_in)) is batch
    assert batch.lcrflg() == next(values_out)
    assert batch.ldtype() == 0
    assert batch.set_ldtype(value=next(values_in)) is batch
    assert batch.ldtype() == next(values_out)
    assert batch.jsaxs() == 0
    assert batch.set_jsaxs(value=next(values_in)) is batch
    assert batch.jsaxs() == next(values_out)
    assert batch.nbscal() == 0
    assert batch.set_nbscal(value=next(values_in)) is batch
    assert batch.nbscal() == next(values_out)
    assert batch.ngonax() == 0
    assert batch.set_ngonax(value=next(values_in)) is batch
    assert batch.ngonax() == next(values_out)
    assert batch.lbmflg() == 0
    assert batch.set_lbmflg(value=next(values_in)) is batch
    assert batch.lbmflg() == next(values_out)
    assert batch.ndet() == 0
    assert batch.set_ndet(value=next(values_in) % 3) is batch
    assert batch.ndet() == next(values_out) % 3
    assert batch.nbsetid() == 0
    assert batch.set_nbsetid(value=next(values_in)) is batch
    assert batch.nbsetid() == next(values_out)
    assert list(batch.cell()) == [0]*6
    assert batch.set_cell(flex.float(range(18,24))) is batch
    assert list(batch.cell()) == list(range(18,24))
    assert list(batch.umat()) == [0]*9
    assert batch.set_umat(flex.float(range(16,25))) is batch
    assert list(batch.umat()) == list(range(16,25))
    assert list(batch.phixyz()) == [0]*6
    assert batch.set_phixyz(flex.float(range(28,34))) is batch
    assert list(batch.phixyz()) == list(range(28,34))
    assert list(batch.crydat()) == [0]*12
    assert batch.set_crydat(flex.float(range(26,38))) is batch
    assert list(batch.crydat()) == list(range(26,38))
    assert list(batch.datum()) == [0]*3
    assert batch.set_datum(flex.float(range(26,29))) is batch
    assert list(batch.datum()) == list(range(26,29))
    assert batch.phistt() == 0
    assert batch.set_phistt(value=next(values_in)) is batch
    assert batch.phistt() == next(values_out)
    assert batch.phiend() == 0
    assert batch.set_phiend(value=next(values_in)) is batch
    assert batch.phiend() == next(values_out)
    assert list(batch.scanax()) == [0]*3
    assert batch.set_scanax(flex.float(range(62,65))) is batch
    assert list(batch.scanax()) == list(range(62,65))
    assert batch.time1() == 0
    assert batch.set_time1(value=next(values_in)) is batch
    assert batch.time1() == next(values_out)
    assert batch.time2() == 0
    assert batch.set_time2(value=next(values_in)) is batch
    assert batch.time2() == next(values_out)
    assert batch.bscale() == 0
    assert batch.set_bscale(value=next(values_in)) is batch
    assert batch.bscale() == next(values_out)
    assert batch.bbfac() == 0
    assert batch.set_bbfac(value=next(values_in)) is batch
    assert batch.bbfac() == next(values_out)
    assert batch.sdbscale() == 0
    assert batch.set_sdbscale(value=next(values_in)) is batch
    assert batch.sdbscale() == next(values_out)
    assert batch.sdbfac() == 0
    assert batch.set_sdbfac(value=next(values_in)) is batch
    assert batch.sdbfac() == next(values_out)
    assert batch.phirange() == 0
    assert batch.set_phirange(value=next(values_in)) is batch
    assert batch.phirange() == next(values_out)
    assert list(batch.e1()) == [0]*3
    assert batch.set_e1(flex.float(range(71,74))) is batch
    assert list(batch.e1()) == list(range(71,74))
    assert list(batch.e2()) == [0]*3
    assert batch.set_e2(flex.float(range(72,75))) is batch
    assert list(batch.e2()) == list(range(72,75))
    assert list(batch.e3()) == [0]*3
    assert batch.set_e3(flex.float(range(73,76))) is batch
    assert list(batch.e3()) == list(range(73,76))
    assert list(batch.source()) == [0]*3
    assert batch.set_source(flex.float(range(74,77))) is batch
    assert list(batch.source()) == list(range(74,77))
    assert list(batch.so()) == [0]*3
    assert batch.set_so(flex.float(range(75,78))) is batch
    assert list(batch.so()) == list(range(75,78))
    assert batch.alambd() == 0
    assert batch.set_alambd(value=next(values_in)) is batch
    assert batch.alambd() == next(values_out)
    assert batch.delamb() == 0
    assert batch.set_delamb(value=next(values_in)) is batch
    assert batch.delamb() == next(values_out)
    assert batch.delcor() == 0
    assert batch.set_delcor(value=next(values_in)) is batch
    assert batch.delcor() == next(values_out)
    assert batch.divhd() == 0
    assert batch.set_divhd(value=next(values_in)) is batch
    assert batch.divhd() == next(values_out)
    assert batch.divvd() == 0
    assert batch.set_divvd(value=next(values_in)) is batch
    assert batch.divvd() == next(values_out)
    assert list(batch.dx()) == [0]*2
    assert batch.set_dx(flex.float(range(84,86))) is batch
    assert list(batch.dx()) == list(range(84,86))
    assert list(batch.theta()) == [0]*2
    assert batch.set_theta(flex.float(range(85,87))) is batch
    assert list(batch.theta()) == list(range(85,87))
    assert list(batch.detlm()) == [0]*8
    assert batch.set_detlm(flex.float(range(86,94))) is batch
    assert list(batch.detlm()) == list(range(86,94))
    if (not verbose): out = StringIO()
    batch.show(out=out)
    if (not verbose and i_batch == 3):
      batch_3_show = out
      assert not show_diff(out.getvalue(), """\
batch number: 4
batch title: Hello MTZHello MTZHello MTZHello MTZHello MTZHello MTZHello MTZHello M
names of the three axes: ['what', 'ever', 'this_is']
type of orientation block: 82
refinement flags for cell: [3, 4, 5, 6, 7, 8]
number of phixyz used (0, 1, or 2): 83
reciprocal axis closest to rotation axis: 84
crystal number: 85
mosaicity model: 0 = isotropic, 1 = anisotropic: 86
type of data: 2D (1), 3D (2), or Laue (3): 87
goniostat scan axis number: 88
number of batch scales & Bfactors (0 if unset): 89
number of goniostat axes: 90
flag for type of beam info: 91
  0: for alambd, delamb; 1: also delcor, divhd, divvd
number of detectors (current maximum 2): 2
dataset id: 93
cell dimensions: [18.0, 19.0, 20.0, 21.0, 22.0, 23.0]
orientation matrix U: [16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0]
  in Fortranic order, i.e. U(1,1), U(2,1) ...
missetting angles at beginning and end of oscillation: [28.0, 29.0, 30.0, 31.0, 32.0, 33.0]
mosaicity: [26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0]
datum values of goniostat axes: [26.0, 27.0, 28.0]
start of phi relative to datum: 94.0
end of phi relative to datum: 95.0
rotation axis in lab frame: [62.0, 63.0, 64.0]
start time: 96.0
stop time: 97.0
batch scale: 98.0
batch temperature factor: 99.0
sd bscale: 100.0
sd bbfac: 101.0
phi range: 102.0
vectors ("Cambridge" laboratory axes) defining ngonax goniostat axes:
  vector 1: [71.0, 72.0, 73.0]
  vector 2: [72.0, 73.0, 74.0]
  vector 3: [73.0, 74.0, 75.0]
idealised source vector: [74.0, 75.0, 76.0]
source vector: [75.0, 76.0, 77.0]
wavelength (A): 103.0
dispersion (deltalambda / lambda): 104.0
correlated component: 105.0
horizontal beam divergence: 106.0
vertical beam divergence: 107.0
xtal to detector distance: [84.0, 85.0]
detector tilt angle: [85.0, 86.0]
min & max values of detector coords (pixels): [86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0]
""")
  mtz_object.write(file_name="tmp_iotbx_mtz_ext.mtz")
  restored = mtz.object(file_name="tmp_iotbx_mtz_ext.mtz")
  assert restored.n_batches() == 10
  if (not verbose): out = StringIO()
  restored.batches()[3].show(out=out)
  if (not verbose):
    assert out.getvalue() == batch_3_show.getvalue()
  for i_trial in range(10):
    perm = flex.random_permutation(size=10)
    for batch,new_num in zip(mtz_object.batches(), perm):
      batch.set_num(value=new_num+1)
    mtz_object.sort_batches()
    assert [batch.num() for batch in mtz_object.batches()] \
        == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  for batch,nbsetid in zip(mtz_object.batches(), [0,0,1,0,1,1,1,0,1,1]):
    batch.set_nbsetid(value=nbsetid)
  for crystal in mtz_object.crystals():
    for dataset in crystal.datasets():
      assert dataset.id() in [0,1]
      assert dataset.n_batches() == [4,6][dataset.id()]
      assert dataset.batches().size() == dataset.n_batches()
      for batch in dataset.batches():
        assert batch.nbsetid() == dataset.id()
      batch = dataset.add_batch()
      assert batch.nbsetid() == dataset.id()
      assert dataset.n_batches() == [5,7][dataset.id()]
  # quick test for delete_reflection
  assert mtz_object.n_reflections() > 3
  mx = mtz_object.extract_miller_indices()[1]
  assert mx in mtz_object.extract_miller_indices()
  mtz_object.delete_reflection(1)
  assert mx not in mtz_object.extract_miller_indices()
  # test for delete_reflections
  isel = flex.size_t((1,0))
  try: mtz_object.delete_reflections(isel)
  except RuntimeError: pass
  else: raise Exception_expected
  isel = flex.size_t((0,2))
  mx = [mtz_object.extract_miller_indices()[i] for i in isel]
  for m in mx:
    assert m in mtz_object.extract_miller_indices()
  mtz_object.delete_reflections(isel)
  for m in mx:
    assert m not in mtz_object.extract_miller_indices()


def exercise():
  if (mtz is None):
    print("Skipping iotbx/mtz/tst_ext.py: ccp4io not available")
    return
  command_line = (option_parser()
    .option(None, "--verbose",
      action="store_true")
    .option(None, "--forever",
      action="store_true",
      help="Infinite loop, for detection of memory leaks")
    .option(None, "--walk",
      action="store",
      type="string",
      metavar="ROOT_DIR",
      help="Find and process all MTZ files under ROOT_DIR")
    .option(None, "--full",
      action="store_true",
      help="Visit all MTZ files")
  ).process(args=sys.argv[1:])
  exercise_read_corrupt()
  exercise_basic()
  exercise_setting_nref_etc()
  exercise_modifiers(verbose=command_line.options.verbose)
  for file_name in command_line.args:
    exercise_extract_any()(file_name=file_name, out=sys.stdout)
  if (command_line.options.walk is not None):
    exercise_walk(
      root_dir=command_line.options.walk,
      full=command_line.options.full,
      verbose=command_line.options.verbose)
  while (command_line.options.forever):
    exercise_basic()
    exercise_modifiers()

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
