from iotbx import mtz
import iotbx.mtz.wrapper
from cctbx.development import debug_utils
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.regression.tst_miller import generate_random_hl
from iotbx.regression.utils import random_f_calc
from libtbx.test_utils import eps_eq
import sys

def to_mtz(miller_array, root_label, column_types=None):
  mtz_object = mtz.wrapper.object()
  mtz_object.set_title("mtz writer test")
  mtz_object.add_history(line="start")
  mtz_object.set_space_group_info(miller_array.space_group_info())
  mtz_object.set_hkl_base(miller_array.unit_cell())
  crystal = mtz_object.add_crystal(
    name="test_crystal",
    project_name="test_project",
    unit_cell=miller_array.unit_cell())
  dataset = crystal.add_dataset(
    name="test_dataset",
    wavelength=1)
  dataset.add_miller_array(
    miller_array=miller_array,
    root_label=root_label,
    column_types=column_types)
  mtz_object.add_history(line="done")
  return dataset

def recycle(miller_array, root_label, column_types=None, verbose=0):
  original_dataset = to_mtz(miller_array, root_label, column_types)
  written = original_dataset.mtz_object()
  if (0 or verbose):
    written.show_summary()
  original_dataset.mtz_object().write(file_name="tmp.mtz")
  restored = mtz.wrapper.object(file_name="tmp.mtz")
  if (0 or verbose):
    restored.show_summary()
  assert restored.title() == written.title()
  assert [line.rstrip() for line in restored.history()] \
      == list(written.history())
  assert restored.space_group_name() == written.space_group_name()
  assert restored.space_group_number() == written.space_group_number()
  assert restored.space_group() == written.space_group()
  assert restored.point_group_name() == written.point_group_name()
  assert restored.lattice_centring_type() == written.lattice_centring_type()
  assert restored.n_batches() == written.n_batches()
  assert restored.n_reflections() == written.n_reflections()
  assert eps_eq(
    restored.max_min_resolution(), written.max_min_resolution(), eps=1.e-5)
  assert restored.n_crystals() == written.n_crystals()
  assert restored.n_active_crystals() == written.n_active_crystals()
  assert restored.n_crystals() == 2
  for rx,wx in zip(restored.crystals(), written.crystals()):
    assert rx.name() == wx.name()
    assert rx.project_name() == wx.project_name()
    assert rx.unit_cell().is_similar_to(wx.unit_cell())
    assert rx.n_datasets() == wx.n_datasets()
    for rd,wd in zip(rx.datasets(), wx.datasets()):
      assert rd.name() == wd.name()
      assert rd.wavelength() == wd.wavelength()
      assert rd.n_columns() == wd.n_columns()
  crystal_symmetry = restored.crystals()[1].crystal_symmetry()
  restored_dataset = restored.crystals()[1].datasets()[0]
  if (not miller_array.anomalous_flag()):
    if (miller_array.sigmas() is None):
      if (miller_array.is_complex()):
        assert restored_dataset.n_columns() == 3+2
        group = restored.extract_complex(
          column_label_ampl=root_label,
          column_label_phi=original_dataset.label_phases(root_label))
      elif (miller_array.is_hendrickson_lattman_array()):
        assert restored_dataset.n_columns() == 3+4
        group = restored.extract_hls(
          column_label_a=root_label+"A",
          column_label_b=root_label+"B",
          column_label_c=root_label+"C",
          column_label_d=root_label+"D")
      else:
        assert restored_dataset.n_columns() == 3+1
        group = restored.extract_reals(
          column_label=root_label)
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=False),
        data=group.data)
    else:
      assert restored_dataset.n_columns() == 3+2
      group = restored.extract_observations(
        column_label_data=root_label,
        column_label_sigmas=original_dataset.label_sigmas(root_label))
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=False),
        data=group.data,
        sigmas=group.sigmas)
  else:
    if (miller_array.sigmas() is None):
      if (miller_array.is_complex()):
        assert restored_dataset.n_columns() == 3+4
        group = restored.extract_complex_anomalous(
          column_label_ampl_plus=original_dataset.label_plus(root_label),
          column_label_phi_plus=original_dataset.label_plus(
            original_dataset.label_phases(root_label)),
          column_label_ampl_minus=original_dataset.label_minus(root_label),
          column_label_phi_minus=original_dataset.label_minus(
            original_dataset.label_phases(root_label)))
      elif (miller_array.is_hendrickson_lattman_array()):
        assert restored_dataset.n_columns() == 3+8
        group = restored.extract_hls_anomalous(
          column_label_a_plus=root_label+"A(+)",
          column_label_b_plus=root_label+"B(+)",
          column_label_c_plus=root_label+"C(+)",
          column_label_d_plus=root_label+"D(+)",
          column_label_a_minus=root_label+"A(-)",
          column_label_b_minus=root_label+"B(-)",
          column_label_c_minus=root_label+"C(-)",
          column_label_d_minus=root_label+"D(-)")
      else:
        assert restored_dataset.n_columns() == 3+2
        group = restored.extract_reals_anomalous(
          column_label_plus=original_dataset.label_plus(root_label),
          column_label_minus=original_dataset.label_minus(root_label))
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=True),
        data=group.data)
    else:
      assert restored_dataset.n_columns() == 3+4
      group = restored.extract_observations_anomalous(
        column_label_data_plus=original_dataset.label_plus(root_label),
        column_label_sigmas_plus=original_dataset.label_plus(
          original_dataset.label_sigmas(root_label)),
        column_label_data_minus=original_dataset.label_minus(root_label),
        column_label_sigmas_minus=original_dataset.label_minus(
          original_dataset.label_sigmas(root_label)))
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=True),
        data=group.data,
        sigmas=group.sigmas)
  verify_miller_arrays(miller_array, r)
  restored_miller_arrays = restored.as_miller_arrays()
  assert len(restored_miller_arrays) == 1
  verify_miller_arrays(miller_array, restored_miller_arrays[0])

def verify_miller_arrays(a1, a2, eps=1.e-5):
  v = a2.adopt_set(a1)
  if (a1.is_bool_array()):
    if (a2.is_integer_array()):
      assert flex.max(flex.abs(a1.data().as_int() - v.data())) == 0
    else:
      assert flex.max(flex.abs(a1.data().as_double() - v.data())) < eps
  elif (a1.is_hendrickson_lattman_array()):
    for i in xrange(4):
      assert flex.max(flex.abs(a1.data().slice(i) - v.data().slice(i))) < eps
  else:
    assert flex.max(flex.abs(a1.data() - v.data())) < eps
  if (v.sigmas() is not None):
    assert flex.max(flex.abs(a1.sigmas() - v.sigmas())) < eps

def exercise(space_group_info, anomalous_flag,
             n_scatterers=8, d_min=2.5, verbose=0):
  f_calc = random_f_calc(
    space_group_info=space_group_info,
    n_scatterers=n_scatterers,
    d_min=d_min,
    anomalous_flag=anomalous_flag,
    verbose=verbose)
  if (f_calc is None): return
  recycle(f_calc, "f_calc", verbose=verbose)
  recycle(abs(f_calc), "f_obs", verbose=verbose)
  if (not anomalous_flag):
    recycle(abs(f_calc), "f_obs", column_types="R", verbose=verbose)
  recycle(miller.array(
    miller_set=f_calc,
    data=flex.abs(f_calc.data()),
    sigmas=flex.abs(f_calc.data())/10), "f_obs", verbose=verbose)
  recycle(f_calc.centric_flags(), "cent", verbose=verbose)
  recycle(generate_random_hl(miller_set=f_calc), "prob", verbose=verbose)

def run_call_back(flags, space_group_info):
  for anomalous_flag in [False, True]:
    exercise(
      space_group_info=space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
