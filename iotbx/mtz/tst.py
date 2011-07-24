import libtbx.load_env
if (libtbx.env.has_module("ccp4io")):
  from iotbx import mtz
else:
  mtz = None
from cctbx.development import debug_utils
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.regression.tst_miller import generate_random_hl
from iotbx.regression.utils import random_f_calc
from libtbx.test_utils import Exception_expected, approx_equal, eps_eq
import math
import sys
import os
op = os.path

def to_mtz(miller_array, column_root_label, column_types=None):
  mtz_object = mtz.object()
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
  assert dataset.add_miller_array(
    miller_array=miller_array,
    column_root_label=column_root_label,
    column_types=column_types) is dataset
  mtz_object.add_history(line="done")
  return dataset

def recycle(miller_array, column_root_label, column_types=None, verbose=0):
  original_dataset = to_mtz(miller_array, column_root_label, column_types)
  label_decorator = mtz.label_decorator()
  written = original_dataset.mtz_object()
  if (0 or verbose):
    written.show_summary()
  original_dataset.mtz_object().write(file_name="tmp.mtz")
  restored = mtz.object(file_name="tmp.mtz")
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
  miller_set = restored.crystals()[1].miller_set()
  assert miller_set.indices().size() == restored.n_reflections()
  crystal_symmetry = restored.crystals()[1].crystal_symmetry()
  restored_dataset = restored.crystals()[1].datasets()[0]
  if (not miller_array.anomalous_flag()):
    if (miller_array.sigmas() is None):
      if (miller_array.is_complex_array()):
        assert restored_dataset.n_columns() == 3+2
        group = restored.extract_complex(
          column_label_ampl=column_root_label,
          column_label_phi=label_decorator.phases(column_root_label))
      elif (miller_array.is_hendrickson_lattman_array()):
        assert restored_dataset.n_columns() == 3+4
        deco = label_decorator.hendrickson_lattman
        group = restored.extract_hendrickson_lattman(
          column_label_a=deco(column_root_label, 0),
          column_label_b=deco(column_root_label, 1),
          column_label_c=deco(column_root_label, 2),
          column_label_d=deco(column_root_label, 3))
      else:
        assert restored_dataset.n_columns() == 3+1
        group = restored.extract_reals(
          column_label=column_root_label)
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=False),
        data=group.data)
    else:
      assert restored_dataset.n_columns() == 3+2
      group = restored.extract_observations(
        column_label_data=column_root_label,
        column_label_sigmas=label_decorator.sigmas(column_root_label))
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=False),
        data=group.data,
        sigmas=group.sigmas)
  else:
    if (miller_array.sigmas() is None):
      if (miller_array.is_complex_array()):
        assert restored_dataset.n_columns() == 3+4
        group = restored.extract_complex_anomalous(
          column_label_ampl_plus=label_decorator.anomalous(
            column_root_label, "+"),
          column_label_phi_plus=label_decorator.phases(
            column_root_label, "+"),
          column_label_ampl_minus=label_decorator.anomalous(
            column_root_label, "-"),
          column_label_phi_minus=label_decorator.phases(
            column_root_label, "-"))
      elif (miller_array.is_hendrickson_lattman_array()):
        assert restored_dataset.n_columns() == 3+8
        deco = label_decorator.hendrickson_lattman
        group = restored.extract_hendrickson_lattman_anomalous(
          column_label_a_plus=deco(column_root_label, 0, "+"),
          column_label_b_plus=deco(column_root_label, 1, "+"),
          column_label_c_plus=deco(column_root_label, 2, "+"),
          column_label_d_plus=deco(column_root_label, 3, "+"),
          column_label_a_minus=deco(column_root_label, 0, "-"),
          column_label_b_minus=deco(column_root_label, 1, "-"),
          column_label_c_minus=deco(column_root_label, 2, "-"),
          column_label_d_minus=deco(column_root_label, 3, "-"))
      else:
        assert restored_dataset.n_columns() == 3+2
        group = restored.extract_reals_anomalous(
          column_label_plus=label_decorator.anomalous(column_root_label, "+"),
          column_label_minus=label_decorator.anomalous(column_root_label, "-"))
      r = miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=group.indices,
          anomalous_flag=True),
        data=group.data)
    else:
      assert restored_dataset.n_columns() == 3+4
      group = restored.extract_observations_anomalous(
        column_label_data_plus=label_decorator.anomalous(
          column_root_label, "+"),
        column_label_sigmas_plus=label_decorator.sigmas(
          column_root_label, "+"),
        column_label_data_minus=label_decorator.anomalous(
          column_root_label, "-"),
        column_label_sigmas_minus=label_decorator.sigmas(
          column_root_label, "-"))
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
  thff = restored_miller_arrays[0].info().type_hints_from_file
  assert thff is not None
  assert miller_array.is_hendrickson_lattman_array() \
      == (thff == "hendrickson_lattman")
  verify_miller_arrays(miller_array, restored_miller_arrays[0])
  mtz_object = miller_array.as_mtz_dataset(
    column_root_label=column_root_label).mtz_object()
  restored_miller_arrays = mtz_object.as_miller_arrays()
  assert len(restored_miller_arrays) == 1
  verify_miller_arrays(miller_array, restored_miller_arrays[0])
  if (   miller_array.is_bool_array()
      or miller_array.is_integer_array()
      or miller_array.is_real_array()):
    cb_op = miller_array.change_of_basis_op_to_niggli_cell()
    mtz_object.change_basis_in_place(cb_op=cb_op)
    cb_array = miller_array.change_basis(cb_op=cb_op)
    assert mtz_object.space_group() == cb_array.space_group()
    for mtz_crystal in mtz_object.crystals():
      assert mtz_crystal.unit_cell().is_similar_to(cb_array.unit_cell())
    restored_miller_arrays = mtz_object.as_miller_arrays()
    assert len(restored_miller_arrays) == 1
    verify_miller_arrays(cb_array, restored_miller_arrays[0])
    mtz_object.change_basis_in_place(cb_op=cb_op.inverse())
    assert mtz_object.space_group() == miller_array.space_group()
    for mtz_crystal in mtz_object.crystals():
      assert mtz_crystal.unit_cell().is_similar_to(miller_array.unit_cell())
    restored_miller_arrays = mtz_object.as_miller_arrays()
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

def exercise_recycle(
      space_group_info, anomalous_flag,
      n_scatterers=8, d_min=2.5, verbose=0):
  f_calc = random_f_calc(
    space_group_info=space_group_info,
    n_scatterers=n_scatterers,
    d_min=d_min,
    anomalous_flag=anomalous_flag,
    verbose=verbose)
  if (f_calc is None): return
  recycle(f_calc, "f_calc", verbose=verbose)
  for column_root_label, column_types in [
        ("f_obs", None),
        ("Ework", "E")]:
    if (anomalous_flag and column_types == "E"): continue
    recycle(
      miller_array=abs(f_calc),
      column_root_label=column_root_label,
      column_types=column_types,
      verbose=verbose)
  if (not anomalous_flag):
    recycle(abs(f_calc), "f_obs", column_types="R", verbose=verbose)
  for column_root_label, column_types in [
        ("f_obs", None),
        ("Ework", "EQ")]:
    if (anomalous_flag and column_types == "EQ"): continue
    recycle(
      miller_array=miller.array(
        miller_set=f_calc,
        data=flex.abs(f_calc.data()),
        sigmas=flex.abs(f_calc.data())/10),
      column_root_label=column_root_label,
      column_types=column_types,
      verbose=verbose)
  recycle(f_calc.centric_flags(), "cent", verbose=verbose)
  recycle(generate_random_hl(miller_set=f_calc), "prob", verbose=verbose)

def run_call_back(flags, space_group_info):
  for anomalous_flag in [False, True]:
    exercise_recycle(
      space_group_info=space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def exercise_miller_array_data_types():
  miller_set = crystal.symmetry(
    unit_cell=(10,10,10,90,90,90),
    space_group_symbol="P1").miller_set(
      indices=flex.miller_index([(1,2,3),(4,5,6)]),
      anomalous_flag=False)
  for data in [
        flex.bool([False,True]),
        flex.int([0,1]),
        flex.size_t([0,1]),
        flex.double([0,1]),
        flex.complex_double([0,1])]:
    miller_array = miller_set.array(data=data)
    if (op.isfile("tmp.mtz")): os.remove("tmp.mtz")
    assert not op.isfile("tmp.mtz")
    miller_array.as_mtz_dataset(column_root_label="DATA").mtz_object().write(
      file_name="tmp.mtz")
    assert op.isfile("tmp.mtz")
    mtz_obj = mtz.object(file_name="tmp.mtz")
    miller_arrays_read_back = mtz_obj.as_miller_arrays()
    assert len(miller_arrays_read_back) == 1
    miller_array_read_back = miller_arrays_read_back[0]
    assert miller_array_read_back.indices().all_eq(miller_array.indices())
    if (miller_array.is_integer_array() or miller_array.is_bool_array()):
      assert miller_array_read_back.data().all_eq(flex.int([0, 1]))
    elif (miller_array.is_real_array()):
      assert miller_array_read_back.data().all_eq(flex.double([0, 1]))
    elif (miller_array.is_complex_array()):
      assert miller_array_read_back.data().all_eq(flex.complex_double([0, 1]))
    else:
      raise RuntimeError("Programming error.")

def exercise_extract_delta_anomalous():
  miller_array_start = miller.set(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10,10,10,90,90,90),
      space_group_symbol="P1"),
    indices=flex.miller_index([(1,2,3),(-1,-2,-3)]),
    anomalous_flag=True).array(
      data=flex.double([3,5]),
      sigmas=flex.double([0.3,0.5]))
  mtz_dataset = miller_array_start.as_mtz_dataset(column_root_label="F")
  mtz_object = mtz_dataset.mtz_object()
  fp = mtz_object.get_column(label="F(+)").extract_values()[0]
  fm = mtz_object.get_column(label="F(-)").extract_values()[0]
  sp = mtz_object.get_column(label="SIGF(+)").extract_values()[0]
  sm = mtz_object.get_column(label="SIGF(-)").extract_values()[0]
  # http://www.ccp4.ac.uk/dist/html/mtzMADmod.html
  # F = 0.5*( F(+) + F(-) )
  # D = F(+) - F(-)
  # SIGD = sqrt( SIGF(+)**2 + SIGF(-)**2 )
  # SIGF = 0.5*SIGD
  f = 0.5 * (fp + fm)
  d = fp - fm
  sigd = math.sqrt(sp**2 + sm**2)
  sigf = 0.5 * sigd
  mtz_dataset.add_column(label="F", type="F").set_reals(
    mtz_reflection_indices=flex.int([0]), data=flex.double([f]))
  mtz_dataset.add_column(label="SIGF", type="Q").set_reals(
    mtz_reflection_indices=flex.int([0]), data=flex.double([sigf]))
  mtz_dataset.add_column(label="D", type="D").set_reals(
    mtz_reflection_indices=flex.int([0]), data=flex.double([d]))
  mtz_dataset.add_column(label="SIGD", type="Q").set_reals(
    mtz_reflection_indices=flex.int([0]), data=flex.double([sigd]))
  miller_arrays = mtz_object.as_miller_arrays()
  assert len(miller_arrays) == 2
  miller_array = miller_arrays[0]
  assert list(miller_array.indices()) == [(1,2,3),(-1,-2,-3)]
  assert approx_equal(miller_array.data(), [3,5])
  assert approx_equal(miller_array.sigmas(), [0.3,0.5])
  miller_array = miller_arrays[1]
  assert list(miller_array.indices()) == [(1,2,3),(-1,-2,-3)]
  # F(+) = F + 0.5*D
  # F(-) = F - 0.5*D
  # SIGF(+) = sqrt( SIGF**2 + 0.25*SIGD**2 )
  # SIGF(-) = SIGF(+)
  fp = f + 0.5 * d
  fm = f - 0.5 * d
  assert approx_equal([fp,fm], miller_array.data())
  sp = math.sqrt(sigf**2 + 0.25 * sigd**2)
  sm = sp
  assert approx_equal([sp,sm], miller_array.sigmas())

def exercise_repair_ccp4i_import_merged_data():
  miller_array_start = miller.set(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10,10,10,90,90,90),
      space_group_symbol="P1"),
    indices=flex.miller_index([(1,2,3),(-1,-2,-3),(1,2,4),(-1,-2,-4)]),
    anomalous_flag=True).array(
      data=flex.double([3,5,6,7]),
      sigmas=flex.double([0.3,0.5,0.6,0.7]))
  mtz_dataset = miller_array_start.as_mtz_dataset(column_root_label="F")
  mtz_object = mtz_dataset.mtz_object()
  selection_valid = flex.bool([False,True])
  for sign in ["+", "-"]:
    for label in ["F(%s)"%sign, "SIGF(%s)"%sign]:
      column = mtz_object.get_column(label=label)
      values = column.extract_values()
      column.set_values(values=values, selection_valid=selection_valid)
    selection_valid = ~selection_valid
  miller_arrays = mtz_object.as_miller_arrays()
  assert len(miller_arrays) == 1
  assert not miller_arrays[0].anomalous_flag()
  assert list(miller_arrays[0].indices()) == [(-1,-2,-3),(1,2,4)]
  assert approx_equal(miller_arrays[0].data(), [5,6])
  assert approx_equal(miller_arrays[0].sigmas(), [0.5,0.6])

def exercise_hl_ab_only(anomalous_flag):
  cs = crystal.symmetry(
    unit_cell=(3.95738, 5.1446, 6.72755, 83, 109, 129),
    space_group_symbol="P1")
  if (not anomalous_flag):
    i = [
      (-1, 0, 1), (-1, 1, 1), (0, -1, 1), (0, 0, 1),
      (0, 0, 2), (0, 1, 0), (0, 1, 1), (1, -1, 0)]
  else:
    i = [
      (-1, 0, 1), (1, 0, -1), (-1, 1, 1), (1, -1, -1),
      (0, -1, 1), (0, 1, -1), (0, 0, 1), (0, 0, -1),
      (0, 0, 2), (0, 0, -2), (0, 1, 0), (0, -1, 0),
      (0, 1, 1), (0, -1, -1), (1, -1, 0), (-1, 1, 0)]
  ms = miller.set(
    crystal_symmetry=cs,
    indices=flex.miller_index(i),
    anomalous_flag=anomalous_flag)
  ma = ms.array(data=flex.size_t_range(ms.indices().size()).as_double()+1)
  mtz_dataset = ma.as_mtz_dataset(column_root_label="HA")
  columns = mtz_dataset.columns()
  if (not anomalous_flag):
    assert columns.size() == 4
    c = columns[3]
    assert c.label() == "HA"
    c.set_type("A")
    values = c.extract_values()
    selection_valid = c.selection_valid()
    c = mtz_dataset.add_column(label="HB", type="A")
    c.set_values(values=-values, selection_valid=selection_valid)
  else:
    assert columns.size() == 5
    for i,l in [(3,"HA(+)"), (4,"HA(-)")]:
      c = columns[i]
      assert c.label() == l
      if (i == 4):
        c.set_label("HB(+)")
      c.set_type("A")
    for i,l in [(3,"HA(-)"), (4,"HB(-)")]:
      c = columns[i]
      values = c.extract_values()
      selection_valid = c.selection_valid()
      c = mtz_dataset.add_column(label=l, type="A")
      c.set_values(values=-values, selection_valid=selection_valid)
  mtz_obj = mtz_dataset.mtz_object()
  mas = mtz_obj.as_miller_arrays()
  assert len(mas) == 1
  assert approx_equal(mas[0].indices(), ma.indices())
  if (not anomalous_flag):
    assert approx_equal(mas[0].data(), [
      (1, -1, 0, 0), (2, -2, 0, 0), (3, -3, 0, 0), (4, -4, 0, 0),
      (5, -5, 0, 0), (6, -6, 0, 0), (7, -7, 0, 0), (8, -8, 0, 0)])
  else:
    assert approx_equal(mas[0].data(), [
      (1, 2, 0, 0), (-1, -2, 0, 0), (3, 4, 0, 0), (-3, -4, 0, 0),
      (5, 6, 0, 0), (-5, -6, 0, 0), (7, 8, 0, 0), (-7, -8, 0, 0),
      (9, 10, 0, 0), (-9, -10, 0, 0), (11, 12, 0, 0), (-11, -12, 0, 0),
      (13, 14, 0, 0), (-13, -14, 0, 0), (15, 16, 0, 0), (-15, -16, 0, 0)])
  #
  columns = mtz_dataset.columns()
  columns[-1].set_type("F")
  try:
    mtz_obj.as_miller_arrays()
  except RuntimeError, e:
    if (not anomalous_flag):
      assert str(e) == 'Invalid MTZ column combination' \
        ' (incomplete Hendrickson-Lattman array),' \
        ' column labels: "HA", "HB" column types: "A", "F"'
    else:
      assert str(e) == 'Invalid MTZ column combination' \
        ' (incomplete Hendrickson-Lattman array),' \
        ' column labels: "HA(-)", "HB(-)" column types: "A", "F"'
  else: raise Exception_expected

def exercise():
  if (mtz is None):
    print "Skipping iotbx/mtz/tst.py: ccp4io not available"
    return
  exercise_extract_delta_anomalous()
  exercise_repair_ccp4i_import_merged_data()
  exercise_miller_array_data_types()
  for anomalous_flag in [False, True]:
    exercise_hl_ab_only(anomalous_flag=anomalous_flag)
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
