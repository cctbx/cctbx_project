from iotbx.reflection_file_utils import reflection_file_server
from iotbx import reflection_file_reader
from iotbx import mtz
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import UserError
from cStringIO import StringIO
import os

def exercise_get_xtal_data():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  input_arrays = [miller_set.array(
    data=flex.random_double(size=miller_set.indices().size()),
    sigmas=flex.random_double(size=miller_set.indices().size())/10)
      .set_observation_type_xray_intensity()
        for i in [0,1]]
  mtz_dataset = input_arrays[0].as_mtz_dataset(column_root_label="F0")
  mtz_dataset.mtz_object().write("tmp1.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp1.mtz")]
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files)
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=None,
    ignore_all_zeros=False,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp1.mtz:F0,SIGF0"
  mtz_dataset.add_miller_array(
    miller_array=input_arrays[1], column_root_label="F1")
  mtz_dataset.mtz_object().write("tmp2.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp2.mtz")]
  err = StringIO()
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=err)
  try:
    f_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=None,
      ignore_all_zeros=True,
      parameter_scope="xray_data")
  except UserError:
    assert err.getvalue() == """\

Multiple equally suitable arrays of observed xray data found.

Possible choices:
  tmp2.mtz:F0,SIGF0
  tmp2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=["F1", "SIGF1"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp2.mtz:F1,SIGF1"
  try:
    f_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=["F1", "SIGF0"],
      ignore_all_zeros=True,
      parameter_scope="xray_data")
  except UserError:
    assert err.getvalue() == """\

No matching array: xray_data.labels=F1 SIGF0

Possible choices:
  tmp2.mtz:F0,SIGF0
  tmp2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=["1"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp2.mtz:F0,SIGF0"
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=["2"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp2.mtz:F1,SIGF1"
  f_obs = reflection_file_srv.get_xray_data(
    file_name="tmp2.mtz",
    labels=["2"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp2.mtz:F1,SIGF1"
  assert len(reflection_file_srv.file_name_miller_arrays) == 1
  f_obs = reflection_file_srv.get_xray_data(
    file_name="tmp1.mtz",
    labels=None,
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(f_obs.info()) == "tmp1.mtz:F0,SIGF0"
  f_obs = reflection_file_srv.get_xray_data(
    file_name=os.path.abspath("tmp1.mtz"),
    labels=["sigf0"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(f_obs.info()) == "tmp1.mtz:F0,SIGF0"
  try:
    f_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=None,
      ignore_all_zeros=True,
      parameter_scope="xray_data")
  except UserError:
    assert err.getvalue() == """\

Multiple equally suitable arrays of observed xray data found.

Possible choices:
  tmp2.mtz:F0,SIGF0
  tmp2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()

def exercise_get_r_free_flags():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(30,31,32,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  n = miller_set.indices().size()
  exercise_flag_arrays = []
  exercise_flag_arrays.append(
    flex.int(list(flex.random_permutation(size=n)%10)))
  exercise_flag_arrays.append(flex.int(xrange(n)))
  exercise_flag_arrays.append(flex.int(n, 0))
  for style in ["ccp4", "cns", "shelx", "bool"]:
    for i_exercise,exercise_flag_array in enumerate(exercise_flag_arrays):
      for reversed in [False, True]:
        if (style == "ccp4"):
          if (reversed): break
          data = exercise_flag_array
          test_flag_value = 3
        else:
          if (not reversed):
            data = (exercise_flag_array == 0)
            test_flag_value = True
          else:
            data = (exercise_flag_array != 0)
            test_flag_value = False
          if (style == "cns"):
            data = data.as_int()
            test_flag_value = int(test_flag_value)
          elif (style == "shelx"):
            data = -data.as_int()
            data.set_selected((data == 0), 1)
            if (not reversed): test_flag_value = -1
            else:              test_flag_value = 1
        input_array = miller_set.array(data=data)
        mtz_dataset = input_array.as_mtz_dataset(
          column_root_label="FreeRflags")
        mtz_dataset.mtz_object().write("tmp.mtz")
        reflection_files = [reflection_file_reader.any_reflection_file(
          file_name="tmp.mtz")]
        err = StringIO()
        reflection_file_srv = reflection_file_server(
          crystal_symmetry=crystal_symmetry,
          force_symmetry=True,
          reflection_files=reflection_files,
          err=err)
        for trial_test_flag_value in [None, test_flag_value]:
          for trial_label in [None, "free", "foo"]:
            try:
              r_free_flags, actual_test_flag_value = \
                reflection_file_srv.get_r_free_flags(
                  file_name=None,
                  label=trial_label,
                  test_flag_value=trial_test_flag_value,
                  disable_suitability_test=False,
                  parameter_scope="r_free_flags")
            except UserError, e:
              if (trial_label != "foo"):
                assert i_exercise > 0
                if (trial_label is None):
                  assert str(e) == """\
No array of R-free flags found.

For manual selection define:
  r_free_flags.label
  r_free_flags.test_flag_value
  r_free_flags.disable_suitability_test=True"""
                else:
                  assert str(e) == \
                      "Not a suitable array of R-free flags:" \
                    + " r_free_flags.label=free\n" \
                    + "To override the suitability test define:" \
                    + " r_free_flags.disable_suitability_test=True"
              else:
                assert str(e) == "No matching array: r_free_flags.label=foo"
                if (i_exercise == 0):
                  assert err.getvalue() == """\

No matching array: r_free_flags.label=foo

Possible choices:
  tmp.mtz:FreeRflags

Please use r_free_flags.label
to specify an unambiguous substring of the target label.

"""
                else:
                  assert err.getvalue() == """\

No matching array: r_free_flags.label=foo

"""
              err = reflection_file_srv.err = StringIO()
            else:
              assert i_exercise == 0
  for second_label in ["test", "foo"]:
    input_array = miller_set.array(data=exercise_flag_arrays[0])
    mtz_dataset = input_array.as_mtz_dataset(
      column_root_label="FreeRflags")
    mtz_dataset.add_miller_array(
      miller_array=input_array,
      column_root_label=second_label)
    mtz_dataset.mtz_object().write("tmp.mtz")
    reflection_files = [reflection_file_reader.any_reflection_file(
      file_name="tmp.mtz")]
    err = StringIO()
    reflection_file_srv = reflection_file_server(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=True,
      reflection_files=reflection_files,
      err=err)
    try:
      reflection_file_srv.get_r_free_flags(
        file_name=None,
        label=None,
        test_flag_value=None,
        disable_suitability_test=False,
        parameter_scope="r_free_flags")
    except UserError, e:
      assert str(e)=="Multiple equally suitable arrays of R-free flags found."
      assert err.getvalue() == """\

Multiple equally suitable arrays of R-free flags found.

Possible choices:
  tmp.mtz:FreeRflags
  tmp.mtz:test

Please use r_free_flags.label
to specify an unambiguous substring of the target label.

"""
      err = reflection_file_srv.err = StringIO()
    else:
      assert str(r_free_flags.info()) == "tmp.mtz:FreeRflags"
  r_free_flags, actual_test_flag_value = \
    reflection_file_srv.get_r_free_flags(
      file_name=None,
      label="FreeRflags",
      test_flag_value=3,
      disable_suitability_test=True,
      parameter_scope="r_free_flags")
  assert r_free_flags.info().label_string() == "FreeRflags"
  assert actual_test_flag_value == 3
  for label,test_flag_value in [(None,3), ("FreeRflags",None)]:
    try:
      reflection_file_srv.get_r_free_flags(
        file_name=None,
        label=label,
        test_flag_value=test_flag_value,
        disable_suitability_test=True,
        parameter_scope="r_free_flags")
    except UserError, e:
      assert str(e) == "r_free_flags.disable_suitability_test=True:" \
        " Suitability test for R-free flags can only be disabled if both" \
        " r_free_flags.label and r_free_flags.test_flag_value are defined."
    else: raise RuntimeError("Exception expected.")

def exercise_get_experimental_phases():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(30,31,32,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  input_array = miller_set.array(
    data=flex.hendrickson_lattman(miller_set.indices().size(), (0,0,0,0)))
  mtz_dataset = input_array.as_mtz_dataset(column_root_label="P")
  mtz_dataset.mtz_object().write("tmp.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp.mtz")]
  err = StringIO()
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=err)
  experimental_phases = reflection_file_srv.get_experimental_phases(
    file_name=None,
    labels=None,
    ignore_all_zeros=False,
    parameter_scope="experimental_phases")
  assert str(experimental_phases.info()) == "tmp.mtz:PA,PB,PC,PD"
  try:
    reflection_file_srv.get_experimental_phases(
      file_name=None,
      labels=None,
      ignore_all_zeros=True,
      parameter_scope="experimental_phases")
  except UserError, e:
    assert str(e) == "No array of experimental phases found."
    assert err.getvalue() == """\

No array of experimental phases found.

"""
  else: raise RuntimeError("Exception expected.")

def exercise():
  exercise_get_xtal_data()
  exercise_get_r_free_flags()
  exercise_get_experimental_phases()
  print "OK"

if (__name__ == "__main__"):
  exercise()
