from __future__ import absolute_import, division, print_function
import libtbx.load_env
from six.moves import range
if (libtbx.env.has_module("ccp4io")):
  from iotbx import reflection_file_reader
  from iotbx.reflection_file_utils import reflection_file_server, \
    guess_r_free_flag_value
  from iotbx import mtz
else:
  mtz = None
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, show_diff
from libtbx.utils import Sorry, null_out
from six.moves import cStringIO as StringIO
import os

def exercise_get_amplitudes_and_get_phases_deg():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  input_arrays = [miller_set.array(
    data=flex.random_double(size=miller_set.indices().size()))
      .set_observation_type_xray_amplitude()
        for i in [0,1]]
  mtz_dataset = input_arrays[0].as_mtz_dataset(column_root_label="F0")
  mtz_dataset.mtz_object().write("tmp_rfu1.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu1.mtz")]
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files)
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=None,
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu1.mtz:F0"
  ampl = reflection_file_srv.get_miller_array(labels="F0")
  assert str(ampl.info()) == "tmp_rfu1.mtz:F0"
  mtz_dataset.add_miller_array(
    miller_array=input_arrays[1], column_root_label="F1")
  mtz_dataset.mtz_object().write("tmp_rfu2.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu2.mtz")]
  err = StringIO()
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=err)
  try:
    reflection_file_srv.get_amplitudes(
      file_name=None,
      labels=None,
      convert_to_amplitudes_if_necessary=True,
      parameter_scope="amplitudes",
      parameter_name="labels")
  except Sorry:
    assert not show_diff(err.getvalue(), """\

Multiple equally suitable arrays of amplitudes found.

Possible choices:
  tmp_rfu2.mtz:F0
  tmp_rfu2.mtz:F1

Please use amplitudes.labels
to specify an unambiguous substring of the target label.

""")
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=["F1"],
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu2.mtz:F1"
  try:
    reflection_file_srv.get_amplitudes(
      file_name=None,
      labels=["F2"],
      convert_to_amplitudes_if_necessary=True,
      parameter_name="labels",
      parameter_scope=None)
  except Sorry:
    assert not show_diff(err.getvalue(), """\

No matching array: labels=F2

Possible choices:
  tmp_rfu2.mtz:F0
  tmp_rfu2.mtz:F1

Please use labels
to specify an unambiguous substring of the target label.

""")
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  assert len(reflection_file_srv.file_name_miller_arrays) == 1
  ampl = reflection_file_srv.get_amplitudes(
    file_name="tmp_rfu1.mtz",
    labels=None,
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(ampl.info()) == "tmp_rfu1.mtz:F0"
  ampl = reflection_file_srv.get_amplitudes(
    file_name=os.path.abspath("tmp_rfu1.mtz"),
    labels=["f0"],
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(ampl.info()) == "tmp_rfu1.mtz:F0"
  try:
    reflection_file_srv.get_amplitudes(
      file_name=None,
      labels=None,
      convert_to_amplitudes_if_necessary=True,
      parameter_scope=None,
      parameter_name=None)
  except Sorry:
    assert not show_diff(err.getvalue(), """\

Multiple equally suitable arrays of amplitudes found.

Possible choices:
  tmp_rfu2.mtz:F0
  tmp_rfu2.mtz:F1

Please specify an unambiguous substring of the target label.

""")
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  #
  mtz_dataset.add_miller_array(
    miller_array=miller_set.array(
      data=flex.polar(
        flex.random_double(size=miller_set.indices().size()),
        flex.random_double(size=miller_set.indices().size()))),
    column_root_label="F2")
  mtz_dataset.add_miller_array(
    miller_array=miller_set.array(
      data=flex.random_double(size=miller_set.indices().size()),
      sigmas=flex.random_double(size=miller_set.indices().size())/10)
        .set_observation_type_xray_intensity(),
    column_root_label="F3")
  mtz_dataset.add_miller_array(
    miller_array=miller_set.array(
      data=flex.hendrickson_lattman(miller_set.indices().size(), (0,0,0,0))),
    column_root_label="P")
  mtz_dataset.mtz_object().write("tmp_rfu3.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu3.mtz")]
  err = StringIO()
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=err)
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=["f2"],
    convert_to_amplitudes_if_necessary=False,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu3.mtz:F2,PHIF2"
  assert ampl.is_complex_array()
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=["f2"],
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu3.mtz:F2"
  assert ampl.is_real_array()
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=["f3"],
    convert_to_amplitudes_if_necessary=False,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu3.mtz:F3,SIGF3"
  assert ampl.is_xray_intensity_array()
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=["f3"],
    convert_to_amplitudes_if_necessary=True,
    parameter_scope="amplitudes",
    parameter_name="labels")
  assert str(ampl.info()) == "tmp_rfu3.mtz:F3,as_amplitude_array"
  assert ampl.is_real_array()
  ampl = reflection_file_srv.get_amplitudes(
    file_name=None,
    labels=None,
    convert_to_amplitudes_if_necessary=False,
    parameter_scope="amplitudes",
    parameter_name="labels",
    return_all_valid_arrays=True,
    strict=True)
  assert (len(ampl) == 2)
  for f in ampl :
    assert (not f.is_xray_intensity_array()) and (not f.is_complex_array())
  #
  phases = reflection_file_srv.get_phases_deg(
    file_name=None,
    labels=["f2"],
    convert_to_phases_if_necessary=False,
    original_phase_units=None,
    parameter_scope="phases",
    parameter_name="labels")
  assert str(phases.info()) == "tmp_rfu3.mtz:F2,PHIF2"
  assert phases.is_complex_array()
  phases = reflection_file_srv.get_phases_deg(
    file_name=None,
    labels=["f2"],
    convert_to_phases_if_necessary=True,
    original_phase_units=None,
    parameter_scope=None,
    parameter_name="labels")
  assert str(phases.info()) == "tmp_rfu3.mtz:PHIF2"
  assert phases.is_real_array()
  assert flex.mean(phases.data()) > 5
  phases = reflection_file_srv.get_phases_deg(
    file_name=None,
    labels=["PA"],
    convert_to_phases_if_necessary=False,
    original_phase_units=None,
    parameter_scope="phases",
    parameter_name="labels")
  assert str(phases.info()) == "tmp_rfu3.mtz:PA,PB,PC,PD"
  phases = reflection_file_srv.get_phases_deg(
    file_name=None,
    labels=["PA"],
    convert_to_phases_if_necessary=True,
    original_phase_units=None,
    parameter_scope="phases",
    parameter_name="labels")
  assert str(phases.info()) \
      == "tmp_rfu3.mtz:PA,PB,PC,PD,converted_to_centroid_phases"
  assert phases.is_real_array()
  for original_phase_units in [None, "deg", "rad"]:
    phases = reflection_file_srv.get_phases_deg(
      file_name=None,
      labels=["F0"],
      convert_to_phases_if_necessary=False,
      original_phase_units=original_phase_units,
      parameter_scope=None,
      parameter_name="labels")
    if (original_phase_units != "rad"):
      assert str(phases.info()) == "tmp_rfu3.mtz:F0"
    else:
      assert str(phases.info()) == "tmp_rfu3.mtz:F0,converted_to_deg"

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
  mtz_dataset.mtz_object().write("tmp_rfu1.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu1.mtz")]
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files)
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=None,
    ignore_all_zeros=False,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp_rfu1.mtz:F0,SIGF0"
  mtz_dataset.add_miller_array(
    miller_array=input_arrays[1], column_root_label="F1")
  mtz_dataset.mtz_object().write("tmp_rfu2.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu2.mtz")]
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
  except Sorry:
    assert err.getvalue() == """\

Multiple equally suitable arrays of observed xray data found.

Possible choices:
  tmp_rfu2.mtz:F0,SIGF0
  tmp_rfu2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  f_obs_list = reflection_file_srv.get_xray_data(
    file_name = None,
    labels = None,
    ignore_all_zeros=True,
    parameter_scope="xray_data",
    return_all_valid_arrays=True,
    minimum_score=1)
  assert len(f_obs_list) == 2
  f_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=["F1", "SIGF1"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert str(f_obs.info()) == "tmp_rfu2.mtz:F1,SIGF1"
  try:
    f_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=["F1", "SIGF0"],
      ignore_all_zeros=True,
      parameter_scope="xray_data")
  except Sorry:
    assert err.getvalue() == """\

No matching array: xray_data.labels=F1 SIGF0

Possible choices:
  tmp_rfu2.mtz:F0,SIGF0
  tmp_rfu2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  assert len(reflection_file_srv.file_name_miller_arrays) == 1
  f_obs = reflection_file_srv.get_xray_data(
    file_name="tmp_rfu1.mtz",
    labels=None,
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(f_obs.info()) == "tmp_rfu1.mtz:F0,SIGF0"
  f_obs = reflection_file_srv.get_xray_data(
    file_name=os.path.abspath("tmp_rfu1.mtz"),
    labels=["sigf0"],
    ignore_all_zeros=True,
    parameter_scope="xray_data")
  assert len(reflection_file_srv.file_name_miller_arrays) == 2
  assert str(f_obs.info()) == "tmp_rfu1.mtz:F0,SIGF0"
  try:
    f_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=None,
      ignore_all_zeros=True,
      parameter_scope="xray_data")
  except Sorry:
    assert err.getvalue() == """\

Multiple equally suitable arrays of observed xray data found.

Possible choices:
  tmp_rfu2.mtz:F0,SIGF0
  tmp_rfu2.mtz:F1,SIGF1

Please use xray_data.labels
to specify an unambiguous substring of the target label.

"""
    err = reflection_file_srv.err = StringIO()
  else:
    raise Exception_expected
  # test preference for anomalous (or merged) data
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=True,
    d_min=3)
  i_obs = miller_set.array(
    data=flex.random_double(size=miller_set.indices().size()),
    sigmas=flex.random_double(size=miller_set.indices().size())/10
      ).set_observation_type_xray_intensity()
  i_mean = i_obs.average_bijvoet_mates()
  mtz_data = i_obs.as_mtz_dataset(column_root_label="I")
  mtz_data.add_miller_array(i_mean, column_root_label="I")
  mtz_data.mtz_object().write("tmp_rfu3.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp_rfu3.mtz")]
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files)
  err = reflection_file_srv.err = StringIO()
  try :
    i_obs = reflection_file_srv.get_xray_data(
      file_name=None,
      labels=None,
      ignore_all_zeros=False,
      parameter_scope="xray_data")
  except Sorry :
    pass
  i_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=None,
    ignore_all_zeros=False,
    parameter_scope="xray_data",
    prefer_anomalous=True)
  assert (i_obs.info().label_string() == "I(+),SIGI(+),I(-),SIGI(-)")
  i_obs = reflection_file_srv.get_xray_data(
    file_name=None,
    labels=None,
    ignore_all_zeros=False,
    parameter_scope="xray_data",
    prefer_anomalous=False)
  assert (i_obs.info().label_string() == "I,SIGI")

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
  exercise_flag_arrays.append(flex.int(range(n)))
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
            except Sorry as e:
              if (trial_label != "foo"):
                assert i_exercise > 0
                if (trial_label is None):
                  assert str(e) == """\
No array of R-free flags found.

For manual selection define:
  r_free_flags.test_flag_value
  r_free_flags.disable_suitability_test=True"""
                else:
                  assert str(e) == \
                      "Not a suitable array of R-free flags:" \
                    + " r_free_flags=free\n" \
                    + "To override the suitability test define:" \
                    + " r_free_flags.disable_suitability_test=True"
              else:
                assert str(e) == "No matching array: r_free_flags=foo"
                if (i_exercise == 0):
                  assert err.getvalue() == """\

No matching array: r_free_flags=foo

Possible choices:
  tmp.mtz:FreeRflags

Please use r_free_flags
to specify an unambiguous substring of the target label.

"""
                else:
                  assert err.getvalue() == """\

No matching array: r_free_flags=foo

"""
              err = reflection_file_srv.err = StringIO()
            else:
              assert i_exercise == 0
              actual_test_flag_value_2 = guess_r_free_flag_value(
                miller_array=r_free_flags,
                test_flag_value=trial_test_flag_value)
              assert (actual_test_flag_value_2 == actual_test_flag_value)
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
    except Sorry as e:
      assert str(e)=="Multiple equally suitable arrays of R-free flags found."
      assert err.getvalue() == """\

Multiple equally suitable arrays of R-free flags found.

Possible choices:
  tmp.mtz:FreeRflags
  tmp.mtz:test

Please use r_free_flags
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
    except Sorry as e:
      assert str(e) == "r_free_flags.disable_suitability_test=True:" \
        " Suitability test for R-free flags can only be disabled if both" \
        " r_free_flags.label and r_free_flags.test_flag_value are defined."
    else: raise Exception_expected
  # test corrupted R-free flags
  r_free_flags = miller_set.generate_r_free_flags()
  int_flags = r_free_flags.data().as_int()
  int_flags[100] = 10000000
  r_free_flags = r_free_flags.customized_copy(data=int_flags)
  mtz_dataset = r_free_flags.as_mtz_dataset(
    column_root_label="TEST")
  mtz_dataset.mtz_object().write("tmp.mtz")
  reflection_files = [reflection_file_reader.any_reflection_file(
    file_name="tmp.mtz")]
  err = StringIO()
  reflection_file_srv = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=err)
  flags, value = reflection_file_srv.get_r_free_flags(
    file_name=None,
    label=None,
    test_flag_value=None,
    disable_suitability_test=False,
    parameter_scope="r_free_flags")
  assert (value == 1)

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
  except Sorry as e:
    assert str(e) == "No array of experimental phases found."
    assert err.getvalue() == """\

No array of experimental phases found.

"""
  else: raise Exception_expected

def exercise_hklf_plus_ins_or_res():
  import textwrap
  from os import path
  import shutil
  from libtbx.test_utils import approx_equal

  try:
    hklf4 = ['   1   0   0    1.11    0.11',
             '   0   1  -1    2.22    0.22',
             '  -1   0   1    4.44    0.44']
    hklf4 = [li + '\n' for li in hklf4]
    ins = textwrap.dedent(
    """\
    TITL 03srv209 in Pbca
    CELL 0.71073 7.35 9.541 12.842 90 90 90
    ZERR 4 0.002 0.002 0.003 0 0 0
    LATT 1
    SYMM 0.5-X,-Y,0.5+Z
    SYMM -X,0.5+Y,0.5-Z
    SYMM 0.5+X,0.5-Y,-Z
    SFAC C H O N
    UNIT 32 40 16 8
    HKLF 4
    """)
    folder = 'temp_data_hklf_plus_ins_res'
    if not path.exists(folder):
      os.mkdir(folder)
    insfn = path.join(folder, '03srv209.ins')
    hklfn = path.join(folder, '03srv209.hkl')
    with open(insfn, 'w') as f:
      f.write(ins)
    with open(hklfn, 'w') as f:
      f.writelines(hklf4)
    rf = reflection_file_reader.any_reflection_file(hklfn+'=hklf+ins/res')
    ma = rf.as_miller_arrays()[0]
    cs = ma.crystal_symmetry()
    assert cs.is_identical_symmetry(
      crystal.symmetry(unit_cell='7.35 9.541 12.842 90 90 90',
                       space_group_symbol='Pbca'))
    assert all(ma.indices() == flex.miller_index([(1,0,0), (0,1,-1), (-1,0,1)]))
    assert approx_equal(ma.data(), [1.11, 2.22, 4.44])
    assert approx_equal(ma.sigmas(), [0.11, 0.22, 0.44])
  finally:
    shutil.rmtree(folder, ignore_errors=True)

def exercise_extract_miller_array_from_file():
  from iotbx import reflection_file_utils as rfu
  from libtbx.test_utils import approx_equal
  log = null_out()
  sorry_counts = 0
  crystal_symmetry = crystal.symmetry(
    unit_cell=(30,31,32,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  size = miller_set.indices().size()
  a1 = miller_set.array(
    data=flex.hendrickson_lattman(size, (1,1,1,1)))
  a2 = miller_set.array(data=flex.double(size, 2))
  a3 = miller_set.array(data=flex.double(size, 3))
  a4 = miller_set.array(data=flex.complex_double(size, 4+4j))
  a5 = miller_set.array(data=flex.complex_double(size, 5+5j))
  #
  mtz_dataset = a1.as_mtz_dataset(column_root_label="A1")
  mtz_dataset.mtz_object().write("tmp.mtz")
  ma = rfu.extract_miller_array_from_file(file_name="tmp.mtz", log=log)
  assert type(ma.data()) == flex.hendrickson_lattman
  #
  mtz_dataset = a5.as_mtz_dataset(column_root_label="A5")
  mtz_dataset.mtz_object().write("tmp.mtz")
  ma = rfu.extract_miller_array_from_file(file_name="tmp.mtz", log=log)
  assert type(ma.data()) == flex.complex_double
  #
  for tp in [None, "complex"]:
    mtz_dataset = a4.as_mtz_dataset(column_root_label="A4")
    mtz_dataset.add_miller_array(
      miller_array=a5, column_root_label="A5")
    mtz_dataset.mtz_object().write("tmp.mtz")
    try:
      rfu.extract_miller_array_from_file(file_name="tmp.mtz",type=tp, log=log)
    except Sorry as e:
      assert ("Multiple choices available." in str(e))
      sorry_counts += 1
  #
  for tp in [None, "real"]:
    mtz_dataset = a2.as_mtz_dataset(column_root_label="A2")
    mtz_dataset.add_miller_array(
      miller_array=a3, column_root_label="A3")
    mtz_dataset.mtz_object().write("tmp.mtz")
    try: rfu.extract_miller_array_from_file(file_name="tmp.mtz",type=tp,log=log)
    except Sorry as e:
      assert ("Multiple choices available." in str(e))
      sorry_counts += 1
  #
  mtz_dataset = a3.as_mtz_dataset(column_root_label="A3")
  mtz_dataset.add_miller_array(
    miller_array=a4, column_root_label="A4")
  mtz_dataset.mtz_object().write("tmp.mtz")
  try: rfu.extract_miller_array_from_file(file_name="tmp.mtz",log=log)
  except Sorry as e:
    assert ("Multiple choices available." in str(e))
    sorry_counts += 1
  #
  mtz_dataset = a4.as_mtz_dataset(column_root_label="A4")
  mtz_dataset.add_miller_array(
    miller_array=a5, column_root_label="A5")
  mtz_dataset.mtz_object().write("tmp.mtz")
  try:
    rfu.extract_miller_array_from_file(file_name="tmp.mtz",type="real",
      log=log)
  except Sorry as e:
    assert str(e)=="No suitable arrays."
    sorry_counts += 1
  #
  mtz_dataset = a2.as_mtz_dataset(column_root_label="A2")
  mtz_dataset.add_miller_array(
    miller_array=a3, column_root_label="A3")
  mtz_dataset.mtz_object().write("tmp.mtz")
  try:
    rfu.extract_miller_array_from_file(file_name="tmp.mtz",type="complex",
      log=log)
  except Sorry as e:
    assert str(e)=="No suitable arrays."
    sorry_counts += 1
  #
  mtz_dataset = a4.as_mtz_dataset(column_root_label="A4")
  mtz_dataset.add_miller_array(
    miller_array=a5, column_root_label="A5")
  mtz_dataset.mtz_object().write("tmp.mtz")
  ma = rfu.extract_miller_array_from_file(file_name="tmp.mtz",label="A5,PHIA5",
    log=log)
  assert approx_equal(ma.data()[0], 5+5j)
  #
  mtz_dataset = a4.as_mtz_dataset(column_root_label="A4")
  mtz_dataset.add_miller_array(
    miller_array=a5, column_root_label="A5")
  mtz_dataset.mtz_object().write("tmp.mtz")
  try:
    rfu.extract_miller_array_from_file(file_name="tmp.mtz",
      label="A5,PHIA5", type="real", log=log)
  except Sorry as e:
    assert str(e)=="No suitable arrays."
    sorry_counts += 1
  #
  assert sorry_counts == 8

def exercise_automation_wrappers():
  from iotbx.reflection_file_utils import process_raw_data, \
    change_space_group, load_f_obs_and_r_free
  from cctbx import sgtbx
  from libtbx.test_utils import approx_equal
  mtz_file = "tmp_iotbx_reflection_file_utils.mtz"
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,90,95,90),
    space_group_symbol="P 2")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=True,
    d_min=1.5)
  n_obs = miller_set.indices().size()
  i_obs = miller_set.array(
    data=flex.random_double(size=n_obs)).set_observation_type_xray_intensity()
  i_obs = i_obs.customized_copy(sigmas=flex.sqrt(i_obs.data()))
  r_free_flags = miller_set.generate_r_free_flags()
  r_free_flags_partial = r_free_flags.select(flex.random_bool(n_obs, 0.9))
  out = StringIO()
  processed = process_raw_data(
    obs=i_obs,
    r_free_flags=None,
    test_flag_value=None,
    log=out)
  assert ("""WARNING: R-free flags not supplied.""" in out.getvalue())
  assert (processed.data_labels() == "F(+),SIGF(+),F(-),SIGF(-)")
  assert (processed.phase_labels() is None)
  assert (processed.flags_are_new())
  out2 = StringIO()
  processed2 = process_raw_data(
    obs=i_obs,
    r_free_flags=r_free_flags_partial,
    test_flag_value=True,
    log=out2)
  assert ("""WARNING: R-free flags are incomplete""" in out2.getvalue())
  assert (not processed2.flags_are_new())
  assert (processed.n_obs() == processed2.n_obs())
  processed.write_mtz_file(mtz_file, title="tst_iotbx", wavelength=0.9792)
  f_obs, r_free = load_f_obs_and_r_free(mtz_file)
  change_space_group(mtz_file, sgtbx.space_group_info("P21"))
  f_obs_new, r_free_new = load_f_obs_and_r_free(mtz_file)
  assert (f_obs_new.size() == f_obs.size() - 4)
  f_obs_new, r_free_new = load_f_obs_and_r_free(mtz_file,
    anomalous_flag=True)
  assert (str(f_obs_new.space_group_info()) == "P 1 21 1")
  assert (approx_equal(f_obs_new.info().wavelength, 0.9792))

def exercise():
  if (mtz is None):
    print("Skipping iotbx/tst_reflection_file_utils.py: ccp4io not available")
    return
  exercise_hklf_plus_ins_or_res()
  exercise_get_amplitudes_and_get_phases_deg()
  exercise_get_xtal_data()
  exercise_get_r_free_flags()
  exercise_get_experimental_phases()
  exercise_extract_miller_array_from_file()
  exercise_automation_wrappers()

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
