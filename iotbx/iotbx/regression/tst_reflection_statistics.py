from cctbx.array_family import flex
from iotbx.command_line import reflection_statistics
from iotbx import mtz
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.development import debug_utils
from cStringIO import StringIO
import sys

def generate_mtz_files(space_group_info, anomalous_flag):
  crystal_symmetry = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    space_group_info=space_group_info)
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=anomalous_flag,
    d_min=1)
  miller_array = miller.array(
    miller_set=miller_set,
    data=flex.random_double(size=miller_set.indices().size()))
  miller_array_p1 = miller_array.expand_to_p1()
  file_names = []
  subgrs = subgroups.subgroups(space_group_info).groups_parent_setting()
  for i_subgroup, subgroup in enumerate(subgrs):
    subgroup_miller_array = miller_array_p1.customized_copy(
      space_group_info=sgtbx.space_group_info(group=subgroup)) \
        .merge_equivalents() \
        .array() \
        .as_reference_setting() \
        .set_observation_type_xray_intensity()
    file_name = "tmp%d.mtz" % i_subgroup
    mtz_object = subgroup_miller_array.as_mtz_dataset(
      column_root_label="FOBS").mtz_object().write(file_name=file_name)
    file_names.append(file_name)
  return file_names

def exercise_reflection_statistics(anomalous_flag, file_names, verbose):
  out = StringIO()
  try:
    sys.stdout = out
    reflection_statistics.run(args=file_names)
  finally:
    sys.stdout = sys.__stdout__
    if (0 or verbose):
      sys.stdout.write(out.getvalue())
  done = {}
  for line in out.getvalue().splitlines():
    if (line.startswith("CC ")):
      type,i,j,cc = line.split()[1:5]
      key = " ".join([type,i,j])
      cc = float(cc)
      if (key not in done):
        if (cc != 1.000):
          raise AssertionError(line.strip())
        done[key] = None
  n = len(file_names)
  expected_number_of_cc_lines = n*(n-1)/2
  if (anomalous_flag):
    expected_number_of_cc_lines *= 2
  assert len(done) == expected_number_of_cc_lines

def run_call_back(flags, space_group_info):
  for anomalous_flag in [False, True]:
    if (flags.Verbose):
      print space_group_info, "anomalous_flag:", anomalous_flag
    if (anomalous_flag and space_group_info.group().is_centric()): continue
    file_names = generate_mtz_files(
      space_group_info=space_group_info,
      anomalous_flag=anomalous_flag)
    exercise_reflection_statistics(
      anomalous_flag=anomalous_flag,
      file_names=file_names,
      verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
