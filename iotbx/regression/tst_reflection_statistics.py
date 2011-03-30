import libtbx.load_env
if (libtbx.env.has_module("ccp4io")):
  from iotbx.command_line import reflection_statistics
  from iotbx import mtz
else:
  mtz = None
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.development import debug_utils
from cStringIO import StringIO
import sys

def exercise_compare_cb_op_as_hkl():
  l = ["k,h,l", "h,k,l"]
  l.sort(reflection_statistics.compare_cb_op_as_hkl)
  assert l == ["h,k,l", "k,h,l"]
  l.sort(reflection_statistics.compare_cb_op_as_hkl)
  assert l == ["h,k,l", "k,h,l"]

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
  miller_arrays = []
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
    miller_arrays.append(
      subgroup_miller_array.f_sq_as_f().expand_to_p1().map_to_asu())
    file_names.append(file_name)
  return miller_arrays, file_names

def exercise_reflection_statistics(
      anomalous_flag,
      miller_arrays,
      file_names,
      verbose):
  out = StringIO()
  try:
    sys.stdout = out
    reflection_statistics.run(
      args=["--lattice_symmetry_max_delta=0.01"]+file_names)
  finally:
    sys.stdout = sys.__stdout__
    if (0 or verbose):
      sys.stdout.write(out.getvalue())
  done = set()
  for line in out.getvalue().splitlines():
    if (line.startswith("CC ")):
      type,i,j,cc,cb_op = line.split()[1:6]
      key = " ".join([type,i,j])
      cc = float(cc)
      if (key not in done):
        if (cc != 1.000):
          raise AssertionError(line.strip())
        done.add(key)
      # check reindexing matrices
      assert type in ["Obs", "Ano"]
      if (type == "Obs"):
        i, j = int(i)-1, int(j)-1
        ma_j = miller_arrays[j].change_basis(cb_op).expand_to_p1().map_to_asu()
        ma_i, ma_j = miller_arrays[i].common_sets(other=ma_j)
        if (float("%6.3f" % ma_i.correlation(ma_j).coefficient()) != cc):
          raise AssertionError(line.strip())
  n = len(file_names)
  expected_number_of_cc_lines = n*(n-1)//2
  if (anomalous_flag):
    expected_number_of_cc_lines *= 2
  assert len(done) == expected_number_of_cc_lines

def run_call_back(flags, space_group_info):
  for anomalous_flag in [False, True]:
    if (flags.Verbose):
      print space_group_info, "anomalous_flag:", anomalous_flag
    if (anomalous_flag and space_group_info.group().is_centric()): continue
    miller_arrays, file_names = generate_mtz_files(
      space_group_info=space_group_info,
      anomalous_flag=anomalous_flag)
    exercise_reflection_statistics(
      anomalous_flag=anomalous_flag,
      miller_arrays=miller_arrays,
      file_names=file_names,
      verbose=flags.Verbose)

def exercise():
  if (mtz is None):
    print \
      "Skipping iotbx/regression/tst_reflection_statistics.py:" \
      " ccp4io not available"
    return
  exercise_compare_cb_op_as_hkl()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
