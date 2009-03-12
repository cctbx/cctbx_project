from cctbx import statistics
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import adptbx
from cctbx.development import debug_utils
from cctbx.development import random_structure
import sys

def exercise_sys_absent_intensity_distribution():
  xs = crystal.symmetry((3,4,5), "F222")
  mi = flex.miller_index(((1,2,3), (1,1,1), (1,2,2), (0,0,4)))
  ms = miller.set(xs, mi)
  data = flex.double((-1,-2,3,4))
  sigmas = flex.double((2,3,4,5))
  f_obs = miller.array(ms, data=data, sigmas=sigmas).set_observation_type_xray_intensity()
  dist = statistics.sys_absent_intensity_distribution(f_obs)
  assert approx_equal(dist.x,(-0.5,0.75))
  assert approx_equal(dist.y,(-1,3))
  assert approx_equal(dist.indices,((1,2,3), (1,2,2)))

def exercise_cumulative_intensity_distribution(space_group_info, anomalous_flag,
             d_min=1.0, reflections_per_bin=200, n_bins=10, verbose=0):
  elements = ("N", "C", "C", "O") * 5
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=True,
    use_u_aniso=False,
    u_iso=adptbx.b_as_u(10)
    ).structure_factors(
        anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct")
  if (0 or verbose):
    structure_factors.xray_structure().show_summary()
  for n_bins in [10,20]:
    f_calc = abs(structure_factors.f_calc())
    f_calc.setup_binner(
      auto_binning=True,
      reflections_per_bin=reflections_per_bin,
      n_bins=n_bins)
    f_obs = miller.array(
      miller_set=f_calc,
      data=f_calc.data()).set_observation_type_xray_amplitude()
    f_obs.use_binner_of(f_calc)
    dist = statistics.cumulative_intensity_distribution(f_obs)
    n_bins_used = f_obs.binner().n_bins_used()
    x_index = {}
    x_index[0.1] = int(n_bins_used/10)
    x_index[0.5] = int(n_bins_used/2)
    x_index[0.9] = int(9*n_bins_used/10)
    if space_group_info.group().is_centric():
      assert 0.23 < dist.y[x_index[0.1]] < 0.32
      assert 0.49 < dist.y[x_index[0.5]] < 0.57
      assert 0.65 < dist.y[x_index[0.9]] < 0.70
    else:
      assert 0.08 < dist.y[x_index[0.1]] < 0.21
      assert 0.37 < dist.y[x_index[0.5]] < 0.49
      assert 0.57 < dist.y[x_index[0.9]] < 0.65

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise_cumulative_intensity_distribution(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  exercise_sys_absent_intensity_distribution()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if __name__ == '__main__':
  run()
