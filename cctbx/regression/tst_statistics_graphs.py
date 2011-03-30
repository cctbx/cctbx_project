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
    f_obs_sq = f_obs.f_as_f_sq()
    f_obs_sq.use_binner_of(f_obs)
    for ma in (f_obs,f_obs_sq):
      if ma.is_xray_amplitude_array():
        dist = statistics.cumulative_intensity_distribution(f_obs=ma)
      else:
        dist = statistics.cumulative_intensity_distribution(f_obs_sq=ma)
      n_bins_used = f_obs.binner().n_bins_used()
      x_index = {}
      x_index[0.1] = n_bins_used//10
      x_index[0.5] = n_bins_used//2
      x_index[0.9] = 9*n_bins_used//10
      if space_group_info.group().is_centric():
        assert 0.23 < dist.y[x_index[0.1]] < 0.32
        assert 0.49 < dist.y[x_index[0.5]] < 0.57
        assert 0.65 < dist.y[x_index[0.9]] < 0.70
      else:
        assert 0.08 < dist.y[x_index[0.1]] < 0.21
        assert 0.37 < dist.y[x_index[0.5]] < 0.49
        assert 0.57 < dist.y[x_index[0.9]] < 0.65

class cumulative_intensity_distribution_python(object):
  # Prototype python version, superseded by faster C++ implementation
  # As described by  Howells, Phillips and Rogers, Acta Cryst. (1950). 3, 210

  def __init__(self, f_obs):
    self.info = f_obs.info()
    n_bins_used = f_obs.binner().n_bins_used()
    data = dict(zip(["%.2f" %(i/float(n_bins_used))
                     for i in range(0,n_bins_used)], [0]*n_bins_used))
    f_obs_sq = f_obs.f_as_f_sq()
    f_obs_sq.use_binner_of(f_obs)
    n_reflections = 0
    self.n_bins = f_obs_sq.binner().n_bins_all()
    self.mean_f_obs_sq = f_obs_sq.mean(use_binning=True)
    for intensity, d_spacing, indices in itertools.izip(
      f_obs_sq.data(), f_obs_sq.d_spacings().data(), f_obs_sq.indices()):
      n_reflections += 1
      i_over_mean_i = intensity/self._get_mean_f_obs_sq(d_spacing)
      rounded_i_over_mean_i = round(i_over_mean_i, 2)
      if i_over_mean_i > rounded_i_over_mean_i:
        rounded_i_over_mean_i += 0.01
      for i in range(n_bins_used,int(rounded_i_over_mean_i*n_bins_used)-1,-1):
        key = "%.2f" %(i/n_bins_used)
        if data.has_key(key):
          data[key] += 1
        else:
          continue

    xy_data = data.items()
    xy_data.sort()
    self.x = [float(x) for x, y in xy_data]
    self.y = [y/n_reflections for x, y in xy_data]

  def _get_mean_f_obs_sq(self, d_spacing):
    for n_bin in xrange(0,self.n_bins):
      if d_spacing >= self.mean_f_obs_sq.binner.bin_d_range(n_bin)[1]:
        break
    return self.mean_f_obs_sq.data[n_bin]

  def xy_plot_info(self):
    r = empty()
    r.title = "Cumulative Intensity Distribution"
    if (self.info != 0):
      r.title += ": " + str(self.info)
    r.x = self.x
    r.y = self.y
    r.xLegend = "z(%)"
    r.yLegend = "N(z)(%)"
    return r

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise_cumulative_intensity_distribution(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  exercise_sys_absent_intensity_distribution()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if __name__ == '__main__':
  run()
