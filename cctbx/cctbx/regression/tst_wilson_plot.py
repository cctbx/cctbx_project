from cctbx import statistics
from cctbx import miller
from cctbx import adptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from scitbx.python_utils import dicts
import sys

def exercise(space_group_info, anomalous_flag,
             d_min=1.0, reflections_per_bin=200, n_bins=10, verbose=0):
  elements = ("N", "C", "C", "O") * 5
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0001,
    anisotropic_flag=00000,
    u_iso=adptbx.b_as_u(10)
    ).structure_factors(
        anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct")
  if (0 or verbose):
    structure_factors.xray_structure().show_summary()
  asu_contents = dicts.with_default_value(0)
  for elem in elements: asu_contents[elem] += 1
  f_calc = abs(structure_factors.f_calc())
  f_calc.setup_binner(
    auto_binning=0001,
    reflections_per_bin=reflections_per_bin,
    n_bins=n_bins)
  if (0 or verbose):
    f_calc.binner().show_summary()
  for k_given in [1,0.1,0.01,10,100]:
    f_obs = miller.array(
      miller_set=f_calc,
      data=f_calc.data()*k_given)
    f_obs.use_binner_of(f_calc)
    wp = statistics.wilson_plot(f_obs, asu_contents)
    if (0 or verbose):
      print "wilson_k, wilson_b:", wp.wilson_k, wp.wilson_b
    assert 0.8 < wp.wilson_k/k_given < 1.2
    assert 9 < wp.wilson_b < 11
    assert wp.xy_plot_info().fit_correlation == wp.fit_correlation

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
