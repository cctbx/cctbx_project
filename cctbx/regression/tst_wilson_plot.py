from __future__ import absolute_import, division, print_function
from cctbx import statistics
from cctbx import miller
from cctbx import adptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from libtbx.test_utils import show_diff
import sys
from cctbx.sgtbx import space_group_info
from libtbx.test_utils import approx_equal

def exercise(space_group_info, anomalous_flag,
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
  asu_contents = dicts.with_default_value(0)
  for elem in elements: asu_contents[elem] += 1
  f_calc = abs(structure_factors.f_calc())
  f_calc.setup_binner(
    auto_binning=True,
    reflections_per_bin=reflections_per_bin,
    n_bins=n_bins)
  if (0 or verbose):
    f_calc.binner().show_summary()
  for k_given in [1,0.1,0.01,10,100]:
    f_obs = miller.array(
      miller_set=f_calc,
      data=f_calc.data()*k_given).set_observation_type_xray_amplitude()
    f_obs.use_binner_of(f_calc)
    wp = statistics.wilson_plot(f_obs, asu_contents, scattering_table="wk1995",
      e_statistics=True)
    if (0 or verbose):
      print("wilson_k, wilson_b:", wp.wilson_k, wp.wilson_b)
      print("space group:", space_group_info.group().type().hall_symbol())
      print("<E^2-1>:", wp.mean_e_sq_minus_1)

    assert 0.8 < wp.wilson_k/k_given < 1.2
    assert 0.64 < wp.wilson_intensity_scale_factor/(k_given*k_given) < 1.44
    assert 9 < wp.wilson_b < 11
    assert wp.xy_plot_info().fit_correlation == wp.fit_correlation
    if space_group_info.group().is_centric():
      assert 0.90 < wp.mean_e_sq_minus_1 < 1.16
      assert 3.15 < wp.percent_e_sq_gt_2 < 6.5
    else:
      assert 0.65 < wp.mean_e_sq_minus_1 < 0.90
      assert 1.0 < wp.percent_e_sq_gt_2 < 3.15
    assert wp.normalised_f_obs.size() == f_obs.size()
  f_obs = f_calc.array(data=flex.double(f_calc.indices().size(), 0))
  f_obs.use_binner_of(f_calc)
  n_bins = f_obs.binner().n_bins_used()
  try:
    statistics.wilson_plot(f_obs, asu_contents, scattering_table="wk1995")
  except RuntimeError as e:
    assert not show_diff(str(e), """\
wilson_plot error: %d empty bins:
  Number of bins: %d
  Number of f_obs > 0: 0
  Number of f_obs <= 0: %d""" % (n_bins, n_bins, f_obs.indices().size()))

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

def run2(n = 100, d_min = 1.0):
  for b in range(10, 110, 10):
    res = []
    for st in ["it1992", "neutron", "electron"]:
      elements = ('C', 'N', 'O', 'H')*n
      group = space_group_info("P1")
      xrs = random_structure.xray_structure(
        space_group_info = group,
        volume_per_atom = 25.,
        general_positions_only = True,
        elements = elements,
        u_iso = adptbx.b_as_u(b),
        min_distance = 1.5)
      xrs.scattering_type_registry(
        table = st,
        d_min = d_min,
        types_without_a_scattering_contribution=["?"])
      fo = abs(xrs.structure_factors(d_min=d_min).f_calc())
      fo.setup_binner(reflections_per_bin=100)
      asu_contents = {}
      for e in elements:
        asu_contents[e] = n
      result = statistics.wilson_plot(f_obs = fo, asu_contents = asu_contents,
        scattering_table = st)
      res.append([st,round(result.wilson_b, 0)])
      assert approx_equal(result.wilson_b, b, 1)
    print(b, res)

if (__name__ == "__main__"):
  run()
  run2()
