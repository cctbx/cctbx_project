from __future__ import division

import math
import sys
from stdlib import random

from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import adptbx, miller, sgtbx
from cctbx.array_family import flex

from iotbx import csv_utils

import libtbx.utils
from libtbx.test_utils import approx_equal

from scitbx.random import variate, normal_distribution, gamma_distribution
from scitbx.math import distributions

from smtbx import absolute_structure

try:
  distributions.students_t_distribution(1)
except RuntimeError, e:
  # XXX Student's t distribution is not supported with GCC 3.2 builds
  if str(e).startswith("Implementation not available in this build."):
    students_t_available = False
    print "Skipping exercise_hooft_analysis() with Student's t distribution."
  else:
    raise RuntimeError(e)
else:
  students_t_available = True


def exercise_hooft_analysis(space_group_info, d_min=1.0,
                            use_students_t_errors=False,
                            debug=False):
  if use_students_t_errors and not students_t_available: return
  elements = ("N", "C", "C", "S") * 5
  xs = random_structure.xray_structure(
    space_group_info,
    elements=elements,
    volume_per_atom=20.,
    min_distance=1.5,
    general_positions_only=True,
    use_u_aniso=False,
    u_iso=adptbx.b_as_u(10),
  )
  xs.set_inelastic_form_factors(1.54, "sasaki")
  k = 0.05 + 10 * flex.random_double()
  fc = xs.structure_factors(
        anomalous_flag=True, d_min=d_min, algorithm="direct").f_calc()
  fo = fc.as_amplitude_array()
  fo.set_observation_type_xray_amplitude()
  if use_students_t_errors:
    nu = random.uniform(1, 20)
    normal_g = variate(normal_distribution())
    gamma_g = variate(gamma_distribution(0.5*nu, 2))
    errors = normal_g(fc.size())/flex.sqrt(2*gamma_g(fc.size()))
  else:
    # use gaussian errors
    g = variate(normal_distribution())
    errors = g(fc.size())
  fo2 = fo.as_intensity_array()
  fo2 = fo2.customized_copy(
    data=(fo2.data()+errors)*k,
    sigmas=flex.double(fc.size(), 1),
  )
  #
  if debug:
    distribution = distributions.normal_distribution()
    observed_deviations = (fo2.data() - k*fc.as_intensity_array().data())
    observed_deviations = observed_deviations.select(
      flex.sort_permutation(observed_deviations))
    expected_deviations = distribution.quantiles(observed_deviations.size())
    csv_utils.writer(
      open('delta_F_npp.csv', 'wb'), (expected_deviations, observed_deviations))
  # first with the correct absolute structure
  analysis = absolute_structure.hooft_analysis(fo2, fc)
  assert approx_equal(analysis.hooft_y, 0, 1e-2)
  assert approx_equal(analysis.p2, 1)
  assert approx_equal(analysis.p3_true, 1)
  assert approx_equal(analysis.p3_false, 0)
  assert approx_equal(analysis.p3_racemic_twin, 0)
  NPP = absolute_structure.bijvoet_differences_probability_plot(analysis)
  if use_students_t_errors:
    tPP = absolute_structure.bijvoet_differences_probability_plot(
      analysis, use_students_t_distribution=True)
    if tPP.distribution.degrees_of_freedom() < 100:
      tPP.correlation.coefficient() > NPP.correlation.coefficient()
  else:
    assert approx_equal(NPP.correlation.coefficient(), 1, 0.005)
  if debug:
    csv_utils.writer(open('npp.csv', 'wb'), (NPP.x,NPP.y))
    if use_students_t_errors:
      csv_utils.writer(open('tpp.csv', 'wb'), (tPP.x,tPP.y))
  assert approx_equal(NPP.fit.y_intercept(), 0)
  # and now with the wrong absolute structure
  xs_i = xs.inverse_hand()
  fc_i = xs_i.structure_factors(
    anomalous_flag=True, d_min=d_min, algorithm="direct").f_calc()
  analysis = absolute_structure.hooft_analysis(fo2, fc_i)
  assert approx_equal(analysis.hooft_y, 1, 1e-2)
  assert approx_equal(analysis.p2, 0)
  assert approx_equal(analysis.p3_true, 0)
  assert approx_equal(analysis.p3_false, 1)
  assert approx_equal(analysis.p3_racemic_twin, 0)
  NPP = absolute_structure.bijvoet_differences_probability_plot(analysis)
  if use_students_t_errors:
    tPP = absolute_structure.bijvoet_differences_probability_plot(
      analysis, use_students_t_distribution=True)
    if tPP.distribution.degrees_of_freedom() < 100:
      assert tPP.correlation.coefficient() > NPP.correlation.coefficient()
  else:
    assert approx_equal(NPP.correlation.coefficient(), 1, 0.002)
    assert approx_equal(NPP.fit.y_intercept(), 0)
  # test for the case of a racemic twin
  fo2_twin = fc.customized_copy(
    data=fc.data()+fc_i.data()).as_intensity_array()
  fo2_twin = fo2_twin.customized_copy(
    data=(errors + fo2_twin.data()) * k,
    sigmas=fo2.sigmas())
  analysis = absolute_structure.hooft_analysis(fo2_twin, fc)
  assert approx_equal(analysis.hooft_y, 0.5, 1e-2)
  assert approx_equal(analysis.p3_true, 0)
  assert approx_equal(analysis.p3_false, 0)
  assert approx_equal(analysis.p3_racemic_twin, 1)
  NPP = absolute_structure.bijvoet_differences_probability_plot(analysis)
  if use_students_t_errors:
    tPP = absolute_structure.bijvoet_differences_probability_plot(
      analysis, use_students_t_distribution=True)
    if tPP.distribution.degrees_of_freedom() < 100:
      assert tPP.correlation.coefficient() > NPP.correlation.coefficient()
  else:
    assert approx_equal(NPP.correlation.coefficient(), 1, 0.002)
    assert approx_equal(NPP.fit.y_intercept(), 0)

def run_call_back(flags, space_group_info):
  if not space_group_info.group().is_centric():
    for use_students_t_errors in (True, False):
      exercise_hooft_analysis(space_group_info,
                              use_students_t_errors=use_students_t_errors,
                              debug=flags.Debug)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if __name__ == '__main__':
  libtbx.utils.show_times_at_exit()
  run()
