from __future__ import division

import math
import sys

from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import adptbx, miller, sgtbx
from cctbx.array_family import flex

import libtbx.utils
from libtbx.test_utils import approx_equal

from scitbx.random import variate, normal_distribution#

from smtbx import absolute_structure


def exercise_hooft_analysis(space_group_info, d_min=1.0):
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
  # use gaussian errors
  g = variate(normal_distribution())
  errors = g(fc.size())
  fo2 = fo.as_intensity_array()
  fo2 = fo2.customized_copy(
    data=(fo2.data()+errors)*k,
    sigmas=flex.double(fc.size(), 1),
  )
  # first with the correct absolute structure
  analysis = absolute_structure.hooft_analysis(fo2, fc)
  assert approx_equal(analysis.hooft_y, 0, 1e-2)
  assert approx_equal(analysis.p2, 1)
  assert approx_equal(analysis.p3_true, 1)
  assert approx_equal(analysis.p3_false, 0)
  assert approx_equal(analysis.p3_racemic_twin, 0)
  NPP = absolute_structure.bijvoet_differences_probability_plot(analysis)
  assert approx_equal(NPP.correlation.coefficient(), 1, 0.04)
  assert approx_equal(NPP.fit.y_intercept(), 0, 0.1)
  assert approx_equal(NPP.correlation.coefficient(), 1, 0.002)
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
  assert approx_equal(NPP.correlation.coefficient(), 1, 0.04)
  assert approx_equal(NPP.fit.y_intercept(), 0, 0.1)
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
  #assert approx_equal(NPP.correlation.coefficient(), 1, 0.11)
  assert approx_equal(NPP.fit.y_intercept(), 0, 0.1)
  assert approx_equal(NPP.correlation.coefficient(), 1, 0.002)
  assert approx_equal(NPP.fit.y_intercept(), 0)

def run_call_back(flags, space_group_info):
  if not space_group_info.group().is_centric():
    exercise_hooft_analysis(space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if __name__ == '__main__':
  libtbx.utils.show_times_at_exit()
  run()
