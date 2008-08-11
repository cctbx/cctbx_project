from cctbx import xray
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

def run():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  f_calc = miller_set.array(
    data=flex.polar(
      flex.random_double(miller_set.size())*10-5,
      flex.random_double(miller_set.size())*10-5))
  obs = miller_set.array(
    data=flex.norm(f_calc.data()) + (flex.random_double(miller_set.size())*2-1),
    sigmas=flex.random_double(miller_set.size()))
  obs.set_observation_type_xray_intensity()

  exercise_shelx_weighting(f_calc, obs)
  exercise_quasi_unit_weighting(obs)

  print 'OK'


def exercise_shelx_weighting(f_calc, obs):
  a,b = 0, 0
  weighting = xray.weighting_schemes.simple_shelx_weighting(a,b)
  weighting_ref = xray.weighting_schemes.pure_statistical_weighting()
  for w in (weighting, weighting_ref):
    w.calculated = f_calc
    w.observed = obs
    w.compute()
  assert approx_equal(weighting.weights, weighting_ref.weights)

  a,b = 10, 100
  weighting = xray.weighting_schemes.shelx_weighting(a,b)
  weighting_ref = xray.weighting_schemes.simple_shelx_weighting(a,b)
  for w in (weighting, weighting_ref):
    w.calculated = f_calc
    w.observed = obs
    w.compute()
  assert approx_equal(weighting.weights, weighting_ref.weights)
  assert weighting.derivatives_wrt_f_c is None
  for w in (weighting, weighting_ref):
    w.computing_derivatives_wrt_f_c = True
    w.compute()
  assert approx_equal(weighting.derivatives_wrt_f_c,
                      weighting_ref.derivatives_wrt_f_c)

  weighting.observed = weighting.observed.discard_sigmas()
  weighting.compute()

def exercise_quasi_unit_weighting(obs):
  w = xray.weighting_schemes.intensity_quasi_unit_weighting()
  w.observed = obs.discard_sigmas()
  w.compute()
  assert approx_equal(list(w.weights), list(1/(4.*obs.data())))

if __name__ == '__main__':
  run()
