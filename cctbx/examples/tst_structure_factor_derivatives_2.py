from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx.examples.structure_factor_derivatives_2 \
  import scatterer_as_list, scatterer_from_list, structure_factors
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

random.seed(0)
flex.set_random_seed(0)

def d_target_d_params_finite(f_obs, xray_structure, eps=1.e-8):
  result = flex.double()
  scatterers = xray_structure.scatterers()
  xray_structure_eps = xray_structure.deep_copy_scatterers()
  scatterers_eps = xray_structure_eps.scatterers()
  for i_scatterer in range(len(scatterers)):
    dx = []
    for ix in range(7):
      ts = []
      for signed_eps in [eps, -eps]:
        si_eps = scatterer_as_list(scatterers[i_scatterer])
        si_eps[ix] += signed_eps
        scatterers_eps[i_scatterer] = scatterer_from_list(si_eps)
        sf = structure_factors(
          xray_structure=xray_structure_eps, miller_set=f_obs)
        sum_target_f = 0
        for obs,f in zip(f_obs.data(), sf.fs()):
          target = least_squares(obs=obs, calc=f)
          sum_target_f += target.f()
        ts.append(sum_target_f)
      result.append((ts[0]-ts[1])/(2*eps))
    scatterers_eps[i_scatterer] = scatterers[i_scatterer]
  return result

def d2_target_d_params_finite(f_obs, xray_structure, eps=1.e-8):
  result = flex.double()
  scatterers = xray_structure.scatterers()
  xray_structure_eps = xray_structure.deep_copy_scatterers()
  scatterers_eps = xray_structure_eps.scatterers()
  for i_scatterer in range(len(scatterers)):
    for ix in range(7):
      gs = []
      for signed_eps in [eps, -eps]:
        si_eps = scatterer_as_list(scatterers[i_scatterer])
        si_eps[ix] += signed_eps
        scatterers_eps[i_scatterer] = scatterer_from_list(si_eps)
        sf = structure_factors(
          xray_structure=xray_structure_eps, miller_set=f_obs)
        dp = sf.d_target_d_params(f_obs=f_obs, target_type=least_squares)
        gs.append(dp)
      result.extend((gs[0]-gs[1])/(2*eps))
    scatterers_eps[i_scatterer] = scatterers[i_scatterer]
  return result

def compare_analytical_and_finite(f_obs, xray_structure, out):
  grads_fin = d_target_d_params_finite(
    f_obs=f_obs, xray_structure=xray_structure)
  print("grads_fin:", list(grads_fin), file=out)
  sf = structure_factors(
    xray_structure=xray_structure, miller_set=f_obs)
  grads_ana = sf.d_target_d_params(f_obs=f_obs, target_type=least_squares)
  print("grads_ana:", list(grads_ana), file=out)
  assert approx_equal(grads_ana, grads_fin)
  curvs_fin = d2_target_d_params_finite(
    f_obs=f_obs, xray_structure=xray_structure)
  print("curvs_fin:", list(curvs_fin), file=out)
  curvs_ana = sf.d2_target_d_params(f_obs=f_obs, target_type=least_squares)
  print("curvs_ana:", list(curvs_ana), file=out)
  assert approx_equal(curvs_ana, curvs_fin, 1.e-5)
  print(file=out)

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  crystal_symmetry = crystal.symmetry(
    unit_cell=(8,9,10,83,97,113),
    space_group_symbol="P1")
  miller_set = miller.set(
    crystal_symmetry=crystal_symmetry,
    indices=flex.miller_index([(1,2,3), (2,3,4), (-1,3,-2)]),
    anomalous_flag=False)
  for n_scatterers in range(2,2+5):
    for i_trial in range(5):
      scatterers = flex.xray_scatterer()
      for i in range(n_scatterers):
        scatterers.append(xray.scatterer(
          site=[random.random() for i in range(3)],
          u=random.random()*0.1,
          occupancy=random.random(),
          scattering_type="const",
          fp=(random.random()-0.5)*2,
          fdp=(random.random()-0.5)*2))
      xray_structure = xray.structure(
        crystal_symmetry=crystal_symmetry,
        scatterers=scatterers)
      sf = structure_factors(
        xray_structure=xray_structure,
        miller_set=miller_set)
      f_calc = miller_set.structure_factors_from_scatterers(
        xray_structure=xray_structure,
        algorithm="direct",
        cos_sin_table=False).f_calc()
      assert approx_equal(sf.fs(), f_calc.data())
      f_obs = miller_set.array(data=flex.abs(sf.fs()))
      compare_analytical_and_finite(
        f_obs=f_obs,
        xray_structure=xray_structure,
        out=out)
      compare_analytical_and_finite(
        f_obs=f_obs.customized_copy(
          data=f_obs.data()*(flex.random_double(size=f_obs.size())+0.5)),
        xray_structure=xray_structure,
        out=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
