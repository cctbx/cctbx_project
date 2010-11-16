from cctbx import miller
from cctbx.examples.structure_factor_derivatives_3 \
  import scatterer_as_list, scatterer_from_list, structure_factors
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
from libtbx.test_utils import approx_equal
import random
from cStringIO import StringIO
import sys

random.seed(0)
flex.set_random_seed(0)

def d_target_d_params_finite(f_obs, xray_structure, eps=1.e-8):
  result = flex.double()
  scatterers = xray_structure.scatterers()
  xray_structure_eps = xray_structure.deep_copy_scatterers()
  scatterers_eps = xray_structure_eps.scatterers()
  for i_scatterer in xrange(len(scatterers)):
    if (not scatterers[i_scatterer].flags.use_u_aniso()):
      np = 7
    else:
      np = 12
    dx = []
    for ix in xrange(np):
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
  for i_scatterer in xrange(len(scatterers)):
    if (not scatterers[i_scatterer].flags.use_u_aniso()):
      np = 7
    else:
      np = 12
    dx = []
    for ix in xrange(np):
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

def compare_derivatives(ana, fin):
  s = max(1, flex.max(flex.abs(ana)))
  assert approx_equal(ana/s, fin/s)

def compare_analytical_and_finite(f_obs, xray_structure, out):
  grads_fin = d_target_d_params_finite(
    f_obs=f_obs, xray_structure=xray_structure)
  print >> out, "grads_fin:", list(grads_fin)
  sf = structure_factors(
    xray_structure=xray_structure, miller_set=f_obs)
  grads_ana = sf.d_target_d_params(f_obs=f_obs, target_type=least_squares)
  print >> out, "grads_ana:", list(grads_ana)
  compare_derivatives(grads_ana, grads_fin)
  curvs_fin = d2_target_d_params_finite(
    f_obs=f_obs, xray_structure=xray_structure)
  print >> out, "curvs_fin:", list(curvs_fin)
  curvs_ana = sf.d2_target_d_params(f_obs=f_obs, target_type=least_squares)
  print >> out, "curvs_ana:", list(curvs_ana)
  compare_derivatives(curvs_ana, curvs_fin)
  print >> out

def exercise(
      space_group_info,
      use_u_aniso,
      anomalous_flag,
      max_n_indices=5,
      verbose=0):
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  for n_scatterers in xrange(3,3+1):
    for i_trial in xrange(1):
      xray_structure = random_structure.xray_structure(
        space_group_info=space_group_info,
        elements=["const"]*n_scatterers,
        volume_per_atom=100,
        general_positions_only=True,
        random_f_prime_d_min=1,
        random_f_double_prime=anomalous_flag,
        use_u_aniso = use_u_aniso,
        use_u_iso = (not use_u_aniso),
        random_u_iso=True,
        random_u_iso_scale=0.3,
        random_occupancy=True)
      xray_structure.show_summary(f=out).show_scatterers(f=out)
      miller_set = miller.build_set(
        crystal_symmetry=xray_structure,
        anomalous_flag=anomalous_flag,
        d_min=max(1, min(xray_structure.unit_cell().parameters()[:3])/2.5))
      n_indices = miller_set.indices().size()
      if (n_indices > max_n_indices):
        miller_set = miller_set.select(
          flex.random_size_t(size=max_n_indices) % n_indices)
      sf = structure_factors(
        xray_structure=xray_structure,
        miller_set=miller_set)
      f_calc = miller_set.structure_factors_from_scatterers(
        xray_structure=xray_structure,
        algorithm="direct",
        cos_sin_table=False).f_calc()
      f_calc.show_summary(f=out)
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

def run_call_back(flags, space_group_info):
  if (flags.isotropic):
    use_u_aniso = [False]
  elif (flags.anisotropic):
    use_u_aniso = [True]
  else:
    use_u_aniso_flags = [False, True]
  for use_u_aniso in use_u_aniso_flags:
    exercise(
      space_group_info=space_group_info,
      use_u_aniso=use_u_aniso,
      anomalous_flag=True,
      verbose=flags.Verbose)

def run(args):
  debug_utils.parse_options_loop_space_groups(args, run_call_back, (
    "isotropic",
    "anisotropic"))

if (__name__ == "__main__"):
  run(sys.argv[1:])
