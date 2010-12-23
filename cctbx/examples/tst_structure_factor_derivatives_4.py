from cctbx import xray
from cctbx import miller
from cctbx.examples.structure_factor_derivatives_4 import structure_factors
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
from iotbx.kriber import strudat
from libtbx.test_utils import approx_equal
import libtbx.load_env
import random
from itertools import count
from cStringIO import StringIO
import sys, os

random.seed(0)
flex.set_random_seed(0)

def scatterer_as_list(self):
  if (self.flags.use_u_iso_only()):
    return list(self.site) + [self.u_iso, self.occupancy, self.fp, self.fdp]
  return list(self.site) + list(self.u_star) \
       + [self.occupancy, self.fp, self.fdp]

def scatterer_from_list(l):
  if (len(l) == 7):
    return xray.scatterer(
      site=l[:3],
      u=l[3],
      occupancy=l[4],
      scattering_type="?",
      fp=l[5],
      fdp=l[6])
  return xray.scatterer(
    site=l[:3],
    u=l[3:9],
    occupancy=l[9],
    scattering_type="?",
    fp=l[10],
    fdp=l[11])

def d_target_d_params_finite(d_order, f_obs, xray_structure, eps=1.e-8):
  assert d_order in [1,2]
  result = flex.double()
  scatterers = xray_structure.scatterers()
  site_symmetry_table = xray_structure.site_symmetry_table()
  xray_structure_eps = xray_structure.deep_copy_scatterers()
  scatterers_eps = xray_structure_eps.scatterers()
  for i_scatterer,scatterer in enumerate(scatterers):
    if (not site_symmetry_table.is_special_position(i_scatterer)):
      site_symmetry_ops = None
      if (not scatterer.flags.use_u_aniso()):
        ips = range(7)
      else:
        ips = range(12)
    else:
      site_symmetry_ops = site_symmetry_table.get(i_scatterer)
      site_constraints = site_symmetry_ops.site_constraints()
      ips = list(site_constraints.independent_indices)
      if (not scatterer.flags.use_u_aniso()):
        ips.extend(range(3,7))
      else:
        adp_constraints = site_symmetry_ops.adp_constraints()
        ips.extend([i+3 for i in adp_constraints.independent_indices])
        ips.extend(range(9,12))
    dx = []
    for ip in ips:
      vs = []
      for signed_eps in [eps, -eps]:
        si_eps = scatterer_as_list(scatterer)
        si_eps[ip] += signed_eps
        if (site_symmetry_ops is not None):
          if (ip < 3):
            all_params = site_constraints.all_params(
              independent_params=site_constraints.independent_params(
              all_params=si_eps[:3]))
            si_eps = list(all_params) + si_eps[3:]
          elif (scatterer.flags.use_u_aniso() and ip < 9):
            all_params = adp_constraints.all_params(
              independent_params=adp_constraints.independent_params(
                all_params=si_eps[3:9]))
            si_eps = si_eps[:3] + list(all_params) + si_eps[9:]
        scatterers_eps[i_scatterer] = scatterer_from_list(si_eps)
        scatterers_eps[i_scatterer].scattering_type = scatterer.scattering_type
        xray_structure_eps.re_apply_symmetry(i_scatterer)
        sf = structure_factors(
          xray_structure=xray_structure_eps, miller_set=f_obs)
        if (d_order == 1):
          sum_target_f = 0
          for obs,f in zip(f_obs.data(), sf.fs()):
            target = least_squares(obs=obs, calc=f)
            sum_target_f += target.f()
          vs.append(sum_target_f)
        else:
          dp = sf.d_target_d_params(f_obs=f_obs, target_type=least_squares)
          vs.append(dp)
      diff = (vs[0]-vs[1])/(2*eps)
      if (d_order == 1):
        result.append(diff)
      else:
        result.extend(diff)
    scatterers_eps[i_scatterer] = scatterer
  return result

def compare_analytical_and_finite(
      f_obs,
      xray_structure,
      gradients_should_be_zero,
      eps,
      out):
  grads_fin = d_target_d_params_finite(
    d_order=1, f_obs=f_obs, xray_structure=xray_structure)
  print >> out, "grads_fin:", list(grads_fin)
  sf = structure_factors(
    xray_structure=xray_structure, miller_set=f_obs)
  grads_ana = sf.d_target_d_params(f_obs=f_obs, target_type=least_squares)
  print >> out, "grads_ana:", list(grads_ana)
  if (gradients_should_be_zero):
    flex.compare_derivatives(grads_ana, flex.double(grads_ana.size(), 0), eps)
  else:
    flex.compare_derivatives(grads_ana, grads_fin, eps)
  curvs_fin = d_target_d_params_finite(
    d_order=2, f_obs=f_obs, xray_structure=xray_structure)
  print >> out, "curvs_fin:", list(curvs_fin)
  curvs_ana = sf.d2_target_d_params(f_obs=f_obs, target_type=least_squares)
  print >> out, "curvs_ana:", list(curvs_ana)
  flex.compare_derivatives(curvs_ana.as_1d(), curvs_fin, eps)
  assert curvs_ana.matrix_is_symmetric(relative_epsilon=1e-10)
  print >> out
  #
  for i_method,curvs_method in enumerate([
        sf.d2_target_d_params_diag,
        sf.d2_target_d_params_diag_cpp]):
    curvs_diag_ana = curvs_method(f_obs=f_obs, target_type=least_squares)
    if (i_method != 0):
      flex.compare_derivatives(grads_ana, curvs_diag_ana.grads, eps=1e-12)
      curvs_diag_ana = curvs_diag_ana.curvs
    assert curvs_diag_ana.size() == curvs_ana.focus()[0]
    flex.compare_derivatives(
      curvs_ana.matrix_diagonal().as_1d(), curvs_diag_ana, eps=1e-12)
  #
  if (gradients_should_be_zero):
    return flex.max(flex.abs(grads_fin))

def exercise(
      xray_structure,
      anomalous_flag,
      max_n_indices,
      out):
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
  noise_fin = compare_analytical_and_finite(
    f_obs=f_obs,
    xray_structure=xray_structure,
    gradients_should_be_zero=True,
    eps=1.e-5,
    out=out)
  compare_analytical_and_finite(
    f_obs=f_obs.customized_copy(
      data=f_obs.data()*(flex.random_double(size=f_obs.size())+0.5)),
    xray_structure=xray_structure,
    gradients_should_be_zero=False,
    eps=max(1.e-5, noise_fin),
    out=out)

zeolite_edi = """\
*EDI
Code: EDI

P-4m2
 6.926 6.926 6.410 90.000 90.000 90.000
SI1 0.2680 0.0000 0.1184 4
SI2 0.0000 0.0000 0.5000 4
----------------------------------------
"""

def run_call_back(flags,
      space_group_info,
      max_n_indices=5,
      anomalous_flag=True):
  if (not flags.Verbose):
    out = StringIO()
  else:
    out = sys.stdout
  if (flags.chunk):
    chunk_n,chunk_i = [int(i) for i in flags.chunk.split(",")]
  else:
    chunk_n = 1
    chunk_i = 0
  if (flags.tag):
    if (flags.tag == "internal"):
      strudat_contents = strudat.read_all_entries(StringIO(zeolite_edi))
      strudat_entries = strudat_contents.entries
    else:
      atlas_file = libtbx.env.find_in_repositories(
        relative_path="phenix_regression/misc/strudat_zeolite_atlas",
        test=os.path.isfile)
      assert atlas_file is not None
      strudat_contents = strudat.read_all_entries(open(atlas_file))
      if (not isinstance(flags.tag, str)):
        strudat_entries = strudat_contents.entries
      else:
        strudat_entries = [strudat_contents.get(tag=flags.tag)]
        assert strudat_entries[0] is not None
  if (flags.isotropic):
    use_u_aniso_flags = [False]
  elif (flags.anisotropic):
    use_u_aniso_flags = [True]
  else:
    use_u_aniso_flags = [False, True]
  if (not flags.tag):
    for n_scatterers in xrange(2,3+1):
      for use_u_aniso in use_u_aniso_flags:
        xray_structure = random_structure.xray_structure(
          space_group_info=space_group_info,
          n_scatterers=n_scatterers,
          elements="random",
          volume_per_atom=100,
          general_positions_only=False,
          random_f_prime_d_min=1,
          random_f_double_prime=anomalous_flag,
          use_u_aniso = use_u_aniso,
          use_u_iso = not(use_u_aniso),
          random_u_iso=True,
          random_u_iso_scale=0.3,
          random_occupancy=True)
        exercise(
          xray_structure=xray_structure,
          anomalous_flag=anomalous_flag,
          max_n_indices=max_n_indices,
          out=out)
        out.flush()
  else:
    i_structure = count()
    for entry in strudat_entries:
      if (i_structure.next() % chunk_n != chunk_i): continue
      print >> sys.stderr, "strudat tag:", entry.tag
      sys.stderr.flush()
      print >> out, "strudat tag:", entry.tag
      out.flush()
      for use_u_aniso in use_u_aniso_flags:
        xray_structure = entry.as_xray_structure()
        xray_structure = random_structure.xray_structure(
          space_group_info=xray_structure.space_group_info(),
          unit_cell=xray_structure.unit_cell(),
          sites_frac=xray_structure.sites_frac(),
          elements="random",
          random_f_prime_d_min=1,
          random_f_double_prime=anomalous_flag,
          use_u_aniso = use_u_aniso,
          use_u_iso = not(use_u_aniso),
          random_u_iso=True,
          random_u_iso_scale=0.3,
          random_occupancy=True)
        exercise(
          xray_structure=xray_structure,
          anomalous_flag=anomalous_flag,
          max_n_indices=max_n_indices,
          out=out)
        out.flush()
  if (flags.tag):
    return False

def run(args):
  debug_utils.parse_options_loop_space_groups(args, run_call_back, (
    "chunk",
    "isotropic",
    "anisotropic",
    "tag"))

if (__name__ == "__main__"):
  run(sys.argv[1:])
