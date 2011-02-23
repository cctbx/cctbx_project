from __future__ import division
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.adptbx import b_as_u
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.lbfgsb
from scitbx import matrix
import libtbx.phil.command_line
from libtbx import easy_pickle
import random
import sys

if (1):
  random.seed(0)
  flex.set_random_seed(0)

class curv_filter(object):

  def __init__(O, curvs, lim_eps):
    c_pos = curvs.select(curvs > 0)
    O.c_lim = flex.max_default(c_pos, default=0) * lim_eps
    if (O.c_lim == 0):
      O.c_lim = 1
      O.c_rms = 1
    else:
      O.c_rms = flex.mean_sq(c_pos)**0.5
    O.n_below_limit = 0
    O.n_above_limit = 0

  def apply(O, some_curvs):
    result = flex.double()
    for c in some_curvs:
      if (c < O.c_lim):
        c = O.c_rms
        O.n_below_limit += 1
      else:
        O.n_above_limit += 1
      result.append(c)
    return result

def get_curvs_work(f_obs, weights, xray_structure, lim_eps):
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=xray_structure,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=False,
    obs_type="F",
    obs=abs(f_calc).data(),
    weights=weights,
    r_free_flags=None,
    f_calc=f_calc.data(),
    derivatives_depth=2,
    scale_factor=1)
  gact = xray_structure.grads_and_curvs_target_simple(
    miller_indices=f_obs.indices(),
    da_db=ls.gradients_work(),
    daa_dbb_dab=ls.hessian_work())
  c_all = gact.curvs
  c_active_site = flex.double()
  c_active_u_iso = flex.double()
  i_all = 0
  sstab = xray_structure.site_symmetry_table()
  for i_sc,sc in enumerate(xray_structure.scatterers()):
    assert sc.flags.use_u_iso()
    assert not sc.flags.use_u_aniso()
    site_symmetry = sstab.get(i_sc)
    if (site_symmetry.is_point_group_1()):
      np = 3
    else:
      np = site_symmetry.site_constraints().n_independent_params()
    c_active_site.extend(c_all[i_all:i_all+np])
    c_active_u_iso.append(c_all[i_all+np])
    np += 4 # u_iso, occ, fp, fdp
    i_all += np
  assert i_all == c_all.size()
  #
  cf_site = curv_filter(curvs=c_active_site, lim_eps=lim_eps)
  cf_u_iso = curv_filter(curvs=c_active_u_iso, lim_eps=lim_eps)
  #
  result = flex.double()
  i_all = 0
  sstab = xray_structure.site_symmetry_table()
  for i_sc,sc in enumerate(xray_structure.scatterers()):
    assert sc.flags.use_u_iso()
    assert not sc.flags.use_u_aniso()
    site_symmetry = sstab.get(i_sc)
    if (site_symmetry.is_point_group_1()):
      np = 3
    else:
      np = site_symmetry.site_constraints().n_independent_params()
    result.extend(cf_site.apply(c_all[i_all:i_all+np]))
    result.extend(cf_u_iso.apply(c_all[i_all+np:i_all+np+1]))
    np += 4 # u_iso, occ, fp, fdp
    i_all += np
  assert i_all == c_all.size()
  print "curv site below limit: %d of %d" % (
    cf_site.n_below_limit,
    cf_site.n_below_limit + cf_site.n_above_limit)
  print "curv u_iso below limit: %d of %d" % (
    cf_u_iso.n_below_limit,
    cf_u_iso.n_below_limit + cf_u_iso.n_above_limit)
  sys.stdout.flush()
  return result

class refinement_stats(libtbx.slots_getstate_setstate):

  __slots__ = ["iter", "nfun", "f", "gnorm", "crmsd", "armsd", "urmsd"]

  def __init__(O, **kw):
    for slot in O.__slots__:
      setattr(O, slot, kw[slot])

class ls_refinement(object):

  def __init__(O,
        f_obs,
        xray_structure,
        params,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None,
        reference_structure=None):
    O.f_obs = f_obs
    O.weights = flex.double(f_obs.data().size(), 1)
    O.xray_structure = xray_structure
    O.params = params
    O.reference_structure = reference_structure
    O.curvs_work = get_curvs_work(
      f_obs=O.f_obs,
      weights=O.weights,
      xray_structure=O.xray_structure,
      lim_eps=O.params.curv_filter_lim_eps)
    if (O.params.curvature_rescaling):
      O.p_as_x = flex.sqrt(O.curvs_work)
    else:
      O.p_as_x = flex.double(O.curvs_work.size(), 1)
    O.pack_parameters()
    O.number_of_function_evaluations = -1
    O.number_of_lbfgs_iterations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.history = []
    O.callback_after_step(minimizer=None)
    if (params.minimizer == "lbfgs"):
      O.minimizer = scitbx.lbfgs.run(
        target_evaluator=O,
        termination_params=lbfgs_termination_params,
        exception_handling_params=lbfgs_exception_handling_params)
    elif (params.minimizer == "lbfgs_raw"):
      O.run_lbfgs_raw()
    elif (params.minimizer == "lbfgsb"):
      O.run_lbfgsb()
    O.f_final, O.g_final = O.compute_functional_and_gradients()

  def pack_parameters(O):
    O.x = flex.double()
    O.l = flex.double()
    O.u = flex.double()
    O.nbd = flex.int()
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      assert sc.flags.use_u_iso()
      assert not sc.flags.use_u_aniso()
      #
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        p = sc.site
      else:
        p = site_symmetry.site_constraints().independent_params(
          all_params=sc.site)
      O.x.extend(flex.double(p))
      O.l.resize(O.x.size(), 0)
      O.u.resize(O.x.size(), 0)
      O.nbd.resize(O.x.size(), 0)
      #
      O.x.append(sc.u_iso)
      O.l.append(0)
      O.u.append(0)
      O.nbd.append(1)
    O.x *= O.p_as_x

  def unpack_parameters(O):
    p = O.x / O.p_as_x
    ix = 0
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        sc.site = tuple(p[ix:ix+3])
        ix += 3
      else:
        constr = site_symmetry.site_constraints()
        np = constr.n_independent_params()
        sc.site = constr.all_params(independent_params=tuple(p[ix:ix+np]))
        ix += np
      sc.u_iso = p[ix]
      ix += 1
    assert ix == O.x.size()

  def adjust_stp(O, stp, csd,
        large_shift_site=0.1,
        large_shift_u_iso_factor=0.5,
        min_large_shift_u_iso=b_as_u(1)):
    print "adjust_stp", stp
    assert csd.size() == O.x.size()
    max_stp = 100
    ix = 0
    uc = O.xray_structure.unit_cell()
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        np = 3
        for i,ucp in enumerate(uc.parameters()[:3]):
          lsx = large_shift_site / ucp \
              * O.p_as_x[ix+i]
          if (abs(stp*csd[ix+i]) > lsx):
            max_stp = min(max_stp, abs(lsx/csd[ix]))
      else:
        constr = site_symmetry.site_constraints()
        np = constr.n_independent_params()
        raise RuntimeError("SPECIAL POSITION LARGE SHIFT NOT IMPLEMENTED.")
      ix += np
      u_iso = O.x[ix]/O.p_as_x[ix]
      lux = max(min_large_shift_u_iso, abs(u_iso) * large_shift_u_iso_factor) \
          * O.p_as_x[ix]
      if (abs(stp*csd[ix]) > lux):
        max_stp = min(max_stp, abs(lux/csd[ix]))
      ix += 1
    assert ix == O.x.size()
    print "max_stp:", max_stp
    sys.stdout.flush()
    if (O.params.use_max_stp):
      return min(max_stp, stp)
    return stp

  def compute_functional_and_gradients(O):
    O.number_of_function_evaluations += 1
    O.unpack_parameters()
    try:
      f_calc = O.f_obs.structure_factors_from_scatterers(
        xray_structure=O.xray_structure,
        algorithm="direct",
        cos_sin_table=False).f_calc()
    except RuntimeError, e:
      print "RuntimeError f_calc:", e
      sys.stdout.flush()
      raise RuntimeError("f_calc")
    ls = xray.targets_least_squares(
      compute_scale_using_all_data=False,
      obs_type="F",
      obs=O.f_obs.data(),
      weights=O.weights,
      r_free_flags=None,
      f_calc=f_calc.data(),
      derivatives_depth=2,
      scale_factor=1)
    gact = O.xray_structure.grads_and_curvs_target_simple(
      miller_indices=O.f_obs.indices(),
      da_db=ls.gradients_work(),
      daa_dbb_dab=ls.hessian_work())
    g_all = gact.grads
    g_active = flex.double()
    i_all = 0
    sstab = O.xray_structure.site_symmetry_table()
    for i_sc,sc in enumerate(O.xray_structure.scatterers()):
      assert sc.flags.use_u_iso()
      assert not sc.flags.use_u_aniso()
      site_symmetry = sstab.get(i_sc)
      if (site_symmetry.is_point_group_1()):
        np = 3
      else:
        np = site_symmetry.site_constraints().n_independent_params()
      g_active.extend(g_all[i_all:i_all+np+1])
      np += 4 # u_iso, occ, fp, fdp
      i_all += np
    assert i_all == g_all.size()
    assert g_active.size() == O.x.size()
    O.f_last = ls.target_work()
    O.g_last = g_active / O.p_as_x
    return O.f_last, O.g_last

  def callback_after_step(O, minimizer, suffix=""):
    O.number_of_lbfgs_iterations += 1
    O.callback_after_step_no_counting(suffix=suffix)

  def callback_after_step_no_counting(O, suffix=""):
    if (O.number_of_lbfgs_iterations % 10 == 0):
      s = "step  fun    f        |g|"
      if (O.reference_structure is not None):
        s += "     cRMSD aRMSD uRMSD"
      print s
    s = "%4d %4d %9.2e %9.2e" % (
      O.number_of_lbfgs_iterations,
      O.number_of_function_evaluations,
      O.f_last, O.g_last.norm())
    s += O.format_rms_info()
    print s+suffix
    O.history.append(refinement_stats(
      iter=O.number_of_lbfgs_iterations,
      nfun=O.number_of_function_evaluations,
      f=O.f_last,
      gnorm=O.g_last.norm(),
      crmsd=O.crmsd,
      armsd=O.armsd,
      urmsd=O.urmsd))
    sys.stdout.flush()

  def get_rms_info(O):
    if (O.reference_structure is None):
      return None
    xs = O.xray_structure
    rs = O.reference_structure
    xf = xs.sites_frac()
    rf = rs.sites_frac()
    # TODO: use scattering power as weights, move to method of xray.structure
    ave_csh = matrix.col((xf-rf).mean())
    ave_csh_perp = matrix.col(xs.space_group_info()
      .subtract_continuous_allowed_origin_shifts(translation_frac=ave_csh))
    caosh_corr = ave_csh - ave_csh_perp
    omx = xs.unit_cell().orthogonalization_matrix()
    O.crmsd = (omx * (rf - xf)).rms_length()
    O.armsd = (omx * (rf - xf + caosh_corr)).rms_length()
    O.urmsd = flex.mean_sq(
        xs.scatterers().extract_u_iso()
      - rs.scatterers().extract_u_iso())**0.5
    return (O.crmsd, O.armsd, O.urmsd)

  def format_rms_info(O):
    s = ""
    info = O.get_rms_info()
    if (info is not None):
      for r in info:
        s += " %5.3f" % r
    return s

  def show_rms_info(O):
    s = O.format_rms_info()
    if (len(s) != 0):
       print " cRMSD aRMSD uRMSD"
       print s
       sys.stdout.flush()

  def run_lbfgs_raw(O):
    lbfgs_impl = [
      scitbx.lbfgs.raw_reference,
      scitbx.lbfgs.raw][O.params.lbfgs_impl_switch]
    diagco = O.params.diagco
    assert diagco in [0,1,2]
    n = O.x.size()
    m = 5
    iprint = O.params.iprint
    eps = 1.0e-5
    xtol = 1.0e-16
    size_w = n*(2*m+1)+2*m
    w = flex.double(size_w)
    diag = None
    diag0 = None
    iflag = 0
    while True:
      if (iflag in [0,1]):
        f, g = O.compute_functional_and_gradients()
        if (iflag == 0):
          if (diagco == 0):
            diag = flex.double(n, -1e20)
          else:
            assert O.curvs_work.size() == O.x.size()
            assert O.curvs_work.all_gt(0)
            if (O.params.curvature_rescaling):
              diag0 = flex.double(n, 1)
            else:
              diag0 = 1 / O.curvs_work
            diag = diag0.deep_copy()
      elif (iflag == 2):
        diag.clear()
        diag.extend(diag0)
      elif (iflag == 100):
        new_stp = O.adjust_stp(
          stp=lbfgs_impl.stp(),
          csd=lbfgs_impl.current_search_direction())
        lbfgs_impl.set_stp(value=new_stp)
      else:
        raise RuntimeError("invalid iflag value: %d" % iflag)
      if (iflag in [0,1]):
        O.show_rms_info()
      iflag = lbfgs_impl(
        n=n, m=m, x=O.x, f=f, g=g, diagco=diagco, diag=diag,
        iprint=iprint, eps=eps, xtol=xtol, w=w, iflag=iflag)
      if (iflag <= 0): break
    if (O.params.lbfgs_impl_switch == 1):
      print "iter, nfun:", lbfgs_impl.iter(), lbfgs_impl.nfun()
      assert O.number_of_function_evaluations == lbfgs_impl.nfun()
      O.number_of_lbfgs_iterations = lbfgs_impl.iter()
    O.callback_after_step_no_counting(suffix=" FINAL")

  def run_lbfgsb(O, iprint=1):
    n = O.x.size()
    minimizer = scitbx.lbfgsb.minimizer(
      n=n,
      m=5,
      l=O.l,
      u=O.u,
      nbd=O.nbd,
      enable_stp_init=True,
      factr=1.0e+7,
      pgtol=1.0e-5,
      iprint=iprint)
    f, g = -1e20, flex.double(n, -1e20)
    while True:
      if (minimizer.process(O.x, f, g)):
        f, g = O.compute_functional_and_gradients()
      elif (minimizer.requests_stp_init()):
        new_stp = O.adjust_stp(
          stp=minimizer.relative_step_length_line_search(),
          csd=minimizer.current_search_direction())
        minimizer.set_relative_step_length_line_search(value=new_stp)
      elif (minimizer.is_terminated()):
        O.callback_after_step_no_counting(suffix=" FINAL")
        break
      else:
        O.callback_after_step(minimizer=None)

def run_refinement(structure_ideal, structure_shake, params, run_id):
  print "Ideal structure:"
  structure_ideal.show_summary().show_scatterers()
  print
  print "Modified structure:"
  structure_shake.show_summary().show_scatterers()
  print
  print "rms difference:", \
    structure_ideal.rms_difference(other=structure_shake)
  print
  print "structure_shake inter-atomic distances:"
  structure_shake.show_distances(distance_cutoff=4)
  print
  f_obs = abs(structure_ideal.structure_factors(
    anomalous_flag=False,
    d_min=1,
    algorithm="direct",
    cos_sin_table=False).f_calc())
  try:
    return ls_refinement(
      f_obs=f_obs,
      xray_structure=structure_shake,
      params=params,
      reference_structure=structure_ideal)
  except RuntimeError, e:
    print "RuntimeError run_id:", run_id
    sys.stdout.flush()
    if (str(e) != "f_calc"):
      raise

def run_call_back(flags, space_group_info, params):
  structure_shake = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "O", "S", "Yb"),
    volume_per_atom=200,
    min_distance=2.0,
    general_positions_only=params.general_positions_only,
    random_u_iso=True)
  structure_ideal = structure_shake.deep_copy_scatterers()
  structure_shake.shake_sites_in_place(rms_difference=params.shake_sites_rmsd)
  structure_shake.shake_adp(spread=params.shake_adp_spread)
  #
  run_id = ""
  if (params.pickle_root_name is not None):
    run_id += params.pickle_root_name + "_"
  run_id += str(space_group_info).replace(" ","").replace("/","_").lower()
  if (params.pickle_root_name is not None):
    pickle_file_name = run_id + "_ideal_shake.pickle"
    print "writing file:", pickle_file_name
    easy_pickle.dump(
      file_name=pickle_file_name,
      obj=(structure_ideal, structure_shake))
    print
    sys.stdout.flush()
  #
  ls_result = run_refinement(
    structure_ideal=structure_ideal,
    structure_shake=structure_shake,
    params=params,
    run_id=run_id)
  if (ls_result is not None and params.pickle_root_name is not None):
    pickle_file_name = run_id + "_ls_history.pickle"
    print "writing file:", pickle_file_name
    easy_pickle.dump(
      file_name=pickle_file_name,
      obj=ls_result.history)
    print
    sys.stdout.flush()

def run(args):
  master_phil = libtbx.phil.parse("""
    general_positions_only = True
      .type = bool
    curvature_rescaling = True
      .type = bool
    minimizer = *lbfgs lbfgs_raw lbfgsb
      .type = choice
      .optional = False
    diagco = 0
      .type = int
    lbfgs_impl_switch = 1
      .type = int
    iprint = 1, 0
      .type = ints(size=2)
    curv_filter_lim_eps = 1e-3
      .type = float
    use_max_stp = True
      .type = bool
    shake_sites_rmsd = 0.5
      .type = float
    shake_adp_spread = 20
      .type = float
    pickle_root_name = None
      .type = str
    unpickle = None
      .type = path
""")
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  phil_objects = []
  remaining_args = []
  for arg in args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print
  params = work_phil.extract()
  if (params.unpickle is None):
    debug_utils.parse_options_loop_space_groups(
      argv=remaining_args, call_back=run_call_back, params=params)
  else:
    structure_ideal, structure_shake = easy_pickle.load(
      file_name=params.unpickle)
    run_refinement(
      structure_ideal=structure_ideal,
      structure_shake=structure_shake,
      params=params,
      run_id=params.unpickle)

if (__name__ == "__main__"):
  from cctbx.regression.tst_refine_xray_curvs import run # for pickle
  run(args=sys.argv[1:])
