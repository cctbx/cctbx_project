from __future__ import absolute_import, division, print_function
import os
import random
import time
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
import libtbx.load_env
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

random.seed(0)
flex.set_random_seed(0)

def get_fmodel(d_min=1.5, target_name="ls_wunit_k1"):
  """
  1yjp (P 1 21 1, 66 atoms) with noise-free simulated amplitudes |Fcalc|
  from direct summation, as the wrappers tests use.
  """
  pdb_file = libtbx.env.find_in_repositories(
    relative_path=os.path.join("iotbx", "regression", "data", "1yjp.pdb"),
    test=os.path.isfile)
  assert pdb_file is not None
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(model_input=pdb_inp)
  xrs = model.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=d_min, algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    xray_structure=xrs,
    target_name=target_name,
    sf_and_grads_accuracy_params=sfg_params)
  return model, fmodel

def shake(fmodel):
  """Move the model off the truth so that the gradients are not zero."""
  xrs = fmodel.xray_structure
  xrs.shake_sites_in_place(mean_distance=0.3)
  xrs.shake_adp(spread=10)
  fmodel.update_xray_structure(update_f_calc=True)

def exercise_joint_gradient_against_finite_differences():
  """
  The joint data term refines Cartesian sites and u_iso together, in block
  order: x = (sites_cart of all atoms | u_iso of all atoms). Its gradient
  must match central finite differences of its target on the selected
  entries and be zero on the unselected ones.
  """
  from mmtbx.refinement import data
  model, fmodel = get_fmodel()
  shake(fmodel)
  n = fmodel.xray_structure.scatterers().size()
  selection = flex.bool(n, True)
  for i in range(0, n, 7): selection[i] = False
  d = data.fs(fmodel=fmodel)
  d.set_refine_sites_and_u_iso(selection=selection)
  x = d.get_x()
  assert x.size() == 4 * n
  assert approx_equal(
    x[:3 * n], fmodel.xray_structure.sites_cart().as_double())
  assert approx_equal(
    x[3 * n:], fmodel.xray_structure.extract_u_iso_or_u_equiv())
  d.update(x=x)
  g = d.gradients()
  assert g.size() == 4 * n
  # unselected atoms get no gradient, in either block
  for i in range(n):
    if selection[i]: continue
    for k in range(3): assert g[3 * i + k] == 0
    assert g[3 * n + i] == 0
  # central differences on a few selected entries of each block
  checked = 0
  for i in range(n):
    if not selection[i] or i % 5 != 1: continue
    for idx, h in [(3 * i, 1.e-4), (3 * i + 1, 1.e-4), (3 * i + 2, 1.e-4),
                   (3 * n + i, 1.e-5)]:
      xp = x.deep_copy(); xp[idx] += h
      xm = x.deep_copy(); xm[idx] -= h
      d.update(x=xp); tp = d.target()
      d.update(x=xm); tm = d.target()
      fd = (tp - tm) / (2 * h)
      assert approx_equal(g[idx], fd, eps=1.e-4 * max(1., abs(fd))), (
        idx, g[idx], fd)
      checked += 1
  assert checked >= 40
  d.update(x=x)

def exercise_joint_gradient_equals_single_modes():
  """
  With all atoms selected the joint gradient is the Cartesian site gradient
  followed by the u_iso gradient of the two existing single-type modes.
  """
  from mmtbx.refinement import data
  model, fmodel = get_fmodel()
  shake(fmodel)
  n = fmodel.xray_structure.scatterers().size()
  d_joint = data.fs(fmodel=fmodel)
  d_joint.set_refine_sites_and_u_iso()
  x = d_joint.get_x()
  d_joint.update(x=x)
  g_joint = d_joint.gradients()
  t_joint = d_joint.target()
  #
  d_sites = data.fs(fmodel=fmodel)
  d_sites.set_refine_sites()
  d_sites.update(x=d_sites.get_x())
  g_sites = d_sites.gradients()
  assert approx_equal(d_sites.target(), t_joint)
  #
  d_u = data.fs(fmodel=fmodel)
  d_u.set_refine_u_iso()
  d_u.update(x=d_u.get_x())
  g_u = d_u.gradients()
  assert approx_equal(d_u.target(), t_joint)
  #
  assert g_sites.size() == 3 * n and g_u.size() == n
  assert approx_equal(g_joint[:3 * n], g_sites)
  assert approx_equal(g_joint[3 * n:], g_u)
  assert flex.max(flex.abs(g_sites)) > 0 and flex.max(flex.abs(g_u)) > 0

def exercise_joint_bounds():
  """
  The restraints manager provides box bounds for the joint parameter
  vector: sites within max_shift of their current values and u_iso
  between u(b_min) and u(b_max) for selected atoms, zero width otherwise.
  Restraint targets are not available in the joint mode.
  """
  from mmtbx.refinement import restraints
  from cctbx import adptbx
  model, fmodel = get_fmodel()
  n = model.size()
  selection = flex.bool(n, True)
  for i in range(0, n, 7): selection[i] = False
  r = restraints.manager(model=model)
  r.set_use_xyz_adp(selection=selection, max_shift=2.0, b_min=1, b_max=200)
  x = r.get_x()
  xrs = model.get_xray_structure()
  assert x.size() == 4 * n
  assert approx_equal(x[:3 * n], xrs.sites_cart().as_double())
  assert approx_equal(x[3 * n:], xrs.extract_u_iso_or_u_equiv())
  lo, up, fl = r.lower_bound(), r.upper_bound(), r.bound_flags()
  assert lo.size() == 4 * n and up.size() == 4 * n and fl.size() == 4 * n
  assert fl.all_eq(2)
  u_min, u_max = adptbx.b_as_u(1), adptbx.b_as_u(200)
  for i in range(n):
    for k in range(3):
      j = 3 * i + k
      if selection[i]:
        assert approx_equal(lo[j], x[j] - 2.0)
        assert approx_equal(up[j], x[j] + 2.0)
      else:
        assert lo[j] == x[j] and up[j] == x[j]
    j = 3 * n + i
    if selection[i]:
      assert approx_equal(lo[j], u_min) and approx_equal(up[j], u_max)
    else:
      assert lo[j] == x[j] and up[j] == x[j]
  assert r.target() is None and r.gradients() is None
  r.update(x=x)
  # restrained joint refinement is not implemented
  try:
    restraints.manager(model=model, use_target=True).set_use_xyz_adp()
  except AssertionError:
    pass
  else:
    raise RuntimeError("use_target=True must be rejected in joint mode")

def joint_calculator(model, fmodel, data_weight=2.0):
  from mmtbx.refinement import data, restraints, calculators
  from cctbx import adptbx
  d = data.fs(fmodel=fmodel)
  r = restraints.manager(model=model)
  return calculators.xyz_adp(
    data        = d,
    restraints  = r,
    data_weight = data_weight,
    max_shift   = 2.0,
    u_min       = adptbx.b_as_u(1),
    u_max       = adptbx.b_as_u(200)).calculator(), d, r

def within_bounds(calc, x):
  return (x >= calc.lower_bound).all_eq(True) and \
         (x <= calc.upper_bound).all_eq(True)

def exercise_xyz_adp_calculator():
  """
  calculators.xyz_adp builds an individual calculator on the joint data
  term and the joint bounds; L-BFGS-B runs on it and lowers the target
  while keeping the variables inside the box.
  """
  from scitbx import minimizers
  model, fmodel = get_fmodel()
  shake(fmodel)
  n = model.size()
  calc, d, r = joint_calculator(model, fmodel)
  assert calc.x.size() == 4 * n
  assert approx_equal(calc.x, d.get_x())
  assert approx_equal(calc.lower_bound, r.lower_bound())
  assert approx_equal(calc.upper_bound, r.upper_bound())
  t, g = calc.target_and_gradients()
  assert approx_equal(t, 2.0 * d.target())
  assert approx_equal(g, 2.0 * d.gradients())
  t_start = t
  minimizers.lbfgs(calculator=calc, max_iterations=5, mode="lbfgsb")
  assert calc.t < t_start, (calc.t, t_start)
  assert within_bounds(calc, calc.x)
  # the model carries the final variables
  assert approx_equal(
    calc.x[:3 * n], fmodel.xray_structure.sites_cart().as_double())
  assert approx_equal(
    calc.x[3 * n:], fmodel.xray_structure.extract_u_iso_or_u_equiv())

def exercise_scaled_calculator():
  """
  calculators.scaled lets the minimizer work on x' = scale * x: it exposes
  the scaled variables and bounds, returns the target of the plain
  calculator at x = x' / scale, and the gradient divided by the scale.
  With a unit scale it reproduces the plain calculator.
  """
  from mmtbx.refinement import calculators
  from scitbx import minimizers
  model, fmodel = get_fmodel()
  shake(fmodel)
  n = model.size()
  calc, d, r = joint_calculator(model, fmodel)
  x0 = calc.x.deep_copy()
  t0, g0 = calc.target_and_gradients()
  g0 = g0.deep_copy()
  # unit scale
  one = calculators.scaled(calculator=calc, scale=flex.double(4 * n, 1.0))
  assert approx_equal(one.x, x0)
  t, g = one.target_and_gradients()
  assert approx_equal(t, t0) and approx_equal(g, g0)
  # one scale per block
  scale = flex.double(3 * n, 0.5).concatenate(flex.double(n, 40.0))
  sc = calculators.scaled(calculator=calc, scale=scale)
  assert approx_equal(sc.x, x0 * scale)
  assert approx_equal(sc.lower_bound, calc.lower_bound * scale)
  assert approx_equal(sc.upper_bound, calc.upper_bound * scale)
  assert sc.bound_flags.all_eq(2)
  t, g = sc.target_and_gradients()
  assert approx_equal(t, t0) and approx_equal(g, g0 / scale)
  # a step in the scaled variables moves the plain variables by 1/scale
  step = flex.double(4 * n, 0.0)
  step[0] = 0.1
  step[3 * n] = 0.1
  xs = sc.x + step
  sc.update(x=xs)
  assert approx_equal(calc.x, xs / scale)
  assert approx_equal(calc.x[0], x0[0] + 0.2)
  assert approx_equal(calc.x[3 * n], x0[3 * n] + 0.0025)
  ts = sc.target()
  calc.update(x=xs / scale)
  assert approx_equal(ts, calc.target())
  # L-BFGS-B runs on the scaled calculator and respects the scaled box
  sc.update(x=x0 * scale)
  minimizers.lbfgs(calculator=sc, max_iterations=5, mode="lbfgsb")
  assert sc.target() < t0
  assert within_bounds(sc, sc.x)
  assert within_bounds(calc, calc.x)
  assert approx_equal(sc.x, calc.x * scale)
  # the scale must be positive
  try:
    calculators.scaled(calculator=calc, scale=flex.double(4 * n, 0.0))
  except AssertionError:
    pass
  else:
    raise RuntimeError("a non-positive scale must be rejected")

def exercise_gauss_newton_weights():
  """
  The Gauss-Newton weights are the work-set reflection weights times the
  expected radial curvature: for ls_wunit_k1, k^2 / sum(Fo^2) on acentric
  and twice that on centric reflections; for ml, (k alpha)^2 / (epsilon
  beta k^2) / n_work where alpha > 0 and beta > 1e-3. Zero on the test set.
  """
  from mmtbx.refinement import curvatures
  for target_name in ["ls_wunit_k1", "ml"]:
    model, fmodel = get_fmodel(target_name=target_name)
    shake(fmodel)
    w = curvatures.gauss_newton_weights(fmodel)
    flags = fmodel.r_free_flags().data()
    assert w.size() == flags.size()
    assert w.select(flags).all_eq(0)
    assert flex.min(w.select(~flags)) >= 0
    assert flex.max(w) > 0
    centric = fmodel.f_obs().centric_flags().data()
    if target_name == "ls_wunit_k1":
      fo = fmodel.f_obs_work().data()
      k = fmodel.scale_k1_w()
      base = k * k / flex.sum(fo * fo)
      expected = flex.double(w.size(), base)
      expected.set_selected(centric, 2 * base)
    else:
      alpha, beta = fmodel.alpha_beta()
      k = fmodel.scale_ml_wrapper()
      a = alpha.data() * k
      sigma = fmodel.epsilons * beta.data() * k * k
      ok = (alpha.data() > 0) & (beta.data() > 1.e-3)
      expected = flex.double(w.size(), 0)
      expected.set_selected(ok, (a * a / sigma).select(ok))
      expected *= 1.0 / (~flags).count(True)
    expected.set_selected(flags, 0)
    assert approx_equal(w, expected)

def exercise_gauss_newton_diagonal():
  """
  The Gauss-Newton diagonal for the Cartesian coordinates and u_iso equals
  sum_h weights_h |dF_h/dp|^2 with dF/dp from central differences of the
  structure factors, in the monoclinic cell of 1yjp.
  """
  from mmtbx.refinement import curvatures
  model, fmodel = get_fmodel()
  shake(fmodel)
  xrs = fmodel.xray_structure
  n = xrs.scatterers().size()
  indices = fmodel.f_model().indices()
  weights = curvatures.gauss_newton_weights(fmodel)
  sites, u_iso = curvatures.gauss_newton_diagonal_sites_u_iso(
    xray_structure=xrs, miller_indices=indices, weights=weights)
  assert sites.size() == 3 * n and u_iso.size() == n
  assert flex.min(sites) > 0 and flex.min(u_iso) > 0
  def f_calc(xs):
    return fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure=xs, algorithm="direct").f_calc().data()
  def gn(plus, minus, h):
    d_f = (f_calc(plus) - f_calc(minus)) * (1.0 / (2 * h))
    return flex.sum(weights * flex.norm(d_f))
  checked = 0
  for i in range(0, n, 7):
    h = 1.e-4
    for k in range(3):
      xp, xm = xrs.deep_copy_scatterers(), xrs.deep_copy_scatterers()
      for xs, sign in [(xp, 1), (xm, -1)]:
        sc = xs.sites_cart()
        site = list(sc[i])
        site[k] += sign * h
        sc[i] = tuple(site)
        xs.set_sites_cart(sc)
      fd = gn(xp, xm, h)
      assert approx_equal(sites[3 * i + k], fd, eps=1.e-5 * fd), (
        i, k, sites[3 * i + k], fd)
    xp, xm = xrs.deep_copy_scatterers(), xrs.deep_copy_scatterers()
    u = xrs.scatterers()[i].u_iso
    xp.scatterers()[i].u_iso = u + h
    xm.scatterers()[i].u_iso = u - h
    fd = gn(xp, xm, h)
    assert approx_equal(u_iso[i], fd, eps=1.e-5 * fd), (i, u_iso[i], fd)
    checked += 1
  assert checked >= 9

def run():
  exercise_joint_gradient_against_finite_differences()
  exercise_joint_gradient_equals_single_modes()
  exercise_joint_bounds()
  exercise_xyz_adp_calculator()
  exercise_scaled_calculator()
  exercise_gauss_newton_weights()
  exercise_gauss_newton_diagonal()
  print("OK")

if __name__ == "__main__":
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time() - t0))
