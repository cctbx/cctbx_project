from __future__ import absolute_import, division, print_function
import time
from cctbx import maptbx
from cctbx import crystal
from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx import matrix

"""
Real-space target (sum of map values at atom centers) and its gradients as
computed by maptbx.target_and_gradients_simple with each method available in
real-space refinement: fd (finite differences of linearly interpolated values),
linear, tricubic (analytical gradients of the interpolant). Quadratic
interpolation is not offered: its target is discontinuous.
"""

METHODS = ["fd", "linear", "tricubic"]

def get_inputs(unit_cell, d_min=2.0):
  cs = crystal.symmetry(unit_cell, "P1")
  sites_frac = flex.vec3_double([
    (0.31, 0.42, 0.53), (0.37, 0.45, 0.50), (0.40, 0.38, 0.58),
    (0.52, 0.61, 0.33), (0.66, 0.25, 0.71)])
  scatterers = flex.xray_scatterer([
    xray.scatterer(label="C%d"%i, site=s, u=0.2, scattering_type="C")
    for i, s in enumerate(sites_frac)])
  xrs = xray.structure(crystal_symmetry=cs, scatterers=scatterers)
  fft_map = xrs.structure_factors(d_min=d_min).f_calc().fft_map(
    resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  return xrs, fft_map.real_map_unpadded()

def tg(unit_cell, map_data, sites_cart, method, selection=None, delta=0.5):
  if(selection is None): selection = flex.bool(sites_cart.size(), True)
  if(method == "fd"):
    return maptbx.target_and_gradients_simple(
      unit_cell   = unit_cell,
      map_target  = map_data,
      sites_cart  = sites_cart,
      delta       = delta,
      selection   = selection)
  return maptbx.target_and_gradients_simple(
    unit_cell     = unit_cell,
    map_target    = map_data,
    sites_cart    = sites_cart,
    selection     = selection,
    interpolation = method)

def shifted(sites_cart, i_site, i_axis, h):
  result = sites_cart.deep_copy()
  s = list(result[i_site])
  s[i_axis] += h
  result[i_site] = s
  return result

def exercise_gradients(unit_cell, methods=METHODS):
  """
  Gradients agree with finite differences of the target of the same method
  (Cartesian shifts). Atoms are displaced from map peaks so that gradients are
  not small.
  """
  xrs, map_data = get_inputs(unit_cell = unit_cell)
  uc = xrs.unit_cell()
  sites_cart = xrs.sites_cart() + flex.vec3_double(
    xrs.scatterers().size(), (0.23, -0.17, 0.31))
  for method in methods:
    # fd gradients are central differences with step delta of the target
    h   = 0.5   if method == "fd" else 1.e-4
    eps = 1.e-6 if method == "fd" else 1.e-3
    g = tg(uc, map_data, sites_cart, method, delta=h).gradients()
    g_max = flex.max(flex.abs(g.as_double()))
    assert g_max > 0.1, [method, g_max]
    for i_site in range(sites_cart.size()):
      for i_axis in [0,1,2]:
        tp = tg(uc, map_data, shifted(sites_cart, i_site, i_axis,  h),
          method, delta=h).target()
        tm = tg(uc, map_data, shifted(sites_cart, i_site, i_axis, -h),
          method, delta=h).target()
        g_fd = (tp-tm)/(2*h)
        assert approx_equal(g[i_site][i_axis], g_fd, eps*max(1., g_max)), [
          unit_cell, method, i_site, i_axis, g[i_site][i_axis], g_fd]

def exercise_selection(unit_cell):
  """
  Atoms that are not selected do not contribute to target and gradients.
  """
  xrs, map_data = get_inputs(unit_cell = unit_cell)
  uc = xrs.unit_cell()
  sites_cart = xrs.sites_cart() + flex.vec3_double(
    xrs.scatterers().size(), (0.23, -0.17, 0.31))
  sel = flex.bool([True, False, True, True, False])
  for method in METHODS:
    o_all = tg(uc, map_data, sites_cart, method)
    o_sel = tg(uc, map_data, sites_cart, method, selection = sel)
    o_ref = tg(uc, map_data, sites_cart.select(sel), method)
    assert approx_equal(o_sel.target(), o_ref.target(), 1.e-9)
    assert abs(o_all.target()-o_sel.target()) > 1.e-3
    g_all, g_sel = o_all.gradients(), o_sel.gradients()
    for i, s in enumerate(sel):
      if(s): assert approx_equal(g_sel[i], g_all[i], 1.e-9)
      else:  assert approx_equal(g_sel[i], (0,0,0), 1.e-12)

def exercise_continuity(unit_cell):
  """
  Target and gradient along a line that crosses several grid planes, sampled
  with step h and h/10. For a continuous function the largest change between
  neighbouring points scales with the step (ratio ~10); if there is a jump it
  does not (ratio ~1). Expected:
    fd:        target and gradient continuous (but the gradient is a finite
               difference, not the exact derivative of the target);
    linear:    target continuous, gradient jumps at grid planes (by up to
               ~40% of its magnitude here);
    tricubic:  target continuous, gradient jumps at grid planes too, several
               times less than with linear.
  """
  xrs, map_data = get_inputs(unit_cell = unit_cell)
  uc = xrs.unit_cell()
  start = xrs.sites_cart()[0]
  direction = (0.61, 0.52, 0.60)
  length = 2.0 # grid step is ~0.5 A
  def max_changes(method, h):
    t, g = flex.double(), flex.double()
    for i in range(int(length/h)):
      site = [start[k] + direction[k]*h*i for k in [0,1,2]]
      o = tg(uc, map_data, flex.vec3_double([site]), method)
      t.append(o.target())
      gi = o.gradients()[0]
      g.append(gi[0]*direction[0]+gi[1]*direction[1]+gi[2]*direction[2])
    return flex.max(flex.abs(t[1:]-t[:-1])), flex.max(flex.abs(g[1:]-g[:-1]))
  expected = { # method: (target continuous, gradient continuous)
    "fd":        (True,  True),
    "linear":    (True,  False),
    "tricubic":  (True,  False)}
  gradient_jump = {}
  for method in METHODS:
    dt1, dg1 = max_changes(method, 1.e-3)
    dt2, dg2 = max_changes(method, 1.e-4)
    gradient_jump[method] = dg2
    for ratio, continuous, what in [
          (dt1/dt2, expected[method][0], "target"),
          (dg1/dg2, expected[method][1], "gradient")]:
      if(continuous): assert ratio > 5, [unit_cell, method, what, ratio]
      else:           assert ratio < 2, [unit_cell, method, what, ratio]
  assert gradient_jump["tricubic"] < 0.2*gradient_jump["linear"], gradient_jump

def exercise_quadratic_is_not_a_choice(unit_cell):
  """
  Quadratic interpolation has a discontinuous target and was removed from the
  methods available for refinement.
  """
  xrs, map_data = get_inputs(unit_cell = unit_cell)
  try:
    tg(xrs.unit_cell(), map_data, xrs.sites_cart(), "quadratic")
  except RuntimeError as e:
    assert str(e).count("Unknown interpolation mode")==1, str(e)
  else: raise RuntimeError("Exception expected")

def exercise_magnification(unit_cell):
  """
  target_and_gradients_simple_magnification: target is the sum of map values
  at K*sites_cart; gradients with respect to elements of K (3x3 matrix or its
  diagonal) agree with finite differences in any cell.
  """
  xrs, map_data = get_inputs(unit_cell = unit_cell)
  uc = xrs.unit_cell()
  sites_cart = xrs.sites_cart() + flex.vec3_double(
    xrs.scatterers().size(), (0.23, -0.17, 0.31))
  def run(k):
    return maptbx.target_and_gradients_simple_magnification(
      unit_cell = uc, map_target = map_data, sites_cart = sites_cart, K = k)
  h = 1.e-5
  for k0, make in [
        ([1.01,0.01,-0.02, 0.015,0.99,0.01, -0.01,0.02,1.02], matrix.sqr),
        ([1.01,0.99,1.02],                                    matrix.col)]:
    o = run(make(k0))
    g = o.gradients()
    assert g.size() == len(k0)
    assert approx_equal(o.target(), tg(uc, map_data,
      flex.vec3_double([make(k0)*matrix.col(s) if len(k0)==9 else
        [k0[i]*s[i] for i in [0,1,2]] for s in sites_cart]),
      "tricubic").target(), 1.e-6)
    g_max = flex.max(flex.abs(g))
    assert g_max > 1
    for i in range(len(k0)):
      kp, km = k0[:], k0[:]
      kp[i] += h
      km[i] -= h
      g_fd = (run(make(kp)).target()-run(make(km)).target())/(2*h)
      assert approx_equal(g[i], g_fd, 1.e-3*g_max), [unit_cell, i, g[i], g_fd]

if(__name__ == "__main__"):
  t0 = time.time()
  for unit_cell in [(20, 22, 25, 90, 90, 90), (20, 22, 25, 80, 105, 70)]:
    exercise_gradients(unit_cell = unit_cell)
    exercise_selection(unit_cell = unit_cell)
    exercise_continuity(unit_cell = unit_cell)
    exercise_quadratic_is_not_a_choice(unit_cell = unit_cell)
    exercise_magnification(unit_cell = unit_cell)
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
