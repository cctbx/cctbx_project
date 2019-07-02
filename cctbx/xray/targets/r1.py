from __future__ import absolute_import, division, print_function
def calc_k(fo, fc):
  "scale factor for (fo-k*fc)**2, only similar to factor for abs(fo-k*fc)"
  from scitbx.array_family import flex
  k_num = flex.sum(fo * fc)
  k_den = flex.sum(fc * fc)
  assert k_den != 0
  assert k_den**2 != 0
  k = k_num / k_den
  k_d_num = fo * k_den - k_num * 2 * fc
  k_d = k_d_num / k_den**2
  k_d2 = -2 * (k_num * k_den + 2 * k_d_num * fc) / k_den**3
  return k, k_d, k_d2

def calc_t(fo, fc, k, k_d, k_d2):
  from scitbx.array_family import flex
  deltas = fo - k * fc
  signs = flex.double(fo.size(), 1)
  signs.set_selected(deltas < 0, -1)
  deltas *= signs
  t_num = flex.sum(deltas)
  t_den = flex.sum(fo)
  assert t_den != 0
  t = t_num / t_den
  if (k_d is None):
    return t
  t_d = -(k_d * flex.sum(signs * fc) + signs * k) / t_den
  if (k_d2 is None):
    return t, t_d
  t_d2 = -(k_d2 * flex.sum(signs * fc) + 2 * signs * k_d) / t_den
  return t, t_d, t_d2

class target(object):

  def __init__(O, f_obs, f_calc=None, f_calc_abs=None, fca_sq_eps=1e-100):
    assert [f_calc, f_calc_abs].count(None) == 1
    if (f_calc is not None):
      from scitbx.array_family import flex
      f_calc_abs = flex.abs(f_calc)
    O.k, O.k_d, O.k_d2 = calc_k(f_obs, f_calc_abs)
    O.t, O.g, O.c = calc_t(f_obs, f_calc_abs, O.k, O.k_d, O.k_d2)
    O.f_calc_gradients = None
    O.f_calc_hessians = None
    if (f_calc is not None):
      fca_sq = f_calc_abs**2
      isel_zero = (fca_sq <= fca_sq_eps).iselection()
      f_calc_abs.set_selected(isel_zero, 1)
      fca_sq.set_selected(isel_zero, 1)
      O.f_calc_gradients = O.g / f_calc_abs * f_calc
      O.f_calc_gradients.set_selected(isel_zero, 0j)
      a = flex.real(f_calc)
      b = flex.imag(f_calc)
      aa, bb, ab = a*a, b*b, a*b
      haa = O.c * aa + O.g * bb / f_calc_abs
      hbb = O.c * bb + O.g * aa / f_calc_abs
      hab = (O.c - O.g / f_calc_abs) * a * b
      O.f_calc_hessians = flex.vec3_double(haa, hbb, hab) / fca_sq
      O.f_calc_hessians.set_selected(isel_zero, (0,0,0))

  def target_work(O):
    return O.t

  def gradients_work(O):
    return O.f_calc_gradients

  def hessians_work(O):
    return O.f_calc_hessians
