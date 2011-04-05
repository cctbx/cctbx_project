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

  def __init__(O, fo, fc):
    O.k, O.k_d, O.k_d2 = calc_k(fo, fc)
    O.t, O.t_d, O.t_d2 = calc_t(fo, fc, O.k, O.k_d, O.k_d2)
