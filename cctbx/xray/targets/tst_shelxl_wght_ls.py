from __future__ import division
from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

def calc_k(f_obs, i_calc):
  num = 0
  den = 0
  for fo,ic in zip(f_obs, i_calc):
    fc = ic**0.5
    num += fo*fc
    den += fc*fc
  assert den != 0
  k = num / den
  return k

def calc_w(wa, wb, i_obs, i_sig, i_calc, k):
  assert i_sig.size() == i_obs.size()
  assert i_calc.size() == i_obs.size()
  weights = flex.double()
  for io,so,ic in zip(i_obs, i_sig, i_calc):
    ik = io / k**2
    sk = so / k**2
    if (ik < 0): ik = 0
    p = (ik + 2 * ic) / 3
    den = sk**2 + (wa*p)**2 + wb*p
    assert den > 1e-8
    w = 1 / den
    weights.append(w)
  return weights

def calc_t(i_obs, i_calc, k, weights):
  t_num = 0
  t_den = 0
  for io,ic,w in zip(i_obs, i_calc, weights):
    delta = io - k**2 * ic
    t_num += w * delta**2
    t_den += w * io**2
  return t_num / t_den

def kwt(f_obs, i_obs, i_sig, i_calc, wa, wb):
  k = calc_k(f_obs, i_calc)
  weights = calc_w(
    wa=wa,
    wb=wb,
    i_obs=i_obs,
    i_sig=i_sig,
    i_calc=i_calc,
    k=k)
  t = calc_t(
    i_obs=i_obs,
    i_calc=i_calc,
    k=k,
    weights=weights)
  result = xray.targets_shelxl_wght_ls_kwt_b_dv(
    f_obs=f_obs,
    i_obs=i_obs,
    i_sig=i_sig,
    ic=i_calc,
    wa=wa,
    wb=wb)
  assert len(result) == 2
  icb, icbd = result
  return t, icb, icbd

def exercise(mt, n_refl):
  f_obs = mt.random_double(size=n_refl)
  i_obs = flex.pow2(f_obs)
  i_sig = mt.random_double(size=i_obs.size())
  i_calc = mt.random_double(size=f_obs.size())
  wa = 1.23
  wb = 2.34
  _, g_ana, c_ana = kwt(
    f_obs=f_obs, i_obs=i_obs, i_sig=i_sig, i_calc=i_calc, wa=wa, wb=wb)
  eps = 1e-6
  g_fin = flex.double()
  c_fin = flex.double()
  for ih in xrange(i_calc.size()):
    fs = []
    gs = []
    c_orig = i_calc[ih]
    for signed_eps in [eps, -eps]:
      i_calc[ih] = c_orig + signed_eps
      t_eps, g_ana_eps, _ = kwt(
        f_obs=f_obs, i_obs=i_obs, i_sig=i_sig, i_calc=i_calc, wa=wa, wb=wb)
      fs.append(t_eps)
      gs.append(g_ana_eps[ih])
    g_fin.append((fs[0]-fs[1])/(2*eps))
    c_fin.append((gs[0]-gs[1])/(2*eps))
    i_calc[ih] = c_orig
  print "g_fin:", numstr(g_fin)
  print "  ana:", numstr(g_ana)
  assert approx_equal(g_ana, g_fin)
  print "c_fin:", numstr(c_fin)
  print "  ana:", numstr(c_ana)
  assert approx_equal(c_ana, c_fin)
  print

def run(args):
  assert len(args) < 3
  arg_vals = [int(arg) for arg in args]
  arg_vals = arg_vals + [3, 2][len(arg_vals):]
  n_refl, n_trials = arg_vals
  assert n_refl > 0
  assert n_trials > 0
  mt = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    exercise(mt, n_refl)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
