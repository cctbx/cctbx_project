from __future__ import division
import cctbx.xray.targets
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

def kwt(f_obs, i_obs, i_sig, f_calc, i_calc, wa, wb):
  if (f_calc is not None):
    i_calc = flex.norm(f_calc)
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
  if (f_calc is not None):
    i_calc = None
  trg = cctbx.xray.targets.shelxl_wght_ls(
    f_obs=f_obs,
    i_obs=i_obs,
    i_sig=i_sig,
    f_calc=f_calc,
    i_calc=i_calc,
    wa=wa,
    wb=wb)
  assert approx_equal(trg.scale_factor, k)
  assert approx_equal(trg.weights, weights)
  assert approx_equal(trg.target, t)
  return trg

def exercise(mt, n_refl):
  f_obs = mt.random_double(size=n_refl)
  i_obs = flex.pow2(f_obs)
  i_sig = mt.random_double(size=i_obs.size())
  f_calc = flex.complex_double(
    mt.random_double(size=f_obs.size()),
    mt.random_double(size=f_obs.size()))
  i_calc = flex.norm(f_calc)
  wa = 1.23
  wb = 2.34
  trg = kwt(
    f_obs=f_obs, i_obs=i_obs, i_sig=i_sig,
    f_calc=f_calc, i_calc=None, wa=wa, wb=wb)
  def check_i_derivs():
    g_ana = trg.i_gradients
    c_ana = trg.i_curvatures
    eps = 1e-6
    g_fin = flex.double()
    c_fin = flex.double()
    for ih in xrange(i_calc.size()):
      fs = []
      gs = []
      c_orig = i_calc[ih]
      for signed_eps in [eps, -eps]:
        i_calc[ih] = c_orig + signed_eps
        trg_eps = kwt(
          f_obs=f_obs, i_obs=i_obs, i_sig=i_sig,
          f_calc=None, i_calc=i_calc, wa=wa, wb=wb)
        fs.append(trg_eps.target)
        gs.append(trg_eps.i_gradients[ih])
      g_fin.append((fs[0]-fs[1])/(2*eps))
      c_fin.append((gs[0]-gs[1])/(2*eps))
      i_calc[ih] = c_orig
    assert approx_equal(g_ana, g_fin)
    assert approx_equal(c_ana, c_fin)
  def check_f_derivs():
    g_ana = trg.f_gradients
    c_ana = trg.f_hessians
    eps = 1e-6
    g_fin = flex.complex_double()
    c_fin = flex.vec3_double()
    for ih in xrange(i_calc.size()):
      c_orig = f_calc[ih]
      g_fin_ab = []
      c_fin_ab = []
      for iab in [0,1]:
        fs = []
        gs = []
        for signed_eps in [eps, -eps]:
          if (iab == 0):
            f_calc[ih] = complex(c_orig.real + signed_eps, c_orig.imag)
          else:
            f_calc[ih] = complex(c_orig.real, c_orig.imag + signed_eps)
          trg_eps = kwt(
            f_obs=f_obs, i_obs=i_obs, i_sig=i_sig,
            f_calc=f_calc, i_calc=None, wa=wa, wb=wb)
          fs.append(trg_eps.target)
          gs.append(trg_eps.f_gradients[ih])
        g_fin_ab.append((fs[0]-fs[1])/(2*eps))
        c_fin_ab.append((gs[0]-gs[1])/(2*eps))
      g_fin.append(complex(*g_fin_ab))
      assert approx_equal(c_fin_ab[0].imag, c_fin_ab[1].real)
      c_fin.append((c_fin_ab[0].real, c_fin_ab[1].imag, c_fin_ab[0].imag))
      f_calc[ih] = c_orig
    assert approx_equal(g_ana, g_fin)
    assert approx_equal(c_ana, c_fin)
  check_i_derivs()
  check_f_derivs()

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
