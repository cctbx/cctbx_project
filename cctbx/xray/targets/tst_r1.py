from __future__ import absolute_import, division, print_function
from cctbx.xray.targets import r1
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import range
from six.moves import zip

def exercise(mt, n_refl, log):
  f_obs = mt.random_double(size=n_refl)
  f_calc = flex.complex_double(
    mt.random_double(size=f_obs.size()),
    mt.random_double(size=f_obs.size()))
  f_calc_abs = flex.abs(f_calc)
  trg = r1.target(f_obs=f_obs, f_calc=f_calc)
  def check_f_calc_abs_derivs():
    eps = 1e-6
    g_fin = flex.double()
    c_fin = flex.double()
    for ih in range(f_calc_abs.size()):
      fs = []
      gs = []
      c_orig = f_calc_abs[ih]
      for signed_eps in [eps, -eps]:
        f_calc_abs[ih] = c_orig + signed_eps
        trg_eps = r1.target(f_obs=f_obs, f_calc_abs=f_calc_abs)
        fs.append(trg_eps.t)
        gs.append(trg_eps.g[ih])
      g_fin.append((fs[0]-fs[1])/(2*eps))
      c_fin.append((gs[0]-gs[1])/(2*eps))
      f_calc_abs[ih] = c_orig
    print("g fin:", numstr(g_fin), file=log)
    print("  ana:", numstr(trg.g), file=log)
    assert approx_equal(trg.g, g_fin)
    print("c fin:", numstr(c_fin), file=log)
    print("  ana:", numstr(trg.c), file=log)
    assert approx_equal(trg.c, c_fin)
  def check_f_calc_derivs():
    eps = 1e-6
    g_fin = flex.complex_double()
    c_fin = flex.vec3_double()
    for ih in range(f_calc.size()):
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
          trg_eps = r1.target(f_obs=f_obs, f_calc=f_calc)
          fs.append(trg_eps.t)
          gs.append(trg_eps.f_calc_gradients[ih])
        g_fin_ab.append((fs[0]-fs[1])/(2*eps))
        c_fin_ab.append((gs[0]-gs[1])/(2*eps))
      g_fin.append(complex(*g_fin_ab))
      assert approx_equal(c_fin_ab[0].imag, c_fin_ab[1].real)
      c_fin.append((c_fin_ab[0].real, c_fin_ab[1].imag, c_fin_ab[0].imag))
      f_calc[ih] = c_orig
    for pn,f,a in zip(
          g_fin.part_names(), g_fin.parts(), trg.f_calc_gradients.parts()):
      print("g fin %s:" % pn, numstr(f), file=log)
      print("  ana %s:" % pn, numstr(a), file=log)
    assert approx_equal(trg.f_calc_gradients, g_fin)
    for pn,f,a in zip(
          ["aa", "bb", "ab"], c_fin.parts(), trg.f_calc_hessians.parts()):
      print("c fin %s:" % pn, numstr(f), file=log)
      print("  ana %s:" % pn, numstr(a), file=log)
    assert approx_equal(trg.f_calc_hessians, c_fin)
  check_f_calc_abs_derivs()
  check_f_calc_derivs()
  #
  f_calc[0] = 0j
  f_calc_abs = flex.abs(f_calc)
  trg = r1.target(f_obs=f_obs, f_calc=f_calc)
  assert trg.gradients_work()[0] == 0j
  assert trg.hessians_work()[0] == (0,0,0)

def run(args):
  assert len(args) < 3
  arg_vals = [int(arg) for arg in args]
  arg_vals = arg_vals + [3, 2][len(arg_vals):]
  n_refl, n_trials = arg_vals
  assert n_refl > 0
  assert n_trials > 0
  if (len(args) == 0):
    from libtbx.utils import null_out
    log = null_out()
  else:
    import sys
    log = sys.stdout
  mt = flex.mersenne_twister(seed=0)
  for i_trial in range(n_trials):
    exercise(mt, n_refl, log)
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
