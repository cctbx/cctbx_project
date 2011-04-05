from cctbx.xray.targets import r1
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise(mt, n_refl, log):
  f_obs = mt.random_double(size=n_refl)
  f_calc = flex.complex_double(
    mt.random_double(size=f_obs.size()),
    mt.random_double(size=f_obs.size()))
  f_calc_abs = flex.abs(f_calc)
  tg = r1.target(f_obs, f_calc_abs)
  def check_f_calc_abs_derivs():
    eps = 1e-6
    g_fin = flex.double()
    c_fin = flex.double()
    for ih in xrange(f_calc_abs.size()):
      fs = []
      gs = []
      c_orig = f_calc_abs[ih]
      for signed_eps in [eps, -eps]:
        f_calc_abs[ih] = c_orig + signed_eps
        tg_eps = r1.target(f_obs, f_calc_abs)
        fs.append(tg_eps.t)
        gs.append(tg_eps.t_d[ih])
      g_fin.append((fs[0]-fs[1])/(2*eps))
      c_fin.append((gs[0]-gs[1])/(2*eps))
      f_calc_abs[ih] = c_orig
    print >> log, "g fin:", numstr(g_fin)
    print >> log, "  ana:", numstr(tg.t_d)
    assert approx_equal(tg.t_d, g_fin)
    print >> log, "c fin:", numstr(c_fin)
    print >> log, "  ana:", numstr(tg.t_d2)
    assert approx_equal(tg.t_d2, c_fin)
  check_f_calc_abs_derivs()

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
  for i_trial in xrange(n_trials):
    exercise(mt, n_refl, log)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
