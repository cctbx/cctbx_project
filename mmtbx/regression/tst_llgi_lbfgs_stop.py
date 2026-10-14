from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import group_args
from six.moves import cStringIO as StringIO
import mmtbx.refinement.minimization as minimization

""" The LLGI stopping rule of mmtbx.refinement.minimization.lbfgs past the
usual iteration limit (llgi_data.lbfgs_stop = free): continue while the
free data target has reached a new minimum within the last lbfgs_patience
iterations, then restore the parameters at that minimum; for other targets
only with adaptive_params. The minimiser is driven through
callback_after_step with a scripted free target. """

class fake_lbfgs_minimizer(object):
  def __init__(self): self.i = 0
  def iter(self): return self.i

class fake_fmodel(object):
  def __init__(self, n_free, llgi_data, target_name="llgi"):
    self.target_name = target_name
    self._flags = flex.bool(100, False)
    self._flags.set_selected(flex.size_t(range(n_free)), True)
    self._llgi_data = llgi_data
  def r_free_flags(self): return group_args(data=lambda: self._flags)
  def llgi_data(self): return self._llgi_data

class fake_fmodels(object):
  def __init__(self, fmodel): self.fmodel = fmodel
  def fmodel_xray(self): return self.fmodel
  def fmodel_neutron(self): return None

def run_rule(free_values, n_min=25, n_max=100, patience=10, stop="free",
             n_free=10, target_name="llgi", adaptive=True):
  """ Feed free_values (one per iteration; x = iteration number) to the
  rule; return (iterations run, x after restore, log). """
  llgi_data = group_args(lbfgs_max_iterations=n_max, lbfgs_stop=stop,
    lbfgs_patience=patience, lbfgs_stop_fraction=0.01)
  o = minimization.lbfgs.__new__(minimization.lbfgs)
  o.fmodels = fake_fmodels(fake_fmodel(n_free, llgi_data, target_name))
  # other targets: the settings come from adaptive_params, if given
  o.adaptive_params = llgi_data if (target_name != "llgi" and adaptive) else None
  o.verbose = 0
  o.refine_xyz = True
  o.log = StringIO()
  params = o._llgi_adaptive_termination(
    group_args(max_iterations=n_min))
  if(o._adaptive is None):              # the usual limit only
    assert params.max_iterations == n_min
    return n_min, None, ""
  assert params.max_iterations == n_max
  m = fake_lbfgs_minimizer()
  o.x = flex.double([0])
  for i, v in enumerate(free_values):
    m.i = i + 1
    o.x[0] = m.i
    if(o._free is not None): o._free.last = v
    o.f = -m.i              # total target always decreasing
    if(o.callback_after_step(m) is True or m.i >= n_max): break
  o.minimizer = m
  o._llgi_restore_best_free(params)
  return m.i, o.x[0], o.log.getvalue()

def exercise_interior_minimum():
  # minimum at 20, before the floor of 25: stop at iteration 30 (10 since
  # the best) and restore iteration 20
  values = [-(i + 1) if i < 20 else -20 + 0.1 * (i - 19) for i in range(100)]
  n, x, log = run_rule(values)
  assert n == 30, n
  assert x == 20, x
  assert "free-target minimum at iteration 20 (restored)" in log, log
  # minimum at 40: continue until 50
  values = [-(i + 1) if i < 40 else -40 + 0.1 * (i - 39) for i in range(100)]
  n, x, log = run_rule(values)
  assert (n, x) == (50, 40), (n, x)

def exercise_still_improving():
  # free target falls throughout: run to the limit, nothing to restore
  n, x, log = run_rule([-i for i in range(100)])
  assert (n, x) == (100, 100), (n, x)
  assert "free-target minimum at iteration 100" in log, log
  assert "(restored)" not in log, log

def exercise_patience():
  values = [-(i + 1) if i < 30 else -30 + 0.1 * (i - 29) for i in range(100)]
  n, x, log = run_rule(values, patience=5)
  assert (n, x) == (35, 30), (n, x)

def exercise_fallbacks():
  # no free set, or lbfgs_stop = progress: the total-target rule, which with
  # a steadily decreasing target runs to the limit; nothing is restored
  for kw in [dict(n_free=0), dict(stop="progress")]:
    n, x, log = run_rule([1.0] * 100, **kw)
    assert (n, x) == (100, 100), (kw, n, x)
    assert "free-target" not in log, log

def exercise_other_targets():
  # ml with adaptive_params: the same rule; without: the usual limit
  values = [-(i + 1) if i < 40 else -40 + 0.1 * (i - 39) for i in range(100)]
  n, x, log = run_rule(values, target_name="ml")
  assert (n, x) == (50, 40), (n, x)
  assert "free-target minimum at iteration 40 (restored)" in log, log
  n, x, log = run_rule(values, target_name="ml", adaptive=False)
  assert (n, x) == (25, None), (n, x)

if(__name__ == "__main__"):
  exercise_interior_minimum()
  exercise_still_improving()
  exercise_patience()
  exercise_fallbacks()
  exercise_other_targets()
  print("OK")
