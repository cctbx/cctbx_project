from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import mmtbx.refinement.minimization as minimization

""" The adaptive stopping rule of mmtbx.refinement.minimization.lbfgs past
the usual iteration limit, switched on by adaptive_params (phenix's
llgi_data.lbfgs_* settings): with lbfgs_stop = free, continue while the
free data target has reached a new minimum within the last lbfgs_patience
iterations, then restore the parameters at that minimum (the starting point
included); with lbfgs_stop = progress, stop when the total target no longer
falls. The rule is driven through callback_after_step with scripted
targets, and exercised in a real minimisation. """

class fake_lbfgs_minimizer(object):
  def __init__(self): self.i = 0
  def iter(self): return self.i

class fake_fmodel(object):
  def __init__(self, n_free, target_name):
    self.target_name = target_name
    self._flags = flex.bool(100, False)
    self._flags.set_selected(flex.size_t(range(n_free)), True)
  def r_free_flags(self): return group_args(data=lambda: self._flags)

class fake_fmodels(object):
  def __init__(self, fmodel): self.fmodel = fmodel
  def fmodel_xray(self): return self.fmodel
  def fmodel_neutron(self): return None

class fake_target_result(object):
  def __init__(self, free): self.free = free
  def target_test(self): return self.free

def run_rule(free_values, n_min=25, n_max=100, patience=10, stop="free",
             n_free=10, target_name="llgi", adaptive=True, free_start=0.,
             f_values=None):
  """ Feed free_values (and total targets f_values; default always
  decreasing) one per iteration, x = iteration number, the starting point
  (x = 0) having free target free_start; return (iterations run, x after
  restore, log). """
  settings = group_args(lbfgs_max_iterations=n_max, lbfgs_stop=stop,
    lbfgs_patience=patience, lbfgs_stop_fraction=0.01)
  o = minimization.lbfgs.__new__(minimization.lbfgs)
  o.fmodels = fake_fmodels(fake_fmodel(n_free, target_name))
  o.adaptive_params = settings if adaptive else None
  o.verbose = 0
  o.refine_xyz = True
  o.log = StringIO()
  params = o._adaptive_termination(group_args(max_iterations=n_min))
  if(o._adaptive is None):              # the usual limit only
    assert params.max_iterations == n_min
    return n_min, None, ""
  assert params.max_iterations == n_max
  o.x = flex.double([0])
  if(o._free is not None):
    o._note_free_target(fake_target_result(free_start))  # the start
  m = fake_lbfgs_minimizer()
  for i, v in enumerate(free_values):
    m.i = i + 1
    o.x[0] = m.i
    if(o._free is not None): o._note_free_target(fake_target_result(v))
    o.f = -m.i if f_values is None else f_values[i]
    if(o.callback_after_step(m) is True or m.i >= n_max): break
  o.minimizer = m
  o._restore_best_free(params)
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
  n, x, log = run_rule([-i for i in range(1, 101)])
  assert (n, x) == (100, 100), (n, x)
  assert "free-target minimum at iteration 100" in log, log
  assert "(restored)" not in log, log

def exercise_start_is_best():
  # the free target only rises: back to the starting point (iteration 0)
  n, x, log = run_rule([0.1 * (i + 1) for i in range(100)])
  assert (n, x) == (25, 0), (n, x)
  assert "free-target minimum at iteration 0 (restored)" in log, log

def exercise_patience():
  values = [-(i + 1) if i < 30 else -30 + 0.1 * (i - 29) for i in range(100)]
  n, x, log = run_rule(values, patience=5)
  assert (n, x) == (35, 30), (n, x)

def exercise_restore_after_exception():
  # scitbx.lbfgs leaves x at the lowest total target after an exception;
  # the best-free iterate is restored even when it is the last one
  o = minimization.lbfgs.__new__(minimization.lbfgs)
  o.refine_xyz, o.log = True, StringIO()
  o._free = group_args(best_x=flex.double([7]), best_iter=7)
  m = fake_lbfgs_minimizer(); m.i = 7
  o.minimizer, o.x = m, flex.double([99])
  o._restore_best_free(group_args(max_iterations=100))
  assert o.x[0] == 7
  assert "iteration 7" in o.log.getvalue()
  assert "(restored)" not in o.log.getvalue()

def exercise_fallbacks():
  # no free set, or lbfgs_stop = progress: the total-target rule, which with
  # a steadily decreasing target runs to the limit; nothing is restored
  for kw in [dict(n_free=0), dict(stop="progress")]:
    n, x, log = run_rule([1.0] * 100, **kw)
    assert (n, x) == (100, 100), (kw, n, x)
    assert "free-target" not in log, log
  # no free target from the functor: also the progress rule
  n, x, log = run_rule([None] * 100, free_start=None)
  assert (n, x) == (100, 100), (n, x)
  # progress: the total target is flat from iteration 30, so at 40 the last
  # 10 iterations have reduced it by less than 1% of the reduction so far
  n, x, log = run_rule([1.0] * 100, stop="progress",
    f_values=[-min(i + 1, 30) for i in range(100)])
  assert (n, x) == (40, 40), (n, x)

def exercise_opt_in():
  # without adaptive_params, the usual limit for every target, llgi too;
  # with them, the same rule for ml
  values = [-(i + 1) if i < 40 else -40 + 0.1 * (i - 39) for i in range(100)]
  for target_name in ("llgi", "ml"):
    n, x, log = run_rule(values, target_name=target_name, adaptive=False)
    assert (n, x) == (25, None), (n, x)
  n, x, log = run_rule(values, target_name="ml")
  assert (n, x) == (50, 40), (n, x)
  assert "free-target minimum at iteration 40 (restored)" in log, log

def exercise_real_minimisation():
  """ Coordinate minimisation (target ml, no restraints) of a shaken
  structure with the free rule: the final free target is the smallest of
  those of the starting point and the accepted iterates, and the
  minimisation runs past the floor of 3 iterations. """
  import scitbx.lbfgs, iotbx.pdb, mmtbx, mmtbx.model
  from mmtbx.regression.llgi_test_utils import build_fmodel
  fmodel = build_fmodel(30, 2.0, seed=2)
  xrs = fmodel.xray_structure.deep_copy_scatterers()
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=flex.split_lines(xrs.as_pdb_file())),
    crystal_symmetry=xrs.crystal_symmetry())
  shaken = xrs.deep_copy_scatterers()
  shaken.shake_sites_in_place(mean_distance=0.3, random_double=None)
  shaken.scatterers().flags_set_grads(state=False)
  shaken.scatterers().flags_set_grad_site(
    iselection=flex.size_t(range(shaken.scatterers().size())))
  fmodel.update_xray_structure(xray_structure=shaken, update_f_calc=True)
  model.set_xray_structure(shaken)
  # free targets as the minimiser sees them (its functor keeps alpha and
  # beta fixed): the starting point, the accepted iterates, and the final
  # evaluation after the restore
  start, accepted = [], []
  class recording_lbfgs(minimization.lbfgs):
    def _note_free_target(self, target_result):
      if(self._free.best_iter is None):
        start.append(target_result.target_test())
      minimization.lbfgs._note_free_target(self, target_result)
    def callback_after_step(self, minimizer):
      if(self._free is not None):
        accepted.append(self._free.target_result.target_test())
      return minimization.lbfgs.callback_after_step(self, minimizer)
  settings = group_args(lbfgs_max_iterations=30, lbfgs_stop="free",
    lbfgs_patience=5, lbfgs_stop_fraction=0.01)
  log = StringIO()
  r = recording_lbfgs(restraints_manager=None,
    fmodels=mmtbx.fmodels(fmodel_xray=fmodel), model=model, refine_xyz=True,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations=3),
    adaptive_params=settings, collect_monitor=False, log=log)
  assert r.minimizer.iter() > 3, r.minimizer.iter()
  assert len(start) == 1 and len(accepted) == r.minimizer.iter()
  free_end = r._free.target_result.target_test()
  assert approx_equal(free_end, min(start + accepted), eps=1e-9)
  assert free_end < start[0]
  assert "free-target minimum at iteration" in log.getvalue()

if(__name__ == "__main__"):
  exercise_interior_minimum()
  exercise_still_improving()
  exercise_start_is_best()
  exercise_patience()
  exercise_restore_after_exception()
  exercise_fallbacks()
  exercise_opt_in()
  exercise_real_minimisation()
  print("OK")
