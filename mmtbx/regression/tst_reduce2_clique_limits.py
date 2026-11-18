"""reduce2: limits of the search of one multi-Mover Clique.

The C++ search cuts a Clique and searches each component again for every
state of the cut, so a large Clique can take hours and tens of GB. Before the
search, Optimizers._CliqueSearchWork() counts its work from the Clique's graph
(state combinations scored plus subsets tried to find cuts); above max_clique_search the Clique is
optimized as singletons. clique_search_seconds is a time backstop per Clique.
Each limited Clique gets a warning in the log and the description file.

Tests:
- test_work_count: _CliqueSearchWork on a path, star, cycle, complete graph and
  two joined cliques against a count that runs the recursion state by state,
  the limit on two Movers, and complete graphs of 18 Movers (counted) and 20
  (abandoned: too large to count).
- test_guard: 6vw1 Thr A 276 + Thr A 445 with max_clique_search=1: the Clique
  is optimized as singletons (the coarse score of the singletons path, as with
  the total time limit passed, not the search's), with a warning in the
  description; with two models, one warning per model.
- test_allowance: a time allowance only warns of a fallback the search shows:
  6vw1's two-Mover Clique (brute force, no time check) finishes with a tiny
  allowance and no warning; 7PRG's three-Mover Cliques fall back and warn,
  one line per alternate pass.
  Verbosity 0 gives the same warnings and models (the C++ fallback line is
  read at every verbosity), also for 1exr with His charges reconciled.
- test_defaults: 6vw1 with the default limits is the same as without limits.
- test_his_rerun: reconcile_his_charges on 1exr His A 107: a rerun whose
  limited Cliques differ from the free run's restores the free result
  (SEARCH_DEADLINE); the same limited Cliques in both runs reconcile; a real
  fallback in pass B of the free run and in pass A of the rerun differ.
- test_his_rerun_state: the same check with his_conformers=state (one state
  across conformers): a rerun that limits other Cliques restores the free
  result (A HD1+HE2, B HD1) instead of setting both to HD1.
"""
import itertools, math
from io import StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers
from mmtbx.regression.tst_reduce2_fresh_scores import pdb_6vw1
from mmtbx.regression.tst_reduce2_clique_polish import pdb_7prg, NEUTRON

def find_cut(nodes, adj):
  """findVertexCut(): the fewest Movers that split nodes, then the most even
  components (spread, then largest), first in combination order. Returns
  (cut, components, subsets tried) or (None, None, tried) without a cut."""
  N = len(nodes)
  tried = 0
  for k in range(N):
    best = None
    for removed in itertools.combinations(nodes, k):
      tried += 1
      keep = [v for v in nodes if v not in removed]
      if len(keep) < 2:
        continue
      comps, left = [], set(keep)
      while left:
        comp, todo = [], [min(left)]
        left.discard(todo[0])
        while todo:
          v = todo.pop()
          comp.append(v)
          for w in sorted(adj[v] & left):
            left.discard(w)
            todo.append(w)
        comps.append(sorted(comp))
      if len(comps) > 1:
        key = (max(map(len, comps)) - min(map(len, comps)), max(map(len, comps)))
        if best is None or key < best[0]:
          best = (key, list(removed), comps)
    if best is not None:
      return best[1], best[2], tried
  return None, None, tried

def stepwise_work(nodes, adj, ns):
  """The work units of _CliqueSearchWork, counted as the C++ recursion runs:
  each state of a cut is visited and its components searched again."""
  N = len(nodes)
  if N <= 2:
    return math.prod(ns[v] for v in nodes)
  cut, comps, tried = find_cut(nodes, adj)
  if cut is None:
    return math.prod(ns[v] for v in nodes) + tried
  w = tried
  for state in itertools.product(*[range(ns[v]) for v in cut]):
    w += 1
    for c in comps:
      w += stepwise_work(c, adj, ns)
  return w

def test_work_count():
  # two Movers: brute force, limited too
  assert Optimizers._CliqueSearchWork(2, [(0, 1)], [11, 11], 121) == 121
  assert Optimizers._CliqueSearchWork(2, [(0, 1)], [11, 11], 120) is None
  # path of 3, two states each: subsets 1 + 3, cut {1}, 2 x (1 + 2 + 2)
  assert Optimizers._CliqueSearchWork(3, [(0, 1), (1, 2)], [2, 2, 2], 1e9) == 14
  # triangle, three states each: no cut, 27 + 2**3 - 1 subsets
  assert Optimizers._CliqueSearchWork(3, [(0, 1), (1, 2), (0, 2)], [3, 3, 3], 1e9) == 34
  k4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]
  graphs = dict(
    path=(5, [(i, i + 1) for i in range(4)]),
    star=(5, [(0, i) for i in range(1, 5)]),
    cycle=(5, [(i, (i + 1) % 5) for i in range(5)]),
    complete=(4, k4),
    joined=(8, k4 + [(i + 4, j + 4) for i, j in k4] + [(3, 4)]))
  for name, (n, edges) in sorted(graphs.items()):
    ns = [2 + (i * 5) % 4 for i in range(n)]
    adj = dict((i, set()) for i in range(n))
    for i, j in edges:
      adj[i].add(j)
      adj[j].add(i)
    expected = stepwise_work(list(range(n)), adj, ns)
    got = Optimizers._CliqueSearchWork(n, edges, ns, 1e12)
    assert got == expected, (name, got, expected)
    # at the limit counted, above it None
    assert Optimizers._CliqueSearchWork(n, edges, ns, expected) == expected, name
    assert Optimizers._CliqueSearchWork(n, edges, ns, expected - 1) is None, name
  # complete graphs: 18 Movers are counted (2**18 states + 2**18 - 1 subsets), 20 are
  # not: their cut search alone visits more than _CLIQUE_WORK_MAX_STEPS Movers
  k = lambda n: [(i, j) for i in range(n) for j in range(i + 1, n)]
  assert Optimizers._CliqueSearchWork(18, k(18), [2] * 18, 1e9) == 2**19 - 1
  assert Optimizers._CliqueSearchWork(20, k(20), [2] * 20, 1e9) == Optimizers._CLIQUE_WORK_ABANDONED
  print("test_work_count OK")

def run(prefix, pdb_text, extra_args=()):
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[
    prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
    "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
    + list(extra_args))
  with open(prefix + ".txt") as f:
    desc = f.read()
  with open(prefix + "_out.pdb") as f:
    return desc, f.read()

def coarse_scores(desc):
  return [l.split("optimized with score ")[1].split()[0] for l in desc.split("\n")
    if l.startswith(" Clique of size ") and "optimized with score" in l]

def limited(desc):
  """The description's lines about Cliques whose search was limited."""
  return [l for l in desc.split("\n") if l.startswith(" Clique of ") and "(model " in l]

FLIPS = ["add_flip_movers=True"]

def test_guard():
  free, _ = run("tst_reduce2_clique_limits_guard", pdb_6vw1, FLIPS)
  assert limited(free) == [], limited(free)
  # the comparisons below need the Clique's line
  assert len(coarse_scores(free)) == 1, coarse_scores(free)
  desc, _ = run("tst_reduce2_clique_limits_guard", pdb_6vw1, FLIPS + ["max_clique_search=1"])
  lines = limited(desc)
  assert len(lines) == 1, lines
  assert lines[0].startswith(" Clique of 2 Movers (model 1, first SingleHydrogenRotator at chain A THR 276 HG1): "
    "optimized as singletons, search work count above the limit of 1; "), lines
  assert "Warning: the search of 1 Clique(s) was limited" in desc
  # the singletons path, as with the total time limit passed
  past, _ = run("tst_reduce2_clique_limits_guard", pdb_6vw1, FLIPS +
    ["treat_cliques_as_singletons_after_seconds=0"])
  assert coarse_scores(desc) == coarse_scores(past), (coarse_scores(desc), coarse_scores(past))
  assert coarse_scores(desc) != coarse_scores(free), (coarse_scores(desc), coarse_scores(free))
  assert "the total time limit had passed" in limited(past)[0], limited(past)
  # two models: the same Clique is limited in each, one line per model
  atoms = [l for l in pdb_6vw1.splitlines() if l.startswith(("ATOM", "HETATM"))]
  two = "\n".join(["MODEL        1"] + atoms + ["ENDMDL", "MODEL        2"] + atoms +
    ["ENDMDL", "END", ""])
  desc, _ = run("tst_reduce2_clique_limits_guard2", two, FLIPS + ["max_clique_search=1"])
  lines = limited(desc)
  assert [l.split(", first")[0] for l in lines] == [" Clique of 2 Movers (model 1",
    " Clique of 2 Movers (model 2"], lines
  assert "Warning: the search of 2 Clique(s) was limited" in desc
  print("test_guard OK")

def test_allowance():
  # brute force of two Movers has no time check: it finishes, no warning
  free, _ = run("tst_reduce2_clique_limits_time", pdb_6vw1, FLIPS)
  desc, _ = run("tst_reduce2_clique_limits_time", pdb_6vw1, FLIPS + ["clique_search_seconds=1e-9"])
  assert limited(desc) == [], limited(desc)
  assert coarse_scores(desc) == coarse_scores(free), (coarse_scores(desc), coarse_scores(free))
  # the same at verbosity 0
  quiet, _ = run("tst_reduce2_clique_limits_time", pdb_6vw1, FLIPS +
    ["clique_search_seconds=1e-9", "verbosity=0"])
  assert limited(quiet) == [], limited(quiet)
  # three-Mover Cliques are cut: out of time, they fall back
  desc, model = run("tst_reduce2_clique_limits_time", pdb_7prg, NEUTRON + ["clique_search_seconds=1e-9"])
  # one line per alternate pass: 2 Cliques in B, 3 in A
  lines = limited(desc)
  assert [l.split(", first")[0] for l in lines] == [" Clique of 3 Movers (model 1, alt B"] * 2 + \
    [" Clique of 3 Movers (model 1, alt A"] * 3, lines
  for l in lines:
    assert l.startswith(" Clique of 3 Movers ") and ("search ran out of its 1e-09 s allowance, "
      "part or all of it optimized as singletons") in l, l
  quiet, model0 = run("tst_reduce2_clique_limits_time", pdb_7prg, NEUTRON +
    ["clique_search_seconds=1e-9", "verbosity=0"])
  assert limited(quiet) == lines, (limited(quiet), lines)
  strip = lambda t: [l for l in t.split("\n") if not l.startswith("REMARK")]
  assert strip(model0) == strip(model)
  # His reruns: the outcome does not depend on verbosity
  from mmtbx.reduce import HisCharges
  from mmtbx.regression import tst_reduce2_his_charges_optimizer as his
  made = {}
  for v in (0, 2):
    atoms, _, m = his.run(his.pdb_1exr, "tst_reduce2_clique_limits_time_his%d" % v, his.ON,
      "verbosity=%d" % v, "clique_search_seconds=1e-9")
    made[v] = (atoms, his.only(m)['code'])
  assert made[0] == made[2] and made[0][1] == HisCharges.RECONCILED, (made[0][1], made[2][1])
  print("test_allowance OK")

def test_defaults():
  desc, model = run("tst_reduce2_clique_limits_default", pdb_6vw1, FLIPS)
  none, model_none = run("tst_reduce2_clique_limits_default", pdb_6vw1, FLIPS +
    ["max_clique_search=None", "clique_search_seconds=None"])
  assert limited(desc) == [] and limited(none) == []
  assert coarse_scores(desc) == coarse_scores(none)
  strip = lambda t: [l for l in t.split("\n") if not l.startswith("REMARK")]
  assert strip(model) == strip(model_none)
  print("test_defaults OK")

def test_his_rerun():
  from mmtbx.reduce import HisCharges
  from mmtbx.regression import tst_reduce2_his_charges_optimizer as his
  orig = Optimizers.Optimizer._OptimizeModel
  def run_with(prefix, runs):
    # a search-limit warning after the given runs (1: the free run, 2: rerun 1)
    calls = []
    def f(self, *a, **k):
      r = orig(self, *a, **k)
      calls.append(1)
      if len(calls) in runs:
        self._warningString += Optimizers._CLIQUE_WARNING + "Clique of 2 Movers (first test): test\n"
      return r
    Optimizers.Optimizer._OptimizeModel = f
    try:
      return his.run(his.pdb_1exr, prefix, his.ON)
    finally:
      Optimizers.Optimizer._OptimizeModel = orig
  off, _, _ = his.run(his.pdb_1exr, "tst_reduce2_clique_limits_his")
  on, _, made = run_with("tst_reduce2_clique_limits_his_r", (2,))
  m = his.only(made)
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.SEARCH_DEADLINE,
    'rerun 1 limited other Cliques than the first run', 1), m
  assert on == off
  on, _, made = run_with("tst_reduce2_clique_limits_his_s", (1, 2))
  m = his.only(made)
  assert m['code'] == HisCharges.RECONCILED and m['reruns'] == 1, m
  assert on != off
  # a real fallback (0 s) in pass B of the free run, in pass A of rerun 1
  origC = Optimizers.Optimizer._OptimizeCliques
  calls = []
  def model(self, *a, **k):
    calls.append(1)
    return orig(self, *a, **k)
  def cliques(self, *a):
    seconds = self._cliqueSeconds
    if (len(calls), self._cliqueAlt) in ((1, 'B'), (2, 'A')):
      self._cliqueSeconds = 0
    try:
      return origC(self, *a)
    finally:
      self._cliqueSeconds = seconds
  Optimizers.Optimizer._OptimizeModel = model
  Optimizers.Optimizer._OptimizeCliques = cliques
  try:
    _, _, made = his.run(his.pdb_1exr, "tst_reduce2_clique_limits_his_a", his.ON)
  finally:
    Optimizers.Optimizer._OptimizeModel = orig
    Optimizers.Optimizer._OptimizeCliques = origC
  m = his.only(made)
  assert (m['code'], m['detail']) == (HisCharges.SEARCH_DEADLINE,
    'rerun 1 limited other Cliques than the first run'), m
  print("test_his_rerun OK")

def test_his_rerun_state():
  from mmtbx.reduce import HisCharges
  from mmtbx.regression import tst_reduce2_his_charges_optimizer as his
  orig = Optimizers.Optimizer._OptimizeModel
  calls = []
  def f(self, *a, **k):
    # a search-limit warning after rerun 1 only
    r = orig(self, *a, **k)
    calls.append(1)
    if len(calls) == 2:
      self._warningString += Optimizers._CLIQUE_WARNING + "Clique of 2 Movers (first test): test\n"
    return r
  Optimizers.Optimizer._OptimizeModel = f
  try:
    atoms, _, made = his.run(his.pdb_1exr, "tst_reduce2_clique_limits_his_st", his.ON, his.STATE)
  finally:
    Optimizers.Optimizer._OptimizeModel = orig
  assert [o._reconcileHisStates for o in made] == [True], made
  m = his.only(made)
  c = his.case(m, '107')
  assert (m['code'], m['detail'], m['reruns']) == (HisCharges.SEARCH_DEADLINE,
    'rerun 1 limited other Cliques than the first run', 1), m
  assert (c['free'], c['final']) == ({'A': 'HD1+HE2', 'B': 'HD1'}, {'A': 'HD1+HE2', 'B': 'HD1'}), c
  assert his.ring_h(atoms, 'A', '107') == [('HD1', 'A'), ('HD1', 'B'), ('HE2', 'A')], \
    his.ring_h(atoms, 'A', '107')
  print("test_his_rerun_state OK")

if __name__ == "__main__":
  test_work_count()
  test_guard()
  test_allowance()
  test_defaults()
  test_his_rerun()
  test_his_rerun_state()
  print("OK")
