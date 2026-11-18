"""reduce2: limits of the search of one multi-Mover Clique.

The C++ search cuts a Clique and searches each component again for every
state of the cut, so a large Clique can take hours and tens of GB. Before the
search, Optimizers._CliqueSearchWork() counts its work from the Clique's graph
(state combinations scored plus subsets tried to find cuts); above maxCliqueWork the Clique is
optimized as singletons. cliqueSeconds is a time backstop per Clique.
Each limited Clique gets a warning.

Tests:
- test_work_count: _CliqueSearchWork on a path, star, cycle, complete graph and
  two joined cliques against a count that runs the recursion state by state,
  the limit on two Movers, and complete graphs of 18 Movers (counted) and 20
  (abandoned: too large to count).
- test_his_rerun: reconcile_his_charges on 1exr His A 107: a rerun whose
  limited Cliques differ from the free run's restores the free result
  (SEARCH_DEADLINE); the same limited Cliques in both runs reconcile; a real
  fallback in pass B of the free run and in pass A of the rerun differ.
"""
import itertools, math
from mmtbx.reduce import Optimizers

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

if __name__ == "__main__":
  test_work_count()
  test_his_rerun()
  print("OK")
