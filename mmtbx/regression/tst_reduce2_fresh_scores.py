"""reduce2: the description's final scores are fresh scores at the final poses.

The search stores each Mover's score when it sets the Mover's state; a Clique
partner that moves afterwards (in its fine step) leaves the stored score stale.
The description prints each Mover's score at its final pose against everything
around it at the end, and each Clique's total is the sum of these.

Test:
- 4wsc Lys C 245 NZ with Gln C 319 (and Phe C 219): the Gln's stored score was
  set before the Lys's fine step turned the NH3; every printed score equals a
  fresh rescoring, and the Clique total their sum.
"""
from __future__ import absolute_import, division, print_function
from io import StringIO
from iotbx.cli_parser import run_program
from libtbx.test_utils import approx_equal
from mmtbx.programs import reduce2
from mmtbx.probe import Helpers
from mmtbx.reduce import Optimizers

# 4wsc Phe C 219, Lys C 245, Gln C 319
lys_gln_str = """\
CRYST1  135.640  259.390  143.310  90.00 105.12  90.00 P 1 21 1
ATOM   9281  N   PHE C 219      -5.510 -12.352 -31.231  1.00 47.10           N
ATOM   9282  CA  PHE C 219      -4.604 -13.105 -32.094  1.00 49.36           C
ATOM   9283  C   PHE C 219      -3.693 -13.996 -31.243  1.00 37.63           C
ATOM   9284  O   PHE C 219      -3.391 -13.690 -30.090  1.00 31.68           O
ATOM   9285  CB  PHE C 219      -3.782 -12.165 -32.969  1.00 50.42           C
ATOM   9286  CG  PHE C 219      -4.555 -11.593 -34.106  1.00 53.72           C
ATOM   9287  CD1 PHE C 219      -4.933 -12.388 -35.169  1.00 55.35           C
ATOM   9288  CD2 PHE C 219      -4.920 -10.256 -34.105  1.00 64.45           C
ATOM   9289  CE1 PHE C 219      -5.657 -11.861 -36.212  1.00 60.53           C
ATOM   9290  CE2 PHE C 219      -5.641  -9.719 -35.148  1.00 68.56           C
ATOM   9291  CZ  PHE C 219      -6.014 -10.524 -36.203  1.00 66.84           C
ATOM   9479  N   LYS C 245      -6.780  -5.665 -31.988  1.00 62.07           N
ATOM   9480  CA  LYS C 245      -7.180  -6.883 -31.275  1.00 52.07           C
ATOM   9481  C   LYS C 245      -6.128  -7.246 -30.233  1.00 47.04           C
ATOM   9482  O   LYS C 245      -4.956  -6.933 -30.426  1.00 47.22           O
ATOM   9483  CB  LYS C 245      -7.366  -8.039 -32.257  1.00 50.07           C
ATOM   9484  CG  LYS C 245      -8.591  -7.921 -33.151  1.00 58.18           C
ATOM   9485  CD  LYS C 245      -8.826  -9.234 -33.875  1.00 62.94           C
ATOM   9486  CE  LYS C 245     -10.159  -9.269 -34.601  1.00 64.88           C
ATOM   9487  NZ  LYS C 245     -10.397 -10.602 -35.241  1.00 63.79           N
ATOM  10021  N   GLN C 319      -6.118 -17.653 -34.053  1.00 55.81           N
ATOM  10022  CA  GLN C 319      -7.173 -16.807 -33.532  1.00 52.62           C
ATOM  10023  C   GLN C 319      -8.283 -17.609 -32.837  1.00 56.71           C
ATOM  10024  O   GLN C 319      -8.580 -18.742 -33.220  1.00 63.19           O
ATOM  10025  CB  GLN C 319      -7.714 -15.979 -34.690  1.00 51.64           C
ATOM  10026  CG  GLN C 319      -8.318 -14.668 -34.292  1.00 67.66           C
ATOM  10027  CD  GLN C 319      -8.460 -13.739 -35.479  1.00 65.76           C
ATOM  10028  OE1 GLN C 319      -9.202 -12.747 -35.429  1.00 62.32           O
ATOM  10029  NE2 GLN C 319      -7.750 -14.061 -36.564  1.00 54.15           N
END
"""

# run() swaps Optimizers.Optimizer for Recorder while reduce2 runs
_Optimizer = Optimizers.Optimizer

class Recorder(_Optimizer):
  """After the fine step (the end of the search without symmetry copies):
  each Mover's stored score and a fresh score at its pose, with its atoms
  checked at the pose its OptimizerC holds."""
  states = {}
  def _OptimizeFine(self, optC, cliques):
    _Optimizer._OptimizeFine(self, optC, cliques)
    dot_scorer = Helpers.createDotScorer(self._extraAtomInfo, self._probePhil)
    # _ScoreMover() looks up each Mover's OptimizerC; the base sets it after this
    self._optimizerOf = dict((x, optC) for x in self._movers)
    for m in self._movers:
      coarse, fine = optC.GetCoarseLocation(m), optC.GetFineLocation(m)
      pr, index = m.CoarsePositions(), coarse
      if fine >= 0:
        pr, index = m.FinePositions(coarse), fine
      for a, xyz in zip(pr.atoms, pr.positions[index]):
        assert approx_equal(a.xyz, xyz, eps=1.e-6, out=None), (
          self._MoverName(m), a.name)
      Recorder.states[self._MoverName(m)] = (optC.GetHighScore(m),
        self._ScoreMover(m, dot_scorer))

def run(prefix, pdb_text, extra_args=()):
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  Recorder.states = {}
  Optimizers.Optimizer = Recorder
  try:
    run_program(program_class=reduce2.Program, logger=StringIO(), args=[
      prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
      "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
      + list(extra_args))
  finally:
    Optimizers.Optimizer = _Optimizer
  with open(prefix + ".txt") as f:
    return f.read()

def assert_fresh(desc, fresh):
  """Every printed final score is the fresh score, and every Clique's final
  total the sum of its Movers' fresh scores, to 2 decimals."""
  printed, groups = {}, []
  for l in desc.split("\n"):
    if "Totals:" in l:
      groups.append((l.split("final score ")[1].strip(), []))
    elif "Singleton Movers:" in l:
      groups.append((None, []))
    elif " final score: " in l:
      name = l.split(" Initial score:")[0].strip()
      printed[name] = l.split(" final score: ")[1].split()[0]
      groups[-1][1].append(name)
  assert sorted(printed) == sorted(fresh), (printed, fresh)
  for name in fresh:
    assert printed[name] == "{:.2f}".format(fresh[name]), (
      name, printed[name], fresh[name])
  for total, names in groups:
    if total is not None:
      assert total == "{:.2f}".format(sum(fresh[n] for n in names)), (
        total, names)

def test_stale_partner():
  desc = run("tst_reduce2_fresh_scores_4wsc", lys_gln_str,
    ["add_flip_movers=True"])
  states = Recorder.states
  gln = "AmideFlip at chain C GLN 319 NE2"
  # stored -11.76, fresh -16.07 (printed was the stored one)
  assert states[gln][0] - states[gln][1] > 1, states[gln]
  assert_fresh(desc, dict((n, s[1]) for n, s in states.items()))

if __name__ == "__main__":
  test_stale_partner()
  print("OK")
