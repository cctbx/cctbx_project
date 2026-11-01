"""reduce2: the Optimizer's excluded atoms for the atoms of Movers, by bond count.

Each atom of a Mover is scored without the atoms within a few bonds of it. The
lists were made with every Mover at its first state and kept an atom only if it
was close enough to touch there, so a rotating H was scored against its own
atoms a few bonds away wherever it came closer than at that state: the three
equivalent poses of an NH3 scored differently, and a hydroxyl H was credited
for contacts with its own atoms that probe2 never counts. The lists, and
placement's list of atoms that cannot be candidate acceptors, are now by bond
count, whatever the distance.

Tests:
- Lys A 11 NZ (1us0, the residue alone): the three equivalent poses of the NH3
  (rel, +120, +240 degrees) score the same;
- CIT A 450 (1us0, the residue alone): the printed score of HO7 equals a
  rescoring at its final pose with the excluded atoms found there (probe2's
  rule);
- Tyr A 48 (1us0, the residue alone): no candidate acceptor of HH is within 4
  bonds of it.
"""
from __future__ import absolute_import, division, print_function
from io import StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.probe import Helpers
from mmtbx.reduce import Movers, Optimizers

# 1us0 Lys A 11, the residue alone
lys_str = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM     80  N   LYS A  11       4.154   9.909  32.720  1.00  6.20           N
ATOM     81  CA  LYS A  11       5.519   9.835  33.152  1.00  6.13           C
ATOM     82  C   LYS A  11       6.078   8.411  32.936  1.00  5.87           C
ATOM     83  O   LYS A  11       5.394   7.416  33.169  1.00  7.52           O
ATOM     84  CB  LYS A  11       5.620  10.187  34.632  1.00  8.59           C
ATOM     85  CG  LYS A  11       5.328  11.681  34.832  1.00 13.51           C
ATOM     86  CD  LYS A  11       5.421  11.992  36.312  1.00 23.56           C
ATOM     87  CE  LYS A  11       5.364  13.488  36.550  1.00 30.54           C
ATOM     88  NZ  LYS A  11       5.446  13.601  38.091  1.00 40.48           N
END
"""

# 1us0 CIT A 450, the residue alone
cit_str = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
HETATM 3054  C1  CIT A 450       6.260 -16.763  24.168  0.51  3.55           C
HETATM 3055  C2  CIT A 450       6.050 -17.259  22.751  0.51  4.50           C
HETATM 3056  C3  CIT A 450       4.553 -17.365  22.401  0.51  4.83           C
HETATM 3057  C4  CIT A 450       4.431 -18.117  21.078  0.51  5.20           C
HETATM 3058  C5  CIT A 450       3.007 -18.440  20.727  0.51  5.80           C
HETATM 3059  C6  CIT A 450       3.988 -15.955  22.237  0.51  4.77           C
HETATM 3060  O1  CIT A 450       5.958 -15.573  24.450  0.51  4.07           O
HETATM 3061  O2  CIT A 450       6.743 -17.602  24.951  0.51  3.60           O
HETATM 3062  O3  CIT A 450       2.178 -17.525  20.672  0.51  5.95           O
HETATM 3063  O4  CIT A 450       2.767 -19.628  20.386  0.51  7.53           O
HETATM 3064  O5  CIT A 450       3.064 -15.660  23.074  0.51  5.33           O
HETATM 3065  O6  CIT A 450       4.423 -15.171  21.396  0.51  5.42           O
HETATM 3066  O7  CIT A 450       3.823 -18.019  23.401  0.51  5.58           O
END
"""

# measure() swaps Optimizers.Optimizer for Recorder while reduce2 runs
_Optimizer = Optimizers.Optimizer

# 1us0 Tyr A 48, the residue alone
tyr_str = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM    382  N   TYR A  48      10.681 -10.951  20.428  1.00  2.77           N
ATOM    383  CA  TYR A  48      11.591 -11.568  21.379  1.00  2.81           C
ATOM    384  C   TYR A  48      11.121 -12.959  21.821  1.00  2.98           C
ATOM    385  O   TYR A  48      11.908 -13.690  22.434  1.00  3.61           O
ATOM    386  CB  TYR A  48      11.786 -10.678  22.610  1.00  2.96           C
ATOM    387  CG  TYR A  48      12.389  -9.325  22.308  1.00  2.86           C
ATOM    388  CD1 TYR A  48      12.993  -9.026  21.089  1.00  2.90           C
ATOM    389  CD2 TYR A  48      12.329  -8.311  23.275  1.00  2.96           C
ATOM    390  CE1 TYR A  48      13.482  -7.753  20.820  1.00  2.88           C
ATOM    391  CE2 TYR A  48      12.825  -7.035  23.021  1.00  3.03           C
ATOM    392  CZ  TYR A  48      13.392  -6.754  21.777  1.00  2.87           C
ATOM    393  OH  TYR A  48      13.834  -5.478  21.524  1.00  3.25           O
END
"""

class Recorder(_Optimizer):
  """After the fine step: NH3 scores at three equivalent poses (the search's
  excluded atoms), and for single-H rotators the stored score and a rescoring
  with the excluded atoms found at the final pose."""
  results = {}

  def _PlaceMovers(self, atoms, rotatableHydrogenIDs, bondedNeighborLists,
                   hParameters, addFlipMovers, alt):
    self._recNeighbors = bondedNeighborLists
    self._recAlt = alt
    return _Optimizer._PlaceMovers(self, atoms, rotatableHydrogenIDs,
      bondedNeighborLists, hParameters, addFlipMovers, alt)

  def _recDots(self, optC, atoms, excluded):
    total = 0.0
    for a in atoms:
      if a in self._deleteMes:
        continue
      r = self._extraAtomInfo.getMappingFor(a).vdwRadius + self._maximumVDWRadius
      total += self._dotScorer.score_dots(a, self._minOccupancy, self._spatialQuery, r,
        self._probeRadius, excluded(a), optC.GetDots(a.i_seq), self._probeDensity,
        False, False).totalScore()
    return total

  def _recPose(self, optC, m):
    coarse, fine = optC.GetCoarseLocation(m), optC.GetFineLocation(m)
    pr, index = m.CoarsePositions(), coarse
    if fine >= 0:
      pr, index = m.FinePositions(coarse), fine
    atoms = [pr.atoms[i] for i in range(len(pr.positions[index]))]
    return atoms, self._preferenceMagnitude * pr.preferenceEnergies[index]

  def _OptimizeFine(self, optC, cliques):
    _Optimizer._OptimizeFine(self, optC, cliques)
    searched = lambda a: self._excludeDict[a.i_seq]
    here = lambda a: Helpers.getAtomsWithinNBonds(a, self._recNeighbors,
      self._extraAtomInfo, self._probeRadius, self._bondedNeighborDepth, 3)
    for m in self._movers:
      atoms, pref = self._recPose(optC, m)
      rec = dict(stored=optC.GetHighScore(m),
                 here=pref + self._recDots(optC, atoms, here))
      if isinstance(m, Movers.MoverNH3Rotator):
        # _posesFor() rotates from the atoms' current (final) positions
        start = [a.xyz for a in atoms]
        rec["equivalent"] = []
        for pose in m._posesFor([0.0, 120.0, 240.0]):
          for a, xyz in zip(atoms, pose):
            self._spatialQuery.remove(a)
            a.xyz = xyz
            self._spatialQuery.add(a)
          rec["equivalent"].append(self._recDots(optC, atoms, searched))
        for a, xyz in zip(atoms, start):
          self._spatialQuery.remove(a)
          a.xyz = xyz
          self._spatialQuery.add(a)
      Recorder.results[(self._recAlt, self._MoverName(m))] = rec

def measure(prefix, pdb_text, extra_args=()):
  """reduce2 add on pdb_text with the Recorder; returns (results, description)."""
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  Recorder.results = {}
  Recorder.acceptors = {}
  rotator = Movers.MoverSingleHydrogenRotator.__init__
  def record(self, atom, bondedNeighborLists, extraAtomInfoMap, hParameters,
             potentialAcceptors=[], *args, **kw):
    Recorder.acceptors[atom.name.strip()] = sorted(
      a.name.strip() for a in potentialAcceptors)
    rotator(self, atom, bondedNeighborLists, extraAtomInfoMap, hParameters,
            potentialAcceptors, *args, **kw)
  Movers.MoverSingleHydrogenRotator.__init__ = record
  Optimizers.Optimizer = Recorder
  try:
    run_program(program_class=reduce2.Program, logger=StringIO(), args=[
      prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
      "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
      + list(extra_args))
  finally:
    Optimizers.Optimizer = _Optimizer
    Movers.MoverSingleHydrogenRotator.__init__ = rotator
  with open(prefix + ".txt") as f:
    description = f.read()
  return dict(Recorder.results), description

def printed_score(description, name):
  lines = [l for l in description.splitlines()
           if l.strip().startswith(name + " Initial score:")]
  assert len(lines) == 1, (name, lines)
  return float(lines[0].split(" final score: ")[1].split()[0])

def test_nh3_equivalent_poses():
  results, _ = measure("tst_reduce2_exclusions_lys", lys_str)
  s = results[("", "NH3Rotator at chain A LYS 11 NZ")]["equivalent"]
  # was 0.18 apart: each HZ excluded the atoms near it at the start only
  assert max(s) - min(s) < 1e-6, s

def test_hydroxyl_own_atoms():
  results, description = measure("tst_reduce2_exclusions_cit", cit_str)
  name = "SingleHydrogenRotator at chain A CIT 450 HO7"
  rec = results[("", name)]
  # was 3.31 printed, 0.00 rescored: HO7 1.89 A from its own O5, 4 bonds away
  assert abs(printed_score(description, name) - rec["here"]) < 0.006, (
    printed_score(description, name), rec)

def test_candidate_acceptors():
  measure("tst_reduce2_exclusions_tyr", tyr_str)
  # was CD1 and CD2, 4 bonds away, for HH placed away from them
  assert not set(Recorder.acceptors["HH"]) & set(
    ["CZ", "CE1", "CE2", "CD1", "CD2"]), Recorder.acceptors

if __name__ == "__main__":
  test_nh3_equivalent_poses()
  test_hydroxyl_own_atoms()
  test_candidate_acceptors()
  print("OK")
