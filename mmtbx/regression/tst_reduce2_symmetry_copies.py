"""reduce2 include_symmetric_copies: symmetry mates within 3.6 A of the model,
from the input's crystal symmetry, as a fixed environment of the Optimizer.

- Ser + Asp, P 1 2 1: Ser OG is 2.65 A from a carboxylate O of the Asp mate
  -x+1,y,-z+1 and nothing else is within 12 A. With copies the mate's O is a
  candidate acceptor of the Ser HG rotator and HG H-bonds to it; without, HG
  sees no acceptor.
- 1XMK His A 354 and Cd A 397, P 21 21 21: NE2 is 2.3 A from the Cd of a
  mate; with copies the His is locked by the ion rule, without it is a flip
  Mover.
- 1R27 Lys A 43 and Gln D 255, I 2 2 2: each faces the other's mate, and the
  two do not settle in 3 passes; the best scoring combination of their
  states is kept, both are listed as unsettled and their reported scores are
  fresh scores with the final copies. Past the time limit no refresh pass
  starts; a pass the time limit overtakes keeps the better of the previous
  and new states.
- Bonds across symmetry map by their operator: two images of a partner
  within 3 A (P 1, a = 5 A), and a metal on a twofold axis whose copy is left
  out as coincident (P 1 2 1), bonded both ways to the copy of its ligand.
"""
from __future__ import absolute_import, division, print_function
import itertools
import time
from io import StringIO
import iotbx.pdb
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from cctbx import crystal, sgtbx
from scitbx import matrix
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers
from mmtbx.probe import Helpers
import mmtbx_probe_ext as probeExt

ser_asp_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   SER A   1      24.099  27.165  12.526  1.00 20.00           N
ATOM      2  CA  SER A   1      22.669  27.246  12.854  1.00 20.00           C
ATOM      3  C   SER A   1      22.230  28.687  12.828  1.00 20.00           C
ATOM      4  O   SER A   1      23.013  29.563  13.105  1.00 20.00           O
ATOM      5  CB  SER A   1      22.432  26.665  14.248  1.00 20.00           C
ATOM      6  OG  SER A   1      22.850  25.299  14.272  1.00 20.00           O
ATOM      7  OXT SER A   1      20.967  28.997  12.497  1.00 20.00           O
ATOM      8  N   ASP B   1       5.600  20.940  11.664  1.00 20.00           N
ATOM      9  CA  ASP B   1       7.082  21.049  11.550  1.00 20.00           C
ATOM     10  C   ASP B   1       7.764  20.210  12.632  1.00 20.00           C
ATOM     11  O   ASP B   1       8.982  20.391  12.840  1.00 20.00           O
ATOM     12  CB  ASP B   1       7.499  22.520  11.683  1.00 20.00           C
ATOM     13  CG  ASP B   1       7.251  23.076  13.088  1.00 20.00           C
ATOM     14  OD1 ASP B   1       7.554  24.272  13.319  1.00 20.00           O
ATOM     15  OD2 ASP B   1       6.766  22.312  13.954  1.00 20.00           O
ATOM     16  OXT ASP B   1       7.061  19.382  13.265  1.00 20.00           O
END
"""

his_cd_str = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21    4
ATOM    548  N   HIS A 354      -1.093  10.746  39.862  1.00  6.48           N
ATOM    549  CA  HIS A 354      -0.970   9.793  40.975  1.00  6.84           C
ATOM    550  C   HIS A 354       0.400   9.117  40.895  1.00  6.70           C
ATOM    551  O   HIS A 354       0.877   8.812  39.808  1.00  7.52           O
ATOM    552  CB  HIS A 354      -2.004   8.651  40.833  1.00  8.27           C
ATOM    553  CG  HIS A 354      -3.476   9.048  40.932  1.00  8.55           C
ATOM    554  ND1 HIS A 354      -4.160   9.785  40.003  1.00  8.74           N
ATOM    555  CD2 HIS A 354      -4.350   8.756  41.946  1.00  8.81           C
ATOM    556  CE1 HIS A 354      -5.447   9.936  40.440  1.00  8.37           C
ATOM    557  NE2 HIS A 354      -5.605   9.331  41.604  1.00  8.44           N
HETATM  666 CD    CD A 397      10.242  12.438  48.077  1.00  9.61          CD
END
"""

lys_gln_str = """\
CRYST1  127.401  298.300  296.853  90.00  90.00  90.00 I 2 2 2
ATOM    139  N   LYS A  43      -9.617 122.553  96.993  1.00 24.63           N
ATOM    140  CA  LYS A  43     -10.111 121.191  97.128  1.00 24.28           C
ATOM    141  C   LYS A  43      -9.013 120.194  96.839  1.00 22.43           C
ATOM    142  O   LYS A  43      -7.997 120.523  96.231  1.00 20.41           O
ATOM    143  CB  LYS A  43     -11.265 120.919  96.160  1.00 26.41           C
ATOM    144  CG  LYS A  43     -10.852 120.973  94.708  1.00 31.17           C
ATOM    145  CD  LYS A  43     -11.847 120.263  93.782  1.00 33.48           C
ATOM    146  CE  LYS A  43     -11.452 120.497  92.322  1.00 34.50           C
ATOM    147  NZ  LYS A  43     -12.293 119.755  91.342  1.00 36.26           N
TER
ATOM  25076  N   GLN D 255      46.498  25.824  58.503  1.00 28.70           N
ATOM  25077  CA  GLN D 255      45.675  26.359  57.422  1.00 28.42           C
ATOM  25078  C   GLN D 255      44.343  25.620  57.320  1.00 26.75           C
ATOM  25079  O   GLN D 255      44.199  24.517  57.846  1.00 25.56           O
ATOM  25080  CB  GLN D 255      46.421  26.245  56.090  1.00 31.27           C
ATOM  25081  CG  GLN D 255      47.742  26.986  56.064  1.00 35.46           C
ATOM  25082  CD  GLN D 255      47.568  28.446  56.394  1.00 38.75           C
ATOM  25083  OE1 GLN D 255      46.767  29.140  55.770  1.00 40.34           O
ATOM  25084  NE2 GLN D 255      48.316  28.925  57.385  1.00 40.98           N
TER
END
"""

cryst1_placeholder = \
  "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"

def run(prefix, pdb_text, extra_args=(), check_copies=False,
        optimizer_class=None):
  """reduce2 on pdb_text; returns (output model text, description, log).
  check_copies: assert at the end of each Optimizer with copies that every
  copy of the last alternate (and model) sits at its original's symmetry
  mate and is present unless the original is deleted.
  optimizer_class: a subclass of Optimizers.Optimizer to run instead."""
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  log = StringIO()
  base = Optimizers.Optimizer
  parent = optimizer_class if optimizer_class is not None else base
  class Checked(parent):
    def __init__(self, *args, **kw):
      parent.__init__(self, *args, **kw)
      if check_copies and kw.get("symmetryCopies") is not None:
        uc = kw["symmetryCopies"].crystal_symmetry.unit_cell()
        deleted = set(self.getHydrogensToDelete())
        for c, o, rt, present in self.getSymmetryCopies():
          assert present == (o not in deleted), o.id_str()
          mate = uc.orthogonalize(rt * uc.fractionalize(o.xyz))
          assert approx_equal(c.xyz, mate, eps=1.e-3, out=None), o.id_str()
  Optimizers.Optimizer = Checked
  try:
    run_program(program_class=reduce2.Program, logger=log, args=[
      prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
      "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
      + list(extra_args))
  finally:
    Optimizers.Optimizer = base
  with open(prefix + "_out.pdb") as f:
    model_text = f.read()
  with open(prefix + ".txt") as f:
    description = f.read()
  return model_text, description, log.getvalue()

def nearest_mate(cs, x, y):
  """Distance from x to the nearest symmetry mate of y (all operators and
  neighboring cells)."""
  uc = cs.unit_cell()
  fy = uc.fractionalize(y)
  best = None
  for op in cs.space_group().all_ops():
    f = op * fy
    for t in itertools.product((-1, 0, 1), repeat=3):
      d = (matrix.col(uc.orthogonalize([f[i] + t[i] for i in range(3)])) -
           matrix.col(x)).length()
      if best is None or d < best:
        best = d
  return best

def hg_to_asp(model_text, model_index=0):
  """Ser HG to the nearest mate of an Asp carboxylate O, in one MODEL."""
  inp = iotbx.pdb.input(lines=model_text.split("\n"), source_info=None)
  m = inp.construct_hierarchy().models()[model_index]
  a = dict(((x.parent().parent().parent().id, x.name.strip()), x)
           for x in m.atoms())
  return min(nearest_mate(inp.crystal_symmetry(), a[("A", "HG")].xyz,
    a[("B", n)].xyz) for n in ("OD1", "OD2"))

def summary_lines(text):
  return [l.strip() for l in text.split("\n")
          if l.strip().startswith("Symmetry copies (fixed environment)")]

def test_hbond_across_contact():
  on, desc_on, log_on = run("tst_reduce2_symmetry_copies_1", ser_asp_str,
    ["include_symmetric_copies=True"], check_copies=True)
  off, desc_off, _ = run("tst_reduce2_symmetry_copies_2", ser_asp_str)
  # the mate's O is a candidate acceptor of the HG rotator
  assert "SER 1 HG with 1 potential nearby acceptors" in desc_on, desc_on
  assert "SER 1 HG with 0 potential nearby acceptors" in desc_off
  assert hg_to_asp(on) < 2.0, hg_to_asp(on)
  assert hg_to_asp(off) > 2.5, hg_to_asp(off)
  lines = summary_lines(log_on)
  assert len(lines) == 1, lines
  # copies changed (H turned by placement), none near a Mover's reach
  assert lines[0] == ("Symmetry copies (fixed environment): model 1, 29 copy "
    "atoms, 2 passes (pass 2: 0 of 2 Cliques re-optimized, 0 Movers changed); "
    "unsettled: none"), lines
  assert summary_lines(desc_on) == lines
  assert summary_lines(desc_off) == []

def test_option_off_identical():
  off, _, _ = run("tst_reduce2_symmetry_copies_3", ser_asp_str,
    ["include_symmetric_copies=False"])
  default, _, _ = run("tst_reduce2_symmetry_copies_4", ser_asp_str)
  assert off == default

def test_no_copies_without_symmetry():
  atoms = "".join(l + "\n" for l in ser_asp_str.split("\n")
                  if l.startswith("ATOM")) + "END\n"
  for i, text in enumerate([cryst1_placeholder + atoms, atoms]):
    prefix = "tst_reduce2_symmetry_copies_%d" % (5 + 2 * i)
    on, desc, log = run(prefix, text, ["include_symmetric_copies=True"])
    off, _, _ = run(prefix + "_off", text)
    assert summary_lines(desc) == [] and summary_lines(log) == []
    assert on == off

def test_ion_lock_through_copy():
  _, desc_on, _ = run("tst_reduce2_symmetry_copies_9", his_cd_str,
    ["include_symmetric_copies=True", "add_flip_movers=True"],
    check_copies=True)
  _, desc_off, _ = run("tst_reduce2_symmetry_copies_10", his_cd_str,
    ["add_flip_movers=True"])
  assert ("Not adding Hydrogen to chain A HIS 354 NE2 and marking as an "
    "acceptor (ionic bond to CD)") in desc_on, desc_on
  assert "Added MoverHisFlip" not in desc_on
  assert "Added MoverHisFlip 1 to chain A HIS 354" in desc_off
  assert "ionic bond" not in desc_off

def test_multi_model():
  # Model 1 has the Asp 6 A away along y: no contact, no copies. Model 2
  # has the contact and its own copies.
  h = iotbx.pdb.input(lines=ser_asp_str.split("\n"), source_info=None
    ).construct_hierarchy()
  m2 = h.models()[0].detached_copy()
  m2.id = "2"
  h.models()[0].id = "1"
  for a in h.models()[0].atoms():
    if a.parent().parent().parent().id == "B":
      a.xyz = (a.xyz[0], a.xyz[1] - 6.0, a.xyz[2])
  h.append_model(m2)
  cs = crystal.symmetry((30, 30, 30, 90, 90, 90), "P 1 2 1")
  out, desc, log = run("tst_reduce2_symmetry_copies_11",
    h.as_pdb_string(crystal_symmetry=cs), ["include_symmetric_copies=True"],
    check_copies=True)
  assert hg_to_asp(out, 0) > 2.5, hg_to_asp(out, 0)
  assert hg_to_asp(out, 1) < 2.0, hg_to_asp(out, 1)
  # a line per model with copies
  lines = summary_lines(log)
  assert len(lines) == 1, lines
  assert lines[0].startswith(
    "Symmetry copies (fixed environment): model 2, 29 copy atoms"), lines

# run() replaces Optimizers.Optimizer while it runs
Optimizer = Optimizers.Optimizer

class FreshScores(Optimizer):
  """Records each Mover's coarse state and a fresh score of its pose with
  the copies as they are after the refresh, by Mover name, and asserts that
  its atoms are at the pose its OptimizerC holds."""
  def _RefreshSymmetryCopies(self, *args):
    Optimizer._RefreshSymmetryCopies(self, *args)
    dot_scorer = Helpers.createDotScorer(self._extraAtomInfo, self._probePhil)
    self.states = {}
    for m in self._movers:
      optC = self._optimizerOf[m]
      coarse = optC.GetCoarseLocation(m)
      fine = optC.GetFineLocation(m)
      pr, index = m.CoarsePositions(), coarse
      if fine >= 0:
        pr, index = m.FinePositions(coarse), fine
      for a, xyz in zip(pr.atoms, pr.positions[index]):
        assert approx_equal(a.xyz, xyz, eps=1.e-6, out=None), (
          self._MoverName(m), a.name, coarse, fine)
      self.states[self._MoverName(m)] = (coarse, self._ScoreMover(m, dot_scorer))
    FreshScores.last = self.states

def assert_fresh_scores(desc, states):
  """Every printed final score is the Mover's fresh score, and every Clique's
  final total the sum of its Movers' fresh scores, to 2 decimals."""
  printed = {}
  groups = []                     # (printed final total, Mover names)
  for l in desc.split("\n"):
    if "Totals:" in l:
      groups.append((l.split("final score ")[1].strip(), []))
    elif "Singleton Movers:" in l:
      groups.append((None, []))
    elif " final score: " in l:
      name = l.split(" Initial score:")[0].strip()
      printed[name] = l.split(" final score: ")[1].split()[0]
      groups[-1][1].append(name)
  assert sorted(printed) == sorted(states), (printed, states)
  for name in states:
    assert printed[name] == "{:.2f}".format(states[name][1]), (
      name, printed[name], states[name])
  for total, names in groups:
    if total is not None:
      assert total == "{:.2f}".format(sum(states[n][1] for n in names)), (
        total, names, states)

def test_unsettled_pair():
  _, desc, log = run("tst_reduce2_symmetry_copies_12", lys_gln_str,
    ["include_symmetric_copies=True", "add_flip_movers=True"],
    check_copies=True, optimizer_class=FreshScores)
  lines = summary_lines(log)
  assert len(lines) == 1, lines
  # only the two Cliques facing the other's copy are optimized again
  assert ("3 passes (pass 2: 2 of 3 Cliques re-optimized, 2 Movers changed; "
    "pass 3: 2 of 3 Cliques re-optimized, 2 Movers changed)") in lines[0], lines
  assert lines[0].endswith("unsettled: NH3Rotator at chain A LYS 43 NZ; "
    "AmideFlip at chain D GLN 255 NE2 (best of 4 combinations, score 4.34)"
    ), lines
  assert summary_lines(desc) == lines
  # Reported scores are fresh scores at the chosen poses with matching
  # copies; the OptimizerC of the Gln's pose had 1.79.
  states = FreshScores.last
  assert approx_equal(states["AmideFlip at chain D GLN 255 NE2"][1], 1.7445,
    eps=1.e-3), states
  assert_fresh_scores(desc, states)

class DeadlineBeforeRefresh(FreshScores):
  """The time limit has passed when the refresh starts."""
  def _RefreshSymmetryCopies(self, singletons, groups, cache, stop, alt):
    FreshScores._RefreshSymmetryCopies(self, singletons, groups, cache, 1.0,
      alt)

class LateDeadline(FreshScores):
  """The time limit expires after the first coarse optimization (in a
  refresh pass, after Initialize())."""
  def _OptimizeCliques(self, optC, singletons, groups, stopOptimizingAt):
    self.searches = getattr(self, "searches", 0) + 1
    if self.searches > 1:
      stopOptimizingAt = 1.0
    return Optimizer._OptimizeCliques(self, optC, singletons,
      groups, stopOptimizingAt)

class ExpiredClock(object):
  """time for Optimizers: past every time limit."""
  def time(self):
    return time.time() + 1.e11

class ExpiresInRefresh(FreshScores):
  """From the first refresh pass's coarse optimization on (after
  Initialize()), Optimizers' clock is past the time limit, for the refresh
  and the coarse solver alike."""
  def _OptimizeCliques(self, optC, singletons, groups, stopOptimizingAt):
    self.searches = getattr(self, "searches", 0) + 1
    if self.searches == 2:
      Optimizers.time = ExpiredClock()
    return Optimizer._OptimizeCliques(self, optC, singletons,
      groups, stopOptimizingAt)

def test_time_limit_in_refresh():
  # The main search leaves Lys 43 NZ in coarse state 7 and Gln 255 flipped.
  # Past the time limit no refresh pass starts: both keep those states, and
  # their scores, from before the copies followed them (Lys 1.55, Gln 1.79),
  # are reported fresh.
  _, desc, log = run("tst_reduce2_symmetry_copies_13", lys_gln_str,
    ["include_symmetric_copies=True", "add_flip_movers=True"],
    check_copies=True, optimizer_class=DeadlineBeforeRefresh)
  states = FreshScores.last
  lys = "NH3Rotator at chain A LYS 43 NZ"
  gln = "AmideFlip at chain D GLN 255 NE2"
  assert states[lys][0] == 7 and states[gln][0] == 1, states
  lines = summary_lines(log)
  assert len(lines) == 1, lines
  assert lines[0].endswith("1 pass (stopped at pass 2 due to time limit); "
    "unsettled: none"), lines
  assert approx_equal(states[lys][1], 0.4229, eps=1.e-3), states
  assert approx_equal(states[gln][1], 0.7241, eps=1.e-3), states
  assert_fresh_scores(desc, states)
  # The time limit expires in the refresh, after Initialize(): its Cliques
  # are solved as singletons, not left in their first state.
  _, desc, log = run("tst_reduce2_symmetry_copies_14", lys_gln_str,
    ["include_symmetric_copies=True", "add_flip_movers=True"],
    check_copies=True, optimizer_class=LateDeadline)
  states = FreshScores.last
  assert desc.count("No time remaining, treating (sub)clique of size 2 as "
    "singletons") == 2, desc
  lines = summary_lines(log)
  assert len(lines) == 1, lines
  assert ("3 passes (pass 2: 2 of 3 Cliques re-optimized, 3 Movers changed; "
    "pass 3: 2 of 3 Cliques re-optimized, 1 Movers changed); unsettled: "
    "NH3Rotator at chain A LYS 43 NZ (best of 2 combinations, score 3.37)"
    ) in lines[0], lines
  assert_fresh_scores(desc, states)
  # The time limit passes after Initialize() of the first refresh pass and
  # stays passed: the pass solves its Cliques as singletons (Lys 0, Gln
  # unflipped: -14.83 in all with matching copies, the main search's states
  # 2.47). The best combination of the previous and new states of the three
  # Movers it changed is kept (Lys 0, Gln flipped, as the unsettled pair
  # without a time limit), and no further pass starts.
  try:
    _, desc, log = run("tst_reduce2_symmetry_copies_15", lys_gln_str,
      ["include_symmetric_copies=True", "add_flip_movers=True"],
      check_copies=True, optimizer_class=ExpiresInRefresh)
  finally:
    Optimizers.time = time
  states = FreshScores.last
  assert desc.count("No time remaining, treating (sub)clique of size 2 as "
    "singletons") == 1, desc
  assert states[lys][0] == 0 and states[gln][0] == 1, states
  assert approx_equal(sum(v[1] for v in states.values()), 4.7688, eps=1.e-3
    ), states
  lines = summary_lines(log)
  assert len(lines) == 1, lines
  assert lines[0].endswith("2 passes (pass 2: 2 of 3 Cliques re-optimized, "
    "3 Movers changed, time limit reached: best of previous and new states "
    "of NH3Rotator at chain A LYS 43 NZ; AmideFlip at chain D GLN 255 NE2; "
    "SingleHydrogenRotator at chain D GLN 255 H (best of 8 combinations, "
    "score 4.34)); unsettled: none"), lines
  assert_fresh_scores(desc, states)

def copies_of(pdb_text, cs):
  """Originals by i_seq, and their copies (copy, original, rt_mx) within
  3.6 A."""
  h = iotbx.pdb.input(lines=pdb_text.split("\n"), source_info=None
    ).construct_hierarchy()
  atoms = h.atoms()
  originals = dict((a.i_seq, a) for a in atoms)
  copies = []
  keep = []
  for e in h.sym_equiv_residue_groups(cs, 3.6):
    keep.append(e.hierarchy)
    for k, c in zip(e.i_seqs, e.hierarchy.atoms()):
      copies.append((c, atoms[k], e.rt_mx))
  return originals, copies, keep, h

def copy_at(copies, name, rt):
  for c, o, r in copies:
    if o.name.strip() == name and r.as_xyz() == sgtbx.rt_mx(rt).as_xyz():
      return c
  return None

def test_bond_partner_by_operator():
  # C1 is bonded to C2 moved by x-1,y,z (2.6 A); C2 moved by x-2,y,z is
  # nearer (2.4 A) but not bonded. C3 brings C1 moved by x-1,y,z into the
  # copies, the copy bonded to C2 moved by x-2,y,z.
  pdb_text = """\
HETATM    1  C1  UNL A   1       2.400  15.000  15.000  1.00 20.00           C
HETATM    2  C2  UNL A   2      10.000  15.000  15.000  1.00 20.00           C
HETATM    3  C3  UNL A   3      -2.600  18.000  15.000  1.00 20.00           C
END
"""
  cs = crystal.symmetry((5, 30, 30, 90, 90, 90), "P 1")
  originals, copies, keep, h = copies_of(pdb_text, cs)
  c1, c2 = originals[0], originals[1]
  neighbors = dict((a, []) for a in originals.values())
  bond = [(c1.i_seq, c2.i_seq, sgtbx.rt_mx("x-1,y,z"))]
  copy_bonds, original_bonds = Optimizers._SymmetryCopyBonds(originals,
    copies, neighbors, bond, cs.unit_cell())
  c2_1 = copy_at(copies, "C2", "x-1,y,z")
  c2_2 = copy_at(copies, "C2", "x-2,y,z")
  c1_1 = copy_at(copies, "C1", "x-1,y,z")
  c1_p1 = copy_at(copies, "C1", "x+1,y,z")
  assert None not in (c2_1, c2_2, c1_1, c1_p1)
  assert approx_equal(c2_1.distance(c1), 2.6)
  assert approx_equal(c2_2.distance(c1), 2.4)
  bonded_to_c1 = [b for a, b in original_bonds if a == c1]
  assert bonded_to_c1 == [c2_1], [b.xyz for b in bonded_to_c1]
  assert [b for a, b in original_bonds if a == c2] == [c1_p1]
  assert copy_bonds[c2_1] == [c1]
  # copy-copy: the copy's operator times the bond's
  assert copy_bonds[c2_2] == [c1_1]
  assert copy_bonds[c1_1] == [c2_2]

def bonded_lists(neighbors, copy_bonds, added_bonds):
  """The Optimizer's bonded-neighbor lists after Mover placement."""
  lists = dict((a, list(b)) for a, b in neighbors.items())
  lists.update(copy_bonds)
  for a, b in added_bonds:
    if b not in lists[a]:
      lists[a].append(b)
  return lists

def test_metal_on_twofold():
  # Ni on the twofold axis -x,y,-z of P 1 2 1, N 2.0 A from it: the copy of
  # N is bonded to the Ni itself, the Ni's own copy being left out, and the
  # Ni to the copy. O, also on the axis, is bonded to the Ni.
  pdb_text = """\
HETATM    1 NI    NI A   1       0.000   5.000   0.000  1.00 20.00          NI
HETATM    2  N   NH3 A   2       1.414   5.000   1.414  1.00 20.00           N
HETATM    3  O   HOH A   3       0.000   7.000   0.000  1.00 20.00           O
END
"""
  cs = crystal.symmetry((20, 20, 20, 90, 90, 90), "P 1 2 1")
  originals, copies, keep, h = copies_of(pdb_text, cs)
  ni, n, o = originals[0], originals[1], originals[2]
  assert [x.name.strip() for c, x, rt in copies] == ["N"], copies
  n2 = copies[0][0]
  probeExt.set_atom_i_seq(n2, 3)
  assert approx_equal(n2.distance(ni), 2.0, eps=1.e-3)
  info = probeExt.ExtraAtomInfoMap([ni, n, o, n2],
    [probeExt.ExtraAtomInfo(r) for r in (1.0, 1.55, 1.4, 1.55)])
  neighbors = {ni: [n, o], n: [ni], o: [ni]}
  # restraints with the Ni-N bond only: the copy's N-Ni bond under -x,y,-z
  # goes to the atom at that site, the Ni
  copy_bonds, original_bonds = Optimizers._SymmetryCopyBonds(originals,
    copies, neighbors, [], cs.unit_cell())
  assert copy_bonds[n2] == [ni], copy_bonds[n2]
  assert original_bonds == [(ni, n2)], original_bonds
  # exclusions both ways, and through the Ni for its other ligand
  lists = bonded_lists(neighbors, copy_bonds, original_bonds)
  def within(a, n_bonds):
    return Helpers.getAtomsWithinNBonds(a, lists, info, 0.25, n_bonds)
  assert n2 in within(ni, 1)
  assert ni in within(n2, 1)
  assert n2 in within(o, 2)
  assert o in within(n2, 2)
  # and with the bond across symmetry as well
  bond = [(ni.i_seq, n.i_seq, sgtbx.rt_mx("-x,y,-z"))]
  copy_bonds, original_bonds = Optimizers._SymmetryCopyBonds(originals,
    copies, neighbors, bond, cs.unit_cell())
  assert copy_bonds[n2] == [ni], copy_bonds[n2]
  assert original_bonds.count((ni, n2)) == 1, original_bonds
  # the reverse bond resolves to the Ni itself, already bonded to N
  assert [b for a, b in original_bonds if a == n] == [ni]

if __name__ == "__main__":
  test_hbond_across_contact()
  test_option_off_identical()
  test_no_copies_without_symmetry()
  test_ion_lock_through_copy()
  test_multi_model()
  test_unsettled_pair()
  test_time_limit_in_refresh()
  test_bond_partner_by_operator()
  test_metal_on_twofold()
  print("OK")
