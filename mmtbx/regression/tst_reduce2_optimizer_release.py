"""reduce2: the Optimizer releases its search data at the end of construction.

After the search the Optimizer drops its OptimizerC (with the dot cache), the
Movers and the dictionaries keyed by them, the interaction graph, the excluded
atoms, the dot scorer and the spatial query; the public getters return what
they returned with these kept.

Tests:
- test_release: Lys A 43 + Gln D 255 (I 2 2 2, symmetry copies, an amide flip)
  and His A 354 + Cd (P 21 21 21, a deleted HE2): after construction the search
  attributes are empty or None and no OptimizerC or Mover is alive; each getter
  returns its known value.
"""
from __future__ import absolute_import, division, print_function
import gc
import weakref
from io import StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers
from mmtbx.regression.tst_reduce2_symmetry_copies import lys_gln_str, his_cd_str

_Optimizer = Optimizers.Optimizer
_OptimizerC = Optimizers.OptimizerC

RELEASED = ["_optimizerOf", "_movers", "_moverInfo", "_initialScore",
  "_finalScore", "_interactionGraph", "_atomMoverLists", "_excludeDict",
  "_dotScorer", "_spatialQuery"]

class Tracked(object):
  """OptimizerC that keeps weak references to itself and the Movers it
  optimizes."""
  refs = []
  def __init__(self, *args):
    self._r = _OptimizerC(*args)
    Tracked.refs.append(weakref.ref(self._r))
  def __getattr__(self, n):
    return getattr(self._r, n)
  def OptimizeSingleMoverCoarse(self, m):
    Tracked.refs.append(weakref.ref(m))
    return self._r.OptimizeSingleMoverCoarse(m)
  def OptimizeCliqueCoarse(self, movers, *args):
    for m in movers:
      Tracked.refs.append(weakref.ref(m))
    return self._r.OptimizeCliqueCoarse(movers, *args)

def residue(a):
  """chain, resname and resseq of the atom's residue (what mmtbx.hydrogens
  reads from a flipped Mover's base atom)."""
  l = a.fetch_labels()
  return (l.chain_id, l.resname, l.resseq.strip())

class Recorder(_Optimizer):
  """Right after construction: the search attributes, the number of live
  search objects and the values of the getters that do not clear what they
  return; the clearing getters' values when reduce2 calls them."""
  built = []
  got = {}
  def __init__(self, *args, **kw):
    _Optimizer.__init__(self, *args, **kw)
    gc.collect()
    Recorder.built.append(dict(
      attributes=dict((n, getattr(self, n)) for n in RELEASED),
      tracked=len(Tracked.refs),
      alive=len([r for r in Tracked.refs if r() is not None])))
    # atoms by id: reduce2 rebuilds the model afterwards
    self._record("delete", [a.id_str() for a in self.getHydrogensToDelete()])
    self._record("amides", [(x.alt, residue(x.baseAtom))
      for x in self.getFlippedAmides()])
    self._record("histidines", [(x.alt, residue(x.baseAtom))
      for x in self.getFlippedHistidines()])
    self._record("copies", [(c.id_str(), o.id_str(), str(rt), present)
      for c, o, rt, present in self.getSymmetryCopies()])
    self._record("summary", self.getSymmetryCopiesSummary())
  def _record(self, name, value):
    Recorder.got.setdefault(name, []).append(value)
    return value
  def getInfo(self):
    return self._record("info", _Optimizer.getInfo(self))
  def getWarnings(self):
    return self._record("warnings", _Optimizer.getWarnings(self))
  def getAtomDump(self):
    return self._record("dump", _Optimizer.getAtomDump(self))

def run(prefix, pdb_text):
  """reduce2 (flips, atom dump) on pdb_text; returns the Optimizer's record
  and the getters' values."""
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  Recorder.built = []
  Recorder.got = {}
  Tracked.refs = []
  Optimizers.Optimizer = Recorder
  Optimizers.OptimizerC = Tracked
  try:
    run_program(program_class=reduce2.Program, logger=StringIO(), args=[
      prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
      "output.description_file_name=" + prefix + ".txt", "overwrite=True",
      "add_flip_movers=True", "print_atom_info=True"])
  finally:
    Optimizers.Optimizer = _Optimizer
    Optimizers.OptimizerC = _OptimizerC
  assert len(Recorder.built) == 1, Recorder.built
  built = Recorder.built[0]
  for n, v in built["attributes"].items():
    assert v is None or len(v) == 0, (n, v)
  assert built["tracked"] > 1 and built["alive"] == 0, built
  return Recorder.got

def report(info):
  """The lines of the info string from BEGIN REPORT to END REPORT."""
  lines = info.splitlines()
  return lines[lines.index(" BEGIN REPORT: Model 0 Alt '':"):
    lines.index(" END REPORT") + 1]

def test_release():
  got = run("tst_reduce2_optimizer_release_lysgln", lys_gln_str)
  info = got["info"][0]
  assert " Inserted 4 Movers" in info, info
  assert report(info) == [
    " BEGIN REPORT: Model 0 Alt '':",
    "  Set of 2 Movers:  Totals: initial score -7.07, final score 1.84",
    "   AmideFlip at chain D GLN 255 NE2 Initial score: -7.39 final score: 1.42"
    " pose Flipped . . AnglesAdjusted .",
    "   SingleHydrogenRotator at chain D GLN 255 H Initial score: 0.32 final"
    " score: 0.42 pose Angle -174.0 deg . .",
    "  Singleton Movers:",
    "   NH3Rotator at chain A LYS 43 NZ Initial score: -8.48 final score: 1.71"
    " pose Angle 177.0 deg . .",
    "   SingleHydrogenRotator at chain A LYS 43 H Initial score: 0.50 final"
    " score: 0.55 pose Angle 97.0 deg . .",
    " END REPORT"], report(info)
  assert got["warnings"][0].count("Overriding radius") == 3, got["warnings"]
  # one line per atom (no copies); the Gln amide flipped
  dump = got["dump"][0].splitlines()
  assert len(dump) == 39, dump
  assert [l for l in dump if " GLN 255 " in l and " OE1 " in l][0] == (
    " D GLN 255  OE1 -  48.291  28.901  57.385  1.40 isAcceptor noDonor"
    " noMetallic"), dump
  assert got["amides"][0] == [("", ("D", "GLN", "255"))], got["amides"]
  assert got["histidines"][0] == [], got
  assert got["delete"][0] == [], got
  copies = got["copies"][0]
  assert len(copies) == 39 and all(present for c, o, rt, present in copies), (
    copies)
  assert got["summary"][0] == ("Symmetry copies (fixed environment): model 1,"
    " 39 copy atoms, 3 passes (pass 2: 2 of 3 Cliques re-optimized, 2 Movers"
    " changed; pass 3: 2 of 3 Cliques re-optimized, 2 Movers changed);"
    " unsettled: NH3Rotator at chain A LYS 43 NZ; AmideFlip at chain D GLN 255"
    " NE2 (best of 4 combinations, score 3.54)\n"), got["summary"]
  got = run("tst_reduce2_optimizer_release_hiscd", his_cd_str)
  assert got["delete"][0] == ['pdb=" HE2 HIS A 354 "'], got["delete"]
  assert len(got["copies"][0]) == 19, got["copies"]
  print("test_release OK")

if __name__ == "__main__":
  test_release()
  print("OK")
