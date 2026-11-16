"""reduce2 with the default add_flip_movers=False: no heavy atom moves, each
His chooses HD1, HE2 or both in its input orientation, and the flips that
would score better are listed (Suggested flips), not made.

- Four isolated His with sulfate O 2.9 A out from ND1 (A), NE2 (B), both (C)
  or CD2 and CE1 (D, only a ring flip gives those two N-H partners): with the
  default A gets HD1, B HE2, C both, D is listed; every heavy atom keeps its
  input coordinates, occupancy and B. add_flip_movers=True flips D.
- One His with a Zn 2.1 A out from CD2, so that only the flipped ring puts an
  N on the Zn: with the default it is not flipped (the ion lock-down sets a
  His only in its input orientation) and is listed with the Zn; with
  add_flip_movers=True it is flipped and the N on the Zn has no H.
- 1PQ7 Asn A 151 with the residues and waters within 6 A of its CG, no
  crystal symmetry: add_flip_movers=True flips it; the default lists it and
  keeps the input orientation, and its amide H where placement puts them.
- The reduce flags in the program description: -build and -flip are
  add_flip_movers=True, -noflip and -nobuild# add preference_magnitude
  (9999 for -noflip): a bias, so no scored flip is made, but a His whose
  metal binds only the flipped ring is flipped, and a His keeps one ring H.
- A later conformer pass that locks a His to its metal withdraws the flip
  an earlier pass suggested; a flip whose scoring fails leaves every atom as
  without scoring, and a His listed for its metal stays listed.
- Scoring the flips changes no atom; the labels name model, insertion code,
  conformer and a metal's symmetry operator; with add_flip_movers=True a His
  locked by set_flip_states gets no BothClash note.
"""
import os
import iotbx.phil
import iotbx.pdb
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from scitbx import matrix
from mmtbx.programs import reduce2
from mmtbx.regression.tst_reduce_switch import _his_fragment_str, his_zn_str

asn_str = """\
ATOM    421  N   LEU A  73      -7.147  13.115   6.192  1.00  5.16           N
ATOM    422  CA  LEU A  73      -6.991  13.138   4.756  1.00  5.01           C
ATOM    423  C   LEU A  73      -5.692  13.843   4.335  1.00  4.83           C
ATOM    424  O   LEU A  73      -5.315  13.691   3.171  1.00  5.26           O
ATOM    425  CB  LEU A  73      -8.197  13.775   4.105  1.00  6.07           C
ATOM    426  CG  LEU A  73      -9.505  13.020   4.347  1.00  6.75           C
ATOM    427  CD1 LEU A  73     -10.658  13.834   3.811  1.00  8.95           C
ATOM    428  CD2 LEU A  73      -9.531  11.659   3.666  1.00  8.07           C
ATOM    429  N   SER A  74      -5.044  14.569   5.224  1.00  5.15           N
ATOM    430  CA  SER A  74      -3.770  15.212   4.950  1.00  5.35           C
ATOM    431  C   SER A  74      -2.630  14.407   5.595  1.00  5.13           C
ATOM    432  O   SER A  74      -2.828  13.702   6.559  1.00  5.92           O
ATOM    433  CB  SER A  74      -3.807  16.619   5.552  1.00  6.31           C
ATOM    434  OG  SER A  74      -2.507  17.189   5.442  1.00  7.61           O
ATOM    977  N   VAL A 150      -2.212  17.583  -1.336  1.00  5.24           N
ATOM    978  CA  VAL A 150      -2.354  17.684   0.104  1.00  5.13           C
ATOM    979  C   VAL A 150      -3.273  16.595   0.643  1.00  5.28           C
ATOM    980  O   VAL A 150      -2.911  15.934   1.648  1.00  6.24           O
ATOM    981  CB  VAL A 150      -2.828  19.100   0.518  1.00  5.30           C
ATOM    982  CG1 VAL A 150      -3.243  19.161   1.970  1.00  6.14           C
ATOM    983  CG2 VAL A 150      -1.703  20.101   0.248  1.00  5.66           C
ATOM    984  N   ASN A 151      -4.420  16.443   0.045  1.00  5.99           N
ATOM    985  CA  ASN A 151      -5.395  15.490   0.541  1.00  5.93           C
ATOM    986  C   ASN A 151      -5.394  14.215  -0.319  1.00  5.25           C
ATOM    987  O   ASN A 151      -5.146  14.178  -1.532  1.00  6.80           O
ATOM    988  CB  ASN A 151      -6.767  16.066   0.619  1.00  7.48           C
ATOM    989  CG  ASN A 151      -6.995  17.020   1.896  1.00 10.49           C
ATOM    990  OD1 ASN A 151      -8.253  17.339   2.045  1.00 13.82           O
ATOM    991  ND2 ASN A 151      -6.225  17.224   2.711  1.00  8.56           N
ATOM    992  N   LEU A 152      -5.759  13.131   0.393  1.00  4.51           N
ATOM    993  CA  LEU A 152      -5.733  11.781  -0.120  1.00  4.26           C
ATOM    994  C   LEU A 152      -6.566  11.654  -1.404  1.00  4.17           C
ATOM    995  O   LEU A 152      -7.730  12.093  -1.454  1.00  4.73           O
ATOM    996  CB  LEU A 152      -6.345  10.891   0.973  1.00  4.25           C
ATOM    997  CG  LEU A 152      -6.361   9.390   0.693  1.00  4.40           C
ATOM    998  CD1 LEU A 152      -4.958   8.849   0.559  1.00  4.58           C
ATOM    999  CD2 LEU A 152      -7.119   8.693   1.822  1.00  5.03           C
HETATM 1675  O   HOH A 715      -9.181  13.655   0.352  1.00  5.89           O
HETATM 1689  O   HOH A 729      -6.579  19.758   3.810  1.00  7.07           O
HETATM 1715  O   HOH A 755     -10.363  16.099   0.509  1.00  9.99           O
HETATM 1749  O   HOH A 789      -6.022  15.899  -3.774  1.00  8.85           O
HETATM 1772  O   HOH A 812      -8.757  15.277   7.247  1.00 14.72           O
HETATM 1789  O   HOH A 829      -9.701  17.817  -1.686  1.00 12.45           O
HETATM 2150  O   HOH A1190     -10.179  17.208   3.200  0.50 19.62           O
END
"""

def atoms(text):
  """{(chain, resseq, name): atom} of a PDB-format model"""
  h = iotbx.pdb.input(lines=text.split("\n"), source_info=None
    ).construct_hierarchy()
  return dict(((a.chain().id, a.parent().parent().resseq.strip(),
    a.name.strip()), a) for a in h.atoms())

def run(prefix, model_str, extra_args=(), log=None):
  """(output atoms, description text); the log goes to log"""
  fn_in, fn_out, fn_txt = prefix + ".pdb", prefix + "_out.pdb", prefix + ".txt"
  with open(fn_in, "w") as f:
    f.write(model_str)
  for fn in (fn_out, fn_txt):
    if os.path.exists(fn): os.remove(fn)
  if log is None:
    log = null_out()
  run_program(program_class=reduce2.Program, logger=log,
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file=" + fn_txt] + list(extra_args))
  return atoms(open(fn_out).read()), open(fn_txt).read()

def all_atoms(text):
  """{(model, chain, resseq, icode, altloc, name): xyz} of a PDB-format
  model, H included"""
  h = iotbx.pdb.input(lines=text.split("\n"), source_info=None
    ).construct_hierarchy()
  return dict(((a.parent().parent().parent().parent().id, a.chain().id,
    a.parent().parent().resseq, a.parent().parent().icode, a.parent().altloc,
    a.name), a.xyz) for a in h.atoms())

def assert_heavy_atoms_unchanged(model_str, out):
  """every heavy atom of the input, with its coordinates, occupancy and B"""
  start = dict((k, a) for k, a in atoms(model_str).items()
    if a.element.strip() not in ("H", "D"))
  assert set(start) <= set(out), set(start) - set(out)
  for k, a in start.items():
    b = out[k]
    assert (a.xyz, a.occ, a.b) == (b.xyz, b.occ, b.b), (k, a.xyz, b.xyz)

def ring_h(out, chain):
  return "+".join(n for n in ("HD1", "HE2") if (chain, "1", n) in out)

def suggested(desc):
  """the residue lines of the Suggested flips section"""
  lines = desc.split("\n")
  i = [j for j, l in enumerate(lines) if l.startswith("Suggested flips")]
  assert len(i) == 1, i
  result = []
  for l in lines[i[0] + 1:]:
    if not l.startswith("  "): break
    if l.startswith(("  chain", "  model")): result.append(l.strip())
  return result

def test_default_is_false():
  params = iotbx.phil.parse(reduce2.master_phil_str,
    process_includes=True).extract()
  assert params.add_flip_movers is False

def test_his_states_by_default():
  s = _his_fragment_str()
  out, desc = run("tst_reduce2_his_states_default_1", s)
  assert_heavy_atoms_unchanged(s, out)
  assert [ring_h(out, c) for c in "ABC"] == ["HD1", "HE2", "HD1+HE2"], desc
  listed = suggested(desc)
  assert len(listed) == 1 and listed[0].startswith("chain D HIS 1: flipped "
    "scores better by "), listed
  assert "Use add_flip_movers=True to apply them." in desc
  out, desc = run("tst_reduce2_his_states_default_2", s,
    ["add_flip_movers=True"])
  start = atoms(s)
  # D flipped: ND1 on the input CD2 site
  assert (matrix.col(out[("D", "1", "ND1")].xyz) -
    matrix.col(start[("D", "1", "CD2")].xyz)).length() < 0.3
  assert "Suggested flips" not in desc

def test_metal_his_listed_not_flipped():
  s = his_zn_str()
  out, desc = run("tst_reduce2_his_states_default_3", s)
  assert_heavy_atoms_unchanged(s, out)
  listed = suggested(desc)
  assert len(listed) == 1, listed
  assert listed[0].startswith("chain A HIS 1: flipped scores better by "), \
    listed
  assert listed[0].endswith("; ZN A 101 binds only the flipped ring"), listed
  out, desc = run("tst_reduce2_his_states_default_4", s,
    ["add_flip_movers=True"])
  start = atoms(s)
  # the ion lock-down flips it: a ring N on the Zn, without its H
  zn = matrix.col(start[("A", "101", "ZN")].xyz)
  d = dict((n, (matrix.col(out[("A", "1", n)].xyz) - zn).length())
    for n in ("ND1", "NE2"))
  assert min(d.values()) < 2.3, d
  near = min(d, key=d.get)
  assert ("A", "1", {"ND1": "HD1", "NE2": "HE2"}[near]) not in out, near

def test_asn_listed_not_flipped():
  start = atoms(asn_str)
  out, desc = run("tst_reduce2_his_states_default_5", asn_str)
  assert_heavy_atoms_unchanged(asn_str, out)
  listed = suggested(desc)
  # the gain over the pose kept (its amide H as placed), not over a newly
  # built amide Mover's own H positions (that gave 6.42)
  assert "chain A ASN 151: flipped scores better by 3.19" in listed, listed
  # scoring the flip leaves the amide H where placement put them
  placed, _ = run("tst_reduce2_his_states_default_7", asn_str,
    ["approach=place"])
  for name in ("HD21", "HD22"):
    k = ("A", "151", name)
    assert out[k].xyz == placed[k].xyz, (name, out[k].xyz, placed[k].xyz)
  # a threshold above that gain: not listed, nothing moves
  out48, desc = run("tst_reduce2_his_states_default_8", asn_str,
    ["non_flip_preference=4.8"])
  assert "Suggested flips (add_flip_movers=False, none made): 0" in desc, desc
  assert_heavy_atoms_unchanged(asn_str, out48)
  for name in ("HD21", "HD22"):
    k = ("A", "151", name)
    assert out48[k].xyz == placed[k].xyz, name
  out, desc = run("tst_reduce2_his_states_default_6", asn_str,
    ["add_flip_movers=True"])
  for name, other in (("ND2", "OD1"), ("OD1", "ND2")):
    assert (matrix.col(out[("A", "151", name)].xyz) -
      matrix.col(start[("A", "151", other)].xyz)).length() < 0.3, name

def test_reduce_flags():
  d = reduce2.Program.description
  for line in (
      "-build: approach=add add_flip_movers=True",
      "-flip: approach=add add_flip_movers=True",
      "-nobuild#: approach=add add_flip_movers=True preference_magnitude=#",
      "-noflip: approach=add add_flip_movers=True preference_magnitude=9999 "
        "(as -nobuild9999)."):
    assert ("    " + line + "\n") in d, line

def test_scoring_leaves_h():
  """Scoring the flips changes no atom: every atom, H included, as in a run
  whose Optimizer scores no flips (the His ring H, deleted or not, too)."""
  from mmtbx.reduce import Optimizers
  for i, s in enumerate((_his_fragment_str(), asn_str)):
    fn = "tst_reduce2_his_states_default_9%d" % i
    run(fn, s)
    with_scores = all_atoms(open(fn + "_out.pdb").read())
    original = Optimizers.Optimizer._SuggestFlips
    Optimizers.Optimizer._SuggestFlips = lambda *args: None
    try:
      run(fn + "n", s)
    finally:
      Optimizers.Optimizer._SuggestFlips = original
    without = all_atoms(open(fn + "n_out.pdb").read())
    assert with_scores == without, (set(with_scores) ^ set(without))

def asn_models_str():
  """asn_str as two MODELs, Asn 151 with insertion code A and its side chain
  in conformers A and B (at the same sites)"""
  body = []
  for l in asn_str.split("\n"):
    if not l.startswith(("ATOM", "HETATM")):
      continue
    if l[17:26] == "ASN A 151":
      l = l[:26] + "A" + l[27:]
      if l[12:16].strip() in ("CB", "CG", "OD1", "ND2"):
        for alt in "AB":
          body.append(l[:16] + alt + l[17:54] + "  0.50" + l[60:])
        continue
    body.append(l)
  lines = []
  for m in (1, 2):
    lines += ["MODEL        %d" % m] + body + ["ENDMDL"]
  return "\n".join(lines + ["END"]) + "\n"

def his_zn_box_str():
  """his_zn_str() with the Zn moved by one edge of its P 1 box: the Zn binds
  only across the box, as the copy at x-1,y,z"""
  out = []
  for l in his_zn_str().split("\n"):
    if l.startswith("HETATM") and l[17:20] == " ZN":
      l = l[:30] + "%8.3f" % (float(l[30:38]) + 40.0) + l[38:]
    out.append(l)
  return "\n".join(out)

def test_label_formats():
  """Model, insertion code and conformer in a Suggested flips line; a metal
  in a symmetry copy with its operator, in the description and the log; a
  metal-only entry whose flip scores worse names the metal as the reason."""
  from libtbx import group_args
  from six.moves import cStringIO as StringIO
  out, desc = run("tst_reduce2_his_states_default_10", asn_models_str())
  listed = suggested(desc)
  for m in (1, 2):
    for alt in "AB":
      w = "model %d chain A ASN 151A conformer %s: flipped scores better by " % (
        m, alt)
      assert [l for l in listed if l.startswith(w)], (w, listed)
  log = StringIO()
  out, desc = run("tst_reduce2_his_states_default_11", his_zn_box_str(),
    log=log)
  listed = suggested(desc)
  assert len(listed) == 1 and listed[0].endswith(
    "; ZN A 101 (symmetry copy x-1,y,z) binds only the flipped ring"), listed
  log = log.getvalue()
  assert ("  chain A HIS 1: ZN A 101 (symmetry copy x-1,y,z) binds only the "
    "flipped ring, CD2 2.10 A from it as modelled, not flipped") in log, log
  # NE2 and CE1 as conformers A and B, CD2 blank: the distance still given
  split = []
  for l in his_zn_str().split("\n"):
    if l.startswith("ATOM") and l[12:16].strip() in ("NE2", "CE1"):
      split += [l[:16] + alt + l[17:54] + "  0.50" + l[60:] for alt in "AB"]
    else:
      split.append(l)
  log = StringIO()
  run("tst_reduce2_his_states_default_18", "\n".join(split), log=log)
  for alt in "AB":
    line = ("  chain A HIS 1 conformer %s: ZN A 101 binds only the flipped ring, "
      "CD2 2.10 A from it as modelled, not flipped" % alt)
    assert line in log.getvalue(), (line, log.getvalue())
  # gain below zero (formatter only: no such case in the fixtures)
  h = iotbx.pdb.input(lines=his_zn_str().split("\n"), source_info=None
    ).construct_hierarchy()
  a = dict((x.name.strip(), x) for x in h.atoms())
  f = group_args(model=1, atom=a["NE2"], kind="HisFlip", gain=-1.5,
    threshold=0.5, chosen="Unflipped HD1Placed HE2NotPlaced",
    flipped="Flipped HD1NotPlaced HE2Placed",
    ion=group_args(atom=a["ZN"], rt_mx=None, ring="CD2", distance=2.1))
  text, metal = reduce2.Program._SuggestedFlipsText(None, [f])
  assert ("  chain A HIS 1: listed for its metal, flipped scores 1.50 worse "
    "(HD1 as chosen, HE2 flipped); ZN A 101 binds only the flipped ring") \
    in text.split("\n"), text
  assert metal == ["  chain A HIS 1: ZN A 101 binds only the flipped ring, "
    "CD2 2.10 A from it as modelled, not flipped"], metal

def test_locked_his_annotation():
  """add_flip_movers=True, His A 1 locked Unflipped (set_flip_states), a
  water 2.4 A out from CE1: its ring H clashes in every protonation state.
  The report line has no BothClash (which compared two of its own
  protonation states as if they were the two orientations)."""
  from scitbx import matrix
  from mmtbx.regression.tst_reduce_switch import _his_atoms
  xyz = dict((n, matrix.col(x)) for n, x in _his_atoms)
  cen = sum([xyz[n] for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
    matrix.col((0, 0, 0))) / 5
  w = xyz["CE1"] + 2.4 * (xyz["CE1"] - cen).normalize()
  fmt = "%-6s%5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %2s"
  lines = ["CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"]
  for i, (n, x) in enumerate(_his_atoms):
    lines.append(fmt % ("ATOM", i + 1, n, "HIS", "A", 1, x[0], x[1], x[2],
      n[0]))
  lines.append(fmt % ("HETATM", 30, "O", "HOH", "A", 301, w[0], w[1], w[2],
    "O"))
  out, desc = run("tst_reduce2_his_states_default_12",
    "\n".join(lines + ["END"]) + "\n",
    ["add_flip_movers=True", "set_flip_states=1 . A HIS 1 Unflipped"])
  pose = [l for l in desc.split("\n") if "HisPlace at chain A HIS 1" in l
    and " pose " in l]
  assert len(pose) == 1, desc
  assert pose[0].endswith(" pose Unflipped HD1Placed HE2NotPlaced . ."), pose

def his_d_zn_str():
  """his_d_str("so4_split") (passes B then A over the blank His D) with a Zn
  as conformer A 2.1 A out from ND1: pass B suggests the flip, pass A locks
  the His unflipped to the Zn and removes HD1."""
  from mmtbx.regression.tst_reduce_switch import his_d_str
  s = his_d_str("so4_split")
  xyz = dict((k[2], matrix.col(a.xyz)) for k, a in atoms(s).items()
    if k[:2] == ("D", "1"))
  cen = sum([xyz[n] for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
    matrix.col((0, 0, 0))) / 5
  zn = xyz["ND1"] + 2.1 * (xyz["ND1"] - cen).normalize()
  lines = [l for l in s.split("\n") if l.strip() and l.strip() != "END"]
  lines.append("HETATM  900 ZN  A ZN D 101    %8.3f%8.3f%8.3f  0.50 10.00"
    "          ZN" % tuple(zn))
  return "\n".join(lines + ["END"]) + "\n"

def test_later_pass_replaces_suggestion():
  """A suggestion of an earlier conformer pass is withdrawn when a later pass
  locks the residue to its metal: none listed, HD1 removed, HE2 kept."""
  out, desc = run("tst_reduce2_his_states_default_15", his_d_zn_str())
  assert ring_h(out, "D") == "HE2", desc
  assert suggested(desc) == [], suggested(desc)
  assert "Suggested flips (add_flip_movers=False, none made): 0" in desc, desc
  assert "flipped scores better by 50.85" not in desc, desc
  assert "Suggested flip of chain D HIS 1 withdrawn" in desc, desc

def _fault_run(prefix, model_str, inject):
  """Run with inject() active around each _SuggestFlips call; assert that
  each call leaves the sites, spatial-query membership, deletions and atom
  info flags of the pass's atoms as before it, and that the output equals a
  run that scores no flips. :return: the description."""
  from mmtbx.reduce import Optimizers
  original = Optimizers.Optimizer._SuggestFlips
  calls = []
  def state(opt):
    st = {}
    for a in opt._atoms:
      if a in opt._environmentAtoms: continue
      found = [n for n in opt._spatialQuery.neighbors(a.xyz, 0.0, 0.001)
        if n.i_seq == a.i_seq]
      e = opt._extraAtomInfo.getMappingFor(a)
      st[a.i_seq] = (tuple(a.xyz), len(found), a in opt._deleteMes,
        e.isAcceptor, e.isDonor, e.isDummyHydrogen, e.isIon, e.vdwRadius,
        e.charge)
    return st
  def checked(opt, *args):
    before = state(opt)
    undo = inject()
    try:
      original(opt, *args)
    finally:
      undo()
    calls.append(len(opt._flipCandidates))
    after = state(opt)
    assert before == after, [k for k in before if before[k] != after[k]]
  Optimizers.Optimizer._SuggestFlips = checked
  try:
    out, desc = run(prefix, model_str)
  finally:
    Optimizers.Optimizer._SuggestFlips = original
  assert sum(calls) > 0, calls
  Optimizers.Optimizer._SuggestFlips = lambda *args: None
  try:
    run(prefix + "n", model_str)
  finally:
    Optimizers.Optimizer._SuggestFlips = original
  got = all_atoms(open(prefix + "_out.pdb").read())
  ref = all_atoms(open(prefix + "n_out.pdb").read())
  assert got == ref, sorted(set(got) ^ set(ref))
  return desc

def test_restore_after_failures():
  """A flip whose scoring fails, after the first flipped state is set, or
  whose amide Mover fails after it re-placed the amide H, leaves the pass's
  atoms and the output as without scoring; a His whose metal binds only the
  flipped ring is listed and named also when its flip is not scored."""
  from mmtbx.reduce import Optimizers, Movers
  def score_fails():
    original = Optimizers.Optimizer._scorePoseAtoms
    n = [0]
    def f(*args):
      n[0] += 1
      if n[0] % 2 == 0:
        raise RuntimeError("injected scoring failure")
      return original(*args)
    Optimizers.Optimizer._scorePoseAtoms = f
    def undo():
      Optimizers.Optimizer._scorePoseAtoms = original
    return undo
  def amide_fails():
    original = Movers.MoverAmideFlip
    def f(*args, **kw):
      original(*args, **kw)
      raise RuntimeError("injected constructor failure")
    Movers.MoverAmideFlip = f
    def undo():
      Movers.MoverAmideFlip = original
    return undo
  for i, s in enumerate((_his_fragment_str(), asn_str, his_zn_str())):
    desc = _fault_run("tst_reduce2_his_states_default_2%d" % i, s, score_fails)
    assert "Suggested flips (add_flip_movers=False, none made): " in desc
  desc = _fault_run("tst_reduce2_his_states_default_23", asn_str, amide_fails)
  assert "Suggested flips (add_flip_movers=False, none made): 0" in desc, desc
  # his_zn_str: the second score (the flipped state) fails
  desc = open("tst_reduce2_his_states_default_22.txt").read()
  assert suggested(desc) == ["chain A HIS 1: listed for its metal, flip not "
    "scored; ZN A 101 binds only the flipped ring"], suggested(desc)

def test_reduce_flags_behaviour():
  """-noflip and -nobuild9999 (add_flip_movers=True preference_magnitude=9999)
  make no scored flip, but the 9999 is a bias: the ion lock still flips the His
  whose Zn binds only the flipped ring, which the default does not. It also
  scales the bias against a His with both ring H: His C, which gets both by
  default, keeps one."""
  s = _his_fragment_str()
  default, _ = run("tst_reduce2_his_states_default_13", s)
  noflip, _ = run("tst_reduce2_his_states_default_14", s,
    ["add_flip_movers=True", "preference_magnitude=9999"])
  assert_heavy_atoms_unchanged(s, noflip)
  assert [ring_h(default, c) for c in "ABC"] == ["HD1", "HE2", "HD1+HE2"]
  assert [ring_h(noflip, c) for c in "ABC"] == ["HD1", "HE2", "HD1"]
  s = his_zn_str()
  start = atoms(s)
  for i, args in enumerate(([], ["add_flip_movers=True",
      "preference_magnitude=9999"])):
    out, _ = run("tst_reduce2_his_states_default_1%d" % (6 + i), s, args)
    d = (matrix.col(out[("A", "1", "ND1")].xyz) -
      matrix.col(start[("A", "1", "CD2")].xyz)).length()
    assert (d < 0.3) == (i == 1), (args, d)

if __name__ == "__main__":
  test_default_is_false()
  test_his_states_by_default()
  test_metal_his_listed_not_flipped()
  test_asn_listed_not_flipped()
  test_reduce_flags()
  test_scoring_leaves_h()
  test_label_formats()
  test_locked_his_annotation()
  test_reduce_flags_behaviour()
  test_later_pass_replaces_suggestion()
  test_restore_after_failures()
  print("OK")
