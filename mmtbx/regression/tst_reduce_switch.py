"""Tests for the reduce1/reduce2 switch in mmtbx.hydrogens (Stage 1 infrastructure).

Covers the in-process reduce2 hydrogen engine (place_and_optimize_hydrogens), the
central add_hydrogens() dispatcher and its USE_OLD_REDUCE switch variable, the legacy
reduce1 path, get_nqh_flips_reduce2() (decision-only NQH flip extraction), and the
clashscore2 delegation refactor.
See docs/superpowers/plans/2026-06-03-reduce2-switch.md.
"""
from __future__ import absolute_import, division, print_function
from mmtbx import hydrogens as reduce_switch
from iotbx.data_manager import DataManager
from libtbx.utils import null_out
from mmtbx.rotamer import nqh
from libtbx import group_args

# Small 3-residue peptide, no hydrogens.
pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       8.000   8.000   8.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       9.000   8.000   8.000  1.00 10.00           C
ATOM      3  C   GLY A   1       9.500   9.400   8.000  1.00 10.00           C
ATOM      4  O   GLY A   1       8.800  10.400   8.000  1.00 10.00           O
ATOM      5  N   SER A   2      10.800   9.500   8.000  1.00 10.00           N
ATOM      6  CA  SER A   2      11.500  10.800   8.000  1.00 10.00           C
ATOM      7  C   SER A   2      13.000  10.600   8.000  1.00 10.00           C
ATOM      8  O   SER A   2      13.500   9.500   8.000  1.00 10.00           O
ATOM      9  CB  SER A   2      11.200  11.700   9.200  1.00 10.00           C
ATOM     10  OG  SER A   2      11.800  12.950   9.000  1.00 10.00           O
ATOM     11  N   ALA A   3      13.700  11.700   8.000  1.00 10.00           N
ATOM     12  CA  ALA A   3      15.150  11.700   8.000  1.00 10.00           C
ATOM     13  C   ALA A   3      15.700  13.100   8.000  1.00 10.00           C
ATOM     14  O   ALA A   3      15.000  14.100   8.000  1.00 10.00           O
ATOM     15  CB  ALA A   3      15.650  10.950   9.230  1.00 10.00           C
TER
END
"""

def _h_count(model):
  es = model.get_hierarchy().atoms().extract_element()
  return sum(1 for e in es if e.strip() in ("H", "D"))

def model_from_str(s):
  dm = DataManager(['model'])
  dm.process_model_str("test.pdb", s)
  m = dm.get_model()
  m.add_crystal_symmetry_if_necessary()
  return m

def test_adds_hydrogens():
  """Engine place_and_optimize_hydrogens() adds hydrogens to an H-less model."""
  model = model_from_str(pdb_str)
  assert _h_count(model) == 0
  out = reduce_switch.place_and_optimize_hydrogens(
    model=model, do_flips=False, nuclear=False, log=null_out())
  assert _h_count(out) > 0, _h_count(out)
  print("test_adds_hydrogens OK")

def test_dispatcher_explicit_reduce2():
  """add_hydrogens(old=False) routes to reduce2 and returns an H-bearing model."""
  out = reduce_switch.add_hydrogens(
    model=model_from_str(pdb_str), old=False, log=null_out())
  assert _h_count(out) > 0
  print("test_dispatcher_explicit_reduce2 OK")

def test_dispatcher_default_uses_switch():
  """add_hydrogens(old=None) honors the USE_OLD_REDUCE switch variable (here -> reduce2)."""
  # old=None must read mmtbx.hydrogens.USE_OLD_REDUCE (the single switch).
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    out = reduce_switch.add_hydrogens(
      model=model_from_str(pdb_str), old=None, log=null_out())
    assert _h_count(out) > 0          # routed to reduce2
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_dispatcher_default_uses_switch OK")

def test_dispatcher_reduce1():
  """add_hydrogens(old=True) routes to the legacy molprobity.reduce path (skipped if absent)."""
  import libtbx.load_env
  if not libtbx.env.has_module(name="reduce"):
    print("reduce module missing, skipping reduce1 path test")
    return
  out = reduce_switch.add_hydrogens(
    model=model_from_str(pdb_str), old=True, log=null_out())
  assert _h_count(out) > 0
  print("test_dispatcher_reduce1 OK")

def _h_restrained_model(model):
  # Produce an H-bearing, restraints-carrying model (what get_nqh_flips_reduce2 expects).
  from mmtbx.hydrogens import reduce_hydrogen
  o = reduce_hydrogen.place_hydrogens(model=model, keep_existing_H=False); o.run()
  m = o.get_model()
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(False)
  p.pdb_interpretation.flip_symmetric_amino_acids = False
  m.get_hierarchy().sort_atoms_in_place(); m.get_hierarchy().atoms().reset_serial()
  m.process(make_restraints=True, pdb_interpretation_params=p)
  return m

def test_get_nqh_flips_requires_h():
  """get_nqh_flips_reduce2() refuses an H-less model -- it never silently adds H."""
  # Design rule: the function must NOT add H; it requires an H-bearing model.
  m = model_from_str(pdb_str)            # no hydrogens
  try:
    reduce_switch.get_nqh_flips_reduce2(m)
  except AssertionError:
    print("test_get_nqh_flips_requires_h OK")
    return
  raise RuntimeError("expected AssertionError for an H-less model")

def test_get_nqh_flips_detects_flip():
  """get_nqh_flips_reduce2() flags a backwards ASN, in flip_selected() key format,
  without mutating the input model (decision-only, atom-preserving)."""
  import os, libtbx.load_env
  p = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/ions/3zu8.pdb", test=os.path.isfile)
  if p is None:
    print("3zu8.pdb not found, skipping flip-detection test"); return
  dm = DataManager(['model']); dm.process_model_file(p)
  m = dm.get_model(); m.add_crystal_symmetry_if_necessary()
  # Force a backwards ASN by swapping OD1/ND2 so reduce2 wants to flip it back.
  h = m.get_hierarchy()
  for rg in h.residue_groups():
    for ag in rg.atom_groups():
      if ag.resname.strip() == "ASN":
        ad = {a.name.strip(): a for a in ag.atoms()}
        if "OD1" in ad and "ND2" in ad:
          ad["OD1"].xyz, ad["ND2"].xyz = ad["ND2"].xyz, ad["OD1"].xyz
  m.set_sites_cart(h.atoms().extract_xyz())
  hm = _h_restrained_model(m)
  n_before = hm.get_number_of_atoms()
  user_mods, atom_notes, score_dict = reduce_switch.get_nqh_flips_reduce2(hm)
  assert hm.get_number_of_atoms() == n_before, "input model must not be mutated"
  assert any("ASN" in k for k in user_mods), user_mods
  assert atom_notes == [] and score_dict == {}
  print("test_get_nqh_flips_detects_flip OK:", user_mods)

def test_clashscore2_delegates():
  """clashscore2.check_and_add_hydrogen still adds H / returns changed=True after
  being refactored to delegate to place_and_optimize_hydrogens (no behavior change)."""
  # Pins clashscore2.check_and_add_hydrogen behavior: an H-less model gets H added
  # and changed=True. Must hold both before and after the delegation refactor.
  from mmtbx.validation import clashscore2
  probe_phil = reduce_switch.default_probe_phil()
  out_model, changed = clashscore2.check_and_add_hydrogen(
    probe_parameters=probe_phil, data_manager_model=model_from_str(pdb_str),
    nuclear=False, keep_hydrogens=True, do_flips=False, log=null_out())
  assert changed is True
  assert _h_count(out_model) > 0
  print("test_clashscore2_delegates OK")

# Same peptide plus an unknown ligand (UNL). mon_lib_query returns None for this
# residue name, so reduce2 records it in no_H_placed_mlq -- which is what the
# missing-restraints handling keys off. (Many made-up 3-letter codes resolve to
# some monomer entry instead and would not exercise this path at all.)
pdb_str_unknown_ligand = pdb_str.replace("TER\nEND\n", "") + """\
HETATM   16  C1  UNL B   1      15.000   5.000   5.000  1.00 20.00           C
HETATM   17  C2  UNL B   1      16.400   5.000   5.000  1.00 20.00           C
HETATM   18  O1  UNL B   1      17.100   6.100   5.000  1.00 20.00           O
TER
END
"""

def test_clashscore2_ignore_missing_restraints():
  """check_and_add_hydrogen aborts on a restraint-less ligand by default, and
  continues best-effort when ignore_missing_restraints=True.

  Pins the plumbing of ignore_missing_restraints through to
  place_and_optimize_hydrogens(raise_on_missing=...). The 2-tuple return shape is
  unchanged on both paths.
  """
  from mmtbx.validation import clashscore2
  from libtbx.utils import Sorry
  probe_phil = reduce_switch.default_probe_phil()

  # Default: the ligand with no restraints aborts the whole run.
  try:
    clashscore2.check_and_add_hydrogen(
      probe_parameters=probe_phil,
      data_manager_model=model_from_str(pdb_str_unknown_ligand),
      nuclear=False, keep_hydrogens=True, do_flips=False, log=null_out())
  except Sorry as e:
    assert "Restraints were not found" in str(e), str(e)
  else:
    raise AssertionError("expected Sorry for a ligand with no restraints when "
                         "ignore_missing_restraints=False")

  # Best-effort: completes, still returns the 2-tuple, and the rest of the model
  # still gets hydrogens.
  out_model, changed = clashscore2.check_and_add_hydrogen(
    probe_parameters=probe_phil,
    data_manager_model=model_from_str(pdb_str_unknown_ligand),
    nuclear=False, keep_hydrogens=True, do_flips=False,
    ignore_missing_restraints=True, log=null_out())
  assert changed is True
  assert _h_count(out_model) > 0
  print("test_clashscore2_ignore_missing_restraints OK")

def _make_all_asn_backwards(hierarchy):
  """Swap OD1/ND2 on every ASN so reduce wants to flip them back. Returns count."""
  n = 0
  for rg in hierarchy.residue_groups():
    for ag in rg.atom_groups():
      if ag.resname.strip() == "ASN":
        ad = {a.name.strip(): a for a in ag.atoms()}
        if "OD1" in ad and "ND2" in ad:
          ad["OD1"].xyz, ad["ND2"].xyz = ad["ND2"].xyz, ad["OD1"].xyz
          n += 1
  return n

def test_flip_nqh_skips_h_less_model_reduce2():
  """Under reduce2 no H is added mid-refinement, so flip_nqh skips an H-less
  model: it logs 'skipping', applies no flip, and leaves the model unchanged."""
  from io import StringIO
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    m = model_from_str(pdb_str)                  # no hydrogens
    assert not m.has_hd()
    log = StringIO(); m.set_log(log)
    n_before = m.get_number_of_atoms()
    nqh.flip(model=m, log=m.log)                 # skips (no raise)
    assert m.get_number_of_atoms() == n_before
    out = log.getvalue().lower()
    assert "skipping" in out and "hydrogen" in out, log.getvalue()
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_flip_nqh_skips_h_less_model_reduce2 OK")

def test_flip_nqh_applies_reduce2_flip():
  """Under reduce2, an H-bearing model with backwards ASN side-chains is flipped
  via the reduce2 decision path, atom-preservingly (no atom added or removed).
  Also exercises the full path get_nqh_flips_reduce2 -> flip_selected, so a
  flip-key format mismatch would surface as 0 flips applied."""
  import os, re, libtbx.load_env
  from io import StringIO
  p = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/ions/3zu8.pdb", test=os.path.isfile)
  if p is None:
    print("3zu8.pdb not found, skipping reduce2 flip-apply test"); return
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    dm = DataManager(['model']); dm.process_model_file(p)
    m = dm.get_model(); m.add_crystal_symmetry_if_necessary()
    assert _make_all_asn_backwards(m.get_hierarchy()) > 0
    m.set_sites_cart(m.get_hierarchy().atoms().extract_xyz())
    hm = _h_restrained_model(m)                   # add H + restraints
    assert hm.get_hd_selection().count(True) > 0
    n_before = hm.get_number_of_atoms()
    log = StringIO(); hm.set_log(log)
    nqh.flip(model=hm, log=hm.log)                # reduce2 decision-only flip
    out = log.getvalue()
    assert hm.get_number_of_atoms() == n_before   # atom-preserving
    mobj = re.search(r"Total number of N/Q/H flips:\s*(\d+)", out)
    assert mobj and int(mobj.group(1)) >= 1, out  # a flip was actually applied
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_flip_nqh_applies_reduce2_flip OK")

def test_flip_nqh_use_ncs_without_groups_reduce2():
  """use_ncs=True on a model with no NCS groups falls back to a whole-model flip
  (get_master_selection is not consulted) and still flips a backwards ASN, atom-
  preservingly. The NCS master-flip + multiply_ncs propagation itself is covered
  end-to-end by real_space_refine tst_73."""
  import os, re, libtbx.load_env
  from io import StringIO
  p = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/ions/3zu8.pdb", test=os.path.isfile)
  if p is None:
    print("3zu8.pdb not found, skipping use_ncs test"); return
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    dm = DataManager(['model']); dm.process_model_file(p)
    m = dm.get_model(); m.add_crystal_symmetry_if_necessary()
    assert _make_all_asn_backwards(m.get_hierarchy()) > 0
    m.set_sites_cart(m.get_hierarchy().atoms().extract_xyz())
    hm = _h_restrained_model(m)
    assert not hm.get_ncs_groups()                 # no NCS -> whole-model fallback
    n_before = hm.get_number_of_atoms()
    log = StringIO(); hm.set_log(log)
    nqh.flip(model=hm, use_ncs=True, log=hm.log)    # safe: falls back to whole model
    assert hm.get_number_of_atoms() == n_before     # atom-preserving
    mo = re.search(r"Total number of N/Q/H flips:\s*(\d+)", log.getvalue())
    assert mo and int(mo.group(1)) >= 1, log.getvalue()
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_flip_nqh_use_ncs_without_groups_reduce2 OK")

def test_clash_none_when_no_h_reduce2():
  """Under reduce2 (USE_OLD_REDUCE=False), geometry clash() is a NON-ADDING
  observer: an H-less model yields score=None (no reduce1 force-add), and
  mp_score()/show() tolerate the None instead of raising."""
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    model = model_from_str(pdb_str)                    # no hydrogens
    model.process(make_restraints=True)
    gs = model.geometry_statistics(use_hydrogens=True) # keep H (there are none)
    c = gs.clash()
    assert c.score is None, c.score
    assert c.clashes is None, c.clashes
    assert gs.mp_score() is None, gs.mp_score()
    gs.show(log=null_out())                            # must not raise on None
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_clash_none_when_no_h_reduce2 OK")

def test_clash_scores_existing_h_reduce2():
  """Under reduce2, an H-bearing model still gets a real (non-None) clashscore:
  Probe runs on the existing H and no reduce is invoked (regression guard for the
  H-present branch of the non-adding observer)."""
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    hmodel = reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(pdb_str), do_flips=False, nuclear=False, log=null_out())
    assert _h_count(hmodel) > 0
    hmodel.process(make_restraints=True)               # engine left no restraints
    gs = hmodel.geometry_statistics(use_hydrogens=True)
    score = gs.clash().score
    assert score is not None and score >= 0, score     # a real clashscore
    assert gs.mp_score() is not None
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_clash_scores_existing_h_reduce2 OK")

def test_rsr_per_cycle_clashscore_keeps_h_reduce2():
  """Under reduce2 the real-space per-cycle stats (rsr_model.update_statistics)
  must report a clashscore for an H-bearing model. The bug: geometry_statistics()
  with default use_hydrogens=None strips H (electron table / riding H), so reduce2
  clash() returns None. The fix passes use_hydrogens=True under reduce2."""
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    hmodel = reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(pdb_str), log=null_out())
    hmodel.process(make_restraints=True)
    hmodel.setup_scattering_dictionaries(scattering_table="electron")
    assert hmodel.has_hd()
    xrs = hmodel.get_xray_structure()
    fft_map = xrs.structure_factors(d_min=2.0).f_calc().fft_map(
      resolution_factor=1./3)
    fft_map.apply_sigma_scaling()
    from mmtbx.refinement.real_space import rsr_model
    rm = rsr_model(model=hmodel, map_data=fft_map.real_map_unpadded(), d_min=2.0)
    rm.update_statistics()
    geo = rm.stats_evaluations[-1].geometry
    assert geo is not None and geo.clash.score is not None, \
      "per-cycle clashscore is None (H stripped before clash)"
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_rsr_per_cycle_clashscore_keeps_h_reduce2 OK")

def test_molprobity_validation_adds_h_reduce2():
  """Standalone comprehensive validation (mmtbx.validation.molprobity.molprobity)
  needs a real clashscore. Under reduce2 clash() is a non-adding observer, so for
  an H-less model clash.clashes is None (-> AttributeError on
  probe_clashscore_manager). Per policy, molprobity adds H itself -- gated on
  use_old_reduce(), on a copy -- so clash() then scores the (added) H. No-op under
  reduce1."""
  import mmtbx.validation.molprobity
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    m = model_from_str(pdb_str)                       # no hydrogens
    assert not m.has_hd()
    v = mmtbx.validation.molprobity.molprobity(model=m, keep_hydrogens=False)
    assert v.clashes is not None                      # H added -> real clash object
    assert v.clashscore() is not None                 # a real clashscore, not None
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_molprobity_validation_adds_h_reduce2 OK")

def test_reduce2_flips_applies_wholesale():
  """USE_REDUCE2_FLIPS on: flip() applies reduce2's flips AND optimized H wholesale
  via the Optimizer (one run), atom-preserving even though reduce2 would otherwise
  delete H (3zu8 has a His near Ni). A flip is actually performed and H are kept."""
  import os, re, libtbx.load_env
  from io import StringIO
  from mmtbx.reduce import Optimizers
  p = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/ions/3zu8.pdb", test=os.path.isfile)
  if p is None:
    print("3zu8.pdb not found, skipping wholesale flip test"); return
  saved_old   = reduce_switch.USE_OLD_REDUCE
  saved_flips = getattr(reduce_switch, "USE_REDUCE2_FLIPS", False)
  try:
    reduce_switch.USE_OLD_REDUCE = False
    dm = DataManager(['model']); dm.process_model_file(p)
    m = dm.get_model(); m.add_crystal_symmetry_if_necessary()
    assert _make_all_asn_backwards(m.get_hierarchy()) > 0
    m.set_sites_cart(m.get_hierarchy().atoms().extract_xyz())
    hm = _h_restrained_model(m)
    n_before = hm.get_number_of_atoms()
    hsel = hm.get_hd_selection()
    # reduce2 WOULD delete H here and wants flips -> atom-preservation is non-trivial
    chk = hm.deep_copy()
    opt = Optimizers.Optimizer(reduce_switch.default_probe_phil(), True, chk,
                               modelIndex=None, altID="")
    assert len(opt.getHydrogensToDelete()) > 0           # reduce2 wants to delete H
    assert len(opt.getFlippedAmides()) + len(opt.getFlippedHistidines()) >= 1
    # switch 1 ON: wholesale flips + H, atom-preserving (deletions NOT applied)
    reduce_switch.USE_REDUCE2_FLIPS = True
    sites_before = hm.get_sites_cart().deep_copy()
    log = StringIO()
    nqh.flip(model=hm, log=log)
    assert hm.get_number_of_atoms() == n_before          # atom-preserving
    assert hm.has_hd()
    mo = re.search(r"Total number of N/Q/H flips:\s*(\d+)", log.getvalue())
    assert mo and int(mo.group(1)) >= 1, log.getvalue()  # a flip was applied
    moved = (hm.get_sites_cart() - sites_before).norms() > 0.1
    assert (moved & ~hsel).count(True) > 0               # heavy atoms flipped
  finally:
    reduce_switch.USE_OLD_REDUCE   = saved_old
    reduce_switch.USE_REDUCE2_FLIPS = saved_flips
  print("test_reduce2_flips_applies_wholesale OK")

# His36 of 1a6m (heavy atoms only).
_his_atoms = [
  ("N",   (18.824, 16.000, -3.674)), ("CA",  (17.931, 16.835, -4.470)),
  ("C",   (18.616, 18.103, -4.877)), ("O",   (18.445, 19.170, -4.297)),
  ("CB",  (16.629, 17.148, -3.681)), ("CG",  (15.962, 15.936, -3.173)),
  ("ND1", (15.130, 15.122, -3.916)), ("CD2", (16.038, 15.336, -1.936)),
  ("CE1", (14.734, 14.111, -3.162)), ("NE2", (15.270, 14.220, -1.964))]

def _his_fragment_str():
  """Four isolated His (chains A-D, 25 A apart), each with SO4 acceptor O 2.9 A
  out from the listed ring atoms. A: ND1, B: NE2, C: ND1+NE2, D: CD2+CE1 (only a
  ring flip gives those two N-H partners)."""
  import math
  from scitbx import matrix
  targets = [("ND1",), ("NE2",), ("ND1", "NE2"), ("CD2", "CE1")]
  xyz = dict((n, matrix.col(x)) for n, x in _his_atoms)
  cen = sum([xyz[n] for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
            matrix.col((0, 0, 0))) / 5
  c, s = math.cos(math.radians(109.47)), math.sin(math.radians(109.47))
  fmt = "%-6s%5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %2s"
  lines = ["CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1"]
  i_seq = 1
  for i, chain in enumerate("ABCD"):
    shift = matrix.col((25.0 * i, 0, 0))
    for n, x in _his_atoms:
      p = matrix.col(x) + shift
      lines.append(fmt % ("ATOM", i_seq, n, "HIS", chain, 1, p[0], p[1], p[2], n[0]))
      i_seq += 1
    for j, t in enumerate(targets[i]):
      u = (xyz[t] - cen).normalize()
      o1 = xyz[t] + 2.9 * u + shift
      sp = o1 + 1.47 * u
      p1 = u.ortho().normalize(); p2 = u.cross(p1)
      atoms = [("S", sp), ("O1", o1)] + [
        ("O%d" % (k + 2), sp + 1.47 * (-u * c + (p1 * math.cos(a) + p2 * math.sin(a)) * s))
        for k, a in enumerate((0, 2 * math.pi / 3, 4 * math.pi / 3))]
      for n, p in atoms:
        lines.append(fmt % ("HETATM", i_seq, n, "SO4", chain, 101 + j,
                            p[0], p[1], p[2], n[0]))
        i_seq += 1
  return "\n".join(lines + ["END"]) + "\n"

def _ring_h_and_heavy_shift(model, ref_model):
  """({chain: 'HD1'|'HE2'|'HD1+HE2'|''}, max heavy-atom shift vs ref_model)."""
  def key(a):
    return (a.chain().id, a.parent().resname, a.parent().parent().resseq, a.name)
  ref = dict((key(a), a.xyz) for a in ref_model.get_hierarchy().atoms())
  ring_h = {}
  shift = 0.0
  for a in model.get_hierarchy().atoms():
    if a.element_is_hydrogen():
      if a.name.strip() in ("HD1", "HE2"):
        ring_h.setdefault(a.chain().id, []).append(a.name.strip())
      continue
    d = sum((x - y) ** 2 for x, y in zip(a.xyz, ref[key(a)])) ** 0.5
    shift = max(shift, d)
  return dict((k, "+".join(sorted(v))) for k, v in ring_h.items()), shift

def test_his_protonation_without_flips():
  """optimize_his_protonation=True with do_flips=False: reduce2 picks each His
  ring protonation (both tautomers and HD1+HE2 reachable) and moves no heavy atom,
  even where a ring flip scores better. Default keeps H on both ring N."""
  s = _his_fragment_str()
  default = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(s), do_flips=False, log=null_out())
  ring_h, shift = _ring_h_and_heavy_shift(default, model_from_str(s))
  assert ring_h == {"A": "HD1+HE2", "B": "HD1+HE2", "C": "HD1+HE2",
                    "D": "HD1+HE2"}, ring_h
  assert shift < 1.e-3, shift
  flipped = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(s), do_flips=True, log=null_out())
  ring_h, shift = _ring_h_and_heavy_shift(flipped, model_from_str(s))
  assert shift > 1.0, shift                        # the fragment favours a flip of D
  his = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(s), do_flips=False, optimize_his_protonation=True,
    log=null_out())
  ring_h, shift = _ring_h_and_heavy_shift(his, model_from_str(s))
  assert ring_h == {"A": "HD1", "B": "HE2", "C": "HD1+HE2", "D": "HD1"}, ring_h
  assert shift < 1.e-3, shift                      # no ring or amide flipped
  print("test_his_protonation_without_flips OK")

def _report():
  return group_args(flips=[], metal_amines_nh2=[], symmetry_copies=[])

def _hg_to_asp(model):
  from mmtbx.regression.tst_reduce2_symmetry_copies import hg_to_asp
  return hg_to_asp(model.model_as_pdb())

def test_symmetry_copies_option():
  """crystal_symmetry: the Optimizer sees the symmetry copies. Ser A 1 OG is
  2.65 A from a carboxylate O of the Asp B 1 mate (P 1 2 1): its HG turns to
  that O with copies and not without (None, the default; or the CRYST1 1 1 1
  placeholder). One report line per model, also in the multi-model path."""
  from mmtbx.regression import tst_reduce2_symmetry_copies as tsc
  ser_asp_str = tsc.ser_asp_str
  default = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(ser_asp_str), log=null_out())
  assert _hg_to_asp(default) > 2.5, _hg_to_asp(default)
  m = model_from_str(ser_asp_str)
  report = _report()
  on = reduce_switch.place_and_optimize_hydrogens(model=m, log=null_out(),
    crystal_symmetry=m.crystal_symmetry(), report=report)
  assert _hg_to_asp(on) < 2.0, _hg_to_asp(on)
  assert len(report.symmetry_copies) == 1, report.symmetry_copies
  assert report.symmetry_copies[0].startswith(
    "Symmetry copies (fixed environment): model 1, 29 copy atoms"), \
    report.symmetry_copies
  # the placeholder cell: no copies, said so
  atoms = "".join(l + "\n" for l in ser_asp_str.splitlines()
    if l.startswith("ATOM"))
  text = tsc.cryst1_placeholder + atoms + "END\n"
  # the input's own symmetry: model_from_str() boxes this model
  import iotbx.pdb
  cs = iotbx.pdb.input(source_info=None, lines=text).crystal_symmetry()
  report = _report()
  reduce_switch.place_and_optimize_hydrogens(model=model_from_str(text),
    log=null_out(), crystal_symmetry=cs, report=report)
  assert report.symmetry_copies == [
    "none: no unit cell, or the CRYST1 1 1 1 placeholder"], report
  # two models: copies per model
  lines = ser_asp_str.splitlines()
  body = [l for l in lines if l.startswith("ATOM")]
  two = "\n".join([lines[0], "MODEL        1"] + body + ["ENDMDL",
    "MODEL        2"] + body + ["ENDMDL", "END"]) + "\n"
  m = model_from_str(two)
  report = _report()
  on = reduce_switch.place_and_optimize_hydrogens(model=m, log=null_out(),
    crystal_symmetry=m.crystal_symmetry(), report=report)
  # each MODEL is named by its id, not the Optimizer's own count
  assert len(report.symmetry_copies) == 2, report.symmetry_copies
  for i, line in enumerate(report.symmetry_copies):
    assert line.startswith("Symmetry copies (fixed environment): model %d, "
      "29 copy atoms" % (i + 1)), report.symmetry_copies
  for i in (0, 1):
    assert tsc.hg_to_asp(on.model_as_pdb(), i) < 2.0, i
  print("test_symmetry_copies_option OK")

def test_metal_amine_option():
  """metal_bound_amines_as_nh2: 1xmk Lys A 359 NZ on Cd A 396 (2.21 A) gets
  two H instead of three, and is reported."""
  from mmtbx.regression.reduce2_metal_amine_utils import pdb_lys_cd
  def nz_h(model):
    return sorted(a.name.strip() for a in model.get_hierarchy().atoms()
      if a.name.strip().startswith("HZ"))
  default = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(pdb_lys_cd), log=null_out())
  assert nz_h(default) == ["HZ1", "HZ2", "HZ3"], nz_h(default)
  report = _report()
  nh2 = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(pdb_lys_cd), log=null_out(),
    metal_bound_amines_as_nh2=True, report=report)
  assert len(nz_h(nh2)) == 2, nz_h(nh2)
  assert len(report.metal_amines_nh2) == 1, report.metal_amines_nh2
  text, contact = report.metal_amines_nh2[0]
  assert text.startswith("NZ  LYS A 359: "), text
  assert "CD    CD A 396" in contact, contact
  print("test_metal_amine_option OK")

def test_new_options_default_off():
  """The two new arguments default to off, and a report changes nothing:
  the same model as without them. The report lists flips (do_flips)."""
  from mmtbx.regression.reduce2_metal_amine_utils import pdb_lys_cd
  from mmtbx.regression.tst_reduce2_symmetry_copies import ser_asp_str
  for s in (pdb_lys_cd, ser_asp_str, _his_fragment_str()):
    plain = reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(s), do_flips=True, log=null_out())
    report = _report()
    explicit = reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(s), do_flips=True, log=null_out(),
      metal_bound_amines_as_nh2=False, crystal_symmetry=None, report=report)
    assert plain.model_as_pdb() == explicit.model_as_pdb()
    assert report.metal_amines_nh2 == [], report
    assert report.symmetry_copies == ["not requested"], report
  # the His fragment favours a ring flip of His D (test above)
  assert [f for f in report.flips if f[1] == "D" and f[6] == "His"], \
    report.flips
  print("test_new_options_default_off OK")

def his_d_str(mode):
  """His D of _his_fragment_str() and its SO4 (a ring flip gives both N-H an
  acceptor). mode 'blank': as it is; 'so4_split': the SO4 as conformer B at
  its site and conformer A 30 A away (passes B then A: the last decides the
  blank His); 'his_split': the His as conformers A and B, the SO4 blank."""
  out = []
  for l in _his_fragment_str().splitlines():
    if l.startswith(("ATOM", "HETATM")) and l[21] != "D": continue
    if not l.startswith(("ATOM", "HETATM")):
      out.append(l); continue
    het = l.startswith("HETATM")
    if (mode == "so4_split" and het) or (mode == "his_split" and not het):
      a = l[:16] + "A" + l[17:54] + "  0.50" + l[60:]
      if het:
        a = a[:30] + "%8.3f" % (float(l[30:38]) + 30.0) + a[38:]
      out += [a, l[:16] + "B" + l[17:54] + "  0.50" + l[60:]]
    else:
      out.append(l)
  return "\n".join(out) + "\n"

def his_zn_str():
  """One His with a Zn 2.1 A out from CD2: only the flipped ring puts an N
  on the Zn (the Optimizer's ion lock-down); an ethanol 4 A out from CE1."""
  from scitbx import matrix
  xyz = dict((n, matrix.col(x)) for n, x in _his_atoms)
  cen = sum([xyz[n] for n in ("CG", "ND1", "CD2", "CE1", "NE2")],
    matrix.col((0, 0, 0))) / 5
  zn = xyz["CD2"] + 2.1 * (xyz["CD2"] - cen).normalize()
  fmt = "%-6s%5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %2s"
  lines = ["CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"]
  for i, (n, x) in enumerate(_his_atoms):
    lines.append(fmt % ("ATOM", i + 1, n, "HIS", "A", 1, x[0], x[1], x[2],
      n[0]))
  lines.append(fmt % ("HETATM", 20, "ZN", " ZN", "A", 101, zn[0], zn[1], zn[2],
    "ZN"))
  e = xyz["CE1"] + 4.0 * (xyz["CE1"] - cen).normalize()
  for k, (n, d) in enumerate((("C1", (0, 0, 0)), ("C2", (1.5, 0, 0)),
      ("O", (-0.5, 1.3, 0)))):
    p = e + matrix.col(d)
    lines.append(fmt % ("HETATM", 21 + k, n, "EOH", "A", 201, p[0], p[1], p[2],
      n[0]))
  return "\n".join(lines + ["END"]) + "\n"

def his_cu_box_str():
  """7T5E His A 1 on Cu A 301 with its glycan, the Cu moved by one edge of a
  P 1 box (CRYST1): the Cu binds the amine only across that box."""
  import iotbx.pdb
  from mmtbx.regression.reduce2_metal_amine_utils import pdb_his_cu_glycan
  h = iotbx.pdb.input(lines=pdb_his_cu_glycan, source_info=None
    ).construct_hierarchy()
  xs = [a.xyz[0] for a in h.atoms()]
  edge = max(xs) - min(xs) + 12.0
  for a in h.atoms():
    if a.parent().resname.strip() == "CU":
      a.xyz = (a.xyz[0] + edge, a.xyz[1], a.xyz[2])
  from cctbx import crystal
  cs = crystal.symmetry((edge, edge, edge, 90, 90, 90), "P 1")
  return h.as_pdb_string(crystal_symmetry=cs)

def test_flips_from_final_model():
  """report.flips are the residues whose final heavy atoms are flipped: a
  blank His flipped only in the first conformer pass is not listed; each
  flipped conformer once; with do_flips=False nothing, except a His whose Zn
  binds only the flipped ring when optimize_his_protonation is False."""
  def flips(s, **kw):
    report = _report()
    reduce_switch.place_and_optimize_hydrogens(model=model_from_str(s),
      log=null_out(), report=report, **kw)
    return report.flips
  key = ("", "D", "   1", " ", "HIS")
  assert flips(his_d_str("blank"), do_flips=True) == [key + ("", "His")]
  assert flips(his_d_str("so4_split"), do_flips=True) == []
  assert flips(his_d_str("his_split"), do_flips=True) == [
    key + ("A", "His"), key + ("B", "His")]
  his = [("", "A", "   1", " ", "HIS", "", "His")]
  assert flips(his_zn_str(), do_flips=True, optimize_his_protonation=True) \
    == his
  assert flips(his_zn_str(), do_flips=False) == his
  assert flips(his_zn_str(), do_flips=False,
    optimize_his_protonation=True) == []
  print("test_flips_from_final_model OK")

def test_ion_lock_without_flips():
  """His A 1 of his_zn_str(), whose Zn binds only the flipped ring, without
  do_flips: the ion lock-down flips it (a ring N on the Zn, without its H)
  unless optimize_his_protonation, which moves no heavy atom (the ring CD2
  and its H stay on the Zn side)."""
  from scitbx import matrix
  def near_zn(**kw):
    m = reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(his_zn_str()), log=null_out(), do_flips=False, **kw)
    atoms = dict((a.name.strip(), a) for a in m.get_hierarchy().atoms()
      if a.parent().resname.strip() in ("HIS", "ZN"))
    zn = matrix.col(atoms["ZN"].xyz)
    return dict((n, (matrix.col(atoms[n].xyz) - zn).length())
      for n in ("ND1", "CD2", "NE2", "CE1", "HD1", "HD2", "HE2") if n in atoms)
  d = near_zn()
  assert d["ND1"] < 2.3 and "HD1" not in d, d
  def key(a):
    return (a.parent().resname, a.name.strip())
  start = dict((key(a), a.xyz) for a in
    model_from_str(his_zn_str()).get_hierarchy().atoms())
  m = reduce_switch.place_and_optimize_hydrogens(
    model=model_from_str(his_zn_str()), log=null_out(), do_flips=False,
    optimize_his_protonation=True)
  for a in m.get_hierarchy().atoms():
    if a.element.strip() != "H":
      # round-off only (1e-15 A)
      d = (matrix.col(a.xyz) - matrix.col(start[key(a)])).length()
      assert d < 1e-6, (key(a), d)
  d = near_zn(optimize_his_protonation=True)
  assert abs(d["CD2"] - 2.1) < 0.01 and d["HD2"] < 1.3, d
  print("test_ion_lock_without_flips OK")

def test_metal_contact_symmetry():
  """input_crystal_symmetry: placement's metal contacts across the model's
  own symmetry by default ('from_model', as before), none with None (a map
  box). The amine of His A 1 binds the Cu only across the P 1 box."""
  def nh2(**kw):
    report = _report()
    reduce_switch.place_and_optimize_hydrogens(
      model=model_from_str(his_cu_box_str()), log=null_out(), report=report,
      metal_bound_amines_as_nh2=True, raise_on_missing=False, **kw)
    return report.metal_amines_nh2
  across = nh2()
  assert len(across) == 1 and across[0][0].startswith("N   HIS A   1: "), \
    across
  assert " sym x-1,y,z" in across[0][1], across
  assert nh2(input_crystal_symmetry=None) == []
  print("test_metal_contact_symmetry OK")

if __name__ == "__main__":
  test_adds_hydrogens()
  test_dispatcher_explicit_reduce2()
  test_dispatcher_default_uses_switch()
  test_dispatcher_reduce1()
  test_get_nqh_flips_requires_h()
  test_get_nqh_flips_detects_flip()
  test_clashscore2_delegates()
  test_clashscore2_ignore_missing_restraints()
  test_clash_none_when_no_h_reduce2()
  test_clash_scores_existing_h_reduce2()
  test_flip_nqh_skips_h_less_model_reduce2()
  test_flip_nqh_applies_reduce2_flip()
  test_flip_nqh_use_ncs_without_groups_reduce2()
  test_rsr_per_cycle_clashscore_keeps_h_reduce2()
  test_molprobity_validation_adds_h_reduce2()
  test_reduce2_flips_applies_wholesale()
  test_his_protonation_without_flips()
  test_symmetry_copies_option()
  test_metal_amine_option()
  test_new_options_default_off()
  test_flips_from_final_model()
  test_ion_lock_without_flips()
  test_metal_contact_symmetry()
  print("OK")
