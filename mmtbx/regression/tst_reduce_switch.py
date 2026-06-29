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

def test_flip_nqh_respects_selection_reduce2():
  """Under reduce2, flip_nqh(selection=...) detects on a SUB-model (model.select)
  and applies/propagates the flip only within the selection."""
  import os, re, libtbx.load_env
  from io import StringIO
  p = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/ions/3zu8.pdb", test=os.path.isfile)
  if p is None:
    print("3zu8.pdb not found, skipping selection test"); return
  saved = reduce_switch.USE_OLD_REDUCE
  try:
    reduce_switch.USE_OLD_REDUCE = False
    dm = DataManager(['model']); dm.process_model_file(p)
    m = dm.get_model(); m.add_crystal_symmetry_if_necessary()
    assert _make_all_asn_backwards(m.get_hierarchy()) > 0
    m.set_sites_cart(m.get_hierarchy().atoms().extract_xyz())
    hm = _h_restrained_model(m)
    sel = hm.selection("chain A")                 # contains the backwards ASN
    assert sel.count(True) > 0
    sites_before = hm.get_sites_cart().deep_copy()
    log = StringIO(); hm.set_log(log)
    nqh.flip(model=hm, selection=sel, log=hm.log)  # reduce2 detects on model.select(sel)
    mo = re.search(r"Total number of N/Q/H flips:\s*(\d+)", log.getvalue())
    assert mo and int(mo.group(1)) >= 1, log.getvalue()
    moved = (hm.get_sites_cart() - sites_before).norms() > 0.1
    assert moved.count(True) > 0                   # the flip propagated to the model
    assert (moved & ~sel).count(True) == 0        # only selection atoms moved
  finally:
    reduce_switch.USE_OLD_REDUCE = saved
  print("test_flip_nqh_respects_selection_reduce2 OK")

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

def test_model_idealization_stat_display_tolerates_none_reduce2():
  """phenix.model_idealization works on an H-less model (it runs remove_hydrogens),
  so under reduce2 every geometry_statistics() snapshot has clash.score=None and
  molprobity_score=None. print_stat_comparison() must render those as the program's
  own 99999 'not available' sentinel instead of crashing on "{:10.2f}".format(None)
  (the contract: all clash()/mp_score() callers must tolerate None)."""
  from io import StringIO
  from libtbx import group_args
  from mmtbx.command_line import model_idealization as mi
  def fake_stats():
    # mimic geometry_statistics().result() under reduce2: clash/molprobity are None,
    # the other validation fields are real numbers (reduce2 does not touch them).
    return group_args(
      molprobity_score = None,
      clash        = group_args(score=None, clashes=None),
      c_beta       = group_args(outliers=0),
      ramachandran = group_args(outliers=0, allowed=0, favored=100),
      rotamer      = group_args(outliers=0),
      omega        = group_args(cis_proline=0, cis_general=0,
                                twisted_proline=0, twisted_general=0),
      cablam       = group_args(outliers=0, disfavored=0, ca_outliers=0))
  # Build the reporter without the heavy __init__ (needs a real model+map); we only
  # exercise the stat-printing method, feeding it the reduce2-shaped stat objects.
  obj = mi.model_idealization.__new__(mi.model_idealization)
  obj.log = StringIO()
  obj.params = group_args(run_minimization_first=False)
  obj.after_cablam_statistics  = None
  obj.init_gm_model_statistics = None
  obj.init_model_statistics    = fake_stats()
  obj.after_ss_idealization    = fake_stats()
  obj.after_loop_idealization  = fake_stats()
  obj.after_rotamer_fixing     = fake_stats()
  obj.final_model_statistics   = fake_stats()
  obj.time_for_init = 0; obj.time_for_run = 0
  obj.get_rmsd_from_start  = lambda: 0.0  # shadow the heavy rmsd calc (unrelated)
  obj.get_rmsd_from_start2 = lambda: 0.0
  obj.print_stat_comparison()             # must NOT raise on the None values
  out = obj.log.getvalue()
  clash_lines = [l for l in out.splitlines() if l.startswith("Clashscore")]
  assert len(clash_lines) == 1, out
  assert "99999" in clash_lines[0], clash_lines[0]   # None -> sentinel, not a crash
  print("test_model_idealization_stat_display_tolerates_none_reduce2 OK")

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

if __name__ == "__main__":
  test_adds_hydrogens()
  test_dispatcher_explicit_reduce2()
  test_dispatcher_default_uses_switch()
  test_dispatcher_reduce1()
  test_get_nqh_flips_requires_h()
  test_get_nqh_flips_detects_flip()
  test_clashscore2_delegates()
  test_clash_none_when_no_h_reduce2()
  test_clash_scores_existing_h_reduce2()
  test_flip_nqh_skips_h_less_model_reduce2()
  test_flip_nqh_applies_reduce2_flip()
  test_flip_nqh_respects_selection_reduce2()
  test_rsr_per_cycle_clashscore_keeps_h_reduce2()
  test_model_idealization_stat_display_tolerates_none_reduce2()
  test_molprobity_validation_adds_h_reduce2()
  print("OK")
