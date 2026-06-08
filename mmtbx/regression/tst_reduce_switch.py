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

if __name__ == "__main__":
  test_adds_hydrogens()
  test_dispatcher_explicit_reduce2()
  test_dispatcher_default_uses_switch()
  test_dispatcher_reduce1()
  test_get_nqh_flips_requires_h()
  test_get_nqh_flips_detects_flip()
  test_clashscore2_delegates()
  print("OK")
