"""
Regression tests for autosol/autobuild bugs found in runs 107 and 109.

Bug 1  (Run 107): inject_user_params adds duplicate wavelength= that crashes PHIL
Bug 2  (Run 109): Directive extraction picks atom_type=S instead of Se
Bug 3  (Run 109): autosol re-runs after autosol+autobuild both succeeded

Fixes:
  1a  Strategy-flag alias awareness in inject_user_params
  1b  Boolean-type-mismatch error learning in bad_inject_params
  2a  Autosol atom_type/mad_ha_add_list deduplication in postprocess_command
  2b  Improved programs.yaml hint for atom_type (covered by YAML check)
  3   _is_program_already_done extended to non-count programs

Run with:
    PYTHONPATH=. python tests/tst_autosol_bugs.py
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
  assert_equal,
  assert_in,
  assert_not_in,
  assert_true,
  assert_false,
  run_tests_with_fail_fast,
)

# =========================================================================
# PHENIX/cctbx Linter "Silencer"
# =========================================================================
(re, assert_equal, assert_in, assert_not_in, assert_true, assert_false,
run_tests_with_fail_fast)


# =========================================================================
# Bug 1a: Strategy-flag alias awareness in inject_user_params
# =========================================================================

def test_bug1a_wavelength_alias_not_duplicated():
  """inject_user_params must NOT re-inject wavelength= when autosol.lambda=
    is already present (strategy_flags maps wavelength → autosol.lambda).

    Run 107: every autosol command had both autosol.lambda=0.9792 AND
    wavelength=0.9792 — PHIL interpreted bare wavelength= as
    autosol.wavelength.added_wavelength (boolean) → crash.
    """
  from agent.command_postprocessor import inject_user_params

  command = ("phenix.autosol autosol.data=p9.sca seq_file=seq.dat "
       "autosol.atom_type=Se mad_ha_add_list=S "
       "autosol.lambda=0.9792 resolution=2.5 nproc=4")

  user_advice = ("run autosol and autobuild at 2.5 A. "
         "use Se and S as anomalous atoms. wavelength is 0.9792")

  result = inject_user_params(
    command, user_advice, program_name="phenix.autosol")

  # The alias check should prevent bare "wavelength=0.9792" from being added
  assert_not_in("wavelength=0.9792", result,
         "inject_user_params should NOT add bare wavelength= when "
         "autosol.lambda= is already present (alias awareness)")

  # autosol.lambda should still be there
  assert_in("autosol.lambda=0.9792", result,
       "autosol.lambda=0.9792 should be preserved")


def test_bug1a_lambda_advice_not_duplicated():
  """When user says 'lambda is 0.9792', the word 'lambda' IS in
    autosol.lambda=0.9792, so the original duplicate check catches it.
    Verify this still works.
    """
  from agent.command_postprocessor import inject_user_params

  command = ("phenix.autosol autosol.data=p9.sca seq_file=seq.dat "
       "autosol.lambda=0.9792 nproc=4")

  user_advice = "run autosol. lambda is 0.9792"

  result = inject_user_params(
    command, user_advice, program_name="phenix.autosol")

  # Count occurrences of lambda — should be exactly one
  lambda_count = result.count("lambda=0.9792")
  assert_equal(lambda_count, 1,
        "lambda=0.9792 should appear exactly once (original dup check)")


def test_bug1a_alias_map_built_correctly():
  """Verify the alias map correctly identifies key → leaf mappings
    from strategy_flags. E.g., wavelength → lambda for autosol.
    """
  from agent.command_postprocessor import _load_prog_allowlist

  allowlist, sf = _load_prog_allowlist("phenix.autosol")

  # Build alias map the same way inject_user_params does
  alias_leaves = {}
  for sfkey, sfdef in sf.items():
    if isinstance(sfdef, dict):
      flag_tpl = sfdef.get('flag', '')
      leaf = flag_tpl.split('=')[0].strip().split('.')[-1].lower()
      if leaf and leaf != '{value}' and leaf != sfkey.lower():
        alias_leaves[sfkey.lower()] = leaf

  assert_in("wavelength", alias_leaves,
       "wavelength should be in alias map")
  assert_equal(alias_leaves["wavelength"], "lambda",
        "wavelength should alias to 'lambda'")


def test_bug1a_wavelength_not_in_allowlist():
  """Bare 'wavelength' must NOT be in the command-line
  allowlist.

  The strategy system maps wavelength -> autosol.lambda.
  The PHIL leaf 'lambda' is the valid bare parameter; bare
  'wavelength' is not a valid PHIL keyword and should be
  blocked by both sanitize_command and inject_user_params.
  """
  from agent.command_postprocessor import _load_prog_allowlist

  allowlist, _ = _load_prog_allowlist("phenix.autosol")

  assert_true("wavelength" not in allowlist,
    "bare 'wavelength' must NOT be in the "
    "allowlist (it is an alias for lambda)")
  assert_in("lambda", allowlist,
    "'lambda' (the PHIL leaf) must be in allowlist")


def test_bug1a_sanitize_strips_bare_wavelength():
  """sanitize_command must strip bare wavelength=xxx from
  autosol commands.

  The LLM or inject_user_params may produce bare
  wavelength= tokens.  Since 'wavelength' is not a valid
  PHIL parameter for autosol (correct form:
  autosol.lambda=xxx), Rule D should strip it.
  """
  from agent.command_postprocessor import sanitize_command

  command = (
    "phenix.autosol autosol.data=p9.sca "
    "seq_file=seq.dat autosol.atom_type=Se "
    "autosol.lambda=0.9792 resolution=2.5 nproc=4 "
    "wavelength=0.9000 wavelength=0.9794 "
    "wavelength=0.9797")

  result = sanitize_command(
    command, program_name="phenix.autosol")

  assert_not_in("wavelength=0.9000", result,
    "bare wavelength=0.9000 must be stripped")
  assert_not_in("wavelength=0.9794", result,
    "bare wavelength=0.9794 must be stripped")
  assert_not_in("wavelength=0.9797", result,
    "bare wavelength=0.9797 must be stripped")
  assert_in("autosol.lambda=0.9792", result,
    "autosol.lambda must be preserved")


def test_bug1a_inject_skips_wavelength_without_lambda():
  """inject_user_params must NOT inject bare wavelength=
  even when autosol.lambda is absent from the command.

  Previous behavior: the alias check only prevented
  injection when autosol.lambda was already present.  Now
  the allowlist itself excludes 'wavelength', so injection
  is blocked regardless.
  """
  from agent.command_postprocessor import inject_user_params

  # Command WITHOUT autosol.lambda -- previously this
  # would let inject_user_params add bare wavelength=
  command = (
    "phenix.autosol autosol.data=p9.sca "
    "seq_file=seq.dat autosol.atom_type=Se "
    "resolution=2.5 nproc=4")

  user_advice = (
    "MAD experiment with peak=0.9000 infl=0.9794 "
    "hrem=0.9797 wavelength=0.9000 "
    "wavelength=0.9794 wavelength=0.9797")

  result = inject_user_params(
    command, user_advice,
    program_name="phenix.autosol")

  assert_not_in("wavelength=0.9000", result,
    "bare wavelength must not be injected "
    "even without autosol.lambda in command")
  assert_not_in("wavelength=0.9794", result,
    "second wavelength must not be injected")
  assert_not_in("wavelength=0.9797", result,
    "third wavelength must not be injected")


# =========================================================================
# Bug 1b: Boolean-type-mismatch error learning
# =========================================================================

def test_bug1b_bool_type_error_blacklists_param():
  """When PHIL reports 'True or False value expected,
    autosol.wavelength.added_wavelength="0.9792" found', the parameter
    path and its components should be blacklisted.
    """
  error_msg = ('One True or False value expected, '
        'autosol.wavelength.added_wavelength="0.9792" found '
        '(command line argument, line 1)')

  # Extract using the same regex as ai_agent.py Fix 1b
  match = re.search(r'(\w+(?:\.\w+)*)="[^"]*"\s*found', error_msg)
  assert_true(match is not None,
        "Regex should match the PHIL error format")

  bad_full = match.group(1)
  assert_equal(bad_full, "autosol.wavelength.added_wavelength",
        "Should extract full PHIL path")

  # Check which components would be blacklisted (len >= 6)
  blacklisted = [part for part in bad_full.split('.') if len(part) >= 6]
  assert_in("autosol", blacklisted,
       "'autosol' (6 chars) should be blacklisted")
  assert_in("wavelength", blacklisted,
       "'wavelength' should be blacklisted")
  assert_in("added_wavelength", blacklisted,
       "'added_wavelength' should be blacklisted")


# =========================================================================
# Bug 2a: Autosol atom_type/mad_ha_add_list deduplication
# =========================================================================

def test_bug2a_duplicate_atom_types_stripped():
  """When atom_type=S and mad_ha_add_list=S (same value), the duplicate
    mad_ha_add_list should be removed.

    Run 109: directive extractor picked atom_type=S instead of Se,
    then the command had both set to S — selenium lost entirely.
    """
  from agent.command_postprocessor import postprocess_command

  command = ("phenix.autosol autosol.data=p9.sca seq_file=seq.dat "
       "autosol.atom_type=S mad_ha_add_list=S "
       "resolution=2.5 autosol.lambda=0.9792 nproc=4")

  result = postprocess_command(
    command, program_name="phenix.autosol",
    user_advice="use Se and S as anomalous atoms")

  assert_not_in("mad_ha_add_list=S", result,
         "Duplicate mad_ha_add_list (same as atom_type) should be removed")
  assert_in("atom_type=S", result,
       "atom_type=S should be preserved")


def test_bug2a_different_atom_types_preserved():
  """When atom_type=Se and mad_ha_add_list=S (different values),
    both should be preserved.
    """
  from agent.command_postprocessor import postprocess_command

  command = ("phenix.autosol autosol.data=p9.sca seq_file=seq.dat "
       "autosol.atom_type=Se mad_ha_add_list=S "
       "resolution=2.5 autosol.lambda=0.9792 nproc=4")

  result = postprocess_command(
    command, program_name="phenix.autosol",
    user_advice="use Se and S as anomalous atoms")

  assert_in("atom_type=Se", result,
       "atom_type=Se should be preserved")
  assert_in("mad_ha_add_list=S", result,
       "mad_ha_add_list=S should be preserved (differs from atom_type)")


# =========================================================================
# Bug 2b: programs.yaml hint quality
# =========================================================================

def test_bug2b_atom_type_hint_mentions_heavier():
  """The atom_type hint should mention using the heavier element."""
  try:
    from libtbx.langchain.knowledge.yaml_loader import get_program
  except ImportError:
    from knowledge.yaml_loader import get_program

  prog_def = get_program("phenix.autosol")
  assert_true(prog_def is not None,
        "phenix.autosol should be in programs.yaml")

  sf = prog_def.get("strategy_flags", {})
  at_hint = sf.get("atom_type", {}).get("hint", "")
  assert_in("heavier", at_hint.lower(),
       "atom_type hint should mention using the heavier element")

  ha_hint = sf.get("additional_atom_types", {}).get("hint", "")
  assert_in("differ", ha_hint.lower(),
       "additional_atom_types hint should say it must differ from atom_type")


# =========================================================================
# Bug 3: _is_program_already_done for non-count programs
# =========================================================================

def test_bug3_autosol_detected_as_done():
  """_is_program_already_done should return True for autosol when
    autosol_done=True in context, even though autosol uses set_flag
    strategy (not run_once).

    Run 109: after autosol+autobuild succeeded, _apply_directives
    re-added autosol from program_settings because _is_program_already_done
    only checked run_once programs.
    """
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  context = {
    "autosol_done": True,
    "autosol_success": True,
    "autobuild_done": True,
  }

  result = engine._is_program_already_done("phenix.autosol", context)
  assert_true(result,
        "autosol should be detected as done (set_flag + program-specific flag)")


def test_bug3_autobuild_detected_as_done():
  """autobuild (set_flag strategy) should also be caught."""
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  context = {"autobuild_done": True}

  result = engine._is_program_already_done("phenix.autobuild", context)
  assert_true(result,
        "autobuild should be detected as done")


def test_bug3_refine_not_blocked():
  """phenix.refine uses count strategy and should NOT be blocked
    by _is_program_already_done (it intentionally repeats).
    """
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  context = {
    "refine_done": True,
    "refine_count": 3,
  }

  result = engine._is_program_already_done("phenix.refine", context)
  assert_false(result,
        "refine (count strategy) should NOT be blocked")


def test_bug3_run_once_still_works():
  """run_once programs (e.g. xtriage) should still be caught
    by the original check.
    """
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  context = {"xtriage_done": True}

  result = engine._is_program_already_done("phenix.xtriage", context)
  assert_true(result,
        "xtriage (run_once) should be detected as done")


def test_bug3_not_done_returns_false():
  """When the done flag is not set, should return False."""
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  context = {"autosol_done": False}

  result = engine._is_program_already_done("phenix.autosol", context)
  assert_false(result,
        "autosol with done=False should NOT be blocked")


def test_bug3_shared_flag_not_blocked():
  """Programs with shared flags (e.g. validation_done) should NOT be
    blocked — only program-specific flags are checked.
    """
  from agent.workflow_engine import WorkflowEngine

  engine = WorkflowEngine()
  # molprobity uses "validation_done" flag which doesn't contain "molprobity"
  context = {"validation_done": True}

  result = engine._is_program_already_done("phenix.molprobity", context)
  assert_false(result,
        "molprobity with shared flag 'validation_done' should NOT be blocked "
        "(flag doesn't contain program short name)")


# =========================================================================
# Combined: end-to-end postprocess_command for autosol
# =========================================================================

def test_e2e_autosol_postprocess_no_wavelength_dup():
  """End-to-end: postprocess_command on a typical autosol command with
    user advice mentioning 'wavelength is 0.9792' should produce a command
    that has autosol.lambda but NOT bare wavelength.
    """
  from agent.command_postprocessor import postprocess_command

  command = ("phenix.autosol autosol.data=p9.sca seq_file=seq.dat "
       "autosol.atom_type=Se mad_ha_add_list=S "
       "autosol.lambda=0.9792 resolution=2.5 nproc=4")

  result = postprocess_command(
    command, program_name="phenix.autosol",
    user_advice="run autosol at 2.5 A. wavelength is 0.9792")

  assert_in("autosol.lambda=0.9792", result,
       "autosol.lambda should survive postprocessing")
  # Bare wavelength= must not appear
  bare = re.search(r'(?<!\.)wavelength=', result)
  assert_true(bare is None,
        "Bare wavelength= must NOT appear in final command "
        "(would crash PHIL as autosol.wavelength.added_wavelength)")


# =========================================================================
# Phase 1: Deterministic atom_type (heavier-atom-wins rule)
# =========================================================================

def test_phase1_swap_S_Se():
  """atom_type=S, mad_ha_add_list=Se → swapped (Z:16 < Z:34)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol autosol.atom_type=S mad_ha_add_list=Se resolution=2.5"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Se", result, "atom_type should become Se (heavier)")
  assert_in("mad_ha_add_list=S", result, "mad_ha_add_list should become S")


def test_phase1_no_swap_Se_S():
  """atom_type=Se, mad_ha_add_list=S → no change (already correct)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol autosol.atom_type=Se mad_ha_add_list=S resolution=2.5"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Se", result, "atom_type=Se should be unchanged")
  assert_in("mad_ha_add_list=S", result, "mad_ha_add_list=S should be unchanged")


def test_phase1_swap_S_Zn():
  """atom_type=S, mad_ha_add_list=Zn → swapped (Z:16 < Z:30)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol atom_type=S mad_ha_add_list=Zn"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Zn", result, "atom_type should become Zn (heavier)")
  assert_in("mad_ha_add_list=S", result, "mad_ha_add_list should become S")


def test_phase1_no_swap_Fe_S():
  """atom_type=Fe, mad_ha_add_list=S → no change (Z:26 > Z:16)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol atom_type=Fe mad_ha_add_list=S"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Fe", result, "atom_type=Fe should be unchanged")
  assert_in("mad_ha_add_list=S", result, "mad_ha_add_list=S should be unchanged")


def test_phase1_unknown_element_no_swap():
  """Unknown element Xx → no swap (do no harm)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol atom_type=Xx mad_ha_add_list=Se"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Xx", result, "Unknown element should not be swapped")
  assert_in("mad_ha_add_list=Se", result, "Se should remain in mad_ha_add_list")


def test_phase1_no_ha_list_no_change():
  """Only atom_type present, no mad_ha_add_list → no change."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol atom_type=Se resolution=2.5"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_equal(cmd, result, "Command should be unchanged without mad_ha_add_list")


def test_phase1_multi_element_swap():
  """atom_type=S, mad_ha_add_list=Se+Zn → swap S with heaviest (Se)."""
  from agent.command_postprocessor import _ensure_primary_scatterer_is_heavier

  cmd = "phenix.autosol atom_type=S mad_ha_add_list=Se+Zn"
  result = _ensure_primary_scatterer_is_heavier(cmd)
  assert_in("atom_type=Se", result,
       "atom_type should become Se (heaviest secondary)")
  assert_in("mad_ha_add_list=S+Zn", result,
       "mad_ha_add_list should have S replacing Se")


def test_phase1_dedup_then_swap():
  """Integration: dedup fires first (S==S → remove ha_list), then swap
    has nothing to do (no mad_ha_add_list left).
    """
  from agent.command_postprocessor import postprocess_command

  cmd = ("phenix.autosol autosol.data=p9.sca autosol.atom_type=S "
     "mad_ha_add_list=S autosol.lambda=0.9792")
  result = postprocess_command(cmd, program_name="phenix.autosol",
                user_advice="use S")
  assert_not_in("mad_ha_add_list", result,
         "Dedup should have removed duplicate mad_ha_add_list")
  assert_in("atom_type=S", result, "atom_type=S should remain")


def test_phase1_e2e_postprocess_swaps():
  """End-to-end: postprocess_command with wrong atom_type order → swapped."""
  from agent.command_postprocessor import postprocess_command

  cmd = ("phenix.autosol autosol.data=p9.sca autosol.atom_type=S "
     "mad_ha_add_list=Se autosol.lambda=0.9792 nproc=4")
  result = postprocess_command(
    cmd, program_name="phenix.autosol",
    user_advice="use Se and S as anomalous atoms. lambda is 0.9792")
  assert_in("atom_type=Se", result,
       "postprocess_command should swap atom_type to Se")
  assert_in("mad_ha_add_list=S", result,
       "postprocess_command should swap mad_ha_add_list to S")


def test_phase1_anomalous_z_table_coverage():
  """Verify the Big 5 + common derivatives are in _ANOMALOUS_Z."""
  from agent.command_postprocessor import _ANOMALOUS_Z

  big5 = {"S": 16, "Se": 34, "Zn": 30, "Fe": 26, "Hg": 80}
  for elem, z in big5.items():
    assert_in(elem, _ANOMALOUS_Z,
         "%s (Big 5) must be in _ANOMALOUS_Z" % elem)
    assert_equal(_ANOMALOUS_Z[elem], z,
          "%s should have Z=%d" % (elem, z))

  # Common heavy-atom derivatives
  for elem in ("Pt", "Au"):
    assert_in(elem, _ANOMALOUS_Z,
         "%s (heavy-atom derivative) must be in _ANOMALOUS_Z" % elem)


# =========================================================================
# Phase 2: Catch-all injection blacklist (streak tracker)
# =========================================================================

class _MockVlog:
  """Minimal mock for vlog."""
  def __init__(self):
    self.messages = []
  def normal(self, msg):
    self.messages.append(msg)


class _MockSession:
  """Minimal mock for session with bad_inject_params tracking."""
  def __init__(self):
    self.data = {}
    self.blacklisted = []  # (program, key) pairs

  def record_bad_inject_param(self, program, key):
    self.blacklisted.append((program, key))


class _MockAgent:
  """Minimal mock for the agent class, hosting _update_inject_fail_streak."""
  def __init__(self):
    self.vlog = _MockVlog()

  # Import the method from ai_agent.py would require PHENIX, so we inline
  # the logic here exactly as implemented.  The real tests validate that
  # the logic in ai_agent.py matches.
  def _update_inject_fail_streak(self, program, error_text, session,
                  is_success=False):
    import re as _re
    if session.data.get("force_retry_program"):
      return
    streaks = session.data.setdefault("inject_fail_streak", {})
    injected = session.data.get("last_injected_params", [])
    if is_success:
      if program in streaks:
        del streaks[program]
      return
    _raw = (error_text or "")[:120].lower()
    fingerprint = _re.sub(r'\d+', '', _raw).strip()
    if not fingerprint:
      return
    streak = streaks.get(program)
    if streak and streak.get("fingerprint") == fingerprint:
      streak["count"] = streak.get("count", 1) + 1
      _existing = set(streak.get("injected", []))
      for p in injected:
        if p not in _existing:
          streak.setdefault("injected", []).append(p)
          _existing.add(p)
    else:
      streaks[program] = {
        "count": 1,
        "fingerprint": fingerprint,
        "injected": list(injected),
      }
      streak = streaks[program]
    if streak["count"] >= 2 and streak.get("injected"):
      _to_blacklist = streak["injected"]
      self.vlog.normal(
        "  [catch-all blacklist] %s failed %dx — blacklisting: %s" %
        (program, streak["count"], ', '.join(_to_blacklist)))
      if hasattr(session, 'record_bad_inject_param'):
        for param in _to_blacklist:
          _key = param.split("=")[0].split(".")[-1]
          session.record_bad_inject_param(program, _key)
          _full = param.split("=")[0]
          if _full != _key:
            session.record_bad_inject_param(program, _full)
      del streaks[program]


def test_phase2_two_failures_triggers_blacklist():
  """Two consecutive failures with same fingerprint → blacklist injected."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["wavelength=0.9792"]

  error = 'True or False value expected, autosol.wavelength.added_wavelength="0.9792" found'

  # Failure 1
  agent._update_inject_fail_streak("phenix.autosol", error, session)
  assert_equal(len(session.blacklisted), 0,
        "No blacklisting after first failure")

  # Failure 2 — should trigger blacklist
  agent._update_inject_fail_streak("phenix.autosol", error, session)
  assert_true(len(session.blacklisted) > 0,
        "Should blacklist after second failure")
  keys = [k for _, k in session.blacklisted]
  assert_in("wavelength", keys,
       "wavelength (leaf of wavelength=0.9792) should be blacklisted")


def test_phase2_success_clears_streak():
  """One failure then success → streak cleared, nothing blacklisted."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["resolution=2.5"]

  error = "Some PHIL error message about resolution"

  # Failure 1
  agent._update_inject_fail_streak("phenix.autosol", error, session)
  assert_equal(len(session.blacklisted), 0)
  assert_in("phenix.autosol", session.data["inject_fail_streak"])

  # Success — clears streak
  agent._update_inject_fail_streak("phenix.autosol", "", session,
                   is_success=True)
  assert_not_in("phenix.autosol", session.data["inject_fail_streak"],
         "Success should clear streak")
  assert_equal(len(session.blacklisted), 0,
        "No blacklisting after success")


def test_phase2_empty_injected_no_blacklist():
  """Two failures but empty injected list → no blacklisting."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = []

  error = "Missing required input file"

  agent._update_inject_fail_streak("phenix.autosol", error, session)
  agent._update_inject_fail_streak("phenix.autosol", error, session)
  assert_equal(len(session.blacklisted), 0,
        "Empty injected list means nothing to blacklist")


def test_phase2_different_fingerprints_reset():
  """Two failures with different fingerprints → streak resets, no blacklist."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["wavelength=0.9792"]

  error1 = "True or False value expected"
  error2 = "Unknown command line parameter: wavelength"

  agent._update_inject_fail_streak("phenix.autosol", error1, session)
  agent._update_inject_fail_streak("phenix.autosol", error2, session)
  assert_equal(len(session.blacklisted), 0,
        "Different fingerprints should not trigger blacklist")
  # Streak should be reset to count=1 with the new fingerprint
  streak = session.data["inject_fail_streak"].get("phenix.autosol", {})
  assert_equal(streak.get("count"), 1,
        "Streak should reset to 1 on different fingerprint")


def test_phase2_force_retry_skips_streak():
  """Force-retry (recovery) cycles don't increment streak."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["wavelength=0.9792"]
  session.data["force_retry_program"] = "phenix.autosol"

  error = "True or False value expected"

  agent._update_inject_fail_streak("phenix.autosol", error, session)
  agent._update_inject_fail_streak("phenix.autosol", error, session)
  assert_equal(len(session.blacklisted), 0,
        "Force-retry cycles should not trigger blacklist")
  assert_equal(session.data.get("inject_fail_streak", {}), {},
        "Streak should not be created during force-retry")


def test_phase2_innocent_bystander_blacklisted():
  """Harmless param (nproc=8) gets blacklisted alongside offender —
    acceptable trade-off to break the loop.
    """
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["wavelength=0.9792", "nproc=8"]

  error = 'True or False value expected, autosol.wavelength.added_wavelength="0.9792" found'

  agent._update_inject_fail_streak("phenix.autosol", error, session)
  agent._update_inject_fail_streak("phenix.autosol", error, session)

  keys = [k for _, k in session.blacklisted]
  assert_in("wavelength", keys, "wavelength should be blacklisted")
  assert_in("nproc", keys, "nproc (innocent bystander) also blacklisted")


def test_phase2_dotted_key_blacklists_both():
  """Dotted injected param blacklists both full path and leaf."""
  agent = _MockAgent()
  session = _MockSession()
  session.data["last_injected_params"] = ["autosol.lambda=0.9792"]

  error = "Some error about lambda"

  agent._update_inject_fail_streak("phenix.autosol", error, session)
  agent._update_inject_fail_streak("phenix.autosol", error, session)

  keys = [k for _, k in session.blacklisted]
  assert_in("lambda", keys, "leaf 'lambda' should be blacklisted")
  assert_in("autosol.lambda", keys,
       "full key 'autosol.lambda' should also be blacklisted")


def test_phase2_return_injected_backward_compat():
  """postprocess_command default returns str (not tuple)."""
  from agent.command_postprocessor import postprocess_command

  result = postprocess_command("phenix.refine model.pdb data.mtz",
                "phenix.refine")
  assert_true(isinstance(result, str),
        "Default return should be str for backward compat")


def test_phase2_return_injected_captures_params():
  """postprocess_command with return_injected=True returns (str, list)."""
  from agent.command_postprocessor import postprocess_command

  cmd = "phenix.autosol autosol.data=p9.sca autosol.atom_type=Se"
  result, injected = postprocess_command(
    cmd, program_name="phenix.autosol",
    user_advice="resolution=2.5", return_injected=True)
  assert_true(isinstance(result, str), "tuple[0] should be str")
  assert_true(isinstance(injected, list), "tuple[1] should be list")
  # resolution=2.5 should have been injected
  assert_true(any("resolution" in p for p in injected),
        "resolution should be in injected list")


# =========================================================================
# Phase 2 supplemental: directive-level atom_type validation
# =========================================================================

def test_directive_atom_type_swap():
  """Session._fix_autosol_atom_type_order swaps S→Se at directive level."""
  try:
    from agent.session import AgentSession
  except ImportError:
    print("  SKIP (session import unavailable)")
    return

  directives = {
    "program_settings": {
      "phenix.autosol": {
        "atom_type": "S",
        "additional_atom_types": "Se",
      }
    }
  }

  logs = []
  result = AgentSession._fix_autosol_atom_type_order(
    directives, log=lambda m: logs.append(m))

  at = result["program_settings"]["phenix.autosol"]["atom_type"]
  ha = result["program_settings"]["phenix.autosol"]["additional_atom_types"]
  assert_equal(at, "Se", "atom_type should be swapped to Se (heavier)")
  assert_equal(ha, "S", "additional_atom_types should be swapped to S")
  assert_true(len(logs) > 0, "Should log the swap")


def test_directive_atom_type_no_swap_correct():
  """No swap when atom_type is already heavier."""
  try:
    from agent.session import AgentSession
  except ImportError:
    print("  SKIP (session import unavailable)")
    return

  directives = {
    "program_settings": {
      "phenix.autosol": {
        "atom_type": "Se",
        "additional_atom_types": "S",
      }
    }
  }

  result = AgentSession._fix_autosol_atom_type_order(directives)
  at = result["program_settings"]["phenix.autosol"]["atom_type"]
  assert_equal(at, "Se", "atom_type=Se should be unchanged (already correct)")


def test_directive_atom_type_no_autosol_noop():
  """No crash when directives don't mention autosol."""
  try:
    from agent.session import AgentSession
  except ImportError:
    print("  SKIP (session import unavailable)")
    return

  directives = {
    "program_settings": {
      "phenix.refine": {"cycles": 5}
    }
  }

  result = AgentSession._fix_autosol_atom_type_order(directives)
  assert_equal(result["program_settings"]["phenix.refine"]["cycles"], 5,
        "Non-autosol directives should be untouched")


# =========================================================================
# Bug 4: rebuild_in_place=False stripped by Rule D (v112.77)
# =========================================================================

def test_bug4_rebuild_in_place_survives_sanitize():
  """rebuild_in_place=False must NOT be stripped by sanitize_command."""
  from agent.command_postprocessor import sanitize_command

  cmd = ("phenix.autobuild data=/path/data.mtz seq_file=/path/seq.dat "
     "model=/path/model.pdb rebuild_in_place=False resolution=2.5 nproc=4")

  result = sanitize_command(cmd, program_name="phenix.autobuild")
  assert_in("rebuild_in_place=False", result,
       "rebuild_in_place=False should survive Rule D")


def test_bug4_rebuild_in_place_survives_postprocess():
  """rebuild_in_place=False must survive full postprocess_command."""
  from agent.command_postprocessor import postprocess_command

  cmd = ("phenix.autobuild data=/path/data.mtz seq_file=/path/seq.dat "
     "rebuild_in_place=False resolution=2.5")

  result = postprocess_command(cmd, program_name="phenix.autobuild")
  assert_in("rebuild_in_place=False", result,
       "rebuild_in_place=False should survive postprocess_command")


def test_bug4_n_cycle_build_max_survives():
  """n_cycle_build_max (added alongside rebuild_in_place) should survive."""
  from agent.command_postprocessor import sanitize_command

  cmd = ("phenix.autobuild data=/path/data.mtz seq_file=/path/seq.dat "
     "n_cycle_build_max=5 nproc=4")

  result = sanitize_command(cmd, program_name="phenix.autobuild")
  assert_in("n_cycle_build_max=5", result,
       "n_cycle_build_max should survive Rule D")


def test_bug4_maps_only_survives():
  """maps_only (added alongside rebuild_in_place) should survive."""
  from agent.command_postprocessor import sanitize_command

  cmd = ("phenix.autobuild data=/path/data.mtz seq_file=/path/seq.dat "
     "maps_only=True nproc=4")

  result = sanitize_command(cmd, program_name="phenix.autobuild")
  assert_in("maps_only=True", result,
       "maps_only should survive Rule D")


def test_bug4_autobuild_allowlist_coverage():
  """Verify autobuild allowlist includes all 6 strategy_flags."""
  from agent.command_postprocessor import _load_prog_allowlist

  allowlist, sf = _load_prog_allowlist("phenix.autobuild")
  expected = {"quick", "nproc", "resolution",
        "rebuild_in_place", "n_cycle_build_max", "maps_only"}
  for key in expected:
    assert_in(key, allowlist,
         "%s must be in autobuild allowlist" % key)


def test_bug4_unknown_bare_still_stripped():
  """Unknown bare params like hallucinated_param=42 should still be stripped
    by Rule D — we only added specific known params.
    """
  from agent.command_postprocessor import sanitize_command

  cmd = ("phenix.autobuild data=/path/data.mtz seq_file=/path/seq.dat "
     "hallucinated_param=42 nproc=4")
  logs = []
  result = sanitize_command(cmd, program_name="phenix.autobuild",
               log=lambda m: logs.append(m))
  assert_not_in("hallucinated_param", result,
         "Unknown bare params should still be stripped by Rule D")
  assert_true(any("hallucinated_param" in l for l in logs),
        "Rule D should log that it stripped hallucinated_param")


def test_bug4_hint_mentions_rebuild_in_place():
  """programs.yaml hints should mention rebuild_in_place recovery."""
  from knowledge.yaml_loader import get_program

  defn = get_program("phenix.autobuild")
  hints = defn.get("hints", [])
  hint_text = " ".join(hints).lower()
  assert_in("rebuild_in_place", hint_text,
       "Autobuild hints should mention rebuild_in_place recovery")


# =========================================================================
# Mock drift detection: verify real _update_inject_fail_streak matches mock
# =========================================================================
# The Phase 2 streak tests use _MockAgent which inlines the streak logic.
# This structural test reads the real ai_agent.py source and verifies key
# logic invariants are present, catching silent drift.
# =========================================================================

def test_phase2_mock_drift_check():
  """Verify ai_agent.py _update_inject_fail_streak has expected structure.

    Reads source code and checks for key logic lines.  If ai_agent.py is
    not found, the test SKIPs (non-PHENIX environment).
    """
  _project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  ai_agent_path = os.path.join(_project_root, "programs", "ai_agent.py")
  if not os.path.isfile(ai_agent_path):
    print("  SKIP (ai_agent.py not found at %s)" % ai_agent_path)
    return

  with open(ai_agent_path) as f:
    source = f.read()

  # Verify the method exists
  assert_in("def _update_inject_fail_streak", source,
       "ai_agent.py must contain _update_inject_fail_streak")

  # Verify N=2 threshold
  assert_in('["count"] >= 2', source,
       "Streak threshold must be >= 2")

  # Verify force_retry exclusion
  assert_in("force_retry_program", source,
       "Must skip streak during force_retry")

  # Verify fingerprint normalization (digits stripped)
  assert_in("sub(r'\\d+'", source,
       "Fingerprint must strip digits")

  # Verify it's called from _record_command_result
  assert_in("_update_inject_fail_streak", source,
       "_record_command_result must call _update_inject_fail_streak")


# ── Bug 5: map_coeffs_mtz empty after refine (GUI mode) ────────────────

def test_bug5_safety_net_handles_empty_best_files():
  """Safety net must fire when best_files is {} (empty dict).

    The condition was `if best_files and not best_files.get(...)` which
    skips {} because empty dicts are falsy.  Fixed to use `is not None`.
    """
  ai_agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(ai_agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(ai_agent_path, 'r') as f:
    source = f.read()

  # The safety net condition must use `is not None`, not truthiness
  assert_in("best_files is not None", source,
       "Safety net must use 'is not None' to handle empty dict")

  # Must NOT have the old falsy-dict-skipping pattern
  # (Check that 'if best_files and not best_files.get("map_coeffs_mtz")'
  #  is not present — the 'is not None' form supersedes it)
  import re
  old_pattern = re.search(
    r'if best_files and not best_files\.get\(["\']map_coeffs_mtz',
    source)
  assert_true(old_pattern is None,
        "Old 'if best_files and' pattern must be replaced with "
        "'best_files is not None'")


def test_bug5_gui_sub_job_returns_output_dir():
  """_execute_sub_job_for_gui must return 4-tuple including output_dir.

    Without this, _execute_command uses os.getcwd() which points to the
    parent agent directory after the OldStyle runner restores CWD.
    """
  ai_agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(ai_agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(ai_agent_path, 'r') as f:
    source = f.read()

  # The main return must include gui_output_dir
  assert_in("return log_text, error_text, executed_command, gui_output_dir",
       source,
       "_execute_sub_job_for_gui must return gui_output_dir as 4th element")


def test_bug5_execute_command_uses_gui_output_dir():
  """_execute_command must use gui_output_dir (not os.getcwd()) in GUI mode.

    The OldStyle runner restores CWD to the parent directory after execution,
    so os.getcwd() points to the wrong place.  output files are in the
    sub-job directory (e.g., sub_03_refine/).
    """
  ai_agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(ai_agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(ai_agent_path, 'r') as f:
    source = f.read()

  # Must unpack 4-tuple from _execute_sub_job_for_gui
  assert_in("gui_output_dir", source,
       "_execute_command must capture gui_output_dir from sub-job")

  # working_dir must prefer gui_output_dir over os.getcwd()
  assert_in("gui_output_dir if gui_output_dir else log_dir",
       source,
       "working_dir must use gui_output_dir when available")


def test_bug5_track_output_files_accepts_working_dir():
  """_track_output_files must accept working_dir param and use it for scanning.

    Without this, it falls back to os.getcwd() which is the parent agent
    directory in GUI mode.
    """
  ai_agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(ai_agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(ai_agent_path, 'r') as f:
    source = f.read()

  # Function signature must accept working_dir
  assert_in("def _track_output_files(self, log_text, session_start_time, "
       "session=None,\n                          cycle=0, working_dir=None)",
       source,
       "_track_output_files must accept working_dir parameter")

  # Must use working_dir for directory scanning, not os.getcwd()
  assert_in("scan_dir = working_dir or os.getcwd()", source,
       "_track_output_files must prefer working_dir over os.getcwd()")

  # Call site must pass working_dir
  assert_in("working_dir=working_dir)", source,
       "_track_output_files call must pass working_dir")


# ── Bug 6: Daily usage limit not raised as Sorry ───────────────────────

def test_bug6_rest_init_raises_sorry_on_daily_limit():
  """rest/__init__.py must raise Sorry on daily_usage_reached."""
  rest_init_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "rest", "__init__.py")
  if not os.path.isfile(rest_init_path):
    print("  SKIP (rest/__init__.py not found)")
    return
  with open(rest_init_path, 'r') as f:
    source = f.read()

  import re
  sorry_raises = re.findall(
    r"daily_usage_reached.*?raise Sorry", source, re.DOTALL)
  assert_true(len(sorry_raises) >= 2,
        "Both daily_usage_reached detection points must raise Sorry "
        "(found %d)" % len(sorry_raises))


def test_bug6_remote_agent_reraises_sorry():
  """RemoteAgent must re-raise Sorry instead of swallowing it.

    The generic `except Exception` around _send_request caught Sorry
    (which inherits from Exception), logged it, and returned None.
    Fix: add `except Sorry: raise` before the generic handler.
    """
  remote_agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "phenix_ai", "remote_agent.py")
  if not os.path.isfile(remote_agent_path):
    print("  SKIP (remote_agent.py not found)")
    return
  with open(remote_agent_path, 'r') as f:
    source = f.read()

  assert_in("except Sorry:", source,
       "RemoteAgent must catch Sorry separately")
  assert_in("raise  # Fatal server errors", source,
       "RemoteAgent must re-raise Sorry")

  # Sorry must be imported
  assert_in("from libtbx.utils import Sorry", source,
       "RemoteAgent must import Sorry")


# ── Bug 7: after_program premature stop on multi-goal requests ──────────

def test_bug7_after_program_not_hard_stop():
  """check_directive_stop must NOT stop when after_program's target completes.

    after_program is now a minimum-run guarantee (enforced in PLAN), not a
    hard stop (was in PERCEIVE).  The LLM decides when all goals are met.
    """
  sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "agent"))
  from perceive_checks import check_directive_stop

  directives = {"stop_conditions": {"after_program": "phenix.map_symmetry"}}
  history = [
    {"program": "phenix.mtriage", "command": "phenix.mtriage ...", "result": "SUCCESS"},
    {"program": "phenix.map_symmetry", "command": "phenix.map_symmetry ...", "result": "SUCCESS"},
  ]

  should_stop, reason = check_directive_stop(directives, history, cycle_number=3)
  assert_true(should_stop is False,
        "after_program must NOT hard-stop in PERCEIVE "
        "(got should_stop=%s, reason=%s)" % (should_stop, reason))


def test_bug7_after_cycle_still_hard_stops():
  """after_cycle must still work as a hard stop (regression guard)."""
  sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "agent"))
  from perceive_checks import check_directive_stop

  directives = {"stop_conditions": {"after_cycle": 2}}
  history = [
    {"program": "phenix.mtriage", "result": "SUCCESS"},
    {"program": "phenix.map_symmetry", "result": "SUCCESS"},
  ]

  should_stop, reason = check_directive_stop(directives, history, cycle_number=3)
  assert_true(should_stop is True,
        "after_cycle must still hard-stop (got should_stop=%s)" % should_stop)
  assert_in("after_cycle", reason,
       "Reason must mention after_cycle")


def test_bug7_metrics_targets_still_hard_stop():
  """r_free_target and map_cc_target must still work as hard stops."""
  sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "agent"))
  from perceive_checks import check_directive_stop

  # r_free target met
  directives = {"stop_conditions": {"r_free_target": 0.25}}
  history = [{"program": "phenix.refine", "result": "SUCCESS"}]
  metrics = {"r_free": 0.22}

  should_stop, reason = check_directive_stop(
    directives, history, cycle_number=2, current_metrics=metrics)
  assert_true(should_stop is True,
        "r_free_target must still hard-stop (got %s)" % should_stop)
  assert_in("R-free", reason, "Reason must mention R-free")

  # map_cc target met
  directives = {"stop_conditions": {"map_cc_target": 0.80}}
  metrics = {"map_cc": 0.85}
  should_stop, reason = check_directive_stop(
    directives, history, cycle_number=2, current_metrics=metrics)
  assert_true(should_stop is True,
        "map_cc_target must still hard-stop (got %s)" % should_stop)
  assert_in("Map CC", reason, "Reason must mention Map CC")


def test_bug7_source_no_after_program_return_true():
  """perceive_checks.py must not contain after_program hard-stop logic."""
  checks_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "perceive_checks.py")
  with open(checks_path, 'r') as f:
    source = f.read()

  # The old pattern that returned True on after_program match must be gone
  assert_true(
    'return True, "Completed %s (directive: after_program)"' not in source,
    "Old after_program hard-stop return must be removed")

  # The comment explaining WHY it was removed must exist
  assert_in("intentionally NOT a hard stop", source,
       "Must document why after_program is not a hard stop")


def test_bug7_session_fallback_no_after_program_stop():
  """session.py fallback check_directive_stop_conditions must not hard-stop
    on after_program either.

    This is a parallel implementation used when directive_extractor can't
    be imported.  Must match perceive_checks.py semantics.
    """
  session_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "session.py")
  with open(session_path, 'r') as f:
    source = f.read()

  # The old pattern that set should_stop=True on after_program match
  assert_true(
    'reason = "Completed %s (directive)" % target_program' not in source,
    "session.py fallback must not hard-stop on after_program")

  # Should have the same "intentionally NOT a hard stop" note
  assert_in("intentionally NOT a hard stop", source,
       "session.py must document after_program is not a hard stop")


def test_bug7_no_dead_last_command_in_perceive_checks():
  """perceive_checks.py should not have unused last_command variable.

    last_command was only used by the after_program block that was removed.
    """
  checks_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "perceive_checks.py")
  with open(checks_path, 'r') as f:
    source = f.read()

  assert_true(
    "last_command" not in source,
    "last_command is dead code after after_program removal")


def test_bug7_directive_extractor_no_after_program_stop():
  """directive_extractor.check_stop_conditions must not hard-stop on
    after_program.

    This is the server-side implementation called by session.py.  Must
    match perceive_checks.py semantics.
    """
  extractor_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "directive_extractor.py")
  with open(extractor_path, 'r') as f:
    source = f.read()

  assert_true(
    'return True, "Completed %s (directive: after_program)" % target_program'
    not in source,
    "directive_extractor must not hard-stop on after_program")

  assert_in("intentionally NOT a hard stop", source,
       "directive_extractor must document why after_program is not a hard stop")


# ===================================================================
# WINDOWS COMPATIBILITY TESTS
# ===================================================================


def test_win_filter_intermediate_normalizes_separators():
  """_filter_intermediate_files must normalize path separators before
    checking temp dir markers.

    On Windows, paths use backslashes (C:\\Users\\...\\TEMP0\\file.pdb).
    The forward-slash markers (/TEMP, /TEMP0/) only match if the path
    is normalized first.
    """
  nodes_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "graph_nodes.py")
  with open(nodes_path, 'r') as f:
    source = f.read()

  # The function must normalize backslashes before matching markers
  assert_in('replace("\\\\", "/")', source,
       "_filter_intermediate_files must normalize \\\\ to / before "
       "checking temp_dir_markers")


def test_win_popen_create_no_window():
  """_run_command_tracked should use CREATE_NO_WINDOW on Windows.

    Source-level check: the Popen call must include creationflags for
    os.name == 'nt' to prevent console windows flashing in GUI mode.
    """
  agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(agent_path, 'r') as f:
    source = f.read()

  assert_in("CREATE_NO_WINDOW", source,
       "Popen must use CREATE_NO_WINDOW on Windows")
  assert_in("os.name == 'nt'", source,
       "CREATE_NO_WINDOW must be conditional on Windows")


def test_win_abort_detection_comment():
  """Abort detection must document Windows return code behavior.

    On Unix, killed processes have return_code < 0.  On Windows,
    taskkill /F produces return_code >= 0.  STOPWIZARD is the reliable
    cross-platform indicator.
    """
  agent_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "programs", "ai_agent.py")
  if not os.path.isfile(agent_path):
    print("  SKIP (ai_agent.py not found)")
    return
  with open(agent_path, 'r') as f:
    source = f.read()

  assert_in("taskkill", source,
       "Must document Windows taskkill behavior")
  assert_in("cross-platform", source,
       "Must note STOPWIZARD as cross-platform indicator")


def test_win_session_utf8_encoding():
  """session.py must use explicit UTF-8 encoding for file I/O.

    On Windows, the default locale encoding (e.g. cp1252) can fail on
    non-ASCII characters in paths or user advice.
    """
  session_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "agent", "session.py")
  with open(session_path, 'r') as f:
    source = f.read()

  # All three open() calls for session JSON should use encoding='utf-8'
  import re
  open_calls = re.findall(r"open\([^)]*session_file[^)]*\)", source)
  for call in open_calls:
    assert_in("encoding='utf-8'", call,
         "session_file open() must specify UTF-8: %s" % call)


def run_all_tests():
  """Run all autosol bugs tests."""
  run_tests_with_fail_fast()


# =========================================================================

if __name__ == "__main__":
  success = run_tests_with_fail_fast()
  if not success:
    sys.exit(1)
  print("\nOK")
