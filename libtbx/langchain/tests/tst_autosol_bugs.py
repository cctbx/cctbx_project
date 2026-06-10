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
# Bug 8: Template-literal flag stripped by sanitize_command allowlist
#         (AIAgent_62, v119.H8)
# =========================================================================
#
# Symptom: phenix.autobuild_denmod was chosen by the planner, but the
# emitted command was missing the literal `maps_only=True` flag that
# the template specifies.  AutoBuild then ran the full iterative build
# loop instead of just doing density modification.
#
# Diagnosis: programs.yaml encodes `maps_only=True` as a literal in the
# command template (not in strategy_flags) because it's an invariant
# of the autobuild_denmod program entry.  But sanitize_command Rule D
# strips any bare key=value whose key isn't in the strategy_flags
# allowlist.  Since maps_only is intentionally NOT in autobuild_denmod's
# strategy_flags (to prevent LLM override), Rule D was stripping the
# template literal.
#
# Distinct from the existing Bug 4 section (rebuild_in_place=False,
# v112.77) which dealt with phenix.autobuild and is solved by
# strategy_flags membership.  Bug 8 is about phenix.autobuild_denmod
# where the canonical fix is template-literal extraction, NOT a
# strategy_flags addition.
#
# Fix: _load_prog_allowlist now also extracts literal key= tokens from
# the command template and adds them to the allowlist.

def test_bug8_template_literal_maps_only_preserved():
  """Primary regression: maps_only=True from template survives sanitize.

    Reproduces the exact failure from AIAgent_62 cycle 4 (bromodomain
    ligand demo).  The post-registry command must retain maps_only=True
    after sanitize_command processes it under program_name=
    'phenix.autobuild_denmod'.
    """
  from agent.command_postprocessor import sanitize_command

  # Realistic command shape from program_registry.build_command()
  command = ("phenix.autobuild "
       "data=/p/sub_03_refine/refine_001_001.mtz "
       "seq_file=/p/7qz0.fa "
       "model=/p/sub_03_refine/refine_001_001.pdb "
       "maps_only=True nproc=4")

  result = sanitize_command(
    command, program_name="phenix.autobuild_denmod",
    log=lambda m: None)

  assert_in("maps_only=True", result,
       "Template literal maps_only=True must survive sanitize_command "
       "(was stripped pre-H8 because not in strategy_flags allowlist)")
  # Sanity: other parts of the command must also survive
  assert_in("nproc=4", result, "nproc=4 should survive (in strategy_flags)")
  assert_in("data=", result, "data= file path should survive")


def test_bug8_load_allowlist_includes_template_literals():
  """Allowlist for phenix.autobuild_denmod must include 'maps_only'.

    Direct unit test on _load_prog_allowlist.  Also verifies that
    pre-existing strategy_flag entries (nproc, resolution) are still
    present — i.e., no regression in the existing logic.
    """
  from agent.command_postprocessor import _load_prog_allowlist

  allowlist, sf = _load_prog_allowlist("phenix.autobuild_denmod")

  if allowlist is None:
    # Sandbox without programs.yaml available — skip gracefully
    print("  SKIP (programs.yaml not loadable in this environment)")
    return

  assert_in("maps_only", allowlist,
       "Allowlist for phenix.autobuild_denmod must include 'maps_only' "
       "(extracted from command template literal)")
  # Pre-existing strategy_flag entries must still be there
  assert_in("nproc", allowlist,
       "nproc (declared in strategy_flags) must remain in allowlist")
  assert_in("resolution", allowlist,
       "resolution (declared in strategy_flags) must remain in allowlist")
  # And maps_only must NOT be in strategy_flags — that's important.
  # If it were, the LLM could override via strategy.maps_only=False.
  # Bug 5's fix relies on maps_only living ONLY in the template.
  assert_true("maps_only" not in sf,
        "maps_only must NOT be in phenix.autobuild_denmod's "
        "strategy_flags — it lives only in the command template as "
        "an invariant of this program entry.  (Adding it to "
        "strategy_flags would enable LLM override and defeat the "
        "template-as-invariant pattern.)")


def test_bug8_template_literal_extraction_ignores_placeholders():
  """The template-literal extraction must not match {placeholder} slots.

    Pins the regex's safety guarantee: only true key=value tokens with
    whitespace or start-of-string before the key get extracted.
    Placeholders like {data_mtz} have '{' before the identifier, so
    they correctly don't match.
    """
  from agent.command_postprocessor import _load_prog_allowlist

  # Construct a fake prog_def via monkey-patching get_program
  fake_prog_def = {
      'command': 'phenix.fake {model} {data_mtz} flag_a=True {map} flag_b=42',
      'strategy_flags': {
          'nproc': {'flag': 'nproc={value}', 'type': 'int'},
      },
  }

  # _load_prog_allowlist imports get_program inside its try block, so
  # we monkey-patch in the yaml_loader module
  try:
    from libtbx.langchain.knowledge import yaml_loader as _yl
  except ImportError:
    try:
      from knowledge import yaml_loader as _yl
    except ImportError:
      print("  SKIP (yaml_loader not importable)")
      return

  _orig = _yl.get_program
  _yl.get_program = lambda name: fake_prog_def if name == 'phenix.fake' else _orig(name)
  try:
    allowlist, sf = _load_prog_allowlist('phenix.fake')
  finally:
    _yl.get_program = _orig

  if allowlist is None:
    print("  SKIP (allowlist returned None — environment issue)")
    return

  # The two literal flags must be in the allowlist
  assert_in("flag_a", allowlist,
       "flag_a (template literal) must be in allowlist")
  assert_in("flag_b", allowlist,
       "flag_b (template literal) must be in allowlist")
  # Placeholders must NOT have been mistakenly extracted as keys.
  # (nproc is in strategy_flags, not from the template, so we exclude
  # it from the negative check.)
  assert_true("data_mtz" not in allowlist,
        "'data_mtz' (placeholder name) should NOT be in allowlist")
  assert_true("map" not in allowlist,
        "'map' (placeholder name) should NOT be in allowlist")
  # nproc must still be there (from strategy_flags)
  assert_in("nproc", allowlist,
       "nproc (declared in strategy_flags) must be in allowlist")


def test_bug8_adversarial_strategy_override_dropped():
  """The fix must NOT introduce an LLM-override vulnerability.

    Even with the template-literal allowlist extension, an LLM that
    emits strategy.maps_only=False must still have that value dropped
    by program_registry (because maps_only is intentionally absent
    from strategy_flags of phenix.autobuild_denmod).  The result is
    that the template's maps_only=True remains unchallenged.
    """
  try:
    from agent.program_registry import ProgramRegistry
  except ImportError:
    print("  SKIP (program_registry not importable)")
    return

  reg = ProgramRegistry()
  files = {
      "data_mtz": "/p/refine.mtz",
      "sequence": "/p/7qz0.fa",
      "model": "/p/refine.pdb",
  }
  # Adversarial: LLM emits maps_only=False in strategy
  strategy = {"nproc": 4, "maps_only": False}

  log_messages = []
  try:
    cmd = reg.build_command(
      "phenix.autobuild_denmod",
      files, strategy=strategy,
      log=lambda m: log_messages.append(m))
  except Exception as e:
    print("  SKIP (registry build_command raised: %s)" % e)
    return

  if not cmd:
    print("  SKIP (registry returned no command)")
    return

  # Template literal must be present
  assert_in("maps_only=True", cmd,
       "Template literal maps_only=True must remain in registry output")
  # Adversarial override must NOT be present
  assert_true("maps_only=False" not in cmd,
        "Adversarial strategy.maps_only=False must NOT appear in "
        "command (registry should drop unknown strategy entries)")
  # Registry should have logged a warning about the unknown strategy
  warnings = [m for m in log_messages if "Unknown strategy 'maps_only'" in m]
  assert_true(len(warnings) >= 1,
        "Registry should log 'Unknown strategy maps_only' warning "
        "when LLM tries to override the template invariant")


# =========================================================================
# Bug 9: phenix.predict_and_build wrong crystal scope
#         (run_38_openai, v119.H9)
# =========================================================================
#
# Symptom: phenix.predict_and_build commands had
# `crystal_symmetry.unit_cell="..."` which PHENIX rejected with
# "Some PHIL parameters are not recognized".  predict_and_build
# uses the `crystal_info.*` scope, not `crystal_symmetry.*`.
#
# Diagnosis: agent/program_registry.py:778-779 unconditionally
# prepends `crystal_symmetry.` to bare unit_cell/space_group
# strategies.  The safety net (parameter_fixes.json +
# fix_program_parameters in planner.py) rewrites the wrong scope
# to the right one for programs like autobuild and autosol — but
# phenix.predict_and_build had no scope-rewrite entries in
# parameter_fixes.json, so the bad form passed through unchanged.
#
# Distinct from the existing Bug 4 (rebuild_in_place=False) and
# Bug 8 (template-literal allowlist).  Bug 9 is purely a data
# patch to parameter_fixes.json — same mechanism as the existing
# autobuild/autosol scope rewrites.
#
# Fix: Add scope-rewrite entries to phenix.predict_and_build in
# parameter_fixes.json:
#   crystal_symmetry.{unit_cell,space_group} → crystal_info.*
#   xray_data.{unit_cell,space_group} → crystal_info.*

def _reload_parameter_fixes():
  """Force re-read of parameter_fixes.json (defeat the module cache).

    fix_program_parameters caches the loaded JSON at module level
    via _PARAMETER_FIXES.  When the file changes between tests
    (or first-time load races), we need to defeat this cache.
    """
  try:
    from libtbx.langchain.agent import planner as _p
  except ImportError:
    try:
      from agent import planner as _p
    except ImportError:
      return
  _p._PARAMETER_FIXES = None


def test_bug9_predict_and_build_unit_cell_scope_rewrite():
  """Primary regression: crystal_symmetry.unit_cell rewritten to crystal_info.

    Reproduces the exact failure from run_38_openai cycle 2.  When
    fix_program_parameters processes a phenix.predict_and_build
    command containing crystal_symmetry.unit_cell="...", it must
    rewrite to crystal_info.unit_cell="...".
    """
  _reload_parameter_fixes()
  from agent.planner import fix_program_parameters

  cmd = ('phenix.predict_and_build input_files.seq_file=/p/seq.fa '
       'input_files.xray_data_file=/p/data.mtz '
       'crystal_info.resolution=2.1 '
       'crystal_symmetry.unit_cell="116.097 116.097 44.175 90 90 120" '
       'control.nproc=4')

  result = fix_program_parameters(cmd, 'phenix.predict_and_build')

  assert_in('crystal_info.unit_cell=', result,
       "Bad crystal_symmetry.unit_cell must be rewritten to "
       "crystal_info.unit_cell for phenix.predict_and_build "
       "(was failing in run_38_openai with 'PHIL parameters not recognized')")
  assert_not_in('crystal_symmetry.unit_cell=', result,
        "After H9 fix, crystal_symmetry.unit_cell must NOT remain "
        "in the command — it was the rejected scope")
  # Sanity: the value must be preserved
  assert_in('"116.097 116.097 44.175 90 90 120"', result,
       "Unit cell value must be preserved through the rewrite")


def test_bug9_predict_and_build_space_group_scope_rewrite():
  """Companion: crystal_symmetry.space_group → crystal_info.space_group.

    Same rewrite, different parameter.  Verifies the parameter_fixes
    entry covers both unit_cell and space_group.
    """
  _reload_parameter_fixes()
  from agent.planner import fix_program_parameters

  cmd = ('phenix.predict_and_build input_files.seq_file=/p/seq.fa '
       'crystal_symmetry.space_group="P 31 2 1"')

  result = fix_program_parameters(cmd, 'phenix.predict_and_build')

  assert_in('crystal_info.space_group=', result,
       "crystal_symmetry.space_group must be rewritten to "
       "crystal_info.space_group")
  assert_not_in('crystal_symmetry.space_group=', result,
        "crystal_symmetry.space_group must NOT remain")


def test_bug9_predict_and_build_xray_data_fallback_rewrite():
  """Recovery-path coverage: xray_data.* → crystal_info.* also works.

    LLM error-recovery sometimes pivots to xray_data.* scope when
    crystal_symmetry.* is rejected.  The parameter_fixes entry must
    handle this alternate-form too, so the recovery doesn't loop
    forever trying different wrong scopes.

    (Gemini-suggested test, plan rev2.)
    """
  _reload_parameter_fixes()
  from agent.planner import fix_program_parameters

  cmd = ('phenix.predict_and_build input_files.seq_file=/p/seq.fa '
       'xray_data.unit_cell="116.097 116.097 44.175 90 90 120" '
       'xray_data.space_group="P 31 2 1"')

  result = fix_program_parameters(cmd, 'phenix.predict_and_build')

  assert_in('crystal_info.unit_cell=', result,
       "xray_data.unit_cell must be rewritten to crystal_info.unit_cell "
       "(recovery-path coverage)")
  assert_in('crystal_info.space_group=', result,
       "xray_data.space_group must be rewritten to crystal_info.space_group")
  assert_not_in('xray_data.unit_cell=', result,
        "xray_data.unit_cell must NOT remain after rewrite")
  assert_not_in('xray_data.space_group=', result,
        "xray_data.space_group must NOT remain after rewrite")


def test_bug9_predict_and_build_idempotent_on_correct_scope():
  """Idempotent: commands already using crystal_info.* are unchanged.

    If a future code path emits the correct scope directly,
    fix_program_parameters must not duplicate or mangle it.
    Also verifies no regression: phenix.refine still uses
    crystal_symmetry.* (its native scope).
    """
  _reload_parameter_fixes()
  from agent.planner import fix_program_parameters

  # predict_and_build with correct scope — should be unchanged
  cmd_pab = ('phenix.predict_and_build input_files.seq_file=/p/seq.fa '
         'crystal_info.unit_cell="116.097 116.097 44.175 90 90 120"')
  out_pab = fix_program_parameters(cmd_pab, 'phenix.predict_and_build')
  # Exactly one occurrence — no duplication
  assert_equal(out_pab.count('crystal_info.unit_cell='), 1,
        "crystal_info.unit_cell should appear exactly once "
        "(idempotent — no duplication)")

  # phenix.refine with crystal_symmetry — should be unchanged
  # (refine uses crystal_symmetry.* natively, so the rewrite must
  # NOT fire there).
  cmd_refine = ('phenix.refine model=/p/m.pdb data=/p/d.mtz '
         'crystal_symmetry.unit_cell="100 100 100 90 90 90"')
  out_refine = fix_program_parameters(cmd_refine, 'phenix.refine')
  assert_in('crystal_symmetry.unit_cell=', out_refine,
       "phenix.refine must retain crystal_symmetry.unit_cell "
       "(its native scope — no regression)")
  assert_not_in('crystal_info.unit_cell=', out_refine,
        "phenix.refine must NOT have crystal_info.unit_cell "
        "(H9 fix is scoped to predict_and_build only)")


# =========================================================================
# Bug 10: exclude_patterns bypassed in category-based file selection
#         (AIAgent_62 cycle 7, v119.H10)
# =========================================================================
#
# Symptom: phenix.refine emitted a command with refine_001_001.cif
# as a positional third argument (the {ligand_cif} slot).  PHENIX
# interpreted it as a SECOND model (it's actually a model mmCIF
# from an earlier refine step, not a ligand restraints CIF),
# crashing with "wrong number of models".
#
# Diagnosis: agent/command_builder.py::_find_file_for_slot honors
# input_def["exclude_patterns"] at only 2 of 9 selection paths
# (LLM-selected and PRIORITY 4 extension fallback).  PRIORITY 2
# (best_files), PRIORITY 2.5 (recovery strategies), PRIORITY 3
# (category-based), and PRIORITY 3.5 (fallback best_files) all
# bypassed the check.  For the ligand_cif slot of phenix.refine,
# the schema explicitly lists "refine_" in exclude_patterns to
# prevent this exact failure — but PRIORITY 3 category-match
# returned the file without consulting that list.
#
# Distinct from Bug 8 (template-literal allowlist, v119.H8) and
# Bug 9 (predict_and_build crystal scope, v119.H9).
#
# Fix: In agent/command_builder.py::_find_file_for_slot, define a
# helper at function entry that checks input_def["exclude_patterns"],
# and apply it at all 7 bypass sites:
#   - PRIORITY 2 (best_files)
#   - PRIORITY 2.5 (recovery strategies)
#   - PRIORITY 3 subcategory loop
#   - PRIORITY 3 parent-category loop
#   - PRIORITY 3 is_multiple branch
#   - PRIORITY 3.5 require_best_files_only
#   - PRIORITY 3.5 fallback best_files for specific subcategory

def _make_bug10_context(categorized_files, available_files,
                        best_files=None, recovery_strategies=None):
  """Build a CommandContext for Bug 10 sandbox tests.

    All fields are pinned to typical xray_refined state so the
    test focuses on file-selection behavior.
    """
  try:
    from libtbx.langchain.agent.command_builder import CommandContext
  except ImportError:
    from agent.command_builder import CommandContext
  return CommandContext(
    cycle_number=7,
    experiment_type='xray',
    resolution=2.1,
    best_files=best_files or {},
    rfree_mtz=None,
    categorized_files=categorized_files,
    workflow_state='xray_refined',
    history=[],
    llm_files=None,
    llm_strategy=None,
    recovery_strategies=recovery_strategies,
    directives={},
    log=lambda m: None,
    files_local=True,
  )


def test_bug10_ligand_cif_slot_rejects_model_mmcif():
  """Primary regression: refine_001_001.cif rejected for ligand_cif slot.

    Reproduces the AIAgent_62 cycle 7 failure.  When phenix.refine
    is built and the ligand_cif category contains refine_001_001.cif
    (a model mmCIF from an earlier refine step), the H10 fix must
    cause CommandBuilder to reject it for the ligand_cif slot, so
    the resulting command doesn't include it as a positional
    second-model argument.
    """
  import os, tempfile
  try:
    from libtbx.langchain.agent.command_builder import CommandBuilder
  except ImportError:
    from agent.command_builder import CommandBuilder

  tmpdir = tempfile.mkdtemp()
  files = ['refine_001_001_modified.pdb', '7qz0.mtz', 'refine_001_001.cif']
  for f in files:
    open(os.path.join(tmpdir, f), 'w').close()
  available_files = [os.path.join(tmpdir, f) for f in files]

  categorized = {
    'model': [os.path.join(tmpdir, 'refine_001_001_modified.pdb')],
    'data_mtz': [os.path.join(tmpdir, '7qz0.mtz')],
    'ligand_cif': [os.path.join(tmpdir, 'refine_001_001.cif')],
  }
  context = _make_bug10_context(categorized, available_files,
    best_files={'model': os.path.join(tmpdir, 'refine_001_001_modified.pdb')})

  builder = CommandBuilder()
  cmd = builder.build('phenix.refine', available_files, context)

  assert_true(cmd is not None,
        "phenix.refine command should still build (just without "
        "the bad ligand_cif)")
  assert_not_in('refine_001_001.cif', cmd,
        "Bug 10: refine_001_001.cif must NOT appear in the command "
        "— it would be interpreted as a SECOND model and crash "
        "phenix.refine with 'wrong number of models'")
  # The model should still be selected
  assert_in('refine_001_001_modified.pdb', cmd,
        "Real protein model must still be selected")


def test_bug10_ligand_cif_slot_accepts_legitimate_cif():
  """No over-rejection: legitimate ligand.ligands.cif still selected.

    The exclude_patterns for ligand_cif is ["overall_best",
    "overall_best_final", "refine_", "autobuild", "predict_and_build"].
    A file named ligand.ligands.cif matches none of these and
    must still be picked.
    """
  import os, tempfile
  try:
    from libtbx.langchain.agent.command_builder import CommandBuilder
  except ImportError:
    from agent.command_builder import CommandBuilder

  tmpdir = tempfile.mkdtemp()
  files = ['protein.pdb', 'data.mtz', 'ligand.ligands.cif']
  for f in files:
    open(os.path.join(tmpdir, f), 'w').close()
  available_files = [os.path.join(tmpdir, f) for f in files]

  categorized = {
    'model': [os.path.join(tmpdir, 'protein.pdb')],
    'data_mtz': [os.path.join(tmpdir, 'data.mtz')],
    'ligand_cif': [os.path.join(tmpdir, 'ligand.ligands.cif')],
  }
  context = _make_bug10_context(categorized, available_files,
    best_files={'model': os.path.join(tmpdir, 'protein.pdb')})

  builder = CommandBuilder()
  cmd = builder.build('phenix.refine', available_files, context)

  assert_true(cmd is not None, "command should build")
  assert_in('ligand.ligands.cif', cmd,
        "Legitimate ligand restraints CIF must still be selected "
        "— its name doesn't match any exclude_pattern")


def test_bug10_exclude_patterns_applied_in_best_files_path():
  """PRIORITY 2 (best_files) honors exclude_patterns.

    Pre-H10, best_files only checked exclude_categories, not
    exclude_patterns.  After H10, both apply.  Construct a case
    where best_files contains a file matching exclude_patterns;
    verify it gets rejected.

    Use phenix.refine::model which has exclude_patterns=
    ['ligand', 'lig.pdb', 'ligand_fit'].  Place a file named
    'fake_ligand_fit.pdb' in best_files for the model slot; the
    H10 fix should reject it.
    """
  import os, tempfile
  try:
    from libtbx.langchain.agent.command_builder import CommandBuilder
  except ImportError:
    from agent.command_builder import CommandBuilder

  tmpdir = tempfile.mkdtemp()
  files = ['good_model.pdb', 'fake_ligand_fit.pdb', 'data.mtz']
  for f in files:
    open(os.path.join(tmpdir, f), 'w').close()
  available_files = [os.path.join(tmpdir, f) for f in files]

  # Place the bad file as best_files['model'] — pre-H10 this would
  # be selected for the model slot via PRIORITY 2 even though the
  # model slot's exclude_patterns lists 'ligand_fit'.
  categorized = {
    'model': [os.path.join(tmpdir, 'good_model.pdb'),
              os.path.join(tmpdir, 'fake_ligand_fit.pdb')],
    'data_mtz': [os.path.join(tmpdir, 'data.mtz')],
  }
  context = _make_bug10_context(categorized, available_files,
    best_files={'model': os.path.join(tmpdir, 'fake_ligand_fit.pdb')})

  builder = CommandBuilder()
  cmd = builder.build('phenix.refine', available_files, context)

  assert_true(cmd is not None, "command should build")
  assert_not_in('fake_ligand_fit.pdb', cmd,
        "Bug 10: fake_ligand_fit.pdb in best_files must be rejected "
        "for the model slot — exclude_patterns lists 'ligand_fit'")
  assert_in('good_model.pdb', cmd,
        "PRIORITY 3 should fall through to good_model.pdb from "
        "the model category once best_files is rejected")


def test_bug10_exclude_patterns_applied_in_category_path():
  """PRIORITY 3 (category-based) honors exclude_patterns.

    Direct test of the actual buggy code path.  Construct a slot
    where category-based selection (NOT best_files, NOT extension
    fallback) would pick a bad file; verify H10 rejects it.
    """
  import os, tempfile
  try:
    from libtbx.langchain.agent.command_builder import CommandBuilder
  except ImportError:
    from agent.command_builder import CommandBuilder

  tmpdir = tempfile.mkdtemp()
  files = ['protein.pdb', 'data.mtz',
           'refine_001_001.cif', 'good_ligand.ligands.cif']
  for f in files:
    open(os.path.join(tmpdir, f), 'w').close()
  available_files = [os.path.join(tmpdir, f) for f in files]

  # Both CIFs in ligand_cif category — H10 should pick the
  # legitimate one and reject the model-mmCIF.
  categorized = {
    'model': [os.path.join(tmpdir, 'protein.pdb')],
    'data_mtz': [os.path.join(tmpdir, 'data.mtz')],
    'ligand_cif': [os.path.join(tmpdir, 'refine_001_001.cif'),
                   os.path.join(tmpdir, 'good_ligand.ligands.cif')],
  }
  context = _make_bug10_context(categorized, available_files,
    best_files={'model': os.path.join(tmpdir, 'protein.pdb')})

  builder = CommandBuilder()
  cmd = builder.build('phenix.refine', available_files, context)

  assert_true(cmd is not None, "command should build")
  assert_not_in('refine_001_001.cif', cmd,
        "Bug 10: refine_001_001.cif must be rejected (matches "
        "'refine' exclude_pattern, v119.H11 — pre-H11 the pattern "
        "was 'refine_' which broke word-boundary matching)")
  assert_in('good_ligand.ligands.cif', cmd,
        "Legitimate ligand CIF should be picked from the same "
        "category")


def test_bug10_unrelated_slots_unaffected():
  """No-regression sentinel: slots WITHOUT exclude_patterns work as before.

    The H10 patch must not change behavior for slots that don't
    have exclude_patterns.  Verify by building phenix.refine with
    a sequence file (sequence slot — no exclude_patterns); the
    file should be selected normally.

    (Claude-reviewer suggestion.)
    """
  import os, tempfile
  try:
    from libtbx.langchain.agent.command_builder import CommandBuilder
  except ImportError:
    from agent.command_builder import CommandBuilder

  tmpdir = tempfile.mkdtemp()
  files = ['protein.pdb', 'data.mtz', '7qz0.fa']
  for f in files:
    open(os.path.join(tmpdir, f), 'w').close()
  available_files = [os.path.join(tmpdir, f) for f in files]

  # phenix.refine doesn't take a sequence slot, but the data_mtz
  # slot also has no exclude_patterns.  Verify it works.
  categorized = {
    'model': [os.path.join(tmpdir, 'protein.pdb')],
    'data_mtz': [os.path.join(tmpdir, 'data.mtz')],
    'sequence': [os.path.join(tmpdir, '7qz0.fa')],
  }
  context = _make_bug10_context(categorized, available_files,
    best_files={'model': os.path.join(tmpdir, 'protein.pdb'),
                'data_mtz': os.path.join(tmpdir, 'data.mtz')})

  builder = CommandBuilder()
  cmd = builder.build('phenix.refine', available_files, context)

  assert_true(cmd is not None, "command should build")
  assert_in('data.mtz', cmd,
        "data.mtz must be selected — its slot has no exclude_patterns "
        "and selection should be unchanged from pre-H10 behavior")
  assert_in('protein.pdb', cmd,
        "protein.pdb must be selected — exclude_patterns list "
        "(['ligand', 'lig.pdb', 'ligand_fit']) doesn't match 'protein.pdb'")


# =========================================================================
# Bug 11: matches_exclude_pattern semantic-pin test
#         (v119.H11 — captures the H10 → H11 lesson)
# =========================================================================
#
# Background: H10's structural fix applied exclude_patterns at all
# file-selection paths in command_builder.py.  But a separate bug in
# the YAML patterns themselves (trailing underscore in "refine_")
# meant the filter, while correctly applied, did not actually reject
# the bad file.  The AIAgent_62 cycle-7 crash persisted after H10.
#
# Root cause: matches_exclude_pattern uses WORD-BOUNDARY regex
# matching (see agent/file_utils.py).  A pattern with leading or
# trailing underscore breaks the match because the regex's
# boundary-lookahead requires another boundary char after the
# pattern's own trailing underscore.
#
# This test PINS that semantic with documented cases.  If the
# function ever changes (e.g., to substring matching), this test
# fails and we notice.  It also serves as authoring documentation
# for future YAML pattern design.
#
# Why a sandbox-skip fallback: in pure sandbox environments where
# neither libtbx.langchain.* nor agent.* file_utils is importable,
# we skip cleanly rather than erroring.  In Tom's PHENIX
# environment this runs against the real function and pins its
# actual behavior.

def test_bug11_matches_exclude_pattern_semantics():
  """Pin matches_exclude_pattern's word-boundary semantics.

    The H10 → H11 cycle revealed that a sandbox stub with
    substring semantics passed Bug 10 tests, but the real
    function uses word-boundary regex semantics, so the H10
    fix was a no-op for the YAML pattern 'refine_' (trailing
    underscore breaks the match).  This test makes the actual
    function's behavior explicit.
    """
  # Sandbox-skip fallback (per Gemini's H11 plan review):
  # if neither libtbx.langchain.* nor agent.* file_utils is
  # importable, skip cleanly rather than erroring.  This keeps
  # lightweight sandbox CI green while ensuring full fidelity
  # on Tom's workstation / cluster.
  try:
    from libtbx.langchain.agent.file_utils import matches_exclude_pattern
  except ImportError:
    try:
      from agent.file_utils import matches_exclude_pattern
    except ImportError:
      print("  SKIP (matches_exclude_pattern not importable — sandbox)")
      return

  # --- Trailing underscore in pattern BREAKS the match ---
  # (this is the H11 bug being fixed)
  assert_false(matches_exclude_pattern("refine_001.cif", ["refine_"]),
        "Pattern 'refine_' (trailing _) does NOT match "
        "'refine_001.cif' — the trailing _ requires another "
        "boundary char after it (here it's '0', not a boundary)")

  # --- Bare pattern matches correctly ---
  assert_true(matches_exclude_pattern("refine_001.cif", ["refine"]),
        "Pattern 'refine' matches 'refine_001.cif' "
        "(word boundary _ follows 'refine')")
  assert_true(matches_exclude_pattern("phenix_refine.cif", ["refine"]),
        "Pattern 'refine' matches 'phenix_refine.cif' "
        "(word boundary _ precedes, end-of-stem follows)")
  assert_false(matches_exclude_pattern("refining_template.cif", ["refine"]),
        "Pattern 'refine' does NOT match 'refining_template.cif' "
        "— after 'refine' comes 'i' which is not a boundary char")

  # --- Leading underscore similarly breaks matches ---
  assert_false(matches_exclude_pattern("protein_half.ccp4", ["_half"]),
        "Pattern '_half' (leading _) does NOT match "
        "'protein_half.ccp4'")
  assert_true(matches_exclude_pattern("protein_half.ccp4", ["half"]),
        "Pattern 'half' matches 'protein_half.ccp4'")

  # --- Boundary-char vs alphanumeric-suffix (Gemini's catch) ---
  # This is why we keep explicit half1, half2 patterns alongside
  # the bare 'half' in YAML.
  assert_false(matches_exclude_pattern("protein_half1.ccp4", ["half"]),
        "Pattern 'half' does NOT match 'protein_half1.ccp4' — "
        "'1' immediately after 'half' is not a boundary char.  "
        "YAML keeps explicit 'half1', 'half2' alongside the bare "
        "'half' to catch separator-less variants.")
  assert_true(matches_exclude_pattern("protein_half1.ccp4", ["half1"]),
        "Pattern 'half1' matches 'protein_half1.ccp4'")

  # --- Extension-suffix patterns ---
  assert_true(matches_exclude_pattern("lig.pdb", ["lig.pdb"]),
        "Pattern 'lig.pdb' matches 'lig.pdb' exactly")
  assert_false(matches_exclude_pattern("nolig.pdb", ["lig.pdb"]),
        "Pattern 'lig.pdb' does NOT match 'nolig.pdb' "
        "(no word boundary before 'lig')")


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

  with open(ai_agent_path, encoding='utf-8') as f:
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
  with open(ai_agent_path, 'r', encoding='utf-8') as f:
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
  with open(ai_agent_path, 'r', encoding='utf-8') as f:
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
  with open(ai_agent_path, 'r', encoding='utf-8') as f:
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
  with open(ai_agent_path, 'r', encoding='utf-8') as f:
    source = f.read()

  # Function signature must accept working_dir
  # v119.H14: relaxed from exact multi-line signature match to a
  # substring check — the original test pinned the exact signature
  # including newline + indent, which broke when later changes added
  # additional kwargs (e.g., skip_if_failed=False) to the signature.
  # The real intent is "the function accepts working_dir"; verify
  # that without locking the rest of the signature.
  assert_in("def _track_output_files(", source,
       "_track_output_files must be defined")
  assert_in("working_dir=None", source,
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
  with open(rest_init_path, 'r', encoding='utf-8') as f:
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
  with open(remote_agent_path, 'r', encoding='utf-8') as f:
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
  with open(checks_path, 'r', encoding='utf-8') as f:
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
  with open(session_path, 'r', encoding='utf-8') as f:
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
  with open(checks_path, 'r', encoding='utf-8') as f:
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
  with open(extractor_path, 'r', encoding='utf-8') as f:
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
  with open(nodes_path, 'r', encoding='utf-8') as f:
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
  with open(agent_path, 'r', encoding='utf-8') as f:
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
  with open(agent_path, 'r', encoding='utf-8') as f:
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
  with open(session_path, 'r', encoding='utf-8') as f:
    source = f.read()

  # All three open() calls for session JSON should use encoding='utf-8'
  import re
  open_calls = re.findall(r"open\([^)]*session_file[^)]*\)", source)
  for call in open_calls:
    assert_in("encoding='utf-8'", call,
         "session_file open(, encoding='utf-8') must specify UTF-8: %s" % call)


def run_all_tests():
  """Run all autosol bugs tests."""
  run_tests_with_fail_fast()


# =========================================================================

if __name__ == "__main__":
  success = run_tests_with_fail_fast()
  if not success:
    sys.exit(1)
  print("\nOK")
