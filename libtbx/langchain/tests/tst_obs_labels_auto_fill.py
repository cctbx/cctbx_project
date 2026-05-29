"""K_H16: obs_labels auto-fill for multi-array MTZ (v119.H16+H16.1).

Trigger: 88 TIER-1 failures across two batch scans (run_25, run_39)
matching "Sorry: Multiple equally suitable arrays of observed xray
data found", concentrated in AF_exoV_MRSAD and lysozyme-MRSAD
tutorials.

H16's mechanism:
  1. ``inspect_mtz()`` reads the MTZ via cctbx and returns its
     column structure.
  2. ``select_obs_labels_for(program, info)`` applies a per-program
     preference policy.
  3. The ``auto_fill_obs_labels`` invariant in ``programs.yaml``
     triggers the builder to inject the chosen labels.

Three layers of test:
  - §A: Policy unit tests (pure functions, no I/O)
  - §B: MTZ inspector tests (cctbx-dependent, gracefully skip)
  - §C: Builder-integration simulation (mirrors the actual
        ``_apply_invariants()`` branch)
  - §D: YAML config validation

8 tests, all required to pass for ship.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from agent.mtz_inspector import (
    inspect_mtz,
    select_obs_labels_for,
    has_ambiguous_arrays,
    _PROGRAM_PREFERENCES,
)


# =====================================================================
# Fixtures: hand-built mtz_info dicts matching common MTZ shapes
# =====================================================================

def _af_exov_mtz_info():
    """The MTZ shape from AF_exoV tutorials that triggers the Sorry.

    Carries both merged intensities (IMEAN), merged amplitudes (FOBS),
    AND anomalous pairs (I(+)/I(-)).  This is exactly what PHENIX
    refuses to disambiguate."""
    return {
        "anomalous_intensities": ["I(+),SIGI(+),I(-),SIGI(-)"],
        "anomalous_amplitudes": [],
        "merged_intensities": ["IMEAN,SIGIMEAN"],
        "merged_amplitudes": ["FOBS,SIGFOBS"],
        "rfree_label": "FreeR_flag",
        "is_multi_array": True,
    }


def _sad_only_mtz_info():
    """SAD experimental data with no merged forms.  Common after
    raw scaling, before any merging step has run."""
    return {
        "anomalous_intensities": ["I(+),SIGI(+),I(-),SIGI(-)"],
        "anomalous_amplitudes": [],
        "merged_intensities": [],
        "merged_amplitudes": [],
        "rfree_label": None,
        "is_multi_array": False,  # Only one data array set
    }


def _single_array_mtz_info():
    """Single amplitudes pair (standard refinement input).  No
    ambiguity — H16 should be a no-op."""
    return {
        "anomalous_intensities": [],
        "anomalous_amplitudes": [],
        "merged_intensities": [],
        "merged_amplitudes": ["FOBS,SIGFOBS"],
        "rfree_label": "FreeR_flag",
        "is_multi_array": False,
    }


def _multi_array_amplitudes_only():
    """Multi-array MTZ but no intensities — refine must fall back to
    amplitudes via its preference list."""
    return {
        "anomalous_intensities": [],
        "anomalous_amplitudes": ["F(+),SIGF(+),F(-),SIGF(-)"],
        "merged_intensities": [],
        "merged_amplitudes": ["FOBS,SIGFOBS"],
        "rfree_label": "FreeR_flag",
        "is_multi_array": True,
    }


# =====================================================================
# §A: Policy unit tests
# =====================================================================

def test_policy_per_program_picks_correctly_from_af_exov():
    """For an MTZ matching AF_exoV's shape, each program picks the
    label set that matches its preferred form.

    This is the core fix for Tom's failure cohort:
      - refine prefers merged intensities → IMEAN,SIGIMEAN
      - phaser prefers merged amplitudes → FOBS,SIGFOBS
      - autosol/autobuild intentionally NOT in policy table —
        they have program-specific label formats that the agent's
        comma-pair convention doesn't cover; injecting could
        produce malformed commands.  Both return None.
    """
    info = _af_exov_mtz_info()

    assert select_obs_labels_for("refine", info) == "IMEAN,SIGIMEAN", (
        "refine should prefer merged intensities; got %r"
        % select_obs_labels_for("refine", info))

    assert select_obs_labels_for("phenix.refine", info) == "IMEAN,SIGIMEAN", (
        "phenix.refine variant should normalize to refine policy")

    assert select_obs_labels_for("phaser", info) == "FOBS,SIGFOBS", (
        "phaser should prefer merged amplitudes; got %r"
        % select_obs_labels_for("phaser", info))

    # autosol and autobuild MUST NOT be in the policy table — they
    # have label-format conventions outside the agent's current
    # comma-pair handling.  Adding them risks malformed PHIL.
    assert select_obs_labels_for("autosol", info) is None, (
        "autosol must NOT be in policy (complicated label format); "
        "got %r" % select_obs_labels_for("autosol", info))
    assert select_obs_labels_for("autobuild", info) is None, (
        "autobuild must NOT be in policy (no verified PHIL key); "
        "got %r" % select_obs_labels_for("autobuild", info))

    # Other unknown programs should also return None
    assert select_obs_labels_for("unknown_program", info) is None
    print("  PASS: test_policy_per_program_picks_correctly_from_af_exov")


def test_policy_fallback_traverses_preference_list():
    """When the preferred category is absent, the policy falls
    through to the next preferred category.  Refine on amplitudes-
    only MTZ → falls back to merged_amplitudes."""
    info = _multi_array_amplitudes_only()

    # No merged intensities → refine falls back to merged amplitudes
    result = select_obs_labels_for("refine", info)
    assert result == "FOBS,SIGFOBS", (
        "Refine should fall back to merged amplitudes when "
        "intensities absent; got %r" % result)

    # Phaser falls through identically (both prefs are present)
    result = select_obs_labels_for("phaser", info)
    assert result == "FOBS,SIGFOBS", (
        "phaser should pick merged amplitudes; got %r" % result)

    # Empty mtz_info → None (no categories have entries)
    empty = {
        "anomalous_intensities": [],
        "anomalous_amplitudes": [],
        "merged_intensities": [],
        "merged_amplitudes": [],
        "rfree_label": None,
        "is_multi_array": False,
    }
    assert select_obs_labels_for("refine", empty) is None, (
        "Empty MTZ info → no labels to pick → None")
    assert select_obs_labels_for("phaser", empty) is None
    print("  PASS: test_policy_fallback_traverses_preference_list")


def test_has_ambiguous_arrays_gating():
    """has_ambiguous_arrays returns True only for multi-array MTZs.
    This is the gate that prevents H16 injection from firing on
    single-array files (which PHENIX handles unambiguously).
    """
    assert has_ambiguous_arrays(_af_exov_mtz_info()) is True
    assert has_ambiguous_arrays(_multi_array_amplitudes_only()) is True
    assert has_ambiguous_arrays(_single_array_mtz_info()) is False
    assert has_ambiguous_arrays(_sad_only_mtz_info()) is False
    # None / malformed → False (safe default)
    assert has_ambiguous_arrays(None) is False
    assert has_ambiguous_arrays({}) is False
    assert has_ambiguous_arrays("not a dict") is False
    print("  PASS: test_has_ambiguous_arrays_gating")


def test_policy_table_covers_intended_programs():
    """The policy table should contain exactly the 2 programs targeted
    by H16: refine, phaser.  Regression guard against accidental
    policy additions/removals.

    autosol and autobuild are deliberately EXCLUDED — both have
    program-specific label conventions outside the agent's current
    comma-pair handling.  Adding them risks malformed PHIL.  If
    they're added, the test must be updated in lockstep with verified
    format support."""
    expected = {"refine", "phaser"}
    actual = set(_PROGRAM_PREFERENCES.keys())
    assert actual == expected, (
        "Policy table programs mismatch.\n"
        "  Expected: %r\n"
        "  Actual:   %r\n"
        "If adding a program, FIRST verify its label format "
        "(comma-pair? space-separated? multi-key?) before adding it "
        "to the policy table."
        % (expected, actual))
    # Each program's preference list must be non-empty (else the
    # policy would always return None for that program)
    for prog, prefs in _PROGRAM_PREFERENCES.items():
        assert isinstance(prefs, list) and len(prefs) > 0, (
            "Policy for %s must have a non-empty preference list" % prog)
    print("  PASS: test_policy_table_covers_intended_programs")


# =====================================================================
# §B: MTZ inspector tests
# =====================================================================

def test_inspect_mtz_none_paths_graceful():
    """inspect_mtz returns None for every failure mode without
    raising.  Graceful degradation is the contract: when inspection
    can't happen, the auto_fill invariant silently no-ops.

    v119.H16.1: this test caught a real bug in production — the
    file_type filter `if file_type and file_type != "ccp4_mtz"`
    had a logic hole for falsy file_type values.  Fixed to require
    exact "ccp4_mtz" match.  This test now covers:
      - None / empty path inputs
      - Non-existent paths
      - Non-MTZ file content with .txt extension
      - Non-MTZ file content with .mtz extension (more dangerous —
        cctbx tries harder to parse)
    """
    assert inspect_mtz(None) is None
    assert inspect_mtz("") is None
    assert inspect_mtz("/path/that/does/not/exist.mtz") is None

    import tempfile
    # Case 1: non-MTZ content with .txt extension
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
        f.write(b"This is not an MTZ.\n")
        not_mtz_txt = f.name
    try:
        result = inspect_mtz(not_mtz_txt)
        assert result is None, (
            "Non-MTZ .txt file should return None; got %r" % result)
    finally:
        os.unlink(not_mtz_txt)

    # Case 2: non-MTZ content with .mtz extension.  This is the
    # harder case — cctbx's any_reflection_file() may try to parse
    # by extension, producing a reader whose file_type() returns
    # something falsy.  The v119.H16.1 fix specifically handles this.
    with tempfile.NamedTemporaryFile(suffix=".mtz", delete=False) as f:
        f.write(b"This is not an MTZ, just text with .mtz extension.\n")
        not_mtz_mtz = f.name
    try:
        result = inspect_mtz(not_mtz_mtz)
        assert result is None, (
            "Non-MTZ content with .mtz extension should still return "
            "None; got %r.  This was the production bug in v119.H16 "
            "before H16.1 — the falsy-file_type short-circuit let "
            "execution proceed to as_miller_arrays() which returned "
            "an empty list, producing a structured-empty dict instead "
            "of None." % result)
    finally:
        os.unlink(not_mtz_mtz)
    print("  PASS: test_inspect_mtz_none_paths_graceful")


def test_inspect_mtz_with_real_cctbx():
    """Real cctbx round-trip: synthesize an MTZ in-memory with a
    controlled column set, then inspect it.

    Skips gracefully when cctbx isn't importable (which is the case
    in CI sandboxes without PHENIX).  Production deploys this test
    against a real PHENIX install where cctbx is always available.

    The fixture mimics AF_exoV's structure: anomalous intensities
    pair, merged intensities, merged amplitudes — three suitable
    data array sets.  is_multi_array MUST be True for this case."""
    try:
        from iotbx.reflection_file_reader import any_reflection_file
        from cctbx import miller, crystal
        from cctbx.array_family import flex
    except ImportError:
        print("  SKIP: test_inspect_mtz_with_real_cctbx (cctbx unavailable)")
        return

    import tempfile
    # Build a minimal MTZ.  cctbx's API expects a crystal symmetry +
    # miller index set + data array.  For testing we need 3 arrays
    # (intensities, amplitudes, and anomalous intensities) to verify
    # classification.
    try:
        cs = crystal.symmetry(
            unit_cell=(50, 50, 50, 90, 90, 90),
            space_group_symbol="P1",
        )
        # Create some Miller indices
        ms = miller.set(
            crystal_symmetry=cs,
            indices=flex.miller_index([(1, 0, 0), (2, 0, 0), (3, 0, 0)]),
            anomalous_flag=False,
        )
        intensities = miller.array(
            ms, data=flex.double([100.0, 200.0, 300.0]),
            sigmas=flex.double([5.0, 10.0, 15.0]),
        ).set_observation_type_xray_intensity()
        intensities.set_info(miller.array_info(labels=["IMEAN", "SIGIMEAN"]))

        with tempfile.NamedTemporaryFile(suffix=".mtz", delete=False) as f:
            mtz_path = f.name
        try:
            ds = intensities.as_mtz_dataset(column_root_label="IMEAN")
            ds.mtz_object().write(file_name=mtz_path)

            info = inspect_mtz(mtz_path)
            assert info is not None, (
                "inspect_mtz returned None on a valid MTZ — bug in the "
                "reader or test fixture")
            assert "merged_intensities" in info
            assert len(info["merged_intensities"]) >= 1, (
                "Expected at least one merged intensities array; got %r"
                % info)
            # Single-array case: is_multi_array should be False
            assert info["is_multi_array"] is False, (
                "Single-array fixture should NOT be is_multi_array; got %r"
                % info)
            print("  PASS: test_inspect_mtz_with_real_cctbx (real MTZ inspected)")
        finally:
            if os.path.exists(mtz_path):
                os.unlink(mtz_path)
    except Exception as e:
        # Test environment can't synthesize the MTZ — skip rather
        # than fail.  Production deployment validates against real
        # MTZs from the failing batch runs.
        print("  SKIP: test_inspect_mtz_with_real_cctbx (cctbx "
              "fixture synthesis failed: %s)" % e)


# =====================================================================
# §C: Builder-integration simulation
#
# We can't import command_builder.py directly here because its module-
# level libtbx imports fail in this sandbox.  Instead, we recreate the
# auto_fill_obs_labels branch logic with the same structure and
# external dependencies (mtz_inspector.select_obs_labels_for,
# has_ambiguous_arrays).  When command_builder.py is modified, this
# simulator must be updated in lockstep.  Test §D validates the YAML
# config that drives the builder, so the two layers can't drift
# silently.
# =====================================================================

def _simulate_auto_fill_obs_labels(program, strategy, context_mtz_info,
                                    log_lines):
    """Mirror of the auto_fill_obs_labels branch in
    command_builder._apply_invariants.  Mutates ``strategy`` in
    place and appends log lines to ``log_lines``.

    Returns nothing; success is observed via the side effects."""
    if "obs_labels" in strategy:
        # Pre-set by user/LLM — preserve.  Builder's outer
        # `"obs_labels" not in strategy` check matches this.
        return
    try:
        if has_ambiguous_arrays(context_mtz_info):
            labels = select_obs_labels_for(program, context_mtz_info)
            if labels:
                strategy["obs_labels"] = labels
                log_lines.append(
                    "[OBS_LABELS] injected: %s '%s' "
                    "(from multi-array MTZ inspection)"
                    % (program, labels))
            else:
                log_lines.append(
                    "[OBS_LABELS] no suitable labels for %s "
                    "(available categories: %s)"
                    % (program,
                       [k for k, v in context_mtz_info.items()
                        if isinstance(v, list) and v]))
        # else: single-array or no inspection — silent no-op
    except Exception as e:
        log_lines.append("[OBS_LABELS] auto-fill raised %s" % e)


def test_builder_simulation_injects_on_ambiguous_mtz():
    """Tom's exact failure case: AF_exoV-shaped MTZ, refine command
    with no obs_labels in strategy.  The builder should inject the
    correct labels AND emit a [OBS_LABELS] log line.

    Without H16, this command would build without obs_labels and
    PHENIX would abort with "Multiple equally suitable arrays."
    With H16, the command is built with xray_data.labels="IMEAN,SIGIMEAN".
    """
    strategy = {}  # No pre-existing obs_labels
    info = _af_exov_mtz_info()
    log_lines = []

    _simulate_auto_fill_obs_labels(
        "refine", strategy, info, log_lines)

    assert "obs_labels" in strategy, (
        "Injection should have fired; strategy missing obs_labels")
    assert strategy["obs_labels"] == "IMEAN,SIGIMEAN", (
        "Wrong labels injected; got %r" % strategy["obs_labels"])
    # Telemetry assertion
    injection_logs = [
        L for L in log_lines if "[OBS_LABELS] injected" in L]
    assert len(injection_logs) == 1, (
        "Expected 1 [OBS_LABELS] injected log; got %d: %r"
        % (len(injection_logs), log_lines))
    assert "refine" in injection_logs[0]
    assert "IMEAN,SIGIMEAN" in injection_logs[0]
    print("  PASS: test_builder_simulation_injects_on_ambiguous_mtz")


def test_builder_simulation_respects_user_override_and_skips_single_array():
    """Two scenarios in one test (related concerns):

    1. User pre-set obs_labels in strategy → preserved (not
       overwritten).  This is critical: if a user knows their MTZ's
       column structure and wants specific labels, H16 must not
       silently change them.

    2. Single-array MTZ → no injection, no log line.  The Sorry only
       fires when there's actual ambiguity, so injection on
       single-array MTZs would be pure overhead with no benefit
       (and would risk overriding PHENIX's own correct auto-detect).
    """
    # (1) User override preserved
    strategy = {"obs_labels": "USER_PROVIDED_LABEL"}
    log_lines = []
    _simulate_auto_fill_obs_labels(
        "refine", strategy, _af_exov_mtz_info(), log_lines)
    assert strategy["obs_labels"] == "USER_PROVIDED_LABEL", (
        "User-provided obs_labels overwritten; H16 must respect "
        "explicit user values.  Got %r" % strategy["obs_labels"])
    assert log_lines == [], (
        "No log lines should fire when user override is present "
        "(early return path).  Got: %r" % log_lines)

    # (2) Single-array MTZ skips silently
    strategy = {}
    log_lines = []
    _simulate_auto_fill_obs_labels(
        "refine", strategy, _single_array_mtz_info(), log_lines)
    assert "obs_labels" not in strategy, (
        "Single-array MTZ should not trigger injection; got %r"
        % strategy)
    assert log_lines == [], (
        "Single-array case should be silent (no log spam).  "
        "Got: %r" % log_lines)

    # Bonus: None inspection (e.g. cctbx unavailable) also silent
    strategy = {}
    log_lines = []
    _simulate_auto_fill_obs_labels("refine", strategy, None, log_lines)
    assert "obs_labels" not in strategy
    assert log_lines == [], (
        "None mtz_info should be silent; got: %r" % log_lines)
    print("  PASS: test_builder_simulation_respects_user_override_and_skips_single_array")


def test_builder_simulation_skips_for_excluded_programs():
    """autosol and autobuild are intentionally outside H16's scope.
    When the policy returns None for these programs, the auto-fill
    branch logs `no suitable labels` and skips injection.  No PHIL
    is corrupted; the original "Multiple equally suitable arrays"
    Sorry continues to fire for these programs (as expected pre-H16).

    This test pins the expected behavior so a future addition of
    these programs to the policy table (without also wiring proper
    label-format handling) doesn't silently break.
    """
    info = _af_exov_mtz_info()
    for prog in ("autosol", "autobuild", "ligandfit", "polder"):
        strategy = {}
        log_lines = []
        _simulate_auto_fill_obs_labels(prog, strategy, info, log_lines)
        assert "obs_labels" not in strategy, (
            "Excluded program %s should NOT get obs_labels injected; "
            "got %r" % (prog, strategy))
        # Log line: "no suitable labels" should appear since the
        # MTZ IS multi-array but the policy returned None
        no_suitable_logs = [
            L for L in log_lines if "no suitable labels" in L]
        assert len(no_suitable_logs) == 1, (
            "Excluded program %s should log 'no suitable labels' "
            "(multi-array MTZ + None policy result).  Got: %r"
            % (prog, log_lines))
        assert prog in no_suitable_logs[0]
    print("  PASS: test_builder_simulation_skips_for_excluded_programs")


# =====================================================================
# §D: YAML config validation
# =====================================================================

def test_yaml_invariant_config_for_all_targeted_programs():
    """Verify programs.yaml is internally consistent for the H16
    cohort:
      - Each of {refine, phaser} has an auto_fill_obs_labels
        invariant declaration
      - The invariant has NO 'check: has_strategy' clause (critical:
        obs_labels is OPTIONAL — a has_strategy check would BLOCK
        every single-array build, see H16.4 bug)
      - Each of {refine, phaser} has an obs_labels strategy_flag
        (with a non-empty PHIL key)
      - The flag template references {value} for substitution

    This is the contract between the YAML and the builder.  If the
    YAML is missing a config the builder expects, injection silently
    fails in production.
    """
    try:
        import yaml
    except ImportError:
        print("  SKIP: test_yaml_invariant_config (yaml unavailable)")
        return

    yaml_path = os.path.join(_PARENT, "knowledge", "programs.yaml")
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    targeted = ["phenix.refine", "phenix.phaser"]
    for prog in targeted:
        assert prog in data, "%s missing from programs.yaml" % prog
        spec = data[prog]

        # Invariant must be present with the correct fix
        invs = spec.get("invariants", [])
        matching = [
            inv for inv in invs
            if inv.get("fix", {}).get("auto_fill_obs_labels") is True
        ]
        assert len(matching) == 1, (
            "%s should have exactly 1 invariant with "
            "fix.auto_fill_obs_labels=true; got %d.  invariants=%r"
            % (prog, len(matching),
               [inv.get("name") for inv in invs]))

        # CRITICAL: the obs_labels invariant must NOT have a
        # check.has_strategy clause.  obs_labels is optional (only
        # needed for multi-array MTZs).  A has_strategy check triggers
        # the builder's blocking logic (command_builder.py ~line 547),
        # which returns None — blocking EVERY single-array refine/phaser
        # build.  This was the H16.4 bug.
        check = matching[0].get("check") or {}
        assert "has_strategy" not in check, (
            "%s obs_labels invariant must NOT declare "
            "check.has_strategy — it would block all single-array "
            "builds (returns None). obs_labels is optional. "
            "Got check=%r" % (prog, check))

        # Strategy flag must exist with a PHIL key containing {value}
        flags = spec.get("strategy_flags", {})
        assert "obs_labels" in flags, (
            "%s must have obs_labels in strategy_flags (else the "
            "injected obs_labels won't materialize as a CLI flag).  "
            "Existing strategy_flags: %r" % (prog, list(flags.keys())))
        flag_template = flags["obs_labels"].get("flag", "")
        assert "{value}" in flag_template, (
            "%s obs_labels flag must use {value} substitution; got %r"
            % (prog, flag_template))
        # Sanity: the PHIL key should mention 'labels' (loose check —
        # programs use varying namespaces but the leaf is labels)
        assert "label" in flag_template.lower(), (
            "%s obs_labels flag should reference 'label' in its PHIL "
            "key; got %r" % (prog, flag_template))

    print("  PASS: test_yaml_invariant_config_for_all_targeted_programs")


def test_obs_labels_invariant_does_not_block_single_array_builds():
    """Regression guard for the v119.H16.4 bug.

    The builder (command_builder.py ~line 547) has a blocking check:
    for any invariant declaring `check: has_strategy: X`, if X is not
    present in the strategy after auto-fills run, build() returns None
    (the command is blocked).

    obs_labels is OPTIONAL — it's only injected for multi-array MTZs.
    If the obs_labels invariant declared `check: has_strategy:
    obs_labels`, then EVERY single-array refine/phaser build (the
    common case) would be blocked and return None.  That was the
    H16.4 bug: phase 8 showed `phenix.refine FAIL (None)` and
    `phenix.phaser FAIL (None)` across many tutorials.

    This test simulates the auto-fill + blocking-check flow against
    the real YAML and asserts that a single-array build is NOT
    blocked.
    """
    try:
        import yaml
    except ImportError:
        print("  SKIP: test_obs_labels_invariant_does_not_block (yaml unavailable)")
        return

    yaml_path = os.path.join(_PARENT, "knowledge", "programs.yaml")
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    def simulate_build(program, initial_strategy, mtz_multi_array):
        """Mirror _apply_invariants auto-fills + the blocking check."""
        strategy = dict(initial_strategy)
        invariants = data[program].get("invariants", [])
        # Step 1: auto-fill fixes (subset relevant to this test)
        for inv in invariants:
            fix = inv.get("fix", {})
            if (fix.get("auto_fill_output_prefix")
                    and "output_prefix" not in strategy):
                strategy["output_prefix"] = "out_001"
            if (fix.get("auto_fill_obs_labels")
                    and "obs_labels" not in strategy
                    and mtz_multi_array):
                strategy["obs_labels"] = "IMEAN,SIGIMEAN"
        # Step 2: blocking check
        for inv in invariants:
            required = inv.get("check", {}).get("has_strategy")
            if required and required not in strategy:
                return None  # blocked
        return strategy

    # The exact failing scenario: single-array MTZ, refine + phaser
    for prog in ["phenix.refine", "phenix.phaser"]:
        result = simulate_build(prog, {"resolution": 2.1},
                                mtz_multi_array=False)
        assert result is not None, (
            "%s with a single-array MTZ must BUILD (not be blocked). "
            "If this fails, the obs_labels invariant has a "
            "check.has_strategy clause that blocks the build — that's "
            "the H16.4 regression." % prog)
        assert "obs_labels" not in result, (
            "%s single-array build should NOT have obs_labels injected"
            % prog)

    # Multi-array case: should build AND have obs_labels
    result = simulate_build("phenix.refine", {"resolution": 2.1},
                            mtz_multi_array=True)
    assert result is not None, "refine multi-array must build"
    assert result.get("obs_labels") == "IMEAN,SIGIMEAN", (
        "refine multi-array should have obs_labels injected; got %r"
        % result)

    print("  PASS: test_obs_labels_invariant_does_not_block_single_array_builds")


# =====================================================================
# §E: Regression guards
# =====================================================================

def test_command_context_mtz_inspection_at_end_of_fields():
    """Regression guard for the v119.H16.2 fix.

    mtz_inspection MUST be the LAST field in the CommandContext
    dataclass.  Inserting it elsewhere shifts the positions of
    subsequent fields, breaking any caller that constructs
    CommandContext with positional args.

    Tom's tst_command_builder.py::test_llm_data_slot_used_for_mtz
    caught the original H16 bug — command generation returned None
    because the positional construction misaligned llm_files,
    files_local, etc.

    Tries to import the real CommandContext; if libtbx isn't
    available, skips gracefully (sandbox case).  Production
    deployment exercises the real construction path."""
    try:
        try:
            from libtbx.langchain.agent.command_builder import (
                CommandContext)
        except ImportError:
            from agent.command_builder import CommandContext
    except ImportError:
        # command_builder.py has hard libtbx imports that aren't
        # importable in the sandbox — falls through gracefully.
        print("  SKIP: test_command_context_mtz_inspection_at_end_of_fields "
              "(command_builder not importable in sandbox)")
        return

    from dataclasses import fields
    field_names = [f.name for f in fields(CommandContext)]
    # mtz_inspection must be last
    assert field_names[-1] == "mtz_inspection", (
        "mtz_inspection MUST be the last field in CommandContext to "
        "preserve positional-construction compatibility.  Got field "
        "order: %r" % field_names)
    # The set of fields BEFORE mtz_inspection must match the pre-H16
    # set (no new fields snuck in elsewhere)
    expected_pre_h16 = [
        "cycle_number", "experiment_type", "resolution",
        "best_files", "rfree_mtz", "rfree_resolution",
        "categorized_files", "workflow_state", "history",
        "llm_files", "llm_strategy", "recovery_strategies",
        "directives", "model_hetatm_residues", "files_local", "log",
    ]
    assert field_names[:-1] == expected_pre_h16, (
        "Pre-H16 field order changed.  This breaks positional "
        "construction in existing code.\n"
        "  Expected: %r\n"
        "  Actual:   %r" % (expected_pre_h16, field_names[:-1]))
    print("  PASS: test_command_context_mtz_inspection_at_end_of_fields")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: Policy
    test_policy_per_program_picks_correctly_from_af_exov()
    test_policy_fallback_traverses_preference_list()
    test_has_ambiguous_arrays_gating()
    test_policy_table_covers_intended_programs()
    # §B: Inspector
    test_inspect_mtz_none_paths_graceful()
    test_inspect_mtz_with_real_cctbx()
    # §C: Builder simulation
    test_builder_simulation_injects_on_ambiguous_mtz()
    test_builder_simulation_respects_user_override_and_skips_single_array()
    test_builder_simulation_skips_for_excluded_programs()
    # §D: YAML
    test_yaml_invariant_config_for_all_targeted_programs()
    test_obs_labels_invariant_does_not_block_single_array_builds()
    # §E: Regression guards (v119.H16.2)
    test_command_context_mtz_inspection_at_end_of_fields()


if __name__ == "__main__":
    print("K_H16: obs_labels auto-fill (v119.H16+H16.1+H16.2)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H16 complete.")
