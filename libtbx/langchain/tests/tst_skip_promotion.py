"""K_H2.1: skip_programs promotion regression tests.

v119.H2.1.  Validates that program_settings[X].skip=truthy gets
promoted to workflow_preferences.skip_programs at the top of
validate_directives.

The promotion fixes a pre-existing gap surfaced by the
skip_programs scenario in tests/llm/tst_directive_extraction.py:
the LLM correctly identifies "don't run mtriage" as
program_settings["phenix.mtriage"].skip=true, but the downstream
consumer reads skip_programs from workflow_preferences.

Test organization (matches plan rev 3 §5):
  §5.1: Promotion mechanics                (4 tests)
  §5.2: Side-effect correctness            (4 tests)
  §5.3: Idempotency / merging              (1 test)
  §5.4: Defensive bail                     (1 test)
  §5.5: End-to-end through validate_directives  (1 test)

Total: 11 tests.

Sandbox-import note: tests that need validate_directives use
_try_import_directive_extractor with graceful skip when the
module's transitive deps are unavailable.  Pure-helper tests
can also use the graceful-skip pattern since the helper lives
in the same module.
"""
from __future__ import absolute_import, division, print_function


def _try_import_directive_extractor():
    """Try to import directive_extractor; return (module, error) tuple.

    Returns (module, None) on success, (None, error_msg) on failure.
    Mirrors the pattern K_H1 uses for its retired-helper tests
    that need imports requiring heavy transitive deps.
    """
    try:
        from libtbx.langchain.agent import directive_extractor
        return directive_extractor, None
    except ImportError:
        pass
    try:
        from agent import directive_extractor
        return directive_extractor, None
    except ImportError as e:
        return None, str(e)


# =====================================================================
# §5.1: Promotion mechanics (4 tests)
# =====================================================================

def test_promotes_program_settings_skip_true_to_skip_programs():
    """Canonical case: program_settings[X].skip=True promotes to skip_programs."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    wf = directives.get("workflow_preferences", {})
    assert wf.get("skip_programs") == ["phenix.mtriage"], (
        "Expected skip_programs == ['phenix.mtriage'], got %r"
        % wf.get("skip_programs"))
    print("  PASS: test_promotes_program_settings_skip_true_to_skip_programs")


def test_promotes_string_truthy_values():
    """All _SKIP_TRUE_VALUES entries promote (Gemini guardrail A regression)."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # Exercise every value in the truthiness list
    for truthy_value in de._SKIP_TRUE_VALUES:
        directives = {
            "program_settings": {
                "phenix.mtriage": {"skip": truthy_value},
            },
        }
        de._promote_skip_settings_to_skip_programs(directives)
        wf = directives.get("workflow_preferences", {})
        assert wf.get("skip_programs") == ["phenix.mtriage"], (
            "Truthy value %r should promote; got skip_programs=%r"
            % (truthy_value, wf.get("skip_programs")))
    print("  PASS: test_promotes_string_truthy_values "
          "(%d truthy values tested)" % len(de._SKIP_TRUE_VALUES))


def test_does_not_promote_falsy_values():
    """Falsy / unrecognized values do NOT promote."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    falsy_cases = [False, "false", "False", "FALSE",
                   "no", "No", "NO",
                   0, "0", "",
                   None, "skip_maybe", "abc"]
    for falsy_value in falsy_cases:
        directives = {
            "program_settings": {
                "phenix.mtriage": {"skip": falsy_value},
            },
        }
        de._promote_skip_settings_to_skip_programs(directives)
        wf = directives.get("workflow_preferences", {})
        skip_progs = wf.get("skip_programs", [])
        assert "phenix.mtriage" not in skip_progs, (
            "Falsy value %r should NOT promote; got skip_programs=%r"
            % (falsy_value, skip_progs))
        # And the skip key should still be there in program_settings
        # (i.e. we didn't pop it for a falsy value).
        ps_skip = directives.get("program_settings", {}).get(
            "phenix.mtriage", {}).get("skip", "MISSING")
        assert ps_skip == falsy_value, (
            "Falsy skip value should remain in program_settings; "
            "value=%r preserved=%r" % (falsy_value, ps_skip))
    print("  PASS: test_does_not_promote_falsy_values "
          "(%d falsy values tested)" % len(falsy_cases))


def test_does_not_promote_unrelated_keys():
    """program_settings without a 'skip' key is unchanged."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.refine": {"resolution": 2.0, "cycles": 5},
            "phenix.autosol": {"atom_type": "Se"},
        },
    }
    import copy
    before = copy.deepcopy(directives)
    de._promote_skip_settings_to_skip_programs(directives)
    assert directives == before, (
        "Helper should be a no-op when no skip keys are present.\n"
        "  before: %r\n  after:  %r" % (before, directives))
    print("  PASS: test_does_not_promote_unrelated_keys")


# =====================================================================
# §5.2: Side-effect correctness (4 tests)
# =====================================================================

def test_removes_skip_flag_from_program_settings_after_promotion():
    """After promotion, the program's skip key is gone (Gemini guardrail B)."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True, "resolution": 3.0},
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    # Skip key is gone from program_settings
    ps = directives.get("program_settings", {}).get("phenix.mtriage", {})
    assert "skip" not in ps, (
        "skip key should be popped from program_settings; "
        "remaining settings=%r" % ps)
    # But the OTHER setting is preserved (covered more in
    # test_preserves_other_settings_when_skip_promoted)
    assert ps.get("resolution") == 3.0
    print("  PASS: test_removes_skip_flag_from_program_settings_after_promotion")


def test_no_unknown_setting_log_fires_for_promoted_skip():
    """Calling validate_directives must not emit "Unknown setting skip".

    Locks in the promotion-before-validation ordering.  If a future
    refactor moves the promotion call to AFTER the per-setting
    validation loop, this test fails immediately.
    """
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
    }
    captured_log_messages = []
    de.validate_directives(directives, log=captured_log_messages.append)
    offenders = [m for m in captured_log_messages
                 if "Unknown setting skip" in m]
    assert not offenders, (
        "validate_directives emitted Unknown-setting log line(s) "
        "for promoted skip value(s); promotion ordering is broken.\n"
        "Offending messages:\n  %s" % "\n  ".join(offenders))
    # Should see the promotion log message
    promotion_msgs = [m for m in captured_log_messages
                      if "Promoted phenix.mtriage" in m]
    assert promotion_msgs, (
        "Expected at least one 'Promoted phenix.mtriage' log line; "
        "captured messages:\n  %s"
        % "\n  ".join(captured_log_messages))
    print("  PASS: test_no_unknown_setting_log_fires_for_promoted_skip")


def test_removes_program_entry_if_only_setting_was_skip():
    """Empty program_settings dict after promotion -> remove entry.

    Gemini guardrail C regression.
    """
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    ps = directives.get("program_settings", {})
    assert "phenix.mtriage" not in ps, (
        "program_settings should not contain phenix.mtriage after "
        "promotion (was empty dict); got %r" % ps)
    print("  PASS: test_removes_program_entry_if_only_setting_was_skip")


def test_preserves_other_settings_when_skip_promoted():
    """Other settings under the same program survive the promotion."""
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True, "resolution": 2.5},
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    # program_settings still contains phenix.mtriage with resolution
    ps = directives.get("program_settings", {}).get("phenix.mtriage")
    assert ps == {"resolution": 2.5}, (
        "Expected program_settings['phenix.mtriage'] == "
        "{'resolution': 2.5}, got %r" % ps)
    # And the program is in skip_programs
    wf = directives.get("workflow_preferences", {})
    assert wf.get("skip_programs") == ["phenix.mtriage"]
    print("  PASS: test_preserves_other_settings_when_skip_promoted")


# =====================================================================
# §5.3: Idempotency / merging (1 test)
# =====================================================================

def test_promotion_merges_with_existing_skip_programs():
    """De-duplication + merging (Gemini guardrail D regression).

    Three sub-cases:
      (a) Merge with existing different entry.
      (b) De-duplicate when same program already in skip_programs.
      (c) Program is still removed from program_settings even when
          its name was already in skip_programs.
    """
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # (a) Merge with existing different entry
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
        "workflow_preferences": {
            "skip_programs": ["phenix.xtriage"],
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    skip_progs = directives["workflow_preferences"]["skip_programs"]
    assert "phenix.xtriage" in skip_progs and "phenix.mtriage" in skip_progs, (
        "(a) Expected both programs in skip_programs, got %r" % skip_progs)
    assert len(skip_progs) == 2, (
        "(a) Expected exactly 2 entries, got %d: %r"
        % (len(skip_progs), skip_progs))

    # (b) De-duplicate when same program already present
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
        "workflow_preferences": {
            "skip_programs": ["phenix.mtriage"],  # already there
        },
    }
    de._promote_skip_settings_to_skip_programs(directives)
    skip_progs = directives["workflow_preferences"]["skip_programs"]
    assert skip_progs == ["phenix.mtriage"], (
        "(b) Expected de-duplicated ['phenix.mtriage'], got %r"
        % skip_progs)
    assert skip_progs.count("phenix.mtriage") == 1, (
        "(b) Duplicate added: skip_programs=%r" % skip_progs)

    # (c) Program still removed from program_settings even when
    # already in skip_programs
    assert "phenix.mtriage" not in directives.get("program_settings", {}), (
        "(c) Empty program_settings entry should be removed even "
        "when name was already in skip_programs; got %r"
        % directives.get("program_settings"))

    print("  PASS: test_promotion_merges_with_existing_skip_programs "
          "(a+b+c sub-cases)")


# =====================================================================
# §5.4: Defensive bail (1 test)
# =====================================================================

def test_does_not_promote_when_workflow_preferences_malformed():
    """Malformed workflow_preferences / skip_programs -> leave input alone.

    Rev 3 finding 4 regression.  If a future refactor "improves"
    the bail to clobber workflow_preferences = {} this test fails.
    """
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # Sub-case 1: workflow_preferences is not a dict
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
        "workflow_preferences": "not a dict",
    }
    captured_log = []
    de._promote_skip_settings_to_skip_programs(
        directives, log=captured_log.append)
    # skip key still present
    assert directives["program_settings"]["phenix.mtriage"]["skip"] is True, (
        "(1) skip key should remain when workflow_preferences malformed; "
        "got %r" % directives["program_settings"])
    # workflow_preferences unchanged
    assert directives["workflow_preferences"] == "not a dict", (
        "(1) workflow_preferences was clobbered; got %r"
        % directives["workflow_preferences"])
    # log message explains the bail
    assert any("skip-promotion skipped" in m for m in captured_log), (
        "(1) Expected 'skip-promotion skipped' log line, got %r"
        % captured_log)

    # Sub-case 2: skip_programs is not a list
    directives = {
        "program_settings": {
            "phenix.mtriage": {"skip": True},
        },
        "workflow_preferences": {
            "skip_programs": "not a list",
        },
    }
    captured_log = []
    de._promote_skip_settings_to_skip_programs(
        directives, log=captured_log.append)
    # skip key still present
    assert directives["program_settings"]["phenix.mtriage"]["skip"] is True, (
        "(2) skip key should remain when skip_programs malformed; "
        "got %r" % directives["program_settings"])
    # skip_programs unchanged
    assert directives["workflow_preferences"]["skip_programs"] == "not a list", (
        "(2) skip_programs was clobbered; got %r"
        % directives["workflow_preferences"])
    # log message explains the bail
    assert any("skip-promotion skipped" in m for m in captured_log), (
        "(2) Expected 'skip-promotion skipped' log line, got %r"
        % captured_log)

    print("  PASS: test_does_not_promote_when_workflow_preferences_malformed "
          "(both bail sub-cases)")


# =====================================================================
# §5.5: End-to-end through validate_directives (1 test)
# =====================================================================

def test_skip_promotion_via_validate_directives():
    """The 'would the failing scenario pass?' test.

    Calls validate_directives (not the helper directly) on the
    exact shape the LLM was producing in the failing
    skip_programs scenario.  This pins the fix to the reported
    failure.
    """
    de, err = _try_import_directive_extractor()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # Shape from the failing scenario's Parsed JSON
    directives = {
        "program_settings": {
            "default": {
                "wavelength": 1.0,
                "resolution": 2.0,
                "space_group": "P 1",
            },
            "phenix.mtriage": {
                "skip": True,
            },
        },
        "stop_conditions": {},
        "file_preferences": {
            "model": "model.pdb",
            "data_mtz": "data.mtz",
        },
        "workflow_preferences": {
            "use_experimental_phasing": False,
            "use_molecular_replacement": False,
        },
        "constraints": [
            "Proceed with refinement and validation.",
        ],
    }
    result = de.validate_directives(directives)
    wf = result.get("workflow_preferences", {})
    skip_progs = wf.get("skip_programs", [])
    assert skip_progs == ["phenix.mtriage"], (
        "End-to-end: validate_directives should produce "
        "skip_programs=['phenix.mtriage']; got %r" % skip_progs)
    # phenix.mtriage should not appear in program_settings of the
    # result (it had only the skip key, which is now promoted).
    ps = result.get("program_settings", {})
    assert "phenix.mtriage" not in ps, (
        "End-to-end: phenix.mtriage should be removed from "
        "program_settings (was empty after promotion); got %r" % ps)
    print("  PASS: test_skip_promotion_via_validate_directives")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §5.1: Promotion mechanics (4)
    test_promotes_program_settings_skip_true_to_skip_programs()
    test_promotes_string_truthy_values()
    test_does_not_promote_falsy_values()
    test_does_not_promote_unrelated_keys()
    # §5.2: Side-effect correctness (4)
    test_removes_skip_flag_from_program_settings_after_promotion()
    test_no_unknown_setting_log_fires_for_promoted_skip()
    test_removes_program_entry_if_only_setting_was_skip()
    test_preserves_other_settings_when_skip_promoted()
    # §5.3: Idempotency / merging (1)
    test_promotion_merges_with_existing_skip_programs()
    # §5.4: Defensive bail (1)
    test_does_not_promote_when_workflow_preferences_malformed()
    # §5.5: End-to-end (1)
    test_skip_promotion_via_validate_directives()


if __name__ == "__main__":
    run_all_tests()
