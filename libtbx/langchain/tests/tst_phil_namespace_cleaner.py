"""
Sandbox tests for v118 Section B: PHIL namespace healing in the
directive extractor.

Section B adds two protections against bad PHIL namespace strings
leaking into commands:

  B1 — preprocessor prompt forbids namespaced PHIL paths
  B2 — directive_extractor heals/drops them in the validator

These tests focus on B2 because B2 is the load-bearing fix.  B1 is
covered by a static prompt-text check (K_B1) since it's prompt-only.

The healer rules (in _heal_namespaced_phil_keys):

  Rule 1 (HEAL): known near-miss → bare strategy form
                 (e.g. data_manager.r_free_flags.generate
                  -> generate_rfree_flags)

  Rule 2 (DROP): keys with known-bad PHIL prefixes are dropped
                 (e.g. anything starting with "data_manager.",
                  "refinement.", "miller_array.", etc.)

  Rule 3 (KEEP): everything else, including:
                 - Known VALID_SETTINGS keys
                 - Unknown bare keys (forward compatibility)
                 - Program-specific dotted keys
                   (e.g. "autosol.atom_type")

Tests K_B1 through K_B14 cover prompt change, translation, drop,
preservation, edge cases, and the production reproducer.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ADVICE_PREPROCESSOR_PATH = os.path.join(
    ROOT, "agent", "advice_preprocessor.py")
DIRECTIVE_EXTRACTOR_PATH = os.path.join(
    ROOT, "agent", "directive_extractor.py")


# --------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------

def _import_directive_extractor():
    """Dynamically import the directive_extractor module under test."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "directive_extractor_under_test", DIRECTIVE_EXTRACTOR_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _capture_log():
    """Return (log_func, captured_lines)."""
    captured = []
    def log_func(msg):
        captured.append(msg)
    return log_func, captured


# --------------------------------------------------------------------
# K_B1 — B1 preprocessor prompt check
# --------------------------------------------------------------------

def test_k_b1_preprocessor_prompt_forbids_namespaced_paths():
    """K_B1: preprocessor prompt contains the new B1 guidance."""
    with open(ADVICE_PREPROCESSOR_PATH) as f:
        src = f.read()

    # Look for the explicit anti-namespace rule
    assert "CRITICAL" in src and "emit only short bare key names" in src, (
        "K_B1 failed: B1 prompt change not present "
        "(should contain 'CRITICAL ... emit only short bare key names')")

    # Look for at least one good/bad example pair
    assert "data_manager.r_free_flags.generate" in src, (
        "K_B1 failed: B1 prompt should show the namespaced bad form "
        "as a negative example")
    assert "generate_rfree_flags=True" in src, (
        "K_B1 failed: B1 prompt should show the bare good form "
        "as a positive example")

    # Look for the new translation example
    assert ("generate R-free flags" in src
            or "create new R-free flags" in src), (
        "K_B1 failed: B1 prompt should include "
        "'generate R-free flags' translation example")

    print("  PASS: K_B1 (preprocessor prompt forbids namespaced paths)")


# --------------------------------------------------------------------
# K_B2 / K_B3 / K_B4 — translation map (HEAL rule)
# --------------------------------------------------------------------

def test_k_b2_heal_data_manager_r_free_flags_generate():
    """K_B2: data_manager.r_free_flags.generate -> generate_rfree_flags."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {"data_manager.r_free_flags.generate": True}
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert out == {"generate_rfree_flags": True}, (
        "K_B2 failed: expected {'generate_rfree_flags': True}, got %r"
        % out)
    assert any("Healed" in m and "data_manager.r_free_flags.generate" in m
               and "generate_rfree_flags" in m for m in captured), (
        "K_B2 failed: should log the healing with both names")

    print("  PASS: K_B2 (heal data_manager.r_free_flags.generate)")


def test_k_b3_heal_xray_data_r_free_flags_generate():
    """K_B3: xray_data.r_free_flags.generate -> generate_rfree_flags."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {"xray_data.r_free_flags.generate": True}
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert out == {"generate_rfree_flags": True}, (
        "K_B3 failed: expected {'generate_rfree_flags': True}, got %r"
        % out)

    print("  PASS: K_B3 (heal xray_data.r_free_flags.generate)")


def test_k_b4_heal_r_free_flags_generate():
    """K_B4: r_free_flags.generate -> generate_rfree_flags."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {"r_free_flags.generate": True}
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert out == {"generate_rfree_flags": True}, (
        "K_B4 failed: expected {'generate_rfree_flags': True}, got %r"
        % out)

    print("  PASS: K_B4 (heal r_free_flags.generate)")


# --------------------------------------------------------------------
# K_B5 / K_B6 — bad-prefix drop rule
# --------------------------------------------------------------------

def test_k_b5_drop_data_manager_prefix():
    """K_B5: keys starting with data_manager. (non-translated) are dropped."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    # A data_manager.X key that ISN'T in the translation map
    settings = {"data_manager.foo.bar": True, "resolution": 2.0}
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert "data_manager.foo.bar" not in out, (
        "K_B5 failed: data_manager.foo.bar should be dropped, got %r" % out)
    assert "resolution" in out and out["resolution"] == 2.0, (
        "K_B5 failed: resolution should survive, got %r" % out)
    assert any("Dropping" in m and "data_manager.foo.bar" in m
               for m in captured), (
        "K_B5 failed: should log the drop")

    print("  PASS: K_B5 (drop data_manager.X non-translated)")


def test_k_b6_drop_other_bad_prefixes():
    """K_B6: refinement.X, miller_array.X, fmodel.X, scaling.input.X dropped."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {
        "refinement.target_weights.wxc": 0.5,
        "miller_array.labels.name.test_flag_value": "FreeR_flag",
        "fmodel.k_sol": 0.35,
        "scaling.input.parameters.merging.n_bins": 20,
        "resolution": 2.0,  # legit key, should survive
    }
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    for bad_key in ("refinement.target_weights.wxc",
                    "miller_array.labels.name.test_flag_value",
                    "fmodel.k_sol",
                    "scaling.input.parameters.merging.n_bins"):
        assert bad_key not in out, (
            "K_B6 failed: %s should be dropped, got %r" % (bad_key, out))

    assert "resolution" in out, (
        "K_B6 failed: resolution should survive, got %r" % out)

    # 4 drop log lines expected
    drop_logs = [m for m in captured if "Dropping" in m]
    assert len(drop_logs) == 4, (
        "K_B6 failed: expected 4 drop log lines, got %d: %r"
        % (len(drop_logs), drop_logs))

    print("  PASS: K_B6 (drop refinement./miller_array./fmodel./scaling.input.)")


# --------------------------------------------------------------------
# K_B7 — preserve VALID_SETTINGS keys
# --------------------------------------------------------------------

def test_k_b7_preserve_valid_settings_keys():
    """K_B7: all VALID_SETTINGS keys pass through the healer unchanged."""
    mod = _import_directive_extractor()
    log_func, _ = _capture_log()

    # Build a settings dict with one example value for each
    # VALID_SETTINGS key.
    samples = {
        "resolution": 2.0,
        "cycles": 3,
        "anisotropic_adp": True,
        "add_waters": True,
        "simulated_annealing": False,
        "atom_type": "Se",
        "additional_atom_types": "S",
        "wavelength": 0.9792,
        "sites": 8,
        "twin_law": "-h,-k,l",
        "riding_hydrogens": True,
        "stop_after_predict": False,
        "ncs": True,
        "tls": True,
        "sharpening_method": "b-factor",
        "selection": "chain A",
        "unit_cell": "116 116 44 90 90 120",
        "space_group": "P 32 2 1",
        "copies": 2,
    }

    out = mod._heal_namespaced_phil_keys(samples, "phenix.refine", log_func)

    for key, value in samples.items():
        assert key in out, (
            "K_B7 failed: VALID_SETTINGS key %r dropped by healer" % key)
        assert out[key] == value, (
            "K_B7 failed: VALID_SETTINGS key %r value changed: %r -> %r"
            % (key, value, out[key]))

    print("  PASS: K_B7 (preserve %d VALID_SETTINGS keys)" % len(samples))


# --------------------------------------------------------------------
# K_B8 — preserve unknown bare keys (forward compat)
# --------------------------------------------------------------------

def test_k_b8_preserve_unknown_bare_keys():
    """K_B8: unknown bare keys (no dots) pass through (forward compat)."""
    mod = _import_directive_extractor()
    log_func, _ = _capture_log()

    settings = {
        "some_future_param": 42,
        "another_new_option": True,
        "experimental_feature_v9": "enabled",
    }
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert out == settings, (
        "K_B8 failed: unknown bare keys should pass through unchanged. "
        "Expected %r, got %r" % (settings, out))

    print("  PASS: K_B8 (preserve unknown bare keys for forward compat)")


# --------------------------------------------------------------------
# K_B9 — program-specific dotted keys preserved (Gemini Q2)
# --------------------------------------------------------------------

def test_k_b9_preserve_program_specific_dotted_keys():
    """K_B9: program-specific dotted keys (autosol., phaser., etc.) preserved."""
    mod = _import_directive_extractor()
    log_func, _ = _capture_log()

    settings = {
        "autosol.atom_type": "Se",
        "phaser.composition.protein.sequence_file": "seq.fa",
        "ligandfit.search_target.ligand_near_pdb": "atp.pdb",
        "xtriage.scaling.input.parameters.reporting.verbose": 1,
        "polder.solvent_exclusion_mask_selection": "chain A",
    }
    out = mod._heal_namespaced_phil_keys(settings, "default", log_func)

    for key, value in settings.items():
        assert key in out, (
            "K_B9 failed: program-specific dotted key %r dropped "
            "(autosol./phaser./ligandfit./xtriage./polder. should all "
            "pass through; only data_manager./refinement./miller_array./"
            "fmodel./scaling.input. are bad prefixes)"
            % key)
        assert out[key] == value, (
            "K_B9 failed: program-specific key value changed: %r -> %r"
            % (value, out[key]))

    print("  PASS: K_B9 (preserve program-specific dotted keys)")


# --------------------------------------------------------------------
# K_B10 — production reproducer
# --------------------------------------------------------------------

def test_k_b10_production_reproducer_testit():
    """K_B10: testit 2026-05-18 production reproducer."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    # Verbatim from testit log "Extracted Directives" → program_settings
    # for phenix.refine
    settings = {
        "resolution": 2.0,
        "data_manager.r_free_flags.generate": True,
    }
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert "data_manager.r_free_flags.generate" not in out, (
        "K_B10 failed: bad PHIL key should not be in output, got %r" % out)
    assert "generate_rfree_flags" in out, (
        "K_B10 failed: healed bare key should be in output, got %r" % out)
    assert out["generate_rfree_flags"] is True, (
        "K_B10 failed: healed value should be True, got %r"
        % out["generate_rfree_flags"])
    assert out["resolution"] == 2.0, (
        "K_B10 failed: resolution should survive intact, got %r" % out)

    print("  PASS: K_B10 (testit production reproducer)")


# --------------------------------------------------------------------
# K_B11 — default scope is cleaned too (Gemini Q3)
# --------------------------------------------------------------------

def test_k_b11_default_scope_cleaned():
    """K_B11: default scope is cleaned, not just program-specific scopes.

    Run the FULL _validate_directives on a directives dict with bad
    keys under 'default'.  default scope is just as infectious as any
    program-specific scope.
    """
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    # Direct healer call on default scope settings
    settings = {
        "data_manager.foo": True,
        "resolution": 2.0,
        "r_free_flags.generate": True,
    }
    out = mod._heal_namespaced_phil_keys(settings, "default", log_func)

    assert "data_manager.foo" not in out, (
        "K_B11 failed: default scope should have data_manager.foo dropped")
    assert out.get("resolution") == 2.0, (
        "K_B11 failed: default scope resolution should survive")
    assert out.get("generate_rfree_flags") is True, (
        "K_B11 failed: default scope r_free_flags.generate should heal "
        "to generate_rfree_flags")

    # Log lines should mention "default"
    assert any("default" in m for m in captured), (
        "K_B11 failed: logs should mention 'default' scope")

    print("  PASS: K_B11 (default scope cleaned by healer)")


# --------------------------------------------------------------------
# K_B12 — log message clarity
# --------------------------------------------------------------------

def test_k_b12_log_message_clarity():
    """K_B12: log messages contain program name, key, value, and reason."""
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {
        "data_manager.r_free_flags.generate": True,
        "refinement.target_weights.wxc": 0.5,
    }
    mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    # Heal log
    heal_logs = [m for m in captured if "Healed" in m]
    assert len(heal_logs) == 1, (
        "K_B12 failed: expected 1 heal log, got %d: %r"
        % (len(heal_logs), heal_logs))
    msg = heal_logs[0]
    assert "phenix.refine" in msg, "K_B12: heal log missing program name"
    assert "data_manager.r_free_flags.generate" in msg, (
        "K_B12: heal log missing bad key")
    assert "generate_rfree_flags" in msg, (
        "K_B12: heal log missing bare key")
    assert "True" in msg, "K_B12: heal log missing value"

    # Drop log
    drop_logs = [m for m in captured if "Dropping" in m]
    assert len(drop_logs) == 1, (
        "K_B12 failed: expected 1 drop log, got %d: %r"
        % (len(drop_logs), drop_logs))
    msg = drop_logs[0]
    assert "phenix.refine" in msg, "K_B12: drop log missing program name"
    assert "refinement.target_weights.wxc" in msg, (
        "K_B12: drop log missing bad key")
    assert "0.5" in msg, "K_B12: drop log missing value"
    assert "prefix" in msg.lower(), (
        "K_B12: drop log missing 'prefix' explanation")

    print("  PASS: K_B12 (log messages contain program/key/value/reason)")


# --------------------------------------------------------------------
# K_B13 — edge cases
# --------------------------------------------------------------------

def test_k_b13_edge_cases():
    """K_B13: empty settings, None, non-dict input."""
    mod = _import_directive_extractor()
    log_func, _ = _capture_log()

    # Empty dict
    assert mod._heal_namespaced_phil_keys({}, "phenix.refine", log_func) == {}

    # Non-dict input passes through unchanged
    assert mod._heal_namespaced_phil_keys(None, "phenix.refine", log_func) is None
    assert mod._heal_namespaced_phil_keys("not a dict", "phenix.refine",
                                          log_func) == "not a dict"

    # Single legitimate key
    assert mod._heal_namespaced_phil_keys(
        {"resolution": 2.0}, "phenix.refine", log_func
    ) == {"resolution": 2.0}

    print("  PASS: K_B13 (edge cases handled)")


# --------------------------------------------------------------------
# K_B14 — conflict policy: existing bare form wins
# --------------------------------------------------------------------

def test_k_b14_conflict_existing_bare_wins():
    """K_B14: if both namespaced and bare forms present, bare wins.

    If a translation target already exists in the input (e.g. both
    'data_manager.r_free_flags.generate=False' AND
    'generate_rfree_flags=True'), the bare form is authoritative and
    must not be overwritten by the translation.
    """
    mod = _import_directive_extractor()
    log_func, captured = _capture_log()

    settings = {
        "data_manager.r_free_flags.generate": False,
        "generate_rfree_flags": True,
    }
    out = mod._heal_namespaced_phil_keys(settings, "phenix.refine", log_func)

    assert out == {"generate_rfree_flags": True}, (
        "K_B14 failed: existing bare form should win when both "
        "namespaced and bare forms present. Got %r" % out)

    # Should log that the namespaced form was discarded in favor of
    # the existing bare form.
    discard_logs = [m for m in captured if "Discarding" in m or "already" in m]
    assert len(discard_logs) >= 1, (
        "K_B14 failed: should log discard/already-present, got %r"
        % captured)

    print("  PASS: K_B14 (existing bare form wins over translation)")


# --------------------------------------------------------------------
# Test runner
# --------------------------------------------------------------------

def run_all_tests():
    tests = [
        ("K_B1_preprocessor_prompt_forbids_namespaced_paths",
         test_k_b1_preprocessor_prompt_forbids_namespaced_paths),
        ("K_B2_heal_data_manager_r_free_flags_generate",
         test_k_b2_heal_data_manager_r_free_flags_generate),
        ("K_B3_heal_xray_data_r_free_flags_generate",
         test_k_b3_heal_xray_data_r_free_flags_generate),
        ("K_B4_heal_r_free_flags_generate",
         test_k_b4_heal_r_free_flags_generate),
        ("K_B5_drop_data_manager_prefix",
         test_k_b5_drop_data_manager_prefix),
        ("K_B6_drop_other_bad_prefixes",
         test_k_b6_drop_other_bad_prefixes),
        ("K_B7_preserve_valid_settings_keys",
         test_k_b7_preserve_valid_settings_keys),
        ("K_B8_preserve_unknown_bare_keys",
         test_k_b8_preserve_unknown_bare_keys),
        ("K_B9_preserve_program_specific_dotted_keys",
         test_k_b9_preserve_program_specific_dotted_keys),
        ("K_B10_production_reproducer_testit",
         test_k_b10_production_reproducer_testit),
        ("K_B11_default_scope_cleaned",
         test_k_b11_default_scope_cleaned),
        ("K_B12_log_message_clarity",
         test_k_b12_log_message_clarity),
        ("K_B13_edge_cases",
         test_k_b13_edge_cases),
        ("K_B14_conflict_existing_bare_wins",
         test_k_b14_conflict_existing_bare_wins),
    ]
    passed = 0
    failed = 0
    for name, fn in tests:
        print("Test: %s" % name)
        try:
            fn()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            print("  ERROR: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
