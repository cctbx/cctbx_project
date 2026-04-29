"""
Phase 9: LLM Perturbation Tests.

Tests that the deterministic pipeline is resilient to LLM output
variations without requiring live LLM calls or recorded sessions.

Items from plan:
  4. Filename perturbation — non-existent file → fallback to category
  5. Program perturbation — wrong program → still buildable
  6. Parameter injection — invalid PHIL → validator strips it
  7. Response truncation — truncated JSON → graceful handling
  8. Empty response — empty string → fallback to rules-based

Run: python tests/tst_phase9_llm_perturbation.py
Produces: findings/phase_9_llm_perturbation.yaml
"""

from __future__ import absolute_import, division, print_function

import json
import os
import sys
import tempfile
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

if 'libtbx' not in sys.modules:
    import types
    _base_dir = os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))
    libtbx = types.ModuleType('libtbx')
    libtbx.__path__ = [_base_dir]
    libtbx.langchain = types.ModuleType('libtbx.langchain')
    libtbx.langchain.__path__ = [_base_dir]
    libtbx.langchain.agent = types.ModuleType(
        'libtbx.langchain.agent')
    libtbx.langchain.agent.__path__ = [
        os.path.join(_base_dir, 'agent')]
    libtbx.langchain.knowledge = types.ModuleType(
        'libtbx.langchain.knowledge')
    libtbx.langchain.knowledge.__path__ = [
        os.path.join(_base_dir, 'knowledge')]
    sys.modules['libtbx'] = libtbx
    sys.modules['libtbx.langchain'] = libtbx.langchain
    sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
    sys.modules['libtbx.langchain.knowledge'] = \
        libtbx.langchain.knowledge
import knowledge.yaml_loader
if 'libtbx.langchain.knowledge.yaml_loader' not in sys.modules:
    sys.modules['libtbx.langchain.knowledge.yaml_loader'] = \
        knowledge.yaml_loader

from agent.graph_nodes import parse_intent_json
from agent.command_builder import CommandBuilder, CommandContext
from agent.phil_validator import validate_phil_strategy
from agent.workflow_state import (
    _categorize_files, _PHASE_COLUMN_CACHE)
import agent.workflow_state as ws

BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
FINDINGS_DIR = os.path.join(BASE_DIR, 'findings')


# =====================================================================
# TEST FIXTURES
# =====================================================================

def make_xray_fixture():
    """Create a standard X-ray fixture: PDB + MTZ.

    Returns (tmpdir, file_paths, categorized_files).
    Caller is responsible for shutil.rmtree(tmpdir).
    """
    tmpdir = tempfile.mkdtemp()
    try:
        pdb = os.path.join(tmpdir, 'model.pdb')
        mtz = os.path.join(tmpdir, 'data.mtz')
        with open(pdb, 'w') as f:
            for i in range(800):
                f.write(
                    "ATOM  %5d  CA  ALA A%4d"
                    "       %7.3f  %7.3f  %7.3f"
                    "  1.00 20.00\n"
                    % (i+1, i+1, i*1.0, 0.0, 0.0))
        with open(mtz, 'w') as f:
            f.write("MTZ_DUMMY\n" * 100)

        orig = ws._mtz_has_phase_columns
        ws._mtz_has_phase_columns = lambda f: False
        try:
            _PHASE_COLUMN_CACHE.clear()
            categorized = _categorize_files([pdb, mtz])
        finally:
            ws._mtz_has_phase_columns = orig

        return tmpdir, [pdb, mtz], categorized
    except Exception:
        shutil.rmtree(tmpdir)
        raise


# =====================================================================
# ITEM 4: FILENAME PERTURBATION
# =====================================================================

def test_nonexistent_llm_file_fallback():
    """LLM suggests a non-existent file. Command builder should
    fall back to category-based selection."""
    print("Test: nonexistent_llm_file_fallback")
    tmpdir, files, categorized = make_xray_fixture()
    try:
        ctx = CommandContext(
            cycle_number=2,
            experiment_type='xray',
            categorized_files=categorized,
            best_files={'model': files[0], 'data_mtz': files[1]},
            # LLM suggests a file that doesn't exist
            llm_files={'model': '/nonexistent/hallucinated.pdb',
                       'data_mtz': files[1]},
            files_local=True,
        )
        builder = CommandBuilder()
        cmd = builder.build('phenix.refine', files, ctx)
        assert cmd is not None, (
            "Command should build despite hallucinated filename")
        assert 'model.pdb' in cmd, (
            "Should fall back to real model.pdb, got: %s"
            % cmd[:80])
        print("  PASSED (fell back to real file)")
    finally:
        shutil.rmtree(tmpdir)


def test_wrong_extension_llm_file():
    """LLM suggests a .mtz where a .pdb is needed."""
    print("Test: wrong_extension_llm_file")
    tmpdir, files, categorized = make_xray_fixture()
    try:
        ctx = CommandContext(
            cycle_number=2,
            experiment_type='xray',
            categorized_files=categorized,
            best_files={'model': files[0], 'data_mtz': files[1]},
            llm_files={'model': files[1]},  # MTZ for model slot
            files_local=True,
        )
        builder = CommandBuilder()
        cmd = builder.build('phenix.refine', files, ctx)
        assert cmd is not None, (
            "Command should build despite wrong extension")
        # Builder should reject MTZ for model slot and fall back
        assert 'model.pdb' in cmd, (
            "Should use real PDB for model, got: %s" % cmd[:80])
        print("  PASSED (rejected wrong extension)")
    finally:
        shutil.rmtree(tmpdir)


def test_empty_llm_files():
    """LLM returns empty files dict."""
    print("Test: empty_llm_files")
    tmpdir, files, categorized = make_xray_fixture()
    try:
        ctx = CommandContext(
            cycle_number=2,
            experiment_type='xray',
            categorized_files=categorized,
            best_files={'model': files[0], 'data_mtz': files[1]},
            llm_files={},  # Empty
            files_local=True,
        )
        builder = CommandBuilder()
        cmd = builder.build('phenix.refine', files, ctx)
        assert cmd is not None, (
            "Command should build with empty llm_files")
        print("  PASSED")
    finally:
        shutil.rmtree(tmpdir)


# =====================================================================
# ITEM 5: PROGRAM PERTURBATION
# =====================================================================

def test_wrong_program_still_builds():
    """Build command for a different valid program than what routing
    suggested. Should still produce a valid command."""
    print("Test: wrong_program_still_builds")
    tmpdir, files, categorized = make_xray_fixture()
    try:
        ctx = CommandContext(
            cycle_number=2,
            experiment_type='xray',
            categorized_files=categorized,
            best_files={'model': files[0], 'data_mtz': files[1]},
            files_local=True,
        )
        builder = CommandBuilder()
        # Routing would suggest refine, but LLM chose xtriage
        cmd = builder.build('phenix.xtriage', files, ctx)
        assert cmd is not None, "xtriage should be buildable"
        assert cmd.startswith('phenix.xtriage'), (
            "Command should start with xtriage")
        print("  PASSED (wrong program still buildable)")
    finally:
        shutil.rmtree(tmpdir)


def test_unknown_program_returns_none():
    """Build command for a non-existent program. Should return None
    gracefully, not crash."""
    print("Test: unknown_program_returns_none")
    tmpdir, files, categorized = make_xray_fixture()
    try:
        ctx = CommandContext(
            cycle_number=2,
            experiment_type='xray',
            categorized_files=categorized,
            best_files={},
            files_local=True,
        )
        builder = CommandBuilder()
        cmd = builder.build('phenix.hallucinated_program',
                            files, ctx)
        # Should return None (unknown program has no YAML entry)
        assert cmd is None, (
            "Unknown program should return None, got: %s"
            % (cmd[:60] if cmd else cmd))
        print("  PASSED (returned None)")
    finally:
        shutil.rmtree(tmpdir)


# =====================================================================
# ITEM 6: PARAMETER INJECTION
# =====================================================================

def test_invalid_phil_stripped():
    """Invalid PHIL parameter should be stripped by validator."""
    print("Test: invalid_phil_stripped")
    strategy = {
        'nproc': '4',
        'xray_data.r_free_flags.generate': 'True',
        'hallucinated_param': 'crash_value',
        'another.fake.param': '42',
    }
    cleaned, stripped = validate_phil_strategy(
        'phenix.refine', strategy)

    stripped_keys = [k for k, v in stripped]
    assert 'hallucinated_param' in stripped_keys, (
        "hallucinated_param should be stripped")
    assert 'another.fake.param' in stripped_keys, (
        "another.fake.param should be stripped")
    assert 'nproc' in cleaned, (
        "nproc should survive (valid strategy_flag)")
    print("  PASSED (stripped: %s)" % stripped_keys)


def test_blocked_param_stripped():
    """Blocked params should be stripped even if they look valid."""
    print("Test: blocked_param_stripped")
    strategy = {
        'resolution': '3.5',
        'mask_atoms': 'True',  # Blocked for resolve_cryo_em
    }
    cleaned, stripped = validate_phil_strategy(
        'phenix.resolve_cryo_em', strategy)

    stripped_keys = [k for k, v in stripped]
    if 'mask_atoms' in stripped_keys:
        print("  PASSED (mask_atoms blocked)")
    elif 'mask_atoms' not in cleaned:
        # Not in cleaned and not in stripped — was never a
        # recognized param. This is acceptable.
        print("  PASSED (mask_atoms not recognized)")
    else:
        raise AssertionError(
            "mask_atoms should be stripped or not in cleaned, "
            "but found in cleaned: %s" % sorted(cleaned.keys()))


def test_empty_strategy_passthrough():
    """Empty strategy should pass through unchanged."""
    print("Test: empty_strategy_passthrough")
    cleaned, stripped = validate_phil_strategy(
        'phenix.refine', {})
    assert cleaned == {}, "Empty strategy should stay empty"
    assert stripped == [], "Nothing to strip"
    print("  PASSED")


def test_none_strategy_passthrough():
    """None strategy should pass through unchanged."""
    print("Test: none_strategy_passthrough")
    cleaned, stripped = validate_phil_strategy(
        'phenix.refine', None)
    assert cleaned is None, "None strategy should stay None"
    assert stripped == [], "Nothing to strip"
    print("  PASSED")


# =====================================================================
# ITEM 7: RESPONSE TRUNCATION
# =====================================================================

def test_truncated_json_handled():
    """Truncated JSON mid-object (unmatched braces) should raise
    ValueError, not crash with an unhandled exception."""
    print("Test: truncated_json_handled")
    truncated = '{"program": "phenix.refine", "files": {"model": '
    try:
        parse_intent_json(truncated)
        raise AssertionError(
            "Should raise ValueError for unmatched braces")
    except (ValueError, json.JSONDecodeError):
        print("  PASSED (raised ValueError)")
    except Exception as e:
        raise AssertionError(
            "Should raise ValueError, got %s: %s"
            % (type(e).__name__, e))


def test_truncated_after_key():
    """JSON truncated after a key (no value, unmatched braces)."""
    print("Test: truncated_after_key")
    truncated = '{"program": "phenix.refine", "strategy":'
    try:
        parse_intent_json(truncated)
        raise AssertionError(
            "Should raise ValueError for unmatched braces")
    except (ValueError, json.JSONDecodeError):
        print("  PASSED (raised ValueError)")
    except Exception as e:
        raise AssertionError(
            "Should raise ValueError, got %s: %s"
            % (type(e).__name__, e))


def test_markdown_wrapped_json():
    """LLM wraps JSON in markdown code block."""
    print("Test: markdown_wrapped_json")
    wrapped = '```json\n{"program": "phenix.refine", ' \
              '"files": {}, "strategy": {}, ' \
              '"reasoning": "test"}\n```'
    intent = parse_intent_json(wrapped)
    assert intent['program'] == 'phenix.refine', (
        "Should unwrap markdown and parse")
    print("  PASSED")


def test_json_with_preamble():
    """LLM adds conversational text before JSON."""
    print("Test: json_with_preamble")
    preamble = 'Here is my analysis:\n' \
               '{"program": "phenix.refine", ' \
               '"files": {}, "strategy": {}, ' \
               '"reasoning": "test"}'
    intent = parse_intent_json(preamble)
    assert intent['program'] == 'phenix.refine', (
        "Should find JSON after preamble")
    print("  PASSED")


# =====================================================================
# ITEM 8: EMPTY RESPONSE
# =====================================================================

def test_empty_string_handled():
    """Empty string should raise ValueError."""
    print("Test: empty_string_handled")
    try:
        parse_intent_json('')
        raise AssertionError("Should have raised ValueError")
    except (ValueError, json.JSONDecodeError):
        print("  PASSED (raised ValueError)")


def test_whitespace_only_handled():
    """Whitespace-only string should raise ValueError."""
    print("Test: whitespace_only_handled")
    try:
        parse_intent_json('   \n\t  ')
        raise AssertionError("Should have raised ValueError")
    except (ValueError, json.JSONDecodeError):
        print("  PASSED (raised ValueError)")


def test_non_json_text_handled():
    """Plain text (no JSON) should raise ValueError."""
    print("Test: non_json_text_handled")
    try:
        parse_intent_json(
            'I think we should run phenix.refine next.')
        raise AssertionError("Should have raised ValueError")
    except (ValueError, json.JSONDecodeError):
        print("  PASSED (raised ValueError)")


def test_intent_type_sanitization():
    """parse_intent_json sanitizes field types even if LLM
    returns wrong types."""
    print("Test: intent_type_sanitization")
    # files as list instead of dict, program as int
    bad_types = '{"program": 42, "files": ["a.pdb"], ' \
                '"strategy": "nproc=4", "reasoning": null}'
    intent = parse_intent_json(bad_types)
    assert isinstance(intent['files'], dict), (
        "files should be coerced to dict")
    assert isinstance(intent['strategy'], dict), (
        "strategy should be coerced to dict")
    assert isinstance(intent['program'], str), (
        "program should be coerced to str")
    assert isinstance(intent['reasoning'], str), (
        "reasoning should be coerced to str")
    print("  PASSED")


# =====================================================================
# RUN ALL
# =====================================================================

def run_all_tests():
    """Run all Phase 9 perturbation tests."""
    print("Phase 9: LLM Perturbation Tests")
    print()

    tests = [
        # Item 4: Filename perturbation
        test_nonexistent_llm_file_fallback,
        test_wrong_extension_llm_file,
        test_empty_llm_files,
        # Item 5: Program perturbation
        test_wrong_program_still_builds,
        test_unknown_program_returns_none,
        # Item 6: Parameter injection
        test_invalid_phil_stripped,
        test_blocked_param_stripped,
        test_empty_strategy_passthrough,
        test_none_strategy_passthrough,
        # Item 7: Response truncation
        test_truncated_json_handled,
        test_truncated_after_key,
        # Format extraction (LLM wrapping)
        test_markdown_wrapped_json,
        test_json_with_preamble,
        # Item 8: Empty/invalid response
        test_empty_string_handled,
        test_whitespace_only_handled,
        test_non_json_text_handled,
        # Type coercion
        test_intent_type_sanitization,
    ]

    passed = 0
    failed = 0
    failures = []

    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            failures.append({
                'test': t.__name__,
                'error': str(e)[:200],
            })
            print("  FAILED: %s" % str(e)[:100])

    print()
    print("  Results: %d passed, %d failed out of %d"
          % (passed, failed, len(tests)))

    status = 'PASS' if failed == 0 else 'FAIL'

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_9_llm_perturbation.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 9: LLM Perturbation Findings\n")
        f.write("phase: 9\n")
        f.write("name: llm_perturbation\n")
        f.write("status: %s\n" % status)
        f.write("total_tests: %d\n" % len(tests))
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        f.write("\ntest_categories:\n")
        f.write("  filename_perturbation: 3 tests\n")
        f.write("  program_perturbation: 2 tests\n")
        f.write("  parameter_injection: 4 tests\n")
        f.write("  response_truncation: 2 tests\n")
        f.write("  format_extraction: 2 tests\n")
        f.write("  empty_invalid_response: 3 tests\n")
        f.write("  type_coercion: 1 test\n")
        f.write("\nnote: \"Playback tests (items 1-3 from plan)\n")
        f.write("  require Run 18 session data. Only perturbation\n")
        f.write("  tests (items 4-8) are implemented here.\"\n")
        if failures:
            f.write("\nfailures:\n")
            for fail in failures:
                f.write("  - test: %s\n" % fail['test'])
                f.write("    error: \"%s\"\n"
                        % fail['error'][:100])

    print("  Findings: %s" % path)
    print("  Phase 9 overall: %s" % status)

    if failed > 0:
        raise AssertionError(
            "Phase 9 LLM perturbation FAILED: %d/%d tests"
            % (failed, len(tests)))

    return {'status': status, 'passed': passed,
            'failed': failed}


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] == 'PASS' else 1)
