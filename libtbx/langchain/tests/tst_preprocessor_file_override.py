"""
Sandbox tests for v118 Section A2: deterministic file-list override
in the preprocessor result.

The override is a post-LLM step in advice_preprocessor.py
(_ensure_file_list_in_processed_advice) that enforces the
architectural invariant: file_list_hint (the disk reality supplied
to the agent) is authoritative over the LLM's section-1 output.

These tests verify the override behavior in isolation, without
invoking any LLM.  They match the K_O1-K_O4 specifications in the
v118 plan §4 Section A.

K_O1: section 1 = "None"; file_list_hint has files
      → override fires, section 1 rewritten

K_O2: section 1 has some files; file_list_hint has more
      → override fires, union written

K_O3: section 1 already complete (matches file_list_hint)
      → no change

K_O4: section 1 = "None"; file_list_hint also empty
      → no change

Additional K_O5-K_O7 boundary tests verify:
- Path normalization (basenames extracted from absolute paths)
- Markdown-formatted section 1
- Section 1 missing entirely (graceful degradation)
"""

import os
import sys

# Path setup for sandbox execution
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Stub libtbx packages for sandbox import
import types
for mod in ("libtbx", "libtbx.langchain", "libtbx.langchain.agent"):
    if mod not in sys.modules:
        sys.modules[mod] = types.ModuleType(mod)

from agent.advice_preprocessor import _ensure_file_list_in_processed_advice


def _make_advice(section1_body, with_other_sections=True):
    """Helper: build a minimal processed-advice string."""
    if with_other_sections:
        return (
            "1. **Input Files Found**: %s\n\n"
            "2. **Experiment Type**: X-ray\n\n"
            "3. **Primary Goal**: Refine the model.\n\n"
            "4. **Key Parameters**:\n  - None\n\n"
            "5. **Program Parameters**: None\n\n"
            "6. **Special Instructions**: None\n\n"
            "7. **Stop Condition**: None\n"
            % section1_body
        )
    return "1. **Input Files Found**: %s\n" % section1_body


def _section1_body(advice):
    """Helper: extract the section-1 body from processed advice."""
    import re
    m = re.search(
        r'1\.\s*\*?\*?Input\s+Files\s+Found\*?\*?\s*:?(.*?)(?=\n\s*2\.\s|\Z)',
        advice,
        re.IGNORECASE | re.DOTALL,
    )
    return m.group(1).strip() if m else None


def test_k_o1_none_with_hint():
    """K_O1: section 1 'None', file_list_hint has files → override fires."""
    advice = _make_advice("None")
    hint = ["nsf-d2_noligand.pdb", "nsf-d2.mtz"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    assert body == "nsf-d2_noligand.pdb, nsf-d2.mtz", (
        "K_O1 failed: expected 'nsf-d2_noligand.pdb, nsf-d2.mtz', got %r"
        % body)
    print("  PASS: K_O1 (None + hint → override)")


def test_k_o2_partial_with_hint():
    """K_O2: section 1 has some files, hint has more → union written."""
    advice = _make_advice("model.pdb")
    hint = ["model.pdb", "data.mtz", "seq.fa"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    # Order: existing first, then additions
    assert body == "model.pdb, data.mtz, seq.fa", (
        "K_O2 failed: expected 'model.pdb, data.mtz, seq.fa', got %r"
        % body)
    print("  PASS: K_O2 (partial + hint → union)")


def test_k_o3_complete_no_change():
    """K_O3: section 1 already complete → no change."""
    advice = _make_advice("model.pdb, data.mtz")
    hint = ["model.pdb", "data.mtz"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    assert out == advice, "K_O3 failed: advice changed when it shouldn't"
    print("  PASS: K_O3 (complete → no change)")


def test_k_o4_none_no_hint():
    """K_O4: section 1 'None', hint empty → no change."""
    advice = _make_advice("None")
    out_empty_list = _ensure_file_list_in_processed_advice(advice, [])
    out_none = _ensure_file_list_in_processed_advice(advice, None)
    assert out_empty_list == advice, "K_O4 failed for empty list"
    assert out_none == advice, "K_O4 failed for None"
    print("  PASS: K_O4 (None + empty hint → no change)")


def test_k_o5_absolute_paths_normalized():
    """K_O5: hint with absolute paths → basenames extracted."""
    advice = _make_advice("None")
    hint = [
        "/Users/terwill/Documents/nsf-d2-ligand_terwill_1/nsf-d2_noligand.pdb",
        "/Users/terwill/Documents/nsf-d2-ligand_terwill_1/nsf-d2.mtz",
    ]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    assert body == "nsf-d2_noligand.pdb, nsf-d2.mtz", (
        "K_O5 failed: expected basenames, got %r" % body)
    print("  PASS: K_O5 (absolute paths normalized to basenames)")


def test_k_o6_non_markdown_section():
    """K_O6: section 1 without markdown bold → still recognized."""
    advice = (
        "1. Input Files Found: None\n\n"
        "2. Experiment Type: X-ray\n"
    )
    hint = ["a.pdb", "b.mtz"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    assert body == "a.pdb, b.mtz", (
        "K_O6 failed: expected 'a.pdb, b.mtz', got %r" % body)
    print("  PASS: K_O6 (non-markdown section recognized)")


def test_k_o7_section_missing_entirely():
    """K_O7: no section 1 at all → graceful degradation, no change."""
    advice = "2. Experiment Type: X-ray\n\n3. Goal: something"
    hint = ["a.pdb"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    # Returns unchanged when section 1 can't be found
    assert out == advice, (
        "K_O7 failed: rewrote advice when section 1 was missing")
    print("  PASS: K_O7 (missing section 1 → graceful no-op)")


def test_k_o8_empty_advice():
    """K_O8: empty/None processed_advice → graceful no-op."""
    hint = ["a.pdb"]
    assert _ensure_file_list_in_processed_advice("", hint) == ""
    assert _ensure_file_list_in_processed_advice(None, hint) is None
    print("  PASS: K_O8 (empty advice → graceful no-op)")


def test_k_o9_duplicate_hint_files():
    """K_O9: hint with duplicates → deduplicated."""
    advice = _make_advice("None")
    hint = ["a.pdb", "b.mtz", "a.pdb"]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    assert body == "a.pdb, b.mtz", (
        "K_O9 failed: expected dedup, got %r" % body)
    print("  PASS: K_O9 (duplicate hints deduplicated)")


def test_k_o10_logging():
    """K_O10: log callable receives override message when override fires."""
    advice = _make_advice("None")
    hint = ["model.pdb", "data.mtz"]
    logs = []
    _ensure_file_list_in_processed_advice(advice, hint, log=logs.append)
    assert any("PREPROCESSOR_OVERRIDE" in m for m in logs), (
        "K_O10 failed: no PREPROCESSOR_OVERRIDE log message; logs=%r"
        % logs)
    print("  PASS: K_O10 (override emits log message)")


def test_k_o11_no_log_when_no_override():
    """K_O11: log not called when no override fires."""
    advice = _make_advice("model.pdb, data.mtz")
    hint = ["model.pdb", "data.mtz"]
    logs = []
    _ensure_file_list_in_processed_advice(advice, hint, log=logs.append)
    # No additions needed → no log
    assert not logs, "K_O11 failed: unexpected log when no override; logs=%r" % logs
    print("  PASS: K_O11 (no log when no override needed)")


def test_k_o12_aiagent_210_reproducer():
    """K_O12: exact production scenario from AIAgent_210."""
    # The actual preprocessor output from the run
    advice = (
        "1. **Input Files Found**: None\n\n"
        "2. **Experiment Type**: X-ray refinement\n\n"
        "3. **Primary Goal**: Refine the provided model "
        "(`nsf-d2_noligand.pdb`) against the reflection data "
        "(`nsf-d2.mtz`).\n\n"
        "4. **Key Parameters**:\n"
        "   - Wavelength: None\n"
        "   - Resolution limit: None\n"
        "   - Number of expected sites: None\n"
        "   - Heavy atom type: None\n"
        "   - Space group: None\n\n"
        "5. **Program Parameters**:\n"
        "None\n\n"
        "6. **Special Instructions**:\n"
        "None\n\n"
        "7. **Stop Condition**: None\n"
    )
    # The actual original_files from the production run
    hint = [
        "/Users/terwill/Documents/nsf-d2-ligand_terwill_1/nsf-d2_noligand.pdb",
        "/Users/terwill/Documents/nsf-d2-ligand_terwill_1/nsf-d2.mtz",
    ]
    out = _ensure_file_list_in_processed_advice(advice, hint)
    body = _section1_body(out)
    assert body == "nsf-d2_noligand.pdb, nsf-d2.mtz", (
        "K_O12 (AIAgent_210 reproducer) failed: got %r" % body)
    # Sanity check: other sections unchanged
    assert "Experiment Type" in out
    assert "Stop Condition" in out
    print("  PASS: K_O12 (AIAgent_210 production reproducer)")


def run_all_tests():
    tests = [
        ("K_O1_none_with_hint", test_k_o1_none_with_hint),
        ("K_O2_partial_with_hint", test_k_o2_partial_with_hint),
        ("K_O3_complete_no_change", test_k_o3_complete_no_change),
        ("K_O4_none_no_hint", test_k_o4_none_no_hint),
        ("K_O5_absolute_paths", test_k_o5_absolute_paths_normalized),
        ("K_O6_non_markdown_section", test_k_o6_non_markdown_section),
        ("K_O7_section_missing", test_k_o7_section_missing_entirely),
        ("K_O8_empty_advice", test_k_o8_empty_advice),
        ("K_O9_duplicate_hints", test_k_o9_duplicate_hint_files),
        ("K_O10_logging", test_k_o10_logging),
        ("K_O11_no_log_no_override", test_k_o11_no_log_when_no_override),
        ("K_O12_aiagent_210_reproducer", test_k_o12_aiagent_210_reproducer),
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
