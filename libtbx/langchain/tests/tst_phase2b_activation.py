"""K_H7: Phase 2B activation tests (v119.H7).

Sandbox-only suite for Phase 2B activation in
phenix_ai/run_ai_analysis.py.  Tests the building blocks the
production code composes, plus the contract invariants the
plan guarantees.

  §A  scan_files_in_advice contract                  4 tests
      (decision-matrix shape: scanner picks files;
      type/sorting/dedup invariants)
  §B  fallback semantics                             3 tests
      (ImportError fallback to LLM extraction; both
      fail = graceful empty result)
  §C  recall metrics                                 3 tests
      (zero-division guards; happy-path math; perfect
      overlap)
  §D  integration smoke                              5 tests
      (mock realism helper; production marker format
      regression; partial-overlap math; production
      decision-structure regression; telemetry
      independence from extracted_files post-H7)

Total: 15 tests.  No libtbx dependency — all tests use either
direct calls to raw_advice_scanner (sandbox-importable as
agent.raw_advice_scanner) or pure-Python recall computation.

All tests PASS in sandbox with no SKIPs.  Under PHENIX they
also PASS (the scanner is importable via libtbx too).
"""
from __future__ import absolute_import, division, print_function

import os
import sys


# =====================================================================
# Import helpers — sandbox + PHENIX both
# =====================================================================

def _try_import_scanner():
    """Import the scanner module, sandbox-tolerant.

    The production code in run_ai_analysis.py uses lazy import:
      try: from libtbx.langchain.agent.raw_advice_scanner import ...
      except ImportError: from agent.raw_advice_scanner import ...

    Tests mirror that pattern for the same robustness.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    try:
        try:
            from libtbx.langchain.agent.raw_advice_scanner import (
                scan_files_in_advice)
            return scan_files_in_advice, None
        except ImportError:
            from agent.raw_advice_scanner import scan_files_in_advice
            return scan_files_in_advice, None
    except ImportError as e:
        return None, str(e)


# =====================================================================
# Mock realism helper (Gemini rev 2.2 recommendation)
# =====================================================================

def _make_mock_scanner_return(files):
    """Canonicalize a list of filenames to match scanner contract.

    Sorts ASCII-ascending case-insensitively and dedups
    case-insensitively (last-wins on casing collision).  Use
    this when constructing mock scanner returns to ensure §A
    invariant assertions don't trigger on fixture bugs rather
    than code regressions.

    Example:
      _make_mock_scanner_return(["z.pdb", "a.mtz", "Z.PDB"])
      # returns: ["a.mtz", "Z.PDB"]  (or ["a.mtz", "z.pdb"];
      # last-wins on case collision is deterministic per scanner
      # docstring at raw_advice_scanner.py:160)
    """
    seen = {}
    for f in files:
        seen[f.lower()] = f
    return sorted(seen.values(), key=str.lower)


# =====================================================================
# Recall computation helper (mirrors production code in
# run_ai_analysis.py post-H7)
# =====================================================================

def _compute_recall_pair(llm_files, regex_files):
    """Compute (scanner_recall_against_llm, llm_recall_against_scanner).

    Mirrors the exact computation in run_ai_analysis.py post-H7.
    Both metrics zero-division-guarded: empty denominator → 1.0
    (mathematically sound for recall).
    """
    # Case-insensitive comparison (matches production normalization)
    llm_keys = set(f.lower() for f in llm_files)
    regex_keys = set(f.lower() for f in regex_files)
    intersection = llm_keys & regex_keys

    scanner_recall = (
        len(intersection) / len(llm_keys) if llm_keys else 1.0)
    llm_recall = (
        len(intersection) / len(regex_keys) if regex_keys else 1.0)

    return scanner_recall, llm_recall


# =====================================================================
# §A: scan_files_in_advice contract (4 tests)
# =====================================================================

def test_scanner_finds_mentioned_files():
    """Scanner returns files mentioned in advice text."""
    scan, err = _try_import_scanner()
    if scan is None:
        print("  SKIP: cannot import raw_advice_scanner (%s)" % err)
        return

    advice = "Please analyze model.pdb and the data file data.mtz."
    result = scan(advice, None)

    # Should contain both files (case may vary by scanner regex)
    result_lower = [f.lower() for f in result]
    assert "model.pdb" in result_lower, (
        "Scanner should find model.pdb in advice; got: %s" % result)
    assert "data.mtz" in result_lower, (
        "Scanner should find data.mtz in advice; got: %s" % result)
    print("  PASS: test_scanner_finds_mentioned_files")


def test_scanner_returns_list_of_strings():
    """Type invariant: scanner returns list of strings."""
    scan, err = _try_import_scanner()
    if scan is None:
        print("  SKIP: cannot import raw_advice_scanner (%s)" % err)
        return

    result = scan("Run model.pdb and refinement.eff.", None)
    assert isinstance(result, list), (
        "Scanner should return a list, got %s" % type(result).__name__)
    assert all(isinstance(f, str) for f in result), (
        "Scanner list should contain only strings; got: %s" % result)
    print("  PASS: test_scanner_returns_list_of_strings")


def test_scanner_sorted_ascii_ascending_case_insensitive():
    """Sorting invariant: ASCII-ascending case-insensitively."""
    scan, err = _try_import_scanner()
    if scan is None:
        print("  SKIP: cannot import raw_advice_scanner (%s)" % err)
        return

    # Mention several files in a deliberately non-sorted order
    advice = "Files: z_file.pdb, a_file.mtz, M_file.eff, b_file.cif"
    result = scan(advice, None)

    lowered = [f.lower() for f in result]
    assert lowered == sorted(lowered), (
        "Scanner output not sorted case-insensitively; got: %s" % result)
    print("  PASS: test_scanner_sorted_ascii_ascending_case_insensitive")


def test_scanner_deduped_case_insensitive():
    """Dedup invariant: no case-insensitive duplicates."""
    scan, err = _try_import_scanner()
    if scan is None:
        print("  SKIP: cannot import raw_advice_scanner (%s)" % err)
        return

    # Mention the same file with different casing
    advice = "Use Model.PDB, then model.pdb, then MODEL.PDB."
    result = scan(advice, None)

    lowered = [f.lower() for f in result]
    assert len(lowered) == len(set(lowered)), (
        "Scanner output has case-insensitive duplicates; got: %s" % result)
    print("  PASS: test_scanner_deduped_case_insensitive")


# =====================================================================
# §B: fallback semantics (3 tests)
# =====================================================================
#
# These tests exercise the decision-matrix from
# run_ai_analysis.py:1009-1067 (the post-H7 extracted_files block)
# by replaying the fallback chain.  We don't call the full
# run_advice_preprocessing function — it's a giant function with
# many arguments.  Instead we test the composed building blocks.


def _simulate_h7_extraction(scanner_fn, llm_fallback_fn,
                            raw_advice, processed_advice):
    """Replay the H7 scanner-first-with-fallback decision logic.

    Mirrors the structure of phenix_ai/run_ai_analysis.py:1009-1067.
    Returns (extracted_files, debug_log) so tests can assert on
    both the outcome and the diagnostic log.
    """
    debug_log = []
    extracted_files = []
    try:
        extracted_files = scanner_fn(raw_advice, None)
        if extracted_files:
            debug_log.append(
                "Extracted files from advice (scanner): %s"
                % extracted_files)
    except Exception as e:
        debug_log.append(
            "raw_advice_scanner failed: %s; "
            "falling back to LLM extraction" % e)
        if processed_advice and processed_advice != raw_advice:
            try:
                extracted_files = llm_fallback_fn(processed_advice)
                if extracted_files:
                    debug_log.append(
                        "Extracted files from advice (LLM fallback): %s"
                        % extracted_files)
            except Exception as e2:
                debug_log.append(
                    "LLM fallback extraction also failed: %s; "
                    "extracted_files remains empty" % e2)
    return extracted_files, debug_log


def test_fallback_scanner_raises_llm_succeeds():
    """Scanner raises; LLM fallback path produces extracted_files."""
    def raising_scanner(advice, hint):
        raise RuntimeError("scanner broken")

    def working_llm(processed):
        return _make_mock_scanner_return(["model.pdb", "data.mtz"])

    extracted, debug_log = _simulate_h7_extraction(
        scanner_fn=raising_scanner,
        llm_fallback_fn=working_llm,
        raw_advice="run model.pdb",
        processed_advice="Input Files Found: model.pdb")

    assert extracted == ["data.mtz", "model.pdb"], (
        "LLM fallback should have populated extracted_files; got: %s"
        % extracted)
    assert any("scanner broken" in msg for msg in debug_log), (
        "debug_log should record scanner failure; got: %s" % debug_log)
    assert any("LLM fallback" in msg for msg in debug_log), (
        "debug_log should record LLM fallback succeeded; got: %s"
        % debug_log)
    print("  PASS: test_fallback_scanner_raises_llm_succeeds")


def test_fallback_both_fail_graceful():
    """Scanner AND LLM both fail; function returns gracefully."""
    def raising_scanner(advice, hint):
        raise RuntimeError("scanner broken")

    def raising_llm(processed):
        raise RuntimeError("llm extraction also broken")

    extracted, debug_log = _simulate_h7_extraction(
        scanner_fn=raising_scanner,
        llm_fallback_fn=raising_llm,
        raw_advice="run model.pdb",
        processed_advice="Input Files Found: model.pdb")

    # Graceful: empty result, NOT a crash
    assert extracted == [], (
        "Both-paths-failed should produce empty extracted_files; "
        "got: %s" % extracted)
    # BOTH errors should be in debug_log (per Gemini rev 2)
    assert any("scanner broken" in msg for msg in debug_log), (
        "debug_log should record scanner error; got: %s" % debug_log)
    assert any("llm extraction also broken" in msg
                for msg in debug_log), (
        "debug_log should record LLM fallback error; got: %s"
        % debug_log)
    print("  PASS: test_fallback_both_fail_graceful")


def test_fallback_scanner_succeeds_no_llm_needed():
    """Scanner returns files; LLM fallback NEVER called."""
    llm_called = [False]

    def working_scanner(advice, hint):
        return _make_mock_scanner_return(["model.pdb"])

    def llm_fn(processed):
        llm_called[0] = True
        return ["should_never_appear.xyz"]

    extracted, debug_log = _simulate_h7_extraction(
        scanner_fn=working_scanner,
        llm_fallback_fn=llm_fn,
        raw_advice="run model.pdb",
        processed_advice="Input Files Found: model.pdb")

    assert extracted == ["model.pdb"]
    assert not llm_called[0], (
        "LLM fallback should NOT be called when scanner succeeds")
    assert "should_never_appear.xyz" not in str(extracted)
    print("  PASS: test_fallback_scanner_succeeds_no_llm_needed")


# =====================================================================
# §C: recall metrics (3 tests)
# =====================================================================


def test_recall_perfect_overlap():
    """Both lists identical → both recalls = 1.0."""
    llm = ["model.pdb", "data.mtz"]
    regex = ["model.pdb", "data.mtz"]
    scanner_recall, llm_recall = _compute_recall_pair(llm, regex)
    assert scanner_recall == 1.0
    assert llm_recall == 1.0
    print("  PASS: test_recall_perfect_overlap")


def test_recall_zero_division_empty_llm():
    """Empty llm_files → scanner_recall_against_llm = 1.0 (zero-div guard)."""
    llm = []
    regex = ["data.mtz", "model.pdb"]
    scanner_recall, llm_recall = _compute_recall_pair(llm, regex)
    # Guard kicks in: 100% of zero targets captured = 1.0
    assert scanner_recall == 1.0, (
        "Empty llm_files should produce scanner_recall=1.0 via "
        "zero-div guard; got %s" % scanner_recall)
    # llm_recall sees intersection=0, |regex|=2 → 0/2 = 0.0
    assert llm_recall == 0.0
    print("  PASS: test_recall_zero_division_empty_llm")


def test_recall_zero_division_empty_regex():
    """Empty regex_files → llm_recall_against_scanner = 1.0 (zero-div guard)."""
    llm = ["model.pdb", "data.mtz"]
    regex = []
    scanner_recall, llm_recall = _compute_recall_pair(llm, regex)
    # scanner_recall sees intersection=0, |llm|=2 → 0/2 = 0.0
    assert scanner_recall == 0.0
    # Guard kicks in for the other direction
    assert llm_recall == 1.0, (
        "Empty regex_files should produce llm_recall=1.0 via "
        "zero-div guard; got %s" % llm_recall)
    print("  PASS: test_recall_zero_division_empty_regex")


# =====================================================================
# §D: integration smoke (3 tests)
# =====================================================================


def test_mock_helper_satisfies_scanner_invariants():
    """The _make_mock_scanner_return helper produces values that
    satisfy the §A invariants — so tests using it as a fixture
    don't accidentally trigger §A defensive assertions."""
    # Pathological input: unsorted, case-mixed, duplicated
    pathological = ["z_file.txt", "a_file.txt", "Z_File.TXT", "B.pdb"]
    canonical = _make_mock_scanner_return(pathological)

    # Verify §A invariants on the helper's output
    assert isinstance(canonical, list)
    assert all(isinstance(f, str) for f in canonical)
    lowered = [f.lower() for f in canonical]
    assert lowered == sorted(lowered), (
        "_make_mock_scanner_return broke sorting invariant; "
        "got: %s" % canonical)
    assert len(lowered) == len(set(lowered)), (
        "_make_mock_scanner_return broke dedup invariant; "
        "got: %s" % canonical)
    print("  PASS: test_mock_helper_satisfies_scanner_invariants")


def test_marker_format_includes_new_recall_keys():
    """The post-H7 production [STEP_1F] marker format includes both
    new recall keys.

    Regression guard: if a future refactor accidentally removes
    one of the recall fields from the production marker string,
    this test catches it without needing a live LLM run.

    Strategy: grep the production source for the marker format
    string to verify both new fields are present.  This is more
    robust than building a parallel marker locally and asserting
    on the local copy (which would only catch typos in the test
    itself, not production-code regressions).
    """
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    src_path = os.path.join(
        repo_root, "phenix_ai", "run_ai_analysis.py")

    if not os.path.isfile(src_path):
        print("  SKIP: cannot locate phenix_ai/run_ai_analysis.py")
        return

    with open(src_path) as f:
        src = f.read()

    # Both new field names must appear in the production marker
    # string.  These are field-name literals embedded in a
    # format-string passed to the % operator.
    assert "scanner_recall_against_llm=" in src, (
        "Production marker missing scanner_recall_against_llm "
        "field literal in run_ai_analysis.py")
    assert "llm_recall_against_scanner=" in src, (
        "Production marker missing llm_recall_against_scanner "
        "field literal in run_ai_analysis.py")

    # The format string should use %.4f for both recall values
    # (per plan rev 2.2; consistent decimal-place formatting for
    # log-parser ergonomics).
    assert "scanner_recall_against_llm=%.4f" in src, (
        "scanner_recall_against_llm not formatted as %%.4f")
    assert "llm_recall_against_scanner=%.4f" in src, (
        "llm_recall_against_scanner not formatted as %%.4f")

    # Pre-H7 fields must still be present (no accidental removal).
    for required in ("preprocessor_mode=%s",
                      "scanner_version=%s",
                      "llm_files=%s", "regex_files=%s",
                      "in_llm_only=%s", "in_regex_only=%s"):
        assert required in src, (
            "Production marker missing pre-H7 field %r" % required)

    # The computation must use the case-insensitive key sets
    # (consistent with how pre-H7 set-difference is computed).
    assert "_intersection_keys = _llm_keys & _regex_keys" in src, (
        "Recall numerator must use _llm_keys & _regex_keys "
        "(case-insensitive intersection); pin avoids drift to a "
        "different intersection definition.")
    print("  PASS: test_marker_format_includes_new_recall_keys")


def test_partial_overlap_realistic_case():
    """Realistic case: scanner finds 4 files, LLM found 3 of them,
    LLM also found one the scanner missed.

    Mirrors the kind of asymmetry the H4.1 corpus measurements
    showed (scanner ≈ 0.9810 recall against LLM).
    """
    llm = ["model.pdb", "data.mtz", "refine.eff", "unique_to_llm.txt"]
    regex = ["model.pdb", "data.mtz", "refine.eff", "extra_from_scan.cif"]

    scanner_recall, llm_recall = _compute_recall_pair(llm, regex)
    # Intersection: 3 files (the three common ones)
    # |llm| = 4 → scanner_recall = 3/4 = 0.75
    # |regex| = 4 → llm_recall = 3/4 = 0.75
    assert abs(scanner_recall - 0.75) < 1e-9, (
        "scanner_recall expected 0.75; got %s" % scanner_recall)
    assert abs(llm_recall - 0.75) < 1e-9, (
        "llm_recall expected 0.75; got %s" % llm_recall)
    print("  PASS: test_partial_overlap_realistic_case")


def test_production_decision_structure_present():
    """Pin the production code's scanner-first-with-fallback
    decision structure as a source-level invariant.

    The K_H7 §B tests use _simulate_h7_extraction to verify the
    decision logic.  But _simulate_h7_extraction is a separate
    code path that mirrors production.  If production diverges
    (refactored to a different structure), §B's coverage no
    longer protects production.

    This test grep's the actual production source to verify the
    expected structure is present.  It catches refactors that
    accidentally drop the fallback path.
    """
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    src_path = os.path.join(
        repo_root, "phenix_ai", "run_ai_analysis.py")

    if not os.path.isfile(src_path):
        print("  SKIP: cannot locate phenix_ai/run_ai_analysis.py")
        return

    with open(src_path) as f:
        src = f.read()

    # Scanner-first call must be present (Q2-strict: None hint)
    assert "extracted_files = scan_files_in_advice(raw_advice, None)" in src, (
        "Production must call scan_files_in_advice(raw_advice, None) "
        "for the H7 primary extraction.  A refactor to file_list "
        "hint or different argument order would break the consumer "
        "contract documented in programs/ai_agent.py:8142.")

    # Fallback to LLM extraction must be present
    assert "extract_files_from_processed_advice" in src, (
        "Production must retain extract_files_from_processed_advice "
        "as the fallback path; dropping it would lose defense-in-depth "
        "against scanner failures.")

    # Scanner import must use libtbx → relative fallback (Gemini rev 2)
    assert "from libtbx.langchain.agent.raw_advice_scanner import" in src
    assert "from agent.raw_advice_scanner import scan_files_in_advice" in src, (
        "Production scanner import must mirror libtbx → relative "
        "fallback pattern.  Bare libtbx import would silently fail "
        "in non-PHENIX environments.")

    # LLM-extractor fallback must ALSO use libtbx → relative
    # (the Gemini rev 2 critical fix)
    assert "from libtbx.langchain.agent.advice_preprocessor import" in src
    assert "from agent.advice_preprocessor import" in src, (
        "Production LLM-fallback import must mirror libtbx → "
        "relative pattern (Gemini rev 2 critical fix); a bare "
        "libtbx import in the fallback would silently fail and "
        "be swallowed by the outer except handler.")

    print("  PASS: test_production_decision_structure_present")


def test_telemetry_independence_from_extracted_files():
    """The [STEP_1F] block must compute LLM extraction INDEPENDENTLY
    of extracted_files, post-H7.

    Pre-H7 the [STEP_1F] block consumed `extracted_files` as the
    LLM's view, because extracted_files WAS the LLM's output.
    Post-H7 extracted_files is scanner-derived, so the [STEP_1F]
    block must call extract_files_from_processed_advice itself
    to obtain the LLM view for comparison.

    Without this, the [STEP_1F] block would compare scanner-output
    to scanner-output and the recall numbers would always be
    inflated (close to 1.0) — defeating the telemetry's purpose
    of detecting scanner-vs-LLM drift.

    This test pins the property by checking that the production
    code in run_ai_analysis.py does NOT use `extracted_files` as
    the input to _step_1f_normalize for llm_normalized.
    """
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    src_path = os.path.join(
        repo_root, "phenix_ai", "run_ai_analysis.py")

    if not os.path.isfile(src_path):
        print("  SKIP: cannot locate phenix_ai/run_ai_analysis.py")
        return

    with open(src_path) as f:
        src = f.read()

    # The post-H7 code must call extract_files_from_processed_advice
    # INSIDE the [STEP_1F] block (in addition to the fallback path)
    # to compute the LLM extraction independently for telemetry.
    # Look for the dedicated _llm_extracted variable that the H7
    # fix introduces.
    assert "_llm_extracted = extract_files_from_processed_advice(" in src, (
        "Post-H7 [STEP_1F] block must extract LLM files "
        "independently for telemetry (not reuse extracted_files, "
        "which is now scanner-derived).  Look for _llm_extracted "
        "assignment in the [STEP_1F] try block.")

    # And the normalization must use _llm_extracted, not extracted_files
    assert "llm_normalized = _step_1f_normalize(_llm_extracted)" in src, (
        "llm_normalized must come from _llm_extracted, not "
        "extracted_files (post-H7 extracted_files is scanner-derived).")

    # Sanity: extracted_files is still used for its own purposes
    # (returned to caller, used in debug_log) — just not as the
    # LLM-comparison axis for telemetry.
    assert "extracted_files = scan_files_in_advice(" in src, (
        "Scanner-first switch should still be present")

    print("  PASS: test_telemetry_independence_from_extracted_files")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: scanner contract (4)
    test_scanner_finds_mentioned_files()
    test_scanner_returns_list_of_strings()
    test_scanner_sorted_ascii_ascending_case_insensitive()
    test_scanner_deduped_case_insensitive()

    # §B: fallback semantics (3)
    test_fallback_scanner_raises_llm_succeeds()
    test_fallback_both_fail_graceful()
    test_fallback_scanner_succeeds_no_llm_needed()

    # §C: recall metrics (3)
    test_recall_perfect_overlap()
    test_recall_zero_division_empty_llm()
    test_recall_zero_division_empty_regex()

    # §D: integration smoke (5)
    test_mock_helper_satisfies_scanner_invariants()
    test_marker_format_includes_new_recall_keys()
    test_partial_overlap_realistic_case()
    test_production_decision_structure_present()
    test_telemetry_independence_from_extracted_files()


if __name__ == "__main__":
    print("K_H7: Phase 2B activation tests (v119.H7)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H7 complete.")
