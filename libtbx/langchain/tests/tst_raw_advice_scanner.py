"""K_H4: Raw advice scanner regression tests.

v119.H4.1 (Step 1F).  Phase 2B prerequisite.

Tests the `scan_files_in_advice` function in
agent/raw_advice_scanner.py — the regex-based file detector that
the metric block in run_advice_preprocessing uses to compute
the regex_recall metric for Phase 2B trigger evaluation.

Test organization:
  §A   Scanner basic detection             (4 tests)
  §B   Scanner edge cases                  (4 tests)
  §C   Scanner pattern robustness          (3 tests)
  §D   file_list_hint UNION + path norm    (3 tests)
  §E   Robustness (never raises)           (1 test)
  §F   ReDoS performance (REQUIRED)        (1 test)
  §G   H4.1: new extensions                (8 tests)
  §H   H4.1: boundary fixes                (3 tests)
  §H2  H4.1: additional context patterns   (3 tests)
  §I   H4.1: golden master corpus          (1 test)

Total: 31 tests (16 from H4 baseline + 15 new in H4.1).

Unlike the v119 H-cluster K-suites, K_H4 does NOT need a
graceful-skip pattern.  The scanner module is stdlib-only (re,
os, time, json), so all 31 tests run and PASS in both sandbox
and PHENIX.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import time


# ---------- Import helper ------------------------------------------

def _try_import_scanner():
    """Import scan_files_in_advice; return (fn, error_msg).

    Returns (callable, None) on success; (None, msg) on failure.
    The scanner is stdlib-only, so import should succeed in any
    environment.  We use the try-import pattern for consistency
    with the v119 K-suite style only.
    """
    try:
        from libtbx.langchain.agent.raw_advice_scanner import (
            scan_files_in_advice)
        return scan_files_in_advice, None
    except ImportError:
        pass
    try:
        # When run from langchain/ directory (sibling import)
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from agent.raw_advice_scanner import scan_files_in_advice
        return scan_files_in_advice, None
    except ImportError as e:
        return None, str(e)


# =====================================================================
# §5.1.1: Scanner basic detection (4)
# =====================================================================

def test_scanner_finds_pdb():
    """A simple advice with a .pdb file finds the PDB."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("refine model.pdb")
    assert result == ["model.pdb"], (
        "Expected ['model.pdb'], got %r" % result)
    print("  PASS: test_scanner_finds_pdb")


def test_scanner_finds_mtz():
    """A simple advice with a .mtz file finds the MTZ."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("using data.mtz")
    assert result == ["data.mtz"], (
        "Expected ['data.mtz'], got %r" % result)
    print("  PASS: test_scanner_finds_mtz")


def test_scanner_finds_multiple_extensions():
    """Advice mentioning multiple file types finds all of them."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("refine model.pdb with data.mtz and density.map")
    expected = sorted(["data.mtz", "density.map", "model.pdb"],
                      key=str.lower)
    assert result == expected, (
        "Expected %r, got %r" % (expected, result))
    print("  PASS: test_scanner_finds_multiple_extensions")


def test_scanner_finds_half_maps():
    """Half-map filenames are caught via .map and .mrc extensions."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("apply half_map_1.map and half1.mrc")
    assert "half_map_1.map" in result, (
        "Expected 'half_map_1.map' in result, got %r" % result)
    assert "half1.mrc" in result, (
        "Expected 'half1.mrc' in result, got %r" % result)
    print("  PASS: test_scanner_finds_half_maps")


# =====================================================================
# §5.1.2: Scanner edge cases (4)
# =====================================================================

def test_scanner_empty_advice():
    """Empty string input returns []."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    assert fn("") == []
    print("  PASS: test_scanner_empty_advice")


def test_scanner_no_files_mentioned():
    """Advice with no filename patterns returns []."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    assert fn("refine and stop") == []
    print("  PASS: test_scanner_no_files_mentioned")


def test_scanner_sorted_and_deduped():
    """Repeated filenames are deduplicated; output is sorted."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("data.mtz, model.pdb, data.mtz")
    assert result == ["data.mtz", "model.pdb"], (
        "Expected ['data.mtz', 'model.pdb'], got %r" % result)
    print("  PASS: test_scanner_sorted_and_deduped")


def test_scanner_case_insensitive_dedup():
    """Same file with different casings collapses to one entry."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("MODEL.PDB and model.pdb")
    assert len(result) == 1, (
        "Expected 1 entry (case-insensitive dedup), got %r" % result)
    # First-seen casing wins
    assert result[0].lower() == "model.pdb"
    print("  PASS: test_scanner_case_insensitive_dedup (%s)"
          % result[0])


# =====================================================================
# §5.1.3: Scanner pattern robustness (3)
# =====================================================================

def test_scanner_handles_path_in_advice():
    """File mentioned with an absolute path: basename is extracted."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("refine /abs/path/to/model.pdb")
    assert result == ["model.pdb"], (
        "Expected ['model.pdb'] (basename), got %r" % result)
    print("  PASS: test_scanner_handles_path_in_advice")


def test_scanner_handles_surrounding_punctuation():
    """Files with surrounding parens, quotes, etc. are caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn('refine (model.pdb) using "data.mtz"')
    assert "model.pdb" in result, (
        "Expected 'model.pdb' in result, got %r" % result)
    assert "data.mtz" in result, (
        "Expected 'data.mtz' in result, got %r" % result)
    print("  PASS: test_scanner_handles_surrounding_punctuation")


def test_scanner_no_false_positive_on_abbreviations():
    """Common abbreviations like 'e.g.' must not be caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    # 'e.g.' has extension 'g' which isn't in the allowlist
    result = fn("e.g., refine and stop")
    assert result == [], (
        "Abbreviation should not be caught; got %r" % result)
    # Also: 'etc.' — extension is empty, doesn't match
    result = fn("refine, etc., and stop")
    assert result == [], (
        "Trailing abbreviation should not be caught; got %r" % result)
    print("  PASS: test_scanner_no_false_positive_on_abbreviations")


# =====================================================================
# §5.1.4: file_list_hint UNION + path normalization (3)
# =====================================================================

def test_scanner_unions_with_hint():
    """file_list_hint is unioned into the result."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("refine model.pdb", file_list_hint=["data.mtz"])
    assert result == ["data.mtz", "model.pdb"], (
        "Expected ['data.mtz', 'model.pdb'], got %r" % result)
    print("  PASS: test_scanner_unions_with_hint")


def test_scanner_basenames_from_hint_paths():
    """Hint paths are reduced to basenames."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("refine", file_list_hint=["/abs/path/to/data.mtz"])
    assert result == ["data.mtz"], (
        "Expected ['data.mtz'] (basename), got %r" % result)
    print("  PASS: test_scanner_basenames_from_hint_paths")


def test_scanner_normalizes_hint_paths():
    """Different path forms of the same file collapse to one entry.

    Per Gemini's path-normalization guardrail: '/data/./refine_1.mtz'
    and 'data//refine_1.mtz' should normalize to the same basename
    via os.path.normpath + basename.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn(
        "",
        file_list_hint=[
            "/data/./refine_1.mtz",
            "data//refine_1.mtz",
            "refine_1.mtz",
        ])
    assert result == ["refine_1.mtz"], (
        "Expected ['refine_1.mtz'] (all paths normalize to same "
        "basename), got %r" % result)
    print("  PASS: test_scanner_normalizes_hint_paths")


# =====================================================================
# §5.1.5: Robustness — never raises (1)
# =====================================================================

def test_scanner_never_raises():
    """Pathological inputs return [] rather than raising."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return

    # raw_advice variants
    for bad_input in [
            None,
            42,
            3.14,
            ["model.pdb"],         # list instead of string
            {"file": "model.pdb"}, # dict
            b"model.pdb",          # bytes
            "\x00\x01\x02",        # control chars
            "\udcff\udcfe",        # invalid unicode surrogates
    ]:
        try:
            result = fn(bad_input)
            assert isinstance(result, list), (
                "Expected list for input %r, got %r"
                % (bad_input, type(result)))
        except Exception as e:
            raise AssertionError(
                "scan_files_in_advice(%r) raised %s: %s"
                % (bad_input, type(e).__name__, e))

    # file_list_hint variants
    try:
        fn("refine model.pdb", file_list_hint=None)
        fn("refine model.pdb", file_list_hint=[])
        fn("refine model.pdb", file_list_hint=[None, "", "data.mtz"])
        fn("refine model.pdb", file_list_hint=["data.mtz", 42, None])
        fn("refine model.pdb", file_list_hint=42)  # not even iterable
    except Exception as e:
        raise AssertionError(
            "scan_files_in_advice with bad hint raised %s: %s"
            % (type(e).__name__, e))

    print("  PASS: test_scanner_never_raises")


# =====================================================================
# §5.1.6: ReDoS performance (REQUIRED — Gemini Q3)
# =====================================================================

def test_scanner_handles_pathological_input():
    """100KB pathological input completes in <50ms.

    Per Gemini's Q3 elevation: this perf test is REQUIRED to
    prevent regex modifications from introducing ReDoS patterns
    that would block the server's hot path.  The current regex
    is non-backtracking by construction (see plan §3.1.1); this
    test asserts that property holds against future changes.

    The pathological input is designed to stress backtracking-
    prone regex engines: dense near-misses, punctuation runs,
    long pseudo-paths, unicode/control chars.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return

    # Build a ~100KB pathological input.
    near_misses = "foo.mtzz " * 5000          # ~45KB
    dense_punct = ",.;:()[];,.,.;" * 2000     # ~30KB
    long_paths = ("/a" * 100 + "/file.pdb ") * 50  # ~15KB
    unicode_chars = "\u0394foo \u03b1\u03b2\x00\x01" * 500
    payload = near_misses + dense_punct + long_paths + unicode_chars

    # Warm up once (regex compile + first-run overhead).
    fn(payload)

    # Time over multiple runs to reduce variance.
    n_runs = 5
    t0 = time.time()
    for _ in range(n_runs):
        result = fn(payload)
    elapsed_ms = (time.time() - t0) / n_runs * 1000.0

    # 50ms ceiling per plan §5.1.6.  In practice should be <10ms
    # on the non-backtracking regex.
    assert elapsed_ms < 50.0, (
        "Scanner took %.2fms on %d-char pathological input "
        "(threshold: 50ms).  Possible ReDoS regression — has the "
        "regex been modified to introduce backtracking?"
        % (elapsed_ms, len(payload)))

    # Sanity: the scanner should still produce a meaningful result
    # (the file.pdb tokens embedded in long_paths should be caught).
    assert "file.pdb" in result, (
        "Expected 'file.pdb' in result from pathological payload; "
        "got %r" % result[:5])

    print("  PASS: test_scanner_handles_pathological_input "
          "(%.2fms on %d chars)" % (elapsed_ms, len(payload)))


# =====================================================================
# Runner
# =====================================================================

# =====================================================================
# §G: H4.1 new extensions (8 tests)
# =====================================================================
# Per-extension confirmation that the expanded allowlist works.
# Selected from telemetry analysis (plan rev 3 §1.1).

def test_scanner_finds_inp_files():
    """`.inp` (PHENIX legacy input param) files caught.

    Top H4 miss class (+72 occurrences across telemetry).
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("rebuild.inp and ligand.inp")
    assert "rebuild.inp" in result and "ligand.inp" in result, (
        "Expected .inp files, got %r" % result)
    print("  PASS: test_scanner_finds_inp_files")


def test_scanner_finds_com_scripts():
    """`.com` script files caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("generate_mad.com generate_mir.com")
    assert "generate_mad.com" in result and "generate_mir.com" in result
    print("  PASS: test_scanner_finds_com_scripts")


def test_scanner_finds_cv_data():
    """`.cv` cross-validation reflection files caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("input data: lysozyme_scale.cv for cross-validation")
    assert "lysozyme_scale.cv" in result, (
        "Expected lysozyme_scale.cv, got %r" % result)
    print("  PASS: test_scanner_finds_cv_data")


def test_scanner_finds_ncs_spec():
    """`.ncs_spec` (8-character extension) caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("apply 3j89_6123.ncs_spec for NCS")
    assert "3j89_6123.ncs_spec" in result, (
        "Expected .ncs_spec, got %r" % result)
    print("  PASS: test_scanner_finds_ncs_spec")


def test_scanner_finds_chimerax_scenes():
    """`.cxs` ChimeraX scene files caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("see helix101-107.cxs in illustrations")
    assert "helix101-107.cxs" in result, (
        "Expected .cxs, got %r" % result)
    print("  PASS: test_scanner_finds_chimerax_scenes")


def test_scanner_finds_xplor_map():
    """`.xplor` map format caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("use if5a.xplor as input map")
    assert "if5a.xplor" in result
    print("  PASS: test_scanner_finds_xplor_map")


def test_scanner_finds_compound_extensions():
    """`.tar.gz` matched as one token via greedy filename body.

    The `[\\w.\\-]+` is greedy so `mr_rosetta.tar.gz` matches with
    the literal `.gz` extension and the captured token includes
    `.tar.gz` as part of the filename body.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("Files: mr_rosetta.tar.gz contains archive")
    assert "mr_rosetta.tar.gz" in result, (
        "Expected compound extension capture, got %r" % result)
    print("  PASS: test_scanner_finds_compound_extensions")


def test_scanner_finds_param_files():
    """`.param` restraint parameter files caught."""
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("Apply curated_ss.param for SS restraints")
    assert "curated_ss.param" in result
    print("  PASS: test_scanner_finds_param_files")


# =====================================================================
# §H: H4.1 boundary fixes (3 tests)
# =====================================================================
# Regression tests for the lookbehind/lookahead boundary change.

def test_scanner_handles_space_adjacent_files():
    """Two filenames separated by ONE space — both caught.

    H4's consuming-boundary regex caught only the first.  H4.1's
    non-consuming lookarounds catch both.  This is the most
    impactful single fix — ~24 misses (.ccp4, .pdb, .seq, .mtz)
    fixed by this change alone.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("Files: half_map_1.ccp4 half_map_2.ccp4 needed")
    assert "half_map_1.ccp4" in result and "half_map_2.ccp4" in result, (
        "Both half-maps should be caught (space-separated); got %r"
        % result)
    # Also: longer sequence
    result = fn("Files: beta.pdb blip.pdb beta.seq blip.seq beta_blip.mtz")
    for expected in ("beta.pdb", "blip.pdb", "beta.seq", "blip.seq"):
        assert expected in result, (
            "Expected %s in space-separated list; got %r"
            % (expected, result))
    print("  PASS: test_scanner_handles_space_adjacent_files")


def test_scanner_handles_markdown_bold():
    """Markdown-bold filenames `**file.pdb**` caught.

    The `**` markers immediately abut the filename with no
    boundary character between.  H4 missed these; H4.1's
    lookarounds handle them correctly because `*` is not in
    `[\\w.\\-]`.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("Select **2cn3.pdb** as the input model")
    assert "2cn3.pdb" in result, (
        "Expected 2cn3.pdb from markdown-bold, got %r" % result)
    print("  PASS: test_scanner_handles_markdown_bold")


def test_scanner_handles_sentence_ending_period():
    """Filename followed by sentence-ending period caught.

    Critical edge case: text like `from foo.eff.\\nNote...` puts
    a period AFTER the file extension as sentence-ending
    punctuation.  H4's lookahead `(?![\\w.\\-])` would fail here
    because `.` IS in that set.  H4.1's lookahead `(?![\\w\\-])`
    EXCLUDES the period, so this case works.

    Tradeoff: filename like `foo.eff.bak` matches `foo.eff` (the
    second dot becomes a non-issue).  Acceptable since `.bak` and
    similar suffix files aren't real workflow inputs.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("adding proper SS from high_resolution_ss.eff.\nNote stats")
    assert "high_resolution_ss.eff" in result, (
        "Expected file with sentence-end period; got %r" % result)
    # Also: end of input
    result = fn("output is foo.pdb.")
    assert "foo.pdb" in result
    print("  PASS: test_scanner_handles_sentence_ending_period")


# =====================================================================
# §H2: H4.1 additional context patterns (3 tests)
# =====================================================================
# Edge cases found during corpus analysis but not part of the
# main fixes.  Documents currently-working behavior so future
# changes don't break it unintentionally.

def test_scanner_handles_markdown_table_cell():
    """Filenames in markdown table cells caught.

    READMEs occasionally use markdown tables like
    `| input.pdb | model file |`.  The `|` characters abut the
    filename but are not in `[\\w.\\-]`, so the lookarounds
    correctly pass them as boundaries.

    Sandbox-confirmed during H4.1 self-review on the golden
    master corpus.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("| input.pdb | model file |")
    assert "input.pdb" in result, (
        "Expected markdown table cell capture, got %r" % result)
    # Also: full row of multiple files
    result = fn("| 1xyz.pdb | 1xyz.mtz | 1xyz.seq |")
    for expected in ("1xyz.pdb", "1xyz.mtz", "1xyz.seq"):
        assert expected in result
    print("  PASS: test_scanner_handles_markdown_table_cell")


def test_scanner_handles_backtick_wrapped_filenames():
    """Filenames wrapped in markdown backticks caught.

    READMEs frequently mention filenames as `` `model.pdb` ``
    using inline code formatting.  Backticks abut the filename
    with no whitespace.  Lookarounds correctly handle this
    because `` ` `` is not in `[\\w.\\-]`.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    result = fn("Use `model.pdb` as the input model")
    assert "model.pdb" in result, (
        "Expected backtick-wrapped capture, got %r" % result)
    print("  PASS: test_scanner_handles_backtick_wrapped_filenames")


def test_scanner_documented_precision_tradeoff_on_bak_suffix():
    """Documented behavior: `foo.pdb.bak` matches `foo.pdb`.

    The lookahead `(?![\\w\\-])` deliberately excludes period
    (`.`) so sentence-ending punctuation doesn't break matches.
    Side effect: when a filename like `foo.pdb` is followed by
    a backup-style suffix `.bak`, `.old`, etc., the regex
    captures the `foo.pdb` portion.

    This is an INTENTIONAL precision tradeoff.  Phase 2B's
    threshold formula measures recall (LLM ∩ regex), not
    precision (regex ∩ LLM-or-real), so this only adds noise
    to `in_regex_only` (which doesn't affect the trigger
    metric).  The LLM also doesn't list `.bak` files as inputs
    so they don't appear in `in_llm_only` either.

    This test pins the behavior so a future change that
    "tightens up" the lookahead is forced to re-confirm the
    tradeoff is still acceptable.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return
    # foo.pdb captured from foo.pdb.bak
    result = fn("old backup: foo.pdb.bak")
    assert "foo.pdb" in result, (
        "Documented behavior changed: foo.pdb.bak no longer "
        "matches foo.pdb (got %r).  If this is intentional, "
        "update the test and re-verify golden master recall."
        % result)
    # And foo.pdb captured from foo.pdb.old
    result = fn("previous run: foo.pdb.old")
    assert "foo.pdb" in result
    print("  PASS: test_scanner_documented_precision_tradeoff_on_bak_suffix")


# =====================================================================
# §I: H4.1 golden master corpus test (1 test)
# =====================================================================
# THE critical test for regression prevention.  Runs the scanner
# against 567 entries extracted from the May 2026 tutorial sweep
# and asserts aggregate recall ≥ 0.95 PLUS per-provider floor ≥
# 0.92.  Future scanner changes can be evaluated in ~3 seconds
# without PHENIX install or re-sweep.

def test_scanner_recall_against_golden_master():
    """Recall on golden master corpus must be ≥ 0.95 aggregate
    AND ≥ 0.92 per provider.

    The fixture `tests/step_1f_corpus.json` contains
    (raw_advice, llm_files) pairs extracted from production
    tutorial runs (run_36_quick_* sweep, May 23 2026).  567
    entries total, ~435 with non-empty llm_files (the rest are
    cases where the LLM extracted zero files).

    This test catches recall regressions caused by accidental
    scanner changes.  Two thresholds:
      - Aggregate: 0.95 (leaves 0.05 margin above Phase 2B's
        0.90 trigger; at H4.1 ship time, measured 0.981)
      - Per-provider: 0.92 (catches asymmetric regressions
        where one provider's recall collapses while others
        compensate; at H4.1 ship time, all three providers
        above 0.96)

    Per-provider floor is essential: an aggregate threshold
    alone could hide a serious regression on one provider that
    other providers happen to compensate for.

    The corpus is loaded via an absolute-path anchor to the test
    file location (per Gemini's rev 3 review §2) so the test
    runs correctly regardless of CWD.
    """
    fn, err = _try_import_scanner()
    if fn is None:
        print("  SKIP: cannot import scanner (%s)" % err)
        return

    # Absolute-path-anchored fixture load (Gemini guardrail).
    test_dir = os.path.dirname(os.path.abspath(__file__))
    corpus_path = os.path.join(test_dir, "step_1f_corpus.json")
    if not os.path.isfile(corpus_path):
        print("  SKIP: golden master corpus not found at %s" % corpus_path)
        return

    import json
    try:
        with open(corpus_path, "r", encoding="utf-8") as f:
            entries = json.load(f)
    except Exception as e:
        print("  SKIP: could not load corpus: %s" % e)
        return

    recalls = []
    by_provider = {}
    for e in entries:
        llm = e.get("llm_files") or []
        if not llm:
            continue
        regex = fn(e.get("raw_advice", ""))
        llm_lower = {f.lower() for f in llm}
        regex_lower = {f.lower() for f in regex}
        r = len(llm_lower & regex_lower) / len(llm_lower)
        recalls.append(r)
        prov = e.get("provider", "unknown")
        by_provider.setdefault(prov, []).append(r)

    aggregate = sum(recalls) / len(recalls) if recalls else 0.0
    per_provider = {
        p: sum(v) / len(v) for p, v in by_provider.items() if v
    }
    provider_summary = ", ".join(
        "%s=%.3f" % (p, per_provider[p]) for p in sorted(per_provider))

    # Aggregate floor: 0.95
    assert aggregate >= 0.95, (
        "Golden master AGGREGATE recall %.4f < 0.95 (N=%d).  "
        "Has the scanner regressed?  Per-provider: %s"
        % (aggregate, len(recalls), provider_summary))

    # Per-provider floor: 0.92 (catches asymmetric regressions)
    for prov, val in sorted(per_provider.items()):
        assert val >= 0.92, (
            "Golden master recall for provider '%s' is %.4f < 0.92 "
            "(aggregate=%.4f, full breakdown: %s).  An asymmetric "
            "regression like this can be hidden behind the "
            "aggregate threshold."
            % (prov, val, aggregate, provider_summary))

    print("  PASS: test_scanner_recall_against_golden_master "
          "(aggregate=%.4f, N=%d, %s)"
          % (aggregate, len(recalls), provider_summary))


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: Basic detection (4)
    test_scanner_finds_pdb()
    test_scanner_finds_mtz()
    test_scanner_finds_multiple_extensions()
    test_scanner_finds_half_maps()
    # §B: Edge cases (4)
    test_scanner_empty_advice()
    test_scanner_no_files_mentioned()
    test_scanner_sorted_and_deduped()
    test_scanner_case_insensitive_dedup()
    # §C: Pattern robustness (3)
    test_scanner_handles_path_in_advice()
    test_scanner_handles_surrounding_punctuation()
    test_scanner_no_false_positive_on_abbreviations()
    # §D: UNION + path normalization (3)
    test_scanner_unions_with_hint()
    test_scanner_basenames_from_hint_paths()
    test_scanner_normalizes_hint_paths()
    # §E: Robustness (1)
    test_scanner_never_raises()
    # §F: ReDoS performance (1, REQUIRED)
    test_scanner_handles_pathological_input()
    # §G: H4.1 new extensions (8)
    test_scanner_finds_inp_files()
    test_scanner_finds_com_scripts()
    test_scanner_finds_cv_data()
    test_scanner_finds_ncs_spec()
    test_scanner_finds_chimerax_scenes()
    test_scanner_finds_xplor_map()
    test_scanner_finds_compound_extensions()
    test_scanner_finds_param_files()
    # §H: H4.1 boundary fixes (3)
    test_scanner_handles_space_adjacent_files()
    test_scanner_handles_markdown_bold()
    test_scanner_handles_sentence_ending_period()
    # §H2: H4.1 additional context patterns (3)
    test_scanner_handles_markdown_table_cell()
    test_scanner_handles_backtick_wrapped_filenames()
    test_scanner_documented_precision_tradeoff_on_bak_suffix()
    # §I: H4.1 golden master corpus (1)
    test_scanner_recall_against_golden_master()


if __name__ == "__main__":
    run_all_tests()
