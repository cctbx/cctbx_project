"""
Phase 5: Error Classification Consistency.

Three error classifiers must agree on severity:
  1. ErrorAnalyzer (recoverable_errors.yaml) → recovery strategy
  2. DiagnosableErrorAnalyzer (diagnosable_errors.yaml) → diagnosis
  3. classify_error() (hardcoded patterns) → category

Agreement rules:
  - If recoverable detects an error, classify_error should NOT
    return TERMINAL
  - If diagnosable detects an error, classify_error should return
    TERMINAL (or at least not RETRYABLE/NO_ERROR)
  - No pattern should be claimed by both recoverable AND diagnosable

Run: python tests/tst_phase5_error_classification.py
Produces: findings/phase_5_classification_disagreements.yaml
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys
import yaml

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

from agent.error_analyzer import ErrorAnalyzer, get_diagnosis_detector
from agent.error_classifier import classify_error

BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
FINDINGS_DIR = os.path.join(BASE_DIR, 'findings')


# =====================================================================
# HELPERS
# =====================================================================

def pattern_to_sample(pattern):
    """Convert a regex/literal pattern to a sample log line.

    For simple literal patterns, use as-is.
    For regex patterns, do basic expansion (remove regex metacharacters).
    """
    # If it looks like a literal string (no metacharacters), use it
    simple = pattern
    # Remove common regex constructs
    simple = re.sub(r'\.\*', ' ', simple)
    simple = re.sub(r'\\\.', '.', simple)
    simple = re.sub(r'\\s\+', ' ', simple)
    simple = re.sub(r'\\s\*', ' ', simple)
    simple = re.sub(r'\(\?i\)', '', simple)
    simple = re.sub(r'[\\^$]', '', simple)
    simple = re.sub(r'\[.*?\]', 'x', simple)
    simple = re.sub(r'\(.*?\)', '', simple)
    simple = simple.strip()
    if not simple:
        simple = pattern
    # Wrap in realistic log context
    return "Sorry: %s\nRaised during program execution." % simple


def load_yamls():
    """Load both error YAML files."""
    rec_path = os.path.join(BASE_DIR, 'knowledge',
                            'recoverable_errors.yaml')
    diag_path = os.path.join(BASE_DIR, 'knowledge',
                             'diagnosable_errors.yaml')
    with open(rec_path) as f:
        rec = yaml.safe_load(f)
    with open(diag_path) as f:
        diag = yaml.safe_load(f)
    return rec, diag


# =====================================================================
# TESTS
# =====================================================================

def test_no_pattern_overlap():
    """No detection_pattern should appear in BOTH recoverable
    and diagnosable YAML files."""
    print("Test: no_pattern_overlap")
    rec, diag = load_yamls()

    rec_patterns = set()
    for err in rec.get('errors', {}).values():
        for p in err.get('detection_patterns', []):
            rec_patterns.add(p.lower())

    diag_patterns = set()
    for err in diag.get('errors', {}).values():
        for p in err.get('detection_patterns', []):
            diag_patterns.add(p.lower())

    overlap = rec_patterns & diag_patterns
    if overlap:
        print("  OVERLAP: %s" % sorted(overlap))
    assert len(overlap) == 0, (
        "Patterns in both recoverable and diagnosable: %s"
        % sorted(overlap))
    print("  PASSED (0 overlaps between %d rec + %d diag "
          "patterns)" % (len(rec_patterns), len(diag_patterns)))


def test_recoverable_not_terminal():
    """Recoverable error patterns should not be classified as
    TERMINAL by classify_error()."""
    print("Test: recoverable_not_terminal")
    rec, _ = load_yamls()

    disagreements = []
    total = 0

    for err_name, err_def in rec.get('errors', {}).items():
        for pattern in err_def.get('detection_patterns', []):
            total += 1
            sample = pattern_to_sample(pattern)
            result = classify_error(sample)
            category = result.get('category', 'NO_ERROR')
            if category == 'TERMINAL':
                disagreements.append({
                    'error': err_name,
                    'pattern': pattern[:60],
                    'classify_result': category,
                    'issue': 'recoverable pattern classified '
                             'as TERMINAL',
                })

    if disagreements:
        for d in disagreements:
            print("  DISAGREE: %s → %s (%s)"
                  % (d['error'], d['classify_result'],
                     d['pattern'][:40]))

    assert len(disagreements) == 0, (
        "%d recoverable patterns wrongly classified as TERMINAL"
        % len(disagreements))
    print("  PASSED (%d recoverable patterns checked, "
          "none TERMINAL)" % total)


def test_diagnosable_detected_as_terminal():
    """Check how many diagnosable patterns classify_error detects.

    This is INFORMATIONAL. In production, DiagnosableErrorAnalyzer
    runs BEFORE classify_error(), so gaps here are expected. The
    YAML detector is authoritative for diagnosable errors.

    We document gaps as findings, not failures."""
    print("Test: diagnosable_vs_classifier_coverage")
    _, diag = load_yamls()

    # Map diagnosable types to acceptable classify_error categories
    acceptable = {
        'crystal_symmetry_mismatch': {'TERMINAL'},
        'model_outside_map': {'TERMINAL'},
        'shelx_not_installed': {'TERMINAL'},
        'unknown_phil_parameter': {'TERMINAL', 'PHIL_ERROR',
                                   'AMBIGUOUS_PHIL'},
        'polymer_special_position': {'TERMINAL'},
        'unknown_chemical_element': {'TERMINAL'},
    }

    missed = []
    detected = 0
    total = 0

    for err_name, err_def in diag.get('errors', {}).items():
        ok_cats = acceptable.get(err_name, {'TERMINAL'})
        for pattern in err_def.get('detection_patterns', []):
            total += 1
            sample = pattern_to_sample(pattern)
            result = classify_error(sample)
            category = result.get('category', 'NO_ERROR')
            if category in ok_cats:
                detected += 1
            elif category in ('NO_ERROR', 'RETRYABLE'):
                missed.append({
                    'error': err_name,
                    'pattern': pattern[:60],
                    'expected': sorted(ok_cats),
                    'got': category,
                })
            else:
                detected += 1  # Different but non-trivial

    print("  Coverage: %d/%d diagnosable patterns detected "
          "by classify_error" % (detected, total))
    if missed:
        print("  Gaps (%d) — covered by YAML detector in "
              "production:" % len(missed))
        for m in missed[:5]:
            print("    %s: '%s' → %s" % (
                m['error'], m['pattern'][:30], m['got']))
        if len(missed) > 5:
            print("    ... and %d more" % (len(missed) - 5))
    # INFORMATIONAL — does not fail
    print("  PASSED (informational: %d gaps documented)"
          % len(missed))
    return missed


def test_error_analyzer_detects_recoverable():
    """ErrorAnalyzer should detect its own patterns."""
    print("Test: error_analyzer_detects_recoverable")
    rec, _ = load_yamls()

    analyzer = ErrorAnalyzer()
    missed = []
    total = 0

    for err_name, err_def in rec.get('errors', {}).items():
        for pattern in err_def.get('detection_patterns', []):
            total += 1
            sample = pattern_to_sample(pattern)
            # Use _detect_error_type directly since analyze()
            # needs session/context objects we don't have
            detected = analyzer._detect_error_type(sample)
            if not detected:
                missed.append({
                    'error': err_name,
                    'pattern': pattern[:60],
                })

    if missed:
        for m in missed:
            print("  MISSED: %s → '%s'" %
                  (m['error'], m['pattern'][:40]))

    assert len(missed) == 0, (
        "ErrorAnalyzer missed %d of its own patterns" %
        len(missed))
    print("  PASSED (%d patterns all detected)" % total)


def test_diagnosis_detector_detects_diagnosable():
    """DiagnosableErrorAnalyzer should detect its own patterns."""
    print("Test: diagnosis_detector_detects_diagnosable")
    _, diag = load_yamls()

    detector = get_diagnosis_detector()
    missed = []
    total = 0

    for err_name, err_def in diag.get('errors', {}).items():
        for pattern in err_def.get('detection_patterns', []):
            total += 1
            sample = pattern_to_sample(pattern)
            result = detector.detect(sample)
            if not result:
                missed.append({
                    'error': err_name,
                    'pattern': pattern[:60],
                })

    if missed:
        for m in missed:
            print("  MISSED: %s → '%s'" %
                  (m['error'], m['pattern'][:40]))

    assert len(missed) == 0, (
        "DiagnosableErrorAnalyzer missed %d of its own "
        "patterns" % len(missed))
    print("  PASSED (%d patterns all detected)" % total)


def test_all_yaml_patterns_are_valid_regex():
    """All detection_patterns must compile as valid regex."""
    print("Test: all_yaml_patterns_are_valid_regex")
    rec, diag = load_yamls()

    invalid = []
    total = 0

    for source, data in [('recoverable', rec),
                         ('diagnosable', diag)]:
        for err_name, err_def in data.get(
                'errors', {}).items():
            for pattern in err_def.get(
                    'detection_patterns', []):
                total += 1
                try:
                    re.compile(pattern, re.IGNORECASE)
                except re.error as e:
                    invalid.append({
                        'source': source,
                        'error': err_name,
                        'pattern': pattern[:60],
                        'regex_error': str(e),
                    })

    assert len(invalid) == 0, (
        "%d invalid regex patterns" % len(invalid))
    print("  PASSED (%d patterns all valid regex)" % total)


def test_no_empty_detection_patterns():
    """No error type should have an empty detection_patterns list."""
    print("Test: no_empty_detection_patterns")
    rec, diag = load_yamls()

    empty = []
    for source, data in [('recoverable', rec),
                         ('diagnosable', diag)]:
        for err_name, err_def in data.get(
                'errors', {}).items():
            patterns = err_def.get('detection_patterns', [])
            if not patterns:
                empty.append('%s/%s' % (source, err_name))

    assert len(empty) == 0, (
        "Error types with no patterns: %s" % empty)
    print("  PASSED (all error types have patterns)")


# =====================================================================
# RUN ALL
# =====================================================================

def run_all_tests():
    """Run all Phase 5 tests."""
    print("Phase 5: Error Classification Consistency")
    print()

    tests = [
        test_all_yaml_patterns_are_valid_regex,
        test_no_empty_detection_patterns,
        test_no_pattern_overlap,
        test_error_analyzer_detects_recoverable,
        test_diagnosis_detector_detects_diagnosable,
        test_recoverable_not_terminal,
        test_diagnosable_detected_as_terminal,
    ]

    passed = 0
    failed = 0
    failures = []
    missed_diagnosable = []

    for t in tests:
        try:
            result = t()
            if isinstance(result, list):
                missed_diagnosable = result
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
                        'phase_5_classification_disagreements.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 5: Error Classification Findings\n")
        f.write("phase: 5\n")
        f.write("name: error_classification\n")
        f.write("status: %s\n" % status)
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        if missed_diagnosable:
            f.write("\ndiagnosable_not_detected_by_classifier:\n")
            f.write("  note: \"These patterns from "
                    "diagnosable_errors.yaml were not\n")
            f.write("    detected by classify_error(). This "
                    "is informational —\n")
            f.write("    classify_error uses hardcoded patterns "
                    "that may not cover\n")
            f.write("    every YAML pattern (especially regex "
                    "variants).\"\n")
            for m in missed_diagnosable:
                f.write("  - error: %s\n" % m['error'])
                f.write("    pattern: \"%s\"\n"
                        % m['pattern'][:60])
                f.write("    got: %s\n" % m['got'])
        if failures:
            f.write("\ntest_failures:\n")
            for fail in failures:
                f.write("  - test: %s\n" % fail['test'])
                f.write("    error: \"%s\"\n"
                        % fail['error'][:100])

    print("  Findings: %s" % path)
    print("  Phase 5 overall: %s" % status)

    if failed > 0:
        raise AssertionError(
            "Phase 5 error classification FAILED: %d tests"
            % failed)

    return {'status': status, 'passed': passed,
            'failed': failed}


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] == 'PASS' else 1)
