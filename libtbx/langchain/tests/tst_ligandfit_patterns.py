#!/usr/bin/env python
"""
Tests for phenix.ligandfit log_parsing patterns in programs.yaml.

Run with: python tests/tst_ligandfit_patterns.py

The previous patterns targeted 'Correlation Coefficient:' and
'Real-space R:', neither of which phenix.ligandfit writes; both real
logs extracted zero metrics. The real summary block is:

    Overall score: 69.04596
    Atoms placed: 20 of 40
    CC to map: 0.777

Fixtures are verbatim excerpts from two real runs, kept long enough to
include the candidate 'SCORE= ... CC=' lines and the raw 'Starting CC'
lines that precede the summary -- the text most likely to steal a match
from an under-anchored pattern. If PHENIX_LIGANDFIT_LOGS points at a
directory of complete .log files, the corpus test also runs against
those; otherwise it skips.
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


EXPECTED_KEYS = set([
    "ligand_cc", "ligand_score", "atoms_placed", "atoms_total",
    "n_ligands", "starting_cc", "starting_local_cc",
])

# ligandfit_5_bromodomain.log. A PARTIAL fit: half the ligand was never
# placed. Note the two different starting CCs (global 0.714, local 0.24)
# and the candidate SCORE=/CC= lines that must not be mistaken for the
# final summary.
LOG_PARTIAL_FIT = """
Ligand file has 40 non-H atoms
Estimated ligand volume: 135.0  A**3
Starting CC of ligand as input to map: 0.714
Starting local CC of ligand as input to map: 0.24
Try:  1 quick
 SCORE=    65.243 CC=     0.692 LIGANDS=       1  LIG=ligand_fit_1_pose_1.pdb
 SCORE=    69.046 CC=     0.777 LIGANDS=       1  LIG=ligand_fit_1_pose_2.pdb
 SCORE=    63.886 CC=     0.563 LIGANDS=       1  LIG=ligand_fit_1_pose_3.pdb
Running refinement with phenix.real_space_refine macro_cycles=2 weight=10
Overall score: 69.04596
Atoms placed: 20 of 40
CC to map: 0.777

This fit is the new best one...
 solution_number  1
 score  69.04596
 cc_overall  0.777
 starting_cc  0.714
 number_of_ligands  1
 cc_list  [0.777]
 fit_score_list  [69.04596]
 atoms_placed_list  [20]
 atoms_total_list  [40]
"""

# ligandfit_nsf-d2-ligand.log. A COMPLETE fit from a very poor start;
# the local starting CC is genuinely negative here.
LOG_COMPLETE_FIT = """
Starting CC of ligand as input to map: 0.069
Starting local CC of ligand as input to map: -0.099
 SCORE=    97.361 CC=     0.697 LIGANDS=       1  LIG=ligand_fit_1_pose_1.pdb
 SCORE=   103.693 CC=     0.844 LIGANDS=       1  LIG=ligand_fit_1_pose_2.pdb
Overall score: 96.98679
Atoms placed: 31 of 31
CC to map: 0.847

This fit is the new best one...
 score  96.98679
 cc_overall  0.847
 starting_cc  0.069
 number_of_ligands  1
"""

# A refine log excerpt: patterns must not fire on another program's log.
LOG_OTHER_PROGRAM = """
Start R-work = 0.3257, R-free = 0.3372
Final R-work = 0.2370, R-free = 0.2603
| r_work= 0.2369 r_free= 0.2601 coordinate error: 0.17 A |
"""


def _extract(log_text):
    """Extract through the real loader, so the test exercises
    programs.yaml rather than duplicating the regexes."""
    from libtbx.langchain.knowledge.metric_patterns import (
        extract_metrics_for_program)
    metrics = extract_metrics_for_program(log_text, "phenix.ligandfit")
    return dict((k, v) for k, v in (metrics or {}).items() if v is not None)


def test_patterns_are_defined():
    """programs.yaml defines exactly the metrics phenix.ligandfit
    writes -- no more, no fewer."""
    print("Test: ligandfit patterns are defined")

    from libtbx.langchain.knowledge.metric_patterns import (
        get_all_metric_patterns)

    patterns = get_all_metric_patterns()
    assert "phenix.ligandfit" in patterns, "no ligandfit patterns"

    defined = set(patterns["phenix.ligandfit"])
    assert defined == EXPECTED_KEYS, \
        "defined %r, expected %r" % (sorted(defined), sorted(EXPECTED_KEYS))

    print("  PASSED")


def test_partial_fit_values():
    """A run that placed only half the ligand."""
    print("Test: partial fit extracts expected values")

    got = _extract(LOG_PARTIAL_FIT)

    assert set(got) == EXPECTED_KEYS, "unexpected key set: %r" % sorted(got)
    assert abs(got["ligand_cc"] - 0.777) < 1e-6, got
    assert abs(got["ligand_score"] - 69.04596) < 1e-6, got
    assert got["atoms_placed"] == 20, got
    assert got["atoms_total"] == 40, got
    assert got["n_ligands"] == 1, got
    assert abs(got["starting_cc"] - 0.714) < 1e-6, got
    assert abs(got["starting_local_cc"] - 0.24) < 1e-6, got

    print("  PASSED")


def test_complete_fit_values():
    """A run that placed every atom from a poor starting point."""
    print("Test: complete fit extracts expected values")

    got = _extract(LOG_COMPLETE_FIT)

    assert set(got) == EXPECTED_KEYS, "unexpected key set: %r" % sorted(got)
    assert abs(got["ligand_cc"] - 0.847) < 1e-6, got
    assert abs(got["ligand_score"] - 96.98679) < 1e-6, got
    assert got["atoms_placed"] == 31, got
    assert got["atoms_total"] == 31, got

    print("  PASSED")


def test_two_starting_ccs_are_distinguished():
    """phenix.ligandfit reports a global and a local starting CC, and
    they can differ sharply (0.714 vs 0.24). The JSON summary's
    starting_cc equals the GLOBAL value; conflating them would misstate
    how much the fit improved."""
    print("Test: global and local starting CC distinguished")

    got = _extract(LOG_PARTIAL_FIT)
    assert abs(got["starting_cc"] - 0.714) < 1e-6, got
    assert abs(got["starting_local_cc"] - 0.24) < 1e-6, got
    assert got["starting_cc"] != got["starting_local_cc"]

    print("  PASSED")


def test_negative_local_starting_cc():
    """The local starting CC is genuinely negative in the nsf-d2 run.
    A character class without a minus sign truncates it silently."""
    print("Test: negative local starting CC is captured")

    got = _extract(LOG_COMPLETE_FIT)
    assert abs(got["starting_local_cc"] + 0.099) < 1e-6, got

    print("  PASSED")


def test_placement_fraction_is_derivable():
    """The tool emits atoms_placed and atoms_total; the placement
    fraction is derived from them by the caller, not extracted."""
    print("Test: placement fraction derivable from extracted metrics")

    partial = _extract(LOG_PARTIAL_FIT)
    complete = _extract(LOG_COMPLETE_FIT)

    assert partial["atoms_placed"] < partial["atoms_total"], partial
    assert complete["atoms_placed"] == complete["atoms_total"], complete

    frac = partial["atoms_placed"] / float(partial["atoms_total"])
    assert abs(frac - 0.5) < 1e-6, frac

    print("  PASSED")


def test_candidate_lines_do_not_steal_the_match():
    """Candidate poses print 'SCORE= ... CC= ...' before the summary.
    The final values must come from the summary block, not from the
    last candidate line."""
    print("Test: candidate SCORE=/CC= lines do not steal the match")

    got = _extract(LOG_PARTIAL_FIT)
    # last candidate CC is 0.563; the summary CC is 0.777
    assert abs(got["ligand_cc"] - 0.777) < 1e-6, got
    # last candidate SCORE is 63.886; the summary score is 69.04596
    assert abs(got["ligand_score"] - 69.04596) < 1e-6, got

    print("  PASSED")


def test_prefixed_phrase_does_not_match():
    """An unanchored pattern would match a sentence merely containing
    the phrase."""
    print("Test: prefixed phrase does not match")

    got = _extract("Best CC to map: 0.9\nPrevious Overall score: 12.0\n")
    assert "ligand_cc" not in got, got
    assert "ligand_score" not in got, got

    print("  PASSED")


def test_malformed_numbers_rejected():
    """The float grammar accepts only syntactically valid numbers."""
    print("Test: malformed numbers rejected")

    for bad in ("CC to map: 1.2.3\n", "Overall score: abc\n",
                "CC to map: 0.5 extra\n"):
        got = _extract(bad)
        assert got == {}, "%r extracted %r" % (bad, got)

    print("  PASSED")


def test_repeated_summary_takes_last():
    """The summary block can repeat -- the log says 'This fit is the new
    best one' -- so the LAST block is the final result. If a future
    change reversed this, every other test would still pass."""
    print("Test: repeated summary takes the last block")

    doubled = (
        "Overall score: 10.0\nAtoms placed: 5 of 40\nCC to map: 0.300\n"
        "This fit is the new best one...\n"
        "Overall score: 69.04596\nAtoms placed: 20 of 40\n"
        "CC to map: 0.777\n")
    got = _extract(doubled)

    assert abs(got["ligand_cc"] - 0.777) < 1e-6, got
    assert abs(got["ligand_score"] - 69.04596) < 1e-6, got
    assert got["atoms_placed"] == 20, got

    print("  PASSED")


def test_cc_primary_beats_fallback():
    """Both 'CC to map:' and 'cc_overall' appear in a real log. Pin which
    one wins, so a loader change cannot silently reverse it."""
    print("Test: CC to map takes precedence over cc_overall")

    got = _extract("CC to map: 0.777\n cc_overall  0.123\n")
    assert abs(got["ligand_cc"] - 0.777) < 1e-6, got

    print("  PASSED")


def test_cc_fallback_used_when_primary_absent():
    """Only reachable if a run omits 'CC to map:'. Not observed in the
    corpus; the fallback is defensive."""
    print("Test: cc_overall fallback when CC to map is absent")

    got = _extract(" cc_overall  0.612\nOverall score: 50.0\n")
    assert abs(got["ligand_cc"] - 0.612) < 1e-6, got

    print("  PASSED")


def test_no_false_positives_on_other_programs():
    """The patterns must not fire on another program's log."""
    print("Test: no false positives on other programs")

    got = _extract(LOG_OTHER_PROGRAM)
    assert got == {}, "ligandfit patterns matched a refine log: %r" % got

    print("  PASSED")


def test_complete_real_logs():
    """Regression against complete real logs, which contain everything
    the fixtures crop out. Skipped unless PHENIX_LIGANDFIT_LOGS points
    at a directory of ligandfit*.log files."""
    print("Test: complete real logs")

    log_dir = os.environ.get("PHENIX_LIGANDFIT_LOGS")
    if not log_dir or not os.path.isdir(log_dir):
        print("  SKIPPED (set PHENIX_LIGANDFIT_LOGS to enable)")
        return

    logs = [f for f in sorted(os.listdir(log_dir))
            if f.startswith("ligandfit") and f.endswith(".log")]
    assert logs, "no ligandfit logs in %s" % log_dir

    for name in logs:
        path = os.path.join(log_dir, name)
        with open(path) as handle:
            text = handle.read()
        got = _extract(text)
        assert set(got) == EXPECTED_KEYS, "%s: got %r" % (name, sorted(got))
        assert 0 < got["atoms_placed"] <= got["atoms_total"], \
            "%s: %r" % (name, got)
        assert -1.0 <= got["ligand_cc"] <= 1.0, "%s: %r" % (name, got)

    print("  PASSED (%d logs)" % len(logs))


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
