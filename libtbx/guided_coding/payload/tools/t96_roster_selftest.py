# Self-test for t96_roster.py, the parser for
# phenix_regression.test_all_parallel / libtbx.run_tests_parallel output.
#
# The suite exits 0 while reporting failures, so exit status is never a
# verdict; the roster of NAMED tests and their outcomes is.  This file
# checks that parse_roster()/compare_rosters() read the log TEXT the way
# parallel.py writes it, and that each known trap (a replay block that
# prints failures a second time, bracketed status words inside captured
# output or inside a command, outcome words in file names, a declared line
# carrying the " [Parallel]" suffix, absent Summary lines) is handled.
#
# compare_rosters() decides blocking by OUTCOME RANK rather than by
# membership in the named keys: passed (rank 2) over ran-but-did-not-pass
# (rank 1) over did-not-run or absent (rank 0), a command blocking when its
# candidate rank falls below its baseline rank, or below 2 when it is new.
#
# Fixtures are synthetic text built here: no log file, no network, no
# libtbx or phenix import.  Runs offline in well under a second, under
# libtbx.python and under plain python3.
#
# Deliberately NOT named tst_*: the regression harness harvests every
# tst* file under the module, and this self-test is not a suite test.

from __future__ import absolute_import, division, print_function

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t96_roster

EQ80 = "=" * 80
STARS80 = "*" * 80

def as_text(lines):
  """Join fixture lines into log text, with the trailing newline a real log has."""
  return "\n".join(lines) + "\n"

def mentions(items, word):
  """Return True if any string in items contains word, ignoring case."""
  return any(word.lower() in item.lower() for item in items)

def names_both(items, first, second):
  """Return True if one string in items contains both numbers as written."""
  return any(str(first) in item and str(second) in item for item in items)

def cmd_for(name):
  """Return the command line parallel.py prints for a regression test file."""
  return 'libtbx.python "/net/x/regression/%s"' % name

def exercise_constants():
  """The eight outcome constants exist with the documented string values."""
  assert t96_roster.PASS == 'pass'
  assert t96_roster.FAIL == 'fail'
  assert t96_roster.WARNING == 'warning'
  assert t96_roster.EXPECTED_FAIL == 'expected_fail'
  assert t96_roster.EXPECTED_UNSTABLE == 'expected_unstable'
  assert t96_roster.SKIPPED == 'skipped'
  assert t96_roster.MISSING == 'missing'
  assert t96_roster.TIMED_OUT == 'timed_out'

def exercise_replay_block_not_double_counted():
  """A failure printed twice (streaming, then replay) is one test, not two."""
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_b = 'libtbx.python "/net/x/regression/tst_b.py"'
  cmd_c = 'libtbx.python "/net/x/regression/tst_c.py"'
  fail_line = '%s [FAIL] 3.4s  (failed after 3 attempts)' % cmd_b
  text = as_text([
    'Running 3 tests on 2 processors:',
    '  ' + cmd_a,
    '  ' + cmd_b,
    '  ' + cmd_c,
    '',
    '%s [OK] 1.2s' % cmd_a,
    fail_line,
    '  Time:  3.40',
    '  Return code: 1',
    '  OKs: 0',
    '%s [OK] 0.9s' % cmd_c,
    EQ80,
    '',
    'Tests finished. Elapsed time: 5.50s',
    '',
    'Error: the following jobs returned non-zero exit codes or suspicious stderr output:',
    '',
    fail_line,
    '  Time:  3.40',
    '  Return code: 1',
    '  OKs: 0',
    '',
    'Please verify these tests manually.',
    '',
    'Summary:',
    '  Tests run                    : 3',
    '  Failures                     : 1',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Retries used                 : 2 extra attempts across 1 tests',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.declared_count == 3
  assert r.nprocs == '2'
  assert sorted(r.declared_commands) == [cmd_a, cmd_b, cmd_c]
  assert len(r.tests) == 3, sorted(r.tests)
  counts = r.counts()
  assert counts[t96_roster.FAIL] == 1, counts
  assert counts[t96_roster.PASS] == 2, counts
  assert sum(counts.values()) == 3, counts
  assert r.names(t96_roster.FAIL) == [cmd_b]
  assert cmd_b in r.replay_commands, r.replay_commands
  assert r.tests[cmd_b].attempt == 3, r.tests[cmd_b].attempt
  assert r.summary['tests_run'] == 3
  assert r.summary['failures'] == 1
  assert r.summary['stderr_output'] == 0
  assert r.reconcile() == [], r.reconcile()

def exercise_indented_status_bracket_is_not_a_result():
  """Status brackets inside a passing test's captured output are not outcomes."""
  cmd = 'libtbx.python "/net/x/regression/tst_recover.py"'
  inner = 'libtbx.python "/net/x/regression/tst_inner.py"'
  text = as_text([
    'Running 1 tests on 1 processors:',
    '  ' + cmd,
    '',
    '%s [OK] 4.1s' % cmd,
    '  Time:  4.10',
    '  Return code: 0',
    '  OKs: 3',
    '  Standard out:',
    '    Starting macro-cycle 1',
    '    [WARNING] RECOVERY LIMIT REACHED',
    '    %s [FAIL] 0.3s' % inner,
    '    OK',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert len(r.tests) == 1, sorted(r.tests)
  assert r.tests[cmd].outcome == t96_roster.PASS
  assert r.counts()[t96_roster.WARNING] == 0, r.counts()
  assert r.counts()[t96_roster.FAIL] == 0, r.counts()
  assert r.summary['warnings'] == 0
  assert not r.not_all_finished

def exercise_outcome_word_in_file_name():
  """A file name containing 'timeout' does not make an outcome."""
  cmd = 'libtbx.python "/net/x/mmtbx/regression/tst_reduce_timeout.py"'
  text = as_text([
    'Running 1 tests on 1 processors:',
    '  ' + cmd,
    '',
    '%s [OK] 8.0s' % cmd,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.tests[cmd].outcome == t96_roster.PASS
  counts = r.counts()
  assert counts[t96_roster.TIMED_OUT] == 0, counts
  for outcome in (t96_roster.PASS, t96_roster.FAIL, t96_roster.WARNING,
                  t96_roster.EXPECTED_FAIL, t96_roster.EXPECTED_UNSTABLE,
                  t96_roster.SKIPPED, t96_roster.MISSING,
                  t96_roster.TIMED_OUT):
    assert outcome in counts, outcome
  assert len(r.tests) == 1, sorted(r.tests)
  assert sum(counts.values()) == 1, counts

def exercise_non_python_command():
  """A bare .csh command is a test like any other, keyed by its full text."""
  cmd_py = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_csh = '/net/x/phenix_regression/cctbx/tst_development.csh'
  text = as_text([
    'Running 2 tests on 1 processors:',
    '  ' + cmd_py,
    '  ' + cmd_csh,
    '',
    '%s [OK] 1.0s' % cmd_py,
    '%s [OK] 2.0s' % cmd_csh,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 2',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert cmd_csh in r.tests, sorted(r.tests)
  assert r.tests[cmd_csh].outcome == t96_roster.PASS
  assert r.tests[cmd_csh].wall_time == 2.0
  assert r.counts()[t96_roster.MISSING] == 0, r.counts()
  assert len(r.tests) == 2, sorted(r.tests)

def exercise_command_containing_bracketed_text():
  """The split is on the trailing status, not on bracketed text in the command."""
  cmd = 'libtbx.python "/net/x/regression/tst_bracket [OK] 1.0s.py"'
  text = as_text([
    'Running 1 tests on 1 processors:',
    '  ' + cmd,
    '',
    '%s [WARNING] 2.5s' % cmd,
    '  Time:  2.50',
    '  Return code: 0',
    '  OKs: 1',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 1',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert len(r.tests) == 1, sorted(r.tests)
  assert cmd in r.tests, sorted(r.tests)
  assert r.tests[cmd].outcome == t96_roster.WARNING
  assert r.tests[cmd].wall_time == 2.5
  assert r.counts()[t96_roster.PASS] == 0, r.counts()
  assert r.counts()[t96_roster.MISSING] == 0, r.counts()

def exercise_preamble_is_ignored():
  """Nothing before the 'Running N tests' header is read, whatever its shape.

  run_tests_parallel.py prints driver chatter such as "Keeping the test
  <command>" before run_command_list ever prints its header
  (libtbx/command_line/run_tests_parallel.py line 261), and a log can hold
  other text in front of the run.  Two preamble shapes change the answer
  unless the parser anchors on the header: a RESULT-SHAPED line at column 0,
  which would enter the roster as a test of its own, and a line of eighty
  '=' characters, which is how parallel.py ends a streaming section (line
  510) and would close this one before it opened, leaving every declared
  test MISSING.
  """
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_b = 'libtbx.python "/net/x/regression/tst_b.py"'
  ghost = 'libtbx.python "/net/x/regression/tst_ghost.py"'
  text = as_text([
    'Keeping the test %s' % ghost,
    '%s [OK] 6.0s' % ghost,
    '  ' + ghost,
    EQ80,
    'Running 2 tests on 1 processors:',
    '  ' + cmd_a,
    '  ' + cmd_b,
    '',
    '%s [OK] 1.0s' % cmd_a,
    '%s [FAIL] 2.0s' % cmd_b,
    '  Time:  2.00',
    '  Return code: 1',
    '  OKs: 0',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 2',
    '  Failures                     : 1',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.declared_count == 2, r.declared_count
  assert sorted(r.declared_commands) == sorted([cmd_a, cmd_b]), \
    sorted(r.declared_commands)
  # The preamble's result-shaped line did not become a test of its own ...
  assert ghost not in r.tests, sorted(r.tests)
  assert sorted(r.tests) == sorted([cmd_a, cmd_b]), sorted(r.tests)
  # ... and the preamble's '=' line did not close the streaming section, so
  # both real result lines were still read.
  assert r.tests[cmd_a].outcome == t96_roster.PASS, r.tests[cmd_a].outcome
  assert r.tests[cmd_b].outcome == t96_roster.FAIL, r.tests[cmd_b].outcome
  assert r.counts()[t96_roster.MISSING] == 0, r.counts()

def exercise_parallel_suffix_is_not_part_of_command():
  """' [Parallel]' marks a declared line; the command is the rest of that line."""
  # parallel.py prints a parallel test's declared line as "  %s [Parallel]" % cmd
  # (libtbx/test_utils/parallel.py).  The suffix is not part of the command, so
  # the declared line has to match the plain command on the later result line.
  cmd_plain = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_par = 'libtbx.python "/net/x/regression/tst_parallel_job.py"'
  text = as_text([
    'Running 2 tests on 2 processors:',
    '  ' + cmd_plain,
    '  %s [Parallel]' % cmd_par,
    '',
    '%s [OK] 1.0s' % cmd_plain,
    '%s [OK] 9.0s' % cmd_par,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 2',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.declared_count == 2
  assert sorted(r.declared_commands) == sorted([cmd_plain, cmd_par]), \
    sorted(r.declared_commands)
  assert sorted(r.tests) == sorted([cmd_plain, cmd_par]), sorted(r.tests)
  assert r.tests[cmd_par].outcome == t96_roster.PASS
  assert r.tests[cmd_par].wall_time == 9.0
  # The whole point of stripping: declared matches result, so nothing is MISSING.
  assert r.counts()[t96_roster.MISSING] == 0, r.counts()
  assert not r.not_all_finished
  assert not mentions(sorted(r.declared_commands), '[Parallel]'), \
    sorted(r.declared_commands)
  assert not mentions(sorted(r.tests), '[Parallel]'), sorted(r.tests)

def exercise_missing_test():
  """A declared test with no result line is MISSING, and that is reported."""
  # 'Tests run' is self.finished, the number of RESULTS (parallel.py line 574
  # printing the value set at line 522), so it is 1 here and not the declared
  # 2: no log ever shows 'Tests run : 2' beside a single result line.  That
  # is exactly why the Summary alone cannot reveal an unfinished run, and why
  # the declared count has to be compared against it.
  cmd_done = 'libtbx.python "/net/x/regression/tst_done.py"'
  cmd_hung = 'libtbx.python "/net/x/regression/tst_hung.py"'
  text = as_text([
    'Running 2 tests on 2 processors:',
    '  ' + cmd_done,
    '  ' + cmd_hung,
    '',
    '%s [OK] 1.0s' % cmd_done,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    STARS80,
    '  WARNING: NOT ALL TESTS FINISHED!',
    STARS80,
    ])
  r = t96_roster.parse_roster(text)
  assert r.tests[cmd_hung].outcome == t96_roster.MISSING
  assert r.tests[cmd_done].outcome == t96_roster.PASS
  assert r.counts()[t96_roster.MISSING] == 1, r.counts()
  assert r.not_all_finished
  assert r.declared_count == 2, r.declared_count
  assert r.summary['tests_run'] == 1, r.summary
  # Every Summary figure agrees with the streamed results now, so the one
  # thing left to report is 2 tests declared against 1 test run.
  discrepancies = r.reconcile()
  assert discrepancies, discrepancies
  assert names_both(discrepancies, 2, 1), discrepancies

def exercise_duplicate_command_skip():
  """'repeated, skipping' is read where parallel.py prints it: before the header."""
  # run_command_list.__init__ filters duplicates out of cmd_list and prints
  # "Test %s repeated, skipping" while doing so (parallel.py line 454), which
  # is BEFORE it prints "Running %d tests on %s processors:" (line 464).  A
  # parser that looked for skip lines only from the header onward would never
  # see one, and every skipped test would silently vanish from the roster.
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_dup = 'libtbx.python "/x/tst_dup.py"'
  text = as_text([
    'Test %s repeated, skipping' % cmd_dup,
    'Running 1 tests on 1 processors:',
    '  ' + cmd_a,
    '',
    '%s [OK] 1.0s' % cmd_a,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert cmd_dup in r.tests, sorted(r.tests)
  assert r.tests[cmd_dup].outcome == t96_roster.SKIPPED
  assert r.names(t96_roster.SKIPPED) == [cmd_dup]
  assert r.tests[cmd_a].outcome == t96_roster.PASS

def exercise_skip_line_does_not_override_a_result():
  """A command with both a skip line and a result line keeps the result."""
  # De-duplication drops only the SECOND copy: the first stays in cmd_list
  # (parallel.py lines 449-454), is printed under the header and runs.  So the
  # ordinary shape of a skip line is a command that ALSO has a result line,
  # and the result is the truth - the test ran.  Letting the skip win would
  # score a passing test as SKIPPED, which is rank 0, and report an unchanged
  # test as blocking.
  #
  # cmd_gone carries a skip line and nothing else.  It is here so that this
  # exercise cannot be satisfied by a parser that simply never reads the skip
  # lines: such a parser gets cmd_twice right by accident, but loses cmd_gone.
  # Both commands are read from the SAME text, so cmd_twice's PASS is an
  # override of a skip line that was seen, not of one that was missed.
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_twice = 'libtbx.python "/net/x/regression/tst_twice.py"'
  cmd_gone = 'libtbx.python "/net/x/regression/tst_gone.py"'
  text = as_text([
    'Test %s repeated, skipping' % cmd_twice,
    'Test %s repeated, skipping' % cmd_gone,
    'Running 2 tests on 1 processors:',
    '  ' + cmd_a,
    '  ' + cmd_twice,
    '',
    '%s [OK] 1.0s' % cmd_a,
    '%s [OK] 5.0s' % cmd_twice,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 2',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  # The result wins: cmd_twice ran, and it passed.
  assert r.tests[cmd_twice].outcome == t96_roster.PASS, \
    r.tests[cmd_twice].outcome
  assert r.tests[cmd_twice].wall_time == 5.0, r.tests[cmd_twice].wall_time
  # The skip lines in this position were read all the same.
  assert cmd_gone in r.tests, sorted(r.tests)
  assert r.tests[cmd_gone].outcome == t96_roster.SKIPPED, \
    r.tests[cmd_gone].outcome
  assert r.names(t96_roster.SKIPPED) == [cmd_gone], r.names(t96_roster.SKIPPED)
  assert r.counts()[t96_roster.SKIPPED] == 1, r.counts()
  assert r.counts()[t96_roster.PASS] == 2, r.counts()
  assert len(r.tests) == 3, sorted(r.tests)
  # The consequence: an unchanged, passing test that happens to be listed
  # twice must not come back as a blocking skip.
  baseline = build_roster([(cmd_a, 'OK', 0), (cmd_twice, 'OK', 0)])
  candidate = build_roster([(cmd_a, 'OK', 0), (cmd_twice, 'OK', 0)],
    skipped=[cmd_twice])
  assert candidate.tests[cmd_twice].outcome == t96_roster.PASS, \
    candidate.tests[cmd_twice].outcome
  diff = t96_roster.compare_rosters(baseline, candidate)
  assert cmd_twice not in diff['newly_skipped'], diff['newly_skipped']
  assert list(diff['blocking']) == [], diff['blocking']

def exercise_reconcile_reports_failure_mismatch():
  """A Summary claiming more failures than were printed is reported, not raised."""
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_b = 'libtbx.python "/net/x/regression/tst_b.py"'
  text = as_text([
    'Running 2 tests on 1 processors:',
    '  ' + cmd_a,
    '  ' + cmd_b,
    '',
    '%s [OK] 1.0s' % cmd_a,
    '%s [FAIL] 2.0s' % cmd_b,
    '  Time:  2.00',
    '  Return code: 1',
    '  OKs: 0',
    EQ80,
    '',
    'Error: the following jobs returned non-zero exit codes or suspicious stderr output:',
    '',
    '%s [FAIL] 2.0s' % cmd_b,
    '  Time:  2.00',
    '  Return code: 1',
    '  OKs: 0',
    '',
    'Summary:',
    '  Tests run                    : 2',
    '  Failures                     : 2',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.counts()[t96_roster.FAIL] == 1, r.counts()
  assert r.summary['failures'] == 2
  discrepancies = r.reconcile()
  assert discrepancies, discrepancies
  assert mentions(discrepancies, 'failure'), discrepancies

def exercise_absent_summary_line_is_not_zero():
  """An omitted Summary line reads as None and is reported, not read as zero."""
  cmd = 'libtbx.python "/net/x/regression/tst_a.py"'
  text = as_text([
    'Running 1 tests on 1 processors:',
    '  ' + cmd,
    '',
    '%s [OK] 1.0s' % cmd,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.counts()[t96_roster.WARNING] == 0, r.counts()
  assert r.summary['warnings'] is None, r.summary
  assert r.summary['retries_used'] is None, r.summary
  assert r.summary['tests_run'] == 1
  discrepancies = r.reconcile()
  assert discrepancies, discrepancies
  assert mentions(discrepancies, 'warning'), discrepancies

def exercise_strip_profile_lines():
  """Leading PROFILE lines go; later ones stay; stripping is the caller's step."""
  cmd = 'libtbx.python "/net/x/regression/tst_a.py"'
  body = [
    'Running 1 tests on 1 processors:',
    '  ' + cmd,
    '',
    '%s [OK] 1.0s' % cmd,
    '  Standard out:',
    '    PROFILE',
    '    OK',
    'PROFILE',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 1',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 0',
    ]
  plain = "\n".join(body)
  assert t96_roster.strip_profile_lines(plain) == plain
  # body carries a PROFILE line indented inside captured output and another at
  # column 0 after the header; equality with body is what proves both survive.
  one = as_text(['PROFILE'] + body)
  kept = t96_roster.strip_profile_lines(one).splitlines()
  assert kept == body, kept[:4]
  two = as_text(['PROFILE', 'PROFILE'] + body)
  assert t96_roster.strip_profile_lines(two).splitlines() == body
  r = t96_roster.parse_roster(t96_roster.strip_profile_lines(two))
  assert r.declared_count == 1
  assert sorted(r.tests) == [cmd], sorted(r.tests)
  assert r.tests[cmd].outcome == t96_roster.PASS
  assert r.counts()[t96_roster.PASS] == 1, r.counts()

def exercise_retry_annotations_and_mode():
  """Retry notes give attempt/attempts_total, and mode names rc and attempts."""
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_b = 'libtbx.python "/net/x/regression/tst_b.py"'
  cmd_c = 'libtbx.python "/net/x/regression/tst_c.py"'
  cmd_d = 'libtbx.python "/net/x/regression/tst_d.py"'
  text = as_text([
    'Running 4 tests on 2 processors:',
    '  ' + cmd_a,
    '  ' + cmd_b,
    '  ' + cmd_c,
    '  ' + cmd_d,
    '',
    '%s [OK] 1.0s' % cmd_a,
    '%s [OK] 12.3s  (passed on attempt 2 of 3)' % cmd_b,
    '%s [FAIL] 4.5s  (failed after 3 attempts)' % cmd_c,
    '  Time:  4.50',
    '  Return code: 1',
    '  OKs: 0',
    '%s [FAIL] 0.5s' % cmd_d,
    '  Time:  0.50',
    '  Return code: 2',
    '  OKs: 0',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 4',
    '  Failures                     : 2',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Retries used                 : 3 extra attempts across 2 tests',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  a = r.tests[cmd_a]
  assert a.attempt == 1 and a.attempts_total == 1, (a.attempt, a.attempts_total)
  assert a.wall_time == 1.0
  assert a.mode is None, a.mode
  b = r.tests[cmd_b]
  assert b.outcome == t96_roster.PASS
  assert b.attempt == 2 and b.attempts_total == 3, (b.attempt, b.attempts_total)
  assert abs(b.wall_time - 12.3) < 1.e-6, b.wall_time
  assert b.mode is None, b.mode
  c = r.tests[cmd_c]
  assert c.outcome == t96_roster.FAIL
  assert c.attempt == 3 and c.attempts_total == 3, (c.attempt, c.attempts_total)
  assert int(c.return_code) == 1, c.return_code
  assert c.mode == 'rc=1; attempt 3 of 3', c.mode
  d = r.tests[cmd_d]
  assert d.attempt == 1, d.attempt
  assert d.mode == 'rc=2', d.mode
  assert r.summary['retries_used'] is not None, r.summary
  assert r.summary['retries_used'].strip() == '3 extra attempts across 2 tests', \
    r.summary['retries_used']

def exercise_summary_parenthetical():
  """The value after the colon is the count; the parenthetical is not."""
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_k = 'libtbx.python "/net/x/regression/tst_known.py"'
  cmd_u = 'libtbx.python "/net/x/regression/tst_unstable.py"'
  text = as_text([
    'Running 3 tests on 1 processors:',
    '  ' + cmd_a,
    '  ' + cmd_k,
    '  ' + cmd_u,
    '',
    '%s [OK] 1.0s' % cmd_a,
    '%s [EXPECTED FAIL] 2.0s' % cmd_k,
    '  Time:  2.00',
    '  Return code: 1',
    '  OKs: 0',
    '%s [EXPECTED UNSTABLE] 3.0s' % cmd_u,
    '  Time:  3.00',
    '  Return code: 1',
    '  OKs: 0',
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 3',
    '  Failures                     : 0',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  8)         : 7',
    '  Known Unstable (  3)         : 1',
    '  Stderr output (discouraged)  : 0',
    ])
  r = t96_roster.parse_roster(text)
  assert r.summary['known_failures'] == 7, r.summary['known_failures']
  assert r.summary['known_unstable'] == 1, r.summary['known_unstable']
  assert r.tests[cmd_k].outcome == t96_roster.EXPECTED_FAIL
  assert r.tests[cmd_u].outcome == t96_roster.EXPECTED_UNSTABLE
  assert r.counts()[t96_roster.FAIL] == 0, r.counts()

def build_roster(entries, skipped=(), missing=()):
  """Build a Roster from (command, status, return_code) triples.

  entries -- one (command, status, return_code) per test that prints a
    result line; status is the bracketed word parallel.py prints, one of
    'OK', 'FAIL', 'WARNING', 'EXPECTED FAIL', 'EXPECTED UNSTABLE'.
  skipped -- commands whose duplicate copy was dropped.  parallel.py prints
    "Test %s repeated, skipping" from the de-duplication loop in
    run_command_list.__init__ (line 454), BEFORE the header (line 464), so
    that is where this puts the line.
  missing -- commands declared under the header that never print a result
    line, which is how a hung test looks; parallel.py then ends with its
    NOT ALL TESTS FINISHED banner.

  'Tests run' is len(entries), never the declared count: parallel.py prints
  self.finished, the number of results (line 574, set at line 522).
  """
  commands = [entry[0] for entry in entries] + list(missing)
  lines = []
  for command in skipped:
    lines.append('Test %s repeated, skipping' % command)
  lines.append('Running %d tests on 1 processors:' % len(commands))
  for command in commands:
    lines.append('  ' + command)
  lines.append('')
  tally = {'FAIL': 0, 'WARNING': 0, 'EXPECTED FAIL': 0, 'EXPECTED UNSTABLE': 0}
  for (command, status, return_code) in entries:
    lines.append('%s [%s] 1.0s' % (command, status))
    if status != 'OK':
      assert status in tally, status
      tally[status] += 1
      lines.append('  Time:  1.00')
      lines.append('  Return code: %d' % return_code)
      lines.append('  OKs: 0')
  lines.extend([
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : %d' % len(entries),
    '  Failures                     : %d' % tally['FAIL'],
    '  Warnings (possible failures) : %d' % tally['WARNING'],
    '  Known Failures (% 3d)         : %d'
      % (tally['EXPECTED FAIL'], tally['EXPECTED FAIL']),
    '  Known Unstable (% 3d)         : %d'
      % (tally['EXPECTED UNSTABLE'], tally['EXPECTED UNSTABLE']),
    '  Stderr output (discouraged)  : 0',
    ])
  if missing:
    lines.extend([STARS80, '  WARNING: NOT ALL TESTS FINISHED!', STARS80])
  return t96_roster.parse_roster(as_text(lines))

def force_outcome(roster, command, outcome):
  """Replace one already-parsed command's outcome, and return the roster.

  Used only for TIMED_OUT.  parallel.py has no timeout and never prints such
  a status, so no log text can make parse_roster produce one; everything else
  about the entry stays exactly as parsed.
  """
  entry = roster.tests[command]
  try:
    entry.outcome = outcome
  except AttributeError:               # an immutable record type
    roster.tests[command] = entry._replace(outcome=outcome)
  assert roster.tests[command].outcome == outcome, \
    roster.tests[command].outcome
  return roster

def exercise_compare_rosters_transitions():
  """Transitions land in the right buckets, and blocking follows the rank rule."""
  cmd_a = 'libtbx.python "/net/x/regression/tst_a.py"'
  cmd_b = 'libtbx.python "/net/x/regression/tst_b.py"'
  cmd_c = 'libtbx.python "/net/x/regression/tst_c.py"'
  cmd_d = 'libtbx.python "/net/x/regression/tst_d.py"'
  cmd_e = 'libtbx.python "/net/x/regression/tst_e.py"'
  cmd_f = 'libtbx.python "/net/x/regression/tst_f.py"'
  baseline = build_roster([
    (cmd_a, 'OK', 0),
    (cmd_b, 'OK', 0),
    (cmd_c, 'FAIL', 1),
    (cmd_d, 'FAIL', 1),
    (cmd_f, 'OK', 0),
    ])
  candidate = build_roster([
    (cmd_a, 'FAIL', 1),
    (cmd_c, 'FAIL', 2),
    (cmd_d, 'OK', 0),
    (cmd_e, 'OK', 0),
    ], skipped=[cmd_f])
  assert baseline.tests[cmd_c].mode == 'rc=1', baseline.tests[cmd_c].mode
  assert candidate.tests[cmd_c].mode == 'rc=2', candidate.tests[cmd_c].mode
  diff = t96_roster.compare_rosters(baseline, candidate)
  assert cmd_a in diff['new_failures'], diff['new_failures']
  assert cmd_a in diff['blocking'], diff['blocking']
  assert cmd_b in diff['disappeared'], diff['disappeared']
  assert cmd_b in diff['blocking'], diff['blocking']
  assert (cmd_c, 'rc=1', 'rc=2') in diff['mode_changed'], diff['mode_changed']
  assert cmd_c not in diff['new_failures'], diff['new_failures']
  assert cmd_c not in diff['blocking'], diff['blocking']
  assert cmd_d in diff['fixed'], diff['fixed']
  assert cmd_d not in diff['blocking'], diff['blocking']
  assert cmd_e in diff['new_tests'], diff['new_tests']
  assert cmd_e not in diff['blocking'], diff['blocking']
  assert cmd_f in diff['newly_skipped'], diff['newly_skipped']
  assert cmd_f in diff['blocking'], diff['blocking']
  assert list(diff['newly_timed_out']) == [], diff['newly_timed_out']
  assert 'newly_missing' in diff, sorted(diff)
  # Every declared command here reports a result, and cmd_b is absent from the
  # candidate roster rather than MISSING in it, so nothing is newly missing.
  assert list(diff['newly_missing']) == [], diff['newly_missing']
  # Rank rule: cmd_a PASS->FAIL (2 to 1), cmd_b PASS->absent (2 to 0) and
  # cmd_f PASS->SKIPPED (2 to 0) block.  cmd_c FAIL->FAIL with a changed
  # return code (1 to 1), cmd_d FAIL->PASS (1 to 2) and cmd_e absent->PASS
  # (new, rank 2) do not.
  assert list(diff['blocking']) == sorted([cmd_a, cmd_b, cmd_f]), \
    diff['blocking']

def exercise_newly_missing_blocks():
  """A test that hangs on the candidate blocks; MISSING on both sides does not."""
  # parallel.py has no timeout, so a hung test prints no result line at all and
  # is only visible as declared-but-unrun.  Without newly_missing such a run
  # reaches no blocking key and blocking comes back empty while a test hangs.
  cmd_ok = 'libtbx.python "/net/x/regression/tst_ok.py"'
  cmd_hung = 'libtbx.python "/net/x/regression/tst_hung.py"'
  cmd_both = 'libtbx.python "/net/x/regression/tst_both.py"'
  cmd_added = 'libtbx.python "/net/x/regression/tst_added.py"'
  baseline = build_roster(
    [(cmd_ok, 'OK', 0), (cmd_hung, 'OK', 0)], missing=[cmd_both])
  candidate = build_roster(
    [(cmd_ok, 'OK', 0)], missing=[cmd_hung, cmd_both, cmd_added])
  # The fixtures say what they are meant to say before anything is compared.
  assert baseline.tests[cmd_hung].outcome == t96_roster.PASS
  assert baseline.tests[cmd_both].outcome == t96_roster.MISSING
  assert cmd_added not in baseline.tests, sorted(baseline.tests)
  assert candidate.tests[cmd_hung].outcome == t96_roster.MISSING
  assert candidate.tests[cmd_both].outcome == t96_roster.MISSING
  assert candidate.tests[cmd_added].outcome == t96_roster.MISSING
  assert candidate.not_all_finished
  diff = t96_roster.compare_rosters(baseline, candidate)
  # Passed, then no result line at all: the hang this key exists to catch.
  assert cmd_hung in diff['newly_missing'], diff['newly_missing']
  assert cmd_hung in diff['blocking'], diff['blocking']
  # MISSING on both sides is rank 0 to rank 0: not a transition, must not block.
  assert cmd_both not in diff['newly_missing'], diff['newly_missing']
  assert cmd_both not in diff['blocking'], diff['blocking']
  # Absent from the baseline roster counts as not MISSING there.
  assert cmd_added in diff['newly_missing'], diff['newly_missing']
  assert cmd_added in diff['blocking'], diff['blocking']
  # Exactly those two, sorted: nothing else is newly missing and nothing else
  # blocks.  Declaration order here is hung before added, so this also says
  # the list is sorted rather than in the order the roster was walked.
  expected = sorted([cmd_added, cmd_hung])
  assert list(diff['newly_missing']) == expected, diff['newly_missing']
  assert list(diff['blocking']) == expected, diff['blocking']

ABSENT = None   # sentinel in the transition table: not in that roster at all

def exercise_blocking_is_an_outcome_rank_rule():
  """blocking is a two-clause rule, not the union of the named keys.

  A command blocks when its candidate outcome ranks BELOW the rank it had to
  reach, OR when it FAILS on the candidate without having failed on the
  baseline.

  Rank 2 is PASS (ran and passed).  Rank 1 is FAIL, WARNING, EXPECTED_FAIL
  and EXPECTED_UNSTABLE (ran, did not pass).  Rank 0 is SKIPPED, MISSING,
  TIMED_OUT and absence from a roster (did not run).  Over the union of both
  rosters' commands the rank a command had to reach is its baseline rank, or
  2 when it is ABSENT FROM THE BASELINE.

  Both clauses are needed, and neither subsumes the other.

  Rank alone is too permissive: FAIL, WARNING, EXPECTED_FAIL and
  EXPECTED_UNSTABLE share rank 1, so coming back [FAIL] from SKIPPED,
  MISSING or TIMED_OUT (rank 0) reads as a rank INCREASE and would not
  block, though the candidate plainly fails.  Ranking FAIL below WARNING
  would not mend that: SKIPPED to FAIL is still a rise.

  new_failures alone is too narrow: a test that passed on the baseline and
  comes back [WARNING], [EXPECTED FAIL] or [EXPECTED UNSTABLE] never FAILS,
  so the second clause cannot see it; and it is in none of new_failures,
  disappeared, newly_skipped, newly_timed_out or newly_missing either, so
  under the old union-of-keys rule it reached no key at all and a regression
  was reported as clean.
  """
  # Every outcome below is produced from real-shaped log text, except
  # TIMED_OUT: parallel.py has no timeout and never prints such a status.
  cmd_p2f = cmd_for('tst_pass_to_fail.py')
  cmd_p2w = cmd_for('tst_pass_to_warning.py')
  cmd_p2ef = cmd_for('tst_pass_to_expected_fail.py')
  cmd_p2eu = cmd_for('tst_pass_to_expected_unstable.py')
  cmd_p2s = cmd_for('tst_pass_to_skipped.py')
  cmd_p2m = cmd_for('tst_pass_to_missing.py')
  cmd_p2t = cmd_for('tst_pass_to_timed_out.py')
  cmd_p2x = cmd_for('tst_pass_to_absent.py')
  cmd_f2m = cmd_for('tst_fail_to_missing.py')
  cmd_f2x = cmd_for('tst_fail_to_absent.py')
  cmd_x2f = cmd_for('tst_absent_to_fail.py')
  cmd_x2m = cmd_for('tst_absent_to_missing.py')
  cmd_f2p = cmd_for('tst_fail_to_pass.py')
  cmd_f2f = cmd_for('tst_fail_to_fail.py')
  cmd_m2m = cmd_for('tst_missing_to_missing.py')
  cmd_x2p = cmd_for('tst_absent_to_pass.py')
  cmd_w2p = cmd_for('tst_warning_to_pass.py')
  cmd_e2e = cmd_for('tst_expected_fail_both.py')
  # Second clause: the candidate FAILS and the baseline did not fail.  Every
  # one of these is a rank RISE or a rank tie, so the rank clause is blind to
  # all six of them.
  cmd_w2f = cmd_for('tst_warning_to_fail.py')
  cmd_ef2f = cmd_for('tst_expected_fail_to_fail.py')
  cmd_eu2f = cmd_for('tst_expected_unstable_to_fail.py')
  cmd_s2f = cmd_for('tst_skipped_to_fail.py')
  cmd_m2f = cmd_for('tst_missing_to_fail.py')
  cmd_t2f = cmd_for('tst_timed_out_to_fail.py')
  baseline = build_roster([
    (cmd_p2f, 'OK', 0),
    (cmd_p2w, 'OK', 0),
    (cmd_p2ef, 'OK', 0),
    (cmd_p2eu, 'OK', 0),
    (cmd_p2s, 'OK', 0),
    (cmd_p2m, 'OK', 0),
    (cmd_p2t, 'OK', 0),
    (cmd_p2x, 'OK', 0),
    (cmd_f2m, 'FAIL', 1),
    (cmd_f2x, 'FAIL', 1),
    (cmd_f2p, 'FAIL', 1),
    (cmd_f2f, 'FAIL', 1),
    (cmd_w2p, 'WARNING', 0),
    (cmd_e2e, 'EXPECTED FAIL', 1),
    (cmd_w2f, 'WARNING', 0),
    (cmd_ef2f, 'EXPECTED FAIL', 1),
    (cmd_eu2f, 'EXPECTED UNSTABLE', 1),
    ], skipped=[cmd_s2f], missing=[cmd_m2m, cmd_m2f, cmd_t2f])
  candidate = build_roster([
    (cmd_p2f, 'FAIL', 1),
    (cmd_p2w, 'WARNING', 0),
    (cmd_p2ef, 'EXPECTED FAIL', 1),
    (cmd_p2eu, 'EXPECTED UNSTABLE', 1),
    (cmd_f2p, 'OK', 0),
    (cmd_f2f, 'FAIL', 1),
    (cmd_w2p, 'OK', 0),
    (cmd_e2e, 'EXPECTED FAIL', 1),
    (cmd_x2f, 'FAIL', 1),
    (cmd_x2p, 'OK', 0),
    (cmd_w2f, 'FAIL', 1),
    (cmd_ef2f, 'FAIL', 1),
    (cmd_eu2f, 'FAIL', 1),
    (cmd_s2f, 'FAIL', 1),
    (cmd_m2f, 'FAIL', 1),
    (cmd_t2f, 'FAIL', 1),
    ], skipped=[cmd_p2s],
    missing=[cmd_p2m, cmd_p2t, cmd_f2m, cmd_m2m, cmd_x2m])
  force_outcome(baseline, cmd_t2f, t96_roster.TIMED_OUT)
  force_outcome(candidate, cmd_p2t, t96_roster.TIMED_OUT)
  # command, baseline outcome, candidate outcome, blocks
  table = [
    (cmd_p2f, t96_roster.PASS, t96_roster.FAIL, True),
    (cmd_p2w, t96_roster.PASS, t96_roster.WARNING, True),
    (cmd_p2ef, t96_roster.PASS, t96_roster.EXPECTED_FAIL, True),
    (cmd_p2eu, t96_roster.PASS, t96_roster.EXPECTED_UNSTABLE, True),
    (cmd_p2s, t96_roster.PASS, t96_roster.SKIPPED, True),
    (cmd_p2m, t96_roster.PASS, t96_roster.MISSING, True),
    (cmd_p2t, t96_roster.PASS, t96_roster.TIMED_OUT, True),
    (cmd_p2x, t96_roster.PASS, ABSENT, True),
    (cmd_f2m, t96_roster.FAIL, t96_roster.MISSING, True),
    (cmd_f2x, t96_roster.FAIL, ABSENT, True),
    (cmd_x2f, ABSENT, t96_roster.FAIL, True),
    (cmd_x2m, ABSENT, t96_roster.MISSING, True),
    (cmd_f2p, t96_roster.FAIL, t96_roster.PASS, False),
    (cmd_f2f, t96_roster.FAIL, t96_roster.FAIL, False),
    (cmd_m2m, t96_roster.MISSING, t96_roster.MISSING, False),
    (cmd_x2p, ABSENT, t96_roster.PASS, False),
    (cmd_w2p, t96_roster.WARNING, t96_roster.PASS, False),
    (cmd_e2e, t96_roster.EXPECTED_FAIL, t96_roster.EXPECTED_FAIL, False),
    # Second clause.  Rank 1 to rank 1 for the first three, rank 0 to rank 1
    # for the last three: not one of them is a rank DROP, so the rank clause
    # alone calls every one of them clean.
    (cmd_w2f, t96_roster.WARNING, t96_roster.FAIL, True),
    (cmd_ef2f, t96_roster.EXPECTED_FAIL, t96_roster.FAIL, True),
    (cmd_eu2f, t96_roster.EXPECTED_UNSTABLE, t96_roster.FAIL, True),
    (cmd_s2f, t96_roster.SKIPPED, t96_roster.FAIL, True),
    (cmd_m2f, t96_roster.MISSING, t96_roster.FAIL, True),
    (cmd_t2f, t96_roster.TIMED_OUT, t96_roster.FAIL, True),
    ]
  # The fixtures say what they are meant to say before anything is compared,
  # so a parse problem cannot be mistaken for a compare_rosters verdict.
  for (command, base_outcome, cand_outcome, blocks) in table:
    for (roster, outcome, side) in ((baseline, base_outcome, 'baseline'),
                                    (candidate, cand_outcome, 'candidate')):
      if outcome is ABSENT:
        assert command not in roster.tests, (side, command)
      else:
        assert command in roster.tests, (side, command, sorted(roster.tests))
        assert roster.tests[command].outcome == outcome, \
          (side, command, roster.tests[command].outcome, outcome)
  diff = t96_roster.compare_rosters(baseline, candidate)
  for (command, base_outcome, cand_outcome, blocks) in table:
    if blocks:
      assert command in diff['blocking'], (command, sorted(diff['blocking']))
    else:
      assert command not in diff['blocking'], \
        (command, sorted(diff['blocking']))
  # Sorted, and nothing beyond the two clauses blocks.
  assert list(diff['blocking']) == sorted(
    [entry[0] for entry in table if entry[3]]), diff['blocking']
  # What tells the two clauses apart, so neither can be dropped in silence.
  # FAIL->FAIL does not block even though the candidate fails: the second
  # clause must not fire when the baseline failed too.
  assert cmd_f2f not in diff['blocking'], diff['blocking']
  # PASS->WARNING does block even though the candidate does not FAIL: the
  # first clause must still fire where the second cannot reach.
  assert cmd_p2w in diff['blocking'], diff['blocking']
  # new_failures is the second clause's source, so all of it blocks.
  for command in diff['new_failures']:
    assert command in diff['blocking'], (command, sorted(diff['blocking']))
  # The named keys keep the meaning and membership they already had.
  assert cmd_p2f in diff['new_failures'], diff['new_failures']
  assert cmd_p2x in diff['disappeared'], diff['disappeared']
  assert cmd_p2s in diff['newly_skipped'], diff['newly_skipped']
  assert cmd_p2t in diff['newly_timed_out'], diff['newly_timed_out']
  assert cmd_p2m in diff['newly_missing'], diff['newly_missing']
  assert cmd_f2m in diff['newly_missing'], diff['newly_missing']
  assert cmd_f2p in diff['fixed'], diff['fixed']
  assert cmd_x2p in diff['new_tests'], diff['new_tests']
  assert 'mode_changed' in diff, sorted(diff)
  # newly_not_passing: in BOTH rosters, PASS on the baseline, not PASS now.
  assert 'newly_not_passing' in diff, sorted(diff)
  for command in (cmd_p2w, cmd_p2ef, cmd_p2f, cmd_p2eu):
    assert command in diff['newly_not_passing'], \
      (command, sorted(diff['newly_not_passing']))
  # Never passed on the baseline, so not a loss of passing.
  assert cmd_f2f not in diff['newly_not_passing'], diff['newly_not_passing']
  assert cmd_x2f not in diff['newly_not_passing'], diff['newly_not_passing']
  # Not in both rosters: absent from the candidate is 'disappeared', not this.
  assert cmd_p2x not in diff['newly_not_passing'], diff['newly_not_passing']
  assert list(diff['newly_not_passing']) == sorted(diff['newly_not_passing']), \
    diff['newly_not_passing']

def exercise_reconcile_declared_count_vs_tests_run():
  """Fewer results than declared tests is a discrepancy naming both numbers.

  parallel.py prints 'Tests run' as self.finished, the length of
  self.results (parallel.py line 574 printing the value set at line 522), so
  that line always equals the streaming result-line count even when tests
  never reported.  It can never disagree with the results, and on its own it
  cannot show an unfinished run.  The declared count from the header is the
  other number available, so declared against run is the comparison that
  notices one.
  """
  finished = [cmd_for('tst_done_%d.py' % i) for i in range(4)]
  unreported = [cmd_for('tst_hung_%d.py' % i) for i in range(3)]
  r = build_roster([(command, 'OK', 0) for command in finished],
    missing=unreported)
  assert r.declared_count == 7, r.declared_count
  assert r.summary['tests_run'] == 4, r.summary
  assert r.counts()[t96_roster.MISSING] == 3, r.counts()
  assert r.not_all_finished
  discrepancies = r.reconcile()
  assert discrepancies, discrepancies
  assert names_both(discrepancies, 7, 4), discrepancies
  # The complete run in exercise_replay_block_not_double_counted declares 3
  # and runs 3 and asserts reconcile() == [], so this check does not fire
  # when the counts agree.

def public_field_names(obj):
  """Return the non-callable public attribute names of an object as a set.

  Methods are excluded: counts() is how a field is reached, not a field.
  """
  names = set()
  for name in dir(obj):
    if name.startswith('_'):
      continue
    try:
      value = getattr(obj, name)
    except AttributeError:             # a property needing state this lacks
      continue
    if callable(value):
      continue
    names.add(name)
  return names

def canonical_key_template(field_origin, names):
  """Return the key spelling FIELD_ORIGIN uses for one group of field names.

  The spelling is discovered, never assumed: every mapping key containing
  one of the group's field names votes for the template got by replacing
  that occurrence with '%s', and the template spelling the most of the group
  wins.  A template written out here would go stale as the mapping can.

  field_origin -- the mapping, read for its keys only.
  names -- the field names of one group, taken off the live objects.
  Returns None when no key spells any name in the group.
  """
  votes = {}
  for key in field_origin:
    for name in names:
      start = key.rfind(name)
      if start < 0:
        continue
      template = key[:start] + '%s' + key[start + len(name):]
      if template.count('%') != 1:     # a stray % would break template % name
        continue
      votes[template] = votes.get(template, 0) + 1
  ranked = sorted(votes.items(), key=lambda item: (-item[1], item[0]))
  if not ranked:
    return None
  return ranked[0][0]

def exercise_field_origin_covers_every_reported_field():
  """FIELD_ORIGIN labels every reported field, and labels nothing that is gone.

  FIELD_ORIGIN is an enumeration, and enumerations go stale: a field added
  later carries no origin label and nothing notices.  So the field list is
  read off the live objects AT RUNTIME, and so is the key spelling -- a
  hand-written list of either in this file would be the same staleness bug
  one level up, passing forever while the reported surface drifts.
  """
  origins = set(t96_roster.ALL_ORIGINS)
  named = [t96_roster.ORIGIN_HARNESS_STATEMENT,
           t96_roster.ORIGIN_HARNESS_FORMATTED,
           t96_roster.ORIGIN_CANDIDATE_CLAIM]
  for origin in named:
    assert isinstance(origin, str), origin
    assert origin in origins, origin
  assert len(set(named)) == 3, named

  # A real-shaped log, so an attribute set only while parsing counts too.
  parsed = build_roster(
    [(cmd_for('tst_pass.py'), 'OK', 0),
     (cmd_for('tst_fail.py'), 'FAIL', 1),
     (cmd_for('tst_warn.py'), 'WARNING', 0),
     (cmd_for('tst_known.py'), 'EXPECTED FAIL', 1),
     (cmd_for('tst_unstable.py'), 'EXPECTED UNSTABLE', 1)],
    skipped=[cmd_for('tst_dup.py')],
    missing=[cmd_for('tst_hung.py')])
  try:
    fresh = t96_roster.Roster()
  except TypeError:                    # a constructor that takes the log text
    fresh = t96_roster.parse_roster('')

  roster_fields = set()
  summary_fields = set()
  counts_fields = set()
  for roster in (fresh, parsed):
    roster_fields |= public_field_names(roster)
    summary_fields |= set(roster.summary.keys())
    counts_fields |= set(roster.counts().keys())
  outcome_fields = set()
  for entry in parsed.tests.values():
    outcome_fields |= public_field_names(entry)
  counts_group = 'counts() key'
  outcome_group = 'TestOutcome attribute'
  groups = [('Roster attribute', roster_fields),
            ('summary key', summary_fields),
            (counts_group, counts_fields),
            (outcome_group, outcome_fields)]
  for (group, names) in groups:
    assert names, group

  field_origin = t96_roster.FIELD_ORIGIN
  templates = {}
  expected = {}
  for (group, names) in groups:
    template = canonical_key_template(field_origin, names)
    assert template is not None, \
      'no FIELD_ORIGIN key spells any %s; the whole group is unlabelled' \
      % group
    templates[group] = template
    for name in names:
      expected[template % name] = (group, name)

  # 1. Every reported field carries a label.
  for key in sorted(expected):
    (group, name) = expected[key]
    assert key in field_origin, \
      '%s %r has no entry in FIELD_ORIGIN (expected key %r)' \
      % (group, name, key)

  for key in sorted(field_origin):
    # 2. No label outlives the field it describes.
    assert key in expected, \
      'FIELD_ORIGIN labels %r, which is not a field of any reported object' \
      % key
    # 3. Every label is one of the known origins.
    assert field_origin[key] in origins, \
      'FIELD_ORIGIN[%r] is %r, which is not in ALL_ORIGINS' \
      % (key, field_origin[key])

  # 4. Nothing contributing to an outcome or a count is a bare candidate
  # claim: captured test output is displayed and never counted.
  outcome_values = set(t96_roster.ALL_OUTCOMES)
  counted = set(templates[counts_group] % name for name in counts_fields)
  for entry in parsed.tests.values():
    for name in public_field_names(entry):
      try:
        is_outcome = getattr(entry, name) in outcome_values
      except TypeError:                # an unhashable field, e.g. a list
        is_outcome = False
      if is_outcome:
        counted.add(templates[outcome_group] % name)
  assert counted, counted
  for key in sorted(counted):
    assert field_origin[key] != t96_roster.ORIGIN_CANDIDATE_CLAIM, \
      '%r contributes to an outcome or a count and cannot be labelled ' \
      'ORIGIN_CANDIDATE_CLAIM' % key

# ---------------------------------------------------------------------------
# Failure modes.  A test that fails on both sides with the same return code
# and attempts can still fail DIFFERENTLY, and the only record of how is its
# captured standard error, which display_result prints four-space indented
# under "  Standard error:" (parallel.py lines 707-709).  The exercises below
# are written from the contract in the plan of change
# 2026-09-13-roster-failure-modes (section 2), not from an implementation.
# ---------------------------------------------------------------------------

# Two frames of one traceback, and what the comparison makes of them.
FM_FRAMES = [('/net/x/regression/tst_mode.py', 30, '<module>'),
             ('/net/x/regression/tst_mode.py', 20, 'run')]
FM_CALL_PATH = [['tst_mode.py', '<module>'], ['tst_mode.py', 'run']]
# Frames of a cleanup that raised while the first exception was handled.
FM_CLEANUP_FRAMES = [('/net/x/regression/tst_mode.py', 40, 'run'),
                     ('/net/x/lib/python3.11/shutil.py', 701,
                      '_rmtree_safe_fd')]
FM_CLEANUP_ERROR = ("FileNotFoundError: [Errno 2] No such file or directory: "
                    "'meta.json.tmp'")
# The lines Python prints between the two tracebacks of a chained exception.
FM_CHAINED_SEPARATOR = [
  '', 'During handling of the above exception, another exception occurred:',
  '']

def result_block(entry):
  """Return the lines display_result prints for one result.

  entry -- (command, status, return_code, attempts, stderr_lines); see
    failure_log for what each field produces.
  """
  (command, status, return_code, attempts, stderr_lines) = entry
  note = ''
  if attempts > 1:
    if return_code == 0:
      note = '  (passed on attempt %d of %d)' % (attempts, attempts)
    else:
      note = '  (failed after %d attempts)' % attempts
  lines = ['%s [%s] 1.0s%s' % (command, status, note)]
  if status != 'OK':
    lines.extend(['  Time:  1.00', '  Return code: %d' % return_code,
                  '  OKs: 0'])
  if stderr_lines:
    lines.append('  Standard error:')
    lines.extend(['    ' + line for line in stderr_lines])
  return lines

def failure_log(entries, replay=None, skipped=(), missing=()):
  """Return suite log text whose results carry Standard error blocks.

  entries -- one (command, status, return_code, attempts, stderr_lines) per
    test that prints a result line.  status is the bracketed word
    parallel.py prints.  attempts > 1 adds the retry note display_result
    prints (parallel.py lines 686-691).  A status other than OK gets the
    Time / Return code / OKs lines (lines 697-701), which save_result asks
    for on any result that is not OK.  stderr_lines None or empty prints no
    Standard error block (line 707); otherwise every line is printed under
    "  Standard error:" with four spaces in front, as line 709 joins them,
    so an empty string prints as four spaces.
  replay -- entries printed a second time after the '=' rule, under the
    replay header.  None means what parallel.py itself prints (lines
    559-570): when any test FAILED, every WARNING result and then every
    FAIL result again, identically.
  skipped, missing -- as for build_roster.
  """
  commands = [entry[0] for entry in entries] + list(missing)
  lines = ['Test %s repeated, skipping' % command for command in skipped]
  lines.append('Running %d tests on 1 processors:' % len(commands))
  lines.extend(['  ' + command for command in commands])
  lines.append('')
  tally = {'OK': 0, 'FAIL': 0, 'WARNING': 0, 'EXPECTED FAIL': 0,
           'EXPECTED UNSTABLE': 0}
  for entry in entries:
    tally[entry[1]] += 1
    lines.extend(result_block(entry))
  if replay is None:
    replay = []
    if tally['FAIL']:
      replay = ([entry for entry in entries if entry[1] == 'WARNING']
                + [entry for entry in entries if entry[1] == 'FAIL'])
  lines.extend([EQ80, '', 'Tests finished. Elapsed time: 1.00s', '', '',
                'Warning: the following are 5 longest jobs:'])
  lines.extend(['  %s: 1.0s' % entry[0] for entry in entries[:5]])
  lines.extend(['Please try to reduce overall runtime - consider splitting '
                'up these tests.', ''])
  if replay:
    lines.extend(['', 'Error: the following jobs returned non-zero exit '
                  'codes or suspicious stderr output:', ''])
    for entry in replay:
      lines.extend(result_block(entry))
    lines.extend(['', 'Please verify these tests manually.', ''])
  retried = [entry for entry in entries if entry[3] > 1]
  lines.extend([
    'Summary:',
    '  Tests run                    : %d' % len(entries),
    '  Failures                     : %d' % tally['FAIL'],
    '  Warnings (possible failures) : %d' % tally['WARNING'],
    '  Known Failures (% 3d)         : %d'
      % (tally['EXPECTED FAIL'], tally['EXPECTED FAIL']),
    '  Known Unstable (% 3d)         : %d'
      % (tally['EXPECTED UNSTABLE'], tally['EXPECTED UNSTABLE']),
    ])
  if retried:
    lines.append('  Retries used                 : %d extra attempts across '
                 '%d tests' % (sum(entry[3] - 1 for entry in retried),
                               len(retried)))
  lines.append('  Stderr output (discouraged)  : %d'
               % len([entry for entry in entries if entry[4]]))
  if missing:
    lines.extend([STARS80, '  WARNING: NOT ALL TESTS FINISHED!', STARS80])
  return as_text(lines)

def traceback_text(frames, exception_line):
  """Return the lines of one Python traceback as the interpreter prints them.

  frames -- (path, line_number, function) per frame, outermost first.  Each
    frame line is followed by a source line and a caret line, which belong
    to the frame run without being frames.
  exception_line -- the last line, or None for a traceback that stops before
    it, which is a truncated one.
  """
  lines = ['Traceback (most recent call last):']
  for (path, line_number, function) in frames:
    lines.append('  File "%s", line %d, in %s' % (path, line_number, function))
    lines.append('    step()')
    lines.append('    ^^^^^^')
  if exception_line is not None:
    lines.append(exception_line)
  return lines

def failure_test(stderr_lines, status='FAIL', return_code=1, attempts=1,
                 replay=None):
  """Parse a one-test log and return that test's TestOutcome.

  stderr_lines, status, return_code, attempts -- as for failure_log.
  replay -- None for the replay copy parallel.py prints itself (identical,
    and only for a FAIL); otherwise a (return_code, stderr_lines) pair
    printed as the replay copy instead.
  """
  command = cmd_for('tst_mode.py')
  replay_entries = None
  if replay is not None:
    replay_entries = [(command, status, replay[0], attempts, replay[1])]
  roster = t96_roster.parse_roster(failure_log(
    [(command, status, return_code, attempts, stderr_lines)],
    replay=replay_entries))
  return roster.tests[command]

def exercise_failure_mode_change_in_sub_step_is_reported():
  """Same return code and attempts, failed in another sub-step: reported.

  The tst_mcp_server.py shape from the 09-13 and 09-12 comparison pairs: the
  same subprocess.TimeoutExpired, hit in one sub-step on the baseline and in
  another on the candidate, reported MODE_CHANGED (0) because the mode was
  built from the return code and attempts alone.  The kind is the literal
  string from the contract rather than a module constant, so that without
  the change the first statement to fail is the lookup of
  'failure_mode_changed' itself.
  """
  cmd = 'libtbx.python "/net/x/phenix/regression/tst_mcp_server.py"'
  baseline_stderr = [
    'Traceback (most recent call last):',
    '  File "/net/x/phenix/regression/tst_mcp_server.py", line 988, '
      'in <module>',
    '    exercise()',
    '  File "/net/x/phenix/regression/tst_mcp_server.py", line 971, '
      'in exercise',
    '    exercise_cli_help_exits_zero()',
    '  File "/net/x/phenix/regression/tst_mcp_server.py", line 941, '
      'in exercise_cli_help_exits_zero',
    '    r = subprocess.run([sys.executable, "-c", probe], env=env,',
    '        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
    '  File "/net/x/conda_base/lib/python3.11/subprocess.py", line 550, '
      'in run',
    '    stdout, stderr = process.communicate(input, timeout=timeout)',
    '  File "/net/x/conda_base/lib/python3.11/subprocess.py", line 1253, '
      'in _check_timeout',
    '    raise TimeoutExpired(',
    "subprocess.TimeoutExpired: Command '['python', '-c', 'probe']' "
      "timed out after 60 seconds",
    ]
  candidate_stderr = [
    line.replace('exercise_cli_help_exits_zero',
                 'exercise_cli_version_exits_zero')
    for line in baseline_stderr]
  assert candidate_stderr != baseline_stderr
  baseline = t96_roster.parse_roster(
    failure_log([(cmd, 'FAIL', 1, 3, baseline_stderr)]))
  candidate = t96_roster.parse_roster(
    failure_log([(cmd, 'FAIL', 1, 3, candidate_stderr)]))
  # The fixtures say what they are meant to say before anything is compared.
  for roster in (baseline, candidate):
    test = roster.tests[cmd]
    assert test.outcome == t96_roster.FAIL, test.outcome
    assert (test.return_code, test.attempt, test.attempts_total) == (1, 3, 3), \
      (test.return_code, test.attempt, test.attempts_total)
  diff = t96_roster.compare_rosters(baseline, candidate)
  # (i) The existing view is unchanged: same return code and attempts, so no
  # mode change, and a both-sides failure does not block.
  assert diff['mode_changed'] == [], diff['mode_changed']
  assert cmd not in diff['blocking'], diff['blocking']
  # (ii) The failure text shows the difference, under the call path.
  baseline_path = [[['tst_mcp_server.py', '<module>'],
                    ['tst_mcp_server.py', 'exercise'],
                    ['tst_mcp_server.py', 'exercise_cli_help_exits_zero'],
                    ['subprocess.py', 'run'],
                    ['subprocess.py', '_check_timeout']]]
  candidate_path = [[['tst_mcp_server.py', '<module>'],
                     ['tst_mcp_server.py', 'exercise'],
                     ['tst_mcp_server.py', 'exercise_cli_version_exits_zero'],
                     ['subprocess.py', 'run'],
                     ['subprocess.py', '_check_timeout']]]
  assert diff['failure_mode_changed'] == [
    (cmd, 'call path', baseline_path, candidate_path)], \
    diff['failure_mode_changed']

def exercise_failure_text_is_the_standard_error_block():
  """failure_text is the four-space block under '  Standard error:', exactly.

  Four characters come off each line and nothing else, so deeper indentation
  survives; the block ends at the first line without four leading spaces;
  trailing blank entries go and internal ones stay.  A '  Standard error:'
  that is itself indented deeper is text, not a new header, and a test that
  prints the header into its own standard output does not get a failure
  text from it.
  """
  cmd_text = cmd_for('tst_text.py')
  cmd_none = cmd_for('tst_no_stderr.py')
  cmd_spoof = cmd_for('tst_stdout_spoof.py')
  cmd_blank = cmd_for('tst_blank_stderr.py')
  cmd_pass = cmd_for('tst_pass.py')
  cmd_dup = cmd_for('tst_dup.py')
  cmd_hung = cmd_for('tst_hung.py')
  text = as_text([
    'Test %s repeated, skipping' % cmd_dup,
    'Running 6 tests on 1 processors:',
    '  ' + cmd_text,
    '  ' + cmd_none,
    '  ' + cmd_spoof,
    '  ' + cmd_blank,
    '  ' + cmd_pass,
    '  ' + cmd_hung,
    '',
    '%s [FAIL] 1.0s' % cmd_text,
    '  Time:  1.00',
    '  Return code: 1',
    '  OKs: 0',
    '  Standard out:',
    '    stdout before the header',
    '  Standard error:',
    '    first line',
    '      six-space line',
    '    ',
    '      Standard error:',
    '    last line',
    '    ',
    '    ',
    '  two-space line',
    '    four-space line after the end of the block',
    '%s [FAIL] 1.0s' % cmd_none,
    '  Time:  1.00',
    '  Return code: 1',
    '  OKs: 0',
    '%s [FAIL] 1.0s' % cmd_spoof,
    '  Time:  1.00',
    '  Return code: 1',
    '  OKs: 0',
    '  Standard out:',
    '    Standard error:',
    '    printed by the test itself',
    '%s [FAIL] 1.0s' % cmd_blank,
    '  Time:  1.00',
    '  Return code: 1',
    '  OKs: 0',
    '  Standard error:',
    '    ',
    '%s [OK] 1.0s' % cmd_pass,
    EQ80,
    '',
    'Summary:',
    '  Tests run                    : 5',
    '  Failures                     : 4',
    '  Warnings (possible failures) : 0',
    '  Known Failures (  0)         : 0',
    '  Known Unstable (  0)         : 0',
    '  Stderr output (discouraged)  : 2',
    STARS80,
    '  WARNING: NOT ALL TESTS FINISHED!',
    STARS80,
    ])
  r = t96_roster.parse_roster(text)
  assert r.tests[cmd_dup].outcome == t96_roster.SKIPPED, r.tests[cmd_dup]
  assert r.tests[cmd_hung].outcome == t96_roster.MISSING, r.tests[cmd_hung]
  got = r.tests[cmd_text].failure_text
  assert got == ['first line', '  six-space line', '', '  Standard error:',
                 'last line'], got
  # No '  Standard error:' line in the detail block: no failure text.
  for command in (cmd_none, cmd_spoof, cmd_pass, cmd_dup, cmd_hung):
    got = r.tests[command].failure_text
    assert got is None, (command, got)
  # A header over nothing but blank lines: the trailing blanks go, and an
  # empty list is left.
  got = r.tests[cmd_blank].failure_text
  assert got == [], got
  # No replay block here, so no copy can disagree.
  for command in (cmd_text, cmd_none, cmd_spoof, cmd_blank, cmd_pass):
    got = r.tests[command].failure_text_conflict
    assert got is False, (command, got)

def exercise_failure_text_conflict_with_replay_copy():
  """A disagreeing replay copy sets failure_text_conflict and changes nothing.

  The streaming result stays the record: its return code, mode and failure
  text are kept whatever the replay copy after the '=' rule says.
  """
  err = traceback_text(FM_FRAMES, 'AssertionError')
  cmd_same = cmd_for('tst_replay_same.py')
  cmd_text = cmd_for('tst_replay_text.py')
  cmd_rc = cmd_for('tst_replay_rc.py')
  cmd_pass = cmd_for('tst_replay_pass.py')
  r = t96_roster.parse_roster(failure_log(
    [(cmd_same, 'FAIL', 1, 3, err),
     (cmd_text, 'FAIL', 1, 3, err),
     (cmd_rc, 'FAIL', 1, 3, err),
     (cmd_pass, 'OK', 0, 1, None)],
    replay=[(cmd_same, 'FAIL', 1, 3, err),
            (cmd_text, 'FAIL', 1, 3, err[:-1] + ['AssertionError: other']),
            (cmd_rc, 'FAIL', 2, 3, err)]))
  assert r.tests[cmd_same].failure_text_conflict is False, \
    r.tests[cmd_same].failure_text_conflict
  assert r.tests[cmd_pass].failure_text_conflict is False, \
    r.tests[cmd_pass].failure_text_conflict
  assert r.tests[cmd_text].failure_text_conflict is True, \
    r.tests[cmd_text].failure_text_conflict
  assert r.tests[cmd_rc].failure_text_conflict is True, \
    r.tests[cmd_rc].failure_text_conflict
  # Replay lines never create or modify anything.
  assert r.tests[cmd_same].failure_text == err, r.tests[cmd_same].failure_text
  assert r.tests[cmd_text].failure_text == err, r.tests[cmd_text].failure_text
  assert r.tests[cmd_rc].return_code == 1, r.tests[cmd_rc].return_code
  assert r.tests[cmd_rc].mode == 'rc=1; attempt 3 of 3', r.tests[cmd_rc].mode

def exercise_failure_mode_kinds():
  """Each layer, differing alone, is reported as its kind with its two values.

  Every pair shares outcome, return code and attempts unless the kind under
  test is one of those, so no earlier check can decide it.
  """
  diff = t96_roster.failure_mode_difference
  changed = t96_roster.CHANGED
  plain = traceback_text(FM_FRAMES, 'AssertionError: expected 1')
  chained = (plain + FM_CHAINED_SEPARATOR
             + traceback_text(FM_CLEANUP_FRAMES, FM_CLEANUP_ERROR))
  # Outcome: the details are the two outcome constants.
  got = diff(failure_test(plain), failure_test(plain, status='EXPECTED FAIL'))
  assert got == (changed, t96_roster.KIND_OUTCOME, t96_roster.FAIL,
                 t96_roster.EXPECTED_FAIL), got
  # Return code, then attempts: [return_code, attempt, attempts_total].
  got = diff(failure_test(plain, attempts=3),
             failure_test(plain, return_code=2, attempts=3))
  assert got == (changed, t96_roster.KIND_RETURN_CODE,
                 [1, 3, 3], [2, 3, 3]), got
  got = diff(failure_test(plain, attempts=3), failure_test(plain))
  assert got == (changed, t96_roster.KIND_RETURN_CODE,
                 [1, 3, 3], [1, 1, 1]), got
  # Exception types, one per traceback: the candidate adds a chained
  # exception (the tst_chat_window_mcp.py shape) ...
  got = diff(failure_test(plain), failure_test(chained))
  assert got == (changed, t96_roster.KIND_EXCEPTION_TYPES,
                 ['AssertionError'],
                 ['AssertionError', 'FileNotFoundError']), got
  # ... a dotted type is one name ...
  got = diff(
    failure_test(traceback_text(FM_FRAMES, 'subprocess.TimeoutExpired: x')),
    failure_test(traceback_text(FM_FRAMES, 'subprocess.CalledProcessError: x')))
  assert got == (changed, t96_roster.KIND_EXCEPTION_TYPES,
                 ['subprocess.TimeoutExpired'],
                 ['subprocess.CalledProcessError']), got
  # ... and an exception line the pattern does not match has type None.
  got = diff(
    failure_test(traceback_text(FM_FRAMES, 'KeyError (no colon)')),
    failure_test(traceback_text(FM_FRAMES, 'KeyError: (no colon)')))
  assert got == (changed, t96_roster.KIND_EXCEPTION_TYPES,
                 [None], ['KeyError']), got
  # Call path, per traceback: [file name, function] of every frame.  Here
  # the second traceback was raised from another function.
  other_cleanup = [FM_CLEANUP_FRAMES[0],
                   ('/net/x/lib/python3.11/shutil.py', 701, '_rmtree_unsafe')]
  got = diff(failure_test(chained), failure_test(
    plain + FM_CHAINED_SEPARATOR
    + traceback_text(other_cleanup, FM_CLEANUP_ERROR)))
  assert got == (changed, t96_roster.KIND_CALL_PATH,
                 [FM_CALL_PATH,
                  [['tst_mode.py', 'run'], ['shutil.py', '_rmtree_safe_fd']]],
                 [FM_CALL_PATH,
                  [['tst_mode.py', 'run'], ['shutil.py', '_rmtree_unsafe']]]), \
    got
  # The file name is what follows the last '/' or '\' ...
  got = diff(
    failure_test(traceback_text([('/net/base/tst_a.py', 5, 'run')],
                                'AssertionError')),
    failure_test(traceback_text([('C:\\cand\\tst_b.py', 5, 'run')],
                                'AssertionError')))
  assert got == (changed, t96_roster.KIND_CALL_PATH,
                 [[['tst_a.py', 'run']]], [[['tst_b.py', 'run']]]), got
  # ... so the same file under another directory differs in no named layer,
  # only in the remaining failure text, reported as the two whole texts.
  in_base = traceback_text([('/net/base/tst_a.py', 5, 'run')],
                           'AssertionError')
  in_cand = traceback_text([('C:\\cand\\tst_a.py', 5, 'run')],
                           'AssertionError')
  got = diff(failure_test(in_base), failure_test(in_cand))
  assert got == (changed, t96_roster.KIND_REMAINING_TEXT,
                 in_base, in_cand), got
  # Exception message: the exception lines, verbatim.
  got = diff(failure_test(plain), failure_test(
    traceback_text(FM_FRAMES, 'AssertionError: expected 2')))
  assert got == (changed, t96_roster.KIND_EXCEPTION_MESSAGE,
                 ['AssertionError: expected 1'],
                 ['AssertionError: expected 2']), got
  # Line numbers, per traceback, as int.
  moved = [FM_FRAMES[0], ('/net/x/regression/tst_mode.py', 21, 'run')]
  got = diff(failure_test(plain), failure_test(
    traceback_text(moved, 'AssertionError: expected 1')))
  assert got == (changed, t96_roster.KIND_LINE_NUMBERS,
                 [[30, 20]], [[30, 21]]), got
  # Other failure text: every entry outside the traceback headers, frame
  # runs and exception lines, in order -- here a line before the first
  # traceback, the chained-exception separator, and a line after the last.
  tail = ['QThread: Destroyed while thread is still running']
  got = diff(
    failure_test(['Sorry: Bond length too long : 16.62'] + chained + tail),
    failure_test(['Sorry: Bond length too long : 17.00'] + chained + tail))
  assert got == (changed, t96_roster.KIND_OTHER_TEXT,
                 ['Sorry: Bond length too long : 16.62']
                 + FM_CHAINED_SEPARATOR + tail,
                 ['Sorry: Bond length too long : 17.00']
                 + FM_CHAINED_SEPARATOR + tail), got
  # Two texts with no traceback at all, differing in an indented line.
  got = diff(failure_test(['DeprecationWarning: invalid escape', "  key = 1"]),
             failure_test(['DeprecationWarning: invalid escape', "  key = 2"]))
  assert got == (changed, t96_roster.KIND_OTHER_TEXT,
                 ['DeprecationWarning: invalid escape', "  key = 1"],
                 ['DeprecationWarning: invalid escape', "  key = 2"]), got
  # Remaining failure text: every named layer is equal and the texts are
  # not.  The failing assert was edited in place, so frames, line numbers
  # and the bare AssertionError are unchanged and only the source line under
  # the innermost frame differs; both whole texts are reported.
  assert_x = (traceback_text(FM_FRAMES, None)[:-2]
              + ['    assert x', '    ^^^^^^^^', 'AssertionError'])
  assert_y = (traceback_text(FM_FRAMES, None)[:-2]
              + ['    assert y', '    ^^^^^^^^', 'AssertionError'])
  got = diff(failure_test(assert_x), failure_test(assert_y))
  assert got == (changed, t96_roster.KIND_REMAINING_TEXT,
                 assert_x, assert_y), got

def exercise_failure_mode_unknown_reasons():
  """What cannot be compared is UNKNOWN with its reason, never the same mode."""
  diff = t96_roster.failure_mode_difference
  unknown = t96_roster.UNKNOWN
  text = traceback_text(FM_FRAMES, 'AssertionError')
  # The streaming and replay copies disagree: in the text on the baseline,
  # in the return code on the candidate.
  got = diff(failure_test(text, replay=(1, text + ['extra line'])),
             failure_test(text))
  assert got == (unknown, t96_roster.REASON_REPLAY_CONFLICT), got
  got = diff(failure_test(text), failure_test(text, replay=(2, text)))
  assert got == (unknown, t96_roster.REASON_REPLAY_CONFLICT), got
  # No failure text: None on both sides, and an empty list counts as none.
  got = diff(failure_test(None), failure_test(None))
  assert got == (unknown, t96_roster.REASON_NO_TEXT_EITHER), got
  got = diff(failure_test(['']), failure_test(None))
  assert got == (unknown, t96_roster.REASON_NO_TEXT_EITHER), got
  got = diff(failure_test(None), failure_test(text))
  assert got == (unknown, t96_roster.REASON_NO_TEXT_BASELINE), got
  got = diff(failure_test(text), failure_test(None))
  assert got == (unknown, t96_roster.REASON_NO_TEXT_CANDIDATE), got
  # A traceback without an exception line: the text ends after the frame
  # run on the baseline; a blank entry follows it on the candidate; and
  # identical truncated texts are still not the same mode.
  cut = traceback_text(FM_FRAMES, None)
  got = diff(failure_test(cut), failure_test(text))
  assert got == (unknown, t96_roster.REASON_TRUNCATED_TRACEBACK), got
  got = diff(failure_test(text),
             failure_test(cut + ['', 'AssertionError']))
  assert got == (unknown, t96_roster.REASON_TRUNCATED_TRACEBACK), got
  got = diff(failure_test(cut), failure_test(cut))
  assert got == (unknown, t96_roster.REASON_TRUNCATED_TRACEBACK), got

def exercise_failure_mode_check_order():
  """When two checks would both apply, the earlier check decides.

  A result of None fails as an assertion showing None, not as a TypeError.
  """
  diff = t96_roster.failure_mode_difference
  changed = t96_roster.CHANGED
  unknown = t96_roster.UNKNOWN
  at_5 = [('/net/x/regression/tst_a.py', 5, 'run')]
  at_6 = [('/net/x/regression/tst_a.py', 6, 'run')]
  in_b = [('/net/x/regression/tst_b.py', 5, 'run')]
  err = traceback_text(at_5, 'AssertionError: one')
  cut = traceback_text(at_5, None)
  # Outcome before return code.
  got = diff(failure_test(err),
             failure_test(err, status='EXPECTED FAIL', return_code=2))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_OUTCOME), got
  # Return code before a replay conflict.
  got = diff(failure_test(err),
             failure_test(err, return_code=2, replay=(2, err + ['extra'])))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_RETURN_CODE), got
  # A replay conflict before no text on either side.
  got = diff(failure_test(None, replay=(1, err)), failure_test(None))
  assert got == (unknown, t96_roster.REASON_REPLAY_CONFLICT), got
  # No text before a truncated traceback.
  got = diff(failure_test(None), failure_test(cut))
  assert got == (unknown, t96_roster.REASON_NO_TEXT_BASELINE), got
  # A truncated traceback before any layer.
  got = diff(failure_test(cut),
             failure_test(traceback_text(in_b, 'KeyError: two')))
  assert got == (unknown, t96_roster.REASON_TRUNCATED_TRACEBACK), got
  # Exception types before call path.
  got = diff(failure_test(err),
             failure_test(traceback_text(in_b, 'KeyError: one')))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_EXCEPTION_TYPES), got
  # Call path before exception message.
  got = diff(failure_test(err),
             failure_test(traceback_text(in_b, 'AssertionError: two')))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_CALL_PATH), got
  # Exception message before line numbers.
  got = diff(failure_test(err),
             failure_test(traceback_text(at_6, 'AssertionError: two')))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_EXCEPTION_MESSAGE), got
  # Line numbers before other failure text.
  got = diff(failure_test(['note one'] + err),
             failure_test(['note two']
                          + traceback_text(at_6, 'AssertionError: one')))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_LINE_NUMBERS), got
  # Other failure text before remaining failure text.  The edited source
  # line alone is remaining text, so both checks apply once the other text
  # differs too.
  edited = err[:2] + ['    assert y'] + err[3:]
  got = diff(failure_test(['note one'] + err),
             failure_test(['note one'] + edited))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_REMAINING_TEXT), got
  got = diff(failure_test(['note one'] + err),
             failure_test(['note two'] + edited))
  assert got is not None and got[:2] == (
    changed, t96_roster.KIND_OTHER_TEXT), got

def exercise_identical_failures_give_no_entry():
  """The same failure text on both sides is the same mode: no entry anywhere."""
  cmd_chained = cmd_for('tst_same_chained.py')
  cmd_plain = cmd_for('tst_same_plain.py')
  chained = (['This plugin does not support propagateSizeHints()']
             + traceback_text(FM_FRAMES, 'AssertionError: expected 1')
             + FM_CHAINED_SEPARATOR
             + traceback_text(FM_CLEANUP_FRAMES, FM_CLEANUP_ERROR)
             + ['QThread: Destroyed while thread is still running'])
  plain = ['phenix.mcp.registry: skipping taam_minus_iam.py: '
           'PyDiscamb not installed.']
  entries = [(cmd_chained, 'FAIL', -6, 3, chained),
             (cmd_plain, 'FAIL', 1, 3, plain)]
  baseline = t96_roster.parse_roster(failure_log(entries))
  candidate = t96_roster.parse_roster(failure_log(entries))
  for command in (cmd_chained, cmd_plain):
    got = t96_roster.failure_mode_difference(baseline.tests[command],
                                             candidate.tests[command])
    assert got is None, (command, got)
  diff = t96_roster.compare_rosters(baseline, candidate)
  assert diff['failure_mode_changed'] == [], diff['failure_mode_changed']
  assert diff['failure_mode_unknown'] == [], diff['failure_mode_unknown']
  assert diff['needs_classification'] == [], diff['needs_classification']

def exercise_failure_mode_population_and_verdict():
  """Only commands ranked 1 on both sides are compared; the verdict holds.

  EXPECTED FAIL, EXPECTED UNSTABLE and WARNING pairs are compared as FAIL
  pairs are.  A PASS pair is not, whatever its standard error, and neither
  is a failure against SKIPPED, MISSING or absence.  The new lists enter
  neither blocking nor mode_changed, which keep their existing values.
  Commands are named so that sorted order is the reverse of log order.
  """
  one = traceback_text(FM_FRAMES, 'AssertionError: one')
  two = traceback_text(FM_FRAMES, 'AssertionError: two')
  cmd_w = cmd_for('tst_a_warning_both.py')
  cmd_ef = cmd_for('tst_b_expected_fail_both.py')
  cmd_eu = cmd_for('tst_c_expected_unstable_both.py')
  cmd_f2e = cmd_for('tst_d_fail_to_expected_fail.py')
  cmd_rc = cmd_for('tst_e_fail_return_code_changed.py')
  cmd_pass = cmd_for('tst_f_pass_both.py')
  cmd_f2s = cmd_for('tst_g_fail_to_skipped.py')
  cmd_m2f = cmd_for('tst_h_missing_to_fail.py')
  cmd_f2x = cmd_for('tst_i_fail_to_absent.py')
  baseline = t96_roster.parse_roster(failure_log([
    (cmd_f2x, 'FAIL', 1, 1, one),
    (cmd_f2s, 'FAIL', 1, 1, one),
    (cmd_pass, 'OK', 0, 1, ['note one']),
    (cmd_rc, 'FAIL', 1, 1, one),
    (cmd_f2e, 'FAIL', 1, 1, one),
    (cmd_eu, 'EXPECTED UNSTABLE', 1, 1, None),
    (cmd_ef, 'EXPECTED FAIL', 1, 1, one),
    (cmd_w, 'WARNING', 0, 1, ['warning one']),
    ], missing=[cmd_m2f]))
  candidate = t96_roster.parse_roster(failure_log([
    (cmd_m2f, 'FAIL', 1, 1, one),
    (cmd_pass, 'OK', 0, 1, ['note two']),
    (cmd_rc, 'FAIL', 2, 1, one),
    (cmd_f2e, 'EXPECTED FAIL', 1, 1, one),
    (cmd_eu, 'EXPECTED UNSTABLE', 1, 1, None),
    (cmd_ef, 'EXPECTED FAIL', 1, 1, two),
    (cmd_w, 'WARNING', 0, 1, ['warning two']),
    ], skipped=[cmd_f2s]))
  # The fixtures say what they are meant to say before anything is compared.
  assert baseline.tests[cmd_m2f].outcome == t96_roster.MISSING
  assert candidate.tests[cmd_f2s].outcome == t96_roster.SKIPPED
  assert cmd_f2x not in candidate.tests, sorted(candidate.tests)
  assert baseline.tests[cmd_pass].outcome == t96_roster.PASS
  assert candidate.tests[cmd_w].outcome == t96_roster.WARNING
  diff = t96_roster.compare_rosters(baseline, candidate)
  assert sorted(diff) == sorted([
    'new_failures', 'fixed', 'mode_changed', 'disappeared', 'newly_skipped',
    'newly_timed_out', 'newly_missing', 'newly_not_passing', 'new_tests',
    'blocking', 'failure_mode_changed', 'failure_mode_unknown',
    'needs_classification']), sorted(diff)
  assert diff['failure_mode_changed'] == [
    (cmd_w, t96_roster.KIND_OTHER_TEXT, ['warning one'], ['warning two']),
    (cmd_ef, t96_roster.KIND_EXCEPTION_MESSAGE,
     ['AssertionError: one'], ['AssertionError: two']),
    (cmd_f2e, t96_roster.KIND_OUTCOME,
     t96_roster.FAIL, t96_roster.EXPECTED_FAIL),
    (cmd_rc, t96_roster.KIND_RETURN_CODE, [1, 1, 1], [2, 1, 1]),
    ], diff['failure_mode_changed']
  assert diff['failure_mode_unknown'] == [
    (cmd_eu, t96_roster.REASON_NO_TEXT_EITHER)], diff['failure_mode_unknown']
  assert diff['needs_classification'] == [
    cmd_w, cmd_ef, cmd_eu, cmd_f2e, cmd_rc], diff['needs_classification']
  # The existing keys keep their values: mode_changed is still FAIL pairs
  # with a different return code or attempts, and blocking is still the rank
  # rule -- the three commands that lost a run or newly fail, and none of
  # the five that need classification.
  assert diff['mode_changed'] == [(cmd_rc, 'rc=1', 'rc=2')], \
    diff['mode_changed']
  assert diff['blocking'] == [cmd_f2s, cmd_m2f, cmd_f2x], diff['blocking']

def exercise_failure_mode_constants():
  """The failure-mode constants, their tuples and the two new origin labels."""
  assert t96_roster.CHANGED == 'changed', t96_roster.CHANGED
  assert t96_roster.UNKNOWN == 'unknown', t96_roster.UNKNOWN
  kinds = [
    ('KIND_OUTCOME', 'outcome'),
    ('KIND_RETURN_CODE', 'return code or attempts'),
    ('KIND_EXCEPTION_TYPES', 'exception types'),
    ('KIND_CALL_PATH', 'call path'),
    ('KIND_EXCEPTION_MESSAGE', 'exception message'),
    ('KIND_LINE_NUMBERS', 'line numbers'),
    ('KIND_OTHER_TEXT', 'other failure text'),
    ('KIND_REMAINING_TEXT', 'remaining failure text'),
    ]
  reasons = [
    ('REASON_REPLAY_CONFLICT',
     'streaming and replay copies of the failure text disagree'),
    ('REASON_NO_TEXT_EITHER', 'no failure text on either side'),
    ('REASON_NO_TEXT_BASELINE', 'no failure text on baseline'),
    ('REASON_NO_TEXT_CANDIDATE', 'no failure text on candidate'),
    ('REASON_TRUNCATED_TRACEBACK', 'traceback without an exception line'),
    ]
  for (name, value) in kinds + reasons:
    got = getattr(t96_roster, name)
    assert got == value, (name, got)
  got = t96_roster.FAILURE_MODE_KINDS
  assert got == tuple(value for (name, value) in kinds), got
  got = t96_roster.UNKNOWN_REASONS
  assert got == tuple(value for (name, value) in reasons), got
  assert callable(t96_roster.failure_mode_difference)
  for key in ('tests[].failure_text', 'tests[].failure_text_conflict'):
    got = t96_roster.FIELD_ORIGIN[key]
    assert got == t96_roster.ORIGIN_CANDIDATE_CLAIM, (key, got)

def run():
  """Run every exercise in turn; print OK if all of them pass."""
  exercise_constants()
  exercise_replay_block_not_double_counted()
  exercise_indented_status_bracket_is_not_a_result()
  exercise_outcome_word_in_file_name()
  exercise_non_python_command()
  exercise_command_containing_bracketed_text()
  exercise_preamble_is_ignored()
  exercise_parallel_suffix_is_not_part_of_command()
  exercise_missing_test()
  exercise_duplicate_command_skip()
  exercise_skip_line_does_not_override_a_result()
  exercise_reconcile_reports_failure_mismatch()
  exercise_reconcile_declared_count_vs_tests_run()
  exercise_absent_summary_line_is_not_zero()
  exercise_strip_profile_lines()
  exercise_retry_annotations_and_mode()
  exercise_summary_parenthetical()
  exercise_compare_rosters_transitions()
  exercise_newly_missing_blocks()
  exercise_blocking_is_an_outcome_rank_rule()
  exercise_field_origin_covers_every_reported_field()
  exercise_failure_mode_change_in_sub_step_is_reported()
  exercise_failure_text_is_the_standard_error_block()
  exercise_failure_text_conflict_with_replay_copy()
  exercise_failure_mode_kinds()
  exercise_failure_mode_unknown_reasons()
  exercise_failure_mode_check_order()
  exercise_identical_failures_give_no_entry()
  exercise_failure_mode_population_and_verdict()
  exercise_failure_mode_constants()
  print("OK")

if __name__ == '__main__':
  run()
