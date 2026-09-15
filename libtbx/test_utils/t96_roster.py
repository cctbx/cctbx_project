# Roster parser for libtbx.run_tests_parallel / phenix_regression.test_all_parallel
# output.
#
# The suite exits 0 even when it reports failures, so a process exit status is
# never a verdict here.  Outcomes come from log CONTENT only: this module turns
# suite output text into a roster of named tests with outcomes, and compares two
# such rosters to find the transitions (new failures, disappearances, ...) that
# block a change.
#
# Everything below is stdlib-only and has no file, process or network access, so
# the module can be dropped unchanged beside libtbx/test_utils/parallel.py.
#
# Trust boundary.  The log is partly candidate-controlled output: a test can
# print "[OK]" at column 0 of its own stdout.  Every value reported here
# therefore has one of three origins, and a consumer must be able to tell which
# of its numbers are the harness speaking and which the code under test
# influenced:
#
#   1. HARNESS STATEMENT -- the suite's own "Summary:" block, and with it the
#      "Running <N> tests on <P> processors:" header, which the harness prints
#      from its own list.  The harness counts these itself; the candidate
#      cannot influence the number, only the events behind it.  The COMMANDS
#      listed under that header are a different matter: they are text out of
#      the tree under test, so they carry a compound label in FIELD_ORIGIN.
#   2. HARNESS-FORMATTED, CANDIDATE-TRIGGERED -- a streaming result line.  The
#      harness chooses the format and the status token; the candidate's exit
#      code and stderr are what select which token appears.  The COMMAND STRING
#      on that line is candidate-influenced text.
#   3. CANDIDATE CLAIM, DISPLAYED AND NEVER COUNTED -- indented captured test
#      output.  This module never derives any number from it: no outcome, no
#      count, nothing in compare_rosters()' blocking and nothing reconcile()
#      reports.  It is READ for exactly one purpose: deciding whether a test
#      that ran and did not pass on BOTH sides (outcome rank 1 in both
#      rosters) failed the same way.  The standard error under a result is
#      kept as TestOutcome.failure_text, and failure_mode_difference()
#      compares two of them; what that finds goes to compare_rosters()'
#      needs_classification, never to blocking.
#
# ORIGIN_* and FIELD_ORIGIN below state this per reported field, so a consumer
# can label a number programmatically rather than by reading this prose.
#
# Two facts make the boundary hold in this parser: a result line is recognised
# at column 0 ONLY (_RESULT_RE begins with \S, and _to_lines splits on '\n'
# alone, so a form feed inside captured output cannot manufacture a column-0
# line), and parallel.py indents captured test output by four spaces.  Candidate
# text therefore cannot present itself as a result line unless the harness's own
# framing is broken.
#
# This is a statement of PROVENANCE, not a security control: the log is a
# trusted-ish artefact produced by the harness, and a candidate that can write
# arbitrary bytes to the harness's own stdout stream at column 0 is outside
# what any parser can distinguish.  That is a real mode of the harness rather
# than a hypothesis: at verbosity EXTRA_VERBOSE (2, "for nightly builds")
# parallel.py's determine_result_status() prints result.error_lines raw and
# unindented.  At the default verbosity (1), captured output reaches
# the log only four-space indented.
#
# One reported field is not column-0 anchored.  not_all_finished is matched on
# stripped content, because parallel.py prints that warning indented, so
# captured output containing exactly that line would raise it too; it can only
# ever be raised, never cleared, so the error runs toward distrusting a run.
# Its FIELD_ORIGIN entry says so rather than calling it a clean harness
# statement.
#
# Limits of the failure text.  parallel.py prints a result's "Standard error:"
# block to sys.stderr, not to the stream its result lines go to, so the block
# is in a log only when the suite's standard error was captured into the same
# file; from a log without it, every comparison that gets past the outcome and
# return-code checks is UNKNOWN.  And parallel.py's run_command() keeps only
# the FINAL attempt's output, so a retried test is compared on its last
# attempt alone.
#
# No origin in this log reports a timeout, which is why a hung test surfaces as
# MISSING; see the TIMED_OUT note below.
#
# Note on TIMED_OUT: parallel.py has no timeout mechanism, so this parser never
# produces that outcome.  A hung test simply never prints a result line and so
# appears as MISSING.  The constant exists only so roster comparison has the
# right shape for a future suite that does report timeouts; no heuristic ever
# infers a timeout from free text (test file names such as tst_reduce_timeout.py
# pass routinely).

from __future__ import absolute_import, division, print_function

import re

# Outcome constants.  These are the values stored in TestOutcome.outcome.
PASS = 'pass'
FAIL = 'fail'
WARNING = 'warning'
EXPECTED_FAIL = 'expected_fail'
EXPECTED_UNSTABLE = 'expected_unstable'
SKIPPED = 'skipped'
MISSING = 'missing'
TIMED_OUT = 'timed_out'

# Every outcome this module knows about, in a stable order.
ALL_OUTCOMES = (PASS, FAIL, WARNING, EXPECTED_FAIL, EXPECTED_UNSTABLE,
                SKIPPED, MISSING, TIMED_OUT)

# Outcome rank: what one run of one test is worth.  compare_rosters() decides
# what blocks from this rank plus one further clause for a test that newly
# FAILS (see compare_rosters), so every outcome pair is covered and no
# transition can slip through for want of a named category.
#
#   rank 2  PASS                                             ran and passed
#   rank 1  FAIL, WARNING, EXPECTED_FAIL, EXPECTED_UNSTABLE  ran, did not pass
#   rank 0  SKIPPED, MISSING, TIMED_OUT                      did not run
#
# A command ABSENT from a roster ranks 0 as well (ABSENT_RANK): a run that
# never happened is worth no more than one that did not run.  An outcome that
# is not in this table ranks 0 too, so an unknown outcome fails safe -- it
# blocks against any baseline that did better instead of passing silently.
ABSENT_RANK = 0

OUTCOME_RANK = {
  PASS: 2,
  FAIL: 1,
  WARNING: 1,
  EXPECTED_FAIL: 1,
  EXPECTED_UNSTABLE: 1,
  SKIPPED: 0,
  MISSING: 0,
  TIMED_OUT: 0,
}

# Origin labels: which kind of log text a reported value came from.  See the
# Trust boundary section in the header comment for what each one means.
ORIGIN_HARNESS_STATEMENT = 'harness statement'
ORIGIN_HARNESS_FORMATTED = 'harness-formatted, candidate-triggered'
ORIGIN_CANDIDATE_CLAIM = 'candidate claim, displayed never counted'

# Two compound labels, for the fields none of the three describes exactly.  They
# are spelled out rather than rounded to the nearest neighbour, so a consumer
# printing a label is never told more than is true.
ORIGIN_HARNESS_STATEMENT_CANDIDATE_TEXT = (
  'harness statement; the command strings in it are candidate-influenced text')
ORIGIN_HARNESS_STATEMENT_UNANCHORED = (
  'harness statement, matched unanchored: captured output could raise it')

# Every origin label, in a stable order.
ALL_ORIGINS = (ORIGIN_HARNESS_STATEMENT,
               ORIGIN_HARNESS_FORMATTED,
               ORIGIN_CANDIDATE_CLAIM,
               ORIGIN_HARNESS_STATEMENT_CANDIDATE_TEXT,
               ORIGIN_HARNESS_STATEMENT_UNANCHORED)

# The origin of every value this module reports.  Keys are Roster attribute
# names, plus 'summary.<key>' for one entry of Roster.summary, 'counts()[<o>]'
# for one value of Roster.counts(), and 'tests[].<attr>' for one attribute of a
# TestOutcome held in Roster.tests.
#
# ORIGIN_CANDIDATE_CLAIM labels exactly two fields, tests[].failure_text and
# tests[].failure_text_conflict, and neither is ever counted: no count and no
# outcome in this table, nothing in compare_rosters()' blocking and nothing
# reconcile() reports is derived from captured test output.  Captured output
# is READ for one purpose only, deciding whether a failure's mode changed
# (failure_mode_difference()).  The one other field captured output could
# touch is not_all_finished, which carries the _UNANCHORED label below.
FIELD_ORIGIN = {
  # Header block: "Running <N> tests on <P> processors:" and the indented
  # command list under it.  The count and the processor string are the
  # harness's own; the commands are text the harness read out of the tree under
  # test, which is why they carry the compound label instead of the plain one.
  'declared_count': ORIGIN_HARNESS_STATEMENT,
  'declared_commands': ORIGIN_HARNESS_STATEMENT_CANDIDATE_TEXT,
  'nprocs': ORIGIN_HARNESS_STATEMENT,

  # Streaming section and the replay block: one harness-formatted line per
  # result, whose status token the candidate's exit code and stderr select.
  'tests': ORIGIN_HARNESS_FORMATTED,
  'replay_commands': ORIGIN_HARNESS_FORMATTED,
  'result_line_count': ORIGIN_HARNESS_FORMATTED,

  # Summary block, which the harness counts itself.
  'summary': ORIGIN_HARNESS_STATEMENT,
  'summary.tests_run': ORIGIN_HARNESS_STATEMENT,
  'summary.failures': ORIGIN_HARNESS_STATEMENT,
  'summary.warnings': ORIGIN_HARNESS_STATEMENT,
  'summary.known_failures': ORIGIN_HARNESS_STATEMENT,
  'summary.known_unstable': ORIGIN_HARNESS_STATEMENT,
  'summary.stderr_output': ORIGIN_HARNESS_STATEMENT,
  'summary.retries_used': ORIGIN_HARNESS_STATEMENT,

  # The suite's own "NOT ALL TESTS FINISHED!" warning, which parallel.py prints
  # INDENTED, so this module matches it on stripped content rather than at
  # column 0; captured test output containing exactly that line would raise it
  # too.  See the Trust boundary section in the header comment.
  'not_all_finished': ORIGIN_HARNESS_STATEMENT_UNANCHORED,

  # Roster.counts(): every value is derived from the streaming section.
  # SKIPPED comes from the harness's column-0 "Test <cmd> repeated, skipping"
  # line, and MISSING is derived here from a declared command that never got a
  # result line; neither is read out of captured test output.  TIMED_OUT is
  # always 0, since no origin in this log reports a timeout.
  'counts()[%s]' % PASS: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % FAIL: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % WARNING: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % EXPECTED_FAIL: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % EXPECTED_UNSTABLE: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % SKIPPED: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % MISSING: ORIGIN_HARNESS_FORMATTED,
  'counts()[%s]' % TIMED_OUT: ORIGIN_HARNESS_FORMATTED,

  # One TestOutcome in Roster.tests.  Every attribute but the last two comes
  # from that test's result line or from the two-space-indented detail block
  # under it (the "Return code:" line), never from the four-space-indented
  # captured output; 'mode' is built here from those values.
  'tests[].command': ORIGIN_HARNESS_FORMATTED,
  'tests[].outcome': ORIGIN_HARNESS_FORMATTED,
  'tests[].wall_time': ORIGIN_HARNESS_FORMATTED,
  'tests[].attempt': ORIGIN_HARNESS_FORMATTED,
  'tests[].attempts_total': ORIGIN_HARNESS_FORMATTED,
  'tests[].return_code': ORIGIN_HARNESS_FORMATTED,
  'tests[].mode': ORIGIN_HARNESS_FORMATTED,
  # The captured standard error under the detail block's "Standard error:"
  # line, and whether the replay block's copy of it (or of the return code)
  # disagrees.  Candidate text, read only by failure_mode_difference().
  'tests[].failure_text': ORIGIN_CANDIDATE_CLAIM,
  'tests[].failure_text_conflict': ORIGIN_CANDIDATE_CLAIM,
}

# Failure-mode comparison (failure_mode_difference()).  Its result is None
# when two failures have the same mode, (CHANGED, kind, baseline_detail,
# candidate_detail) when they differ, and (UNKNOWN, reason) when the log does
# not hold enough to tell; UNKNOWN is never folded into "same".  A kind names
# the first layer that differs, most significant first.  The last kind,
# KIND_REMAINING_TEXT, compares the two failure texts whole (KIND_OTHER_TEXT
# holds only the entries outside every traceback), so two failure texts that
# differ in ANY entry never compare as the same mode: the named layers before
# it only choose the label, and cannot hide a difference.
CHANGED = 'changed'
UNKNOWN = 'unknown'

KIND_OUTCOME = 'outcome'
KIND_RETURN_CODE = 'return code or attempts'
KIND_EXCEPTION_TYPES = 'exception types'
KIND_CALL_PATH = 'call path'
KIND_EXCEPTION_MESSAGE = 'exception message'
KIND_LINE_NUMBERS = 'line numbers'
KIND_OTHER_TEXT = 'other failure text'
KIND_REMAINING_TEXT = 'remaining failure text'

# Every kind, in the order the layers are checked.
FAILURE_MODE_KINDS = (KIND_OUTCOME,
                      KIND_RETURN_CODE,
                      KIND_EXCEPTION_TYPES,
                      KIND_CALL_PATH,
                      KIND_EXCEPTION_MESSAGE,
                      KIND_LINE_NUMBERS,
                      KIND_OTHER_TEXT,
                      KIND_REMAINING_TEXT)

REASON_REPLAY_CONFLICT = (
  'streaming and replay copies of the failure text disagree')
REASON_NO_TEXT_EITHER = 'no failure text on either side'
REASON_NO_TEXT_BASELINE = 'no failure text on baseline'
REASON_NO_TEXT_CANDIDATE = 'no failure text on candidate'
REASON_TRUNCATED_TRACEBACK = 'traceback without an exception line'

# Every reason for UNKNOWN, in the order the checks are made.
UNKNOWN_REASONS = (REASON_REPLAY_CONFLICT,
                   REASON_NO_TEXT_EITHER,
                   REASON_NO_TEXT_BASELINE,
                   REASON_NO_TEXT_CANDIDATE,
                   REASON_TRUNCATED_TRACEBACK)

# Status token printed by parallel.py -> outcome constant.
_STATUS_TO_OUTCOME = {
  'OK': PASS,
  'WARNING': WARNING,
  'FAIL': FAIL,
  'EXPECTED FAIL': EXPECTED_FAIL,
  'EXPECTED UNSTABLE': EXPECTED_UNSTABLE,
}

# "Running %d tests on %s processors:" at column 0.
_HEADER_RE = re.compile(r'^Running (\d+) tests on (.+) processors:$')

# One declared-roster entry, printed as "  %s" % cmd: exactly two spaces then
# the command.  Anything else (blank line, deeper indent) ends the block.
_DECLARED_RE = re.compile(r'^  (\S.*)$')

# One streaming result line, printed at column 0 as
#   "%s [%s] %.1fs%s" % (command, status, wall_time, retry_note)
# The command is NON-GREEDY and the remainder is END-ANCHORED, so the split is
# fixed by the status/time/note tail rather than by a scan for '['.  A command
# that itself contains something like " [OK] 1.0s" therefore still parses
# correctly: the short split fails to reach the end of the line and the engine
# backtracks to the real one.  The leading \S keeps indented lines (captured
# test output such as "    [WARNING] RECOVERY LIMIT REACHED") from ever
# matching.
_RESULT_RE = re.compile(
  r'^(?P<command>\S.*?)'
  r' \[(?P<status>OK|WARNING|FAIL|EXPECTED FAIL|EXPECTED UNSTABLE)\]'
  r' (?P<wall_time>-?\d+\.\d+)s'
  r'(?:'
  r'  \(passed on attempt (?P<attempt>\d+) of (?P<attempts_total>\d+)\)'
  r'|'
  r'  \(failed after (?P<failed_attempts>\d+) attempts\)'
  r')?$')

# "Test %s repeated, skipping" at column 0, printed for a duplicate command.
_SKIP_RE = re.compile(r'^Test (.+) repeated, skipping$')

# "  Return code: %s" in a result's detail block.  The detail block is indented
# by exactly two spaces; captured test output is indented by four, so a
# "Return code:" line inside a test's own output cannot match.
_RETURN_CODE_RE = re.compile(r'^  Return code: (.*)$')

# "  Standard error:" in a result's detail block, matched as the WHOLE line.
# parallel.py prints each captured line under it as "    " + line, so the
# failure text is the following lines that begin with exactly that prefix.
_STDERR_HEADER = '  Standard error:'
_CAPTURED_INDENT = '    '

# Structure of a failure text, used only to compare failure modes: a traceback
# header, the frame lines of its run, and the exception type at the start of
# its exception line.
_TRACEBACK_HEADER = 'Traceback (most recent call last):'
_FRAME_RE = re.compile(
  r'^  File "(?P<path>[^"]*)", line (?P<line>\d+), in (?P<func>.*)$')
_EXCEPTION_TYPE_RE = re.compile(r'^([A-Za-z_][\w.]*)(?::|$)')

# Column-0 line of eighty (or more) '=' characters: end of the streaming
# section.
_TERMINATOR_RE = re.compile(r'^={80,}$')

# Column-0 line introducing the replay block, in which every warning and every
# failure is printed a SECOND time in the identical format.  Used only as a
# fallback end-of-streaming marker for truncated logs that lack the '=' rule.
_REPLAY_HEADER = ('Error: the following jobs returned non-zero exit codes or '
                  'suspicious stderr output:')

_SUMMARY_HEADER = 'Summary:'

_NOT_ALL_FINISHED = 'WARNING: NOT ALL TESTS FINISHED!'

# Indented "label : value" line inside the Summary block.
_SUMMARY_ITEM_RE = re.compile(r'^\s+([^:]+?)\s*:\s*(.*)$')

# Normalised summary label (parenthetical removed, lower case) -> summary key.
# The parenthetical in e.g. "Known Failures (  8)" is the size of the configured
# expected-failure list, NOT the outcome count, so it is discarded.
_SUMMARY_LABELS = {
  'tests run': 'tests_run',
  'failures': 'failures',
  'warnings': 'warnings',
  'known failures': 'known_failures',
  'known unstable': 'known_unstable',
  'retries used': 'retries_used',
  'stderr output': 'stderr_output',
}

# Summary keys holding an integer; 'retries used' is kept as its raw string.
_SUMMARY_INT_KEYS = ('tests_run', 'failures', 'warnings', 'known_failures',
                     'known_unstable', 'stderr_output')

_PAREN_RE = re.compile(r'\([^)]*\)')


def strip_profile_lines(text):
  """Remove the login profile's stray leading PROFILE lines from text.

  Every remote ``bash -lc`` invocation on the group's server emits a line
  reading exactly ``PROFILE`` on stdout, and ``bash -l`` cannot be dropped
  because the PHENIX environment needs it.  Only LEADING PROFILE lines are
  removed, stopping at the first line that is not one, so a PROFILE line
  inside captured test output survives.

  This is deliberately separate from parsing and is never called implicitly
  by :func:`parse_roster`; callers apply it to captured command stdout.

  Parameters
  ----------
  text : str
      Captured stdout, possibly starting with one or more PROFILE lines.

  Returns
  -------
  str
      The text with leading PROFILE lines removed, or the identical object
      when there are none.
  """
  lines = text.split('\n')
  n_strip = 0
  while n_strip < len(lines) and lines[n_strip].strip() == 'PROFILE':
    n_strip += 1
  if n_strip == 0:
    return text
  return '\n'.join(lines[n_strip:])


def _build_mode(outcome, return_code, attempt, attempts_total):
  """Return the deterministic mode string for one outcome, or None.

  Parameters
  ----------
  outcome : str
      One of the outcome constants.
  return_code : int or None
      Return code seen in the detail block, if any.
  attempt : int
      Attempt on which the result was produced.
  attempts_total : int
      Number of attempts configured/used for that test.

  Returns
  -------
  str or None
      ``None`` for a pass or when nothing stable is known, otherwise a short
      string such as ``'rc=1'`` or ``'rc=1; attempt 2 of 3'``.  Never contains
      a wall time or a path, so it is stable across runs.
  """
  if outcome == PASS:
    return None
  parts = []
  if return_code is not None:
    parts.append('rc=%s' % return_code)
  if attempt > 1:
    parts.append('attempt %d of %d' % (attempt, attempts_total))
  if not parts:
    return None
  return '; '.join(parts)


class TestOutcome(object):
  """The recorded outcome of one named test in a suite run."""

  def __init__(self, command, outcome, wall_time=None, attempt=1,
               attempts_total=1, return_code=None, mode=None,
               failure_text=None, failure_text_conflict=False):
    """Store one test's outcome, computing ``mode`` when it is not supplied.

    Parameters
    ----------
    command : str
        The test command exactly as the suite printed it.
    outcome : str
        One of the module's outcome constants.
    wall_time : float or None, optional
        Seconds reported on the result line; ``None`` when there was none.
    attempt : int, optional
        Attempt that produced this result.  Default ``1``.
    attempts_total : int, optional
        Attempts configured/used for this test.  Default ``1``.
    return_code : int or None, optional
        Return code from the detail block, when one was printed.
    mode : str or None, optional
        Explicit mode string.  When omitted it is derived from ``outcome``,
        ``return_code``, ``attempt`` and ``attempts_total``; see
        :func:`_build_mode`.  The failure text never enters it.
    failure_text : list of str or None, optional
        The captured standard error under the result, one entry per line with
        its four-space indent removed and trailing empty entries dropped;
        ``None`` (the default) when the result printed no ``Standard error:``
        line; :func:`parse_roster` leaves it ``None`` for SKIPPED and MISSING
        entries.  Candidate text: read only by
        :func:`failure_mode_difference`, never for an outcome or a count.
    failure_text_conflict : bool, optional
        ``True`` when the replay block's copy of this result gives a different
        return code or a different failure text.  Default ``False``.
    """
    self.command = command
    self.outcome = outcome
    self.wall_time = wall_time
    self.attempt = attempt
    self.attempts_total = attempts_total
    self.return_code = return_code
    if mode is None:
      mode = _build_mode(outcome, return_code, attempt, attempts_total)
    self.mode = mode
    self.failure_text = failure_text
    self.failure_text_conflict = failure_text_conflict

  def __repr__(self):
    """Return a short unambiguous representation for debugging."""
    return '<TestOutcome %s %s mode=%s>' % (
      self.outcome, self.command, self.mode)


class Roster(object):
  """The set of named tests and outcomes parsed from one suite log.

  Every attribute has a stated origin in :data:`FIELD_ORIGIN`.  ``summary`` is
  the harness counting for itself; ``declared_count``, ``declared_commands``
  and ``nprocs`` come from the harness's header block, whose command STRINGS
  are text out of the tree under test; ``tests``, ``replay_commands`` and
  ``result_line_count`` come from harness-formatted result lines whose status
  token the candidate's exit code and stderr select.  ``not_all_finished`` is
  the one field matched on stripped content rather than at column 0, so
  captured test output containing exactly that line would raise it too.  Apart
  from that flag, and from ``failure_text`` and ``failure_text_conflict`` on
  each entry of ``tests`` (captured standard error, read only to compare
  failure modes), no attribute here is derived from captured test output.
  """

  def __init__(self):
    """Create an empty roster; :func:`parse_roster` fills it in."""
    self.declared_count = None
    self.declared_commands = []
    self.nprocs = None
    self.tests = {}
    self.summary = {
      'tests_run': None,
      'failures': None,
      'warnings': None,
      'known_failures': None,
      'known_unstable': None,
      'stderr_output': None,
      'retries_used': None,
    }
    self.not_all_finished = False
    self.replay_commands = []
    # Number of result lines seen in the streaming section.  This is not
    # len(self.tests): the roster also holds skipped and missing entries, and a
    # repeated command collapses in the dict.  None means "not recorded", in
    # which case reconcile() derives it from the outcomes.
    self.result_line_count = None

  def counts(self):
    """Return a dict mapping every outcome constant to how many tests have it.

    Origin: harness-formatted, candidate-triggered (:data:`FIELD_ORIGIN` keys
    ``counts()[<outcome>]``).  Every value is derived from the streaming
    section rather than from the suite's ``Summary:`` block, and never from
    captured test output: ``SKIPPED`` comes from the harness's column-0
    ``repeated, skipping`` line and ``MISSING`` is derived here from a declared
    command that never got a result line.  ``TIMED_OUT`` is always 0.

    Returns
    -------
    dict
        Every outcome constant is a key, including the ones with zero tests.
    """
    result = {}
    for outcome in ALL_OUTCOMES:
      result[outcome] = 0
    for test in self.tests.values():
      if test.outcome not in result:
        result[test.outcome] = 0
      result[test.outcome] += 1
    return result

  def names(self, outcome):
    """Return the sorted commands whose outcome is ``outcome``.

    Parameters
    ----------
    outcome : str
        One of the outcome constants.

    Returns
    -------
    list of str
        Sorted commands; empty when no test has that outcome.
    """
    return sorted([command for command, test in self.tests.items()
                   if test.outcome == outcome])

  def _result_lines_seen(self):
    """Return how many streaming result lines this roster was built from."""
    if self.result_line_count is not None:
      return self.result_line_count
    return len([test for test in self.tests.values()
                if test.outcome not in (SKIPPED, MISSING)])

  def reconcile(self):
    """Return human-readable discrepancies between the roster and its summary.

    What this compares: the ``Summary:`` block is a harness statement and the
    streaming result lines are harness-formatted, so this is a consistency
    check WITHIN the harness's own reporting.  An empty list means the
    harness's two accounts of the same run agree; it is not evidence about the
    candidate's behaviour, and neither account is a candidate claim.

    Nothing is raised for a mismatch; a mismatch is reported.  An empty list
    means the roster reconciles.  A summary value of ``None`` means the suite
    never printed that line, which is reported as "summary line absent" and is
    never treated as a zero.  (A line that was printed but whose value is not
    an integer also leaves ``None`` and so is reported the same way; the two
    cases are not distinguished.)

    A run that died halfway is reported as well: when the header's declared
    count and the summary's ``Tests run`` are both known and ``Tests run`` is
    the smaller, the shortfall is reported as declared tests that never
    reported a result.  The check is needed because ``parallel.py`` sets
    ``self.finished = len(self.results)`` and prints exactly that as ``Tests
    run``, so in any real log that number agrees with the streaming result
    lines even when the suite died with tests still unreported; only the
    declared count sees the gap.  It is also safe on a complete run: the
    header prints ``len(cmd_list) + len(parallel_list)``, which is the very
    quantity ``parallel.py`` compares ``self.finished`` against for its own
    ``NOT ALL TESTS FINISHED!`` warning.  Nothing is reported when the two are
    equal, when ``Tests run`` is the larger, or when either is unknown -- an
    absent header or summary line is already reported on its own.

    Returns
    -------
    list of str
        One string per discrepancy, in a stable order.  Empty when the roster
        is internally consistent.
    """
    problems = []
    counts = self.counts()
    if self.declared_count is None:
      problems.append(
        'header line absent: "Running <N> tests on <P> processors:" was not '
        'found, so the declared roster is unknown')
    elif self.declared_count != len(self.declared_commands):
      problems.append(
        'declared count is %d but %d declared commands were listed'
        % (self.declared_count, len(self.declared_commands)))
    tests_run = self.summary.get('tests_run')
    if (self.declared_count is not None and tests_run is not None
        and tests_run < self.declared_count):
      problems.append(
        'summary "Tests run" is %d but %d tests were declared: %d declared '
        'test(s) never reported a result, so the run did not finish'
        % (tests_run, self.declared_count, self.declared_count - tests_run))
    checks = [
      ('tests_run', 'Tests run', None, self._result_lines_seen()),
      ('failures', 'Failures', FAIL, counts[FAIL]),
      ('warnings', 'Warnings', WARNING, counts[WARNING]),
      ('known_failures', 'Known Failures', EXPECTED_FAIL,
       counts[EXPECTED_FAIL]),
      ('known_unstable', 'Known Unstable', EXPECTED_UNSTABLE,
       counts[EXPECTED_UNSTABLE]),
    ]
    for key, label, outcome, observed in checks:
      value = self.summary.get(key)
      if value is None:
        problems.append('summary line absent: "%s"' % label)
      elif value != observed:
        if outcome is None:
          problems.append(
            'summary "%s" is %d but %d streaming result lines were parsed'
            % (label, value, observed))
        else:
          problems.append(
            'summary "%s" is %d but %d tests have outcome %r'
            % (label, value, observed, outcome))
    return problems

  def __repr__(self):
    """Return a short unambiguous representation for debugging."""
    return '<Roster declared=%s tests=%d counts=%s>' % (
      self.declared_count, len(self.tests), self.counts())


def _to_lines(text):
  """Return ``text`` as a list of lines, accepting a str or a list of lines.

  Splitting is on '\\n' only (not str.splitlines), so an exotic character such
  as a form feed inside captured test output cannot manufacture a new
  column-0 line.
  """
  if hasattr(text, 'splitlines'):
    return [line.rstrip('\r') for line in text.split('\n')]
  return [line.rstrip('\r\n') for line in text]


def _parse_return_code(raw):
  """Return the int value of a printed return code, or None when it is not one.

  ``parallel.py`` prints ``result.return_code`` with ``%s``, so the text can be
  ``None`` as well as an integer.
  """
  raw = raw.strip()
  if re.match(r'^-?\d+$', raw):
    return int(raw)
  return None


def _find_return_code(lines, start):
  """Return the return code from a result's detail block, or None.

  Scans forward from ``start`` over the indented detail block, stopping at the
  first non-blank column-0 line (the next result line, the terminator, ...).
  The first ``  Return code:`` line wins.
  """
  index = start
  while index < len(lines):
    line = lines[index]
    if line.strip() and not line.startswith(' '):
      return None
    match = _RETURN_CODE_RE.match(line)
    if match:
      return _parse_return_code(match.group(1))
    index += 1
  return None


def _detail_block_end(lines, start):
  """Return the index just past the detail block that begins at ``start``.

  The detail block is the span :func:`_find_return_code` scans: every line
  from ``start`` up to the first non-blank column-0 line (the next result
  line, the terminator, ...), blank lines included.  Returns ``len(lines)``
  when the text ends first.
  """
  index = start
  while index < len(lines):
    line = lines[index]
    if line.strip() and not line.startswith(' '):
      return index
    index += 1
  return index


def _find_failure_text(lines, start):
  """Return the captured standard error in a result's detail block, or None.

  Parameters
  ----------
  lines : list of str
      The whole log, one entry per line.
  start : int
      Index of the first line after the result line.

  Returns
  -------
  list of str or None
      ``None`` when the detail block holds no line exactly
      ``  Standard error:``.  Otherwise, reading from just after the FIRST
      such line and stopping at the first line that does not begin with four
      spaces, every line read with exactly its first four characters
      removed, then with trailing ``''`` entries removed; this can be empty.
  """
  end = _detail_block_end(lines, start)
  index = start
  while index < end and lines[index] != _STDERR_HEADER:
    index += 1
  if index >= end:
    return None
  text = []
  index += 1
  while index < end and lines[index].startswith(_CAPTURED_INDENT):
    text.append(lines[index][len(_CAPTURED_INDENT):])
    index += 1
  while text and text[-1] == '':
    text.pop()
  return text


def _parse_result_line(line):
  """Return (command, outcome, wall_time, attempt, attempts_total) or None.

  ``line`` matches only at column 0 and only when the whole line ends with the
  ``[STATUS] <float>s`` tail and an optional known retry note.
  """
  match = _RESULT_RE.match(line)
  if not match:
    return None
  command = match.group('command')
  outcome = _STATUS_TO_OUTCOME[match.group('status')]
  wall_time = float(match.group('wall_time'))
  attempt = 1
  attempts_total = 1
  if match.group('attempt') is not None:
    attempt = int(match.group('attempt'))
    attempts_total = int(match.group('attempts_total'))
  elif match.group('failed_attempts') is not None:
    # "(failed after N attempts)" reports only N; the test used all of them.
    attempt = int(match.group('failed_attempts'))
    attempts_total = attempt
  return (command, outcome, wall_time, attempt, attempts_total)


def _parse_summary(lines, start, roster):
  """Fill ``roster.summary`` from the indented block after a ``Summary:`` line.

  Reading stops at the first non-blank column-0 line.  Labels are matched by
  their leading words with any parenthetical removed, and the value is
  whatever follows the colon.
  """
  index = start
  while index < len(lines):
    line = lines[index]
    if line.strip() and not line.startswith(' '):
      return
    match = _SUMMARY_ITEM_RE.match(line)
    if match:
      label = _PAREN_RE.sub(' ', match.group(1))
      label = ' '.join(label.split()).lower()
      key = _SUMMARY_LABELS.get(label)
      if key is not None:
        value = match.group(2).strip()
        if key in _SUMMARY_INT_KEYS:
          if re.match(r'^-?\d+$', value):
            roster.summary[key] = int(value)
        else:
          roster.summary[key] = value
    index += 1


def parse_roster(text):
  """Parse suite output text into a :class:`Roster` of named tests.

  Content only: the suite exits 0 even when it reports failures, so this
  function never consults, and must never be given, a process exit status.  It
  performs no file, process or network access.

  Origins of what it returns: the header block and the ``Summary:`` block are
  harness statements, the streaming and replay result lines are
  harness-formatted with a status token that the candidate's exit code and
  stderr select, and indented captured test output is a candidate claim from
  which no outcome and no count here is ever derived.  The one exception is
  the ``not_all_finished`` flag, matched on stripped content rather than at
  column 0, which captured output could therefore raise as well.  Captured
  standard error IS read, into each result's ``failure_text`` and
  ``failure_text_conflict``, for one purpose only: comparing failure modes
  (:func:`failure_mode_difference`).
  :data:`FIELD_ORIGIN` states this per field; the header comment's Trust
  boundary section states the boundary itself and its limits.

  Parsing rules, all anchored on what ``parallel.py`` prints:

  * Everything before the column-0 ``Running <N> tests on <P> processors:``
    header is ignored (the preamble can hold column-0 text that looks like a
    command).  When the header is absent, ``declared_count`` and ``nprocs``
    stay ``None``, ``declared_commands`` stays empty, result lines are still
    parsed, and :meth:`Roster.reconcile` reports the absent header.
  * Declared commands are the following ``"  %s"`` lines, with a trailing
    ``" [Parallel]"`` stripped; the block ends at the first line of another
    shape.
  * The streaming section runs to the first column-0 line of eighty or more
    ``=`` characters.  Result lines are counted ONLY there, because the replay
    block after it prints every warning and every failure a second time, with
    the same retry annotation; in a real 2425-test baseline that is 18 FAIL
    lines for 9 failures.  The replay block is recorded in
    ``replay_commands`` as corroboration and never creates or modifies an
    outcome.  When no ``=`` rule is present (a truncated log) the replay
    header line is used as the end of the streaming section instead.
  * Result lines match at column 0 only and are end-anchored on the
    ``[STATUS] <float>s`` tail, so neither an indented ``[WARNING]`` inside
    captured test output nor a bracketed fragment inside a command can be
    mistaken for one.  No outcome is ever inferred from free text: a test
    named ``tst_reduce_timeout.py`` passes like any other.
  * ``Return code: <n>`` in a result's two-space-indented detail block fills
    ``return_code``.
  * A line exactly ``  Standard error:`` in a streaming result's detail block
    (the first one, if there are several) fills ``failure_text``: the lines
    after it that begin with four spaces, each with exactly those four
    characters removed, stopping at the first line that does not begin with
    four spaces, with trailing empty entries removed.  Without such a line
    ``failure_text`` is ``None``; for SKIPPED and MISSING it is always
    ``None``.  Only the final attempt's standard error is ever in the log,
    and only when the suite's standard error was captured into it.
  * A replay-block result line for a command that has a streaming result
    sets that test's ``failure_text_conflict`` to ``True`` when its own
    detail block gives a different ``return_code`` or a different failure
    text; otherwise it stays ``False``.  Nothing else about the test changes.
  * A retry note fills ``attempt`` and ``attempts_total``.  ``(passed on
    attempt N of M)`` gives both; ``(failed after N attempts)`` reports only
    N, and both are then set to N, since a test that failed after N attempts
    used all of them.
  * ``Test <command> repeated, skipping`` yields outcome ``SKIPPED``.  That
    line is printed while the command list is being de-duplicated, i.e. before
    the header, so it is looked for from the start of the text up to the end
    of the streaming section.  A command that also has a result line keeps the
    result: it did run.
  * A declared command with neither a result line nor a skip line gets outcome
    ``MISSING``.
  * ``TIMED_OUT`` is never produced: no origin in this log reports a timeout.
    ``parallel.py`` has no timeout mechanism, so a hung test prints no result
    line and appears as ``MISSING``.  The constant exists only so roster
    comparison has the right shape for a future suite that does report
    timeouts.

  Parameters
  ----------
  text : str or list of str
      Suite output.  A list of lines is accepted as well as one string.

  Returns
  -------
  Roster
      The parsed roster.  Never raises for an inconsistent log; use
      :meth:`Roster.reconcile` to report inconsistencies.
  """
  lines = _to_lines(text)
  roster = Roster()

  # The suite prints this indented, so compare stripped content.
  for line in lines:
    if line.strip() == _NOT_ALL_FINISHED:
      roster.not_all_finished = True
      break

  # 1. Header.  Everything before it is preamble and is ignored.
  header_index = None
  for index, line in enumerate(lines):
    match = _HEADER_RE.match(line)
    if match:
      header_index = index
      roster.declared_count = int(match.group(1))
      roster.nprocs = match.group(2)
      break

  # 2. Declared roster.
  if header_index is None:
    stream_start = 0
  else:
    index = header_index + 1
    while index < len(lines):
      match = _DECLARED_RE.match(lines[index])
      if not match:
        break
      command = match.group(1)
      if command.endswith(' [Parallel]'):
        command = command[:-len(' [Parallel]')]
      roster.declared_commands.append(command)
      index += 1
    stream_start = index

  # 3. End of the streaming section: the '=' rule, or failing that the replay
  # header, or failing that the end of the text.
  stream_end = len(lines)
  for index in range(stream_start, len(lines)):
    if _TERMINATOR_RE.match(lines[index]):
      stream_end = index
      break
  else:
    for index in range(stream_start, len(lines)):
      if lines[index] == _REPLAY_HEADER:
        stream_end = index
        break

  # 4. Streaming results.
  result_line_count = 0
  streamed = set()
  for index in range(stream_start, stream_end):
    parsed = _parse_result_line(lines[index])
    if parsed is None:
      continue
    command, outcome, wall_time, attempt, attempts_total = parsed
    result_line_count += 1
    streamed.add(command)
    return_code = _find_return_code(lines, index + 1)
    roster.tests[command] = TestOutcome(
      command=command,
      outcome=outcome,
      wall_time=wall_time,
      attempt=attempt,
      attempts_total=attempts_total,
      return_code=return_code,
      failure_text=_find_failure_text(lines, index + 1))
  roster.result_line_count = result_line_count

  # 5. Skipped duplicates.  These are printed before the header, so scan from
  # the start of the text; a command that also produced a result line ran, and
  # the result wins.
  for index in range(0, stream_end):
    match = _SKIP_RE.match(lines[index])
    if not match:
      continue
    command = match.group(1)
    if command not in roster.tests:
      roster.tests[command] = TestOutcome(command=command, outcome=SKIPPED)

  # 6. Declared but never reported.
  for command in roster.declared_commands:
    if command not in roster.tests:
      roster.tests[command] = TestOutcome(command=command, outcome=MISSING)

  # 7. Tail: replay block (corroboration only) and the Summary block.  A
  # replay copy that disagrees with its streaming result about the return
  # code or the failure text marks that test's failure_text_conflict, and
  # changes nothing else.
  summary_seen = False
  for index in range(stream_end, len(lines)):
    line = lines[index]
    parsed = _parse_result_line(line)
    if parsed is not None:
      command = parsed[0]
      roster.replay_commands.append(command)
      if command in streamed:
        test = roster.tests[command]
        if (_find_return_code(lines, index + 1) != test.return_code
            or _find_failure_text(lines, index + 1) != test.failure_text):
          test.failure_text_conflict = True
      continue
    if not summary_seen and line == _SUMMARY_HEADER:
      summary_seen = True
      _parse_summary(lines, index + 1, roster)
  return roster


def _rank_of(roster, command):
  """Return the outcome rank of ``command`` in ``roster``.

  Parameters
  ----------
  roster : Roster
      The roster to look ``command`` up in.
  command : str
      A test command.

  Returns
  -------
  int
      The :data:`OUTCOME_RANK` of the recorded outcome; :data:`ABSENT_RANK`
      when the command is not in the roster at all, and also when its outcome
      is not in :data:`OUTCOME_RANK`, since an unknown outcome is never
      credited with a run.
  """
  test = roster.tests.get(command)
  if test is None:
    return ABSENT_RANK
  return OUTCOME_RANK.get(test.outcome, ABSENT_RANK)


def _line_number(raw):
  """Return a frame's line number as an int, or the digit string itself.

  ``raw`` matched ``\\d+``, so ``int()`` fails only on a string past the
  interpreter's integer-conversion digit limit (4300 digits by default in
  current Python 3 releases, including 3.10.9 and 3.11 here).  The text comes
  from candidate output and must never raise, so it is then kept verbatim,
  where a difference in it is still seen.
  """
  try:
    return int(raw)
  except ValueError:
    return raw


def _frame_file_name(path):
  """Return the part of a frame path after its last ``/`` or ``\\``.

  The path is split as TEXT, not as a filesystem path, so a log written on
  one platform reads the same on another.
  """
  return path[max(path.rfind('/'), path.rfind('\\')) + 1:]


def _split_failure_text(text):
  """Split a failure text into its tracebacks and its other text.

  A traceback starts at an entry exactly ``Traceback (most recent call
  last):``.  Its frame run is the following entries that begin with a space;
  of those, a frame is one matching ``  File "<path>", line <n>, in <func>``,
  and the rest (source lines, caret markers, ``[Previous line repeated
  ...]``) belong to the run without being frames.  The exception line is the
  first entry after the run, whatever it holds (even another traceback
  header); when the text ends first, or that entry is ``''``, the traceback
  is truncated and has no exception line (the ``''`` entry is then other
  text).  Every entry that is not a traceback header, not
  in a frame run and not an exception line is other text, in order.

  Parameters
  ----------
  text : list of str
      A failure text, as in :attr:`TestOutcome.failure_text`.

  Returns
  -------
  tuple
      ``(tracebacks, other, truncated)``: ``tracebacks`` is a list with one
      ``(frames, exception_line)`` per traceback in order, where ``frames``
      is a list of ``(path, line_number, function)`` and ``exception_line``
      is ``None`` for a truncated traceback; ``other`` is the list of other
      entries; ``truncated`` is ``True`` when any traceback is truncated.
  """
  tracebacks = []
  other = []
  truncated = False
  index = 0
  while index < len(text):
    if text[index] != _TRACEBACK_HEADER:
      other.append(text[index])
      index += 1
      continue
    frames = []
    index += 1
    while index < len(text) and text[index].startswith(' '):
      match = _FRAME_RE.match(text[index])
      if match:
        frames.append((match.group('path'),
                       _line_number(match.group('line')),
                       match.group('func')))
      index += 1
    if index >= len(text) or text[index] == '':
      truncated = True
      tracebacks.append((frames, None))
    else:
      tracebacks.append((frames, text[index]))
      index += 1
  return (tracebacks, other, truncated)


def _failure_text_layers(text):
  """Return whether a failure text is truncated, and its comparison layers.

  Parameters
  ----------
  text : list of str
      A failure text, as in :attr:`TestOutcome.failure_text`.

  Returns
  -------
  tuple
      ``(truncated, layers)``.  ``layers`` is a list of ``(kind, value)``
      pairs in comparison order: :data:`KIND_EXCEPTION_TYPES` (one exception
      type or ``None`` per traceback), :data:`KIND_CALL_PATH` (per traceback,
      a list of ``[file name, function]`` per frame),
      :data:`KIND_EXCEPTION_MESSAGE` (the exception lines verbatim),
      :data:`KIND_LINE_NUMBERS` (per traceback, the list of frame line
      numbers) and :data:`KIND_OTHER_TEXT` (the other entries verbatim).
      These layers do not cover every entry (source lines, caret markers and
      frame-path directories, for example, are in none), so they only label
      a difference; :func:`failure_mode_difference` compares the whole texts
      after them.
  """
  tracebacks, other, truncated = _split_failure_text(text)
  exception_types = []
  call_paths = []
  exception_lines = []
  line_numbers = []
  for frames, exception_line in tracebacks:
    exception_type = None
    if exception_line is not None:
      match = _EXCEPTION_TYPE_RE.match(exception_line)
      if match:
        exception_type = match.group(1)
    exception_types.append(exception_type)
    call_paths.append([[_frame_file_name(path), function]
                       for (path, line, function) in frames])
    exception_lines.append(exception_line)
    line_numbers.append([line for (path, line, function) in frames])
  layers = [(KIND_EXCEPTION_TYPES, exception_types),
            (KIND_CALL_PATH, call_paths),
            (KIND_EXCEPTION_MESSAGE, exception_lines),
            (KIND_LINE_NUMBERS, line_numbers),
            (KIND_OTHER_TEXT, other)]
  return (truncated, layers)


def failure_mode_difference(baseline_test, candidate_test):
  """Return how two results of one test differ in failure mode, or None.

  The checks are made in this order, and the first that applies decides:

  1. The outcomes differ: ``(CHANGED, KIND_OUTCOME, baseline outcome,
     candidate outcome)``.
  2. ``[return_code, attempt, attempts_total]`` differ: ``(CHANGED,
     KIND_RETURN_CODE, baseline list, candidate list)``.
  3. Either side has ``failure_text_conflict``: ``(UNKNOWN,
     REASON_REPLAY_CONFLICT)``.
  4. A side has no failure text (``None`` or empty): ``(UNKNOWN,
     REASON_NO_TEXT_EITHER)``, ``REASON_NO_TEXT_BASELINE`` or
     ``REASON_NO_TEXT_CANDIDATE``.
  5. Either failure text has a traceback without an exception line:
     ``(UNKNOWN, REASON_TRUNCATED_TRACEBACK)``.
  6. The first differing layer, in the order exception types, call path,
     exception message, line numbers, other failure text, gives ``(CHANGED,
     kind, baseline value, candidate value)``; each value is a list with one
     element per traceback, in order (other failure text: one per entry).
     See :func:`_split_failure_text` for how a failure text is read.
  7. The two failure texts, as whole lists, differ: ``(CHANGED,
     KIND_REMAINING_TEXT, baseline failure text, candidate failure text)``,
     both lists verbatim.  Reached only when every earlier check and layer is
     equal, that is when the texts differ only where no named layer looks:
     for example a source line, a caret marker, ``[Previous line repeated
     ...]``, the directory part of a frame's path, or a line number spelled
     differently with the same value (``05`` and ``5``).
  8. Otherwise ``None``: the same mode.

  What a result carries.  A CHANGED result names the FIRST check or layer
  above that differs and carries that one's two values only; later checks
  and layers may differ as well and appear nowhere in it (a changed call
  path is reported, and a changed exception message after it is not).  An
  UNKNOWN result carries a reason and no values.  So classifying a reported
  test means reading both tests in full -- outcome, return code, attempts
  and the whole ``failure_text`` of each, and for a replay conflict both
  copies in the log -- never only what the result carries.  If a side has
  no failure text (check 4), there is nothing to read for that side; the
  only place left to look is the log itself, for example the test's
  printed output, if the log kept it.  UNKNOWN is never treated as the
  same mode.

  Limits, stated so that ``None`` is not read as more than it is:

  * The failure text is the captured standard error of the FINAL attempt
    only; ``parallel.py``'s ``run_command()`` discards earlier attempts.
  * ``parallel.py`` prints the ``Standard error:`` block to ``sys.stderr``,
    so a log that did not capture standard error has no failure text, and
    every comparison from it that reaches check 4 is UNKNOWN.
  * Two failure texts that differ in any entry never compare as the same
    mode; the named layers only choose the label.  ``None`` therefore means
    the two ``failure_text`` lists are identical as :func:`parse_roster`
    read them, and reading changes the text beyond removing the harness's
    four-space indent: ``_to_lines`` removes every trailing carriage return
    from each log line (every trailing carriage return and newline when the
    log is given as a list of lines), and trailing empty entries are
    dropped, so standard errors that differ only in those ways are not told
    apart.  Only the first ``Standard error:`` block in a result's detail
    block is read, up to its first line without the four-space indent.
    Standard output is not read at all.

  Parameters
  ----------
  baseline_test : TestOutcome
      The test's result in the baseline run.
  candidate_test : TestOutcome
      The same test's result in the candidate run.

  Returns
  -------
  tuple or None
      ``None``, a 4-tuple starting with :data:`CHANGED` whose second element
      is one of :data:`FAILURE_MODE_KINDS`, or a 2-tuple starting with
      :data:`UNKNOWN` whose second element is one of
      :data:`UNKNOWN_REASONS`.
  """
  if baseline_test.outcome != candidate_test.outcome:
    return (CHANGED, KIND_OUTCOME, baseline_test.outcome,
            candidate_test.outcome)
  baseline_run = [baseline_test.return_code, baseline_test.attempt,
                  baseline_test.attempts_total]
  candidate_run = [candidate_test.return_code, candidate_test.attempt,
                   candidate_test.attempts_total]
  if baseline_run != candidate_run:
    return (CHANGED, KIND_RETURN_CODE, baseline_run, candidate_run)
  if (baseline_test.failure_text_conflict
      or candidate_test.failure_text_conflict):
    return (UNKNOWN, REASON_REPLAY_CONFLICT)
  if not baseline_test.failure_text and not candidate_test.failure_text:
    return (UNKNOWN, REASON_NO_TEXT_EITHER)
  if not baseline_test.failure_text:
    return (UNKNOWN, REASON_NO_TEXT_BASELINE)
  if not candidate_test.failure_text:
    return (UNKNOWN, REASON_NO_TEXT_CANDIDATE)
  baseline_truncated, baseline_layers = _failure_text_layers(
    baseline_test.failure_text)
  candidate_truncated, candidate_layers = _failure_text_layers(
    candidate_test.failure_text)
  if baseline_truncated or candidate_truncated:
    return (UNKNOWN, REASON_TRUNCATED_TRACEBACK)
  for (kind, baseline_value), (_, candidate_value) in zip(baseline_layers,
                                                          candidate_layers):
    if baseline_value != candidate_value:
      return (CHANGED, kind, baseline_value, candidate_value)
  if baseline_test.failure_text != candidate_test.failure_text:
    return (CHANGED, KIND_REMAINING_TEXT, list(baseline_test.failure_text),
            list(candidate_test.failure_text))
  return None


def compare_rosters(baseline, candidate):
  """Return the roster transitions between a baseline run and a candidate run.

  ``blocking`` is DEFINED by one rule over every command in either roster,
  not by a union of the named keys below: a command blocks when its candidate
  outcome ranks BELOW the rank it had to reach, or when it FAILS on the
  candidate without having failed on the baseline.

  The rank is 2 (PASS: ran and passed), 1 (FAIL, WARNING, EXPECTED_FAIL,
  EXPECTED_UNSTABLE: ran, did not pass) or 0 (SKIPPED, MISSING, TIMED_OUT:
  did not run), and a command ABSENT from a roster ranks 0 as well; see
  :data:`OUTCOME_RANK`.  The rank a command had to reach is its BASELINE rank,
  except that a command ABSENT FROM THE BASELINE had to reach PASS, so a newly
  added or newly enabled test that is not passing is a new failure, as it was
  before.  The second clause is ``new_failures`` below: it catches a candidate
  FAIL whose baseline outcome ranks no higher -- WARNING, EXPECTED FAIL,
  EXPECTED UNSTABLE, SKIPPED, MISSING or TIMED_OUT to FAIL -- which rank alone
  cannot see and which no ordering invented among the rank-1 outcomes should
  have to settle.

  So PASS to anything else blocks, including PASS to WARNING, to EXPECTED FAIL
  and to EXPECTED UNSTABLE -- all three of which ``parallel.py`` can produce
  for a test that passed before -- as well as PASS to SKIPPED, MISSING,
  TIMED_OUT or absent.  Anything to FAIL blocks unless the baseline FAILED
  too.  FAIL to MISSING and FAIL to absent block, because coverage was lost.
  FAIL to FAIL, FAIL to PASS, absent-from-baseline to PASS, and a command that
  did not run on EITHER side (MISSING to SKIPPED and the like) do not block:
  measured against the baseline, nothing was lost.  Neither clause enumerates
  blocking outcomes, so no transition goes unnoticed for want of a named
  category.

  Origin: every command this returns -- in the sorted lists, and as the first
  element of each tuple in ``mode_changed``, ``failure_mode_changed`` and
  ``failure_mode_unknown`` -- is a key of the rosters' ``tests``, that is a
  command string the harness printed on a result line, on a ``repeated,
  skipping`` line or in the declared roster.  Membership of every key except
  the three failure-mode keys (``failure_mode_changed``,
  ``failure_mode_unknown``, ``needs_classification``) rests only on what the
  harness printed: on outcomes, on presence in each roster and, for
  ``mode_changed``, on the return code and attempts from a result line and
  its detail block.  For those keys both sides are harness-formatted
  accounts, so a difference between them is evidence about the two RUNS and
  not a claim made by either candidate.  Membership of the three
  failure-mode keys, and what their tuples report, can also rest on captured
  standard error, text written by the test process: on ``failure_text``, and
  on ``failure_text_conflict``, which compares a result's return code and
  ``failure_text`` with the second copy of that result printed at the end
  of the log.  Many results never get a second copy -- every EXPECTED FAIL
  and EXPECTED UNSTABLE result, and WARNING results in a run with no FAIL
  -- and for those this flag is always ``False``.  None of the three
  enters ``blocking``.

  Parameters
  ----------
  baseline : Roster
      Roster parsed from the run of the unmodified code.
  candidate : Roster
      Roster parsed from the run of the changed code.

  Returns
  -------
  dict
      ``blocking`` is the sorted list of commands selected by the rule above.
      Of the other keys only ``new_failures`` enters that rule, as its second
      clause; as a list it is DIAGNOSTIC, like every other key but
      ``blocking``.  The keys below report most changes between the two
      runs.  One kind is missing: a test that starts passing is reported
      only if it was a plain FAIL before (the ``fixed`` key); a test that
      starts passing after being EXPECTED FAIL, MISSING or SKIPPED is not
      reported anywhere.  The three failure-mode keys report tests that
      failed on both sides but failed in a different way.  They can only
      see what is in the captured error text of the final attempt: two
      failures that differ only in ordinary printed output, only in
      trailing carriage returns or blank lines, or only in an earlier
      attempt look the same to these keys (the Limits of
      :func:`failure_mode_difference` state these bounds).
      The established keys keep their established meaning:

      ``new_failures`` (not failing on baseline, or absent there, and FAIL on
      candidate), ``fixed`` (FAIL on baseline, PASS on candidate),
      ``newly_not_passing`` (present in BOTH rosters, PASS on baseline and not
      PASS on candidate -- this is what makes a PASS to WARNING, EXPECTED FAIL
      or EXPECTED UNSTABLE transition legible in a report instead of only
      counted in ``blocking``), ``disappeared`` (in baseline, absent from
      candidate), ``newly_skipped``, ``newly_timed_out``, ``newly_missing``
      (MISSING on candidate and not MISSING on baseline, including when the
      command is absent from the baseline roster) and ``new_tests`` (in
      candidate, not baseline) are sorted lists of commands.  ``mode_changed``
      is a list of ``(command, baseline_mode, candidate_mode)`` tuples, sorted
      by command, for tests that FAILED on both sides with a different mode.

      Three keys compare HOW a test failed, over every command present in
      BOTH rosters whose outcome ranks 1 (FAIL, WARNING, EXPECTED_FAIL,
      EXPECTED_UNSTABLE) on BOTH sides, using
      :func:`failure_mode_difference`, which reads the captured standard
      error.  ``failure_mode_changed`` is a list of ``(command, kind,
      baseline_detail, candidate_detail)`` tuples for the CHANGED results,
      ``failure_mode_unknown`` a list of ``(command, reason)`` tuples for
      the UNKNOWN ones, both sorted by command, and ``needs_classification``
      the sorted list of every command in either.  A ``failure_mode_changed``
      tuple names only the FIRST check or layer that differs and carries only
      that one's two values, though later ones may differ as well, and a
      ``failure_mode_unknown`` tuple carries only a reason; so classifying a
      command in ``needs_classification`` means reading both tests in full
      -- outcome, return code, attempts and the whole ``failure_text`` of
      each, and for a replay conflict both copies in the log -- never only
      its tuple; and if a side has no failure text, there is nothing to
      read for that side, and the only place left to look is the log
      itself (for example the test's printed output, if the log kept it).
      None of the three enters ``blocking``.
      A roster comparison is NOT a clean result while
      ``needs_classification`` is non-empty and its commands have not yet
      been classified: an empty ``blocking`` says nothing was lost,
      not that every failure failed the same way.

      Every command in ``new_failures`` blocks, by the rule's second clause.
      The other diagnostic keys may hold a command that does not block:
      ``new_tests`` holds a newly added test that passes, and ``newly_missing``
      holds a command that was already SKIPPED on the baseline, which lost no
      coverage.  Read ``blocking`` for the verdict and the other keys for
      detail.

      ``candidate.not_all_finished`` is still worth checking as independent
      corroboration: it comes from the suite's own warning line rather than
      from the roster.  That line is matched on stripped content rather than
      at column 0, since the suite prints it indented, so captured test output
      containing exactly that line would raise it too; it can only be raised,
      never cleared.  See its :data:`FIELD_ORIGIN` entry.
  """
  new_failures = []
  fixed = []
  mode_changed = []
  disappeared = []
  newly_skipped = []
  newly_timed_out = []
  newly_missing = []
  newly_not_passing = []
  new_tests = []

  for command, test in candidate.tests.items():
    before = baseline.tests.get(command)
    if before is None:
      new_tests.append(command)
    if test.outcome == FAIL and (before is None or before.outcome != FAIL):
      new_failures.append(command)
    if (test.outcome == PASS and before is not None
        and before.outcome == FAIL):
      fixed.append(command)
    if (test.outcome == FAIL and before is not None
        and before.outcome == FAIL and before.mode != test.mode):
      mode_changed.append((command, before.mode, test.mode))
    if (test.outcome == SKIPPED
        and (before is None or before.outcome != SKIPPED)):
      newly_skipped.append(command)
    if (test.outcome == TIMED_OUT
        and (before is None or before.outcome != TIMED_OUT)):
      newly_timed_out.append(command)
    if (test.outcome == MISSING
        and (before is None or before.outcome != MISSING)):
      newly_missing.append(command)
    if (test.outcome != PASS and before is not None
        and before.outcome == PASS):
      newly_not_passing.append(command)

  for command in baseline.tests:
    if command not in candidate.tests:
      disappeared.append(command)

  new_failures = sorted(new_failures)
  disappeared = sorted(disappeared)
  newly_skipped = sorted(newly_skipped)
  newly_timed_out = sorted(newly_timed_out)
  newly_missing = sorted(newly_missing)
  # The rule, over the union of both rosters.  This is the DEFINITION of
  # blocking.  Clause one: a command must reach the rank it reached on the
  # baseline, and a command the baseline never had must reach PASS.  Clause
  # two: a command that FAILS on the candidate without having failed on the
  # baseline blocks whatever its baseline rank was, which is how WARNING,
  # EXPECTED FAIL, EXPECTED UNSTABLE, SKIPPED, MISSING and TIMED_OUT to FAIL
  # block without inventing an ordering among the rank-1 outcomes.
  new_failure_set = set(new_failures)
  blocking = []
  for command in set(baseline.tests) | set(candidate.tests):
    if command in baseline.tests:
      required_rank = _rank_of(baseline, command)
    else:
      required_rank = OUTCOME_RANK[PASS]
    if (_rank_of(candidate, command) < required_rank
        or command in new_failure_set):
      blocking.append(command)
  blocking = sorted(blocking)

  # How each both-sides rank-1 test failed.  Diagnostic only: nothing here
  # touches blocking or any key above.
  failure_mode_changed = []
  failure_mode_unknown = []
  for command in sorted(set(baseline.tests) & set(candidate.tests)):
    if (_rank_of(baseline, command) != OUTCOME_RANK[FAIL]
        or _rank_of(candidate, command) != OUTCOME_RANK[FAIL]):
      continue
    difference = failure_mode_difference(baseline.tests[command],
                                         candidate.tests[command])
    if difference is None:
      continue
    if difference[0] == CHANGED:
      failure_mode_changed.append((command,) + tuple(difference[1:]))
    else:
      failure_mode_unknown.append((command, difference[1]))
  needs_classification = sorted(
    set([item[0] for item in failure_mode_changed])
    | set([item[0] for item in failure_mode_unknown]))
  return {
    'new_failures': new_failures,
    'fixed': sorted(fixed),
    'mode_changed': sorted(mode_changed, key=lambda item: item[0]),
    'disappeared': disappeared,
    'newly_skipped': newly_skipped,
    'newly_timed_out': newly_timed_out,
    'newly_missing': newly_missing,
    'newly_not_passing': sorted(newly_not_passing),
    'new_tests': sorted(new_tests),
    'blocking': blocking,
    'failure_mode_changed': failure_mode_changed,
    'failure_mode_unknown': failure_mode_unknown,
    'needs_classification': needs_classification,
  }
