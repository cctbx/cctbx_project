"""
Error Classifier for PHENIX AI Agent (Fix 3, v115).

Classifies program failures from log text into categories that
drive the tiered failure response:

  TERMINAL     → immediate pivot, no retry
  PHIL_ERROR   → self-correctable by stripping/fixing params
  LABEL_ERROR  → self-correctable by fixing data labels
  RETRYABLE    → generic retryable (one self-correction attempt)

Each classification includes extracted details (bad parameter
names, error messages, suggested fixes) that the self-correction
prompt can use.

Usage:
    from agent.error_classifier import classify_error

    result = classify_error(log_text, program="phenix.autobuild")
    # result = {
    #   "category": "PHIL_ERROR",
    #   "error_message": "some phil parameters are not recognized",
    #   "bad_params": ["obs_labels"],
    #   "suggestion": "Remove obs_labels from the command",
    #   "is_terminal": False,
    # }
"""

from __future__ import absolute_import, division, print_function

import re
import logging

logger = logging.getLogger(__name__)

# ================================================================
# Error categories
# ================================================================

TERMINAL = "TERMINAL"
PHIL_ERROR = "PHIL_ERROR"
LABEL_ERROR = "LABEL_ERROR"
AMBIGUOUS_PHIL = "AMBIGUOUS_PHIL"
RETRYABLE = "RETRYABLE"
NO_ERROR = "NO_ERROR"

# ================================================================
# Terminal error patterns
# ================================================================

# These errors cannot be resolved by parameter changes.
# The agent should immediately pivot to a different program.

# (A Python traceback is terminal too, unless it ends in a Sorry: see
# _last_traceback.)

_TERMINAL_PATTERNS = [
  # Data/model structural errors
  (re.compile(
    r'polymer\s+crosses\s+special\s+position',
    re.IGNORECASE),
   "Model has polymer crossing special position "
   "(data/model error)"),

  # Fatal space group / symmetry errors
  (re.compile(
    r'(?:cannot|unable to)\s+determine\s+'
    r'(?:space\s+group|symmetry)',
    re.IGNORECASE),
   "Cannot determine space group or symmetry"),

  # No data at all
  (re.compile(
    r'no\s+(?:reflection|diffraction)\s+data\s+found',
    re.IGNORECASE),
   "No reflection data found in file"),

  # License / installation errors
  (re.compile(
    r'(?:license|licence)\s+(?:expired|invalid|not found)',
    re.IGNORECASE),
   "License error (installation issue)"),

  # Memory / resource errors
  (re.compile(
    r'(?:MemoryError|out\s+of\s+memory|cannot\s+allocate)',
    re.IGNORECASE),
   "Out of memory (resource limit)"),

  # Segmentation fault
  (re.compile(
    r'(?:segmentation\s+fault|SIGSEGV|core\s+dumped)',
    re.IGNORECASE),
   "Segmentation fault (program crash)"),

  # Kill signal (whole words: not "unskilled" or "killed_atoms")
  (re.compile(
    r'\b(?:killed|SIGKILL|signal\s+9)\b',
    re.IGNORECASE),
   "Process killed (resource limit or timeout)"),
]

# ================================================================
# The step's own error, and tracebacks
# ================================================================

# The agent writes each step's log in labelled parts; the part after the
# line that starts with this label is the step's own error.
_ERROR_TEXT_LABEL = "ERROR TEXT:"
_STEP_LOG_LABELS = (
  "WORKING DIRECTORY:", "COMMAND THAT WAS RUN:", "LOG TEXT:",
  "ERROR TEXT:", "MTZ LABEL INFO:")

_TRACEBACK_RE = re.compile(
  r'Traceback \(most recent call last\)', re.IGNORECASE)

# The last line of a traceback that reports a Sorry (bad input), with or
# without its module: "Sorry: ...", "libtbx.utils.Sorry: ...",
# "exceptions.Sorry_No_array_of_the_required_type: ..."
# (a Sorry with no message is printed as the bare class name)
_SORRY_EXCEPTION_RE = re.compile(r'(?:[\w.]+\.)?Sorry(?:_\w+)?(?::|$)')

# ================================================================
# PHIL error patterns
# ================================================================

_PHIL_UNRECOGNIZED_RE = re.compile(
  r'(?:sorry[:\s]*)?'
  r'(?:some\s+)?phil\s+parameters?\s+'
  r'(?:are\s+)?not\s+recognized'
  r'(?:\s+by\s+(\w+))?',
  re.IGNORECASE,
)

# Extract specific bad parameter names from PHIL error text.
# PHENIX prints lines like:
#   Unrecognized: obs_labels
#   unrecognized PHIL parameter: "autobuild.input.xray_data.obs_labels"
_PHIL_BAD_PARAM_RE = re.compile(
  r'(?:unrecognized|unknown|invalid)'
  r'(?:\s+phil)?\s*(?:parameter)?'
  r'[:\s]+["\']?([a-zA-Z_][a-zA-Z0-9_.]*)["\']?',
  re.IGNORECASE,
)

# ================================================================
# Ambiguous PHIL patterns
# ================================================================

_PHIL_AMBIGUOUS_RE = re.compile(
  r'(?:sorry[:\s]*)?'
  r'(?:ambiguous|multiple\s+(?:definitions?|matches?))'
  r'(?:\s+for)?\s+'
  r'(?:phil\s+)?parameter\s+'
  r'["\']?([a-zA-Z_]\w*)["\']?',
  re.IGNORECASE,
)

# ================================================================
# Label / column errors
# ================================================================

_LABEL_ERROR_PATTERNS = [
  re.compile(
    r'(?:no\s+array\s+of\s+)?R-free\s+flags?\s+'
    r'(?:found|not\s+found)',
    re.IGNORECASE),
  re.compile(
    r'(?:cannot|could\s+not)\s+(?:find|determine)'
    r'\s+(?:data\s+)?(?:labels?|columns?|arrays?)',
    re.IGNORECASE),
  re.compile(
    r'ambiguous\s+(?:data\s+)?(?:labels?|columns?)',
    re.IGNORECASE),
  re.compile(
    r'multiple\s+(?:data\s+)?arrays?\s+(?:found|present)',
    re.IGNORECASE),
  re.compile(
    r'input\s+labels?\s+(?:error|not\s+found|ambiguous)',
    re.IGNORECASE),
]

# ================================================================
# Sorry pattern (generic PHENIX error)
# ================================================================

# An error line: "Sorry:" at the start of a line (spaces before it
# allowed). Not the word inside a line or a name ("raise_sorry = False").
_SORRY_RE = re.compile(
  r'^[ \t]*Sorry:[ \t]*(.+?)[ \t]*$',
  re.MULTILINE,
)


# ================================================================
# Main classifier
# ================================================================

def classify_error(log_text, program=None):
  """Classify a program failure from its log output.

  Args:
      log_text: Full log text from the failed program.
          Can be empty/None if program produced no output.
      program: Program name (e.g., "phenix.autobuild").
          Optional; used for context in suggestions.

  Returns:
      dict with keys:
        category: One of TERMINAL, PHIL_ERROR, AMBIGUOUS_PHIL,
            LABEL_ERROR, RETRYABLE, NO_ERROR.
        error_message: Human-readable error summary.
        is_terminal: bool — True if pivot immediately.
        bad_params: list of str — param names to strip
            (PHIL_ERROR only).
        ambiguous_param: str — the ambiguous param name
            (AMBIGUOUS_PHIL only).
        suggestion: str — actionable fix suggestion for
            self-correction prompt.

  Never raises.
  """
  if not log_text:
    return _result(NO_ERROR, "No log output to analyze")

  # The step's own error decides whether a failure is final, not
  # everything the program printed (a wizard's output may quote the
  # errors of sub-steps it recovered from). Bad parameter names are
  # still collected from the whole text.
  whole_text = log_text
  output_text = None
  parts = _step_log_parts(log_text)
  if parts is not None:
    output_text, log_text = parts

  # --- 1. Terminal errors (check first — they override everything) ---
  # A traceback is a crash unless it ends in a Sorry (bad input raised
  # inside the agent's own process). Such a traceback is replaced by its
  # Sorry, so that the text of the traceback itself (source lines, an
  # earlier exception of a chain) cannot decide.
  found = _last_traceback(log_text)
  if found is not None:
    start, end, exception_line = found
    m = _SORRY_EXCEPTION_RE.match(exception_line)
    if m is None:
      return _result(
        TERMINAL, "Program crashed with Python traceback",
        is_terminal=True,
        suggestion="Pivot to a different program. "
          "This error cannot be resolved by "
          "parameter changes.")
    message = exception_line[m.end():].strip() or "(no message)"
    log_text = "%sSorry: %s\n%s" % (
      log_text[:start], message, log_text[end:])
  for pattern, description in _TERMINAL_PATTERNS:
    if pattern.search(log_text):
      return _result(
        TERMINAL, description,
        is_terminal=True,
        suggestion="Pivot to a different program. "
          "This error cannot be resolved by "
          "parameter changes.")

  # --- 2-4. Errors the agent can correct itself ---
  correctable = _classify_correctable(log_text, whole_text, program)
  if correctable is not None:
    return correctable

  # --- 5. Generic Sorry error (retryable) ---
  sorry = _extract_sorry(log_text)
  if sorry:
    return _result(
      RETRYABLE,
      sorry,
      suggestion="Review the error and adjust "
        "parameters or try a different approach.")

  # A step log whose own error says nothing specific (for example "The
  # program exited with code 1"): the program may have printed a
  # parameter or label complaint in its output. That is not a decision
  # to give up on the program, so the output may be read for it.
  if output_text is not None:
    correctable = _classify_correctable(output_text, whole_text, program)
    if correctable is not None:
      return correctable

  # --- 6. Check for generic failure indicators ---
  if re.search(r'(?:FAILED|ERROR|Error|failed)', log_text):
    return _result(
      RETRYABLE,
      "Program reported an error",
      suggestion="Review the log and adjust parameters.")

  # --- 7. A step log with an error part that no rule above knows ---
  # The planner classifies only a step that did not end in SUCCESS, so
  # this is a failure, and it must not become NO_ERROR: the planner
  # counts failures only for other categories (should_pivot). (The agent
  # writes the label for any text on the error stream, warnings of a
  # step that worked included; such a step is not classified.)
  if parts is not None:
    first_line = "The step failed with no error message"
    if log_text.strip():
      first_line = log_text.strip().splitlines()[0].strip()
    if len(first_line) > 200:
      first_line = first_line[:200] + "..."
    return _result(
      RETRYABLE,
      first_line,
      suggestion="Review the error and adjust "
        "parameters or try a different approach.")

  # No error detected
  return _result(NO_ERROR, "No error detected in log")


def _classify_correctable(log_text, whole_text, program):
  """Return the result for an error the agent can correct itself (an
  unrecognized or ambiguous PHIL parameter, a data label error) found in
  log_text, or None. Bad parameter names come from whole_text."""
  # --- 2. Unrecognized PHIL parameters ---
  m = _PHIL_UNRECOGNIZED_RE.search(log_text)
  if m:
    # Extract specific bad param names
    bad_params = _extract_bad_params(whole_text)
    target = m.group(1) or (
      program.replace("phenix.", "") if program else "program")
    return _result(
      PHIL_ERROR,
      "Unrecognized PHIL parameter(s) for %s" % target,
      bad_params=bad_params,
      suggestion="Remove the unrecognized parameter(s): "
        "%s" % ", ".join(bad_params) if bad_params
        else "Remove parameters not recognized by %s"
        % target)

  # --- 3. Ambiguous PHIL parameters ---
  m = _PHIL_AMBIGUOUS_RE.search(log_text)
  if m:
    ambiguous = m.group(1)
    return _result(
      AMBIGUOUS_PHIL,
      "Ambiguous PHIL parameter: %s" % ambiguous,
      ambiguous_param=ambiguous,
      suggestion="Use the fully qualified parameter "
        "name instead of '%s'" % ambiguous)

  # --- 4. Label / column errors ---
  for pattern in _LABEL_ERROR_PATTERNS:
    if pattern.search(log_text):
      error_msg = _extract_sorry(log_text)
      return _result(
        LABEL_ERROR,
        error_msg or "Data label/column error",
        suggestion="Check data labels. "
          "May need to specify obs_labels or "
          "generate R-free flags.")
  return None


def _result(category, error_message, is_terminal=None,
            bad_params=None, ambiguous_param=None,
            suggestion=None):
  """Build a classification result dict."""
  if is_terminal is None:
    is_terminal = (category == TERMINAL)
  return {
    "category": category,
    "error_message": error_message,
    "is_terminal": is_terminal,
    "bad_params": bad_params or [],
    "ambiguous_param": ambiguous_param or "",
    "suggestion": suggestion or "",
  }


def _step_log_parts(log_text):
  """Split the text of a step log: (what the program printed, the step's
  own error). None for a text that is not a step log.

  The agent writes a step log in labelled parts (_STEP_LOG_LABELS). The
  step's own error is what follows the last line that starts with
  "ERROR TEXT:", up to the next label; what the program printed is the
  text before that line.
  """
  lines = log_text.splitlines(True)
  start = None
  for i, line in enumerate(lines):
    if line.startswith(_ERROR_TEXT_LABEL):
      start = i
  if start is None:
    return None
  error_lines = [lines[start][len(_ERROR_TEXT_LABEL):]]
  for line in lines[start + 1:]:
    if line.startswith(_STEP_LOG_LABELS):
      break
    error_lines.append(line)
  return "".join(lines[:start]), "".join(error_lines)


def _last_traceback(text):
  """Find the last Python traceback in text.

  Return None if there is none, or (start, end, exception_line): start
  is where the first traceback of text begins, end is just after the
  line that names the exception of the last one, and exception_line is
  that line ("" if the traceback was cut off before it).

  A traceback whose exception is a Sorry (_SORRY_EXCEPTION_RE) reports
  bad input, raised inside the agent's own process and logged with the
  traceback of the agent's launcher; it is not a crash of the program.
  """
  first = None
  last = None
  for m in _TRACEBACK_RE.finditer(text):
    if first is None:
      first = m
    last = m
  if last is None:
    return None
  position = last.start()
  header = True
  for line in text[last.start():].splitlines(True):
    position += len(line)
    if header: # the "Traceback (most recent call last):" line
      header = False
    elif line.strip() and not line[0].isspace():
      return first.start(), position, line.strip()
  return first.start(), len(text), ""


def _extract_bad_params(log_text):
  """Extract unrecognized parameter names from log text.

  Returns:
      list of str — parameter names found in error messages.
  """
  params = []
  for m in _PHIL_BAD_PARAM_RE.finditer(log_text):
    name = m.group(1).strip()
    # Filter out common false positives
    if name and name not in ("sorry", "Sorry", "Error",
                              "error", "by", "the"):
      # Extract the leaf name (last component)
      leaf = name.rsplit(".", 1)[-1]
      if leaf not in params:
        params.append(leaf)
  return params


def _extract_sorry(log_text):
  """Extract the first Sorry: error message from log text.

  Returns:
      str or None
  """
  m = _SORRY_RE.search(log_text)
  if m:
    msg = m.group(1).strip()
    # Truncate long messages
    if len(msg) > 200:
      msg = msg[:200] + "..."
    return msg
  return None


# ================================================================
# Failure tracking helpers
# ================================================================

def count_consecutive_failures(history, program):
  """Count consecutive recent failures of a specific program.

  Walks backward through history from the most recent entry.
  Stops counting when it finds a different program or a
  success.

  Args:
      history: List of history dicts (oldest first).
      program: Program name to count failures for.

  Returns:
      int — number of consecutive recent failures of this
      program (0 if last entry is not a failure of this
      program).
  """
  if not history or not program:
    return 0

  count = 0
  prog_short = program.replace("phenix.", "")

  for h in reversed(history):
    if not isinstance(h, dict):
      continue
    h_prog = (h.get("program") or "").replace("phenix.", "")
    h_result = (h.get("result") or h.get("status") or "")

    if h_prog != prog_short:
      break  # Different program — stop counting
    if str(h_result).startswith("SUCCESS"):
      break  # This program succeeded — stop counting
    # This program failed
    count += 1

  return count



# Per-program session block thresholds (P4 fix).
# refine/rsr are expected to run many times; use a higher threshold.
# All other programs use the default.
_SESSION_BLOCK_THRESHOLD = {
    "refine": 6,
    "real_space_refine": 6,
    "__default__": 4,
}


def count_total_failures(history, program):
  """Count all failures of a program in the session, regardless of order.

  Used by the session_blocked_programs mechanism (P4 fix).
  Catches the case where a program fails with different parameters each
  time (Fix 5 allows each because the command string differs, but the
  aggregate still crosses the session threshold).

  Args:
      history: List of history dicts (oldest first).
      program: Program name to count.

  Returns:
      int -- total failures of this program in the session.
  """
  if not history or not program:
    return 0
  prog_short = program.replace("phenix.", "")
  count = 0
  for h in history:
    if not isinstance(h, dict):
      continue
    h_prog = (h.get("program") or "").replace("phenix.", "")
    if h_prog != prog_short:
      continue
    h_result = str(h.get("result") or h.get("status") or "")
    if not h_result.startswith("SUCCESS"):
      count += 1
  return count


def session_block_threshold(program):
  """Return the session-block threshold for a given program.

  Args:
      program: Program name (with or without phenix. prefix).

  Returns:
      int -- number of total failures before session block fires.
  """
  prog_short = program.replace("phenix.", "")
  return _SESSION_BLOCK_THRESHOLD.get(
      prog_short, _SESSION_BLOCK_THRESHOLD["__default__"])


def should_pivot(history, program, error_classification):
  """Decide whether to pivot away from a failed program.

  Implements the tiered failure response from the PLAN:
    - Terminal error → always pivot
    - First failure → allow self-correction (return False)
    - Second failure of same program → pivot

  Args:
      history: List of history dicts.
      program: The program that just failed.
      error_classification: Dict from classify_error().

  Returns:
      (should_pivot, reason) tuple.
      should_pivot: bool.
      reason: str explanation for logging.
  """
  if not error_classification:
    return False, ""

  cat = error_classification.get("category", NO_ERROR)

  if cat == NO_ERROR:
    return False, ""

  if cat == TERMINAL:
    return True, (
      "Terminal error: %s"
      % error_classification.get("error_message", ""))

  # Count consecutive failures of this program
  failures = count_consecutive_failures(history, program)

  if failures >= 2:
    return True, (
      "%s failed %d consecutive times -- pivoting "
      "to different program"
      % (program, failures))

  # Check total session failures (P4 fix).
  # This catches the case where a program fails with different parameters
  # each time (consecutive count resets, but total keeps climbing).
  total = count_total_failures(history, program)
  threshold = session_block_threshold(program)
  if total >= threshold:
    return True, (
      "%s failed %d times total (session threshold=%d) "
      "-- session block"
      % (program, total, threshold))

  # First/early failure -- allow self-correction attempt
  return False, (
    "%s failed once (%s) -- allowing self-correction"
    % (program, cat))
