"""
Tests for Fix 3: Error classifier + failure tracking (v115).

Verifies error classification from program log text and the
tiered failure response (terminal → pivot, first failure →
self-correct, second failure → pivot).

Run with: python3 tests/tst_error_classifier.py
"""

from __future__ import absolute_import, division, print_function
import os
import sys

PROJECT_ROOT = os.path.dirname(
  os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

try:
  from tests.tst_utils import (
    assert_equal, assert_true, assert_false,
    assert_in, run_tests_with_fail_fast,
  )
except ImportError:
  def assert_equal(a, b, msg=""):
    assert a == b, "%s: %r != %r" % (msg, a, b)
  def assert_true(val, msg=""):
    assert val, msg
  def assert_false(val, msg=""):
    assert not val, msg
  def assert_in(needle, haystack, msg=""):
    assert needle in haystack, "%s: %r not in %r" % (
      msg, needle, haystack)
  def run_tests_with_fail_fast():
    g = globals()
    tests = sorted(k for k in g if k.startswith("test_"))
    for name in tests:
      print("  Running %s..." % name)
      g[name]()
      print("  PASS: %s" % name)
    print("\nAll %d tests passed." % len(tests))


from agent.error_classifier import (
  classify_error, count_consecutive_failures,
  should_pivot,
  TERMINAL, PHIL_ERROR, AMBIGUOUS_PHIL,
  LABEL_ERROR, RETRYABLE, NO_ERROR,
)


# ===========================================================
# Terminal error tests
# ===========================================================

def test_terminal_traceback():
  """Python traceback is terminal."""
  log = (
    "Running phenix.autobuild...\n"
    "Traceback (most recent call last):\n"
    "  File \"autobuild.py\", line 42\n"
    "    raise RuntimeError('something broke')\n"
    "RuntimeError: something broke"
  )
  r = classify_error(log, "phenix.autobuild")
  assert_equal(r["category"], TERMINAL,
    "Traceback should be TERMINAL")
  assert_true(r["is_terminal"],
    "is_terminal should be True")


def test_terminal_polymer_crosses():
  """'polymer crosses special position' is terminal."""
  log = (
    "Sorry: polymer crosses special position — "
    "cannot refine this model"
  )
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], TERMINAL,
    "Polymer crosses special position is TERMINAL")


def test_terminal_segfault():
  """Segmentation fault is terminal."""
  log = "Segmentation fault (core dumped)"
  r = classify_error(log, "phenix.autobuild")
  assert_equal(r["category"], TERMINAL,
    "Segfault should be TERMINAL")


def test_terminal_oom():
  """Out of memory is terminal."""
  log = "MemoryError: unable to allocate array"
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], TERMINAL,
    "OOM should be TERMINAL")


# ===========================================================
# PHIL error tests
# ===========================================================

def test_phil_unrecognized():
  """Unrecognized PHIL parameter detected."""
  log = (
    "Sorry: some phil parameters are not recognized "
    "by autobuild\n"
    "  Unrecognized: obs_labels\n"
    "  Unrecognized: xray_data.obs_labels"
  )
  r = classify_error(log, "phenix.autobuild")
  assert_equal(r["category"], PHIL_ERROR,
    "Unrecognized PHIL should be PHIL_ERROR")
  assert_false(r["is_terminal"],
    "PHIL_ERROR is not terminal")
  assert_in("obs_labels", r["bad_params"],
    "Should extract obs_labels as bad param")


def test_phil_multiple_bad_params():
  """Multiple unrecognized params extracted."""
  log = (
    "Sorry: some phil parameters are not recognized\n"
    "  Unrecognized: fake_param\n"
    "  Unrecognized: another_bad"
  )
  r = classify_error(log)
  assert_equal(r["category"], PHIL_ERROR)
  assert_true(len(r["bad_params"]) >= 2,
    "Should extract 2+ bad params, got %d"
    % len(r["bad_params"]))


def test_plan_case_autobuild_obs_labels():
  """PLAN: AF_exoV autobuild failed with obs_labels."""
  log = (
    "phenix.autobuild data=exoV_data.mtz "
    "obs_labels=I(+) resolution=3.0\n\n"
    "Sorry: some phil parameters are not recognized "
    "by autobuild\n"
    "  Unrecognized: \"autobuild.input.xray_data"
    ".obs_labels\""
  )
  r = classify_error(log, "phenix.autobuild")
  assert_equal(r["category"], PHIL_ERROR)
  assert_true(
    any("obs_labels" in p for p in r["bad_params"]),
    "Should extract obs_labels")


# ===========================================================
# Ambiguous PHIL tests
# ===========================================================

def test_phil_ambiguous():
  """Ambiguous PHIL parameter detected."""
  log = (
    "Sorry: ambiguous parameter rmsd\n"
    "  Matches: process_predicted_model.rmsd\n"
    "  Matches: process_predicted_model.maximum_rmsd"
  )
  r = classify_error(log, "phenix.process_predicted_model")
  assert_equal(r["category"], AMBIGUOUS_PHIL,
    "Ambiguous PHIL should be AMBIGUOUS_PHIL")
  assert_equal(r["ambiguous_param"], "rmsd",
    "Should extract 'rmsd' as ambiguous param")


def test_plan_case_ppm_rmsd():
  """PLAN: AF_exoV_MRSAD PPM failed with ambiguous rmsd."""
  log = (
    "Sorry: multiple definitions for phil parameter "
    "rmsd\n"
    "  maximum_rmsd and rmsd both match"
  )
  r = classify_error(
    log, "phenix.process_predicted_model")
  assert_equal(r["category"], AMBIGUOUS_PHIL)
  assert_equal(r["ambiguous_param"], "rmsd")


# ===========================================================
# Label error tests
# ===========================================================

def test_label_error_rfree():
  """Missing R-free flags detected."""
  log = "Sorry: No array of R-free flags found in data"
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], LABEL_ERROR,
    "Missing R-free should be LABEL_ERROR")
  assert_false(r["is_terminal"],
    "Label errors are not terminal")


def test_label_error_ambiguous():
  """Ambiguous data labels detected."""
  log = "Sorry: ambiguous data labels in input MTZ"
  r = classify_error(log)
  assert_equal(r["category"], LABEL_ERROR)


# ===========================================================
# Retryable / no-error tests
# ===========================================================

def test_generic_sorry():
  """Generic Sorry error is retryable."""
  log = "Sorry: unable to process this file correctly"
  r = classify_error(log)
  assert_equal(r["category"], RETRYABLE,
    "Generic Sorry should be RETRYABLE")


def test_no_error():
  """Successful log has no error."""
  log = (
    "phenix.refine completed successfully\n"
    "R-free = 0.2500\n"
    "Final R-work = 0.2100"
  )
  r = classify_error(log)
  assert_equal(r["category"], NO_ERROR,
    "Successful log should be NO_ERROR")


def test_empty_log():
  """Empty log returns NO_ERROR."""
  r = classify_error("")
  assert_equal(r["category"], NO_ERROR)
  r2 = classify_error(None)
  assert_equal(r2["category"], NO_ERROR)


# ===========================================================
# Priority tests (terminal overrides PHIL)
# ===========================================================

def test_traceback_overrides_phil():
  """Traceback + PHIL error → TERMINAL wins."""
  log = (
    "Sorry: some phil parameters are not recognized\n"
    "Traceback (most recent call last):\n"
    "  File 'x.py', line 1\n"
    "RuntimeError: crash"
  )
  r = classify_error(log)
  assert_equal(r["category"], TERMINAL,
    "Terminal should override PHIL_ERROR")


# ===========================================================
# count_consecutive_failures tests
# ===========================================================

def test_count_zero_no_failures():
  """No failures returns 0."""
  history = [
    {"program": "phenix.autobuild",
     "result": "SUCCESS: OK"},
  ]
  assert_equal(
    count_consecutive_failures(
      history, "phenix.autobuild"),
    0, "Success should count as 0 failures")


def test_count_one_failure():
  """One failure returns 1."""
  history = [
    {"program": "phenix.autobuild",
     "result": "SUCCESS: OK"},
    {"program": "phenix.autobuild",
     "result": "FAILED: Sorry: something"},
  ]
  assert_equal(
    count_consecutive_failures(
      history, "phenix.autobuild"),
    1, "One failure should count as 1")


def test_count_multiple_failures():
  """Multiple consecutive failures counted."""
  history = [
    {"program": "phenix.xtriage",
     "result": "SUCCESS: OK"},
    {"program": "autobuild",
     "result": "FAILED: error 1"},
    {"program": "autobuild",
     "result": "FAILED: error 2"},
    {"program": "autobuild",
     "result": "FAILED: error 3"},
  ]
  assert_equal(
    count_consecutive_failures(
      history, "phenix.autobuild"),
    3, "Three failures should count as 3")


def test_count_different_program_breaks():
  """Different program in between resets count."""
  history = [
    {"program": "autobuild",
     "result": "FAILED: error"},
    {"program": "refine",
     "result": "SUCCESS: OK"},
    {"program": "autobuild",
     "result": "FAILED: error"},
  ]
  assert_equal(
    count_consecutive_failures(
      history, "phenix.autobuild"),
    1, "Only last run counts — refine in between")


def test_count_empty_history():
  """Empty history returns 0."""
  assert_equal(
    count_consecutive_failures([], "phenix.autobuild"),
    0)


# ===========================================================
# should_pivot tests
# ===========================================================

def test_pivot_terminal():
  """Terminal errors always pivot."""
  classification = {
    "category": TERMINAL,
    "error_message": "Traceback",
    "is_terminal": True,
  }
  pivot, reason = should_pivot(
    [], "phenix.autobuild", classification)
  assert_true(pivot, "Terminal should always pivot")
  assert_in("Terminal", reason,
    "Reason should mention terminal")


def test_pivot_first_failure():
  """First failure allows self-correction."""
  classification = {
    "category": PHIL_ERROR,
    "error_message": "bad params",
    "is_terminal": False,
  }
  # No prior failures
  history = [
    {"program": "xtriage", "result": "SUCCESS: OK"},
  ]
  pivot, reason = should_pivot(
    history, "phenix.autobuild", classification)
  assert_false(pivot,
    "First failure should NOT pivot")
  assert_in("self-correction", reason,
    "Reason should mention self-correction")


def test_pivot_second_failure():
  """Second consecutive failure forces pivot."""
  classification = {
    "category": PHIL_ERROR,
    "error_message": "bad params",
    "is_terminal": False,
  }
  # Two prior failures
  history = [
    {"program": "autobuild",
     "result": "FAILED: error 1"},
    {"program": "autobuild",
     "result": "FAILED: error 2"},
  ]
  pivot, reason = should_pivot(
    history, "phenix.autobuild", classification)
  assert_true(pivot,
    "Second failure should pivot")
  assert_in("2 consecutive", reason,
    "Reason should mention consecutive failures")


def test_pivot_no_error():
  """No error → no pivot."""
  classification = {
    "category": NO_ERROR,
    "error_message": "",
    "is_terminal": False,
  }
  pivot, reason = should_pivot(
    [], "phenix.autobuild", classification)
  assert_false(pivot, "No error should not pivot")


# ===========================================================
# Entry point
# ===========================================================

# ===========================================================
# The step's own error decides (the ERROR TEXT part of a step log)
# ===========================================================

def _step_log(log_text, error_text):
  """A step log as the agent writes it."""
  return (
    "WORKING DIRECTORY: /tmp/run\n"
    "COMMAND THAT WAS RUN: phenix.autosol data.sca\n"
    "LOG TEXT: %s\n"
    "ERROR TEXT: %s\n" % (log_text, error_text))


def test_phrase_in_output_is_not_terminal():
  """A terminal phrase in what the program printed is not the step's
  error (AutoSol quotes a refinement message and carries on)."""
  log = _step_log(
    "NOTE: ERROR MESSAGE FROM REFINEMENT:\n"
    "Polymer crosses special position element:\n"
    "NOTE: Refinement against original data did not work..."
    "using working refinement\n"
    "Traceback (most recent call last):\n"
    "  File 'sub.py', line 1\n"
    "RuntimeError: in a nested job, recovered\n"
    "3 waters killed\n",
    "Sorry: Not enough heavy-atom sites found")
  r = classify_error(log, "phenix.autosol")
  assert_equal(r["category"], RETRYABLE,
    "Phrases in the output must not make the failure terminal")
  assert_in("Not enough heavy-atom sites", r["error_message"])


def test_terminal_in_error_text():
  """The same phrase in the step's own error is terminal."""
  log = _step_log(
    "Starting refinement\n",
    "Sorry: Polymer crosses special position element:")
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], TERMINAL,
    "A terminal phrase in the step's error is TERMINAL")


def test_traceback_ending_in_sorry_is_not_a_crash():
  """An input error raised inside the agent's process is logged with
  the traceback of the agent's launcher, ending in a Sorry."""
  log = (
    "WORKING DIRECTORY: /tmp/run\n"
    "COMMAND THAT WAS RUN: phenix.xtriage data.mtz\n"
    "ERROR TEXT: No reflection arrays available.\n"
    "Traceback (most recent call last):\n"
    "  File \"ai_agent.py\", line 6881, in _try_native_execution\n"
    "    launcher()\n"
    "exceptions.Sorry_No_array_of_the_required_type: "
    "No reflection arrays available.\n")
  r = classify_error(log, "phenix.xtriage")
  assert_false(r["is_terminal"],
    "A traceback that ends in a Sorry is not a crash")
  assert_equal(r["category"], RETRYABLE,
    "The Sorry of the traceback is the error")
  assert_equal(r["error_message"], "No reflection arrays available.")
  for last_line in ["Sorry: no data", "libtbx.utils.Sorry: no data",
                    "libtbx.utils.Sorry"]:
    r = classify_error(
      "Traceback (most recent call last):\n"
      "  File 'x.py', line 1\n" + last_line)
    assert_false(r["is_terminal"], last_line)


def test_plain_error_text_is_still_an_error():
  """A failed step whose own error has no known word is RETRYABLE, so
  that the failure is counted; the words in its output do not matter."""
  log = _step_log(
    "Refining\nR-free = 0.25\n",
    "The program exited with code 1")
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], RETRYABLE,
    "A plain error text must still be an error")
  assert_equal(r["error_message"], "The program exited with code 1")
  # a text with no ERROR TEXT part and no error is still NO_ERROR
  r = classify_error("Refining\nR-free = 0.25\n")
  assert_equal(r["category"], NO_ERROR)


def test_last_error_text_part_is_the_steps_own():
  """An 'ERROR TEXT:' line inside the output (the log of a nested
  job) is not the step's own error; the last one is."""
  log = _step_log(
    "nested log follows\n"
    "ERROR TEXT: Segmentation fault in a nested job\n"
    "recovered\n",
    "Sorry: Not enough sites")
  r = classify_error(log, "phenix.autosol")
  assert_equal(r["category"], RETRYABLE)
  assert_equal(r["error_message"], "Not enough sites")


def test_error_part_ends_at_the_next_label():
  """Label information written after the error is not the error."""
  log = _step_log(
    "Refining\n",
    "Sorry: Model and data do not match\n"
    "MTZ LABEL INFO: multiple data arrays found; "
    "no reflection data found in column FOBS")
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], RETRYABLE,
    "Label information must not decide the class")
  assert_equal(r["error_message"], "Model and data do not match")


def test_traceback_text_does_not_decide():
  """A traceback that ends in a Sorry is classed by that Sorry, not
  by the source lines it shows."""
  traceback = (
    "Traceback (most recent call last):\n"
    "  File \"example.py\", line 3, in run\n"
    "    raise Sorry(\"No data\" if no_data else\n"
    "                \"Polymer crosses special position\")\n")
  r = classify_error(traceback + "libtbx.utils.Sorry: No data\n")
  assert_equal(r["category"], RETRYABLE,
    "A phrase in a traceback's source lines must not be terminal")
  assert_equal(r["error_message"], "No data")
  # the same phrase as the message of the Sorry is terminal
  r = classify_error(
    traceback + "libtbx.utils.Sorry: Polymer crosses special position\n")
  assert_equal(r["category"], TERMINAL)
  # a Sorry with no message is still an error
  r = classify_error(traceback + "libtbx.utils.Sorry\n")
  assert_equal(r["category"], RETRYABLE,
    "A Sorry with no message must be RETRYABLE")


def test_step_log_with_empty_error_part_is_an_error():
  """The agent writes the label only for a failed step, so a step log
  with nothing after it is still a failure to count."""
  log = _step_log("All done\n", "")
  r = classify_error(log, "phenix.refine")
  assert_equal(r["category"], RETRYABLE,
    "An empty error part must not be NO_ERROR")


def test_correctable_error_found_in_output():
  """With an error part that says nothing specific, a label or
  parameter complaint in the program's output is still found."""
  log = _step_log(
    "Reading data\nmultiple data arrays found in data.mtz\n",
    "The program exited with code 1")
  r = classify_error(log, "phenix.polder")
  assert_equal(r["category"], LABEL_ERROR,
    "A label complaint in the output should give LABEL_ERROR")
  # but a specific Sorry in the error part comes first
  log = _step_log(
    "Reading data\nmultiple data arrays found in data.mtz\n",
    "Sorry: Selection is empty")
  r = classify_error(log, "phenix.polder")
  assert_equal(r["category"], RETRYABLE)
  assert_equal(r["error_message"], "Selection is empty")


def test_killed_as_a_whole_word():
  """'killed' inside another word or name is not a kill."""
  for text in ["12 unskilled guesses", "killed_atoms=3"]:
    r = classify_error(text)
    assert_false(r["is_terminal"], text)
  for text in ["Killed", "Killed: 9", "received SIGKILL"]:
    r = classify_error(text)
    assert_equal(r["category"], TERMINAL, text)


def test_sorry_as_an_error_line():
  """A parameter name or the ordinary word is not a Sorry error."""
  for text in ["    raise_sorry = False\n    background = True",
               "we are sorry to say this took a while"]:
    r = classify_error(text)
    assert_equal(r["category"], NO_ERROR, text)
  # only as Phenix prints it
  r = classify_error("SORRY: no data")
  assert_equal(r["category"], NO_ERROR)
  r = classify_error("Done\n  Sorry: no data found in file\n")
  assert_equal(r["category"], RETRYABLE)
  assert_equal(r["error_message"], "no data found in file")


def run_all_tests():
  run_tests_with_fail_fast()


if __name__ == "__main__":
  run_all_tests()
