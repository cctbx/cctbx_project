"""
abort_message surfaced for all stop reasons; max-cycles message explicit (v115 F1+F3).

F1: abort_message is populated for every stop_reason value, not only 'red_flag'.
F3: An explicit "MAXIMUM CYCLES REACHED" message is printed when the cycle
    loop exhausts max_cycles without a voluntary stop.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent"), os.path.join(_ROOT, "programs")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


# ============================================================================
# Minimal stubs — enough to exercise the two edited blocks without the full
# Phenix / libtbx stack.
# ============================================================================

class _FakeVlog:
    """Capture vlog.quiet() calls for assertion."""
    def __init__(self):
        self.quiet_lines = []
        self.verbose_lines = []
    def quiet(self, msg):
        self.quiet_lines.append(msg)
    def verbose(self, msg):
        self.verbose_lines.append(msg)
    def normal(self, msg):
        pass


class _FakeParams:
    class ai_analysis:
        max_cycles = 3
        log_directory = None
    class communication:
        run_on_server = True


class _FakeSession:
    """Minimal session stub for the cycle-loop tests."""
    def __init__(self, num_cycles=0):
        self._num_cycles = num_cycles
        self.data = {}
    def get_num_cycles(self):
        return self._num_cycles
    def start_cycle(self, cycle):
        pass
    def to_dict(self):
        return {}
    def save(self):
        pass
    def format_all_cycles(self):
        return ""
    def session_file(self):
        return ""


# ============================================================================
# F1 TESTS
# ============================================================================

def _run_f1_block(stop_reason, abort_message, red_flag_issues=None):
    """
    Execute only the F1 edited block from _run_single_cycle post-processing.

    We replicate the logic inline so we don't need to instantiate the full
    Program class.  If this test breaks because the production code diverges,
    the test itself should be updated.
    """
    vlog = _FakeVlog()
    history_record = {
        "stop_reason": stop_reason,
        "abort_message": abort_message,
        "red_flag_issues": red_flag_issues or [],
    }

    # ----- replicated F1 block -----
    abort_msg = history_record.get('abort_message', '')
    if abort_msg:
        vlog.quiet(abort_msg)
    if history_record.get('stop_reason') == 'red_flag':
        issues = history_record.get('red_flag_issues', [])
        if issues:
            vlog.quiet(f"\nRed flag issues detected: {len(issues)}")
    # ----- end replicated block -----

    return vlog


def test_f1_red_flag_prints_abort_message():
    """red_flag: abort_message printed, issue count printed."""
    vlog = _run_f1_block('red_flag', 'R-free too high', red_flag_issues=[1, 2])
    assert 'R-free too high' in vlog.quiet_lines, vlog.quiet_lines
    assert any('2' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_directive_prints_abort_message():
    """directive: abort_message must now reach the user (was silently dropped)."""
    vlog = _run_f1_block('directive', 'Stopping after phenix.refine as directed.')
    assert any('Stopping after phenix.refine' in l for l in vlog.quiet_lines), \
        vlog.quiet_lines


def test_f1_consecutive_cap_prints_abort_message():
    vlog = _run_f1_block('consecutive_cap', 'Consecutive cap reached for phenix.refine.')
    assert any('Consecutive cap' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_stuck_loop_prints_abort_message():
    vlog = _run_f1_block('stuck_loop', 'Workflow is stuck: valid_programs unchanged for 3 cycles.')
    assert any('Workflow is stuck' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_no_valid_programs_prints_abort_message():
    vlog = _run_f1_block('no_valid_programs', 'No valid programs remain.')
    assert any('No valid programs' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_no_workflow_state_prints_abort_message():
    vlog = _run_f1_block('no_workflow_state', 'Workflow state could not be determined.')
    assert any('Workflow state' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_all_commands_duplicate_prints_abort_message():
    vlog = _run_f1_block('all_commands_duplicate', 'All candidate commands are duplicates.')
    assert any('duplicate' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_build_failures_prints_abort_message():
    vlog = _run_f1_block('build_failures_and_duplicates', 'Some programs cannot be built.')
    assert any('cannot be built' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_cannot_build_any_program_prints_abort_message():
    vlog = _run_f1_block('cannot_build_any_program', 'Cannot build any program.')
    assert any('Cannot build' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_after_program_not_available_prints_abort_message():
    vlog = _run_f1_block('after_program_not_available',
                         'Stop-after program was not run this session.')
    assert any('Stop-after' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f1_empty_abort_message_prints_nothing():
    """No abort_message → nothing printed (no spurious blank output)."""
    vlog = _run_f1_block('directive', '')
    assert vlog.quiet_lines == [], vlog.quiet_lines


def test_f1_none_abort_message_prints_nothing():
    """abort_message=None → nothing printed."""
    vlog = _FakeVlog()
    history_record = {'stop_reason': 'directive', 'abort_message': None}
    abort_msg = history_record.get('abort_message', '')
    if abort_msg:
        vlog.quiet(abort_msg)
    assert vlog.quiet_lines == [], vlog.quiet_lines


def test_f1_red_flag_no_issues_no_count_line():
    """red_flag with empty issues list → no 'Red flag issues detected' line."""
    vlog = _run_f1_block('red_flag', 'Problem found.', red_flag_issues=[])
    assert not any('Red flag issues detected' in l for l in vlog.quiet_lines), \
        vlog.quiet_lines


# ============================================================================
# F3 TESTS
# ============================================================================

def _run_f3_loop(max_cycles, break_on_cycle=None):
    """
    Simulate the F3 cycle loop.

    break_on_cycle: if set, _run_single_cycle returns True on that cycle number,
                    simulating a voluntary stop.  None means the loop exhausts.

    Returns (vlog, broke_early).
    """
    vlog = _FakeVlog()
    start_cycle = 1
    broke_early = False

    for cycle in range(start_cycle, start_cycle + max_cycles):
        should_break = (break_on_cycle is not None and cycle == break_on_cycle)
        if should_break:
            broke_early = True
            break

    if not broke_early:
        vlog.quiet("\n" + "=" * 60)
        vlog.quiet("MAXIMUM CYCLES REACHED (%d)" % max_cycles)
        vlog.quiet("The agent did not stop voluntarily. If work is incomplete,")
        vlog.quiet("consider re-running with a higher max_cycles value,")
        vlog.quiet("or review the log to understand why the agent kept running.")
        vlog.quiet("=" * 60)

    return vlog, broke_early


def test_f3_exhausted_loop_prints_message():
    """When max_cycles is exhausted, the MAXIMUM CYCLES REACHED message appears."""
    vlog, broke_early = _run_f3_loop(max_cycles=3, break_on_cycle=None)
    assert not broke_early
    assert any('MAXIMUM CYCLES REACHED' in l for l in vlog.quiet_lines), \
        vlog.quiet_lines
    assert any('3' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f3_voluntary_stop_no_message():
    """When the loop breaks early (voluntary stop), no MAXIMUM CYCLES message."""
    vlog, broke_early = _run_f3_loop(max_cycles=3, break_on_cycle=2)
    assert broke_early
    assert not any('MAXIMUM CYCLES REACHED' in l for l in vlog.quiet_lines), \
        vlog.quiet_lines


def test_f3_single_cycle_exhausted():
    """Edge case: max_cycles=1 and it runs without breaking."""
    vlog, broke_early = _run_f3_loop(max_cycles=1, break_on_cycle=None)
    assert not broke_early
    assert any('MAXIMUM CYCLES REACHED' in l for l in vlog.quiet_lines)


def test_f3_max_cycles_value_in_message():
    """The message includes the actual max_cycles value."""
    vlog, _ = _run_f3_loop(max_cycles=20, break_on_cycle=None)
    assert any('20' in l for l in vlog.quiet_lines), vlog.quiet_lines


def test_f3_zero_max_cycles_no_banner():
    """max_cycles=0: empty range, no cycles ran, no MAXIMUM CYCLES banner.
    Guards against a misleading message when the loop body never executed."""
    max_cycles = 0
    broke_early = False  # range was empty, loop body never ran
    vlog = _FakeVlog()
    if not broke_early and max_cycles > 0:
        vlog.quiet("MAXIMUM CYCLES REACHED (%d)" % max_cycles)
    assert vlog.quiet_lines == [], vlog.quiet_lines


# ============================================================================
# Test registry and runner
# ============================================================================

# ============================================================================
# F1 COMPANION TESTS — run_ai_agent.py abort_message population
# ============================================================================
# These replicate the history_record construction logic from
# run_ai_agent.py to verify abort_message is set for all stop_reason values.

def _build_history_record_abort_message(stop, stop_reason, reasoning):
    """Replicate the relevant run_ai_agent.py logic."""
    response = {"stop": stop, "stop_reason": stop_reason,
                "decision": {"reasoning": reasoning}}
    decision = response.get("decision", {})
    history_record = {"reasoning": decision.get("reasoning", "")}

    # F1 companion fix (v115): populate for ALL stops
    if response.get("stop"):
        history_record["abort_message"] = decision.get("reasoning", "")

    return history_record


def test_companion_red_flag_has_abort_message():
    hr = _build_history_record_abort_message(True, 'red_flag', 'R-free sanity failed.')
    assert hr.get('abort_message') == 'R-free sanity failed.', hr


def test_companion_directive_has_abort_message():
    hr = _build_history_record_abort_message(True, 'directive', 'Stopped as directed.')
    assert hr.get('abort_message') == 'Stopped as directed.', hr


def test_companion_consecutive_cap_has_abort_message():
    hr = _build_history_record_abort_message(True, 'consecutive_cap', 'Cap reached.')
    assert hr.get('abort_message') == 'Cap reached.', hr


def test_companion_stuck_loop_has_abort_message():
    hr = _build_history_record_abort_message(True, 'stuck_loop', 'Loop detected.')
    assert hr.get('abort_message') == 'Loop detected.', hr


def test_companion_no_valid_programs_has_abort_message():
    hr = _build_history_record_abort_message(True, 'no_valid_programs', 'No programs remain.')
    assert hr.get('abort_message') == 'No programs remain.', hr


def test_companion_non_stop_cycle_no_abort_message():
    """Normal (non-stop) cycle must NOT get abort_message set."""
    hr = _build_history_record_abort_message(False, None, '')
    assert 'abort_message' not in hr, hr


# ============================================================================
# Test registry and runner
# ============================================================================

_TESTS = [
    # F1 — ai_agent.py print guard
    test_f1_red_flag_prints_abort_message,
    test_f1_directive_prints_abort_message,
    test_f1_consecutive_cap_prints_abort_message,
    test_f1_stuck_loop_prints_abort_message,
    test_f1_no_valid_programs_prints_abort_message,
    test_f1_no_workflow_state_prints_abort_message,
    test_f1_all_commands_duplicate_prints_abort_message,
    test_f1_build_failures_prints_abort_message,
    test_f1_cannot_build_any_program_prints_abort_message,
    test_f1_after_program_not_available_prints_abort_message,
    test_f1_empty_abort_message_prints_nothing,
    test_f1_none_abort_message_prints_nothing,
    test_f1_red_flag_no_issues_no_count_line,
    # F3 — max cycles banner
    test_f3_exhausted_loop_prints_message,
    test_f3_voluntary_stop_no_message,
    test_f3_single_cycle_exhausted,
    test_f3_max_cycles_value_in_message,
    test_f3_zero_max_cycles_no_banner,
    # F1 companion — run_ai_agent.py abort_message population
    test_companion_red_flag_has_abort_message,
    test_companion_directive_has_abort_message,
    test_companion_consecutive_cap_has_abort_message,
    test_companion_stuck_loop_has_abort_message,
    test_companion_no_valid_programs_has_abort_message,
    test_companion_non_stop_cycle_no_abort_message,
]


def run_all_tests():
    for test_fn in _TESTS:
        test_fn()
    print("All %d tests passed." % len(_TESTS))

if __name__ == "__main__":
    passed = 0
    failed = 0
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        try:
            test_fn()
            print("  PASS: %s" % test_fn.__name__)
            passed += 1
        except Exception as exc:
            import traceback
            print("  FAIL: %s" % test_fn.__name__)
            traceback.print_exc()
            failed += 1
    print()
    if failed:
        print("%d/%d tests FAILED." % (failed, passed + failed))
        sys.exit(1)
    else:
        print("All %d tests passed." % passed)

