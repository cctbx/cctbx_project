"""
Phase 4: History Flag Consistency.

Tests that flag names written by _analyze_history / build_context
match what detect_step / get_valid_programs read.

Bug B1 (last_program) found and fixed in this session.

Run: python tests/tst_phase4_history_flags.py
Produces: findings/phase_4_flag_mismatches.yaml
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys

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

FINDINGS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'findings')


def assert_true(condition, msg=""):
    if not condition:
        raise AssertionError(msg or "Assertion failed")


def assert_in(item, container, msg=""):
    if item not in container:
        raise AssertionError(
            "%s: %r not in %r" % (msg, item, container)
            if msg else "%r not in container" % item)


def assert_equal(a, b, msg=""):
    if a != b:
        raise AssertionError(
            "%s: %r != %r" % (msg, a, b)
            if msg else "%r != %r" % (a, b))


# =====================================================================
# CHECK 1: All flags read by detect_step/get_valid_programs are
#           written by build_context
# =====================================================================

def test_all_read_flags_are_written():
    """Every flag read from context must be written by build_context."""
    print("Test: all_read_flags_are_written")

    with open(os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))),
            'agent', 'workflow_engine.py')) as f:
        src = f.read()

    # Extract flags from dict literal in build_context
    ctx_start = src.find('context = {')
    ctx_end = src.find('\n        }', ctx_start) + 10
    ctx_block = src[ctx_start:ctx_end]
    literal_flags = set()
    for m in re.finditer(r'"(\w+)":', ctx_block):
        literal_flags.add(m.group(1))

    # Flags from auto-include (done_flag_map)
    from knowledge.program_registration import (
        get_program_done_flag_map)
    auto_flags = set(get_program_done_flag_map().values())

    # Post-literal context assignments
    post_flags = set()
    after_literal = src[ctx_end:]
    for m in re.finditer(r'context\["(\w+)"\]\s*=',
                         after_literal):
        post_flags.add(m.group(1))

    all_written = literal_flags | auto_flags | post_flags

    # All flags read
    all_read = set()
    for m in re.finditer(r'context\.get\(["\'](\w+)["\']',
                         src):
        all_read.add(m.group(1))
    for m in re.finditer(r'context\[["\'](\w+)["\']\]', src):
        all_read.add(m.group(1))

    # Remove false positives (comments, strings)
    false_positives = {'has_a', 'has_b'}  # Only in docstring
    all_read -= false_positives

    read_not_written = all_read - all_written

    if read_not_written:
        print("  WARNING: %d flags read but not written:"
              % len(read_not_written))
        for f in sorted(read_not_written):
            print("    %s" % f)
    else:
        print("  All read flags are written")

    # All gaps should now be resolved (B1 fix added last_program)
    known_gaps = set()
    unknown_gaps = read_not_written - known_gaps
    assert_true(len(unknown_gaps) == 0,
        "Unexpected flags read but not written: %s"
        % sorted(unknown_gaps))
    print("  PASSED (all read flags are written)")


# =====================================================================
# CHECK 2: last_program IS in context (B1 fix verified)
# =====================================================================

def test_last_program_in_context():
    """Verify last_program is transferred from history_info to context.
    (B1 fix: was missing before v115.08.)"""
    print("Test: last_program_in_context")

    from agent.workflow_state import _analyze_history

    # Create minimal history with a successful refine
    history = [
        {"program": "phenix.xtriage",
         "result": "SUCCESS: OK", "cycle": 1},
        {"program": "phenix.refine",
         "result": "SUCCESS: R-work=0.22 R-free=0.26",
         "cycle": 2},
    ]
    info = _analyze_history(history)

    # last_program set in history_info
    assert_equal(info.get("last_program"), "phenix.refine",
        "last_program should be set in history_info")

    # Verify build_context copies it
    from agent.workflow_engine import WorkflowEngine
    engine = WorkflowEngine()
    context = engine.build_context(
        files={"data_mtz": [], "model": [], "pdb": []},
        history_info=info)

    assert_equal(context.get("last_program"), "phenix.refine",
        "last_program should be in context (B1 fix)")
    print("  PASSED")


# =====================================================================
# CHECK 3: done_tracking YAML flags all exist in program_registration
# =====================================================================

def test_yaml_done_flags_registered():
    """Every done_tracking.flag in programs.yaml should be in
    get_program_done_flag_map()."""
    print("Test: yaml_done_flags_registered")
    import yaml

    with open(os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))),
            'knowledge', 'programs.yaml')) as f:
        progs = yaml.safe_load(f)

    from knowledge.program_registration import (
        get_program_done_flag_map)
    registered = set(get_program_done_flag_map().values())

    yaml_flags = {}
    for prog_name, prog in progs.items():
        if not isinstance(prog, dict):
            continue
        dt = prog.get('done_tracking', {})
        flag = dt.get('flag')
        if flag:
            yaml_flags[prog_name] = flag

    missing = []
    for prog, flag in yaml_flags.items():
        if flag not in registered:
            missing.append((prog, flag))

    if missing:
        print("  YAML done flags not in registration:")
        for prog, flag in missing:
            print("    %s → %s" % (prog, flag))
    assert_equal(len(missing), 0,
        "All YAML done_tracking flags should be registered")
    print("  PASSED (%d YAML done flags all registered)"
          % len(yaml_flags))


# =====================================================================
# CHECK 4: Post-probe correction clears validation_done
# =====================================================================

def test_post_probe_clears_validation_done():
    """When model_vs_data runs as placement probe (before refine),
    validation_done should be cleared by post-probe correction."""
    print("Test: post_probe_clears_validation_done")

    from agent.workflow_state import _analyze_history

    # Scenario: xtriage → model_vs_data (probe), no refine
    history = [
        {"program": "phenix.xtriage",
         "result": "SUCCESS: OK", "cycle": 1},
        {"program": "phenix.model_vs_data",
         "result": "SUCCESS: R-work=0.386 R-free=0.386",
         "cycle": 2},
    ]
    info = _analyze_history(history)

    assert_true(info.get("placement_probed"),
        "placement_probed should be True after model_vs_data")
    assert_true(not info.get("validation_done"),
        "validation_done should be cleared (no refine before "
        "model_vs_data)")
    print("  PASSED")


def test_post_probe_keeps_validation_after_refine():
    """When model_vs_data runs AFTER refine, validation_done stays."""
    print("Test: post_probe_keeps_validation_after_refine")

    from agent.workflow_state import _analyze_history

    # Scenario: xtriage → refine → model_vs_data (real validation)
    history = [
        {"program": "phenix.xtriage",
         "result": "SUCCESS: OK", "cycle": 1},
        {"program": "phenix.refine",
         "result": "SUCCESS: R-work=0.22 R-free=0.26",
         "cycle": 2},
        {"program": "phenix.model_vs_data",
         "result": "SUCCESS: R-work=0.22 R-free=0.26",
         "cycle": 3},
    ]
    info = _analyze_history(history)

    assert_true(info.get("validation_done"),
        "validation_done should stay True (model_vs_data ran "
        "after refine)")
    print("  PASSED")


# =====================================================================
# CHECK 5: Program → done_flag mapping is consistent
# =====================================================================

def test_done_flag_set_by_correct_program():
    """Each done flag should only be set by its registered program."""
    print("Test: done_flag_set_by_correct_program")

    from agent.workflow_state import _analyze_history

    # Test each major program sets its own done flag
    test_cases = [
        ("phenix.xtriage", "xtriage_done"),
        ("phenix.refine", "refine_done"),
        ("phenix.phaser", "phaser_done"),
        ("phenix.autobuild", "autobuild_done"),
        ("phenix.real_space_refine", "rsr_done"),
    ]

    for prog, expected_flag in test_cases:
        history = [
            {"program": prog,
             "result": "SUCCESS: OK",
             "cycle": 1},
        ]
        info = _analyze_history(history)
        assert_true(info.get(expected_flag),
            "%s should set %s" % (prog, expected_flag))

    print("  PASSED (%d program→flag mappings verified)"
          % len(test_cases))


# =====================================================================
# CHECK 6: refine_count increments correctly
# =====================================================================

def test_refine_count_increments():
    """refine_count should count successful refine cycles."""
    print("Test: refine_count_increments")

    from agent.workflow_state import _analyze_history

    history = [
        {"program": "phenix.xtriage",
         "result": "SUCCESS: OK", "cycle": 1},
        {"program": "phenix.refine",
         "result": "SUCCESS: R-work=0.25 R-free=0.30",
         "cycle": 2},
        {"program": "phenix.refine",
         "result": "SUCCESS: R-work=0.22 R-free=0.26",
         "cycle": 3},
        {"program": "phenix.refine",
         "result": "FAILED: crash",
         "cycle": 4},
    ]
    info = _analyze_history(history)

    assert_equal(info.get("refine_count", 0), 2,
        "refine_count should be 2 (two successes, one failure)")
    print("  PASSED")


# =====================================================================
# CHECK 7: Dead flags (written but never read)
# =====================================================================

def test_dead_flags_documented():
    """Flags written to context but never read are dead code.
    Document them but don't fail — they're harmless."""
    print("Test: dead_flags_documented")

    base = os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))

    # Scan both files that read context flags
    all_read = set()
    for module in ['agent/workflow_engine.py',
                    'agent/graph_nodes.py']:
        with open(os.path.join(base, module)) as f:
            src = f.read()
        for m in re.finditer(
                r'context\.get\(["\'](\w+)["\']', src):
            all_read.add(m.group(1))
        for m in re.finditer(
                r'context\[["\'](\w+)["\']\]', src):
            all_read.add(m.group(1))

    # Extract all written flags from workflow_engine.py
    with open(os.path.join(base,
              'agent/workflow_engine.py')) as f:
        src = f.read()
    ctx_start = src.find('context = {')
    ctx_end = src.find('\n        }', ctx_start) + 10
    ctx_block = src[ctx_start:ctx_end]
    literal_flags = set()
    for m in re.finditer(r'"(\w+)":', ctx_block):
        literal_flags.add(m.group(1))

    from knowledge.program_registration import (
        get_program_done_flag_map)
    auto_flags = set(get_program_done_flag_map().values())
    post_flags = set()
    after_literal = src[ctx_end:]
    for m in re.finditer(r'context\["(\w+)"\]\s*=',
                         after_literal):
        post_flags.add(m.group(1))
    all_written = literal_flags | auto_flags | post_flags

    dead = all_written - all_read
    print("  Dead flags (written, never read in "
          "workflow_engine or graph_nodes): %d" % len(dead))
    for f in sorted(dead):
        print("    %s" % f)
    # Note: _check_conditions reads flags dynamically from
    # YAML (context.get("has_" + val)), which the regex
    # can't see. Currently no "dead" flags are read this way.
    print("  PASSED (informational — %d dead flags)" % len(dead))


# =====================================================================
# RUN ALL TESTS
# =====================================================================

def run_all_tests():
    print("Phase 4: History Flag Consistency")
    print()

    tests = [
        test_all_read_flags_are_written,
        test_last_program_in_context,
        test_yaml_done_flags_registered,
        test_post_probe_clears_validation_done,
        test_post_probe_keeps_validation_after_refine,
        test_done_flag_set_by_correct_program,
        test_refine_count_increments,
        test_dead_flags_documented,
    ]

    passed = 0
    failed = 0
    failures = []

    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            failures.append({
                'test': t.__name__,
                'error': str(e)[:200],
            })
            print("  FAILED: %s" % str(e)[:100])

    print()
    print("  Results: %d passed, %d failed out of %d" %
          (passed, failed, len(tests)))

    status = 'PASS' if failed == 0 else 'FAIL'

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_4_flag_mismatches.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 4: History Flag Consistency Findings\n")
        f.write("phase: 4\n")
        f.write("name: history_flags\n")
        f.write("status: %s\n" % status)
        f.write("total_tests: %d\n" % len(tests))
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        f.write("\nbugs:\n")
        f.write("  - id: BUG_last_program\n")
        f.write("    status: FIXED\n")
        f.write("    severity: LOW\n")
        f.write("    description: \"last_program was missing "
                "from build_context().\n")
        f.write("      Fixed by adding: last_program = "
                "history_info.get('last_program')\"\n")
        if failures:
            f.write("\ntest_failures:\n")
            for fail in failures:
                f.write("  - test: %s\n" % fail['test'])
                f.write("    error: \"%s\"\n"
                        % fail['error'][:100])

    print("  Findings: %s" % path)
    print("  Phase 4 overall: %s" % status)

    result = {'status': status, 'passed': passed,
              'failed': failed}
    if failed > 0:
        raise AssertionError(
            "Phase 4 history flags FAILED: %d/%d tests failed"
            % (failed, passed + failed))
    return result


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] == 'PASS' else 1)
