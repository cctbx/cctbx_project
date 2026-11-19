"""
Test runner return-value handling and _read_sequence encoding.

Scope: pure logic plus one subprocess.  No network, no keys.

Covers:
    - H1 tests/run_all_tests.py _result_failed truth table: None, True, 0,
        {'status': 'PASS'}, {'status': 'PARTIAL', 'failed': 2},
        {'failed': 0} and {} pass; False, 1, 2, {'status': 'FAIL'} and
        {'failed': 1} fail.  The runner is loaded from its path with
        importlib (main() is not run).  ({'status': 'PASS', 'failed': 3}
        moved to H8: it now fails.)
    - H2 run_test_module(name, f) returns success False when f returns False,
        returns 1, raises AssertionError or returns {'status': 'FAIL'}; True when
        f returns None or 0.
    - H3 agent/graph_nodes.py _read_sequence reads a FASTA file under
        -X warn_default_encoding -W error::EncodingWarning (subprocess) and
        returns 'MKVLLA'; without an explicit encoding the warning-as-error is
        swallowed and '' is returned.
    - H4 run_all_tests.py references tst_claude5_readiness and
        tst_provider_reasoning_models, and both expose a callable
        run_all_tests.
    - H5 run_all_tests.py _module_failure_counters truth table, using a
        synthetic module registered in sys.modules: FAIL=2 and
        FAILURES=['x'] are flagged; FAIL=0, FAILURES=[], FAIL=True,
        FAIL='FAIL' and no counters are not.
    - H6 run_test_module returns success False when the run function
        returns None but its module has FAIL=1; True when FAIL=0.
    - H7 tst_directive_layer_diagnostics.run_all_tests returns False (not
        True) when phenix's ai_agent.py is not found (AI_AGENT_PATH set to
        None, then restored).
    - H8 _result_failed dict cases: status ERROR, FAILED, 'failed' (any
        case) fail; PASS or OK with a nonzero 'failed' count fail; PARTIAL
        with 'failed': 2, PASS with 'failed': 0 and bare PASS pass.
    - H9 run_test_module with run functions calling sys.exit(1),
        sys.exit(0), sys.exit(None) and sys.exit("message"): no SystemExit
        escapes and success is False, True, True, False.
    - H10 AST scan of every production .py file under libtbx/langchain
        (all subdirectories except tests/): each open() / io.open() call
        whose mode is a literal text mode (absent, or without 'b') passes
        encoding=.  Non-literal modes and **kwargs calls are skipped.
        tst_file_encoding's production scan covers only agent/, knowledge/
        and utils/ (and is line-regex based); H10 also checks that it
        reached analysis/, core/ and log_extraction/, where 8 such calls
        were missed.
"""
from __future__ import absolute_import, division, print_function

import contextlib
import importlib.util
import io
import os
import shutil
import subprocess
import ast
import sys
import tempfile
import types

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_LANGCHAIN_DIR = os.path.dirname(_TESTS_DIR)
# Tree root that contains the libtbx package this test belongs to.
_TREE_ROOT = os.path.dirname(os.path.dirname(_LANGCHAIN_DIR))
_RUNNER_PATH = os.path.join(_TESTS_DIR, "run_all_tests.py")


def _load_module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_runner_cache = []


def _runner():
    if not _runner_cache:
        _runner_cache.append(
            _load_module_from_path("_tst_runner_under_test", _RUNNER_PATH))
    return _runner_cache[0]


def test_result_failed_truth_table():
    result_failed = _runner()._result_failed
    cases = [
        (None, False),
        (True, False),
        (False, True),
        (0, False),
        (1, True),
        (2, True),
        ({"status": "FAIL"}, True),
        ({"status": "PASS"}, False),
        ({"status": "PARTIAL", "failed": 2}, False),
        ({"failed": 1}, True),
        ({"failed": 0}, False),
        ({}, False),
    ]
    for value, expected in cases:
        got = result_failed(value)
        assert bool(got) == expected, (
            "_result_failed(%r) = %r, expected %r" % (value, got, expected))


def _raise_assertion():
    raise AssertionError("deliberate failure")


def test_run_test_module_return_values():
    run_test_module = _runner().run_test_module
    cases = [
        ("returns_false", lambda: False, False),
        ("returns_one", lambda: 1, False),
        ("raises_assertion", _raise_assertion, False),
        ("returns_status_fail", lambda: {"status": "FAIL"}, False),
        ("returns_none", lambda: None, True),
        ("returns_zero", lambda: 0, True),
    ]
    for name, func, expected in cases:
        with contextlib.redirect_stdout(io.StringIO()):
            outcome = run_test_module(name, func)
        success = outcome[0]
        assert success is expected, (
            "run_test_module(%s) success = %r, expected %r"
            % (name, success, expected))


_SUBPROCESS_CODE = """
import sys
sys.path.insert(0, sys.argv[1])
import libtbx.langchain.agent.graph_nodes as graph_nodes
sys.stdout.write("RESULT=" + repr(graph_nodes._read_sequence(sys.argv[2])))
"""


def test_read_sequence_explicit_encoding():
    tmpdir = tempfile.mkdtemp(prefix="tst_read_sequence_")
    try:
        fasta = os.path.join(tmpdir, "seq.fa")
        with open(fasta, "w", encoding="utf-8") as f:
            f.write(">chain A\nMKV\nLLA\n")
        cmd = [sys.executable, "-X", "warn_default_encoding",
                      "-W", "error::EncodingWarning",
                      "-c", _SUBPROCESS_CODE, _TREE_ROOT, fasta]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                                                    stderr=subprocess.PIPE, cwd=tmpdir)
        out = proc.stdout.decode("utf-8", "replace")
        err = proc.stderr.decode("utf-8", "replace")
        assert proc.returncode == 0, (
            "subprocess failed (rc=%d):\n%s" % (proc.returncode, err[-2000:]))
        assert "RESULT='MKVLLA'" in out, (
            "_read_sequence under EncodingWarning-as-error gave %r" % out)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_runner_registers_new_modules():
    with open(_RUNNER_PATH, "r", encoding="utf-8") as f:
        source = f.read()
    for mod in ("tst_claude5_readiness", "tst_provider_reasoning_models"):
        assert mod in source, "run_all_tests.py does not reference %s" % mod
        module = _load_module_from_path(
            "_tst_loaded_" + mod, os.path.join(_TESTS_DIR, mod + ".py"))
        run = getattr(module, "run_all_tests", None)
        assert callable(run), "%s has no callable run_all_tests" % mod


_SYNTHETIC_MODULE = "_tst_synthetic_counter_module"


@contextlib.contextmanager
def _synthetic_module(_run_source="def run():\n    return None\n",
                      **globals_):
    """Register a throwaway module holding globals_ and a run function.

    _run_source defines run() (default: returns None).

    Yields the module's run function (whose __module__ names the module);
    the module is removed from sys.modules afterwards.
    """
    module = types.ModuleType(_SYNTHETIC_MODULE)
    for name, value in globals_.items():
        setattr(module, name, value)
    exec(_run_source, module.__dict__)
    assert module.run.__module__ == _SYNTHETIC_MODULE
    saved = sys.modules.get(_SYNTHETIC_MODULE)
    sys.modules[_SYNTHETIC_MODULE] = module
    try:
        yield module.run
    finally:
        if saved is None:
            sys.modules.pop(_SYNTHETIC_MODULE, None)
        else:
            sys.modules[_SYNTHETIC_MODULE] = saved


def test_module_failure_counters_truth_table():
    counters_fn = getattr(_runner(), "_module_failure_counters", None)
    assert callable(counters_fn), (
        "run_all_tests.py has no _module_failure_counters")
    cases = [
        ({"FAIL": 2}, True),
        ({"FAIL": 0}, False),
        ({"FAILURES": ["x"]}, True),
        ({"FAILURES": []}, False),
        ({"FAIL": True}, False),
        ({"FAIL": "FAIL"}, False),
        ({}, False),
    ]
    for globals_, expected in cases:
        with _synthetic_module(**globals_) as run:
            got = counters_fn(run)
        assert bool(got) == expected, (
            "_module_failure_counters with %r = %r, expected flagged=%r"
            % (globals_, got, expected))
    assert _SYNTHETIC_MODULE not in sys.modules


def test_run_test_module_counts_module_failures():
    run_test_module = _runner().run_test_module
    for fail, expected in ((1, False), (0, True)):
        with _synthetic_module(FAIL=fail) as run:
            with contextlib.redirect_stdout(io.StringIO()):
                outcome = run_test_module("counter_fail_%d" % fail, run)
        assert outcome[0] is expected, (
            "run_test_module with FAIL=%d and None returned: success = %r, "
            "expected %r" % (fail, outcome[0], expected))
    assert _SYNTHETIC_MODULE not in sys.modules


def test_directive_layer_diagnostics_missing_ai_agent_fails():
    module = _load_module_from_path(
        "_tst_loaded_tst_directive_layer_diagnostics",
        os.path.join(_TESTS_DIR, "tst_directive_layer_diagnostics.py"))
    saved = (module.AI_AGENT_PATH, module._AI_AGENT_HOW)
    module.AI_AGENT_PATH = None
    module._AI_AGENT_HOW = ["(test: forced missing)"]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            result = module.run_all_tests()
    finally:
        module.AI_AGENT_PATH, module._AI_AGENT_HOW = saved
    assert result is False, (
        "tst_directive_layer_diagnostics.run_all_tests() without ai_agent.py "
        "returned %r, expected False" % (result,))


def test_result_failed_dict_status_cases():
    result_failed = _runner()._result_failed
    cases = [
        ({"status": "ERROR"}, True),
        ({"status": "FAILED"}, True),
        ({"status": "failed"}, True),
        ({"status": "PASS", "failed": 3}, True),
        ({"status": "OK", "failed": 1}, True),
        ({"status": "PARTIAL", "failed": 2}, False),
        ({"status": "PASS", "failed": 0}, False),
        ({"status": "PASS"}, False),
    ]
    wrong = []
    for value, expected in cases:
        got = result_failed(value)
        if bool(got) != expected:
            wrong.append("_result_failed(%r) = %r, expected %r"
                         % (value, got, expected))
    assert not wrong, "\n".join(wrong)


def test_run_test_module_catches_system_exit():
    run_test_module = _runner().run_test_module
    cases = [
        ("exit_1", "1", False),
        ("exit_0", "0", True),
        ("exit_none", "None", True),
        ("exit_message", "'message'", False),
    ]
    got = []
    for name, arg, _ in cases:
        source = "import sys\ndef run():\n    sys.exit(%s)\n" % arg
        with _synthetic_module(_run_source=source) as run:
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    outcome = run_test_module(name, run)
            except SystemExit as e:
                raise AssertionError(
                    "run_test_module(%s) let SystemExit(%r) escape"
                    % (name, e.code))
        got.append(outcome[0])
    expected = [case[2] for case in cases]
    assert got == expected, (
        "run_test_module success for sys.exit(1), (0), (None), ('message') "
        "= %r, expected %r" % (got, expected))
    source = "import sys\ndef run():\n    sys.exit(0)\n"
    with _synthetic_module(FAIL=1, _run_source=source) as run:
        with contextlib.redirect_stdout(io.StringIO()):
            outcome = run_test_module("exit_0_with_fail_counter", run)
    assert outcome[0] is False, (
        "run_test_module success for sys.exit(0) with FAIL=1 = %r, "
        "expected False" % (outcome[0],))
    source = "def run():\n    raise KeyboardInterrupt\n"
    with _synthetic_module(_run_source=source) as run:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                run_test_module("keyboard_interrupt", run)
        except KeyboardInterrupt:
            pass
        else:
            raise AssertionError(
                "run_test_module swallowed KeyboardInterrupt")
    assert _SYNTHETIC_MODULE not in sys.modules


def _is_open_call(node):
    func = node.func
    if isinstance(func, ast.Name):
        return func.id == "open"
    return (isinstance(func, ast.Attribute) and func.attr == "open"
            and isinstance(func.value, ast.Name) and func.value.id == "io")


def _text_open_without_encoding(node):
    """True if node is open()/io.open() in literal text mode, no encoding."""
    if not _is_open_call(node):
        return False
    if any(k.arg is None for k in node.keywords):
        return False  # **kwargs: cannot tell
    mode = node.args[1] if len(node.args) > 1 else None
    for k in node.keywords:
        if k.arg == "mode":
            mode = k.value
    if mode is not None:
        if not (isinstance(mode, ast.Constant)
                and isinstance(mode.value, str)):
            return False  # non-literal mode: skipped
        if "b" in mode.value:
            return False
    if len(node.args) > 3:
        return False  # encoding passed positionally
    return not any(k.arg == "encoding" for k in node.keywords)


def test_production_text_opens_pass_encoding():
    offenders = []
    scanned_dirs = set()
    n_files = 0
    for root, dirs, files in os.walk(_LANGCHAIN_DIR):
        rel_root = os.path.relpath(root, _LANGCHAIN_DIR)
        if rel_root == ".":
            dirs[:] = [d for d in dirs if d != "tests"]
        dirs[:] = sorted(d for d in dirs if d != "__pycache__")
        for fname in sorted(files):
            if not fname.endswith(".py") or fname.startswith("._"):
                continue
            path = os.path.join(root, fname)
            rel = os.path.relpath(path, _LANGCHAIN_DIR)
            with open(path, "rb") as f:
                tree = ast.parse(f.read(), filename=path)
            n_files += 1
            scanned_dirs.add(rel.split(os.sep)[0] if os.sep in rel else ".")
            for node in ast.walk(tree):
                if (isinstance(node, ast.Call)
                        and _text_open_without_encoding(node)):
                    offenders.append("%s:%d" % (rel, node.lineno))
    # The scan must reach the directories tst_file_encoding does not.
    missing = {"agent", "analysis", "core", "log_extraction"} - scanned_dirs
    assert not missing, "scan did not reach %s" % sorted(missing)
    assert n_files >= 50, "only %d production files scanned" % n_files
    assert not offenders, (
        "%d text-mode open() calls without encoding= in production code:"
        "\n  %s" % (len(offenders), "\n  ".join(offenders)))


_TESTS = [
    test_result_failed_truth_table,
    test_run_test_module_return_values,
    test_read_sequence_explicit_encoding,
    test_runner_registers_new_modules,
    test_module_failure_counters_truth_table,
    test_run_test_module_counts_module_failures,
    test_directive_layer_diagnostics_missing_ai_agent_fails,
    test_result_failed_dict_status_cases,
    test_run_test_module_catches_system_exit,
    test_production_text_opens_pass_encoding,
]


def run_all_tests():
    """Run every test in _TESTS, fail-fast."""
    for test_fn in _TESTS:
        print("  Running %s..." % test_fn.__name__)
        test_fn()
        print("  PASS: %s" % test_fn.__name__)


if __name__ == "__main__":
    run_all_tests()
    print()
    print("All %d tests passed." % len(_TESTS))
    print("OK")
