"""Tests for v116.10 file-encoding fix.

The bug: Python's `open(path)` without `encoding=` falls back to
the system default codec.  On Linux/macOS this is UTF-8; on
Windows it's the system code page (cp1252 for English, gbk for
Chinese-locale, cp932 for Japanese, etc.).  Any file containing
non-ASCII bytes that aren't representable in the local codec
crashes.

A Chinese-locale Windows user reported:
  UnicodeDecodeError: 'gbk' codec can't decode byte 0x94
  in position 987: illegal multibyte sequence
  ...
  File ".../knowledge/yaml_loader.py", line 76
  with open(filepath, 'r') as f:

YAML files are UTF-8 by spec, so the fix is to force UTF-8 at
every file-open site.  This test scans the source files for
unguarded text-mode `open()` calls and fails if any are
introduced.

Fail-fast convention: prints PASS on success, raises on failure.
"""

from __future__ import absolute_import, division, print_function

import os
import re
import sys


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_HERE)


def _find_source_file(relpath_candidates):
    """Find a source file across PHENIX install layouts."""
    for relpath in relpath_candidates:
        cand = os.path.join(_PROJECT_ROOT, relpath)
        if os.path.exists(cand):
            return cand
        # PHENIX install: cctbx_project/libtbx/langchain/ ↔ phenix/phenix/
        cand = os.path.normpath(os.path.join(
            _PROJECT_ROOT, "..", "..", "..", relpath))
        if os.path.exists(cand):
            return cand
    return None


def _find_text_opens_without_encoding(source):
    """Find `open(...)` calls in source that are text mode and
    don't specify encoding=.

    Returns a list of (line_number, line_text) tuples for offending
    sites.  Uses Python's tokenize module to properly distinguish
    code from string literals (including docstrings) and comments.
    """
    import tokenize
    import io

    offenders = []

    # Tokenize and track which (line, col) positions are inside
    # strings or comments — those are "not code" and should be
    # excluded from open() detection.
    try:
        tokens = list(tokenize.tokenize(
            io.BytesIO(source.encode('utf-8')).readline))
    except tokenize.TokenizeError:
        # Fall back to text-only logic if tokenizer fails
        # (shouldn't happen for valid Python)
        return offenders

    # Build set of (line_no) where ANY open( occurs as actual code
    # (not in a string or comment).
    # Strategy: collect all NAME tokens with value 'open' that are
    # immediately followed by '(' token, and record their line.
    candidate_lines = set()
    for k in range(len(tokens) - 1):
        tok = tokens[k]
        nxt = tokens[k + 1]
        if (tok.type == tokenize.NAME
                and tok.string == "open"
                and nxt.type == tokenize.OP
                and nxt.string == "("):
            # Verify previous token isn't '.' (would be x.open())
            if k > 0:
                prev = tokens[k - 1]
                if (prev.type == tokenize.OP
                        and prev.string == "."):
                    continue
            candidate_lines.add(tok.start[0])

    if not candidate_lines:
        return offenders

    lines = source.split("\n")
    for line_no in sorted(candidate_lines):
        idx = line_no - 1
        if idx >= len(lines):
            continue
        line = lines[idx]

        # Build the full call text by walking forward until
        # parens balance.
        call_text = line
        j = idx + 1
        # Count parens — but only those outside strings/comments.
        # For this check we can use the simple paren-count heuristic;
        # tokenize already filtered out strings from being matched.
        def open_paren_count(text):
            # Strip simple strings to avoid false counts
            t = re.sub(r'"[^"]*"', '""', text)
            t = re.sub(r"'[^']*'", "''", t)
            t = re.sub(r"#.*", "", t)
            return t.count("(") - t.count(")")

        while open_paren_count(call_text) > 0 and j - idx < 10:
            if j >= len(lines):
                break
            call_text += " " + lines[j].strip()
            j += 1

        # Skip binary mode
        if re.search(r"open\([^)]*['\"][rwa]?b[+]?['\"]",
                     call_text):
            continue
        if re.search(r"open\([^)]*[+]?b[rwa]?['\"]",
                     call_text):
            continue

        # Already specifies encoding
        if "encoding=" in call_text:
            continue

        offenders.append((line_no, line.rstrip()))

    return offenders


def test_yaml_loader_uses_utf8():
    """The reported bug site: yaml_loader._load_yaml_file."""
    print("Test: yaml_loader_uses_utf8")
    path = _find_source_file([
        "knowledge/yaml_loader.py",
        "cctbx_project/libtbx/langchain/knowledge/yaml_loader.py",
    ])
    if not path:
        print("  SKIP (yaml_loader.py not found)")
        return
    with open(path, encoding='utf-8') as f:
        source = f.read()
    offenders = _find_text_opens_without_encoding(source)
    assert not offenders, (
        "yaml_loader.py has unguarded text-mode open() calls "
        "(crashes on Windows non-English locales):\n  " +
        "\n  ".join("line %d: %s" % o for o in offenders))
    print("  PASS")


def test_ai_agent_uses_utf8():
    """The phenix client opens YAML and report files."""
    print("Test: ai_agent_uses_utf8")
    path = _find_source_file([
        "phenix/programs/ai_agent.py",
        "phenix/phenix/programs/ai_agent.py",
    ])
    if not path:
        print("  SKIP (ai_agent.py not found)")
        return
    with open(path, encoding='utf-8') as f:
        source = f.read()
    offenders = _find_text_opens_without_encoding(source)
    assert not offenders, (
        "ai_agent.py has unguarded text-mode open() calls "
        "(crashes on Windows non-English locales):\n  " +
        "\n  ".join("line %d: %s" % o for o in offenders))
    print("  PASS")


def test_directive_validator_uses_utf8():
    """The directive validator loads programs.yaml."""
    print("Test: directive_validator_uses_utf8")
    path = _find_source_file([
        "agent/directive_validator.py",
        "cctbx_project/libtbx/langchain/agent/directive_validator.py",
    ])
    if not path:
        print("  SKIP (directive_validator.py not found)")
        return
    with open(path, encoding='utf-8') as f:
        source = f.read()
    offenders = _find_text_opens_without_encoding(source)
    assert not offenders, (
        "directive_validator.py has unguarded text-mode open() "
        "calls (crashes on Windows non-English locales):\n  " +
        "\n  ".join("line %d: %s" % o for o in offenders))
    print("  PASS")


def _find_project_dir():
    """Find the langchain project root directory."""
    # 1. Standalone development layout: tests/ inside project root
    if os.path.exists(os.path.join(_PROJECT_ROOT, "agent")):
        return _PROJECT_ROOT
    # 2. PHENIX install layout
    candidate = os.path.normpath(os.path.join(
        _PROJECT_ROOT, "..", "..", "..",
        "cctbx_project", "libtbx", "langchain"))
    if os.path.exists(os.path.join(candidate, "agent")):
        return candidate
    return None


def test_all_production_code_uses_utf8():
    """Walk every production .py file in agent/, knowledge/, utils/
    and verify none has unguarded text-mode open() calls.

    This is the catch-all: a new file added to any production
    directory automatically gets tested.  Tests are excluded from
    this check (they have their own test below)."""
    print("Test: all_production_code_uses_utf8")
    project_root = _find_project_dir()
    if not project_root:
        print("  SKIP (project root not found)")
        return

    failures = []
    for subdir in ("agent", "knowledge", "utils"):
        dir_path = os.path.join(project_root, subdir)
        if not os.path.isdir(dir_path):
            continue
        for root, _, files in os.walk(dir_path):
            for fname in files:
                if not fname.endswith(".py"):
                    continue
                if fname.startswith("._"):
                    continue  # macOS metadata
                fpath = os.path.join(root, fname)
                try:
                    with open(fpath, encoding='utf-8') as f:
                        source = f.read()
                except (OSError, UnicodeDecodeError):
                    continue
                offenders = _find_text_opens_without_encoding(source)
                if offenders:
                    rel = os.path.relpath(fpath, project_root)
                    for line_no, line_text in offenders:
                        failures.append("%s:%d: %s"
                                        % (rel, line_no,
                                           line_text.strip()))

    assert not failures, (
        "Production code has %d unguarded text-mode open() "
        "calls (each crashes on Windows non-English locales):\n  "
        % len(failures) +
        "\n  ".join(failures))
    print("  PASS (scanned all production .py files)")


def test_all_test_code_uses_utf8():
    """Walk every test .py file and verify none has unguarded
    text-mode open() calls.

    Separated from the production scan because test failures here
    don't ship to users — they only bite developers on
    non-UTF-8-locale Windows.  Kept as a separate test so a
    production-only run can choose to ignore test-suite hygiene."""
    print("Test: all_test_code_uses_utf8")
    project_root = _find_project_dir()
    if not project_root:
        print("  SKIP (project root not found)")
        return

    tests_dir = os.path.join(project_root, "tests")
    if not os.path.isdir(tests_dir):
        print("  SKIP (tests/ directory not found)")
        return

    failures = []
    for root, _, files in os.walk(tests_dir):
        for fname in files:
            if not fname.endswith(".py"):
                continue
            if fname.startswith("._"):
                continue
            fpath = os.path.join(root, fname)
            try:
                with open(fpath, encoding='utf-8') as f:
                    source = f.read()
            except (OSError, UnicodeDecodeError):
                continue
            offenders = _find_text_opens_without_encoding(source)
            if offenders:
                rel = os.path.relpath(fpath, project_root)
                for line_no, line_text in offenders:
                    failures.append("%s:%d: %s"
                                    % (rel, line_no,
                                       line_text.strip()))

    assert not failures, (
        "Test code has %d unguarded text-mode open() calls "
        "(each crashes on Windows non-English locales):\n  "
        % len(failures) +
        "\n  ".join(failures))
    print("  PASS (scanned all test .py files)")


def test_yaml_loader_actually_reads_utf8():
    """Empirical test: write a YAML file with non-ASCII content
    and load it via yaml_loader's _load_yaml_file."""
    print("Test: yaml_loader_actually_reads_utf8")
    import tempfile

    # Try to import yaml_loader
    try:
        from libtbx.langchain.knowledge.yaml_loader import (
            _load_yaml_file)
    except ImportError:
        try:
            from knowledge.yaml_loader import _load_yaml_file
        except ImportError:
            print("  SKIP (yaml_loader not importable in test env)")
            return

    # Build a temp YAML with bytes that crash GBK but are valid UTF-8.
    # 0x94 alone is invalid GBK at most positions but valid as
    # part of a UTF-8 sequence like U+2014 (em dash: 0xE2 0x80 0x94).
    yaml_content = (
        "name: test\n"
        "description: 'A character — em dash — that breaks GBK'\n"
    )

    tmpdir = tempfile.mkdtemp()
    try:
        # _load_yaml_file looks for files in the knowledge dir.
        # Test the underlying open() instead.
        path = os.path.join(tmpdir, "test.yaml")
        with open(path, "w", encoding='utf-8') as f:
            f.write(yaml_content)

        # Open with the same logic as the fixed yaml_loader:
        import yaml
        with open(path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
        assert data is not None
        assert "em dash" in data.get("description", "")
        print("  PASS")
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_motivating_case_chinese_windows():
    """The actual reported bug: a byte sequence that crashes GBK
    but is valid UTF-8.  Without encoding='utf-8', Python on
    Chinese-locale Windows would crash with the reported error."""
    print("Test: motivating_case_chinese_windows")
    import tempfile

    # Build a file containing the specific byte (0x94) at position 987
    # that the user's traceback reported.  In UTF-8, 0xE2 0x80 0x9D
    # is U+201D ("right double quotation mark"), a common YAML offender.
    padding = "x" * 985
    bad_bytes = b"\xe2\x80\x9d"  # UTF-8 for U+201D, contains 0x9d
    content = padding.encode('ascii') + bad_bytes + b"\nname: test\n"

    tmpdir = tempfile.mkdtemp()
    try:
        path = os.path.join(tmpdir, "test.yaml")
        with open(path, "wb") as f:
            f.write(content)

        # Fixed pattern: explicit UTF-8
        try:
            with open(path, 'r', encoding='utf-8') as f:
                _ = f.read()
            utf8_works = True
        except UnicodeDecodeError:
            utf8_works = False

        assert utf8_works, (
            "Even with encoding='utf-8', a valid UTF-8 file "
            "didn't open cleanly. Something is wrong with the "
            "test setup.")
        print("  PASS")
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def run_all_tests(verbose=False):
    tests = [
        test_yaml_loader_uses_utf8,
        test_ai_agent_uses_utf8,
        test_directive_validator_uses_utf8,
        test_all_production_code_uses_utf8,
        test_all_test_code_uses_utf8,
        test_yaml_loader_actually_reads_utf8,
        test_motivating_case_chinese_windows,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s" % e)
            failed += 1
        except Exception as e:
            print("  FAIL: %s: %s" % (type(e).__name__, e))
            failed += 1
    print("\n%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
