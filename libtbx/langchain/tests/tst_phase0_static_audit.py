"""
Phase 0: Automated Static Audit Gate.

Checks specific patterns that cause production bugs in the AI agent.
This is NOT a style linter — it targets the exact failure modes
observed in sessions 12a and 13.

Run: python tests/tst_phase0_static_audit.py
Exit code 0 = PASS, non-zero = FAIL (blocks downstream phases).

Produces: findings/phase_0_static_audit.yaml
"""

from __future__ import absolute_import, division, print_function

import ast
import os
import re
import sys

AGENT_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'agent')
TESTS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'tests')
FINDINGS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'findings')


def get_py_files(directory):
    """Get all .py files in a directory (non-recursive).
    Excludes macOS resource forks (._*) and __pycache__."""
    if not os.path.isdir(directory):
        return []
    return sorted([
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith('.py') and not f.startswith('__')
        and not f.startswith('._')
    ])


# =====================================================================
# CHECK 1: ast.parse() every .py file
# =====================================================================

def check_parse(files):
    """Verify every .py file is syntactically valid."""
    failures = []
    for path in files:
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                ast.parse(f.read(), filename=path)
        except SyntaxError as e:
            failures.append({
                'file': os.path.basename(path),
                'line': e.lineno,
                'error': str(e),
            })
    return failures


# =====================================================================
# CHECK 2: No bare except:
# =====================================================================

def check_bare_except(files):
    """Find bare 'except:' without a named exception type.

    Bare excepts swallow KeyboardInterrupt, SystemExit, and
    generator cleanup — all of which cause silent hangs in the
    agent's long-running process.
    """
    failures = []
    pattern = re.compile(r'^\s*except\s*:\s*(#.*)?$')
    for path in files:
        try:
            with open(path) as f:
                for i, line in enumerate(f, 1):
                    if pattern.match(line):
                        failures.append({
                            'file': os.path.basename(path),
                            'line': i,
                            'content': line.strip(),
                        })
        except Exception:
            pass  # Binary or unreadable file
    return failures


# =====================================================================
# CHECK 3: Exception swallowing has comments
# =====================================================================

def check_exception_comments(files):
    """Find 'except Exception: pass' without explanatory comment.

    Swallowing exceptions silently makes debugging impossible.
    Every 'except Exception: pass' must have a comment on the
    same line or the line immediately before.
    """
    failures = []
    # Match: except Exception: pass  or  except (Exception): pass
    # Also: except (OSError, TypeError): pass
    pass_pattern = re.compile(
        r'^\s*except\s+[\w(,\s)]+:\s*$')
    for path in files:
        try:
            with open(path) as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                stripped = line.strip()
                # Look for "except SomeError:" followed by "pass" on next line
                if pass_pattern.match(line):
                    # Check if next line is just "pass"
                    if i + 1 < len(lines):
                        next_stripped = lines[i + 1].strip()
                        if next_stripped == 'pass':
                            # Need a comment on except line, pass line,
                            # or the line before except
                            has_comment = (
                                '#' in line or
                                '#' in lines[i + 1] or
                                (i > 0 and '#' in lines[i - 1])
                            )
                            if not has_comment:
                                failures.append({
                                    'file': os.path.basename(path),
                                    'line': i + 1,
                                    'content': stripped,
                                    'note': 'except...pass without '
                                            'explanatory comment',
                                })
                # Also match single-line: except Exception: pass
                if re.match(r'^\s*except\s+[\w(,\s)]+:\s*pass\s*$',
                            stripped):
                    has_comment = (
                        '#' in line or
                        (i > 0 and '#' in lines[i - 1])
                    )
                    if not has_comment:
                        failures.append({
                            'file': os.path.basename(path),
                            'line': i + 1,
                            'content': stripped,
                            'note': 'single-line except...pass '
                                    'without comment',
                        })
        except Exception:
            pass  # Binary or unreadable file
    return failures


# =====================================================================
# CHECK 4: No unguarded print() in agent/ modules
# =====================================================================

def check_unguarded_print(agent_files):
    """Find print() calls not behind an env-var gate.

    print() in agent/ writes to stdout, which in production is the
    log file. Unguarded prints pollute logs and can cause I/O
    blocking. Only print() behind PHENIX_AGENT_DIAG_* checks
    are allowed.

    Excludes:
    - print(..., file=...) — writes to StringIO, not stdout
    - Lines inside comments or strings
    """
    failures = []
    print_pattern = re.compile(r'^\s*print\s*\(')
    gate_pattern = re.compile(
        r'PHENIX_AGENT_DIAG|if\s+_diag|if\s+diag')
    # print() with file= argument is not stdout
    file_arg_pattern = re.compile(r'file\s*=')

    for path in agent_files:
        try:
            with open(path) as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                if print_pattern.match(line):
                    # Skip print(..., file=...) — not stdout
                    if file_arg_pattern.search(line):
                        continue
                    # Check if gated: look up to 5 lines back for
                    # an if _diag or env-var check
                    gated = False
                    for j in range(max(0, i - 5), i):
                        if gate_pattern.search(lines[j]):
                            gated = True
                            break
                    if not gated:
                        failures.append({
                            'file': os.path.basename(path),
                            'line': i + 1,
                            'content': line.strip()[:80],
                        })
        except Exception:
            pass  # Binary or unreadable file
    return failures


# =====================================================================
# CHECK 5: Import fallbacks for iotbx
# =====================================================================

def check_import_fallbacks(agent_files):
    """Find iotbx imports without try/except ImportError fallback.

    The agent runs in test environments without iotbx. Every iotbx
    import must either:
    - Be inside a try: block (primary import)
    - Be inside an except ImportError: block (IS the fallback)
    - Be in a function documented as "May raise" (caller handles it)
    """
    # Accepted exclusions: imports protected at caller level.
    # Each entry: (basename, line_number, reason)
    _ACCEPTED = {
        ('validation_inspector.py', 673):
            'Caller run_validation() wraps in except Exception',
        ('validation_inspector.py', 816):
            'Inside _find_diff_peaks_inner (May raise)',
        ('validation_inspector.py', 74):
            'Inside _run_validation_inner (May raise)',
    }
    failures = []
    iotbx_pattern = re.compile(
        r'^\s*(from\s+iotbx|import\s+iotbx)')

    for path in agent_files:
        try:
            with open(path) as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                if iotbx_pattern.match(line):
                    indent = len(line) - len(line.lstrip())
                    protected = False
                    # Scan upward for try: or except ImportError:
                    for j in range(i - 1, max(0, i - 15) - 1, -1):
                        prev = lines[j]
                        prev_stripped = prev.strip()
                        prev_indent = len(prev) - len(prev.lstrip())
                        if prev_stripped.startswith('try:'):
                            if prev_indent <= indent:
                                protected = True
                                break
                        if prev_stripped.startswith('except ImportError'):
                            if prev_indent <= indent:
                                # This IS a fallback import
                                protected = True
                                break
                        if prev_stripped.startswith('except'):
                            if prev_indent <= indent:
                                protected = True
                                break
                        # Check for "May raise" in nearby docstring
                        if 'May raise' in prev_stripped:
                            protected = True
                            break
                        # Stop at function def or class def
                        if (prev_stripped.startswith('def ') or
                                prev_stripped.startswith('class ')):
                            # Check if the function docstring says
                            # "May raise"
                            for k in range(j + 1,
                                           min(j + 5, len(lines))):
                                if 'May raise' in lines[k]:
                                    protected = True
                                    break
                            break
                    if not protected:
                        bname = os.path.basename(path)
                        key = (bname, i + 1)
                        if key in _ACCEPTED:
                            continue  # Accepted exclusion
                        failures.append({
                            'file': bname,
                            'line': i + 1,
                            'content': line.strip(),
                            'note': 'iotbx import without '
                                    'try/except fallback',
                        })
        except Exception:
            pass  # Binary or unreadable file
    return failures


# =====================================================================
# RUN ALL CHECKS
# =====================================================================

def run_audit():
    """Run all static audit checks. Returns (pass, findings_dict)."""
    agent_files = get_py_files(AGENT_DIR)
    test_files = get_py_files(TESTS_DIR)
    all_files = agent_files + test_files

    print("Phase 0: Static Audit")
    print("  Agent files: %d" % len(agent_files))
    print("  Test files:  %d" % len(test_files))
    print()

    findings = {
        'phase': 0,
        'name': 'static_audit',
        'checks': {},
    }

    all_passed = True

    # (check_id, description, check_fn, blocking)
    # blocking=True means failure blocks downstream phases
    # blocking=False means informational (pre-existing debt)
    checks = [
        ('parse', 'ast.parse() all .py files',
         lambda: check_parse(all_files), True),
        ('bare_except', 'No bare except:',
         lambda: check_bare_except(agent_files), True),
        ('exception_comments', 'Exception swallowing has comments',
         lambda: check_exception_comments(agent_files), False),
        ('unguarded_print', 'No unguarded print() in agent/',
         lambda: check_unguarded_print(agent_files), False),
        ('import_fallbacks', 'iotbx imports have fallbacks',
         lambda: check_import_fallbacks(agent_files), True),
    ]

    for check_id, description, check_fn, blocking in checks:
        failures = check_fn()
        if failures:
            status = 'FAIL' if blocking else 'WARN'
        else:
            status = 'PASS'
        if failures and blocking:
            all_passed = False
        findings['checks'][check_id] = {
            'description': description,
            'status': status,
            'blocking': blocking,
            'failure_count': len(failures),
            'failures': failures[:20],  # Cap at 20 to avoid noise
        }
        if not failures:
            marker = "  PASS"
        elif blocking:
            marker = "  FAIL (%d) — BLOCKING" % len(failures)
        else:
            marker = "  WARN (%d) — informational" % len(failures)
        print("  Check: %-40s %s" % (description, marker))
        if failures:
            for f in failures[:5]:
                print("    %s:%s %s" % (
                    f.get('file', '?'),
                    f.get('line', '?'),
                    f.get('content', f.get('error', ''))[:60]))
            if len(failures) > 5:
                print("    ... and %d more" % (len(failures) - 5))

    findings['status'] = 'PASS' if all_passed else 'FAIL'
    print()
    print("Phase 0 overall: %s" % findings['status'])

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    findings_path = os.path.join(FINDINGS_DIR,
                                 'phase_0_static_audit.yaml')
    # Write as simple YAML-like format (no pyyaml dependency)
    with open(findings_path, 'w') as f:
        f.write("# Phase 0: Static Audit Findings\n")
        f.write("# Generated by tst_phase0_static_audit.py\n")
        f.write("phase: 0\n")
        f.write("name: static_audit\n")
        f.write("status: %s\n" % findings['status'])
        f.write("checks:\n")
        for check_id, check_data in findings['checks'].items():
            f.write("  %s:\n" % check_id)
            f.write("    status: %s\n" % check_data['status'])
            f.write("    failure_count: %d\n"
                    % check_data['failure_count'])
            if check_data['failures']:
                f.write("    failures:\n")
                for fail in check_data['failures']:
                    f.write("      - file: %s\n"
                            % fail.get('file', '?'))
                    f.write("        line: %s\n"
                            % fail.get('line', '?'))
    print("Findings written to: %s" % findings_path)

    return all_passed, findings


if __name__ == "__main__":
    passed, _ = run_audit()
    sys.exit(0 if passed else 1)


def run_all_tests():
    """Entry point for run_all_tests.py integration."""
    passed, findings = run_audit()
    if not passed:
        raise AssertionError(
            "Phase 0 static audit FAILED: %d blocking checks failed"
            % sum(1 for c in findings.get('checks', {}).values()
                  if c.get('status') == 'FAIL'
                  and c.get('blocking')))

