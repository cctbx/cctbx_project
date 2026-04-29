#!/usr/bin/env python
"""
Automatic documentation generator for safety checks and validations.

This script scans the codebase and generates a comprehensive table of:
1. Sanity checks (pre-execution validation)
2. Post-LLM processing (directive corrections)
3. File validation (input/output verification)
4. Workflow validation (state machine checks)
5. Command validation (security and correctness)

Run with: python agent/generate_safety_docs.py > docs/SAFETY_CHECKS.md
"""

import os
import re
from collections import defaultdict

# Categories of safety checks
CATEGORIES = {
    "sanity_check": "Sanity Checks (Pre-Execution)",
    "directive_validation": "Directive Validation (Post-LLM)",
    "file_validation": "File Validation",
    "workflow_validation": "Workflow State Validation",
    "command_validation": "Command Building Validation",
    "input_validation": "Input Validation",
    "post_processing": "Post-Processing Corrections",
}

# Files to scan
SCAN_DIRS = ["agent", "phenix_ai", "programs", "knowledge"]


def find_safety_checks(base_dir):
    """
    Scan code for safety check functions and classes.

    Returns:
        dict: category -> list of {name, file, line, description, conditions}
    """
    checks = defaultdict(list)

    # Patterns to look for
    patterns = [
        # Function definitions
        (r'def (_?(?:check|validate|fix|verify|sanity|abort|red_flag)\w*)', 'function'),
        (r'def (_?(?:is_valid|has_error|is_safe|is_intermediate)\w*)', 'function'),
        (r'def (_fix_\w+)', 'function'),  # All _fix_ functions are post-processing

        # Class definitions
        (r'class ((?:Sanity|Validation|Check)\w+)', 'class'),

        # Specific check implementations
        (r'# ((?:SAFETY|POST-PROCESS|VALIDATION|CHECK|FIX):.*)', 'comment'),

        # Sanity issues
        (r'SanityIssue\s*\(\s*severity\s*=\s*["\'](\w+)["\'].*?code\s*=\s*["\']([\w_]+)["\']', 'sanity_issue'),
    ]

    for scan_dir in SCAN_DIRS:
        dir_path = os.path.join(base_dir, scan_dir)
        if not os.path.isdir(dir_path):
            continue

        for root, dirs, files in os.walk(dir_path):
            # Skip __pycache__ and test directories
            dirs[:] = [d for d in dirs if d not in ('__pycache__', 'tests', 'test')]

            for filename in files:
                if not filename.endswith('.py'):
                    continue

                filepath = os.path.join(root, filename)
                relative_path = os.path.relpath(filepath, base_dir)

                try:
                    with open(filepath, 'r') as f:
                        content = f.read()
                        lines = content.split('\n')
                except Exception:
                    continue

                # Search for patterns
                for pattern, pattern_type in patterns:
                    for match in re.finditer(pattern, content, re.MULTILINE | re.DOTALL):
                        # Find line number
                        pos = match.start()
                        line_num = content[:pos].count('\n') + 1

                        # Extract info based on pattern type
                        if pattern_type == 'function':
                            name = match.group(1)
                            desc = _extract_docstring(lines, line_num - 1)
                            category = _categorize_check(name, filepath, desc)

                            checks[category].append({
                                'name': name,
                                'file': relative_path,
                                'line': line_num,
                                'description': desc or "(No docstring)",
                                'type': 'function'
                            })

                        elif pattern_type == 'class':
                            name = match.group(1)
                            desc = _extract_docstring(lines, line_num - 1)

                            checks['sanity_check'].append({
                                'name': name,
                                'file': relative_path,
                                'line': line_num,
                                'description': desc or "(No docstring)",
                                'type': 'class'
                            })

                        elif pattern_type == 'sanity_issue':
                            severity = match.group(1)
                            code = match.group(2)

                            checks['sanity_check'].append({
                                'name': code,
                                'file': relative_path,
                                'line': line_num,
                                'description': f"Severity: {severity}",
                                'type': 'sanity_issue'
                            })

    return dict(checks)


def _extract_docstring(lines, start_line):
    """Extract docstring from function/class definition."""
    # Look for docstring in the lines following the definition
    for i in range(start_line + 1, min(start_line + 10, len(lines))):
        line = lines[i].strip()
        if line.startswith('"""') or line.startswith("'''"):
            # Single-line docstring
            if line.count('"""') >= 2 or line.count("'''") >= 2:
                return line.strip('"\' ')
            # Multi-line docstring - get first line
            quote = '"""' if '"""' in line else "'''"
            doc = line[3:]
            for j in range(i + 1, min(i + 20, len(lines))):
                if quote in lines[j]:
                    break
                doc += " " + lines[j].strip()
            return doc.strip()
        elif line and not line.startswith('#'):
            # No docstring
            break
    return None


def _categorize_check(name, filepath, description):
    """Categorize a check based on its name and location."""
    name_lower = name.lower()
    filepath_lower = filepath.lower()
    desc_lower = (description or "").lower()

    # Sanity checks
    if 'sanity' in name_lower or 'sanity' in filepath_lower:
        return 'sanity_check'
    if 'red_flag' in name_lower or 'abort' in name_lower:
        return 'sanity_check'

    # Directive validation
    if 'directive' in filepath_lower:
        if 'fix' in name_lower:
            return 'post_processing'
        return 'directive_validation'

    # Workflow validation
    if 'workflow' in filepath_lower or 'workflow' in name_lower:
        return 'workflow_validation'

    # Command validation
    if 'command' in filepath_lower or 'build' in filepath_lower:
        return 'command_validation'

    # File validation
    if 'file' in name_lower or 'path' in name_lower:
        return 'file_validation'

    # Input validation
    if 'input' in name_lower or 'param' in name_lower:
        return 'input_validation'

    # Post-processing
    if 'fix' in name_lower or 'correct' in name_lower:
        return 'post_processing'

    # Default to input validation
    return 'input_validation'


def generate_markdown(checks):
    """Generate markdown documentation from checks."""
    lines = [
        "# Safety Checks and Validations",
        "",
        "This document is auto-generated by `generate_safety_docs.py`.",
        "",
        "## Overview",
        "",
        "The PHENIX AI Agent includes multiple layers of safety checks:",
        "",
    ]

    # Summary table
    lines.append("| Category | Count |")
    lines.append("|----------|-------|")
    for cat_key, cat_name in CATEGORIES.items():
        count = len(checks.get(cat_key, []))
        if count > 0:
            lines.append(f"| {cat_name} | {count} |")
    lines.append("")

    # Detailed sections
    for cat_key, cat_name in CATEGORIES.items():
        cat_checks = checks.get(cat_key, [])
        if not cat_checks:
            continue

        lines.append(f"## {cat_name}")
        lines.append("")

        # Sort by file then name
        cat_checks.sort(key=lambda x: (x['file'], x['name']))

        # Group by file
        by_file = defaultdict(list)
        for check in cat_checks:
            by_file[check['file']].append(check)

        for filepath, file_checks in sorted(by_file.items()):
            lines.append(f"### `{filepath}`")
            lines.append("")
            lines.append("| Name | Line | Type | Description |")
            lines.append("|------|------|------|-------------|")

            for check in file_checks:
                desc = check['description'][:80] + "..." if len(check['description']) > 80 else check['description']
                desc = desc.replace("|", "\\|").replace("\n", " ")
                lines.append(f"| `{check['name']}` | {check['line']} | {check['type']} | {desc} |")

            lines.append("")

    return "\n".join(lines)


def generate_sanity_checks_table(base_dir):
    """
    Generate a specific table of all SanityChecker checks from the YAML or code.
    """
    sanity_file = os.path.join(base_dir, "agent", "sanity_checker.py")
    if not os.path.exists(sanity_file):
        return ""

    with open(sanity_file, 'r') as f:
        content = f.read()

    # Find all SanityIssue creations
    pattern = r'SanityIssue\s*\(\s*severity\s*=\s*["\'](\w+)["\'],\s*code\s*=\s*["\']([\w_]+)["\'],\s*message\s*=\s*["\']([^"\']+)["\']'

    issues = []
    for match in re.finditer(pattern, content, re.MULTILINE | re.DOTALL):
        issues.append({
            'severity': match.group(1),
            'code': match.group(2),
            'message': match.group(3)
        })

    if not issues:
        return ""

    lines = [
        "## Sanity Check Issues",
        "",
        "These are the specific issues that can be detected by the SanityChecker:",
        "",
        "| Code | Severity | Message |",
        "|------|----------|---------|",
    ]

    for issue in sorted(issues, key=lambda x: (x['severity'], x['code'])):
        msg = issue['message'][:60] + "..." if len(issue['message']) > 60 else issue['message']
        lines.append(f"| `{issue['code']}` | {issue['severity']} | {msg} |")

    lines.append("")
    return "\n".join(lines)


def main():
    # Find base directory (script is in agent/, project root is one level up)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)  # Up from agent/ to project root

    # Find all checks
    checks = find_safety_checks(base_dir)

    # Generate markdown
    md = generate_markdown(checks)

    # Add sanity checks table
    sanity_table = generate_sanity_checks_table(base_dir)
    if sanity_table:
        md = md.replace("## Sanity Checks (Pre-Execution)",
                        sanity_table + "\n## Sanity Check Functions")

    print(md)


if __name__ == "__main__":
    main()
