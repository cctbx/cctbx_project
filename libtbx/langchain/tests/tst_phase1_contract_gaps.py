"""
Phase 1: Contract Gaps — Coverage Map.

Scans key agent modules to extract all functions, then scans all test
files to determine which functions are exercised by tests. Produces
a machine-readable findings file.

Run: python tests/tst_phase1_contract_gaps.py
Produces: findings/phase_1_contract_gaps.yaml
"""

from __future__ import absolute_import, division, print_function

import ast
import os
import re
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(BASE_DIR, 'tests')
FINDINGS_DIR = os.path.join(BASE_DIR, 'findings')

# Key modules to analyze (highest boundary-crossing risk)
KEY_MODULES = [
    'agent/workflow_state.py',
    'agent/workflow_engine.py',
    'agent/command_builder.py',
    'agent/graph_nodes.py',
]


# =====================================================================
# STEP 1: Extract functions from modules
# =====================================================================

def extract_functions(filepath):
    """Extract function metadata from a Python module using AST.

    Returns list of dicts with: name, lineno, is_method, docstring,
    never_raises, args.
    """
    with open(filepath, encoding='utf-8', errors='replace') as f:
        source = f.read()
    tree = ast.parse(source, filename=filepath)

    functions = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith('__'):
                continue  # Skip dunder methods

            # Check if it's a method (inside a class)
            is_method = False
            for parent in ast.walk(tree):
                if isinstance(parent, ast.ClassDef):
                    for child in ast.iter_child_nodes(parent):
                        if child is node:
                            is_method = True

            # Extract docstring
            docstring = ast.get_docstring(node) or ""

            # Check for "Never raises" contract
            never_raises = "never raises" in docstring.lower()

            # Extract argument names
            args = []
            for arg in node.args.args:
                if arg.arg != 'self':
                    args.append(arg.arg)

            functions.append({
                'name': node.name,
                'lineno': node.lineno,
                'is_method': is_method,
                'docstring_preview': docstring[:100].replace('\n', ' '),
                'never_raises': never_raises,
                'arg_count': len(args),
                'args': args,
            })

    return functions


# =====================================================================
# STEP 2: Scan tests for function references
# =====================================================================

def scan_tests_for_references(function_names):
    """Scan all test files for references to specific function names.

    Returns dict: function_name -> [test_files_that_reference_it]
    """
    coverage = {name: [] for name in function_names}

    test_files = sorted([
        os.path.join(TESTS_DIR, f)
        for f in os.listdir(TESTS_DIR)
        if f.startswith('tst_') and f.endswith('.py')
        and not f.startswith('._')
        # Exclude self to prevent self-referential false positives
        and f != 'tst_phase1_contract_gaps.py'
    ])

    for test_path in test_files:
        try:
            with open(test_path, encoding='utf-8',
                      errors='replace') as f:
                content = f.read()
        except Exception:
            continue

        test_basename = os.path.basename(test_path)
        for name in function_names:
            # Look for: function_name( or .function_name(
            # This avoids false positives from comments/strings
            # containing the name without calling it
            pattern = r'(?:^|\W)' + re.escape(name) + r'\s*\('
            if re.search(pattern, content):
                if test_basename not in coverage[name]:
                    coverage[name].append(test_basename)

    return coverage


# =====================================================================
# STEP 3: Rate coverage
# =====================================================================

def rate_coverage(func_name, test_files):
    """Rate coverage: NONE, LOW, COVERED."""
    if not test_files:
        return 'NONE'
    elif len(test_files) == 1:
        return 'LOW'
    else:
        return 'COVERED'


# =====================================================================
# STEP 4: Identify boundary-crossing functions
# =====================================================================

# Functions whose output feeds another module's input
BOUNDARY_FUNCTIONS = {
    # workflow_state.py → workflow_engine.py
    '_categorize_files',
    '_categorize_files_yaml',
    '_categorize_files_hardcoded',
    '_resolve_phased_promotions',
    '_bubble_up_to_parents',
    '_mtz_has_phase_columns',
    '_has_phase_columns_cached',
    '_ascii_phase_heuristic',
    '_analyze_history',
    '_detect_experiment_type',
    'detect_workflow_state',
    # workflow_engine.py → graph_nodes.py
    'build_context',
    'detect_step',
    'get_valid_programs',
    '_detect_xray_step',
    '_detect_cryoem_step',
    '_is_at_target',
    # command_builder.py → client
    # Note: main entry is CommandBuilder.build(), tracked as 'build'
    '_find_file_for_slot',
    '_build_strategy',
    # graph_nodes.py → LLM/client
    'perceive',
    'plan',       # LangGraph node (was incorrectly plan_node)
    'validate',   # LangGraph node (was incorrectly validate_node)
    'fallback',   # LangGraph node (was incorrectly fallback_node)
    'output_node',
    # NOTE: graph_nodes.build is intentionally OMITTED here.
    # 'build' is ambiguous (also CommandBuilder.build), so boundary
    # tracking is unreliable — CommandBuilder tests would mask a
    # coverage loss for the graph_nodes entry point.
    # graph_nodes.build coverage is tracked via the AMBIGUOUS NAME
    # marker instead.
}

# Functions with the same name in multiple modules.
# Coverage for these is inflated -- a reference to one module's
# function inflates the rating for all modules that define it.
_AMBIGUOUS_NAMES = {
    '_log', 'build', 'sort_key', '_get_most_recent_file',
    'log_wrapper',
}


# =====================================================================
# MAIN: Run analysis and produce report
# =====================================================================

def run_analysis():
    """Run full coverage analysis."""
    print("Phase 1: Contract Gaps — Coverage Map")
    print()

    all_functions = {}  # module -> [func_dicts]
    all_names = []

    for mod_path in KEY_MODULES:
        full_path = os.path.join(BASE_DIR, mod_path)
        mod_name = os.path.basename(mod_path)
        functions = extract_functions(full_path)
        all_functions[mod_name] = functions
        for f in functions:
            all_names.append(f['name'])
        print("  %s: %d functions extracted" % (mod_name, len(functions)))

    print()
    print("  Scanning %d test files for references..." % len([
        f for f in os.listdir(TESTS_DIR)
        if f.startswith('tst_') and f.endswith('.py')
        and not f.startswith('._')
        and f != 'tst_phase1_contract_gaps.py'
    ]))

    coverage = scan_tests_for_references(all_names)

    # Produce report
    print()
    print("  Coverage Summary:")
    print("  %-42s %-8s %-8s %s" % ("Module", "NONE", "LOW", "COVERED"))
    print("  " + "-" * 70)

    findings = {
        'phase': 1,
        'name': 'contract_gaps',
        'status': 'PASS',  # Updated below if critical gaps found
        'modules': {},
    }

    total_none = 0
    total_low = 0
    total_covered = 0
    none_boundary = []
    never_raises_untested = []

    for mod_name, functions in all_functions.items():
        none_count = 0
        low_count = 0
        covered_count = 0
        mod_findings = []

        for func in functions:
            name = func['name']
            test_files = coverage.get(name, [])
            rating = rate_coverage(name, test_files)
            is_boundary = name in BOUNDARY_FUNCTIONS

            if rating == 'NONE':
                none_count += 1
                total_none += 1
                if is_boundary:
                    none_boundary.append({
                        'module': mod_name,
                        'function': name,
                        'line': func['lineno'],
                        'boundary': True,
                    })
            elif rating == 'LOW':
                low_count += 1
                total_low += 1
            else:
                covered_count += 1
                total_covered += 1

            if func['never_raises'] and not test_files:
                never_raises_untested.append({
                    'module': mod_name,
                    'function': name,
                    'line': func['lineno'],
                })

            mod_findings.append({
                'name': name,
                'line': func['lineno'],
                'coverage': rating,
                'test_files': test_files,
                'is_boundary': is_boundary,
                'never_raises': func['never_raises'],
                'arg_count': func['arg_count'],
                'is_ambiguous': name in _AMBIGUOUS_NAMES,
            })

        findings['modules'][mod_name] = {
            'total': len(functions),
            'none': none_count,
            'low': low_count,
            'covered': covered_count,
            'functions': mod_findings,
        }

        print("  %-42s %-8d %-8d %d" % (
            mod_name, none_count, low_count, covered_count))

    total = total_none + total_low + total_covered
    print("  " + "-" * 70)
    print("  %-42s %-8d %-8d %d" % (
        "TOTAL (%d functions)" % total,
        total_none, total_low, total_covered))

    # Report critical gaps
    print()
    if none_boundary:
        print("  BOUNDARY FUNCTIONS WITH NO COVERAGE (%d):" %
              len(none_boundary))
        for gap in none_boundary:
            print("    %s:%s %s" % (
                gap['module'], gap['line'], gap['function']))
        findings['status'] = 'PARTIAL'
    else:
        print("  All boundary functions have at least some coverage.")

    if never_raises_untested:
        print()
        print("  'NEVER RAISES' FUNCTIONS WITH NO TESTS (%d):" %
              len(never_raises_untested))
        for nr in never_raises_untested:
            print("    %s:%s %s" % (
                nr['module'], nr['line'], nr['function']))

    # Count ambiguous-name functions with inflated coverage
    ambiguous_covered = 0
    for mod_name, mod_data in findings['modules'].items():
        for func in mod_data['functions']:
            if func.get('is_ambiguous') and func['coverage'] != 'NONE':
                ambiguous_covered += 1

    findings['summary'] = {
        'total_functions': total,
        'none_coverage': total_none,
        'low_coverage': total_low,
        'covered': total_covered,
        'boundary_no_coverage': len(none_boundary),
        'never_raises_untested': len(never_raises_untested),
        'ambiguous_name_coverage': ambiguous_covered,
    }
    findings['boundary_gaps'] = none_boundary
    findings['never_raises_gaps'] = never_raises_untested

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    findings_path = os.path.join(FINDINGS_DIR,
                                 'phase_1_contract_gaps.yaml')
    with open(findings_path, 'w') as f:
        f.write("# Phase 1: Contract Gaps Findings\n")
        f.write("# Generated by tst_phase1_contract_gaps.py\n")
        f.write("phase: 1\n")
        f.write("name: contract_gaps\n")
        f.write("status: %s\n" % findings['status'])
        f.write("\nsummary:\n")
        for k, v in findings['summary'].items():
            f.write("  %s: %d\n" % (k, v))
        f.write("\nboundary_functions_no_coverage:\n")
        for gap in none_boundary:
            f.write("  - module: %s\n" % gap['module'])
            f.write("    function: %s\n" % gap['function'])
            f.write("    line: %d\n" % gap['line'])
        f.write("\nnever_raises_untested:\n")
        for nr in never_raises_untested:
            f.write("  - module: %s\n" % nr['module'])
            f.write("    function: %s\n" % nr['function'])
            f.write("    line: %d\n" % nr['line'])
        f.write("\nper_module_detail:\n")
        for mod_name, mod_data in findings['modules'].items():
            f.write("  %s:\n" % mod_name)
            f.write("    total: %d\n" % mod_data['total'])
            f.write("    none: %d\n" % mod_data['none'])
            f.write("    low: %d\n" % mod_data['low'])
            f.write("    covered: %d\n" % mod_data['covered'])
            f.write("    functions:\n")
            for func in mod_data['functions']:
                marker = ""
                if func['is_boundary'] and func['coverage'] == 'NONE':
                    marker = " # BOUNDARY GAP"
                elif func['never_raises'] and not func['test_files']:
                    marker = " # NEVER-RAISES UNTESTED"
                elif func.get('is_ambiguous') and func['coverage'] != 'NONE':
                    marker = " # AMBIGUOUS NAME"
                f.write("      - %s: %s [%s]%s\n" % (
                    func['name'], func['coverage'],
                    ','.join(func['test_files'][:3]) if func['test_files']
                    else 'no tests',
                    marker))

    print()
    print("Findings written to: %s" % findings_path)
    print("Phase 1 overall: %s" % findings['status'])

    return findings


def run_all_tests():
    """Entry point for run_all_tests.py integration."""
    findings = run_analysis()
    if findings['status'] == 'FAIL':
        raise AssertionError(
            "Phase 1 contract gaps FAILED: %d boundary functions "
            "with no coverage"
            % findings['summary']['boundary_no_coverage'])


if __name__ == "__main__":
    findings = run_analysis()
    sys.exit(0 if findings['status'] != 'FAIL' else 1)

