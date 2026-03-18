"""
Phase 7: Routing Simulation.

For each tutorial with known files, simulate 3 workflow cycles
using the real production routing pipeline (detect_step +
get_valid_programs). Verifies:
  - The routing pipeline doesn't crash on real file combinations
  - The workflow doesn't stop prematurely (STOP before cycle 3)
  - The workflow doesn't get stuck in one step for 2+ cycles

Does NOT verify that specific programs are chosen — only that
the routing engine produces a viable program sequence.

Only mock: _mtz_has_phase_columns (iotbx unavailable in tests).

Run: python tests/tst_phase7_routing_simulation.py
Produces: findings/phase_7_routing_failures.yaml
"""

from __future__ import absolute_import, division, print_function

import json
import os
import re
import sys
import tempfile
import shutil

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

from agent.workflow_state import (
    _categorize_files, _analyze_history, _PHASE_COLUMN_CACHE)
from agent.workflow_engine import WorkflowEngine
import agent.workflow_state as ws

BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
FINDINGS_DIR = os.path.join(BASE_DIR, 'findings')


def load_tutorials():
    """Extract tutorial file lists from report_solve.json."""
    report_path = os.path.join(BASE_DIR, 'report_solve.json')
    with open(report_path) as f:
        solve = json.load(f)['tutorials']

    tutorials = {}
    for name, modes in sorted(solve.items()):
        mode_data = list(modes.values())[0]
        advice = mode_data.get('advice', '')

        m = re.search(
            r'\*\*Input Files Found\*\*:\s*(.+?)(?:\n|$)',
            advice)
        if not m:
            continue
        files = [f.strip().rstrip(',')
                 for f in m.group(1).split(',')]
        files = [f for f in files if '.' in f and len(f) > 2]
        if not files:
            continue

        exp_type = 'xray'
        exts = set(f.rsplit('.', 1)[1].lower()
                   for f in files if '.' in f)
        if exts & {'ccp4', 'mrc', 'map'}:
            exp_type = 'cryoem'

        tutorials[name] = {
            'files': files,
            'experiment_type': exp_type,
        }

    return tutorials


def create_tutorial_files(tmpdir, file_list):
    """Create temp files with realistic content/sizes."""
    paths = []
    for name in file_list:
        path = os.path.join(tmpdir, name)
        parent = os.path.dirname(path)
        if parent and not os.path.exists(parent):
            os.makedirs(parent)
        ext = os.path.splitext(name)[1].lower()
        with open(path, 'w') as f:
            if ext == '.pdb':
                for i in range(800):
                    f.write(
                        "ATOM  %5d  CA  ALA A%4d"
                        "       %7.3f  %7.3f  %7.3f"
                        "  1.00 20.00\n"
                        % (i+1, i+1, i*1.0, 0.0, 0.0))
            elif ext in ('.seq', '.dat', '.fasta'):
                f.write(">sequence\n")
                f.write("MVLSPADKTNVKAAWGKVGA" * 10 + "\n")
            elif ext in ('.hkl', '.sca', '.xplor'):
                if 'phases' in name.lower():
                    f.write("COLUMN_LABELS H K L PHIB FOM "
                            "FP SIGFP\n")
                else:
                    f.write("COLUMN_LABELS H K L FP SIGFP\n")
                f.write("1 0 0 120.5 3.2\n")
            elif ext == '.mtz':
                f.write("MTZ_DUMMY\n" * 100)
            elif ext in ('.ccp4', '.mrc', '.map', '.sit'):
                f.write("MAP " + "\x00" * 200 + "\n")
            elif ext == '.cv':
                f.write("COLUMN_LABELS H K L FREE\n")
            elif ext == '.cif':
                f.write("data_ligand\n_chem_comp.id LIG\n")
            elif ext == '.json':
                f.write("{}\n")
            elif ext == '.ncs_spec':
                f.write("new_operator\nrota_matrix 1 0 0\n")
            elif ext in ('.param', '.eff'):
                f.write("refinement {\n}\n")
            else:
                f.write("dummy content for %s\n" % ext)
        paths.append(path)
    return paths


def simulate_result(program, cycle):
    """Generate a realistic SUCCESS result string."""
    p = program.lower()
    if 'xtriage' in p:
        return "SUCCESS: OK"
    elif 'mtriage' in p:
        return "SUCCESS: Resolution=3.5"
    elif 'real_space' in p:
        return "SUCCESS: map_cc=0.85"
    elif 'refine' in p:
        r = 0.30 - 0.03 * cycle
        return "SUCCESS: R-work=%.4f R-free=%.4f" % (
            r - 0.04, r)
    elif 'phaser' in p:
        return "SUCCESS: TFZ=12.5"
    elif 'autosol' in p:
        return "SUCCESS: R=0.35"
    elif 'autobuild' in p:
        return "SUCCESS: R-work=0.22 R-free=0.26"
    elif 'model_vs_data' in p:
        return "SUCCESS: R-work=0.38 R-free=0.40"
    elif 'molprobity' in p:
        return "SUCCESS: clashscore=5.2"
    elif 'resolve_cryo_em' in p:
        return "SUCCESS: OK"
    elif 'dock_in_map' in p:
        return "SUCCESS: map_cc=0.80"
    elif 'map_to_model' in p:
        return "SUCCESS: R=0.30"
    return "SUCCESS: OK"


def simulate_tutorial(name, tut_data, engine):
    """Simulate 3 cycles. Returns list of cycle results."""
    tmpdir = tempfile.mkdtemp()
    results = []

    try:
        file_paths = create_tutorial_files(
            tmpdir, tut_data['files'])
        exp_type = tut_data['experiment_type']

        _PHASE_COLUMN_CACHE.clear()
        categorized = _categorize_files(file_paths)
        history = []

        for cycle in range(1, 4):
            history_info = _analyze_history(history)
            context = engine.build_context(
                files=categorized,
                history_info=history_info)
            step_info = engine.detect_step(
                exp_type, context)
            valid = engine.get_valid_programs(
                exp_type, step_info, context)

            results.append({
                'cycle': cycle,
                'step': step_info.get("step", "?"),
                'reason': step_info.get("reason", "?")[:80],
                'valid_programs': valid[:5],
                'has_stop': 'STOP' in valid,
            })

            if valid == ['STOP']:
                break

            prog = valid[0] if valid[0] != 'STOP' else (
                valid[1] if len(valid) > 1 else 'STOP')
            if prog == 'STOP':
                break

            history.append({
                'program': prog,
                'result': simulate_result(prog, cycle),
                'cycle': cycle,
            })

    finally:
        shutil.rmtree(tmpdir)

    return results


def run_all_simulations():
    """Run routing simulation for all tutorials."""
    print("Phase 7: Routing Simulation")
    print()

    tutorials = load_tutorials()
    print("  Tutorials with file data: %d" % len(tutorials))
    if len(tutorials) == 0:
        raise AssertionError(
            "Phase 7: no tutorials extracted from "
            "report_solve.json — data source may be broken")
    print()

    orig_check = ws._mtz_has_phase_columns
    ws._mtz_has_phase_columns = lambda f: False

    engine = WorkflowEngine()

    # Tutorials expected to get stuck due to test-env limitations.
    # 7 cryo-EM: map files fail _is_valid_file (invalid headers)
    # 1 cryo-EM: map only, no model (groel_map_symmetry)
    # 1 X-ray: single .sca file, no model/sequence (p9-xtriage)
    _EXPECTED_STUCK = {
        'apo-ferritin',  # cryo-EM: invalid map headers in test
        'bgal_denmod', 'curli_alphafold',
        'emd_6123_map_to_model', 'groel_map_symmetry',
        'ion_channel_denmod', 'model-building-scripting',
        'p9-xtriage', 'real-space-refine-5ljv',
        'real-space-refine-6crz',
        'rotavirus-autosharpen',       # cryo-EM: invalid map
        'rotavirus-model-building',    # cryo-EM: invalid map
    }

    passed = 0
    warnings = 0
    expected_stuck = 0
    unexpected_stuck = 0
    issues = []

    try:
        for name, tut_data in sorted(tutorials.items()):
            try:
                results = simulate_tutorial(
                    name, tut_data, engine)

                cycles_run = len(results)
                last_step = results[-1]['step']
                last_valid = results[-1]['valid_programs']

                premature_stop = (
                    cycles_run < 3 and
                    last_valid == ['STOP'] and
                    last_step not in ('complete',))

                steps = [r['step'] for r in results]
                stuck_in_step = (
                    len(set(steps)) == 1 and
                    cycles_run >= 2 and
                    steps[0] != 'complete')

                if premature_stop:
                    is_expected = name in _EXPECTED_STUCK
                    status = "STOP@%d" % cycles_run
                    if is_expected:
                        status += "(exp)"
                        expected_stuck += 1
                    else:
                        unexpected_stuck += 1
                    issues.append({
                        'tutorial': name,
                        'issue': 'premature_stop',
                        'cycle': cycles_run,
                        'step': last_step,
                        'reason': results[-1]['reason'],
                        'valid': str(last_valid),
                        'expected': is_expected,
                    })
                elif stuck_in_step:
                    status = "STUCK:%s" % steps[0]
                    warnings += 1
                    issues.append({
                        'tutorial': name,
                        'issue': 'stuck_in_step',
                        'step': steps[0],
                        'cycles': str(cycles_run),
                        'expected': name in _EXPECTED_STUCK,
                    })
                else:
                    status = "OK"
                    passed += 1

                cycle_summary = " -> ".join(
                    "%s[%s]" % (
                        r['step'][:12],
                        r['valid_programs'][0][:20]
                        if r['valid_programs'] else '?')
                    for r in results)
                print("  %-35s %-10s %s" % (
                    name, status, cycle_summary))

            except Exception as e:
                print("  %-35s ERROR: %s" % (
                    name, str(e)[:60]))
                issues.append({
                    'tutorial': name,
                    'issue': 'exception',
                    'error': str(e)[:200],
                    'expected': False,
                })
                unexpected_stuck += 1

    finally:
        ws._mtz_has_phase_columns = orig_check
        _PHASE_COLUMN_CACHE.clear()

    total = passed + warnings + expected_stuck + unexpected_stuck
    print()
    print("  Results: %d OK, %d warnings, %d expected stuck, "
          "%d unexpected stuck/error out of %d"
          % (passed, warnings, expected_stuck,
             unexpected_stuck, total))

    if unexpected_stuck > 0:
        status = 'FAIL'
    elif expected_stuck > 0 or warnings > 0:
        status = 'PARTIAL'
    else:
        status = 'PASS'

    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_7_routing_failures.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 7: Routing Simulation Findings\n")
        f.write("phase: 7\n")
        f.write("name: routing_simulation\n")
        f.write("status: %s\n" % status)
        f.write("tutorials_tested: %d\n" % total)
        f.write("ok: %d\n" % passed)
        f.write("warnings: %d\n" % warnings)
        f.write("expected_stuck: %d\n" % expected_stuck)
        f.write("unexpected_stuck: %d\n" % unexpected_stuck)
        f.write("q3_tutorials: NOT_TESTED (no file data)\n")
        if issues:
            f.write("\nissues:\n")
            for iss in issues:
                f.write("  - tutorial: %s\n"
                        % iss['tutorial'])
                f.write("    type: %s\n" % iss['issue'])
                if iss.get('expected'):
                    f.write("    assessment: EXPECTED "
                            "(test-env limitation)\n")
                for k, v in iss.items():
                    if k not in ('tutorial', 'issue',
                                 'expected'):
                        f.write("    %s: %s\n" % (k, v))

    print("  Findings: %s" % path)
    print("  Phase 7 overall: %s" % status)
    return {'status': status, 'passed': passed,
            'unexpected_stuck': unexpected_stuck,
            'issues': issues}


def run_all_tests():
    """Entry point for run_all_tests.py integration."""
    result = run_all_simulations()
    if result['status'] == 'FAIL':
        unexpected = [i for i in result['issues']
                      if not i.get('expected')]
        raise AssertionError(
            "Phase 7 routing simulation FAILED: %d unexpected "
            "stuck tutorials: %s"
            % (result['unexpected_stuck'],
               [i['tutorial'] for i in unexpected]))


if __name__ == "__main__":
    result = run_all_simulations()
    sys.exit(0 if result['status'] != 'FAIL' else 1)
