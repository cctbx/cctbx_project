"""
Phase 8: Command Building Simulation.

For each tutorial x expected program, construct a CommandContext
and call CommandBuilder.build(). Verifies:
  - build() returns a non-None command (unless expect_none)
  - the command string starts with the program name
  - build() doesn't crash on real file combinations

Does NOT verify specific file selection or parameter correctness.

Bug class: unbuildable commands, wrong program prefix, crashes.

Run: python tests/tst_phase8_command_building.py
Produces: findings/phase_8_command_failures.yaml
"""

from __future__ import absolute_import, division, print_function

import os
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
    _categorize_files, _PHASE_COLUMN_CACHE)
from agent.command_builder import CommandBuilder, CommandContext
import agent.workflow_state as ws

BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
FINDINGS_DIR = os.path.join(BASE_DIR, 'findings')


# =====================================================================
# TUTORIAL DEFINITIONS
# =====================================================================
# Each tutorial: files, experiment_type, and (program, cycle, extras)
# tuples to test. extras may contain best_files, resolution, etc.
# needed to make the command buildable.

TUTORIALS = {
    'hipip-refine': {
        'type': 'xray',
        'files': ['hipip.pdb', 'hipip.mtz'],
        'commands': [
            ('phenix.xtriage', 1, {}),
            ('phenix.model_vs_data', 2, {}),
            ('phenix.refine', 3, {}),
        ],
    },
    'rab3a-refine': {
        'type': 'xray',
        'files': ['rab3a.pdb', 'rab3a_scale.hkl',
                  'rab3a_phases.hkl'],
        'commands': [
            ('phenix.xtriage', 1, {}),
            ('phenix.refine', 3, {}),
        ],
    },
    'nsf-d2-refine': {
        'type': 'xray',
        'files': ['nsf-d2_start.pdb', 'nsf-d2_phases.hkl'],
        'commands': [
            ('phenix.xtriage', 1, {}),
            ('phenix.refine', 3, {}),
        ],
    },
    'beta-blip-phaser': {
        'type': 'xray',
        'files': ['beta.pdb', 'beta_blip.mtz'],
        'commands': [
            ('phenix.xtriage', 1, {}),
            # Phaser with model in 'model' (not search_model)
            # — tests the fallback_categories fix (B2)
            ('phenix.phaser', 2, {}),
        ],
    },
    'gene-5-mad': {
        'type': 'xray',
        'files': ['sequence.dat', 'high.sca'],
        'commands': [
            ('phenix.xtriage', 1, {}),
        ],
    },
    'rsr-5ljv': {
        'type': 'cryoem',
        'files': ['model.pdb', 'map.ccp4'],
        'commands': [
            ('phenix.mtriage', 1, {}),
            # RSR needs resolution (set by mtriage)
            ('phenix.real_space_refine', 2,
             {'resolution': 3.5}),
        ],
        'real_map': True,  # Needs valid CCP4 header
    },
    'apoferritin-build': {
        'type': 'cryoem',
        'files': ['apoferritin.seq', 'apoferritin.ccp4',
                  'apoferritin.pdb'],
        'commands': [
            ('phenix.mtriage', 1, {}),
        ],
        'real_map': True,
    },
}

# Edge case tests (from plan: empty categories, long filenames,
# special characters, missing-on-disk)
EDGE_CASES = {
    'empty_data_mtz': {
        'type': 'xray',
        'files': ['model.pdb'],  # No MTZ at all
        'commands': [
            # xtriage needs data_mtz — should return None
            ('phenix.xtriage', 1, {'expect_none': True}),
        ],
    },
    'long_filename': {
        'type': 'xray',
        'files': [
            'this_is_a_very_long_protein_model_name.pdb',
            'this_is_a_very_long_data_filename.mtz',
        ],
        'commands': [
            ('phenix.refine', 2, {}),
        ],
    },
}


# =====================================================================
# HELPERS
# =====================================================================

def create_pdb(path, n_atoms=800):
    """Create a realistic PDB file."""
    with open(path, 'w') as f:
        for i in range(n_atoms):
            f.write(
                "ATOM  %5d  CA  ALA A%4d"
                "       %7.3f  %7.3f  %7.3f"
                "  1.00 20.00\n"
                % (i+1, i+1, i*1.0, 0.0, 0.0))


def create_ccp4_map(path):
    """Create a minimal valid CCP4/MRC map file."""
    with open(path, 'wb') as f:
        header = bytearray(1024)
        header[0:4] = (10).to_bytes(4, 'little')  # NX
        header[4:8] = (10).to_bytes(4, 'little')  # NY
        header[8:12] = (10).to_bytes(4, 'little')  # NZ
        header[12:16] = (2).to_bytes(4, 'little')  # MODE=2
        header[208:212] = b'MAP '  # Magic at byte 208
        header[212:216] = b'\x44\x41\x00\x00'
        f.write(header)
        f.write(b'\x00' * (10*10*10*4))


def create_tutorial_files(tmpdir, file_list, real_map=False):
    """Create temp files for a tutorial."""
    paths = []
    for name in file_list:
        path = os.path.join(tmpdir, name)
        ext = os.path.splitext(name)[1].lower()
        if ext == '.pdb':
            create_pdb(path)
        elif ext in ('.ccp4', '.mrc', '.map') and real_map:
            create_ccp4_map(path)
        elif ext in ('.ccp4', '.mrc', '.map'):
            with open(path, 'w') as f:
                f.write("MAP " + "\x00" * 200 + "\n")
        elif ext in ('.seq', '.dat', '.fasta'):
            with open(path, 'w') as f:
                f.write(">sequence\n")
                f.write("MVLSPADKTNVKAAWGKVGA" * 10 + "\n")
        elif ext in ('.hkl', '.sca'):
            with open(path, 'w') as f:
                if 'phases' in name.lower():
                    f.write("COLUMN_LABELS H K L PHIB FOM "
                            "FP SIGFP\n")
                else:
                    f.write("COLUMN_LABELS H K L FP SIGFP\n")
                f.write("1 0 0 120.5 3.2\n")
        elif ext == '.mtz':
            with open(path, 'w') as f:
                f.write("MTZ_DUMMY\n" * 100)
        elif ext == '.cif':
            with open(path, 'w') as f:
                f.write("data_ligand\n_chem_comp.id LIG\n")
        else:
            with open(path, 'w') as f:
                f.write("dummy\n")
        paths.append(path)
    return paths


# =====================================================================
# RUN SIMULATIONS
# =====================================================================

def run_all_simulations():
    """Build commands for all tutorials x programs."""
    print("Phase 8: Command Building Simulation")
    print()

    orig_check = ws._mtz_has_phase_columns
    ws._mtz_has_phase_columns = lambda f: False

    builder = CommandBuilder()

    passed = 0
    failed = 0
    issues = []
    total_commands = 0

    all_tutorials = dict(TUTORIALS)
    all_tutorials.update(EDGE_CASES)

    try:
        for tut_name, tut_def in sorted(all_tutorials.items()):
            tmpdir = tempfile.mkdtemp()
            try:
                real_map = tut_def.get('real_map', False)
                file_paths = create_tutorial_files(
                    tmpdir, tut_def['files'], real_map)

                _PHASE_COLUMN_CACHE.clear()
                categorized = _categorize_files(file_paths)

                # Build best_files from categorized output,
                # using the same keys that CommandBuilder
                # looks up (best_category from programs.yaml)
                best = {}
                for cat in ('model', 'refined'):
                    for f in categorized.get(cat, []):
                        if isinstance(f, str):
                            best.setdefault('model', f)
                            break
                for cat in ('data_mtz', 'original_data_mtz'):
                    for f in categorized.get(cat, []):
                        if isinstance(f, str):
                            best.setdefault('data_mtz', f)
                            break
                for cat in ('full_map', 'optimized_full_map'):
                    for f in categorized.get(cat, []):
                        if isinstance(f, str):
                            best.setdefault('full_map', f)
                            break

                for prog, cycle, extras in tut_def['commands']:
                    total_commands += 1
                    expect_none = extras.get(
                        'expect_none', False)

                    try:
                        ctx = CommandContext(
                            cycle_number=cycle,
                            experiment_type=tut_def['type'],
                            categorized_files=categorized,
                            best_files=best,
                            resolution=extras.get(
                                'resolution'),
                            files_local=True,
                        )

                        cmd = builder.build(
                            prog, file_paths, ctx)

                        if expect_none:
                            if cmd is None:
                                passed += 1
                                print(
                                    "  %-20s %-28s "
                                    "PASS (expected None)"
                                    % (tut_name, prog))
                            else:
                                failed += 1
                                issues.append({
                                    'tutorial': tut_name,
                                    'program': prog,
                                    'issue':
                                        'expected_none_got_cmd',
                                    'command': cmd[:60],
                                })
                                print(
                                    "  %-20s %-28s "
                                    "FAIL (expected None, "
                                    "got cmd)"
                                    % (tut_name, prog))
                        elif cmd:
                            # Verify command starts with
                            # program
                            if cmd.startswith(prog):
                                passed += 1
                                print(
                                    "  %-20s %-28s "
                                    "PASS → %s"
                                    % (tut_name, prog,
                                       cmd[:50]))
                            else:
                                failed += 1
                                issues.append({
                                    'tutorial': tut_name,
                                    'program': prog,
                                    'issue': 'wrong_prefix',
                                    'command': cmd[:60],
                                })
                                print(
                                    "  %-20s %-28s "
                                    "FAIL (wrong prefix)"
                                    % (tut_name, prog))
                        else:
                            failed += 1
                            issues.append({
                                'tutorial': tut_name,
                                'program': prog,
                                'issue':
                                    'build_returned_none',
                            })
                            print(
                                "  %-20s %-28s "
                                "FAIL (None)"
                                % (tut_name, prog))

                    except Exception as e:
                        failed += 1
                        issues.append({
                            'tutorial': tut_name,
                            'program': prog,
                            'issue': 'exception',
                            'error': str(e)[:200],
                        })
                        print(
                            "  %-20s %-28s "
                            "ERROR: %s"
                            % (tut_name, prog,
                               str(e)[:50]))

            except Exception as e:
                # Categorization or setup failed for this
                # tutorial — count all its commands as failed
                n_cmds = len(tut_def.get('commands', []))
                failed += n_cmds
                total_commands += n_cmds
                issues.append({
                    'tutorial': tut_name,
                    'program': '(setup)',
                    'issue': 'tutorial_exception',
                    'error': str(e)[:200],
                })
                print("  %-20s %-28s ERROR: %s"
                      % (tut_name, '(all commands)',
                         str(e)[:50]))
            finally:
                shutil.rmtree(tmpdir)

    finally:
        ws._mtz_has_phase_columns = orig_check
        _PHASE_COLUMN_CACHE.clear()

    print()
    print("  Results: %d passed, %d failed out of %d"
          % (passed, failed, total_commands))

    status = 'PASS' if failed == 0 else 'FAIL'

    # Write findings
    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_8_command_failures.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 8: Command Building Findings\n")
        f.write("phase: 8\n")
        f.write("name: command_building\n")
        f.write("status: %s\n" % status)
        f.write("total_commands: %d\n" % total_commands)
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        if issues:
            f.write("\nissues:\n")
            for iss in issues:
                f.write("  - tutorial: %s\n"
                        % iss['tutorial'])
                f.write("    program: %s\n"
                        % iss['program'])
                f.write("    issue: %s\n"
                        % iss['issue'])
                if 'command' in iss:
                    f.write("    command: \"%s\"\n"
                            % iss['command'])

    print("  Findings: %s" % path)
    print("  Phase 8 overall: %s" % status)
    return {'status': status, 'passed': passed,
            'failed': failed, 'issues': issues}


def run_all_tests():
    """Entry point for run_all_tests.py integration."""
    result = run_all_simulations()
    if result['status'] == 'FAIL':
        raise AssertionError(
            "Phase 8 command building FAILED: %d/%d commands "
            "failed" % (result['failed'],
                        result['passed'] + result['failed']))


if __name__ == "__main__":
    result = run_all_simulations()
    sys.exit(0 if result['status'] == 'PASS' else 1)
