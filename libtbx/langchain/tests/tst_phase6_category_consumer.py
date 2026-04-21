"""
Phase 6: Category-Consumer Alignment.

For every program × tutorial combination, verify that required
input slots can find matching files in the categorized output.

Bug class found by this approach: Fix 5 (refine doesn't check
phased_data_mtz).

Run: python tests/tst_phase6_category_consumer.py
Produces: findings/phase_6_missing_inputs.yaml
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

import yaml
from agent.workflow_state import (
    _categorize_files, _PHASE_COLUMN_CACHE)
import agent.workflow_state as ws

FINDINGS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'findings')
BASE_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))


# =====================================================================
# TEST DATA: tutorials with their expected first programs
# =====================================================================

TUTORIALS = {
    'rab3a-refine': {
        'files': ['rab3a.pdb', 'rab3a_scale.hkl',
                  'rab3a_phases.hkl'],
        'expected_programs': ['phenix.xtriage', 'phenix.refine'],
    },
    'nsf-d2-refine': {
        'files': ['nsf-d2_start.pdb', 'nsf-d2_reference.pdb',
                  'nsf-d2_phases.hkl', 'nsf-d2.cv'],
        'expected_programs': ['phenix.xtriage', 'phenix.refine',
                              'phenix.autobuild'],
    },
    'hipip-refine': {
        'files': ['hipip.mtz', 'hipip.pdb', '1IUA.pdb'],
        'expected_programs': ['phenix.xtriage', 'phenix.refine'],
    },
    'lysozyme-MRSAD': {
        'files': ['1fkq_prot.pdb', 'hewl.seq',
                  'lyso2001_scala1.mtz'],
        'expected_programs': ['phenix.xtriage', 'phenix.phaser'],
    },
    'gene-5-mad': {
        'files': ['sequence.dat', 'high.sca', 'peak.sca',
                  'infl.sca'],
        'expected_programs': ['phenix.xtriage', 'phenix.autosol'],
    },
    'real-space-refine-5ljv': {
        'files': ['model.pdb', 'map.ccp4'],
        'expected_programs': ['phenix.mtriage',
                              'phenix.real_space_refine'],
    },
    'apoferritin_model_building': {
        'files': ['apoferritin.seq',
                  'apoferritin_denmod_box.ccp4',
                  'apoferritin_chainA_rsr.pdb'],
        'expected_programs': ['phenix.mtriage'],
    },
    'beta-blip': {
        'files': ['beta.pdb', 'blip.pdb', 'beta.seq',
                  'blip.seq', 'beta_blip_P3221.mtz'],
        'expected_programs': ['phenix.xtriage', 'phenix.phaser'],
    },
}


def create_tutorial_files(tmpdir, file_list):
    """Create temp files with realistic content."""
    paths = []
    for name in file_list:
        path = os.path.join(tmpdir, name)
        ext = os.path.splitext(name)[1].lower()
        with open(path, 'w') as f:
            if ext == '.pdb':
                for i in range(800):
                    f.write(
                        "ATOM  %5d  CA  ALA A%4d"
                        "       %7.3f  %7.3f  %7.3f"
                        "  1.00 20.00\n"
                        % (i+1, i+1, i*1.0, 0.0, 0.0))
            elif ext in ('.seq', '.dat'):
                f.write(">sequence\n")
                f.write("MVLSPADKTNVKAAWGKVGA" * 10 + "\n")
            elif ext in ('.hkl', '.sca'):
                if 'phases' in name.lower():
                    f.write("COLUMN_LABELS H K L PHIB FOM "
                            "FP SIGFP\n")
                else:
                    f.write("COLUMN_LABELS H K L FP SIGFP\n")
                f.write("1 0 0 120.5 3.2\n")
            elif ext == '.mtz':
                f.write("MTZ_DUMMY\n" * 100)
            elif ext in ('.ccp4', '.mrc', '.map'):
                f.write("MAP " + "\x00" * 100 + "\n")
            elif ext == '.cv':
                f.write("COLUMN_LABELS H K L FREE\n")
            else:
                f.write("dummy\n")
        paths.append(path)
    return paths


# =====================================================================
# CHECK 1: All referenced categories exist
# =====================================================================

def test_all_categories_exist():
    """Every category in input_priorities exists in categorizer."""
    print("Test: all_categories_exist")

    with open(os.path.join(BASE_DIR,
                           'knowledge', 'programs.yaml')) as f:
        progs = yaml.safe_load(f)

    # Get a representative categorizer output
    tmpdir = tempfile.mkdtemp()
    try:
        # Create files that cover all extensions
        test_files = []
        for name in ['test.mtz', 'test.pdb', 'test.seq',
                     'test.ccp4', 'test.hkl', 'test.sca']:
            path = os.path.join(tmpdir, name)
            with open(path, 'w') as f:
                if name.endswith('.pdb'):
                    for i in range(800):
                        f.write("ATOM  %5d  CA  ALA A%4d"
                                "       0.0   0.0   0.0"
                                "  1.00 20.00\n" % (i+1, i+1))
                else:
                    f.write("dummy\n")
            test_files.append(path)

        _PHASE_COLUMN_CACHE.clear()
        result = _categorize_files(test_files)
        produced = set(result.keys())

        # Check all input_priorities categories
        referenced = set()
        for prog_name, prog in progs.items():
            if not isinstance(prog, dict):
                continue
            ip = prog.get('input_priorities', {})
            for slot, slot_def in ip.items():
                for cat in slot_def.get('categories', []):
                    referenced.add(cat)
                for cat in slot_def.get(
                        'fallback_categories', []):
                    referenced.add(cat)

        missing = referenced - produced
        if missing:
            print("  Categories referenced but not produced: %s"
                  % sorted(missing))
        assert len(missing) == 0, (
            "Missing categories: %s" % sorted(missing))
    finally:
        shutil.rmtree(tmpdir)
    print("  PASSED (%d categories all exist)" % len(referenced))


# =====================================================================
# CHECK 2: Per-tutorial, expected programs can find required inputs
# =====================================================================

def test_tutorial_programs_have_inputs():
    """For each tutorial × expected program, verify required
    input slots have files in the categorized output.

    Checks BOTH primary categories AND fallback_categories,
    matching what the command builder actually does."""
    print("Test: tutorial_programs_have_inputs")

    with open(os.path.join(BASE_DIR,
                           'knowledge', 'programs.yaml')) as f:
        progs = yaml.safe_load(f)

    orig_check = ws._mtz_has_phase_columns
    ws._mtz_has_phase_columns = lambda f: False
    _PHASE_COLUMN_CACHE.clear()

    # Test-env limitations that produce expected missing inputs
    _EXPECTED_MISSING = {
        # CCP4 maps have invalid headers in test env
        ('real-space-refine-5ljv', 'phenix.real_space_refine',
         'map'),
    }

    missing_inputs = []
    total_checks = 0
    errors = 0

    try:
        for tut_name, tut_def in sorted(TUTORIALS.items()):
            tmpdir = tempfile.mkdtemp()
            try:
                files = create_tutorial_files(
                    tmpdir, tut_def['files'])
                _PHASE_COLUMN_CACHE.clear()
                categorized = _categorize_files(files)

                for prog_name in tut_def['expected_programs']:
                    prog = progs.get(prog_name, {})
                    ip = prog.get('input_priorities', {})
                    required = prog.get('inputs', {}).get(
                        'required', {})

                    for slot_name, slot_def in ip.items():
                        # Only check slots that are required
                        is_required = slot_name in required
                        if not is_required:
                            continue

                        # Check primary categories
                        categories = slot_def.get(
                            'categories', [])
                        # Also check fallback_categories
                        # (command builder uses these when
                        # primary categories are empty)
                        fallback = slot_def.get(
                            'fallback_categories', [])
                        all_cats = categories + fallback

                        total_checks += 1

                        # Check if any category has files
                        found = False
                        for cat in all_cats:
                            cat_files = categorized.get(cat, [])
                            if isinstance(cat_files, list):
                                if cat_files:
                                    found = True
                                    break

                        if not found:
                            missing_inputs.append({
                                'tutorial': tut_name,
                                'program': prog_name,
                                'slot': slot_name,
                                'categories': categories,
                                'fallback': fallback,
                                'expected': (tut_name, prog_name,
                                             slot_name)
                                            in _EXPECTED_MISSING,
                            })
            except Exception as e:
                print("    ERROR in %s: %s" % (
                    tut_name, str(e)[:60]))
                errors += 1
            finally:
                shutil.rmtree(tmpdir)
    finally:
        ws._mtz_has_phase_columns = orig_check
        _PHASE_COLUMN_CACHE.clear()

    print("  Checked %d program×slot combinations" % total_checks)

    unexpected = [mi for mi in missing_inputs
                  if not mi['expected']]
    expected = [mi for mi in missing_inputs
                if mi['expected']]

    if unexpected:
        print("  UNEXPECTED missing inputs (%d):" %
              len(unexpected))
        for mi in unexpected:
            print("    %s / %s / %s → cats %s fallback %s"
                  % (mi['tutorial'], mi['program'],
                     mi['slot'], mi['categories'],
                     mi['fallback']))
    if expected:
        print("  Expected missing (test-env): %d" %
              len(expected))
    if errors:
        print("  Errors: %d tutorials crashed" % errors)

    if unexpected or errors:
        raise AssertionError(
            "%d unexpected missing inputs, %d errors"
            % (len(unexpected), errors))

    print("  PASSED (%d checks, %d expected missing)"
          % (total_checks, len(expected)))
    return missing_inputs


# =====================================================================
# CHECK 3: phased_data_mtz is in autobuild but NOT refine
# =====================================================================

def test_phased_data_mtz_not_in_refine():
    """phenix.refine input_priorities should NOT include
    phased_data_mtz (it only checks data_mtz/original_data_mtz).
    This is the Fix 5 pattern — documented, not a bug."""
    print("Test: phased_data_mtz_not_in_refine")

    with open(os.path.join(BASE_DIR,
                           'knowledge', 'programs.yaml')) as f:
        progs = yaml.safe_load(f)

    refine = progs.get('phenix.refine', {})
    refine_cats = refine.get('input_priorities', {}).get(
        'data_mtz', {}).get('categories', [])

    autobuild = progs.get('phenix.autobuild', {})
    autobuild_cats = autobuild.get('input_priorities', {}).get(
        'data_mtz', {}).get('categories', [])

    assert 'phased_data_mtz' not in refine_cats, (
        "phenix.refine should NOT check phased_data_mtz — "
        "this was the Fix 5 finding")
    assert 'phased_data_mtz' in autobuild_cats, (
        "phenix.autobuild should check phased_data_mtz")
    print("  PASSED (refine: %s, autobuild: %s)" % (
        refine_cats, autobuild_cats))


# =====================================================================
# RUN ALL
# =====================================================================

def run_all_tests():
    print("Phase 6: Category-Consumer Alignment")
    print()

    tests = [
        test_all_categories_exist,
        test_tutorial_programs_have_inputs,
        test_phased_data_mtz_not_in_refine,
    ]

    passed = 0
    failed = 0
    failures = []
    missing_inputs = []

    for t in tests:
        try:
            result = t()
            if isinstance(result, list):
                missing_inputs = result
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

    os.makedirs(FINDINGS_DIR, exist_ok=True)
    path = os.path.join(FINDINGS_DIR,
                        'phase_6_missing_inputs.yaml')
    with open(path, 'w') as f:
        f.write("# Phase 6: Category-Consumer Findings\n")
        f.write("phase: 6\n")
        f.write("name: category_consumer\n")
        f.write("status: %s\n" % status)
        f.write("passed: %d\n" % passed)
        f.write("failed: %d\n" % failed)
        if missing_inputs:
            f.write("\nmissing_inputs:\n")
            for mi in missing_inputs:
                f.write("  - tutorial: %s\n" % mi['tutorial'])
                f.write("    program: %s\n" % mi['program'])
                f.write("    slot: %s\n" % mi['slot'])
                f.write("    categories: %s\n"
                        % mi['categories'])
                if mi.get('fallback'):
                    f.write("    fallback_categories: %s\n"
                            % mi['fallback'])
                if mi.get('expected'):
                    f.write("    assessment: EXPECTED "
                            "(test-env limitation)\n")
        if failures:
            f.write("\ntest_failures:\n")
            for fail in failures:
                f.write("  - test: %s\n" % fail['test'])
                f.write("    error: \"%s\"\n"
                        % fail['error'][:100])

    print("  Findings: %s" % path)
    print("  Phase 6 overall: %s" % status)
    result = {'status': status}
    if status == 'FAIL':
        raise AssertionError(
            "Phase 6 category-consumer FAILED: %d test failures"
            % failed)
    return result


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] == 'PASS' else 1)
