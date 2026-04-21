"""
Phase 2: Path Consistency — YAML vs Hardcoded Categorization.

For 10 representative tutorials, calls _categorize_files twice:
  Path A: with real YAML rules (production path)
  Path B: with _load_category_rules monkeypatched to return None
          (forces hardcoded fallback path)

Both paths run the full post-processing pipeline (bubble-up,
orphan promotion, phased promotion, ignored formats). Diffs
the output categories to detect path-specific bugs.

Run: python tests/tst_phase2_path_consistency.py
Produces: findings/phase_2_path_divergences.yaml
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import tempfile
import shutil

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

# Mock libtbx imports
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
    _categorize_files,
    _load_category_rules,
)
import agent.workflow_state as ws


# =====================================================================
# TEST DATA: 10 representative tutorials
# =====================================================================

TUTORIALS = {
    'rab3a-refine': {
        'type': 'xray',
        'files': ['rab3a.pdb', 'rab3a_scale.hkl',
                  'rab3a_phases.hkl'],
        'note': 'Two .hkl files — phases bug target',
    },
    'nsf-d2-refine': {
        'type': 'xray',
        'files': ['nsf-d2_start.pdb', 'nsf-d2_reference.pdb',
                  'nsf-d2_phases.hkl', 'nsf-d2.cv'],
        'note': 'Single .hkl + .cv — phases bug target',
    },
    'hipip-refine': {
        'type': 'xray',
        'files': ['hipip.mtz', 'hipip.pdb', '1IUA.pdb'],
        'note': '.mtz + .pdb — simple refine',
    },
    'lysozyme-MRSAD': {
        'type': 'xray',
        'files': ['1fkq_prot.pdb', 'hewl.seq',
                  'lyso2001_scala1.mtz'],
        'note': '.mtz + .pdb + .seq — MR-SAD',
    },
    'p9-build': {
        'type': 'xray',
        'files': ['seq.dat', 'p9_hires.mtz',
                  'p9_data_and_phases.mtz'],
        'note': 'Autobuild — two MTZ files',
    },
    'gene-5-mad': {
        'type': 'xray',
        'files': ['sequence.dat', 'high.sca', 'peak.sca',
                  'infl.sca'],
        'note': '.sca + .dat — MAD',
    },
    'porin-twin': {
        'type': 'xray',
        'files': ['porin.mtz', 'porin.pdb', 'porin.cv'],
        'note': '.mtz + .pdb + .cv — twinned',
    },
    'beta-blip': {
        'type': 'xray',
        'files': ['beta.pdb', 'blip.pdb', 'beta.seq',
                  'blip.seq', 'beta_blip_P3221.mtz'],
        'note': '.mtz + 2 .pdb + 2 .seq — MR',
    },
    'real-space-refine-5ljv': {
        'type': 'cryoem',
        'files': ['model.pdb', 'map.ccp4'],
        'note': 'Cryo-EM RSR — .ccp4 + .pdb',
    },
    'apoferritin_model_building': {
        'type': 'cryoem',
        'files': ['apoferritin.seq',
                  'apoferritin_denmod_box.ccp4',
                  'apoferritin_chainA_rsr.pdb'],
        'note': 'Cryo-EM model building — .ccp4 + .pdb + .seq',
    },
}

# Key categories to compare between paths
KEY_CATEGORIES = [
    'data_mtz', 'phased_data_mtz', 'map_coeffs_mtz',
    'model', 'search_model', 'pdb',
    'ligand_cif', 'ligand_pdb',
    'map', 'full_map', 'half_map',
    'sequence',
    'ignored_formats',
]


# =====================================================================
# HELPERS
# =====================================================================

def create_tutorial_files(tmpdir, file_list):
    """Create temp files with realistic content and sizes."""
    paths = []
    for name in file_list:
        path = os.path.join(tmpdir, name)
        ext = os.path.splitext(name)[1].lower()
        with open(path, 'w') as f:
            if ext == '.pdb':
                # Protein-sized PDB (>10KB to avoid ligand_pdb)
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
                    f.write("1 0 0 45.2 0.89 120.5 3.2\n")
                else:
                    f.write("COLUMN_LABELS H K L FP SIGFP\n")
                    f.write("1 0 0 120.5 3.2\n")
            elif ext == '.mtz':
                # Binary-ish (iotbx won't parse in tests anyway)
                f.write("MTZ_DUMMY\n" * 100)
            elif ext in ('.ccp4', '.mrc', '.map'):
                # Map files need minimal content
                f.write("MAP " + "\x00" * 100 + "\n")
            elif ext == '.cv':
                f.write("COLUMN_LABELS H K L FREE\n")
                f.write("1 0 0 1\n")
            elif ext == '.cif':
                f.write("data_ligand\n_chem_comp.id LIG\n")
            else:
                f.write("dummy content\n")
        paths.append(path)
    return paths


def extract_basenames(file_dict, category):
    """Get sorted basenames from a category."""
    files = file_dict.get(category, [])
    if not isinstance(files, list):
        return []
    result = []
    for f in files:
        if isinstance(f, str):
            result.append(os.path.basename(f))
        elif isinstance(f, dict):
            # ignored_formats entries are dicts with 'file' key
            name = f.get('file', '')
            if name:
                result.append(os.path.basename(str(name)))
    return sorted(result)


def compare_paths(yaml_result, hard_result, categories):
    """Compare two categorization results across key categories.

    Returns list of divergences.
    """
    divergences = []
    for cat in categories:
        yaml_files = extract_basenames(yaml_result, cat)
        hard_files = extract_basenames(hard_result, cat)

        if yaml_files != hard_files:
            yaml_only = set(yaml_files) - set(hard_files)
            hard_only = set(hard_files) - set(yaml_files)
            divergences.append({
                'category': cat,
                'yaml_files': yaml_files,
                'hard_files': hard_files,
                'yaml_only': sorted(yaml_only),
                'hard_only': sorted(hard_only),
            })
    return divergences


# =====================================================================
# TESTS
# =====================================================================

def run_all_tests():
    """Run path consistency tests for all 10 tutorials."""
    print("Phase 2: Path Consistency — YAML vs Hardcoded")
    print()

    # Verify YAML rules are available
    rules = _load_category_rules()
    if not rules:
        print("  ERROR: Cannot load file_categories.yaml — "
              "cannot run YAML path tests")
        raise AssertionError(
            "Phase 2: file_categories.yaml not loadable")

    # Monkeypatch _mtz_has_phase_columns for test environment
    orig_check = ws._mtz_has_phase_columns
    orig_load_rules = ws._load_category_rules
    def mock_check(f):
        return False  # iotbx not available
    ws._mtz_has_phase_columns = mock_check

    passed = 0
    failed = 0
    errors = 0
    all_divergences = {}

    try:
        for tutorial_name, tutorial_def in sorted(
                TUTORIALS.items()):
            print("  Test: %-30s" % tutorial_name, end="")

            tmpdir = tempfile.mkdtemp()
            try:
                files = create_tutorial_files(
                    tmpdir, tutorial_def['files'])

                # Path A: YAML — normal _categorize_files (uses
                # YAML rules via _load_category_rules)
                ws._PHASE_COLUMN_CACHE.clear()
                ws._load_category_rules = orig_load_rules
                yaml_result = _categorize_files(files)

                # Path B: Hardcoded — force _load_category_rules
                # to return None so _categorize_files falls back
                # to _categorize_files_hardcoded
                ws._PHASE_COLUMN_CACHE.clear()
                ws._load_category_rules = lambda: None
                hard_result = _categorize_files(files)

                divergences = compare_paths(
                    yaml_result, hard_result, KEY_CATEGORIES)

                if divergences:
                    print("DIVERGENCE (%d categories)" %
                          len(divergences))
                    failed += 1
                    all_divergences[tutorial_name] = divergences
                    for d in divergences:
                        print("    %s: yaml=%s hard=%s" % (
                            d['category'],
                            d['yaml_files'],
                            d['hard_files']))
                else:
                    print("PASS")
                    passed += 1

            except Exception as e:
                print("ERROR: %s" % str(e)[:60])
                errors += 1
            finally:
                shutil.rmtree(tmpdir)

    finally:
        ws._mtz_has_phase_columns = orig_check
        ws._load_category_rules = orig_load_rules
        ws._PHASE_COLUMN_CACHE.clear()

    # Known expected divergences (YAML is more correct in both cases).
    # These are whitelisted — the test fails on any NEW divergence.
    _EXPECTED_DIVERGENCES = {
        ('nsf-d2-refine', 'pdb'):
            'EXPECTED: YAML puts nsf-d2_reference.pdb in '
            'intermediate (*reference* pattern). Hardcoded '
            'lacks this pattern.',
        ('p9-build', 'data_mtz'):
            'EXPECTED: YAML promotes p9_data_and_phases.mtz '
            'to phased_data_mtz (*phases* pattern). Hardcoded '
            'relies on content check (mocked False).',
        ('p9-build', 'phased_data_mtz'):
            'EXPECTED: Same root cause as data_mtz divergence.',
    }

    # Count unexpected divergences (not in whitelist)
    unexpected = []
    for tut, divs in all_divergences.items():
        for d in divs:
            key = (tut, d['category'])
            if key not in _EXPECTED_DIVERGENCES:
                unexpected.append(
                    '%s/%s' % (tut, d['category']))

    if unexpected or errors > 0:
        status = 'FAIL'
    elif failed > 0:
        status = 'PARTIAL'  # Only expected divergences
    else:
        status = 'PASS'

    print()
    print("  Results: %d passed, %d divergent, %d errors"
          % (passed, failed, errors))
    if unexpected:
        print("  UNEXPECTED divergences: %s" % unexpected)
    elif failed > 0:
        print("  (all divergences are whitelisted/expected)")

    # Write findings
    findings_dir = os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), 'findings')
    os.makedirs(findings_dir, exist_ok=True)
    findings_path = os.path.join(
        findings_dir, 'phase_2_path_divergences.yaml')
    with open(findings_path, 'w') as f:
        f.write("# Phase 2: Path Consistency Findings\n")
        f.write("phase: 2\n")
        f.write("name: path_consistency\n")
        f.write("status: %s\n" % status)
        f.write("tutorials_tested: %d\n" % len(TUTORIALS))
        f.write("passed: %d\n" % passed)
        f.write("divergent: %d\n" % failed)
        if all_divergences:
            f.write("\ndivergences:\n")
            for tut, divs in sorted(all_divergences.items()):
                f.write("  %s:\n" % tut)
                for d in divs:
                    f.write("    - category: %s\n"
                            % d['category'])
                    f.write("      yaml_only: %s\n"
                            % d['yaml_only'])
                    f.write("      hard_only: %s\n"
                            % d['hard_only'])
                    note = _EXPECTED_DIVERGENCES.get(
                        (tut, d['category']))
                    if note:
                        f.write("      assessment: %s\n"
                                % note)

    print("  Findings: %s" % findings_path)
    print("  Phase 2 overall: %s" % status)

    result = {
        'status': status,
        'passed': passed,
        'failed': failed,
        'divergences': all_divergences,
    }
    if status == 'FAIL':
        parts = []
        if unexpected:
            parts.append("unexpected divergences: %s" % unexpected)
        if errors > 0:
            parts.append("%d tutorial errors" % errors)
        raise AssertionError(
            "Phase 2 path consistency FAILED: %s"
            % '; '.join(parts))
    return result


if __name__ == "__main__":
    result = run_all_tests()
    sys.exit(0 if result['status'] != 'FAIL' else 1)
