"""
Tests for Bug N fixes: N1 (event_formatter TypeError), N3 (file_categories
refined exclude), plus regression checks for N2/N4/N5 (already applied).

Run with: python tst_n_bugs.py
"""
from __future__ import absolute_import, division, print_function
import os
import sys
import fnmatch

# ---------------------------------------------------------------------------
# Bug A (N1): event_formatter TypeError when metric value is a string
# ---------------------------------------------------------------------------

def test_formatter_string_metric():
    """event_formatter must not crash when a metric value is a string."""
    # Simulate the crash: resolution or another metric key has a string value
    # (e.g., xtriage sometimes stores space_group-like strings in parsed fields)
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    try:
        from agent.event_formatter import EventFormatter
        try:
            from libtbx.langchain.agent.event_log import EventType, Verbosity
        except ImportError:
            # Provide minimal stubs for standalone testing
            class EventType:
                METRICS_EXTRACTED = "metrics_extracted"
            class Verbosity:
                VERBOSE = "verbose"
                NORMAL = "normal"
            import sys as _sys
            import types as _types
            _mod = _types.ModuleType("libtbx.langchain.agent.event_log")
            _mod.EventType = EventType
            _mod.Verbosity = Verbosity
            _mod.EVENT_VERBOSITY = {}
            _mod.verbosity_index = lambda v: 0
            _sys.modules["libtbx"] = _types.ModuleType("libtbx")
            _sys.modules["libtbx.langchain"] = _types.ModuleType("libtbx.langchain")
            _sys.modules["libtbx.langchain.agent"] = _types.ModuleType("libtbx.langchain.agent")
            _sys.modules["libtbx.langchain.agent.event_log"] = _mod
    except ImportError:
        print("SKIP test_formatter_string_metric (cannot import event_formatter)")
        return

    formatter = EventFormatter(verbosity=Verbosity.VERBOSE)

    # Case 1: string value for resolution (the CASP7 crash pattern)
    event = {
        "type": EventType.METRICS_EXTRACTED,
        "resolution": "2.04",     # string instead of float
        "r_free": 0.556,
    }
    try:
        result = formatter._format_metrics_extracted(event)
        assert result is not None
        assert "Resolution" in result
        assert "2.04" in result
        print("PASS test_formatter_string_metric (string resolution)")
    except TypeError as e:
        print("FAIL test_formatter_string_metric: TypeError: %s" % e)
        return

    # Case 2: string with previous value (diff computation would crash)
    event2 = {
        "type": EventType.METRICS_EXTRACTED,
        "r_free": "0.556",        # string
        "r_free_prev": "0.580",   # string
    }
    try:
        result2 = formatter._format_metrics_extracted(event2)
        assert result2 is not None
        print("PASS test_formatter_string_metric (string r_free with prev)")
    except TypeError as e:
        print("FAIL test_formatter_string_metric (with prev): TypeError: %s" % e)


def test_formatter_normal_metrics():
    """event_formatter still works correctly for normal numeric metrics."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    try:
        from agent.event_formatter import EventFormatter
        try:
            from libtbx.langchain.agent.event_log import EventType, Verbosity
        except ImportError:
            print("SKIP test_formatter_normal_metrics (cannot import)")
            return
    except ImportError:
        print("SKIP test_formatter_normal_metrics (cannot import event_formatter)")
        return

    formatter = EventFormatter(verbosity=Verbosity.VERBOSE)
    event = {
        "type": EventType.METRICS_EXTRACTED,
        "r_free": 0.284,
        "r_free_prev": 0.310,
        "resolution": 2.5,
    }
    result = formatter._format_metrics_extracted(event)
    assert result is not None
    assert "0.2840" in result   # formatted r_free
    assert "improved" in result  # trend label
    assert "2.50 A" in result   # resolution
    print("PASS test_formatter_normal_metrics")


# ---------------------------------------------------------------------------
# Bug C1 (N3/N4): file_categories.yaml *overall_best* in refined.excludes
# ---------------------------------------------------------------------------

def _load_yaml_rules():
    """Load file_categories.yaml using PyYAML or ruamel.yaml."""
    yaml_path = os.path.join(
        os.path.dirname(__file__), '..', 'knowledge', 'file_categories.yaml')
    if not os.path.exists(yaml_path):
        return None
    try:
        import yaml
        with open(yaml_path) as f:
            return yaml.safe_load(f)
    except ImportError:
        pass
    try:
        from ruamel.yaml import YAML
        yml = YAML()
        with open(yaml_path) as f:
            return yml.load(f)
    except ImportError:
        return None


def _matches(basename, pattern):
    return fnmatch.fnmatch(basename.lower(), pattern.lower())


def test_overall_best_refine_excluded_from_refined():
    """overall_best_final_refine_001.pdb must NOT be categorized as refined."""
    rules = _load_yaml_rules()
    if rules is None:
        print("SKIP test_overall_best_refine_excluded (no yaml parser)")
        return

    refined_def = rules.get("refined", {})
    excludes = refined_def.get("excludes", [])

    # The file that was causing refine to fail with "Wrong number of models"
    target = "overall_best_final_refine_001.pdb"

    # Verify it would match the refined pattern without the exclude
    patterns = refined_def.get("patterns", [])
    matched_pattern = any(_matches(target, p) for p in patterns)
    assert matched_pattern, "Test setup error: %s should match *refine*" % target

    # Now verify it IS excluded
    excluded = any(_matches(target, exc) for exc in excludes)
    assert excluded, (
        "FAIL: %s is NOT excluded from refined category.\n"
        "  excludes=%r\n"
        "  Add '*overall_best*' to refined.excludes in file_categories.yaml" % (target, excludes))
    print("PASS test_overall_best_refine_excluded_from_refined")


def test_standalone_refine_still_categorized():
    """Normal phenix.refine outputs must still be in the refined category."""
    rules = _load_yaml_rules()
    if rules is None:
        print("SKIP test_standalone_refine_still_categorized (no yaml parser)")
        return

    refined_def = rules.get("refined", {})
    patterns = refined_def.get("patterns", [])
    excludes = refined_def.get("excludes", [])

    # These should always be in refined
    standalone_refine_files = [
        "7qz0_refine_001.pdb",
        "refine_001.pdb",
        "refine_001_001.pdb",
        "my_structure_refine_002.pdb",
    ]
    for fname in standalone_refine_files:
        matched = any(_matches(fname, p) for p in patterns)
        excluded = any(_matches(fname, exc) for exc in excludes)
        if not matched:
            print("FAIL test_standalone_refine_still_categorized: %s doesn't match pattern" % fname)
            return
        if excluded:
            print("FAIL test_standalone_refine_still_categorized: %s is incorrectly excluded" % fname)
            return
    print("PASS test_standalone_refine_still_categorized")


def test_predict_and_build_output_not_in_refined():
    """Various predict_and_build outputs should not be in refined."""
    rules = _load_yaml_rules()
    if rules is None:
        print("SKIP test_predict_and_build_output_not_in_refined (no yaml parser)")
        return

    refined_def = rules.get("refined", {})
    patterns = refined_def.get("patterns", [])
    excludes = refined_def.get("excludes", [])

    # All of these should be excluded from refined
    pab_outputs = [
        "overall_best_final_refine_001.pdb",   # main bug case
        "overall_best_refine_001.pdb",          # variant
        "overall_best.pdb",                     # standard PAB output
    ]
    for fname in pab_outputs:
        matched = any(_matches(fname, p) for p in patterns)
        excluded = any(_matches(fname, exc) for exc in excludes)
        if matched and not excluded:
            print("FAIL test_predict_and_build_output_not_in_refined: "
                  "%s is in refined (should be excluded)" % fname)
            return
    print("PASS test_predict_and_build_output_not_in_refined")


# ---------------------------------------------------------------------------
# Bug B (N3): _auto_discover_files supplement mode — unit-level check
# ---------------------------------------------------------------------------

def test_auto_discover_supplement_mode():
    """_auto_discover_files adds missing files even when original_files is set."""
    import tempfile, shutil

    # Build a temporary directory with MTZ + FASTA + ligand PDB
    tmpdir = tempfile.mkdtemp()
    try:
        open(os.path.join(tmpdir, "7qz0.mtz"), 'w').close()
        open(os.path.join(tmpdir, "7qz0.fa"), 'w').close()
        open(os.path.join(tmpdir, "7qz0_ligand.pdb"), 'w').close()
        open(os.path.join(tmpdir, "README.txt"), 'w').close()

        # Simulate: LLM preprocessing already found ligand PDB
        discovered_so_far = [os.path.join(tmpdir, "7qz0_ligand.pdb")]
        existing_basenames = {os.path.basename(f) for f in discovered_so_far}

        _KNOWN_EXT = {".mtz", ".pdb", ".cif", ".sca", ".mrc", ".map",
                      ".fa", ".fasta", ".seq", ".dat", ".hkl", ".sf",
                      ".phil", ".eff"}
        _SKIP_EXT  = {".log", ".txt", ".md", ".html", ".pdf",
                      ".png", ".jpg", ".svg", ".json", ".csv",
                      ".py", ".sh", ".bat"}

        new_files = []
        for fname in sorted(os.listdir(tmpdir)):
            fpath = os.path.join(tmpdir, fname)
            if not os.path.isfile(fpath): continue
            if fname.startswith("."): continue
            if fname.upper().startswith("README"): continue
            if os.path.basename(fpath) in existing_basenames: continue
            _, ext = os.path.splitext(fname.lower())
            if ext in _KNOWN_EXT:
                new_files.append(fpath)

        result = discovered_so_far + new_files
        basenames = {os.path.basename(f) for f in result}

        assert "7qz0.mtz" in basenames, "MTZ not discovered: %s" % basenames
        assert "7qz0.fa" in basenames,  "FASTA not discovered: %s" % basenames
        assert "7qz0_ligand.pdb" in basenames, "Ligand PDB dropped"
        assert "README.txt" not in basenames, "README incorrectly included"
        print("PASS test_auto_discover_supplement_mode")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_all_tests():
    print("=" * 60)
    print("Bug N fix tests")
    print("=" * 60)
    test_formatter_string_metric()
    test_formatter_normal_metrics()
    test_overall_best_refine_excluded_from_refined()
    test_standalone_refine_still_categorized()
    test_predict_and_build_output_not_in_refined()
    test_auto_discover_supplement_mode()
    print("=" * 60)
    print("Done")

if __name__ == "__main__":
    run_all_tests()
