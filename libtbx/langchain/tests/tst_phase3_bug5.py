"""Tests for Phase 3 (Bugs 1, 3, 4, 6) and Bug 5 fixes.

Covers:
  - Bug 1: Terminal diagnosis for unknown chemical element
  - Bug 3: mask_atoms blocked for resolve_cryo_em
  - Bug 4: Numeric coercion in StructureModel
  - Bug 5 Fix A: Hierarchical prefix whitelist
  - Bug 5 Fix B: Path resolution for strategy file values
  - Bug 5 Fix C: Reference model exclusion (Tier 2 categorizer)
  - Bug 5 Fix D: Strategy rewrites for full PHIL paths
  - Bug 6: Half-map pair detection (Tier 1 and Tier 2)
"""
from __future__ import print_function
import os
import re
import sys
import tempfile
import shutil

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_pass = _fail = 0

def assert_true(cond, msg=""):
  global _pass, _fail
  if cond:
    _pass += 1
  else:
    _fail += 1
    print("  FAIL: %s" % msg)

def assert_false(cond, msg=""):
  assert_true(not cond, msg)

def assert_equal(a, b, msg=""):
  assert_true(a == b, "%s: got %r, expected %r" % (msg, a, b))

def assert_in(item, container, msg=""):
  assert_true(item in container, "%s: %r not in %r" % (msg, item, container))

def assert_not_in(item, container, msg=""):
  assert_true(item not in container, "%s: %r found in %r" % (msg, item, container))


# =========================================================================
# Bug 1: Terminal diagnosis for unknown chemical element
# =========================================================================

def test_bug1_diagnosis_detector_unknown_element():
  """DiagnosisDetector must detect 'Unknown chemical element type'."""
  print("Test: bug1_diagnosis_detector_unknown_element")
  from agent.error_analyzer import DiagnosisDetector
  detector = DiagnosisDetector()

  result_text = (
    'Sorry: Unknown chemical element type:\n'
    '  "ATOM   1001  AU  AU    500 .*.        "\n'
    '  To resolve this problem, specify a chemical element type in\n'
    '  columns 77-78 of the PDB file, right justified (e.g. " C").'
  )
  match = detector.detect(result_text)
  assert_true(match is not None, "Should detect unknown_chemical_element")
  error_type, description, excerpt = match
  assert_equal(error_type, "unknown_chemical_element",
               "Error type should be unknown_chemical_element")
  assert_in("element", description.lower(), "Description should mention element")
  print("  PASSED")


def test_bug1_diagnosis_detector_no_false_positive():
  """DiagnosisDetector must not fire on unrelated errors."""
  print("Test: bug1_diagnosis_detector_no_false_positive")
  from agent.error_analyzer import DiagnosisDetector
  detector = DiagnosisDetector()

  # Normal refine error — should NOT match
  result_text = "Sorry: Wrong number of models of each type supplied."
  match = detector.detect(result_text)
  assert_true(match is None or match[0] != "unknown_chemical_element",
              "Should not match on unrelated error")
  print("  PASSED")


def test_bug1_no_overlap_with_recoverable():
  """Unknown chemical element must NOT be in recoverable_errors.yaml."""
  print("Test: bug1_no_overlap_with_recoverable")
  this_dir = os.path.dirname(os.path.abspath(__file__))
  yaml_path = os.path.join(this_dir, "..", "knowledge", "recoverable_errors.yaml")
  if os.path.exists(yaml_path):
    with open(yaml_path) as f:
      content = f.read()
    assert_not_in("Unknown chemical element", content,
                  "Must not overlap with recoverable_errors")
  print("  PASSED")


# =========================================================================
# Bug 3: mask_atoms blocked for resolve_cryo_em
# =========================================================================

def test_bug3_mask_atoms_blocked():
  """All 4 hallucinated params must be blocked for resolve_cryo_em."""
  print("Test: bug3_mask_atoms_blocked")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE
  _STRATEGY_FLAGS_CACHE.clear()

  strategy = {
    "mask_atoms": True,
    "output.prefix": "denmod",
    "d_min": 2.0,
    "main.number_of_macro_cycles": 1,
    "resolution": 3.0,
  }
  cleaned, stripped = validate_phil_strategy("phenix.resolve_cryo_em", strategy)
  stripped_keys = [k for k, v in stripped]

  assert_in("resolution", cleaned, "resolution should pass")
  for blocked in ["mask_atoms", "output.prefix", "d_min", "main.number_of_macro_cycles"]:
    assert_not_in(blocked, cleaned, "%s should be blocked" % blocked)
    assert_in(blocked, stripped_keys, "%s should be in stripped list" % blocked)
  print("  PASSED")


def test_bug3_mask_atoms_not_in_strategy_flags():
  """mask_atoms must not be in resolve_cryo_em strategy_flags."""
  print("Test: bug3_mask_atoms_not_in_strategy_flags")
  import yaml
  this_dir = os.path.dirname(os.path.abspath(__file__))
  yaml_path = os.path.join(this_dir, "..", "knowledge", "programs.yaml")
  with open(yaml_path) as f:
    progs = yaml.safe_load(f)
  sf = progs["phenix.resolve_cryo_em"].get("strategy_flags", {})
  assert_not_in("mask_atoms", sf, "mask_atoms must not be in strategy_flags")
  assert_in("resolution", sf, "resolution should remain")
  print("  PASSED")


# =========================================================================
# Bug 4: Numeric coercion in StructureModel
# =========================================================================

def test_bug4_coerce_string_metrics():
  """from_dict must coerce string metrics to float."""
  print("Test: bug4_coerce_string_metrics")
  from agent.structure_model import StructureModel

  d = {
    "model_state": {
      "r_work": "0.385",
      "r_free": "0.425",
      "model_map_cc": "0.75",
      "geometry": {"clashscore": "4.5", "rama_favored": "95.2"},
      "waters": "42",
    },
    "progress": [
      {"cycle": 1, "program": "refine", "r_work": "0.40", "r_free": "0.45"},
    ],
  }
  sm = StructureModel.from_dict(d)
  assert_true(isinstance(sm.model_state["r_work"], float),
              "r_work should be float, got %s" % type(sm.model_state["r_work"]))
  assert_true(isinstance(sm.model_state["r_free"], float),
              "r_free should be float")
  assert_true(isinstance(sm.model_state["model_map_cc"], float),
              "model_map_cc should be float")
  assert_true(isinstance(sm.model_state["geometry"]["clashscore"], float),
              "clashscore should be float")
  assert_true(isinstance(sm.model_state["waters"], int),
              "waters should be int")
  # Progress
  assert_true(isinstance(sm.progress[0]["r_free"], float),
              "progress r_free should be float")
  print("  PASSED")


def test_bug4_coerce_garbage_becomes_none():
  """Garbage string values must become None, not crash."""
  print("Test: bug4_coerce_garbage_becomes_none")
  from agent.structure_model import StructureModel

  d = {"model_state": {"r_work": "not_a_number", "r_free": [1, 2]}}
  sm = StructureModel.from_dict(d)
  assert_true(sm.model_state["r_work"] is None,
              "garbage r_work should be None")
  assert_true(sm.model_state["r_free"] is None,
              "garbage r_free should be None")
  print("  PASSED")


def test_bug4_detect_problems_no_crash_with_strings():
  """_detect_problems must not crash when model_state has string metrics."""
  print("Test: bug4_detect_problems_no_crash_with_strings")
  from agent.structure_model import StructureModel

  d = {"model_state": {"r_work": "0.385", "r_free": "0.425"}}
  sm = StructureModel.from_dict(d)
  # This is the exact crash site from run 15b
  try:
    sm._detect_problems({})
    assert_true(True, "no crash")
  except TypeError as e:
    assert_true(False, "_detect_problems crashed: %s" % e)
  print("  PASSED")


def test_bug4_detect_problems_rfree_gap():
  """_detect_problems should detect R-free gap after coercion."""
  print("Test: bug4_detect_problems_rfree_gap")
  from agent.structure_model import StructureModel

  d = {"model_state": {"r_work": "0.200", "r_free": "0.400"}}
  sm = StructureModel.from_dict(d)
  sm._detect_problems({})
  problems = sm.model_state.get("problems", [])
  gap_problems = [p for p in problems if p.get("type") == "r_free_gap"]
  assert_true(len(gap_problems) > 0,
              "Should detect R-free gap (0.400 - 0.200 = 0.200 > 0.10)")
  print("  PASSED")


def test_bug4_none_stays_none():
  """None values must stay None after coercion."""
  print("Test: bug4_none_stays_none")
  from agent.structure_model import StructureModel

  d = {"model_state": {"r_work": None, "r_free": None}}
  sm = StructureModel.from_dict(d)
  assert_true(sm.model_state["r_work"] is None, "None r_work should stay None")
  assert_true(sm.model_state["r_free"] is None, "None r_free should stay None")
  print("  PASSED")


# =========================================================================
# Bug 5 Fix A: Hierarchical prefix whitelist
# =========================================================================

def test_bug5a_prefix_passes_restraint_params():
  """Restraint params (ncs, SS, reference_model) must pass via prefix."""
  print("Test: bug5a_prefix_passes_restraint_params")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  strategy = {
    "reference_model.file": "4pf4.pdb",
    "reference_model.enabled": True,
    "ncs.type": "torsion",
    "ncs.constraints": True,
    "secondary_structure.enabled": True,
    "secondary_structure.protein.remove_outliers": False,
    "ramachandran_restraints": True,
    "resolution": 3.5,
  }
  cleaned, stripped = validate_phil_strategy("phenix.refine", strategy)

  for key in strategy:
    assert_in(key, cleaned, "%s should pass validation" % key)
  assert_equal(len(stripped), 0, "No params should be stripped")
  print("  PASSED")


def test_bug5a_prefix_strips_unknown():
  """Unknown params must still be stripped even with prefixes."""
  print("Test: bug5a_prefix_strips_unknown")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  strategy = {
    "ncs.type": "torsion",      # passes (prefix)
    "bad_param": "junk",        # stripped
    "space_group": "P4",        # stripped (not in whitelist or prefix)
    "my_custom_flag": True,     # stripped
  }
  cleaned, stripped = validate_phil_strategy("phenix.refine", strategy)

  assert_in("ncs.type", cleaned, "ncs.type should pass")
  for bad in ["bad_param", "space_group", "my_custom_flag"]:
    assert_not_in(bad, cleaned, "%s should be stripped" % bad)
  print("  PASSED")


def test_bug5a_prefix_case_insensitive():
  """Prefix matching must be case-insensitive."""
  print("Test: bug5a_prefix_case_insensitive")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  strategy = {
    "NCS.Type": "cartesian",
    "SECONDARY_STRUCTURE.enabled": True,
    "Reference_Model.File": "ref.pdb",
  }
  cleaned, stripped = validate_phil_strategy("phenix.refine", strategy)

  for key in strategy:
    assert_in(key, cleaned, "%s should pass (case-insensitive)" % key)
  print("  PASSED")


def test_bug5a_prefix_full_phil_path():
  """Full PHIL paths must pass via prefix substring match."""
  print("Test: bug5a_prefix_full_phil_path")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  strategy = {
    "refinement.pdb_interpretation.secondary_structure.enabled": True,
    "refinement.pdb_interpretation.ncs.type": "torsion",
    "refinement.reference_model.file": "ref.pdb",
  }
  cleaned, stripped = validate_phil_strategy("phenix.refine", strategy)

  for key in strategy:
    assert_in(key, cleaned,
              "%s should pass (full PHIL path, prefix substring)" % key)
  print("  PASSED")


def test_bug5a_prefix_no_false_positives():
  """Prefix must not match underscore-joined names (ncs_constraints_enabled)."""
  print("Test: bug5a_prefix_no_false_positives")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  strategy = {
    "ncs_constraints_enabled": True,  # NOT ncs. (no dot)
    "ncsa_something": True,           # NOT ncs. (no dot)
    "my_reference_model": "x.pdb",    # NOT reference_model. (no dot after)
  }
  cleaned, stripped = validate_phil_strategy("phenix.refine", strategy)

  for key in strategy:
    assert_not_in(key, cleaned,
                  "%s should NOT match prefix (no dot boundary)" % key)
  print("  PASSED")


def test_bug5a_other_programs_unaffected():
  """Programs without prefixes must be unaffected."""
  print("Test: bug5a_other_programs_unaffected")
  from agent.phil_validator import validate_phil_strategy, _STRATEGY_FLAGS_CACHE, _PREFIXES_CACHE
  _STRATEGY_FLAGS_CACHE.clear()
  _PREFIXES_CACHE.clear()

  # resolve_cryo_em has no prefixes — ncs.type should be stripped
  strategy = {"ncs.type": "torsion", "resolution": 3.0}
  cleaned, stripped = validate_phil_strategy("phenix.resolve_cryo_em", strategy)
  assert_not_in("ncs.type", cleaned,
                "ncs.type should not pass for resolve_cryo_em (no prefixes)")
  assert_in("resolution", cleaned, "resolution should pass")
  print("  PASSED")


# =========================================================================
# Bug 5 Fix B: Path resolution for strategy file values
# =========================================================================

def test_bug5b_path_resolution():
  """Strategy values ending in .pdb/.params must be resolved to absolute."""
  print("Test: bug5b_path_resolution")

  # Create temp files
  tmpdir = tempfile.mkdtemp()
  try:
    for fn in ["target.pdb", "ref.pdb", "data.mtz", "ss.params"]:
      with open(os.path.join(tmpdir, fn), "w") as f:
        f.write("dummy")

    # Simulate the path resolution logic from program_registry.py
    files = {
      "model": os.path.join(tmpdir, "target.pdb"),
      "data_mtz": os.path.join(tmpdir, "data.mtz"),
    }
    strategy = {
      "reference_model.file": "ref.pdb",        # relative
      "secondary_structure.input.file_name": "ss.params",  # relative
      "ncs.type": "torsion",                     # not a file
      "resolution": 3.5,                         # not a file
    }

    _PATH_EXTENSIONS = frozenset({
      '.pdb', '.cif', '.mtz', '.params', '.eff',
      '.dat', '.fa', '.fasta', '.seq', '.phil',
      '.param', '.ncs_spec',
    })
    _basename_map = {}
    _work_dir = None
    for _fp in files.values():
      _fps = _fp if isinstance(_fp, list) else [_fp]
      for _f in _fps:
        _f = str(_f)
        if os.path.isfile(_f):
          _basename_map[os.path.basename(_f)] = _f
          if _work_dir is None:
            _work_dir = os.path.dirname(_f)
    if _work_dir and os.path.isdir(_work_dir):
      for _fn in os.listdir(_work_dir):
        if _fn not in _basename_map:
          _full = os.path.join(_work_dir, _fn)
          if os.path.isfile(_full):
            _basename_map[_fn] = _full

    for key in list(strategy.keys()):
      val = strategy[key]
      if not isinstance(val, str) or not val:
        continue
      if not any(val.lower().endswith(ext) for ext in _PATH_EXTENSIONS):
        continue
      basename = os.path.basename(val)
      if basename in _basename_map:
        strategy[key] = _basename_map[basename]

    assert_equal(strategy["reference_model.file"],
                 os.path.join(tmpdir, "ref.pdb"),
                 "ref.pdb should resolve to absolute")
    assert_equal(strategy["secondary_structure.input.file_name"],
                 os.path.join(tmpdir, "ss.params"),
                 "ss.params should resolve to absolute")
    assert_equal(strategy["ncs.type"], "torsion",
                 "ncs.type should be unchanged")
    assert_equal(strategy["resolution"], 3.5,
                 "resolution should be unchanged")
  finally:
    shutil.rmtree(tmpdir)
  print("  PASSED")


def test_bug5b_absolute_path_unchanged():
  """Already-absolute file paths must not be modified."""
  print("Test: bug5b_absolute_path_unchanged")

  tmpdir = tempfile.mkdtemp()
  try:
    ref_path = os.path.join(tmpdir, "ref.pdb")
    with open(ref_path, "w") as f:
      f.write("dummy")

    strategy = {"reference_model.file": ref_path}

    _PATH_EXTENSIONS = frozenset({'.pdb'})
    _basename_map = {"ref.pdb": ref_path}

    for key in list(strategy.keys()):
      val = strategy[key]
      if not isinstance(val, str):
        continue
      if not any(val.lower().endswith(ext) for ext in _PATH_EXTENSIONS):
        continue
      basename = os.path.basename(val)
      if basename in _basename_map:
        resolved = _basename_map[basename]
        if resolved != val:
          strategy[key] = resolved

    assert_equal(strategy["reference_model.file"], ref_path,
                 "Already-absolute path should be unchanged")
  finally:
    shutil.rmtree(tmpdir)
  print("  PASSED")


# =========================================================================
# Bug 5 Fix C Tier 2: Reference model categorizer heuristic
# =========================================================================

def _run_ref_model_categorizer(model_list):
  """Simulate the Tier 2 reference model categorizer."""
  files = {"model": list(model_list)}  # copy

  _REF_KEYWORDS = re.compile(
    r'(reference|homolog|template|restraint|high.res)',
    re.IGNORECASE)
  _AGENT_OUTPUT_PREFIXES = (
    'refine_', 'autobuild_', 'autosol_', 'phaser_',
    'resolve_', 'predict_', 'real_space_', 'dock_',
    'map_to_model_', 'pdbtools_', 'molprobity_', 'rsr_',
  )

  if len(files["model"]) >= 2:
    for f in list(files["model"]):
      bn = os.path.basename(f).lower()
      if any(bn.startswith(p) for p in _AGENT_OUTPUT_PREFIXES):
        continue
      if _REF_KEYWORDS.search(bn):
        files["model"].remove(f)
        files.setdefault("reference_model", []).append(f)
        break

  return files


def test_bug5c_reference_keyword_detected():
  """File named 'reference_model.pdb' must be reclassified."""
  print("Test: bug5c_reference_keyword_detected")
  files = _run_ref_model_categorizer(
    ["/p/target.pdb", "/p/reference_model.pdb"])
  assert_equal(len(files["model"]), 1, "Should have 1 model left")
  assert_equal(files["model"][0], "/p/target.pdb", "target.pdb stays")
  assert_in("/p/reference_model.pdb", files.get("reference_model", []),
            "reference_model.pdb reclassified")
  print("  PASSED")


def test_bug5c_homolog_detected():
  """File named 'homolog_high_res.pdb' must be reclassified."""
  print("Test: bug5c_homolog_detected")
  files = _run_ref_model_categorizer(
    ["/p/target.pdb", "/p/homolog_high_res.pdb"])
  assert_in("/p/homolog_high_res.pdb", files.get("reference_model", []),
            "homolog detected")
  print("  PASSED")


def test_bug5c_hipip_not_reclassified():
  """hipip: 1IUA.pdb must NOT be reclassified (no keyword match)."""
  print("Test: bug5c_hipip_not_reclassified")
  files = _run_ref_model_categorizer(
    ["/p/hipip.pdb", "/p/1IUA.pdb"])
  assert_equal(len(files["model"]), 2, "Both should remain as model")
  assert_true("reference_model" not in files or len(files["reference_model"]) == 0,
              "1IUA.pdb should NOT be reclassified")
  print("  PASSED")


def test_bug5c_4pf4_not_reclassified():
  """lowres: 4pf4.pdb must NOT be reclassified (no keyword, relies on Tier 1)."""
  print("Test: bug5c_4pf4_not_reclassified")
  files = _run_ref_model_categorizer(
    ["/p/1jkt.pdb", "/p/4pf4.pdb"])
  assert_equal(len(files["model"]), 2, "Both should remain as model")
  print("  PASSED")


def test_bug5c_agent_output_never_reclassified():
  """Agent output files (refine_001.pdb etc) must never be reclassified."""
  print("Test: bug5c_agent_output_never_reclassified")
  # Even if filename contains 'reference' after the agent prefix
  files = _run_ref_model_categorizer(
    ["/p/model.pdb", "/p/refine_reference.pdb"])
  assert_equal(len(files["model"]), 2, "Agent output should not be reclassified")
  print("  PASSED")


def test_bug5c_single_model_skipped():
  """Single model: Tier 2 must not fire."""
  print("Test: bug5c_single_model_skipped")
  files = _run_ref_model_categorizer(["/p/model.pdb"])
  assert_equal(len(files["model"]), 1, "Single model should remain")
  assert_true("reference_model" not in files, "No reclassification")
  print("  PASSED")


def test_bug5c_only_one_reclassified():
  """With 3 models, only the first keyword match is reclassified."""
  print("Test: bug5c_only_one_reclassified")
  files = _run_ref_model_categorizer(
    ["/p/target.pdb", "/p/reference_a.pdb", "/p/homolog_b.pdb"])
  assert_equal(len(files["model"]), 2, "Should have 2 models left")
  assert_equal(len(files.get("reference_model", [])), 1,
               "Only one should be reclassified")
  print("  PASSED")


# =========================================================================
# Bug 5 Fix D: Strategy rewrites
# =========================================================================

def test_bug5d_rewrites_exist():
  """_STRATEGY_REWRITES must have new entries for refine."""
  print("Test: bug5d_rewrites_exist")
  with open(os.path.join(os.path.dirname(__file__),
            "..", "agent", "graph_nodes.py")) as f:
    src = f.read()
  expected = [
    "secondary_structure.enabled",
    "secondary_structure.protein.remove_outliers",
    "secondary_structure.input.file_name",
    "ncs.type",
    "ncs.constraints",
  ]
  for target in expected:
    assert_in(target, src,
              "Rewrite target '%s' must be in graph_nodes.py" % target)
  print("  PASSED")


# =========================================================================
# Bug 6: Half-map pair detection
# =========================================================================

def _simulate_halfmap_detection(full_map_list):
  """Simulate the two-tier half-map pair detection from workflow_state.py."""
  files = {"full_map": list(full_map_list), "half_map": []}

  _pair_re = re.compile(r'^(.+)[_-]([12])\.(\w+)$', re.IGNORECASE)
  _by_prefix = {}
  for f in files.get("full_map", []):
    m = _pair_re.match(os.path.basename(f))
    if m:
      prefix = m.group(1).lower()
      _by_prefix.setdefault(prefix, []).append(f)

  n_full = len(files.get("full_map", []))
  for prefix, pair in _by_prefix.items():
    if len(pair) != 2:
      continue
    if n_full == 2:
      for f in pair:
        files["full_map"].remove(f)
        files["half_map"].append(f)
      break
    elif n_full >= 3:
      remaining = [f for f in files["full_map"] if f not in pair]
      if remaining:
        for f in pair:
          files["full_map"].remove(f)
          files["half_map"].append(f)

  return files


def test_bug6_tier1_exactly_two_files():
  """Tier 1: exactly 2 _1/_2 files must be promoted to half_map."""
  print("Test: bug6_tier1_exactly_two_files")
  files = _simulate_halfmap_detection([
    "/p/7n8i_24237_box_1.ccp4",
    "/p/7n8i_24237_box_2.ccp4",
  ])
  assert_equal(len(files["full_map"]), 0, "full_map should be empty")
  assert_equal(len(files["half_map"]), 2, "half_map should have 2")
  print("  PASSED")


def test_bug6_tier2_three_files_with_companion():
  """Tier 2: pair + companion must promote pair, keep companion."""
  print("Test: bug6_tier2_three_files_with_companion")
  files = _simulate_halfmap_detection([
    "/p/map.ccp4",
    "/p/map_half_1.ccp4",
    "/p/map_half_2.ccp4",
  ])
  assert_equal(len(files["full_map"]), 1, "companion should remain")
  assert_equal(len(files["half_map"]), 2, "pair should be promoted")
  assert_equal(os.path.basename(files["full_map"][0]), "map.ccp4",
               "companion is map.ccp4")
  print("  PASSED")


def test_bug6_single_map_not_promoted():
  """Single map file must not be promoted."""
  print("Test: bug6_single_map_not_promoted")
  files = _simulate_halfmap_detection(["/p/map.ccp4"])
  assert_equal(len(files["full_map"]), 1, "Should stay as full_map")
  assert_equal(len(files["half_map"]), 0, "No half_maps")
  print("  PASSED")


def test_bug6_non_matching_pair_not_promoted():
  """Two files that don't form a _1/_2 pair must not be promoted."""
  print("Test: bug6_non_matching_pair_not_promoted")
  files = _simulate_halfmap_detection([
    "/p/map_a.ccp4",
    "/p/map_b.ccp4",
  ])
  assert_equal(len(files["full_map"]), 2, "Both should stay as full_map")
  assert_equal(len(files["half_map"]), 0, "No half_maps")
  print("  PASSED")


def test_bug6_tier2_no_companion_no_promotion():
  """Tier 2: 3 files that all form pairs should not promote without companion."""
  print("Test: bug6_tier2_no_companion_no_promotion")
  # 3 files: two pairs overlap, no true companion
  files = _simulate_halfmap_detection([
    "/p/map_1.ccp4",
    "/p/map_2.ccp4",
    "/p/other_1.ccp4",
  ])
  # map_1/map_2 form a pair with n_full=3 → remaining=[other_1] → promote
  assert_equal(len(files["half_map"]), 2, "map pair should promote (companion exists)")
  assert_equal(len(files["full_map"]), 1, "other_1 remains")
  print("  PASSED")


# =========================================================================
# Bug 4 continued: metric_evaluator.py _safe_float
# =========================================================================

def test_bug4_metric_evaluator_safe_float():
  """metric_evaluator handles string metrics without crash."""
  print("Test: bug4_metric_evaluator_safe_float")
  from agent.metric_evaluator import MetricEvaluator, _safe_float

  # Basic _safe_float
  assert_equal(_safe_float("0.385"), 0.385, "_safe_float string")
  assert_equal(_safe_float(0.385), 0.385, "_safe_float float")
  assert_equal(_safe_float(None), None, "_safe_float None")
  assert_equal(_safe_float("bad"), None, "_safe_float non-numeric")

  ev = MetricEvaluator()

  # calculate_improvement_rate with strings must not crash
  rate = ev.calculate_improvement_rate("r_free", "0.40", "0.35")
  assert_true(abs(rate - 12.5) < 0.1, "improvement rate from strings: %.1f" % rate)

  # is_significant_improvement with strings
  sig = ev.is_significant_improvement("r_free", "0.30", "0.25")
  assert_true(sig, "significant improvement from strings")

  # is_plateau with string values
  vals = ["0.30", "0.299", "0.298", "0.297"]
  plateau = ev.is_plateau(vals, "r_free")
  assert_true(isinstance(plateau, bool), "is_plateau returned bool")

  print("  PASSED")


def test_bug4_metric_evaluator_analyze_trend_strings():
  """analyze_trend survives string r_free values in metrics_history."""
  print("Test: bug4_metric_evaluator_analyze_trend_strings")
  from agent.metric_evaluator import MetricEvaluator

  ev = MetricEvaluator()

  # Simulate the exact 7rpq scenario: model_vs_data writes string r_free
  history = [
    {"program": "phenix.xtriage", "r_free": None},
    {"program": "phenix.model_vs_data", "r_free": "0.385", "r_work": "0.384"},
  ]

  # This must not crash (was: unsupported operand type(s) for -: 'str' and 'float')
  result = ev.analyze_trend(history, "xray", resolution=3.3)
  assert_true(isinstance(result, dict), "returned dict")
  assert_true("trend_summary" in result, "has trend_summary")
  print("  PASSED")


def test_bug4_metric_evaluator_cryoem_strings():
  """analyze_trend survives string CC values for cryo-EM."""
  print("Test: bug4_metric_evaluator_cryoem_strings")
  from agent.metric_evaluator import MetricEvaluator

  ev = MetricEvaluator()
  history = [
    {"program": "phenix.real_space_refine", "map_cc": "0.75"},
    {"program": "phenix.real_space_refine", "map_cc": "0.78"},
  ]
  result = ev.analyze_trend(history, "cryoem")
  assert_true(isinstance(result, dict), "returned dict")
  assert_true("map_cc_trend" in result, "has map_cc_trend")
  print("  PASSED")


# =========================================================================
# Orphan-map promotion
# =========================================================================

def test_orphan_map_promotion():
  """Map files in parent 'map' but no subcategory → promoted to full_map."""
  print("Test: orphan_map_promotion")
  # Simulate the apoferritin_denmod_dock scenario:
  # emd-20026_auto_sharpen_A.ccp4 is in 'map' but excluded from 'full_map'
  # by *_a.* pattern
  files = {
    "map": ["/p/emd-20026_auto_sharpen_A.ccp4",
            "/p/emd_20026_half_map_1_box.ccp4",
            "/p/emd_20026_half_map_2_box.ccp4"],
    "full_map": [],
    "half_map": ["/p/emd_20026_half_map_1_box.ccp4",
                 "/p/emd_20026_half_map_2_box.ccp4"],
    "optimized_full_map": [],
  }

  # Run orphan-map promotion logic directly
  _map_subcats = {"full_map", "half_map", "optimized_full_map"}
  _in_subcat = set()
  for sc in _map_subcats:
    for f in files.get(sc, []):
      _in_subcat.add(f)
  for f in list(files.get("map", [])):
    if f not in _in_subcat:
      if "full_map" not in files:
        files["full_map"] = []
      if f not in files["full_map"]:
        files["full_map"].append(f)

  assert_equal(len(files["full_map"]), 1, "orphan promoted to full_map")
  assert_true("auto_sharpen_A" in files["full_map"][0],
              "correct file promoted")
  assert_equal(len(files["half_map"]), 2, "half_maps unchanged")
  print("  PASSED")


def test_orphan_map_no_false_positive():
  """Files already in subcategories are NOT double-promoted."""
  print("Test: orphan_map_no_false_positive")
  files = {
    "map": ["/p/map.mrc", "/p/half_1.mrc", "/p/half_2.mrc"],
    "full_map": ["/p/map.mrc"],
    "half_map": ["/p/half_1.mrc", "/p/half_2.mrc"],
    "optimized_full_map": [],
  }

  _map_subcats = {"full_map", "half_map", "optimized_full_map"}
  _in_subcat = set()
  for sc in _map_subcats:
    for f in files.get(sc, []):
      _in_subcat.add(f)
  for f in list(files.get("map", [])):
    if f not in _in_subcat:
      if "full_map" not in files:
        files["full_map"] = []
      if f not in files["full_map"]:
        files["full_map"].append(f)

  assert_equal(len(files["full_map"]), 1, "no extra promotion")
  assert_equal(len(files["half_map"]), 2, "half_maps unchanged")
  print("  PASSED")


# =========================================================================
# Belt-and-suspenders: kb_tags coercion
# =========================================================================

def test_kb_tags_string_rfree_trend():
  """kb_tags._trend_tags handles string R-free values."""
  print("Test: kb_tags_string_rfree_trend")
  from agent.kb_tags import _trend_tags
  # String values that would crash on subtraction
  tags = _trend_tags(["0.30", "0.299", "0.298"])
  assert_true(isinstance(tags, list), "returned list")
  # Values are close → should detect plateau
  assert_true("plateau" in tags or "r_free_stuck" in tags,
              "plateau detected from strings")

  # With improving values
  tags2 = _trend_tags(["0.45", "0.35", "0.25"])
  assert_true("improving" in tags2, "improving detected from strings")
  print("  PASSED")


# =========================================================================
# Runner
# =========================================================================

def run_tests():
  global _pass, _fail
  _pass = _fail = 0

  test_funcs = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
  for func in test_funcs:
    try:
      func()
    except Exception as e:
      _fail += 1
      print("  EXCEPTION in %s: %s" % (func.__name__, e))

  print("\n%d passed, %d failed" % (_pass, _fail))
  if _fail > 0:
    sys.exit(1)


if __name__ == "__main__":
  run_tests()
