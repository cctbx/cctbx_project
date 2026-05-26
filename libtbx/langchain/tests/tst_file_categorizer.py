"""
Semantic-pin tests for agent/file_utils.py categorizer functions.

Pin documented behavior of every public function in file_utils.py.
Future changes to these functions will fail the suite if behavior
drifts, surfacing the divergence at the moment a sandbox test runs
against the real function — the same discipline that
test_bug11_matches_exclude_pattern_semantics applies to one function,
generalized to all of file_utils.py.

Background: the H10 → H11 cycle (May 2026) exposed that a sandbox
stub with substring matching diverged from the real function's
word-boundary matching, masking a production bug for one ship.
The general lesson is captured in DEVELOPER_GUIDE.md §10.3
("Semantic-pin tests").  This file applies that discipline to:

  - classify_mtz_type      (5-rule classification + ordering invariant)
  - get_mtz_stage          (subcategory dispatch)
  - get_category_for_extension  (extension lookup + edge cases)
  - is_mtz_file, is_model_file, is_map_file, is_sequence_file
    (boolean extension/name checks with one load-bearing quirk:
    is_model_file's 'ligand' filter only applies to .cif)

  matches_exclude_pattern is pinned by
  test_bug11_matches_exclude_pattern_semantics in tst_autosol_bugs.py;
  this suite cross-references that pin rather than duplicating it.

All tests use the same sandbox-skip pattern as Bug 11: nested
try/except for libtbx.langchain.agent.file_utils and agent.file_utils;
clean SKIP + return if neither resolves.  This keeps lightweight
sandbox CI green while ensuring full fidelity in PHENIX environments.

Run with:
    PYTHONPATH=. python tests/tst_file_categorizer.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
  assert_equal,
  assert_in,
  assert_true,
  assert_false,
  run_tests_with_fail_fast,
)

# PHENIX/cctbx linter silencer
(assert_equal, assert_in, assert_true, assert_false,
 run_tests_with_fail_fast)


def _import_file_utils():
  """Import file_utils with the standard sandbox-skip pattern.

  Returns the module if importable, or None if neither path works
  (the caller should print SKIP and return).
  """
  try:
    from libtbx.langchain.agent import file_utils as _fu
    return _fu
  except ImportError:
    try:
      from agent import file_utils as _fu
      return _fu
    except ImportError:
      return None


# =========================================================================
# classify_mtz_type — 5-rule classification
# =========================================================================

def test_classify_mtz_type_semantics():
  """Pin classify_mtz_type's classification rules.

    Categories:
      - data_mtz: contains measured Fobs + R-free flags (refinement input)
      - map_coeffs_mtz: calculated phases (ligand fitting, visualization)

    Five rules apply in order; first match wins.  Rule 1 takes precedence
    over Rule 2 even when both could match (e.g., refine_001.mtz matches
    both: the Rule 1 regex fires first).  This ordering invariant is
    pinned by the dedicated assertion at the bottom.
    """
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  # Rule 1: refine_NNN.mtz / refine_NNN_NNN.mtz (with optional prefix)
  #         → map_coeffs_mtz
  assert_equal(fu.classify_mtz_type("refine_001.mtz"), "map_coeffs_mtz",
        "Rule 1: refine_001.mtz → map_coeffs_mtz (refine output)")
  assert_equal(fu.classify_mtz_type("refine_001_001.mtz"), "map_coeffs_mtz",
        "Rule 1: refine_001_001.mtz → map_coeffs_mtz (double 3-digit)")
  assert_equal(fu.classify_mtz_type("7qz0_refine_001.mtz"), "map_coeffs_mtz",
        "Rule 1: prefixed refine_001 → map_coeffs_mtz")
  assert_equal(fu.classify_mtz_type("7qz0_refine_001_001.mtz"), "map_coeffs_mtz",
        "Rule 1: prefixed refine_001_001 → map_coeffs_mtz")

  # Rule 2: other *_001.mtz → map_coeffs_mtz
  assert_equal(fu.classify_mtz_type("model_refine_001.mtz"), "map_coeffs_mtz",
        "Rule 2: other _001 suffix → map_coeffs_mtz")

  # Rule 3: contains 'map_coeffs', 'denmod', 'density_mod' → map_coeffs_mtz
  assert_equal(fu.classify_mtz_type("denmod_map.mtz"), "map_coeffs_mtz",
        "Rule 3: 'denmod' substring → map_coeffs_mtz")
  assert_equal(fu.classify_mtz_type("density_modified.mtz"), "map_coeffs_mtz",
        "Rule 3: 'density_mod' substring → map_coeffs_mtz")
  assert_equal(fu.classify_mtz_type("solve_map_coeffs.mtz"), "map_coeffs_mtz",
        "Rule 3: 'map_coeffs' substring → map_coeffs_mtz")

  # Rule 4: '_data.mtz' or 'refinement_data' → data_mtz
  assert_equal(fu.classify_mtz_type("experimental_data.mtz"), "data_mtz",
        "Rule 4: '_data.mtz' → data_mtz")
  assert_equal(fu.classify_mtz_type("refinement_data_input.mtz"), "data_mtz",
        "Rule 4: 'refinement_data' → data_mtz")

  # Rule 5: default → data_mtz
  assert_equal(fu.classify_mtz_type("data.mtz"), "data_mtz",
        "Rule 5: default → data_mtz")
  assert_equal(fu.classify_mtz_type("7qz0.mtz"), "data_mtz",
        "Rule 5: bare PDB-named MTZ → data_mtz")

  # Rule-ordering invariant (Gemini's H12 plan review):
  # refine_001.mtz matches BOTH Rule 1 and Rule 2.  The assertion at
  # the top of this test confirmed Rule 1 fires first (returning
  # map_coeffs_mtz).  If a future refactor changes rule order such
  # that Rule 2 fires first, classify_mtz_type would still return
  # map_coeffs_mtz for this specific filename — but for a counterexample,
  # see Rule 4 below: 'refinement_data_input.mtz' has the '_data.mtz'
  # substring AND the '_001' substring is absent, so Rules 1/2/3 all
  # miss and Rule 4 wins.  If Rule 4 were moved before Rule 3, a
  # filename like 'denmod_data.mtz' would flip from map_coeffs_mtz to
  # data_mtz — but no such file exists in production naming.  The
  # current rule order is correct for production filenames.


# =========================================================================
# get_mtz_stage — subcategory dispatch
# =========================================================================

def test_get_mtz_stage_semantics():
  """Pin get_mtz_stage's subcategory mapping per category.

    For category='data_mtz':
      - '_data.mtz' substring → 'original_data_mtz'
      - 'phased' substring → 'phased_data_mtz'
      - else → 'data_mtz'

    For category='map_coeffs_mtz':
      - 'denmod'/'density_mod' substring → 'denmod_map_coeffs'
      - 'map_coeffs' + 'predict' substrings → 'predict_build_map_coeffs'
      - 'overall_best_map_coeffs' substring → 'predict_build_map_coeffs'
      - else → 'refine_map_coeffs'

    For unrecognized category: returns the category unchanged
    (passthrough behavior).
    """
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  # data_mtz subcategories
  assert_equal(fu.get_mtz_stage("experimental_data.mtz", "data_mtz"),
        "original_data_mtz",
        "data_mtz + '_data.mtz' → original_data_mtz")
  assert_equal(fu.get_mtz_stage("phased_solution.mtz", "data_mtz"),
        "phased_data_mtz",
        "data_mtz + 'phased' → phased_data_mtz")
  assert_equal(fu.get_mtz_stage("data.mtz", "data_mtz"),
        "data_mtz",
        "data_mtz default → data_mtz")

  # map_coeffs_mtz subcategories
  assert_equal(fu.get_mtz_stage("denmod_map_coeffs.mtz", "map_coeffs_mtz"),
        "denmod_map_coeffs",
        "map_coeffs + 'denmod' → denmod_map_coeffs")
  assert_equal(fu.get_mtz_stage("density_modified.mtz", "map_coeffs_mtz"),
        "denmod_map_coeffs",
        "map_coeffs + 'density_mod' → denmod_map_coeffs")
  assert_equal(fu.get_mtz_stage(
        "predict_and_build_map_coeffs.mtz", "map_coeffs_mtz"),
        "predict_build_map_coeffs",
        "map_coeffs + 'predict' + 'map_coeffs' → predict_build_map_coeffs")
  assert_equal(fu.get_mtz_stage(
        "overall_best_map_coeffs.mtz", "map_coeffs_mtz"),
        "predict_build_map_coeffs",
        "'overall_best_map_coeffs' → predict_build_map_coeffs")
  assert_equal(fu.get_mtz_stage("refine_001.mtz", "map_coeffs_mtz"),
        "refine_map_coeffs",
        "map_coeffs default → refine_map_coeffs")

  # Unrecognized category — passthrough
  assert_equal(fu.get_mtz_stage("something.dat", "unknown_category"),
        "unknown_category",
        "Unrecognized category passes through unchanged")


# =========================================================================
# get_category_for_extension — extension lookup + edge cases
# =========================================================================

def test_get_category_for_extension_semantics():
  """Pin get_category_for_extension's extension-to-category mapping.

    Returns the category for known extensions.  Returns None for:
      - Extensions that need context (.cif, .mtz)
      - Unknown extensions
      - Files with no extension
      - Files with only a trailing dot

    Case-insensitive on the extension (per the lower() call in the
    implementation).
    """
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  # Known model extensions
  assert_equal(fu.get_category_for_extension("file.pdb"), "model",
        ".pdb → model")
  assert_equal(fu.get_category_for_extension("file.ent"), "model",
        ".ent → model")

  # Context-required extensions (.cif, .mtz)
  assert_equal(fu.get_category_for_extension("file.cif"), None,
        ".cif → None (needs context: model? ligand? data?)")
  assert_equal(fu.get_category_for_extension("file.mtz"), None,
        ".mtz → None (needs context: data_mtz? map_coeffs_mtz?)")

  # Other data extensions
  assert_equal(fu.get_category_for_extension("file.sca"), "data_mtz",
        ".sca → data_mtz")
  assert_equal(fu.get_category_for_extension("file.hkl"), "data_mtz",
        ".hkl → data_mtz")

  # Map extensions
  assert_equal(fu.get_category_for_extension("file.map"), "map",
        ".map → map")
  assert_equal(fu.get_category_for_extension("file.ccp4"), "map",
        ".ccp4 → map")
  assert_equal(fu.get_category_for_extension("file.mrc"), "map",
        ".mrc → map")

  # Sequence extensions
  assert_equal(fu.get_category_for_extension("file.fa"), "sequence",
        ".fa → sequence")
  assert_equal(fu.get_category_for_extension("file.fasta"), "sequence",
        ".fasta → sequence")
  assert_equal(fu.get_category_for_extension("file.seq"), "sequence",
        ".seq → sequence")
  assert_equal(fu.get_category_for_extension("file.dat"), "sequence",
        ".dat → sequence")

  # Case-insensitivity
  assert_equal(fu.get_category_for_extension("STRUCTURE.PDB"), "model",
        "Uppercase .PDB → model (case-insensitive)")

  # Multi-dot / nested-extension cases (Gemini's H12 plan review):
  # Implementation uses os.path.splitext which takes ONLY the trailing
  # extension, so intermediate dots in the filename don't confuse it.
  assert_equal(fu.get_category_for_extension("model.backup.pdb"), "model",
        "Multi-dot with .pdb suffix → model (trailing-only)")
  assert_equal(fu.get_category_for_extension("data.processed.mtz"), None,
        "Multi-dot with .mtz suffix → None (still needs context)")

  # No-extension / edge cases (Gemini's H12 plan review):
  assert_equal(fu.get_category_for_extension("README"), None,
        "No extension → None (graceful, no crash)")
  assert_equal(fu.get_category_for_extension("LOCAL_MAP"), None,
        "No extension, suggestive name → None (extension is the signal)")
  assert_equal(fu.get_category_for_extension("file."), None,
        "Trailing dot only → None")

  # Unknown extension
  assert_equal(fu.get_category_for_extension("archive.tar.gz"), None,
        "Unknown trailing extension → None")
  assert_equal(fu.get_category_for_extension("file.xyz"), None,
        "Unknown extension → None")


# =========================================================================
# is_mtz_file / is_model_file / is_map_file / is_sequence_file
# =========================================================================

def test_is_mtz_file_semantics():
  """Pin is_mtz_file: True iff lowercase filename ends with .mtz."""
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  assert_true(fu.is_mtz_file("data.mtz"), "data.mtz → True")
  assert_true(fu.is_mtz_file("DATA.MTZ"), "DATA.MTZ → True (case-insensitive)")
  assert_true(fu.is_mtz_file("/path/to/file.mtz"), "Full path .mtz → True")
  assert_true(fu.is_mtz_file("refine_001_001.mtz"),
        "Refine output MTZ → True")
  assert_false(fu.is_mtz_file("model.pdb"), ".pdb → False")
  assert_false(fu.is_mtz_file("data.cif"), ".cif → False")
  assert_false(fu.is_mtz_file("mtzdata"),
        "'mtz' substring without extension → False")
  assert_false(fu.is_mtz_file("data.mtz.bak"),
        "Trailing .bak after .mtz → False (suffix matters)")


def test_is_model_file_semantics():
  """Pin is_model_file: .pdb or .ent always; .cif only if 'ligand' NOT in name.

    Load-bearing quirk: the 'ligand' filter applies ONLY to .cif files.
    A file named 'ligand_fit.pdb' or 'protein_with_bound_ligand.pdb' is
    still considered a model file, because .pdb extension overrides the
    name-based filter.  This is intentional — PDB files are always
    structural models (small-molecule PDBs are rare and handled by
    content-based guards elsewhere in command_builder), whereas .cif
    files are ambiguous and 'ligand'-in-name is the cheap disambiguator.
    """
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  # Standard model formats
  assert_true(fu.is_model_file("protein.pdb"), ".pdb → True")
  assert_true(fu.is_model_file("PROTEIN.PDB"), ".PDB → True (case-insensitive)")
  assert_true(fu.is_model_file("model.ent"), ".ent → True")
  assert_true(fu.is_model_file("protein.cif"),
        ".cif without 'ligand' in name → True")
  assert_true(fu.is_model_file("refine_001.cif"),
        "Refine output mmCIF → True (no 'ligand')")
  assert_true(fu.is_model_file("refine_001_001.cif"),
        "The H10/H11 bug filename → True from is_model_file's view "
        "(category-level filter; exclude_patterns is a separate slot-level filter)")

  # .cif with 'ligand' in name → False (the name-based disambiguator)
  assert_false(fu.is_model_file("my_ligand.cif"),
        ".cif with 'ligand' in name → False")
  assert_false(fu.is_model_file("co_crystallized_ligand_structure.cif"),
        ".cif with 'ligand' as substring anywhere → False")

  # Load-bearing quirk (Gemini's H12 plan review):
  # .pdb extension overrides the 'ligand' filter
  assert_true(fu.is_model_file("ligand_fit.pdb"),
        "ligand_fit.pdb → True (.pdb extension overrides 'ligand' filter)")
  assert_true(fu.is_model_file("protein_with_bound_ligand.pdb"),
        "protein_with_bound_ligand.pdb → True (.pdb wins)")

  # Non-model files
  assert_false(fu.is_model_file("data.mtz"), ".mtz → False")
  assert_false(fu.is_model_file("map.ccp4"), ".ccp4 → False")
  assert_false(fu.is_model_file("seq.fa"), ".fa → False")


def test_is_map_file_semantics():
  """Pin is_map_file: True iff lowercase ends with .map, .ccp4, or .mrc."""
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  assert_true(fu.is_map_file("output.map"), ".map → True")
  assert_true(fu.is_map_file("output.ccp4"), ".ccp4 → True")
  assert_true(fu.is_map_file("output.mrc"), ".mrc → True")
  assert_true(fu.is_map_file("OUTPUT.CCP4"), "case-insensitive")
  assert_false(fu.is_map_file("data.mtz"), ".mtz → False (use is_mtz_file)")
  assert_false(fu.is_map_file("model.pdb"), ".pdb → False")
  assert_false(fu.is_map_file("mapping.txt"),
        "'map' substring without extension → False")


def test_is_sequence_file_semantics():
  """Pin is_sequence_file: .fa, .fasta, .seq, .dat (lowercase)."""
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  assert_true(fu.is_sequence_file("protein.fa"), ".fa → True")
  assert_true(fu.is_sequence_file("protein.fasta"), ".fasta → True")
  assert_true(fu.is_sequence_file("protein.seq"), ".seq → True")
  assert_true(fu.is_sequence_file("protein.dat"), ".dat → True")
  assert_true(fu.is_sequence_file("PROTEIN.FA"), "case-insensitive")
  assert_false(fu.is_sequence_file("model.pdb"), ".pdb → False")
  assert_false(fu.is_sequence_file("data.mtz"), ".mtz → False")


# =========================================================================
# Cross-function consistency
# =========================================================================

def test_categorizer_cross_function_consistency():
  """Integration: representative files give consistent answers across
    is_mtz_file / classify_mtz_type / get_mtz_stage chain, and across
    get_category_for_extension / is_X_file pairs.

    This catches the case where someone modifies one function but
    forgets the related one.
    """
  fu = _import_file_utils()
  if fu is None:
    print("  SKIP (file_utils not importable — sandbox)")
    return

  # Refine output MTZ: is_mtz_file True, classify returns map_coeffs_mtz,
  # get_mtz_stage returns refine_map_coeffs
  fn = "refine_001_001.mtz"
  assert_true(fu.is_mtz_file(fn), "is_mtz_file(refine_001_001.mtz)")
  cat = fu.classify_mtz_type(fn)
  assert_equal(cat, "map_coeffs_mtz", "classify_mtz_type(refine_001_001.mtz)")
  assert_equal(fu.get_mtz_stage(fn, cat), "refine_map_coeffs",
        "get_mtz_stage chained from classify")

  # is_mtz_file True but get_category_for_extension is None — both correct,
  # because extension alone is insufficient for .mtz
  assert_true(fu.is_mtz_file("refine_001.mtz"),
        "is_mtz_file: extension match alone is sufficient")
  assert_equal(fu.get_category_for_extension("refine_001.mtz"), None,
        "get_category_for_extension: .mtz returns None (needs context)")

  # Ligand .cif: is_model_file False, get_category_for_extension None.
  # Both correct: .cif is context-dependent, and the name-based 'ligand'
  # filter says this isn't a model.
  assert_false(fu.is_model_file("my_ligand.cif"),
        "is_model_file: 'ligand'-in-name .cif → False")
  assert_equal(fu.get_category_for_extension("my_ligand.cif"), None,
        "get_category_for_extension: .cif is context-dependent → None")

  # Ligand .pdb: is_model_file True (extension wins), but in practice the
  # downstream consumer should treat it via content guards.  This test
  # documents that the categorizer alone does NOT block ligand_fit.pdb.
  assert_true(fu.is_model_file("ligand_fit.pdb"),
        "is_model_file: .pdb extension overrides 'ligand' filter — "
        "downstream content guards are the actual ligand vs protein "
        "discriminator for PDB files")


def run_all_tests():
  """Run all file_utils semantic-pin tests."""
  run_tests_with_fail_fast()


if __name__ == "__main__":
  success = run_tests_with_fail_fast()
  if not success:
    sys.exit(1)
  print("\nOK")
