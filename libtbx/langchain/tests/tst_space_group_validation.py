"""K_H14_ITEM_3: space_group validation extensions (v119.H14).

The pre-H14 validator at agent/directive_extractor.py:2608 had two
gaps when handling LLM-emitted space_group values:

  1. The _SYMMETRY_SENTINELS frozenset did not include the
     "Not explicitly mentioned" family — Tom's qwen2.5:72b
     xtriage tutorial run emitted "Not explicitly mentio"
     (truncated mid-phrase by LLM length cap) which slipped
     through.  The H13.1 ship's plan-side override prevented
     downstream harm, but the value was still in directives.

  2. The negative checks (sentinel + starts-with-letter +
     length<=25) let multi-word phrases through if they
     happened to start with a letter — e.g., "Solve the
     structure" passes all three checks.

H14 fix:
  a. Extend _SYMMETRY_SENTINELS with "not explicitly *"
     variants (including truncated forms).
  b. Add a positive Hermann-Mauguin shape check
     (_looks_like_space_group) and call it from the validator.

This file pins:
  - The new sentinels are recognised
  - _looks_like_space_group accepts standard space group symbols
  - _looks_like_space_group rejects prose phrases
  - validate_directives drops invalid values via either path
  - validate_directives preserves valid values
  - Existing sentinel behavior is unchanged (regression guard)

Total: 24 tests.
"""
from __future__ import absolute_import, division, print_function

import os
import sys


# =====================================================================
# Import helpers — sandbox-friendly
# =====================================================================

def _try_import_de():
    """Import the directive_extractor module."""
    try:
        from libtbx.langchain.agent import directive_extractor as de
        return de, None
    except ImportError:
        pass
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        if parent not in sys.path:
            sys.path.insert(0, parent)
        from agent import directive_extractor as de
        return de, None
    except ImportError as e:
        return None, str(e)


def _validate_sg(de, value):
    """Helper: run validate_directives on a single space_group value,
    return the surviving value (or sentinel 'DROPPED' if removed)."""
    directives = {
        "program_settings": {
            "phenix.refine": {"space_group": value}
        }
    }
    result = de.validate_directives(directives, log=lambda m: None)
    refine_settings = (
        result.get("program_settings", {}).get("phenix.refine", {}))
    return refine_settings.get("space_group", "DROPPED")


# =====================================================================
# §A: New sentinels (3)
# =====================================================================

def test_sentinel_not_explicitly_mentioned():
    """The full phrase Tom's LLMs emit."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert de._is_symmetry_sentinel("Not explicitly mentioned")
    assert de._is_symmetry_sentinel("not explicitly mentioned")  # case-insensitive
    assert de._is_symmetry_sentinel("  Not explicitly mentioned  ")  # stripped
    print("  PASS: test_sentinel_not_explicitly_mentioned")


def test_sentinel_truncated_forms():
    """LLM output sometimes hits length caps mid-phrase.  The
    truncated forms are explicitly added to the sentinel set."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # Tom's actual case from H13.1 verification
    assert de._is_symmetry_sentinel("Not explicitly mentio"), (
        "The truncated form Tom saw must be recognised")
    # Other plausible truncations
    assert de._is_symmetry_sentinel("Not explicitly state")
    assert de._is_symmetry_sentinel("Not explicitly giv")
    assert de._is_symmetry_sentinel("Not explicitly specif")
    print("  PASS: test_sentinel_truncated_forms")


def test_sentinel_other_explicit_variants():
    """Cover the "Not explicitly {stated,given,specified,defined,
    listed,noted}" family."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert de._is_symmetry_sentinel("Not explicitly stated")
    assert de._is_symmetry_sentinel("Not explicitly given")
    assert de._is_symmetry_sentinel("Not explicitly specified")
    assert de._is_symmetry_sentinel("Not explicitly defined")
    assert de._is_symmetry_sentinel("Not explicitly listed")
    assert de._is_symmetry_sentinel("Not explicitly noted")
    print("  PASS: test_sentinel_other_explicit_variants")


# =====================================================================
# §B: Existing sentinels still work (regression guard) (1)
# =====================================================================

def test_existing_sentinels_unchanged():
    """The H14 extension is ADDITIVE — pre-H14 sentinels must
    still be recognised."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # Pre-H14 sentinels — these must still work
    pre_h14 = ["none", "not mentioned", "not specified", "not given",
               "not provided", "not applicable", "unknown", "n/a",
               "null", "auto", "default"]
    for s in pre_h14:
        assert de._is_symmetry_sentinel(s), (
            "Pre-H14 sentinel %r must still be recognised" % s)
    print("  PASS: test_existing_sentinels_unchanged")


# =====================================================================
# §C: _looks_like_space_group positive shape check (4)
# =====================================================================

def test_hm_form_accepts_standard_symbols():
    """Common Hermann-Mauguin space group symbols must be accepted."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    for sg in ["P 1 21 1", "P21", "P212121", "P 21 21 21",
               "P-1", "P1", "F222", "I 4", "C 2 2 21", "I4"]:
        assert de._looks_like_space_group(sg), (
            "Standard HM symbol %r must be accepted" % sg)
    print("  PASS: test_hm_form_accepts_standard_symbols")


def test_hm_form_accepts_all_230_space_groups():
    """The 230 official Hermann-Mauguin space-group symbols
    (International Tables for Crystallography, Vol. A) must ALL
    match the H14 shape check.

    v119.H14 initially shipped a too-strict regex that rejected
    159/230 valid groups (any symbol containing m, c, n, d, a, b,
    or '/' after the lattice letter — i.e., all the monoclinic
    slash forms, orthorhombic mirror/glide groups, tetragonal
    mirror forms, all the cubic high-symmetry groups like Pm-3m,
    Im-3m, Fd-3m, etc.).  The H14 review fixed it to admit the
    full 230.  This test pins that fix.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # All 230 space groups, standard short Hermann-Mauguin symbols.
    # Source: International Tables for Crystallography Vol A.
    all_230 = [
        # Triclinic (1-2)
        "P1", "P-1",
        # Monoclinic (3-15)
        "P2", "P21", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m", "P21/m",
        "C2/m", "P2/c", "P21/c", "C2/c",
        # Orthorhombic (16-74)
        "P222", "P2221", "P21212", "P212121", "C2221", "C222",
        "F222", "I222", "I212121",
        "Pmm2", "Pmc21", "Pcc2", "Pma2", "Pca21", "Pnc2", "Pmn21",
        "Pba2", "Pna21", "Pnn2",
        "Cmm2", "Cmc21", "Ccc2", "Amm2", "Aem2", "Ama2", "Aea2",
        "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2",
        "Pmmm", "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna",
        "Pcca", "Pbam", "Pccn", "Pbcm", "Pnnm", "Pmmn", "Pbcn",
        "Pbca", "Pnma",
        "Cmcm", "Cmce", "Cmmm", "Cccm", "Cmme", "Ccce",
        "Fmmm", "Fddd", "Immm", "Ibam", "Ibca", "Imma",
        # Tetragonal (75-142)
        "P4", "P41", "P42", "P43", "I4", "I41",
        "P-4", "I-4",
        "P4/m", "P42/m", "P4/n", "P42/n", "I4/m", "I41/a",
        "P422", "P4212", "P4122", "P41212", "P4222", "P42212",
        "P4322", "P43212", "I422", "I4122",
        "P4mm", "P4bm", "P42cm", "P42nm", "P4cc", "P4nc", "P42mc",
        "P42bc", "I4mm", "I4cm", "I41md", "I41cd",
        "P-42m", "P-42c", "P-421m", "P-421c", "P-4m2", "P-4c2",
        "P-4b2", "P-4n2", "I-4m2", "I-4c2", "I-42m", "I-42d",
        "P4/mmm", "P4/mcc", "P4/nbm", "P4/nnc", "P4/mbm", "P4/mnc",
        "P4/nmm", "P4/ncc", "P42/mmc", "P42/mcm", "P42/nbc",
        "P42/nnm", "P42/mbc", "P42/mnm", "P42/nmc", "P42/ncm",
        "I4/mmm", "I4/mcm", "I41/amd", "I41/acd",
        # Trigonal (143-167)
        "P3", "P31", "P32", "R3",
        "P-3", "R-3",
        "P312", "P321", "P3112", "P3121", "P3212", "P3221", "R32",
        "P3m1", "P31m", "P3c1", "P31c", "R3m", "R3c",
        "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c",
        # Hexagonal (168-194)
        "P6", "P61", "P65", "P62", "P64", "P63",
        "P-6",
        "P6/m", "P63/m",
        "P622", "P6122", "P6522", "P6222", "P6422", "P6322",
        "P6mm", "P6cc", "P63cm", "P63mc",
        "P-6m2", "P-6c2", "P-62m", "P-62c",
        "P6/mmm", "P6/mcc", "P63/mcm", "P63/mmc",
        # Cubic (195-230)
        "P23", "F23", "I23", "P213", "I213",
        "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3", "Ia-3",
        "P432", "P4232", "F432", "F4132", "I432", "P4332",
        "P4132", "I4132",
        "P-43m", "F-43m", "I-43m", "P-43n", "F-43c", "I-43d",
        "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m",
        "Fm-3m", "Fm-3c", "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d",
    ]
    assert len(all_230) == 230, (
        "Expected exactly 230 space groups, got %d" % len(all_230))

    failures = [sg for sg in all_230
                if not de._looks_like_space_group(sg)]
    assert not failures, (
        "All 230 official Hermann-Mauguin symbols must be accepted "
        "by _looks_like_space_group.  Rejected: %r" % failures)
    print("  PASS: test_hm_form_accepts_all_230_space_groups "
          "(230/230)")


def test_hm_form_accepts_spaced_variants():
    """The whitespace-separated and PERCEIVE-emitted forms must
    also be accepted (e.g., 'P 1 21 1', 'P 21 21 21 (No. 19)')."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    spaced = [
        "P 1 21 1", "P 21 21 21", "P 32 2 1", "C 2 2 21",
        "P 41 21 2", "I 41 22", "P 6 2 2", "R 3 2",
        "P 21/c", "P 21/n",
    ]
    for sg in spaced:
        assert de._looks_like_space_group(sg), (
            "Spaced HM form %r must be accepted" % sg)
    print("  PASS: test_hm_form_accepts_spaced_variants")


def test_hm_form_accepts_with_parenthetical_number():
    """Forms like 'P 1 21 1 (No. 4)' (with the IUCr number annotation)
    must be accepted."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert de._looks_like_space_group("P 1 21 1 (No. 4)")
    assert de._looks_like_space_group("P21 (No. 4)")
    assert de._looks_like_space_group("P 21 21 21 (No. 19)")
    print("  PASS: test_hm_form_accepts_with_parenthetical_number")


def test_hm_form_rejects_prose():
    """Multi-word English phrases must be rejected.  This catches
    the pre-H14 gap where phrases starting with a letter passed
    all the negative checks."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    for phrase in ["Solve the structure", "garbage value 42",
                   "Not explicitly mentioned",
                   "From data file", "see header"]:
        assert not de._looks_like_space_group(phrase), (
            "Prose %r must be rejected by HM shape check" % phrase)
    print("  PASS: test_hm_form_rejects_prose")


def test_hm_form_rejects_empty_and_none():
    """Edge cases: empty string, None, non-strings."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert not de._looks_like_space_group("")
    assert not de._looks_like_space_group(None)
    assert not de._looks_like_space_group(42)  # int — not string
    assert not de._looks_like_space_group([])
    print("  PASS: test_hm_form_rejects_empty_and_none")


# =====================================================================
# §D: End-to-end validate_directives behavior (5)
# =====================================================================

def test_e2e_drops_not_explicitly_mentioned():
    """The Tom-case value must be dropped by validate_directives."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert _validate_sg(de, "Not explicitly mentioned") == "DROPPED"
    assert _validate_sg(de, "Not explicitly mentio") == "DROPPED"
    print("  PASS: test_e2e_drops_not_explicitly_mentioned")


def test_e2e_drops_prose_phrases():
    """Phrases that start with a letter but aren't space groups
    must be dropped."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert _validate_sg(de, "Solve the structure") == "DROPPED"
    assert _validate_sg(de, "see header") == "DROPPED"
    print("  PASS: test_e2e_drops_prose_phrases")


def test_e2e_keeps_valid_hermann_mauguin():
    """Valid HM symbols must survive validation unchanged."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    for sg in ["P 1 21 1", "P21", "P212121", "P-1", "F222", "C 2 2 21"]:
        result = _validate_sg(de, sg)
        assert result == sg, (
            "Valid space group %r must survive validation.  Got: %r"
            % (sg, result))
    print("  PASS: test_e2e_keeps_valid_hermann_mauguin")


def test_e2e_keeps_hm_with_parenthetical():
    """The PERCEIVE-emitted form with (No. N) annotation must
    survive."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    assert _validate_sg(de, "P 1 21 1 (No. 4)") == "P 1 21 1 (No. 4)"
    print("  PASS: test_e2e_keeps_hm_with_parenthetical")


def test_e2e_other_settings_preserved():
    """When space_group is dropped, OTHER settings in the same
    program_settings dict must be preserved."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    directives = {
        "program_settings": {
            "phenix.refine": {
                "space_group": "Not explicitly mentioned",
                "resolution": "2.3",
                "atom_type": "Se",
            }
        }
    }
    result = de.validate_directives(directives, log=lambda m: None)
    refine = result.get("program_settings", {}).get("phenix.refine", {})
    assert "space_group" not in refine, "space_group should be dropped"
    assert refine.get("resolution") == 2.3, (
        "resolution should be preserved (and coerced to float)")
    assert refine.get("atom_type") == "Se", (
        "atom_type should be preserved")
    print("  PASS: test_e2e_other_settings_preserved")


# =====================================================================
# §E: Equivalence-class coverage (4)
#
# Per Tom (H14 review): "P 1 21 1" and P21 and p21 are all the same
# space group (#4).  The validator must accept all common variants.
# =====================================================================

def test_equivalence_case_insensitive():
    """Lattice letter case doesn't matter: p21 ≡ P21 ≡ P 21.
    PHENIX accepts both; the validator must too."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # All these are space group #4
    variants = ["P21", "p21", "P 21", "p 21", "P21 ", "p 21 "]
    for v in variants:
        assert de._looks_like_space_group(v.strip()), (
            "Case/whitespace variant %r of P21 must be accepted" % v)
        # And end-to-end through validate_directives
        result = _validate_sg(de, v)
        assert result != "DROPPED", (
            "Case/whitespace variant %r of P21 should NOT be dropped. "
            "Got: %r" % (v, result))
    # Lowercase lattice letters across the alphabet
    lowercase_lattices = ["p1", "p-1", "c2", "i4", "f222", "r3",
                          "p 21 21 21", "i 41 22", "c 2 2 21",
                          "p 32 2 1"]
    for v in lowercase_lattices:
        assert de._looks_like_space_group(v), (
            "Lowercase lattice variant %r must be accepted" % v)
    print("  PASS: test_equivalence_case_insensitive")


def test_equivalence_spacing_variants():
    """Whitespace between symmetry elements doesn't matter.
    P21 ≡ P 21; P212121 ≡ P 21 21 21 ≡ P  21  21  21."""
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # P 21 21 21 (#19) spacing variants
    p212121_variants = [
        "P212121",
        "P 21 21 21",
        "P  21  21  21",   # double spaces
        "P21 21 21",       # mixed: no space after first, space after others
        "P 212121",        # space after lattice only
    ]
    for v in p212121_variants:
        assert de._looks_like_space_group(v), (
            "P212121 spacing variant %r must be accepted" % v)
        result = _validate_sg(de, v)
        assert result != "DROPPED", (
            "P212121 variant %r should NOT be dropped. "
            "Got: %r" % (v, result))
    # Slash forms with spacing variants
    p21c_variants = [
        "P21/c",
        "P 21/c",
        "P21 /c",          # space before slash
        "P21/ c",          # space after slash
        "p21/c",           # lowercase
        "P 21 / c",        # spaces around slash
    ]
    for v in p21c_variants:
        assert de._looks_like_space_group(v), (
            "P21/c spacing variant %r must be accepted" % v)
    print("  PASS: test_equivalence_spacing_variants")


def test_equivalence_one_placeholder_forms():
    """"P 1 21 1" (full monoclinic form with "1" placeholders for the
    a-axis and c-axis trivial symmetry) is the same space group as
    P21 (short form).  Same for orthorhombic etc.

    The validator must accept both forms — it doesn't normalize, just
    checks the shape looks like a space group.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # Space group #4 (monoclinic): P21 / P 21 / P 1 21 1
    sg4_variants = ["P21", "P 21", "P 1 21 1", "p 1 21 1"]
    for v in sg4_variants:
        assert de._looks_like_space_group(v), (
            "Space group #4 variant %r must be accepted" % v)
    # Space group #145 (trigonal): P32 / P 32 / P 3 2 (with trailing 1s)
    sg145_variants = ["P32", "P 32"]
    for v in sg145_variants:
        assert de._looks_like_space_group(v), (
            "Space group #145 variant %r must be accepted" % v)
    # Space group #4 with parenthetical number — both compact and full
    annotated = ["P21 (No. 4)", "P 1 21 1 (No. 4)",
                 "p21 (no. 4)", "P21 (no 4)"]
    for v in annotated:
        assert de._looks_like_space_group(v), (
            "Annotated variant %r must be accepted" % v)
    print("  PASS: test_equivalence_one_placeholder_forms")


def test_equivalence_setting_variants_for_sg14():
    """Space group #14 has multiple legitimate settings depending on
    cell choice: P21/c (b-unique standard), P21/n (alternative),
    P21/a (older convention).  All are valid Hermann-Mauguin symbols.
    The validator accepts all without trying to canonicalize.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    sg14_settings = [
        "P21/c",       # standard b-unique setting
        "P21/n",       # alternative setting
        "P21/a",       # older convention
        "P 1 21/c 1",  # full b-unique form
        "P 1 21/n 1",  # full alternative form
        "P 1 21/a 1",  # full older form
    ]
    for v in sg14_settings:
        assert de._looks_like_space_group(v), (
            "Space group #14 setting %r must be accepted" % v)
    print("  PASS: test_equivalence_setting_variants_for_sg14")


# =====================================================================
# §F: Alternative cell/origin settings via colon suffix (per Gemini)
#
# International Tables Vol A uses colon-prefixed suffixes to denote
# alternative origin choices and cell settings.  These appear in real
# PDB/mmCIF metadata and cctbx tool output, so the validator must
# accept them.  Pre-Gemini-review H14 rejected all of these.
# =====================================================================

def test_setting_rhombohedral_axes_specifier():
    """Rhombohedral space groups (146, 148, 155, 160, 161, 166, 167)
    have two equivalent settings:
      :H — hexagonal axes (default in PHENIX, PDB)
      :R — rhombohedral axes
    Software needs to know which axes the user chose.  Both forms
    must be accepted by the validator.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    rhombohedral = [
        "R3:H", "R3:R", "R3:h", "R3:r",       # compact forms
        "R 3 :H", "R 3 :R",                    # spaced forms
        "R32:H", "R32:R",                      # space group 155
        "R-3:H", "R-3:R",                      # space group 148
        "R3m:H", "R3m:R",                      # space group 160
        "R-3m:H", "R-3m:R",                    # space group 166
        "R3c:H",                               # space group 161
        "R-3c:H",                              # space group 167
    ]
    for v in rhombohedral:
        assert de._looks_like_space_group(v), (
            "Rhombohedral axis-spec %r must be accepted" % v)
    print("  PASS: test_setting_rhombohedral_axes_specifier")


def test_setting_origin_choice_specifier():
    """Centrosymmetric groups with two origin choices (mostly cubic
    and some tetragonal/orthorhombic) use :1 vs :2 suffix.
      :1 — origin on a higher-symmetry point
      :2 — origin on an inversion center
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    origin_choice = [
        # Common cases
        "P 4/n :1", "P 4/n :2",
        "P4/n:1", "P4/n:2",
        "Pnnn:1", "Pnnn:2",
        "Pccn:1", "Pccn:2",
        # Cubic
        "Fd-3m:1", "Fd-3m:2",
        "Pn-3:1", "Pn-3:2",
        "Pn-3m:1", "Pn-3m:2",
        "Fd-3:1", "Fd-3:2",
        "Fd-3c:1", "Fd-3c:2",
    ]
    for v in origin_choice:
        assert de._looks_like_space_group(v), (
            "Origin-choice spec %r must be accepted" % v)
    print("  PASS: test_setting_origin_choice_specifier")


def test_setting_unique_axis_specifier():
    """Monoclinic groups can be specified with a unique-axis suffix:
      :a, :b, :c — which axis carries the 2-fold (or screw)
    PHENIX default is :b but :c is also common.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    unique_axis = [
        "P21/c:b", "P21/c:c", "P21/c:a",
        "P2/c:b", "P2/c:c",
        "P21:b", "P21:c",
        "C2/c:b", "C2/c:c",
    ]
    for v in unique_axis:
        assert de._looks_like_space_group(v), (
            "Unique-axis spec %r must be accepted" % v)
    print("  PASS: test_setting_unique_axis_specifier")


# =====================================================================
# §G: Known limitation — Risk B from Gemini review
#
# The regex is permissive within the HM alphabet.  Short English words
# that happen to start with a lattice letter (P, F, I, C, R, H, A, B)
# and use only alphabet characters (m, c, n, d, a, b, e, h, r, digits,
# space, /, _, -, :) will inadvertently match.
#
# This test documents the limitation rather than enforcing it as a
# spec — the limitation is acceptable because: (1) no realistic LLM
# emits these as space_group values, (2) the pre-H14 negative checks
# already let them through, (3) closing this gap fully would require
# enumerating the 230 valid HM symbols, which is beyond the scope of
# a directive sanity check.
#
# If production surfaces such hallucinations, the right response is
# probably to enumerate the 230 symbols rather than chase regex
# tightening.
# =====================================================================

def test_hm_form_known_limitation_short_words():
    """DOCUMENTED LIMITATION (Gemini H14 review, Risk B):

    Short English words starting with a Bravais lattice letter and
    using only HM-alphabet characters DO match the regex.  This
    test pins the limitation so future contributors don't think
    "Panda accepts" is a bug.

    The acceptable upstream filters are:
      1. The sentinel set (_SYMMETRY_SENTINELS), which catches
         common LLM placeholder phrases.
      2. The LLM's own output discipline — it would not emit
         'Panda' as a space group value.

    If production batches surface short-word hallucinations as
    space_group values, the right next step is to enumerate the
    230 official symbols (e.g., via cctbx.sgtbx) rather than
    tighten the regex.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return

    # These short English words DO match (limitation, not bug).
    # If a future tightening of the regex drops them, that's also
    # fine — but it must NOT drop any of the 230 official symbols.
    leak_examples = [
        "Panda",   # P + anda
        "Fame",    # F + ame
        "Bad",     # B + ad
        "Cab",     # C + ab
        "Bed",     # B + ed
        "Acme",    # A + cme
        "Pad",     # P + ad
        "Bee",     # B + ee
        "Ham",     # H + am
    ]
    leaks = [w for w in leak_examples if de._looks_like_space_group(w)]
    # We expect MOST to leak; this test pins the behavior rather than
    # asserts the limit goes away.  The assertion is that the
    # limitation is documented and recognised.
    print("  PASS: test_hm_form_known_limitation_short_words "
          "(%d/%d documented leaks confirmed)" %
          (len(leaks), len(leak_examples)))


def test_hm_form_rejects_words_with_non_alphabet_chars():
    """Words containing characters OUTSIDE the HM alphabet must be
    rejected even if they start with a lattice letter.

    The HM alphabet (post-Gemini-review) is:
      0-9 m c n d a b e h r [space] / _ - :

    Letters NOT in the alphabet: f, g, i, j, k, l, o, p, q, s, t,
    u, v, w, x, y, z.  Words containing any of these letters after
    the lattice letter must be rejected.
    """
    de, err = _try_import_de()
    if de is None:
        print("  SKIP: cannot import directive_extractor (%s)" % err)
        return
    # These should ALL be rejected because they contain at least one
    # non-alphabet character after the lattice letter.
    must_reject = [
        "Phaser",       # P + h_a_s_e_r — 's' not in alphabet (but h/r are now)
        "Pizza",        # P + i — 'i' not in alphabet
        "Cake",         # C + a_k_e — 'k' not in alphabet
        "From data",    # F + r_o_m — 'o' not in alphabet
        "Process",      # P + r_o — 'o' not in alphabet
        "Place",        # P + l — 'l' not in alphabet
        "Solve",        # S — not a lattice letter
        "Identify",     # I + d_e_n_t — 't' not in alphabet
        "ratio",        # r/R + a_t — 't' not in alphabet
        "see PDB",      # s — not a lattice letter
        "above",        # a — not a lattice letter
    ]
    failures = [w for w in must_reject if de._looks_like_space_group(w)]
    assert not failures, (
        "These words contain non-alphabet characters and must be "
        "rejected: %r" % failures)
    # Special: "Phaser" — after Gemini fix, h IS in alphabet
    # so trace: P + h_a_s_e_r → 's' is not in alphabet → reject. Good.
    print("  PASS: test_hm_form_rejects_words_with_non_alphabet_chars")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: new sentinels (3)
    test_sentinel_not_explicitly_mentioned()
    test_sentinel_truncated_forms()
    test_sentinel_other_explicit_variants()
    # §B: existing sentinels regression (1)
    test_existing_sentinels_unchanged()
    # §C: HM shape check (6)
    test_hm_form_accepts_standard_symbols()
    test_hm_form_accepts_all_230_space_groups()
    test_hm_form_accepts_spaced_variants()
    test_hm_form_accepts_with_parenthetical_number()
    test_hm_form_rejects_prose()
    test_hm_form_rejects_empty_and_none()
    # §D: end-to-end (5)
    test_e2e_drops_not_explicitly_mentioned()
    test_e2e_drops_prose_phrases()
    test_e2e_keeps_valid_hermann_mauguin()
    test_e2e_keeps_hm_with_parenthetical()
    test_e2e_other_settings_preserved()
    # §E: equivalence-class coverage (4)
    test_equivalence_case_insensitive()
    test_equivalence_spacing_variants()
    test_equivalence_one_placeholder_forms()
    test_equivalence_setting_variants_for_sg14()
    # §F: alternative cell/origin settings (3 — Gemini review)
    test_setting_rhombohedral_axes_specifier()
    test_setting_origin_choice_specifier()
    test_setting_unique_axis_specifier()
    # §G: limitation documentation (2 — Gemini review)
    test_hm_form_known_limitation_short_words()
    test_hm_form_rejects_words_with_non_alphabet_chars()


if __name__ == "__main__":
    print("K_H14_ITEM_3: space_group validation extensions (v119.H14)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H14_ITEM_3 complete.")
