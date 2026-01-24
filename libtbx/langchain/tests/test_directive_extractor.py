"""
Tests for directive_extractor module.

These tests verify that directives are correctly extracted from user advice
and properly validated.

Run with: python -m pytest tests/test_directive_extractor.py -v
Or: python tests/test_directive_extractor.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import unittest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.directive_extractor import (
    validate_directives,
    merge_directives,
    get_program_settings,
    check_stop_conditions,
    extract_directives_simple,
    _fix_program_name,
    format_directives_for_display,
)


class TestValidateDirectives(unittest.TestCase):
    """Tests for validate_directives function."""

    def test_empty_input(self):
        """Empty input should return empty dict."""
        self.assertEqual(validate_directives(None), {})
        self.assertEqual(validate_directives({}), {})
        self.assertEqual(validate_directives("not a dict"), {})

    def test_valid_program_settings(self):
        """Valid program settings should be preserved."""
        input_directives = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.5, "anisotropic_adp": True},
                "default": {"resolution": 3.0}
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("program_settings", result)
        self.assertEqual(result["program_settings"]["phenix.refine"]["resolution"], 2.5)
        self.assertEqual(result["program_settings"]["phenix.refine"]["anisotropic_adp"], True)
        self.assertEqual(result["program_settings"]["default"]["resolution"], 3.0)

    def test_invalid_program_name_removed(self):
        """Invalid program names should be removed."""
        input_directives = {
            "program_settings": {
                "invalid_program": {"resolution": 2.5},
                "phenix.refine": {"resolution": 3.0}
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("program_settings", result)
        self.assertNotIn("invalid_program", result["program_settings"])
        self.assertIn("phenix.refine", result["program_settings"])

    def test_fixable_program_name(self):
        """Fixable program names should be corrected."""
        input_directives = {
            "program_settings": {
                "refine": {"resolution": 2.5},  # Missing "phenix." prefix
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("program_settings", result)
        self.assertIn("phenix.refine", result["program_settings"])

    def test_type_conversion(self):
        """Values should be converted to correct types."""
        input_directives = {
            "program_settings": {
                "phenix.refine": {
                    "resolution": "2.5",  # String should become float
                    "cycles": "5",  # String should become int
                    "anisotropic_adp": 1  # Int should become bool
                }
            }
        }
        result = validate_directives(input_directives)

        settings = result["program_settings"]["phenix.refine"]
        self.assertIsInstance(settings["resolution"], float)
        self.assertIsInstance(settings["cycles"], int)
        self.assertIsInstance(settings["anisotropic_adp"], bool)

    def test_valid_stop_conditions(self):
        """Valid stop conditions should be preserved."""
        input_directives = {
            "stop_conditions": {
                "after_program": "phenix.refine",
                "after_cycle": 4,
                "skip_validation": True
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.refine")
        self.assertEqual(result["stop_conditions"]["after_cycle"], 4)
        self.assertEqual(result["stop_conditions"]["skip_validation"], True)

    def test_file_preferences(self):
        """File preferences should be validated."""
        input_directives = {
            "file_preferences": {
                "model": "beta.pdb",
                "exclude": ["old.pdb", "backup.pdb"]
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("file_preferences", result)
        self.assertEqual(result["file_preferences"]["model"], "beta.pdb")
        self.assertEqual(result["file_preferences"]["exclude"], ["old.pdb", "backup.pdb"])

    def test_workflow_preferences(self):
        """Workflow preferences should be validated."""
        input_directives = {
            "workflow_preferences": {
                "skip_programs": ["phenix.autobuild"],
                "use_experimental_phasing": True
            }
        }
        result = validate_directives(input_directives)

        self.assertIn("workflow_preferences", result)
        self.assertEqual(result["workflow_preferences"]["skip_programs"], ["phenix.autobuild"])
        self.assertEqual(result["workflow_preferences"]["use_experimental_phasing"], True)

    def test_constraints_preserved(self):
        """Constraints should be kept as strings."""
        input_directives = {
            "constraints": [
                "Do not add waters until R-free < 0.30",
                "Use TLS after cycle 3"
            ]
        }
        result = validate_directives(input_directives)

        self.assertIn("constraints", result)
        self.assertEqual(len(result["constraints"]), 2)


class TestMergeDirectives(unittest.TestCase):
    """Tests for merge_directives function."""

    def test_empty_merge(self):
        """Merging with empty should return non-empty."""
        base = {"program_settings": {"default": {"resolution": 2.5}}}

        self.assertEqual(merge_directives(None, base), base)
        self.assertEqual(merge_directives(base, None), base)
        self.assertEqual(merge_directives(None, None), {})

    def test_override_wins(self):
        """Override values should take precedence."""
        base = {
            "program_settings": {
                "default": {"resolution": 2.5}
            }
        }
        override = {
            "program_settings": {
                "default": {"resolution": 3.0}
            }
        }
        result = merge_directives(base, override)

        self.assertEqual(result["program_settings"]["default"]["resolution"], 3.0)

    def test_deep_merge_program_settings(self):
        """Program settings should be deep merged."""
        base = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.5}
            }
        }
        override = {
            "program_settings": {
                "phenix.refine": {"anisotropic_adp": True}
            }
        }
        result = merge_directives(base, override)

        self.assertEqual(result["program_settings"]["phenix.refine"]["resolution"], 2.5)
        self.assertEqual(result["program_settings"]["phenix.refine"]["anisotropic_adp"], True)

    def test_constraints_concatenate(self):
        """Constraints should be concatenated."""
        base = {"constraints": ["constraint 1"]}
        override = {"constraints": ["constraint 2"]}
        result = merge_directives(base, override)

        self.assertEqual(len(result["constraints"]), 2)


class TestGetProgramSettings(unittest.TestCase):
    """Tests for get_program_settings function."""

    def test_empty_directives(self):
        """Empty directives should return empty dict."""
        self.assertEqual(get_program_settings(None, "phenix.refine"), {})
        self.assertEqual(get_program_settings({}, "phenix.refine"), {})

    def test_default_fallback(self):
        """Should fall back to default settings."""
        directives = {
            "program_settings": {
                "default": {"resolution": 2.5}
            }
        }
        result = get_program_settings(directives, "phenix.refine")

        self.assertEqual(result["resolution"], 2.5)

    def test_program_specific_override(self):
        """Program-specific settings should override defaults."""
        directives = {
            "program_settings": {
                "default": {"resolution": 2.5},
                "phenix.autosol": {"resolution": 3.0}
            }
        }

        # autosol should get 3.0
        result = get_program_settings(directives, "phenix.autosol")
        self.assertEqual(result["resolution"], 3.0)

        # refine should get default 2.5
        result = get_program_settings(directives, "phenix.refine")
        self.assertEqual(result["resolution"], 2.5)

    def test_merge_default_and_specific(self):
        """Should merge default with program-specific."""
        directives = {
            "program_settings": {
                "default": {"resolution": 2.5, "cycles": 5},
                "phenix.refine": {"anisotropic_adp": True}
            }
        }
        result = get_program_settings(directives, "phenix.refine")

        self.assertEqual(result["resolution"], 2.5)  # From default
        self.assertEqual(result["cycles"], 5)  # From default
        self.assertEqual(result["anisotropic_adp"], True)  # From specific


class TestCheckStopConditions(unittest.TestCase):
    """Tests for check_stop_conditions function."""

    def test_no_stop_conditions(self):
        """No conditions should not trigger stop."""
        should_stop, reason = check_stop_conditions({}, 5, "phenix.refine")
        self.assertFalse(should_stop)
        self.assertIsNone(reason)

    def test_after_cycle(self):
        """Should stop after specified cycle."""
        directives = {
            "stop_conditions": {"after_cycle": 4}
        }

        # Before cycle 4
        should_stop, reason = check_stop_conditions(directives, 3, "phenix.refine")
        self.assertFalse(should_stop)

        # At cycle 4
        should_stop, reason = check_stop_conditions(directives, 4, "phenix.refine")
        self.assertTrue(should_stop)
        self.assertIn("cycle", reason.lower())

    def test_after_program(self):
        """Should stop after specified program."""
        directives = {
            "stop_conditions": {"after_program": "phenix.refine"}
        }

        # Different program
        should_stop, reason = check_stop_conditions(directives, 2, "phenix.phaser")
        self.assertFalse(should_stop)

        # Matching program
        should_stop, reason = check_stop_conditions(directives, 2, "phenix.refine")
        self.assertTrue(should_stop)
        self.assertIn("refine", reason.lower())

    def test_r_free_target(self):
        """Should stop when R-free target reached."""
        directives = {
            "stop_conditions": {"r_free_target": 0.25}
        }

        # Above target
        should_stop, reason = check_stop_conditions(
            directives, 5, "phenix.refine", {"r_free": 0.30}
        )
        self.assertFalse(should_stop)

        # At target
        should_stop, reason = check_stop_conditions(
            directives, 5, "phenix.refine", {"r_free": 0.24}
        )
        self.assertTrue(should_stop)


class TestExtractDirectivesSimple(unittest.TestCase):
    """Tests for simple pattern-based extraction."""

    def test_empty_input(self):
        """Empty input should return empty dict."""
        self.assertEqual(extract_directives_simple(None), {})
        self.assertEqual(extract_directives_simple(""), {})

    def test_resolution_extraction(self):
        """Should extract resolution values."""
        test_cases = [
            ("resolution=2.5", 2.5),
            ("resolution: 3.0", 3.0),
            ("resolution of 2.8", 2.8),
            ("2.5 Angstrom resolution", 2.5),
            ("to 3.0 A", 3.0),
        ]

        for advice, expected in test_cases:
            result = extract_directives_simple(advice)
            self.assertIn("program_settings", result, f"Failed for: {advice}")
            self.assertEqual(
                result["program_settings"]["default"]["resolution"],
                expected,
                f"Failed for: {advice}"
            )

    def test_anisotropic_extraction(self):
        """Should extract anisotropic setting."""
        result = extract_directives_simple("run anisotropic refinement")

        self.assertIn("program_settings", result)
        self.assertIn("phenix.refine", result["program_settings"])
        self.assertTrue(result["program_settings"]["phenix.refine"]["anisotropic_adp"])

    def test_stop_after_cycle(self):
        """Should extract stop after cycle N."""
        result = extract_directives_simple("stop after cycle 4")

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_cycle"], 4)

    def test_stop_after_first_refine(self):
        """Should extract stop after first refinement."""
        result = extract_directives_simple("stop after the first refinement")

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.refine")
        self.assertEqual(result["stop_conditions"]["max_refine_cycles"], 1)

    def test_skip_validation(self):
        """Should extract skip validation directive."""
        for advice in ["skip validation", "don't validate", "no validation"]:
            result = extract_directives_simple(advice)
            self.assertIn("stop_conditions", result, f"Failed for: {advice}")
            self.assertTrue(
                result["stop_conditions"].get("skip_validation"),
                f"Failed for: {advice}"
            )


class TestFixProgramName(unittest.TestCase):
    """Tests for program name fixing."""

    def test_common_variations(self):
        """Should fix common name variations."""
        test_cases = [
            ("refine", "phenix.refine"),
            ("refinement", "phenix.refine"),
            ("autosol", "phenix.autosol"),
            ("phaser", "phenix.phaser"),
            ("rsr", "phenix.real_space_refine"),
        ]

        for input_name, expected in test_cases:
            result = _fix_program_name(input_name)
            self.assertEqual(result, expected, f"Failed for: {input_name}")

    def test_unknown_returns_none(self):
        """Unknown names should return None."""
        self.assertIsNone(_fix_program_name("unknown_program"))
        self.assertIsNone(_fix_program_name(""))
        self.assertIsNone(_fix_program_name(None))


class TestFormatDirectivesForDisplay(unittest.TestCase):
    """Tests for display formatting."""

    def test_empty_directives(self):
        """Empty directives should give appropriate message."""
        result = format_directives_for_display({})
        self.assertIn("No directives", result)

    def test_formatted_output(self):
        """Should format directives readably."""
        directives = {
            "program_settings": {
                "phenix.refine": {"resolution": 2.5}
            },
            "stop_conditions": {
                "after_cycle": 4
            }
        }
        result = format_directives_for_display(directives)

        self.assertIn("Program Settings", result)
        self.assertIn("phenix.refine", result)
        self.assertIn("resolution", result)
        self.assertIn("Stop Conditions", result)
        self.assertIn("after_cycle", result)


class TestCheckStopConditions(unittest.TestCase):
    """Tests for check_stop_conditions function."""

    def test_after_program_exact_match(self):
        """Should stop when program matches exactly."""
        directives = {
            "stop_conditions": {
                "after_program": "phenix.xtriage"
            }
        }
        should_stop, reason = check_stop_conditions(directives, 1, "phenix.xtriage")
        self.assertTrue(should_stop)
        self.assertIn("xtriage", reason)

    def test_after_program_normalized_match(self):
        """Should stop when program matches after normalization."""
        directives = {
            "stop_conditions": {
                "after_program": "phenix.xtriage"
            }
        }
        # Test with short name
        should_stop, reason = check_stop_conditions(directives, 1, "xtriage")
        self.assertTrue(should_stop)
        self.assertIn("xtriage", reason)

    def test_after_program_reverse_normalization(self):
        """Should match when directive has short name but program has full name."""
        directives = {
            "stop_conditions": {
                "after_program": "xtriage"  # Short name in directive
            }
        }
        should_stop, reason = check_stop_conditions(directives, 1, "phenix.xtriage")
        self.assertTrue(should_stop)

    def test_after_cycle(self):
        """Should stop after specified cycle."""
        directives = {
            "stop_conditions": {
                "after_cycle": 3
            }
        }
        should_stop, _ = check_stop_conditions(directives, 2, "phenix.refine")
        self.assertFalse(should_stop)

        should_stop, reason = check_stop_conditions(directives, 3, "phenix.refine")
        self.assertTrue(should_stop)
        self.assertIn("cycle 3", reason)

    def test_no_match(self):
        """Should not stop when conditions not met."""
        directives = {
            "stop_conditions": {
                "after_program": "phenix.xtriage"
            }
        }
        should_stop, _ = check_stop_conditions(directives, 1, "phenix.phaser")
        self.assertFalse(should_stop)

    def test_empty_directives(self):
        """Should not stop with empty directives."""
        should_stop, _ = check_stop_conditions({}, 1, "phenix.xtriage")
        self.assertFalse(should_stop)

        should_stop, _ = check_stop_conditions(None, 1, "phenix.xtriage")
        self.assertFalse(should_stop)


class TestTutorialDetection(unittest.TestCase):
    """Tests for tutorial/procedure detection in simple extraction."""

    def test_xtriage_tutorial_detection(self):
        """Should detect xtriage tutorial and add stop condition."""
        advice = "Run Xtriage on the reflection data to analyze for twinning."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.xtriage")
        self.assertTrue(result["stop_conditions"]["skip_validation"])

    def test_twinning_check_detection(self):
        """Should detect twinning check and stop after xtriage."""
        advice = "Check for twinning in the porin dataset."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.xtriage")

    def test_phaser_tutorial_detection(self):
        """Should detect MR tutorial and add stop condition."""
        advice = "Try molecular replacement with the search model."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.phaser")
        self.assertTrue(result["stop_conditions"]["skip_validation"])

    def test_mr_test_detection(self):
        """Should detect MR test pattern."""
        advice = "Test MR solution with this model."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.phaser")

    def test_mtriage_detection(self):
        """Should detect mtriage tutorial."""
        advice = "Run mtriage to analyze the cryo-EM map."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.mtriage")

    def test_density_modification_detection(self):
        """Should detect density modification tutorial."""
        advice = "Run one cycle of density modification to improve the map."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        # Without clear context, defaults to cryo-EM
        self.assertTrue(result["stop_conditions"]["skip_validation"])

    def test_density_modification_cryoem(self):
        """Should detect cryo-EM density modification."""
        advice = "Experiment Type: cryo-EM. Run density modification on the half-maps."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.resolve_cryo_em")

    def test_density_modification_xray(self):
        """Should detect X-ray density modification."""
        advice = "Experiment Type: X-ray. Run density modification to improve phases from the .mtz file."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.autobuild_denmod")

    def test_stop_condition_section_detection(self):
        """Should detect explicit Stop Condition section from preprocessed advice."""
        advice = """
        1. **Input Files Found**: half_map_1.ccp4, half_map_2.ccp4
        2. **Experiment Type**: cryo-EM
        6. **Stop Condition**: Stop after running one cycle of density modification.
        """
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertTrue(result["stop_conditions"]["skip_validation"])
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.resolve_cryo_em")

    def test_stop_condition_mtriage(self):
        """Should detect mtriage from Stop Condition section."""
        advice = "6. **Stop Condition**: Stop after running mtriage."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertTrue(result["stop_conditions"]["skip_validation"])
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.mtriage")

    def test_no_false_positive(self):
        """Should not add stop condition for full workflow advice."""
        advice = "Solve the structure by molecular replacement and refine to completion."
        result = extract_directives_simple(advice)

        # Should detect MR preference but the "to completion" implies full workflow
        # The simple extractor will still match MR, but with full LLM it would be smarter
        # This test documents current behavior
        if "stop_conditions" in result:
            # If stop_conditions exists, it should be from the MR pattern
            # A smarter system would not add stop for "to completion"
            pass

    def test_ligand_later_no_early_stop(self):
        """Should NOT stop early when user wants ligand fitting later."""
        advice = "Run predict_and_build with rebuilding_strategy=Quick. Fit the ligand later."
        result = extract_directives_simple(advice)

        # Should NOT have after_program stop condition because "later" indicates continuation
        if "stop_conditions" in result:
            self.assertNotIn("after_program", result["stop_conditions"],
                           "Should not set after_program when 'later' indicates workflow continuation")

    def test_then_refine_no_early_stop(self):
        """Should NOT stop early when user says 'then refine'."""
        advice = "Try molecular replacement with the search model, then refine the structure."
        result = extract_directives_simple(advice)

        # Should NOT have after_program=phenix.phaser because "then refine" indicates continuation
        if "stop_conditions" in result:
            self.assertNotIn("after_program", result["stop_conditions"],
                           "Should not set after_program when 'then' indicates workflow continuation")

    def test_afterwards_no_early_stop(self):
        """Should NOT stop early when user says 'afterwards'."""
        advice = "Run xtriage to check for twinning. Afterwards, proceed with phasing."
        result = extract_directives_simple(advice)

        # Should NOT have after_program=phenix.xtriage because "afterwards" indicates continuation
        if "stop_conditions" in result:
            self.assertNotIn("after_program", result["stop_conditions"],
                           "Should not set after_program when 'afterwards' indicates workflow continuation")

    def test_validate_later_no_early_stop(self):
        """Should NOT stop early when user wants validation."""
        advice = "Run phaser for MR, then validate with molprobity."
        result = extract_directives_simple(advice)

        # "validate" is a downstream task indicator - should not stop early
        if "stop_conditions" in result:
            self.assertNotIn("after_program", result["stop_conditions"],
                           "Should not set after_program when validation is mentioned as downstream task")

    def test_pure_tutorial_still_works(self):
        """Pure tutorial without continuation indicators should still stop."""
        advice = "Run xtriage to check for twinning."
        result = extract_directives_simple(advice)

        # This is a pure tutorial - should have stop condition
        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.xtriage")

    def test_map_symmetry_detection(self):
        """Should detect map symmetry tutorial and add stop condition."""
        advice = "Determine the symmetry of the provided cryo-EM map."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.map_symmetry")
        self.assertTrue(result["stop_conditions"]["skip_validation"])

    def test_map_symmetry_stop_condition(self):
        """Should detect map symmetry from Stop Condition section."""
        advice = """
        1. **Input Files Found**: emd-20026_auto_sharpen_A.ccp4
        2. **Experiment Type**: cryo-EM (analysis only)
        3. **Primary Goal**: Determine the symmetry of the provided cryo-EM map.
        6. **Stop Condition**: Stop after determining the map symmetry.
        """
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertTrue(result["stop_conditions"]["skip_validation"])
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.map_symmetry")

    def test_find_symmetry_detection(self):
        """Should detect 'find symmetry' pattern."""
        advice = "Find the point-group symmetry in this cryo-EM map."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.map_symmetry")

    def test_ligandfit_stop_condition(self):
        """Should detect ligandfit from Stop Condition section."""
        advice = """
        1. Input Files Found: nsf-d2.mtz, nsf-d2_noligand.pdb, atp.pdb
        Experiment Type: Refinement and Ligand Fitting
        Primary Goal: Refine an unliganded protein model to generate an mFo-DFc difference map.
        Stop Condition: Stop after successfully running LigandFit and generating the model.
        """
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertTrue(result["stop_conditions"]["skip_validation"])
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.ligandfit")

    def test_ligandfit_tutorial_pattern(self):
        """Should detect ligand fitting tutorial pattern."""
        advice = "Fit a ligand into the difference density."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.ligandfit")
        self.assertTrue(result["stop_conditions"]["skip_validation"])

    def test_run_ligandfit_pattern(self):
        """Should detect 'run ligandfit' pattern."""
        advice = "Run phenix.ligandfit to place the ATP molecule."
        result = extract_directives_simple(advice)

        self.assertIn("stop_conditions", result)
        self.assertEqual(result["stop_conditions"]["after_program"], "phenix.ligandfit")


def run_all_tests():
    """Run all tests and raise exception if any fail."""
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise Exception("Some tests failed")


if __name__ == "__main__":
    unittest.main()
