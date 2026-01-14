#!/usr/bin/env python
"""
Test script for the dry run system.

This tests the DryRunManager independently of the full agent,
verifying that it correctly loads scenarios and returns dummy data.
"""

from __future__ import absolute_import, division, print_function
import os
import tempfile
import shutil

from libtbx.langchain.agent.dry_run_manager import DryRunManager, list_scenarios, get_scenario_info


def test_list_scenarios():
    """Test listing available scenarios."""
    print("Test: list_scenarios")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    scenarios = list_scenarios(scenarios_dir)
    print(f"  Found {len(scenarios)} scenario(s): {scenarios}")

    assert len(scenarios) >= 1, "Should have at least one scenario"
    assert "xray_basic" in scenarios, "Should have xray_basic scenario"

    print("  PASSED")


def test_load_scenario():
    """Test loading a scenario."""
    print("Test: load_scenario")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    manager = DryRunManager("xray_basic", scenarios_dir)

    assert manager.scenario is not None, "Should load scenario"
    assert manager.scenario.get('name') == 'xray_basic', "Should have correct name"
    assert manager.scenario.get('experiment_type') == 'xray', "Should be xray type"

    steps = manager.scenario.get('steps', [])
    assert len(steps) >= 5, "Should have at least 5 steps"

    print(f"  Scenario has {len(steps)} steps")
    print("  PASSED")


def test_get_initial_files():
    """Test getting initial input files."""
    print("Test: get_initial_files")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    manager = DryRunManager("xray_basic", scenarios_dir)
    files = manager.get_initial_files()

    assert len(files) >= 2, "Should have at least 2 input files"

    basenames = [os.path.basename(f) for f in files]
    assert "data.mtz" in basenames, "Should have data.mtz"
    assert "sequence.fa" in basenames, "Should have sequence.fa"

    print(f"  Initial files: {basenames}")
    print("  PASSED")


def test_get_dummy_result_xtriage():
    """Test getting dummy result for xtriage."""
    print("Test: get_dummy_result_xtriage")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    manager = DryRunManager("xray_basic", scenarios_dir)

    # Create temp working directory
    work_dir = tempfile.mkdtemp()

    try:
        log_text, error_text, output_files = manager.get_dummy_result(
            "phenix.xtriage data.mtz",
            work_dir
        )

        assert len(log_text) > 100, "Should have substantial log text"
        assert "Resolution" in log_text, "Log should mention resolution"
        assert error_text == "", "Should have no error"

        print(f"  Log length: {len(log_text)} chars")
        print("  PASSED")

    finally:
        shutil.rmtree(work_dir)


def test_get_dummy_result_refine_multiple():
    """Test getting dummy results for multiple refine calls."""
    print("Test: get_dummy_result_refine_multiple")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    manager = DryRunManager("xray_basic", scenarios_dir)
    work_dir = tempfile.mkdtemp()

    try:
        # Skip xtriage and predict_and_build
        manager.get_dummy_result("phenix.xtriage data.mtz", work_dir)
        manager.get_dummy_result("phenix.predict_and_build seq.fa", work_dir)

        # Now test multiple refine calls
        r_free_values = []

        for i in range(3):
            log_text, error_text, output_files = manager.get_dummy_result(
                "phenix.refine model.pdb data.mtz",
                work_dir
            )

            # Extract R-free from log
            for line in log_text.split('\n'):
                if 'Final R-free' in line:
                    # Parse "Final R-free: 0.325"
                    parts = line.split(':')
                    if len(parts) >= 2:
                        try:
                            r_free = float(parts[-1].strip())
                            r_free_values.append(r_free)
                        except ValueError:
                            pass

            assert len(output_files) == 1, f"Refine {i+1} should produce 1 output"
            assert output_files[0].endswith('.pdb'), "Output should be PDB"

        print(f"  R-free progression: {r_free_values}")

        # Verify R-free is decreasing
        assert len(r_free_values) == 3, "Should have 3 R-free values"
        assert r_free_values[0] > r_free_values[1] > r_free_values[2], \
            "R-free should decrease with each cycle"

        print("  PASSED")

    finally:
        shutil.rmtree(work_dir)


def test_scenario_info():
    """Test getting scenario information."""
    print("Test: scenario_info")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    info = get_scenario_info("xray_basic", scenarios_dir)

    assert info['name'] == 'xray_basic', "Should have correct name"
    assert info['experiment_type'] == 'xray', "Should be xray"
    assert info['num_steps'] >= 5, "Should have at least 5 steps"

    programs = info['programs']
    assert 'phenix.xtriage' in programs, "Should include xtriage"
    assert 'phenix.refine' in programs, "Should include refine"
    assert 'phenix.molprobity' in programs, "Should include molprobity"

    print(f"  Programs: {programs}")
    print("  PASSED")


def test_full_scenario_walkthrough():
    """Test walking through entire scenario."""
    print("Test: full_scenario_walkthrough")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    manager = DryRunManager("xray_basic", scenarios_dir)
    work_dir = tempfile.mkdtemp()

    try:
        steps = manager.scenario.get('steps', [])
        programs_run = []

        for step in steps:
            program = step['program']
            log_text, error_text, output_files = manager.get_dummy_result(
                f"{program} dummy_args",
                work_dir
            )

            programs_run.append(program)
            print(f"    {program}: log={len(log_text)} chars, outputs={len(output_files)}")

        print(f"  Completed {len(programs_run)} steps")

        # Verify we ran all expected programs
        assert 'phenix.xtriage' in programs_run
        assert 'phenix.molprobity' in programs_run
        assert programs_run.count('phenix.refine') == 3

        print("  PASSED")

    finally:
        shutil.rmtree(work_dir)


def test_cryoem_scenario():
    """Test the cryo-EM scenario."""
    print("Test: cryoem_scenario")

    scenarios_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scenarios"
    )

    # Check scenario exists
    scenarios = list_scenarios(scenarios_dir)
    assert "cryoem_basic" in scenarios, "Should have cryoem_basic scenario"

    manager = DryRunManager("cryoem_basic", scenarios_dir)
    work_dir = tempfile.mkdtemp()

    try:
        # Verify initial files
        initial_files = manager.get_initial_files()
        basenames = [os.path.basename(f) for f in initial_files]
        assert "map.mrc" in basenames, "Should have map.mrc"
        assert "sequence.fa" in basenames, "Should have sequence.fa"

        # Walk through scenario
        steps = manager.scenario.get('steps', [])
        programs_run = []
        cc_values = []

        for step in steps:
            program = step['program']
            log_text, error_text, output_files = manager.get_dummy_result(
                f"{program} dummy_args",
                work_dir
            )
            programs_run.append(program)

            # Extract CC from real_space_refine logs
            if 'real_space_refine' in program:
                for line in log_text.split('\n'):
                    if 'Final map-model CC' in line:
                        parts = line.split(':')
                        if len(parts) >= 2:
                            try:
                                cc = float(parts[-1].strip())
                                cc_values.append(cc)
                            except ValueError:
                                pass

        print(f"  Programs: {programs_run}")
        print(f"  CC progression: {cc_values}")

        # Verify workflow
        assert 'phenix.mtriage' in programs_run, "Should run mtriage"
        assert 'phenix.predict_and_build' in programs_run, "Should run predict_and_build"
        assert programs_run.count('phenix.real_space_refine') == 3, "Should run 3 RSR cycles"
        assert 'phenix.molprobity' in programs_run, "Should run molprobity"

        # Verify CC is increasing
        assert len(cc_values) == 3, "Should have 3 CC values"
        assert cc_values[0] < cc_values[1] < cc_values[2], "CC should increase"
        assert cc_values[-1] >= 0.75, "Final CC should be >= 0.75"

        print("  PASSED")

    finally:
        shutil.rmtree(work_dir)


def run_all_tests():
    """Run all dry run tests."""
    print("=" * 60)
    print("DRY RUN MANAGER TESTS")
    print("=" * 60)
    print()

    test_list_scenarios()
    test_load_scenario()
    test_get_initial_files()
    test_get_dummy_result_xtriage()
    test_get_dummy_result_refine_multiple()
    test_scenario_info()
    test_full_scenario_walkthrough()
    test_cryoem_scenario()

    print()
    print("=" * 60)
    print("ALL DRY RUN TESTS PASSED!")
    print("=" * 60)


if __name__ == '__main__':
    run_all_tests()
