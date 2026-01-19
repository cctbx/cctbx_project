"""
Tests for YAML-driven configuration.

Tests:
1. YAML files load correctly
2. Program registry works with YAML
3. Metrics can be extracted from logs
4. Workflow phases are defined correctly
"""

from __future__ import absolute_import, division, print_function


def test_yaml_loading():
    """Test that all YAML files load correctly."""
    print("Test: yaml_loading")

    from libtbx.langchain.knowledge.yaml_loader import (
        load_programs, load_workflows, load_metrics, validate_yaml_files
    )

    # Validate files
    is_valid, errors = validate_yaml_files()
    assert is_valid, "YAML validation failed: %s" % errors

    # Check programs loaded
    programs = load_programs()
    assert len(programs) > 0, "No programs loaded"
    assert "phenix.refine" in programs, "phenix.refine not found"
    assert "phenix.xtriage" in programs, "phenix.xtriage not found"

    # Check workflows loaded
    workflows = load_workflows()
    assert "xray" in workflows, "xray workflow not found"
    assert "cryoem" in workflows, "cryoem workflow not found"

    # Check metrics loaded
    metrics = load_metrics()
    assert "r_free" in metrics, "r_free metric not found"
    assert "map_cc" in metrics, "map_cc metric not found"

    print("  PASSED")


def test_program_registry():
    """Test ProgramRegistry with YAML backend."""
    print("Test: program_registry")

    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry(use_yaml=True)

    # Check YAML is being used
    assert registry.use_yaml, "Registry not using YAML"

    # Check program info
    refine = registry.get_program("phenix.refine")
    assert refine is not None, "phenix.refine not found"
    assert "description" in refine, "No description for phenix.refine"

    # Check required inputs
    required = registry.get_required_inputs("phenix.refine")
    assert "model" in required, "model not in required inputs"
    assert "mtz" in required, "mtz not in required inputs"

    # Check command building
    files = {"model": "test.pdb", "mtz": "test.mtz"}
    cmd = registry.build_command("phenix.refine", files)
    assert "phenix.refine" in cmd, "Command doesn't start with program name"
    assert "test.pdb" in cmd, "Model file not in command"
    assert "test.mtz" in cmd, "MTZ file not in command"

    # Check strategy flags
    files = {"model": "test.pdb", "mtz": "test.mtz"}
    strategy = {"generate_rfree_flags": True}
    cmd = registry.build_command("phenix.refine", files, strategy)
    assert "r_free_flags.generate" in cmd, "Strategy flag not in command"

    print("  PASSED")


def test_experiment_type_filtering():
    """Test filtering programs by experiment type."""
    print("Test: experiment_type_filtering")

    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry(use_yaml=True)

    # X-ray programs
    xray_progs = registry.get_programs_for_experiment("xray")
    assert "phenix.refine" in xray_progs, "phenix.refine not in xray programs"
    assert "phenix.xtriage" in xray_progs, "phenix.xtriage not in xray programs"
    assert "phenix.real_space_refine" not in xray_progs, "real_space_refine in xray programs"

    # Cryo-EM programs
    cryoem_progs = registry.get_programs_for_experiment("cryoem")
    assert "phenix.real_space_refine" in cryoem_progs, "real_space_refine not in cryoem"
    assert "phenix.mtriage" in cryoem_progs, "mtriage not in cryoem"
    assert "phenix.refine" not in cryoem_progs, "phenix.refine in cryoem programs"

    # Both types
    assert "phenix.predict_and_build" in xray_progs, "predict_and_build not in xray"
    assert "phenix.predict_and_build" in cryoem_progs, "predict_and_build not in cryoem"

    print("  PASSED")


def test_metric_thresholds():
    """Test resolution-dependent metric thresholds."""
    print("Test: metric_thresholds")

    from libtbx.langchain.knowledge.yaml_loader import (
        get_metric_threshold, get_target_r_free, is_metric_acceptable
    )

    # Test R-free thresholds
    rfree_good_15 = get_metric_threshold("r_free", "good", resolution=1.5)
    rfree_good_30 = get_metric_threshold("r_free", "good", resolution=3.0)
    assert rfree_good_15 < rfree_good_30, "R-free threshold not resolution-dependent"

    # Test target R-free
    target_20 = get_target_r_free(2.0)
    target_30 = get_target_r_free(3.0)
    assert target_20 < target_30, "Target R-free not resolution-dependent"
    assert 0.20 < target_20 < 0.30, "Target R-free at 2.0A out of range"

    # Test is_acceptable
    assert is_metric_acceptable("r_free", 0.25, resolution=2.0), "0.25 should be acceptable at 2.0A"
    assert not is_metric_acceptable("r_free", 0.40, resolution=2.0), "0.40 should not be acceptable"

    print("  PASSED")


def test_log_pattern_extraction():
    """Test extracting metrics from log text."""
    print("Test: log_pattern_extraction")

    from libtbx.langchain.knowledge.yaml_loader import extract_metric_from_log

    # Test R-free extraction
    refine_log = """
    Final R-work = 0.1856, R-free = 0.2234
    Bonds RMSD: 0.012
    """

    r_free = extract_metric_from_log(refine_log, "r_free", "phenix.refine")
    assert r_free is not None, "Failed to extract R-free"
    assert abs(r_free - 0.2234) < 0.001, "R-free value incorrect: %s" % r_free

    # Test TFZ extraction (should get last value)
    phaser_log = """
    SOLU SET  RFZ=4.6 TFZ=8.3 PAK=2 LLG=94 TFZ==9.8 LLG=544 TFZ==24.0
    """

    tfz = extract_metric_from_log(phaser_log, "tfz", "phenix.phaser")
    assert tfz is not None, "Failed to extract TFZ"
    assert abs(tfz - 24.0) < 0.1, "TFZ value incorrect (should be last): %s" % tfz

    print("  PASSED")


def test_measurability_parsing():
    """Test anomalous measurability parsing from xtriage logs."""
    print("Test: measurability_parsing")

    from phenix.phenix_ai.log_parsers import parse_xtriage_log

    # Test with actual anomalous data (should have measurability)
    xtriage_log_anomalous = """
    Anomalous signal analysis:

    Measurability of anomalous signal:
    | Resolution | Count  | Measurability |
    |  50.0-5.0  | [23/30] | 0.183 |
    |  5.0-4.0   | [45/52] | 0.142 |
    |  4.0-3.5   | [38/48] | 0.098 |
    |  3.5-3.0   | [22/40] | 0.065 |

    The anomalous signal appears to be significant.
    """

    result_anom = parse_xtriage_log(xtriage_log_anomalous)

    assert "anomalous_measurability" in result_anom, \
        "Should have anomalous_measurability for anomalous data"
    assert result_anom["anomalous_measurability"] > 0.05, \
        "Measurability should be > 0.05 for anomalous data"
    assert result_anom["anomalous_measurability"] < 0.3, \
        "Measurability should be reasonable (< 0.3)"

    # Test with non-anomalous data (should NOT have measurability)
    xtriage_log_normal = """
    Summary of data:
      Resolution: 2.5 Angstroms
      Completeness: 98.5%

    Data completeness table:
    | Resolution | Count  | Completeness |
    |  50.0-5.0  | [100/100] | 1.000 |
    |  5.0-3.0   | [500/510] | 0.980 |
    |  3.0-2.5   | [300/320] | 0.938 |

    Twinning analysis:
      No twinning detected.
    """

    result_normal = parse_xtriage_log(xtriage_log_normal)

    # Should NOT have measurability (or it should be None/absent)
    meas = result_normal.get("anomalous_measurability")
    assert meas is None or meas == 0, \
        "Non-anomalous data should not have measurability, got: %s" % meas

    print("  PASSED")


def test_measurability_sanity_check():
    """Test that measurability values are sanity-checked (0-0.5 range)."""
    print("Test: measurability_sanity_check")

    from phenix.phenix_ai.log_parsers import parse_xtriage_log

    # Log with a table that looks like measurability but has 1.0 values
    # (this is actually a completeness table, not measurability)
    fake_meas_log = """
    Data statistics:

    Completeness analysis:
    | Resolution | Count  | Completeness |
    |  50.0-5.0  | [100/100] | 1.000 |
    |  5.0-3.0   | [500/500] | 1.000 |
    """

    result = parse_xtriage_log(fake_meas_log)

    # Should NOT extract 1.0 as measurability (sanity check should reject it)
    meas = result.get("anomalous_measurability")
    assert meas is None or meas < 0.5, \
        "Should not extract 1.0 as measurability (sanity check), got: %s" % meas

    print("  PASSED")


def test_workflow_phases():
    """Test workflow phase definitions."""
    print("Test: workflow_phases")

    from libtbx.langchain.knowledge.yaml_loader import get_workflow_phases, get_phase_programs

    # X-ray phases
    xray_phases = get_workflow_phases("xray")
    assert "analyze" in xray_phases, "analyze phase not found"
    assert "refine" in xray_phases, "refine phase not found"
    assert "validate" in xray_phases, "validate phase not found"

    # Phase programs
    analyze_progs = get_phase_programs("xray", "analyze")
    assert len(analyze_progs) > 0, "No programs for analyze phase"

    # Cryo-EM phases
    cryoem_phases = get_workflow_phases("cryoem")
    assert "analyze" in cryoem_phases, "analyze phase not found in cryoem"

    print("  PASSED")


def test_command_with_strategy():
    """Test command building with various strategy options."""
    print("Test: command_with_strategy")

    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry(use_yaml=True)

    # Test predict_and_build with stop_after_predict
    files = {"sequence": "seq.fa", "mtz": "data.mtz"}
    strategy = {"stop_after_predict": True, "resolution": 2.5}
    cmd = registry.build_command("phenix.predict_and_build", files, strategy)

    assert "stop_after_predict" in cmd, "stop_after_predict not in command"
    assert "resolution" in cmd or "2.5" in cmd, "resolution not in command"

    # Test refine with twin law
    files = {"model": "model.pdb", "mtz": "data.mtz"}
    strategy = {"twin_law": "-h,-k,l"}
    cmd = registry.build_command("phenix.refine", files, strategy)

    assert "twin_law" in cmd or "twinning" in cmd, "twin_law not in command"

    print("  PASSED")


def test_workflow_engine():
    """Test WorkflowEngine for phase detection."""
    print("Test: workflow_engine")

    from libtbx.langchain.agent.workflow_engine import WorkflowEngine

    engine = WorkflowEngine()

    # Test X-ray initial state
    files = {"mtz": ["data.mtz"], "sequence": ["seq.fa"]}
    history = {}
    context = engine.build_context(files, history)
    phase = engine.detect_phase("xray", context)

    assert phase["phase"] == "analyze", "X-ray initial should be analyze phase"

    programs = engine.get_valid_programs("xray", phase, context)
    assert "phenix.xtriage" in programs, "xtriage should be valid in analyze phase"

    # Test X-ray after xtriage
    history = {"xtriage_done": True}
    context = engine.build_context(files, history)
    phase = engine.detect_phase("xray", context)

    assert phase["phase"] == "obtain_model", "After xtriage should be obtain_model"

    programs = engine.get_valid_programs("xray", phase, context)
    assert "phenix.predict_and_build" in programs, "predict_and_build should be valid"

    # Test cryo-EM initial state
    files = {"full_map": ["map.mrc"], "sequence": ["seq.fa"]}
    history = {}
    context = engine.build_context(files, history)
    phase = engine.detect_phase("cryoem", context)

    assert phase["phase"] == "analyze", "Cryo-EM initial should be analyze"

    programs = engine.get_valid_programs("cryoem", phase, context)
    assert "phenix.mtriage" in programs, "mtriage should be valid"

    print("  PASSED")


def test_workflow_engine_refined_state():
    """Test WorkflowEngine for refined/validation states."""
    print("Test: workflow_engine_refined_state")

    from libtbx.langchain.agent.workflow_engine import WorkflowEngine

    engine = WorkflowEngine()

    # Test X-ray refined state with good R-free
    files = {
        "mtz": ["data.mtz"],
        "sequence": ["seq.fa"],
        "pdb": ["model.pdb"],
        "refined": ["model_refine_001.pdb"],
        "phaser_output": ["PHASER.1.pdb"],
    }
    history = {
        "xtriage_done": True,
        "phaser_done": True,
        "refine_done": True,
        "refine_count": 3,
    }
    analysis = {"r_free": 0.24, "resolution": 2.0}

    context = engine.build_context(files, history, analysis)
    phase = engine.detect_phase("xray", context)

    # Should need validation since R-free is good
    assert phase["phase"] in ["validate", "refine"], "Should be validate or refine phase"

    # Test cryo-EM with half-maps only
    files = {
        "half_map": ["half1.mrc", "half2.mrc"],
        "sequence": ["seq.fa"],
        "pdb": ["model.pdb"],
    }
    history = {
        "mtriage_done": True,
        "predict_full_done": True,  # Model placed
    }

    context = engine.build_context(files, history)
    phase = engine.detect_phase("cryoem", context)

    # Should need to optimize map since only half-maps
    assert phase["phase"] == "optimize_map", "Should need map optimization: %s" % phase

    print("  PASSED")


def test_yaml_workflow_state_detection():
    """Test that YAML engine correctly detects workflow states."""
    print("Test: yaml_workflow_state_detection")

    from libtbx.langchain.agent.workflow_state import detect_workflow_state

    # Test case 1: X-ray initial
    available_files = ["data.mtz", "seq.fa"]
    history = []

    result = detect_workflow_state(history, available_files, use_yaml_engine=True)

    assert result["experiment_type"] == "xray", \
        "Expected xray, got %s" % result["experiment_type"]
    assert "phenix.xtriage" in result["valid_programs"], \
        "xtriage not in valid_programs: %s" % result["valid_programs"]

    # Test case 2: Cryo-EM initial
    available_files = ["map.mrc", "seq.fa"]
    history = []

    result = detect_workflow_state(history, available_files, use_yaml_engine=True)

    assert result["experiment_type"] == "cryoem", \
        "Expected cryoem, got %s" % result["experiment_type"]
    assert "phenix.mtriage" in result["valid_programs"], \
        "mtriage not in valid_programs: %s" % result["valid_programs"]

    # Test case 3: X-ray after xtriage - should offer model building options
    available_files = ["data.mtz", "seq.fa"]
    history = [{"program": "phenix.xtriage", "result": "SUCCESS"}]

    result = detect_workflow_state(history, available_files, use_yaml_engine=True)

    assert "phenix.predict_and_build" in result["valid_programs"], \
        "predict_and_build should be valid after xtriage"

    print("  PASSED")


def test_metric_evaluator():
    """Test MetricEvaluator for quality assessment."""
    print("Test: metric_evaluator")

    from libtbx.langchain.agent.metric_evaluator import MetricEvaluator

    evaluator = MetricEvaluator()

    # Test threshold access
    target = evaluator.get_target("r_free", resolution=2.0)
    assert target is not None, "R-free target should not be None"
    assert 0.20 < target < 0.35, "R-free target out of range: %s" % target

    # Test quality levels
    level = evaluator.get_quality_level("r_free", 0.22, resolution=2.0)
    assert level in ["good", "acceptable"], "0.22 should be good/acceptable at 2.0A"

    level = evaluator.get_quality_level("r_free", 0.40, resolution=2.0)
    assert level == "poor", "0.40 should be poor at 2.0A"

    # Test quality assessment
    metrics = {"r_free": 0.24, "clashscore": 3.5}
    quality = evaluator.assess_quality(metrics, resolution=2.0)
    assert quality["overall"] in ["good", "acceptable"], "Quality should be acceptable"
    assert "r_free" in quality["details"], "R-free should be in details"

    # Test improvement detection
    assert evaluator.is_significant_improvement("r_free", 0.30, 0.28), \
        "0.30->0.28 should be significant"
    assert not evaluator.is_significant_improvement("r_free", 0.30, 0.299), \
        "0.30->0.299 should not be significant"

    print("  PASSED")


def test_metric_evaluator_trend():
    """Test MetricEvaluator for trend analysis."""
    print("Test: metric_evaluator_trend")

    from libtbx.langchain.agent.metric_evaluator import MetricEvaluator

    evaluator = MetricEvaluator()

    # Test improving trend
    history = [
        {"program": "phenix.refine", "r_free": 0.35},
        {"program": "phenix.refine", "r_free": 0.30},
        {"program": "phenix.refine", "r_free": 0.26},
    ]
    trend = evaluator.analyze_trend(history, "xray", resolution=2.0)

    assert trend["recommendation"] == "continue", \
        "Improving trend should recommend continue: %s" % trend
    assert trend["improvement_rate"] > 0, \
        "Improvement rate should be positive: %s" % trend["improvement_rate"]

    # Test success case
    history = [
        {"program": "phenix.refine", "r_free": 0.30},
        {"program": "phenix.refine", "r_free": 0.25},
        {"program": "phenix.refine", "r_free": 0.22},
        {"program": "phenix.molprobity"},  # Validation done
    ]
    trend = evaluator.analyze_trend(history, "xray", resolution=2.5)

    # At 2.5A, target is ~0.30, so 0.22 is well below
    assert trend["should_stop"] or "TARGET" in trend.get("trend_summary", ""), \
        "Should recommend stop or show target reached"

    print("  PASSED")


def test_metrics_analyzer_yaml_mode():
    """Test that metrics_analyzer can use YAML evaluator."""
    print("Test: metrics_analyzer_yaml_mode")

    from libtbx.langchain.agent.metrics_analyzer import analyze_metrics_trend

    # Test with YAML evaluator
    history = [
        {"program": "phenix.refine", "r_free": 0.35},
        {"program": "phenix.refine", "r_free": 0.30},
    ]

    # Hardcoded mode
    hardcoded = analyze_metrics_trend(history, resolution=2.0, use_yaml_evaluator=False)

    # YAML mode
    yaml_result = analyze_metrics_trend(history, resolution=2.0, use_yaml_evaluator=True)

    # Both should have similar structure
    assert "should_stop" in hardcoded, "Missing should_stop in hardcoded"
    assert "should_stop" in yaml_result, "Missing should_stop in YAML"
    assert "recommendation" in hardcoded, "Missing recommendation in hardcoded"
    assert "recommendation" in yaml_result, "Missing recommendation in YAML"

    # Both should agree on basic assessment
    assert hardcoded["should_stop"] == yaml_result["should_stop"], \
        "Should_stop differs: %s vs %s" % (hardcoded["should_stop"], yaml_result["should_stop"])

    print("  PASSED")


def test_rules_selector():
    """Test RulesSelector for program selection."""
    print("Test: rules_selector")

    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()

    # Test X-ray initial state - should select xtriage
    workflow_state = {
        "state": "analyze",
        "experiment_type": "xray",
        "valid_programs": ["phenix.xtriage"],
        "resolution": 2.0,
    }
    files = {"mtz": ["data.mtz"], "sequence": ["seq.fa"]}

    intent = selector.select_next_action(workflow_state, files)

    assert intent["program"] == "phenix.xtriage", \
        "Should select xtriage: %s" % intent["program"]
    assert not intent["stop"], "Should not stop"
    assert intent["selection_method"] == "rules", "Should be rules-based"

    # Test refinement state with good metrics - should prefer validation
    workflow_state = {
        "state": "refine",
        "experiment_type": "xray",
        "valid_programs": ["phenix.refine", "phenix.molprobity", "STOP"],
        "resolution": 2.0,
    }
    files = {"mtz": ["data.mtz"], "pdb": ["model.pdb"]}
    metrics_trend = {
        "r_free_trend": [0.30, 0.27, 0.25],
        "should_stop": False,
    }

    intent = selector.select_next_action(workflow_state, files, metrics_trend)

    # Should prefer molprobity since R-free is good
    assert intent["program"] == "phenix.molprobity", \
        "Should select molprobity for validation: %s" % intent["program"]

    print("  PASSED")


def test_rules_selector_stop():
    """Test RulesSelector stop conditions."""
    print("Test: rules_selector_stop")

    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()

    # Test explicit stop
    workflow_state = {
        "state": "complete",
        "valid_programs": ["STOP"],
    }
    metrics_trend = {"should_stop": True, "reason": "Target reached"}

    intent = selector.select_next_action(workflow_state, {}, metrics_trend)

    assert intent["stop"], "Should stop"
    assert intent["program"] is None, "Should have no program"
    assert "Target reached" in intent.get("stop_reason", ""), \
        "Should include stop reason"

    print("  PASSED")


# =============================================================================
# STEP 1 TESTS: INVARIANTS
# =============================================================================

def test_invariants_yaml_structure():
    """Test that invariants are properly defined in programs.yaml."""
    print("Test: invariants_yaml_structure")

    from libtbx.langchain.knowledge.yaml_loader import load_programs

    programs = load_programs()

    # Check predict_and_build has invariants
    pab = programs.get("phenix.predict_and_build", {})
    assert "invariants" in pab, "predict_and_build should have invariants"

    invariants = pab["invariants"]
    assert isinstance(invariants, list), "invariants should be a list"
    assert len(invariants) >= 2, "predict_and_build should have at least 2 invariants"

    # Check invariant structure
    for inv in invariants:
        assert "name" in inv, "invariant should have 'name'"
        assert "check" in inv, "invariant should have 'check'"
        assert "fix" in inv, "invariant should have 'fix'"

    # Check real_space_refine has resolution invariant
    rsr = programs.get("phenix.real_space_refine", {})
    assert "invariants" in rsr, "real_space_refine should have invariants"

    inv_names = [i.get("name") for i in rsr["invariants"]]
    assert "requires_resolution" in inv_names, "RSR should have requires_resolution invariant"

    print("  PASSED")


def test_invariants_validate_and_fix():
    """Test that invariants are checked and fixed by TemplateBuilder."""
    print("Test: invariants_validate_and_fix")

    from libtbx.langchain.agent.template_builder import TemplateBuilder

    builder = TemplateBuilder()

    # Test predict_and_build: no data file -> should set stop_after_predict
    files = {"sequence": "test.seq"}
    strategy = {}
    context = {}

    files2, strategy2, warnings = builder.validate_and_fix(
        "phenix.predict_and_build", files, strategy, context=context
    )

    assert strategy2.get("stop_after_predict") == True, \
        "Should set stop_after_predict=True when no data file"
    assert len(warnings) > 0, "Should have warnings about fix applied"

    # Test real_space_refine: no resolution -> should auto-fill from context
    files = {"model": "test.pdb", "map": "test.mrc"}
    strategy = {}
    context = {"resolution": 3.5}

    files2, strategy2, warnings = builder.validate_and_fix(
        "phenix.real_space_refine", files, strategy, context=context
    )

    assert "resolution" in strategy2, "Should have resolution after fix"
    assert abs(strategy2["resolution"] - 3.5) < 0.01, "Resolution should be 3.5"

    # Test when resolution already set - should not override
    files = {"model": "test.pdb", "map": "test.mrc"}
    strategy = {"resolution": 2.8, "output_prefix": "rsr_001"}  # Both required
    context = {"resolution": 3.5}

    files3, strategy3, warnings3 = builder.validate_and_fix(
        "phenix.real_space_refine", files, strategy, context=context
    )

    assert strategy3["resolution"] == 2.8, "Should not override existing resolution"
    assert strategy3["output_prefix"] == "rsr_001", "Should not override existing output_prefix"
    assert len(warnings3) == 0, "Should have no warnings when all invariants satisfied"

    print("  PASSED")


def test_invariants_check_types():
    """Test different check types in invariants."""
    print("Test: invariants_check_types")

    from libtbx.langchain.agent.template_builder import TemplateBuilder

    builder = TemplateBuilder()

    # Test has_file check
    files_with_mtz = {"mtz": "test.mtz", "sequence": "test.seq"}
    files_no_mtz = {"sequence": "test.seq"}
    strategy = {}

    # With MTZ - should pass data_or_predict_only check
    f1, s1, w1 = builder.validate_and_fix(
        "phenix.predict_and_build", files_with_mtz, strategy, context={}
    )
    assert s1.get("stop_after_predict") != True, \
        "With MTZ file, should NOT set stop_after_predict"

    # Without MTZ - should fail check and apply fix
    f2, s2, w2 = builder.validate_and_fix(
        "phenix.predict_and_build", files_no_mtz, strategy, context={}
    )
    assert s2.get("stop_after_predict") == True, \
        "Without MTZ file, SHOULD set stop_after_predict"

    # Test has_strategy check
    strategy_with_res = {"resolution": 2.5}
    strategy_no_res = {}

    f3, s3, w3 = builder.validate_and_fix(
        "phenix.real_space_refine",
        {"model": "t.pdb", "map": "t.mrc"},
        strategy_with_res,
        context={"resolution": 3.0}
    )
    assert s3["resolution"] == 2.5, "Should keep existing resolution"

    print("  PASSED")


# =============================================================================
# STEP 2 TESTS: INPUT PRIORITIES
# =============================================================================

def test_input_priorities_yaml_structure():
    """Test that input_priorities are properly defined in programs.yaml."""
    print("Test: input_priorities_yaml_structure")

    from libtbx.langchain.knowledge.yaml_loader import load_programs

    programs = load_programs()

    # Check phenix.refine has input_priorities
    refine = programs.get("phenix.refine", {})
    assert "input_priorities" in refine, "phenix.refine should have input_priorities"

    priorities = refine["input_priorities"]
    assert "model" in priorities, "Should have model priorities"
    assert "mtz" in priorities, "Should have mtz priorities"

    model_prio = priorities["model"]
    assert "categories" in model_prio, "model should have categories"
    assert "exclude_categories" in model_prio, "model should have exclude_categories"
    assert "predicted" in model_prio["exclude_categories"], \
        "Should exclude predicted for refine"

    # Check real_space_refine has different priorities
    rsr = programs.get("phenix.real_space_refine", {})
    assert "input_priorities" in rsr, "real_space_refine should have input_priorities"

    rsr_model = rsr["input_priorities"]["model"]
    assert "rsr_output" in rsr_model["categories"], "RSR should prefer rsr_output"
    assert "processed_predicted" in rsr_model["exclude_categories"], \
        "RSR should exclude processed_predicted"

    print("  PASSED")


def test_input_priorities_registry():
    """Test ProgramRegistry.get_input_priorities() method."""
    print("Test: input_priorities_registry")

    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry(use_yaml=True)

    # Test phenix.refine priorities
    refine_model = registry.get_input_priorities("phenix.refine", "model")
    assert "categories" in refine_model, "Should have categories"
    assert len(refine_model["categories"]) > 0, "Should have at least one category"
    assert "with_ligand" in refine_model["categories"], "Should include with_ligand"

    refine_mtz = registry.get_input_priorities("phenix.refine", "mtz")
    assert "refined_mtz" in refine_mtz["categories"], "Should prefer refined_mtz"

    # Test real_space_refine priorities
    rsr_model = registry.get_input_priorities("phenix.real_space_refine", "model")
    assert "rsr_output" in rsr_model["categories"], "RSR should prefer rsr_output"

    # Test program without priorities
    xtriage = registry.get_input_priorities("phenix.xtriage", "model")
    assert xtriage["categories"] == [], "xtriage has no model input priorities"

    print("  PASSED")


def test_input_priorities_selection():
    """Test that input priorities affect file selection in RulesSelector."""
    print("Test: input_priorities_selection")

    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()

    # Create categorized files
    categorized_files = {
        "pdb": ["model.pdb"],
        "refined": ["refined_001.pdb"],
        "rsr_output": ["model_real_space_refined.pdb"],
        "predicted": ["predicted_model.pdb"],
        "processed_predicted": ["processed_model.pdb"],
        "mtz": ["data.mtz"],
        "refined_mtz": ["refine_001.mtz"],
        "sequence": ["seq.fa"],
        "full_map": ["map.mrc"],
    }

    workflow_state = {
        "state": "refine",
        "experiment_type": "xray",
        "valid_programs": ["phenix.refine"],
        "categorized_files": categorized_files,
    }

    # Test X-ray refinement - should prefer refined model
    intent = selector.select_next_action(workflow_state, {}, {})

    if intent.get("files", {}).get("model"):
        model = intent["files"]["model"]
        assert "refined" in model.lower() or model == "refined_001.pdb", \
            "Should select refined model for X-ray refine, got: %s" % model

    # Test cryo-EM refinement - should prefer rsr_output
    workflow_state_cryoem = {
        "state": "refine",
        "experiment_type": "cryoem",
        "valid_programs": ["phenix.real_space_refine"],
        "categorized_files": categorized_files,
        "resolution": 3.0,
    }

    intent2 = selector.select_next_action(workflow_state_cryoem, {}, {})

    if intent2.get("files", {}).get("model"):
        model2 = intent2["files"]["model"]
        assert "real_space_refined" in model2.lower() or "rsr" in model2.lower(), \
            "Should select RSR output for cryo-EM refine, got: %s" % model2

    print("  PASSED")


# =============================================================================
# STEP 3 TESTS: USER ADVICE KEYWORDS AND PRIORITY_WHEN
# =============================================================================

def test_user_advice_keywords_yaml_structure():
    """Test that user_advice_keywords are properly defined in programs.yaml."""
    print("Test: user_advice_keywords_yaml_structure")

    from libtbx.langchain.knowledge.yaml_loader import load_programs

    programs = load_programs()

    # Check autosol has keywords
    autosol = programs.get("phenix.autosol", {})
    assert "user_advice_keywords" in autosol, "autosol should have user_advice_keywords"

    keywords = autosol["user_advice_keywords"]
    assert isinstance(keywords, list), "keywords should be a list"
    assert "experimental phasing" in keywords, "Should include 'experimental phasing'"
    assert "sad phasing" in keywords, "Should include 'sad phasing'"
    assert "anomalous" in keywords, "Should include 'anomalous'"

    # Check phaser has keywords
    phaser = programs.get("phenix.phaser", {})
    assert "user_advice_keywords" in phaser, "phaser should have user_advice_keywords"
    assert "molecular replacement" in phaser["user_advice_keywords"]

    # Check predict_and_build has keywords
    pab = programs.get("phenix.predict_and_build", {})
    assert "user_advice_keywords" in pab, "predict_and_build should have keywords"
    assert "alphafold" in pab["user_advice_keywords"]

    print("  PASSED")


def test_user_advice_keywords_registry():
    """Test ProgramRegistry.get_user_advice_keywords() method."""
    print("Test: user_advice_keywords_registry")

    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry(use_yaml=True)

    # Test autosol keywords
    autosol_kw = registry.get_user_advice_keywords("phenix.autosol")
    assert len(autosol_kw) > 0, "autosol should have keywords"
    assert "anomalous" in autosol_kw, "Should include 'anomalous'"

    # Test phaser keywords
    phaser_kw = registry.get_user_advice_keywords("phenix.phaser")
    assert "molecular replacement" in phaser_kw, "Should include 'molecular replacement'"

    # Test program without keywords
    xtriage_kw = registry.get_user_advice_keywords("phenix.xtriage")
    assert xtriage_kw == [], "xtriage should have no keywords"

    print("  PASSED")


def test_user_advice_keyword_matching():
    """Test that user advice keywords affect program selection."""
    print("Test: user_advice_keyword_matching")

    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()

    workflow_state = {
        "state": "obtain_model",
        "experiment_type": "xray",
        "valid_programs": ["phenix.predict_and_build", "phenix.phaser", "phenix.autosol"],
        "categorized_files": {
            "mtz": ["data.mtz"],
            "sequence": ["seq.fa"],
        },
    }

    # Test "experimental phasing" advice - should select autosol
    intent1 = selector.select_next_action(
        workflow_state, {}, {},
        user_advice="Use experimental phasing for this dataset"
    )
    assert intent1["program"] == "phenix.autosol", \
        "Should select autosol for 'experimental phasing', got: %s" % intent1["program"]

    # Test "molecular replacement" advice - should select phaser
    intent2 = selector.select_next_action(
        workflow_state, {}, {},
        user_advice="Try molecular replacement"
    )
    assert intent2["program"] == "phenix.phaser", \
        "Should select phaser for 'molecular replacement', got: %s" % intent2["program"]

    # Test "alphafold" advice - should select predict_and_build
    intent3 = selector.select_next_action(
        workflow_state, {}, {},
        user_advice="Use alphafold prediction"
    )
    assert intent3["program"] == "phenix.predict_and_build", \
        "Should select predict_and_build for 'alphafold', got: %s" % intent3["program"]

    print("  PASSED")


def test_priority_when_yaml_structure():
    """Test that priority_when is properly defined in workflows.yaml."""
    print("Test: priority_when_yaml_structure")

    from libtbx.langchain.knowledge.yaml_loader import load_workflows

    workflows = load_workflows()

    # Check xray workflow
    xray = workflows.get("xray", {})
    phases = xray.get("phases", {})
    obtain_model = phases.get("obtain_model", {})
    programs = obtain_model.get("programs", [])

    # Find autosol entry
    autosol_entry = None
    for p in programs:
        if isinstance(p, dict) and p.get("program") == "phenix.autosol":
            autosol_entry = p
            break

    assert autosol_entry is not None, "autosol should be in obtain_model phase"
    assert "priority_when" in autosol_entry, "autosol should have priority_when"
    assert autosol_entry["priority_when"] == "strong_anomalous", \
        "autosol priority_when should be 'strong_anomalous'"

    print("  PASSED")


def test_priority_when_workflow_engine():
    """Test that priority_when affects program_priorities in workflow state."""
    print("Test: priority_when_workflow_engine")

    from libtbx.langchain.agent.workflow_engine import WorkflowEngine

    engine = WorkflowEngine()

    # Build context with strong anomalous signal
    # Need xtriage_done=True to be in obtain_model phase where autosol is valid
    files = {
        "mtz": ["data.mtz"],
        "sequence": ["seq.fa"],
    }
    history_info = {
        "strong_anomalous": True,
        "anomalous_measurability": 0.15,
        "xtriage_done": True,  # Must have completed analyze phase
    }

    workflow_state = engine.get_workflow_state("xray", files, history_info, {})

    # Check that program_priorities includes autosol
    priorities = workflow_state.get("program_priorities", [])
    assert "phenix.autosol" in priorities, \
        "autosol should be in program_priorities when strong_anomalous=True"

    # Test without strong anomalous
    history_info_no_anom = {
        "strong_anomalous": False,
        "anomalous_measurability": 0.02,
        "xtriage_done": True,
    }

    workflow_state2 = engine.get_workflow_state("xray", files, history_info_no_anom, {})
    priorities2 = workflow_state2.get("program_priorities", [])
    assert "phenix.autosol" not in priorities2, \
        "autosol should NOT be in program_priorities when strong_anomalous=False"

    print("  PASSED")


def test_priority_when_rules_selector():
    """Test that program_priorities affect program selection."""
    print("Test: priority_when_rules_selector")

    from libtbx.langchain.agent.rules_selector import RulesSelector

    selector = RulesSelector()

    # Workflow state with autosol as priority (strong anomalous detected)
    workflow_state_with_priority = {
        "state": "obtain_model",
        "experiment_type": "xray",
        "valid_programs": ["phenix.predict_and_build", "phenix.phaser", "phenix.autosol"],
        "program_priorities": ["phenix.autosol"],  # From priority_when
        "categorized_files": {
            "mtz": ["data.mtz"],
            "sequence": ["seq.fa"],
        },
    }

    intent = selector.select_next_action(workflow_state_with_priority, {}, {})

    assert intent["program"] == "phenix.autosol", \
        "Should select autosol when in program_priorities, got: %s" % intent["program"]

    # Without priority - should select predict_and_build (default preferred)
    workflow_state_no_priority = {
        "state": "obtain_model",
        "experiment_type": "xray",
        "valid_programs": ["phenix.predict_and_build", "phenix.phaser", "phenix.autosol"],
        "program_priorities": [],  # No priority
        "categorized_files": {
            "mtz": ["data.mtz"],
            "sequence": ["seq.fa"],
        },
    }

    intent2 = selector.select_next_action(workflow_state_no_priority, {}, {})

    # Should be predict_and_build (has sequence, default priority)
    assert intent2["program"] in ["phenix.predict_and_build", "phenix.phaser"], \
        "Should select predict_and_build or phaser without priority, got: %s" % intent2["program"]

    print("  PASSED")


def run_all_tests():
    """Run all YAML infrastructure tests."""
    print("=" * 60)
    print("YAML INFRASTRUCTURE TESTS")
    print("=" * 60)

    test_yaml_loading()
    test_program_registry()
    test_experiment_type_filtering()
    test_metric_thresholds()
    test_log_pattern_extraction()
    test_measurability_parsing()
    test_measurability_sanity_check()
    test_workflow_phases()
    test_command_with_strategy()
    test_workflow_engine()
    test_workflow_engine_refined_state()
    test_yaml_workflow_state_detection()
    test_metric_evaluator()
    test_metric_evaluator_trend()
    test_metrics_analyzer_yaml_mode()
    test_rules_selector()
    test_rules_selector_stop()

    # Step 1: Invariants tests
    print()
    print("-" * 60)
    print("STEP 1: INVARIANTS TESTS")
    print("-" * 60)
    test_invariants_yaml_structure()
    test_invariants_validate_and_fix()
    test_invariants_check_types()

    # Step 2: Input priorities tests
    print()
    print("-" * 60)
    print("STEP 2: INPUT PRIORITIES TESTS")
    print("-" * 60)
    test_input_priorities_yaml_structure()
    test_input_priorities_registry()
    test_input_priorities_selection()

    # Step 3: User advice keywords and priority_when tests
    print()
    print("-" * 60)
    print("STEP 3: USER ADVICE & PRIORITY_WHEN TESTS")
    print("-" * 60)
    test_user_advice_keywords_yaml_structure()
    test_user_advice_keywords_registry()
    test_user_advice_keyword_matching()
    test_priority_when_yaml_structure()
    test_priority_when_workflow_engine()
    test_priority_when_rules_selector()

    print("=" * 60)
    print("ALL YAML TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
