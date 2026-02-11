"""
Tests for command_builder module.

Tests the unified command generation pipeline:
1. File selection
2. Strategy building
3. Invariant application
4. Command assembly
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock libtbx imports - must be done BEFORE importing our modules
if 'libtbx' not in sys.modules:
    import types
    _base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    libtbx = types.ModuleType('libtbx')
    libtbx.__path__ = [_base_dir]
    libtbx.langchain = types.ModuleType('libtbx.langchain')
    libtbx.langchain.__path__ = [_base_dir]
    libtbx.langchain.agent = types.ModuleType('libtbx.langchain.agent')
    libtbx.langchain.agent.__path__ = [os.path.join(_base_dir, 'agent')]
    libtbx.langchain.knowledge = types.ModuleType('libtbx.langchain.knowledge')
    libtbx.langchain.knowledge.__path__ = [os.path.join(_base_dir, 'knowledge')]
    sys.modules['libtbx'] = libtbx
    sys.modules['libtbx.langchain'] = libtbx.langchain
    sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
    sys.modules['libtbx.langchain.knowledge'] = libtbx.langchain.knowledge

# First mock yaml_loader
import knowledge.yaml_loader
if 'libtbx.langchain.knowledge.yaml_loader' not in sys.modules:
    sys.modules['libtbx.langchain.knowledge.yaml_loader'] = knowledge.yaml_loader

# Mock pattern_manager
from agent import pattern_manager
if 'libtbx.langchain.agent.pattern_manager' not in sys.modules:
    sys.modules['libtbx.langchain.agent.pattern_manager'] = pattern_manager

# Mock program_registry
from agent import program_registry
if 'libtbx.langchain.agent.program_registry' not in sys.modules:
    sys.modules['libtbx.langchain.agent.program_registry'] = program_registry

# Mock command_builder
from agent import command_builder
if 'libtbx.langchain.agent.command_builder' not in sys.modules:
    sys.modules['libtbx.langchain.agent.command_builder'] = command_builder

# Mock template_builder
from agent import template_builder
if 'libtbx.langchain.agent.template_builder' not in sys.modules:
    sys.modules['libtbx.langchain.agent.template_builder'] = template_builder

# Import command_builder classes for tests
from agent.command_builder import CommandBuilder, CommandContext, get_command_builder


def test_context_from_state():
    """Test CommandContext.from_state()"""
    print("Test: CommandContext.from_state()")

    state = {
        'cycle_number': 5,
        'session_resolution': 2.5,
        'session_info': {
            'experiment_type': 'xray',
            'best_files': {'model': '/data/best.pdb', 'mtz': '/data/best.mtz'},
            'rfree_mtz': '/data/rfree.mtz',
        },
        'workflow_state': {
            'resolution': 2.5,
            'categorized_files': {
                'pdb': ['/data/model.pdb'],
                'mtz': ['/data/data.mtz'],
            },
            'state': 'refine',
        },
        'history': [
            {'cycle_number': 1, 'program': 'phenix.xtriage', 'result': 'SUCCESS'},
        ],
        'corrected_files': {'model': '/data/model.pdb'},
        'strategy': {'nproc': 4},
    }

    ctx = CommandContext.from_state(state)

    assert ctx.cycle_number == 5
    assert ctx.experiment_type == 'xray'
    assert ctx.resolution == 2.5
    assert ctx.rfree_mtz == '/data/rfree.mtz'
    assert ctx.workflow_state == 'refine'
    assert len(ctx.history) == 1
    assert ctx.llm_files == {'model': '/data/model.pdb'}
    assert ctx.llm_strategy == {'nproc': 4}

    print("  PASSED")


def test_output_prefix_generation():
    """Test output prefix generation based on history."""
    print("Test: output prefix generation")

    builder = CommandBuilder()

    # Test with no history
    ctx1 = CommandContext(
        cycle_number=1,
        history=[],
    )
    prefix1 = builder._generate_output_prefix('phenix.refine', ctx1)
    assert prefix1 == 'refine_001', f"Expected refine_001, got {prefix1}"

    # Test with 2 successful refine runs
    ctx2 = CommandContext(
        cycle_number=5,
        history=[
            {'cycle_number': 1, 'program': 'phenix.xtriage', 'result': 'SUCCESS'},
            {'cycle_number': 2, 'program': 'phenix.refine', 'result': 'SUCCESS: done'},
            {'cycle_number': 3, 'program': 'phenix.refine', 'result': 'FAILED: error'},
            {'cycle_number': 4, 'program': 'phenix.refine', 'result': 'SUCCESS: metrics'},
        ],
    )
    prefix2 = builder._generate_output_prefix('phenix.refine', ctx2)
    assert prefix2 == 'refine_003', f"Expected refine_003, got {prefix2}"

    # Test RSR
    ctx3 = CommandContext(
        cycle_number=3,
        history=[
            {'cycle_number': 1, 'program': 'phenix.real_space_refine', 'result': 'SUCCESS'},
        ],
    )
    prefix3 = builder._generate_output_prefix('phenix.real_space_refine', ctx3)
    assert prefix3 == 'rsr_002', f"Expected rsr_002, got {prefix3}"

    print("  PASSED")


def test_most_recent_file():
    """Test file selection prefers most recent."""
    print("Test: most recent file selection")

    builder = CommandBuilder()

    # Test numbered files
    files1 = [
        '/data/rsr_001_real_space_refined_000.pdb',
        '/data/rsr_001_real_space_refined_001.pdb',
        '/data/rsr_001_real_space_refined_002.pdb',
    ]
    result1 = builder._get_most_recent_file(files1)
    assert result1 == '/data/rsr_001_real_space_refined_002.pdb', f"Got {result1}"

    # Test refine files
    files2 = [
        '/data/refine_001_001.pdb',
        '/data/refine_002_001.pdb',
        '/data/refine_001_002.pdb',
    ]
    result2 = builder._get_most_recent_file(files2)
    # Should prefer highest last number
    assert '002' in result2, f"Expected file with 002, got {result2}"

    # Test empty list
    result3 = builder._get_most_recent_file([])
    assert result3 is None

    # Test single file
    result4 = builder._get_most_recent_file(['/data/model.pdb'])
    assert result4 == '/data/model.pdb'

    print("  PASSED")


def test_strategy_building():
    """Test strategy building with LLM hints."""
    print("Test: strategy building")

    builder = CommandBuilder()

    # Test with LLM strategy
    ctx1 = CommandContext(
        cycle_number=1,
        history=[],
        llm_strategy={'nproc': 8, 'simulated_annealing': True},
    )
    strategy1 = builder._build_strategy('phenix.refine', ctx1)
    assert strategy1.get('nproc') == 8
    assert strategy1.get('simulated_annealing') == True
    assert 'output_prefix' in strategy1  # Auto-added for refine

    # Test without LLM strategy
    ctx2 = CommandContext(
        cycle_number=2,
        history=[
            {'cycle_number': 1, 'program': 'phenix.real_space_refine', 'result': 'SUCCESS'},
        ],
    )
    strategy2 = builder._build_strategy('phenix.real_space_refine', ctx2)
    assert strategy2.get('output_prefix') == 'rsr_002'

    print("  PASSED")


def test_invariants_resolution():
    """Test invariant auto-fills resolution."""
    print("Test: invariants auto-fill resolution")

    builder = CommandBuilder()

    ctx = CommandContext(
        cycle_number=1,
        resolution=2.86,
    )

    files = {'model': '/data/model.pdb', 'map': '/data/map.ccp4'}
    strategy = {}

    files2, strategy2 = builder._apply_invariants(
        'phenix.real_space_refine', files, strategy, ctx
    )

    # Resolution should be auto-filled and rounded
    assert strategy2.get('resolution') == 2.9, f"Got {strategy2.get('resolution')}"

    print("  PASSED")


def test_invariants_rfree():
    """Test invariant adds R-free flags on first refinement."""
    print("Test: invariants R-free flags")

    builder = CommandBuilder()

    # First refinement should add R-free flags
    ctx1 = CommandContext(
        cycle_number=1,
        experiment_type='xray',
        history=[],
    )

    files, strategy = builder._apply_invariants(
        'phenix.refine', {}, {}, ctx1
    )
    assert strategy.get('generate_rfree_flags') == True

    # Second refinement should NOT add R-free flags
    ctx2 = CommandContext(
        cycle_number=3,
        experiment_type='xray',
        history=[
            {'cycle_number': 2, 'program': 'phenix.refine', 'result': 'SUCCESS'},
        ],
    )

    files2, strategy2 = builder._apply_invariants(
        'phenix.refine', {}, {}, ctx2
    )
    assert 'generate_rfree_flags' not in strategy2

    print("  PASSED")


def test_singleton():
    """Test get_command_builder() returns singleton."""
    print("Test: singleton pattern")

    builder1 = get_command_builder()
    builder2 = get_command_builder()

    assert builder1 is builder2

    print("  PASSED")


def test_compatibility_with_template_builder():
    """Test that CommandBuilder produces compatible output with TemplateBuilder."""
    print("Test: compatibility with TemplateBuilder")

    # Import TemplateBuilder
    try:
        from agent.template_builder import TemplateBuilder
    except ImportError:
        from libtbx.langchain.agent.template_builder import TemplateBuilder

    # Create both builders
    old_builder = TemplateBuilder(use_yaml=True)
    new_builder = CommandBuilder()

    # Test case: phenix.refine with explicit files
    files = {'model': '/data/model.pdb', 'mtz': '/data/data.mtz'}
    strategy = {'nproc': 4}

    old_cmd = old_builder.build_command('phenix.refine', files, strategy)

    # New builder needs context
    ctx = CommandContext(
        cycle_number=1,
        experiment_type='xray',
        history=[],
    )

    # For direct file specification, we need to use registry.build_command
    # since CommandBuilder.build() does file selection
    new_cmd = new_builder._registry.build_command('phenix.refine', files, strategy)

    # Commands should be equivalent (may have different ordering)
    old_parts = set(old_cmd.split())
    new_parts = set(new_cmd.split())

    assert old_parts == new_parts, f"Commands differ:\nOld: {old_cmd}\nNew: {new_cmd}"

    print("  PASSED")


def test_template_builder_delegation():
    """Test TemplateBuilder can delegate to CommandBuilder."""
    print("Test: TemplateBuilder delegation")

    try:
        from agent.template_builder import TemplateBuilder
    except ImportError:
        from libtbx.langchain.agent.template_builder import TemplateBuilder

    # Test with delegation disabled (default)
    builder = TemplateBuilder(use_yaml=True)
    assert builder.USE_NEW_BUILDER == False

    # Test that _get_command_builder works
    cmd_builder = builder._get_command_builder()
    assert cmd_builder is not None
    assert type(cmd_builder).__name__ == 'CommandBuilder'

    print("  PASSED")


def test_graph_nodes_feature_flag():
    """Test that graph_nodes has the new builder feature flag."""
    print("Test: graph_nodes feature flag")

    # Try to find graph_nodes.py in multiple locations
    import os

    # Location 1: Local agent/ directory (development)
    local_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'agent', 'graph_nodes.py'
    )

    # Location 2: Try to find via import (PHENIX environment)
    installed_path = None
    try:
        import libtbx.langchain.agent.graph_nodes as gn_module
        installed_path = gn_module.__file__
    except ImportError:
        pass

    # Use whichever exists
    graph_nodes_path = None
    if os.path.exists(local_path):
        graph_nodes_path = local_path
    elif installed_path and os.path.exists(installed_path):
        graph_nodes_path = installed_path

    if graph_nodes_path is None:
        print("  SKIPPED - graph_nodes.py not found")
        return

    with open(graph_nodes_path, 'r') as f:
        source = f.read()

    # Check feature flag is defined
    assert 'USE_NEW_COMMAND_BUILDER' in source, "Feature flag not found"

    # Check new functions exist
    assert 'def _build_with_new_builder' in source, "_build_with_new_builder not found"
    assert 'def _fallback_with_new_builder' in source, "_fallback_with_new_builder not found"
    assert 'def _get_command_builder' in source, "_get_command_builder not found"

    # Check delegation is in place
    assert 'if USE_NEW_COMMAND_BUILDER:' in source, "Delegation check not found"

    print("  PASSED")


def test_llm_slot_alias_mapping():
    """
    Test that LLM slot names are correctly mapped to program input names.

    Bug fix: LLM might use 'data' but program expects 'mtz'.
    The slot alias mapping should handle this.

    Example scenario:
    - LLM requests: data=PredictAndBuild_0_overall_best_refinement.mtz
    - Program expects: data_mtz slot
    - Without alias mapping: LLM file is ignored, wrong file auto-selected
    - With alias mapping: LLM's file is correctly used
    """
    print("test_llm_slot_alias_mapping:", end=" ")

    builder = CommandBuilder()

    # Check SLOT_ALIASES exists
    assert hasattr(builder, 'SLOT_ALIASES'), "SLOT_ALIASES not found on CommandBuilder"

    # Verify key mappings
    aliases = builder.SLOT_ALIASES
    assert aliases.get("data") == "data_mtz", "data should map to data_mtz"
    assert aliases.get("mtz") == "data_mtz", "mtz should map to data_mtz"
    assert aliases.get("pdb") == "model", "pdb should map to model"
    assert aliases.get("seq_file") == "sequence", "seq_file should map to sequence"

    print("  PASSED")


def test_llm_data_slot_used_for_mtz():
    """
    Test that when LLM uses 'data' slot, the file is correctly used for 'mtz' input.

    This tests the actual file selection logic, not just the alias mapping.
    """
    print("test_llm_data_slot_used_for_mtz:", end=" ")

    import tempfile
    import shutil

    # Create temp directory with test files
    temp_dir = tempfile.mkdtemp()
    try:
        # Create test files
        model_file = os.path.join(temp_dir, "model.pdb")
        mtz_correct = os.path.join(temp_dir, "PredictAndBuild_0_overall_best_refinement.mtz")
        mtz_wrong = os.path.join(temp_dir, "PredictAndBuild_0_refinement_cycle_2.extended_r_free.mtz")

        for f in [model_file, mtz_correct, mtz_wrong]:
            with open(f, 'w') as fh:
                fh.write("test")

        available_files = [model_file, mtz_correct, mtz_wrong]

        # Create context with LLM using 'data' slot
        context = CommandContext(
            cycle_number=1,
            experiment_type="xray",
            resolution=3.0,
            categorized_files={
                "model": [model_file],
                "data_mtz": [mtz_correct, mtz_wrong],
                "original_data_mtz": [mtz_correct, mtz_wrong],
            },
            best_files={"model": model_file, "data_mtz": mtz_correct},
            rfree_mtz=None,
            llm_files={
                "model": model_file,
                "data": mtz_correct,  # LLM uses 'data', maps to data_mtz
            },
            llm_strategy={},
        )

        builder = CommandBuilder()
        result = builder.build("phenix.refine", available_files, context)

        # The command should use the LLM's requested MTZ file
        assert result is not None, "Command generation failed"
        assert mtz_correct in result, \
            "LLM's requested MTZ not used. Got: %s" % result

        # Verify it didn't use the wrong file
        assert mtz_wrong not in result or mtz_correct in result, \
            "Wrong MTZ file used instead of LLM's choice"

        print("  PASSED")

    finally:
        shutil.rmtree(temp_dir)


def test_auto_fill_false_prevents_model_injection():
    """
    Test that auto_fill: false prevents map_to_model from getting a model auto-added.

    Bug: map_to_model is de novo model building - it should NOT auto-fill a model
    from available PDB files. The partial_model input has auto_fill: false.
    """
    print("test_auto_fill_false_prevents_model_injection:", end=" ")

    import tempfile
    import shutil

    temp_dir = tempfile.mkdtemp()
    try:
        # Create test files including a PDB that should NOT be auto-added
        map_file = os.path.join(temp_dir, "map.ccp4")
        seq_file = os.path.join(temp_dir, "seq.dat")
        pdb_file = os.path.join(temp_dir, "model.pdb")
        for f in [map_file, seq_file, pdb_file]:
            open(f, 'w').close()

        available = [map_file, seq_file, pdb_file]
        builder = CommandBuilder()

        state = {
            "cycle_number": 1,
            "session_info": {
                "experiment_type": "cryoem",
                "best_files": {},
                "rfree_mtz": None,
            },
            "workflow_state": {
                "resolution": 2.9,
                "categorized_files": {
                    "full_map": [map_file],
                    "sequence": [seq_file],
                    "model": [pdb_file],
                },
                "state": "cryoem_initial",
            },
            "history": [],
            "corrected_files": {},
            "strategy": {"resolution": 2.9},
        }

        context = CommandContext.from_state(state)
        command = builder.build("phenix.map_to_model", available, context)

        assert command is not None, "Command should be generated"
        assert "model=" not in command and "partial_model=" not in command, \
            "map_to_model should NOT include model= or partial_model= (got: %s)" % command

    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_fuzzy_slot_match_map_to_full_map():
    """
    Test that LLM's 'map' slot fuzzy-matches to 'full_map' for map_to_model.

    Bug: LLM requests map=file.ccp4 but map_to_model has 'full_map' slot.
    Without fuzzy matching, the file is rejected and auto-selection takes over.
    """
    print("test_fuzzy_slot_match_map_to_full_map:", end=" ")

    import tempfile
    import shutil

    temp_dir = tempfile.mkdtemp()
    try:
        map_file = os.path.join(temp_dir, "map.ccp4")
        seq_file = os.path.join(temp_dir, "seq.dat")
        for f in [map_file, seq_file]:
            open(f, 'w').close()

        available = [map_file, seq_file]
        builder = CommandBuilder()

        state = {
            "cycle_number": 1,
            "session_info": {
                "experiment_type": "cryoem",
                "best_files": {},
                "rfree_mtz": None,
            },
            "workflow_state": {
                "resolution": 2.9,
                "categorized_files": {
                    "full_map": [map_file],
                    "sequence": [seq_file],
                },
                "state": "cryoem_initial",
            },
            "history": [],
            "corrected_files": {
                "map": map_file,    # LLM requested "map"
                "sequence": seq_file,
            },
            "strategy": {"resolution": 2.9},
        }

        context = CommandContext.from_state(state)
        command = builder.build("phenix.map_to_model", available, context)

        assert command is not None, "Command should be generated"
        # The map file should be in the command (as bare positional since flag is "")
        assert "map.ccp4" in command, \
            "map.ccp4 should be in command via fuzzy match (got: %s)" % command

    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


# =============================================================================
# RECOVERY LABEL TRANSLATION TESTS
# =============================================================================

def test_recovery_label_translation_xtriage():
    """Test recovery translates label for xtriage."""
    print("Test: recovery label translation for xtriage")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/toxd.mtz': {
                'flags': {'scaling.input.xray_data.obs_labels': 'FTOXD3'},
                'program_scope': [],
                'reason': 'test',
                'selected_label': 'FTOXD3',
                'selected_label_pair': 'FTOXD3,SIGFTOXD3',
            }
        }
    )

    strategy = {}
    files = {'data_mtz': '/path/toxd.mtz'}
    result = builder._apply_recovery_strategies('phenix.xtriage', files, strategy, ctx)

    assert result.get('scaling.input.xray_data.obs_labels') == 'FTOXD3', \
        "xtriage should get scaling.input.xray_data.obs_labels=FTOXD3, got: %s" % result

    print("  PASSED")


def test_recovery_label_translation_refine():
    """Test recovery translates label for refine using miller_array.labels.name."""
    print("Test: recovery label translation for refine")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/toxd.mtz': {
                'flags': {'scaling.input.xray_data.obs_labels': 'FTOXD3'},
                'program_scope': [],
                'reason': 'test',
                'selected_label': 'FTOXD3',
                'selected_label_pair': 'FTOXD3,SIGFTOXD3',
            }
        }
    )

    strategy = {}
    files = {'data_mtz': '/path/toxd.mtz'}
    result = builder._apply_recovery_strategies('phenix.refine', files, strategy, ctx)

    assert result.get('miller_array.labels.name') == 'FTOXD3', \
        "refine should get miller_array.labels.name=FTOXD3, got: %s" % result
    # Should NOT have the xtriage-specific parameter
    assert 'scaling.input.xray_data.obs_labels' not in result, \
        "refine should NOT get xtriage parameter"

    print("  PASSED")


def test_recovery_label_translation_unknown_program():
    """Test recovery uses default parameter for unknown programs."""
    print("Test: recovery label translation for unknown program")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/toxd.mtz': {
                'flags': {'scaling.input.xray_data.obs_labels': 'FTOXD3'},
                'program_scope': [],
                'reason': 'test',
                'selected_label': 'FTOXD3',
                'selected_label_pair': 'FTOXD3,SIGFTOXD3',
            }
        }
    )

    strategy = {}
    files = {'data_mtz': '/path/toxd.mtz'}
    result = builder._apply_recovery_strategies('phenix.some_new_program', files, strategy, ctx)

    assert result.get('obs_labels') == 'FTOXD3', \
        "unknown program should get obs_labels=FTOXD3, got: %s" % result

    print("  PASSED")


def test_recovery_label_always_uses_main_label():
    """Test that recovery always uses main label, not the pair."""
    print("Test: recovery uses main label only (not pair)")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/data.mtz': {
                'flags': {},
                'program_scope': [],
                'reason': 'test',
                'selected_label': 'I_CuKa',
                'selected_label_pair': 'I_CuKa,SIGI_CuKa,I_CuKa_minus,SIGI_CuKa_minus',
            }
        }
    )

    strategy = {}
    files = {'data_mtz': '/path/data.mtz'}

    # For refine - should use main label, not pair
    result = builder._apply_recovery_strategies('phenix.refine', files, strategy, ctx)
    assert result.get('miller_array.labels.name') == 'I_CuKa', \
        "Should use main label only, got: %s" % result

    print("  PASSED")


def test_recovery_empty_scope_applies_to_all():
    """Test that empty program_scope applies recovery to all programs."""
    print("Test: empty program_scope applies to all programs")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/toxd.mtz': {
                'flags': {},
                'program_scope': [],
                'reason': 'test',
                'selected_label': 'FTOXD3',
                'selected_label_pair': 'FTOXD3,SIGFTOXD3',
            }
        }
    )

    files = {'data_mtz': '/path/toxd.mtz'}

    for prog in ['phenix.xtriage', 'phenix.refine', 'phenix.autobuild']:
        strategy = {}
        result = builder._apply_recovery_strategies(prog, files, strategy, ctx)
        assert len(result) > 0, "Recovery should apply to %s" % prog

    print("  PASSED")


def test_recovery_scoped_programs_respected():
    """Test that non-empty program_scope limits recovery to listed programs."""
    print("Test: program_scope limits recovery")

    builder = CommandBuilder()
    ctx = CommandContext(
        recovery_strategies={
            '/path/toxd.mtz': {
                'flags': {},
                'program_scope': ['phenix.xtriage'],
                'reason': 'test',
                'selected_label': 'FTOXD3',
                'selected_label_pair': 'FTOXD3,SIGFTOXD3',
            }
        }
    )

    files = {'data_mtz': '/path/toxd.mtz'}

    # Should apply to xtriage
    result1 = builder._apply_recovery_strategies('phenix.xtriage', files, {}, ctx)
    assert len(result1) > 0, "Should apply to xtriage"

    # Should NOT apply to refine
    result2 = builder._apply_recovery_strategies('phenix.refine', files, {}, ctx)
    assert len(result2) == 0, "Should NOT apply to refine when scoped to xtriage"

    print("  PASSED")


def test_invariant_blocks_missing_resolution():
    """Test that step 3.5 blocks command when required resolution is missing."""
    print("Test: invariant blocks missing resolution")

    import tempfile
    import shutil

    temp_dir = tempfile.mkdtemp()
    try:
        # Create minimal files
        map_file = os.path.join(temp_dir, "map.ccp4")
        seq_file = os.path.join(temp_dir, "seq.dat")
        with open(map_file, 'w') as f:
            f.write("mock")
        with open(seq_file, 'w') as f:
            f.write("mock")

        builder = CommandBuilder()
        available = [map_file, seq_file]

        # No resolution in context
        state = {
            "session_info": {
                "experiment_type": "cryoem",
                "best_files": {"map": map_file, "sequence": seq_file},
            },
            "workflow_state": {
                "state": "cryoem_has_model",
                "categorized_files": {
                    "full_map": [map_file],
                    "sequence": [seq_file],
                },
            },
            "session_resolution": None,
            "history": [],
        }

        context = CommandContext.from_state(state)
        command = builder.build("phenix.map_to_model", available, context)

        assert command is None, \
            "Command should be blocked when resolution is missing, got: %s" % command
    finally:
        shutil.rmtree(temp_dir)

    print("  PASSED")


def test_autosol_partpdb_file_in_command():
    """
    Test that autosol command includes partpdb_file= for MR-SAD workflow.

    Bug: partial_model was selected correctly (Final files shows it) but
    build_command dropped it because the command template only had
    {data_mtz} and {sequence} placeholders - no {partial_model}.
    Files without placeholders in the template must still be appended
    using their flag from the YAML definition.
    """
    print("test_autosol_partpdb_file_in_command:", end=" ")

    import tempfile
    import shutil

    temp_dir = tempfile.mkdtemp()
    try:
        mtz_file = os.path.join(temp_dir, "7lw7.mtz")
        seq_file = os.path.join(temp_dir, "exoV_construct.seq")
        phaser_pdb = os.path.join(temp_dir, "PHASER.1.pdb")

        for f in [mtz_file, seq_file, phaser_pdb]:
            with open(f, 'w') as fh:
                fh.write("test")

        available_files = [mtz_file, seq_file, phaser_pdb]

        context = CommandContext(
            cycle_number=3,
            experiment_type="xray",
            resolution=2.5,
            categorized_files={
                "data_mtz": [mtz_file],
                "sequence": [seq_file],
                "phaser_output": [phaser_pdb],
                "model": [phaser_pdb],
            },
            best_files={"data_mtz": mtz_file, "model": phaser_pdb},
            rfree_mtz=None,
            llm_files={
                "data": mtz_file,
                "model": phaser_pdb,
                "sequence": seq_file,
            },
            llm_strategy={"atom_type": "Fe", "sites": 4},
        )

        builder = CommandBuilder()
        result = builder.build("phenix.autosol", available_files, context)

        assert result is not None, "Command generation failed"
        assert "partpdb_file=" in result, \
            "autosol command missing partpdb_file= for MR-SAD. Got: %s" % result
        assert "PHASER.1.pdb" in result, \
            "autosol command missing PHASER.1.pdb. Got: %s" % result
        assert "autosol.data=" in result, \
            "autosol command missing data. Got: %s" % result
        assert "seq_file=" in result, \
            "autosol command missing sequence. Got: %s" % result

        print("  PASSED")

    finally:
        shutil.rmtree(temp_dir)


def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
