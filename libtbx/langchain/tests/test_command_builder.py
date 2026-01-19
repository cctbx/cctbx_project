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
import tempfile

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock libtbx imports - must be done BEFORE importing our modules
import types
libtbx = types.ModuleType('libtbx')
libtbx.langchain = types.ModuleType('libtbx.langchain')
libtbx.langchain.agent = types.ModuleType('libtbx.langchain.agent')
libtbx.langchain.knowledge = types.ModuleType('libtbx.langchain.knowledge')
sys.modules['libtbx'] = libtbx
sys.modules['libtbx.langchain'] = libtbx.langchain
sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
sys.modules['libtbx.langchain.knowledge'] = libtbx.langchain.knowledge

# First mock yaml_loader
from knowledge.yaml_loader import get_program
import knowledge.yaml_loader
libtbx.langchain.knowledge.yaml_loader = knowledge.yaml_loader
sys.modules['libtbx.langchain.knowledge.yaml_loader'] = knowledge.yaml_loader

# Mock program_registry
from agent import program_registry
libtbx.langchain.agent.program_registry = program_registry
sys.modules['libtbx.langchain.agent.program_registry'] = program_registry

# Mock command_builder
from agent import command_builder
libtbx.langchain.agent.command_builder = command_builder
sys.modules['libtbx.langchain.agent.command_builder'] = command_builder

# Mock template_builder
from agent import template_builder
libtbx.langchain.agent.template_builder = template_builder
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


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Running command_builder tests")
    print("=" * 60)
    
    test_context_from_state()
    test_output_prefix_generation()
    test_most_recent_file()
    test_strategy_building()
    test_invariants_resolution()
    test_invariants_rfree()
    test_singleton()
    test_compatibility_with_template_builder()
    test_template_builder_delegation()
    test_graph_nodes_feature_flag()
    
    print("=" * 60)
    print("All tests PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
