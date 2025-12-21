"""
Comprehensive unit tests for the libtbx.langchain modules.

Tests cover:
- Core module (LLM setup, types)
- Analysis module (summarization, state extraction)
- RAG module (document loading, retrieval)
- Validation module (command validation)
- Knowledge module (programs, prompts)
- Agent module (memory, planning)
- Utils module (text processing)
- Backward compatibility (langchain_tools shim)

Run with: phenix.python tst_langchain_tools.py
Or: python -m pytest tst_langchain_tools.py -v
"""
from __future__ import division
import unittest
import os
import json
import tempfile
from unittest.mock import MagicMock, patch, AsyncMock

from langchain_core.documents import Document
from langchain_core.prompts import PromptTemplate


class TestCoreModule(unittest.TestCase):
    """Tests for libtbx.langchain.core"""

    def test_imports(self):
        """Test that core module imports work."""
        from libtbx.langchain.core import get_llm_and_embeddings
        self.assertIsNotNone(get_llm_and_embeddings)

    def test_types_agent_plan(self):
        """Test AgentPlan dataclass."""
        from libtbx.langchain.core.types import AgentPlan
        plan = AgentPlan(
            program="phenix.refine",
            strategy="Run refinement",
            commands=["phenix.refine model.pdb data.mtz"],
            reasoning="Model needs refinement"
        )
        self.assertEqual(plan.program, "phenix.refine")
        self.assertEqual(len(plan.commands), 1)

    def test_types_validation_result(self):
        """Test ValidationResult dataclass."""
        from libtbx.langchain.core.types import ValidationResult
        result = ValidationResult(
            is_valid=False,
            error_message="Invalid parameter",
            original_command="phenix.refine bad_param=True",
            fixed_command="phenix.refine good_param=True"
        )
        self.assertFalse(result.is_valid)
        self.assertIn("Invalid", result.error_message)


class TestAnalysisModule(unittest.TestCase):
    """Tests for libtbx.langchain.analysis"""

    def test_imports(self):
        """Test that analysis module imports work."""
        from libtbx.langchain.analysis import (
            summarize_log_text,
            analyze_log_summary,
            extract_project_state_updates,
            get_log_info,
        )
        self.assertIsNotNone(summarize_log_text)
        self.assertIsNotNone(analyze_log_summary)

    def test_get_log_map_prompt(self):
        """Test log map prompt content."""
        from libtbx.langchain.analysis import get_log_map_prompt
        prompt = get_log_map_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("data extraction bot", prompt.template)
        self.assertIn("overall_best", prompt.template)

    def test_get_log_combine_prompt(self):
        """Test log combine prompt content."""
        from libtbx.langchain.analysis import get_log_combine_prompt
        prompt = get_log_combine_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("expert research scientist", prompt.template)

    def test_get_log_analysis_prompt(self):
        """Test log analysis prompt content."""
        from libtbx.langchain.analysis import get_log_analysis_prompt
        prompt = get_log_analysis_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("crystallography and cryo-EM", prompt.template)
        self.assertIn("Phenix power-user", prompt.template)

    def test_get_chunk_size_google(self):
        """Test chunk size for Google provider."""
        from libtbx.langchain.analysis import get_chunk_size
        chunk_size, chunk_overlap = get_chunk_size('google')
        self.assertEqual(chunk_size, 750000)
        self.assertEqual(chunk_overlap, 50000)

    def test_get_chunk_size_openai(self):
        """Test chunk size for OpenAI provider."""
        from libtbx.langchain.analysis import get_chunk_size
        chunk_size, chunk_overlap = get_chunk_size('openai')
        self.assertEqual(chunk_size, 100000)
        self.assertEqual(chunk_overlap, 10000)


class TestRAGModule(unittest.TestCase):
    """Tests for libtbx.langchain.rag"""

    def test_imports(self):
        """Test that RAG module imports work."""
        from libtbx.langchain.rag import (
            load_all_docs_from_folder,
            load_specific_docs,
            create_and_persist_db,
            load_persistent_db,
            create_reranking_retriever,
        )
        self.assertIsNotNone(load_all_docs_from_folder)
        self.assertIsNotNone(load_persistent_db)

    def test_phenix_html_loader(self):
        """Test custom HTML loader."""
        from libtbx.langchain.rag import PhenixHTMLLoader

        # Create a temporary HTML file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False) as f:
            f.write("<html><body><p>Test content</p><script>ignore</script></body></html>")
            temp_path = f.name

        try:
            loader = PhenixHTMLLoader(temp_path)
            docs = loader.load()
            self.assertEqual(len(docs), 1)
            self.assertIn("Test content", docs[0].page_content)
            self.assertNotIn("ignore", docs[0].page_content)  # Script should be removed
        finally:
            os.unlink(temp_path)

    def test_custom_chunker(self):
        """Test document chunker with keyword detection."""
        from libtbx.langchain.rag.document_loader import _custom_chunker

        docs = [
            Document(
                page_content="Some text before. List of all available keywords is here. Some text after.",
                metadata={'source': 'test.txt'}
            ),
            Document(
                page_content="Another document without the keyword.",
                metadata={'source': 'test2.txt'}
            )
        ]
        chunks = _custom_chunker(docs)

        # Should have chunks from both documents
        self.assertGreater(len(chunks), 0)
        # Each chunk should have context header
        for chunk in chunks:
            self.assertIn("** DOCUMENT:", chunk.page_content)


class TestValidationModule(unittest.TestCase):
    """Tests for libtbx.langchain.validation"""

    def test_imports(self):
        """Test that validation module imports work."""
        from libtbx.langchain.validation import (
            validate_phenix_command,
            fix_command_syntax,
            mtz_has_rfree_flags,
            get_validator,
        )
        self.assertIsNotNone(validate_phenix_command)
        self.assertIsNotNone(get_validator)

    def test_get_validator_refine(self):
        """Test getting validator for phenix.refine."""
        from libtbx.langchain.validation import get_validator
        validator = get_validator('phenix.refine')
        self.assertIsNotNone(validator)
        self.assertEqual(validator.program_name, 'phenix.refine')

    def test_get_validator_unknown(self):
        """Test getting validator for unknown program returns None."""
        from libtbx.langchain.validation import get_validator
        validator = get_validator('phenix.nonexistent')
        self.assertIsNone(validator)

    def test_get_all_validators(self):
        """Test getting all registered validators."""
        from libtbx.langchain.validation import get_all_validators
        validators = get_all_validators()
        self.assertIsInstance(validators, dict)
        self.assertIn('phenix.refine', validators)

    def test_refine_validator_common_errors(self):
        """Test RefineValidator has common error patterns."""
        from libtbx.langchain.validation import RefineValidator
        validator = RefineValidator()
        errors = validator.get_common_errors()
        self.assertIsInstance(errors, dict)
        # Should have twin_law error pattern
        self.assertTrue(any('twin_law' in key for key in errors.keys()))


class TestKnowledgeModule(unittest.TestCase):
    """Tests for libtbx.langchain.knowledge"""

    def test_imports(self):
        """Test that knowledge module imports work."""
        from libtbx.langchain.knowledge import (
            get_phenix_program_list,
            get_keywords_as_phil_string,
            get_strategic_planning_prompt,
            get_command_writer_prompt,
        )
        self.assertIsNotNone(get_phenix_program_list)
        self.assertIsNotNone(get_strategic_planning_prompt)

    def test_get_phenix_program_list(self):
        """Test that program list contains expected programs."""
        from libtbx.langchain.knowledge import get_phenix_program_list
        programs = get_phenix_program_list()
        self.assertIsInstance(programs, list)
        self.assertGreater(len(programs), 0)
        # All should start with phenix.
        for prog in programs:
            self.assertTrue(prog.startswith('phenix.'))

    def test_get_strategic_planning_prompt(self):
        """Test strategic planning prompt content."""
        from libtbx.langchain.knowledge import get_strategic_planning_prompt
        prompt = get_strategic_planning_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("Lead Crystallographer", prompt.template)
        self.assertIn("FORBIDDEN PROGRAMS", prompt.template)

    def test_get_command_writer_prompt(self):
        """Test command writer prompt content."""
        from libtbx.langchain.knowledge import get_command_writer_prompt
        prompt = get_command_writer_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("Command-Line Expert", prompt.template)

    def test_get_docs_query_prompt(self):
        """Test docs query prompt content."""
        from libtbx.langchain.knowledge import get_docs_query_prompt
        prompt = get_docs_query_prompt()
        self.assertIsInstance(prompt, PromptTemplate)
        self.assertIn("expert assistant for the Phenix software suite", prompt.template)


class TestAgentModule(unittest.TestCase):
    """Tests for libtbx.langchain.agent"""

    def test_imports(self):
        """Test that agent module imports work."""
        from libtbx.langchain.agent import (
            generate_next_move,
            get_run_history,
            load_learned_memory,
            save_learned_memory,
            get_memory_file_path,
        )
        self.assertIsNotNone(generate_next_move)
        self.assertIsNotNone(load_learned_memory)

    def test_get_memory_file_path(self):
        """Test memory file path generation."""
        from libtbx.langchain.agent import get_memory_file_path
        path = get_memory_file_path("./test_db")
        self.assertIn("phenix_learned_memory.json", path)

    def test_load_save_learned_memory(self):
        """Test loading and saving learned memory."""
        from libtbx.langchain.agent import load_learned_memory, save_learned_memory

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_path = f.name

        try:
            # Test saving
            test_memory = {"phenix.refine": ["Use full PHIL paths"]}
            save_learned_memory(test_memory, temp_path)

            # Test loading
            loaded = load_learned_memory(temp_path)
            self.assertEqual(loaded, test_memory)
        finally:
            os.unlink(temp_path)

    def test_load_nonexistent_memory(self):
        """Test loading memory from nonexistent file returns empty dict."""
        from libtbx.langchain.agent import load_learned_memory
        result = load_learned_memory("/nonexistent/path/memory.json")
        self.assertEqual(result, {})

    def test_get_run_history_empty_dir(self):
        """Test get_run_history with nonexistent directory."""
        from libtbx.langchain.agent import get_run_history
        result = get_run_history("/nonexistent/directory")
        self.assertEqual(result, [])

    def test_extract_output_files(self):
        """Test extracting output files from summary."""
        from libtbx.langchain.agent import extract_output_files

        summary = """
        **Key Output Files:**
        model_refined.pdb, data_free.mtz
        **Next Section**
        """
        files = extract_output_files(summary)
        self.assertIn("model_refined.pdb", files)
        self.assertIn("data_free.mtz", files)


class TestUtilsModule(unittest.TestCase):
    """Tests for libtbx.langchain.utils"""

    def test_imports(self):
        """Test that utils module imports work."""
        from libtbx.langchain.utils import (
            find_text_block,
            get_processed_log_dict,
            query_docs,
        )
        self.assertIsNotNone(find_text_block)
        self.assertIsNotNone(query_docs)

    def test_find_text_block(self):
        """Test extracting text blocks from logs."""
        from libtbx.langchain.utils import find_text_block

        log_text = """
**1. Some header**
Some text
**2. Another header**
More text
And more
**3. The header we want**
This is the text to extract.
Line 2 of the text.
**4. Final header**
End of text.
        """
        extracted = find_text_block(log_text, "**3", end_text="**")
        self.assertIn("This is the text to extract", extracted)
        self.assertIn("Line 2 of the text", extracted)

    def test_find_text_block_not_found(self):
        """Test find_text_block when target not found."""
        from libtbx.langchain.utils import find_text_block
        result = find_text_block("Some text", "**NOTFOUND**")
        self.assertEqual(result, "")

    def test_find_text_block_none_input(self):
        """Test find_text_block with None input."""
        from libtbx.langchain.utils import find_text_block
        result = find_text_block(None, "**3")
        self.assertEqual(result, "")

    def test_get_processed_log_dict(self):
        """Test extracting program info from log."""
        from libtbx.langchain.utils import get_processed_log_dict

        log_text = """
**3. The Phenix Program**
phenix.real_space_refine
**8. Key Output Files and Evaluation:**
This is the summary.
        """
        processed_dict = get_processed_log_dict(log_text)
        self.assertIn("phenix_program", processed_dict)
        self.assertIn("phenix.real_space_refine", processed_dict['phenix_program'])


class TestBackwardCompatibility(unittest.TestCase):
    """Tests for backward compatibility with langchain_tools shim."""

    def test_import_as_lct(self):
        """Test importing as langchain_tools (legacy style)."""
        from libtbx.langchain import langchain_tools as lct
        self.assertIsNotNone(lct)

    def test_all_functions_available(self):
        """Test that all expected functions are available via shim."""
        from libtbx.langchain import langchain_tools as lct

        # Core
        self.assertTrue(hasattr(lct, 'get_llm_and_embeddings'))

        # Analysis
        self.assertTrue(hasattr(lct, 'summarize_log_text'))
        self.assertTrue(hasattr(lct, 'analyze_log_summary'))
        self.assertTrue(hasattr(lct, 'extract_project_state_updates'))
        self.assertTrue(hasattr(lct, 'get_log_info'))

        # RAG
        self.assertTrue(hasattr(lct, 'load_all_docs_from_folder'))
        self.assertTrue(hasattr(lct, 'create_and_persist_db'))
        self.assertTrue(hasattr(lct, 'load_persistent_db'))
        self.assertTrue(hasattr(lct, 'create_reranking_retriever'))

        # Validation
        self.assertTrue(hasattr(lct, 'validate_phenix_command'))
        self.assertTrue(hasattr(lct, 'mtz_has_rfree_flags'))

        # Knowledge
        self.assertTrue(hasattr(lct, 'get_phenix_program_list'))
        self.assertTrue(hasattr(lct, 'get_strategic_planning_prompt'))

        # Agent
        self.assertTrue(hasattr(lct, 'generate_next_move'))
        self.assertTrue(hasattr(lct, 'get_run_history'))

        # Utils
        self.assertTrue(hasattr(lct, 'find_text_block'))
        self.assertTrue(hasattr(lct, 'query_docs'))

    def test_functions_are_same(self):
        """Test that shim functions point to the same implementations."""
        from libtbx.langchain import langchain_tools as lct
        from libtbx.langchain.core import get_llm_and_embeddings
        from libtbx.langchain.analysis import summarize_log_text
        from libtbx.langchain.agent import generate_next_move

        self.assertIs(lct.get_llm_and_embeddings, get_llm_and_embeddings)
        self.assertIs(lct.summarize_log_text, summarize_log_text)
        self.assertIs(lct.generate_next_move, generate_next_move)


class TestIntegration(unittest.TestCase):
    """Integration tests that verify modules work together."""

    def test_prompts_have_required_variables(self):
        """Test that all prompts have their required input variables."""
        from libtbx.langchain.analysis import (
            get_log_map_prompt,
            get_log_combine_prompt,
            get_log_analysis_prompt,
        )
        from libtbx.langchain.knowledge import (
            get_strategic_planning_prompt,
            get_command_writer_prompt,
        )

        # Log map prompt needs 'text'
        prompt = get_log_map_prompt()
        self.assertIn('text', prompt.input_variables)

        # Log combine prompt needs 'context'
        prompt = get_log_combine_prompt()
        self.assertIn('context', prompt.input_variables)

        # Analysis prompt needs 'context' and 'log_summary'
        prompt = get_log_analysis_prompt()
        self.assertIn('context', prompt.input_variables)
        self.assertIn('log_summary', prompt.input_variables)

        # Strategic planning prompt
        prompt = get_strategic_planning_prompt()
        self.assertIn('num_runs', prompt.input_variables)
        self.assertIn('history', prompt.input_variables)

        # Command writer prompt
        prompt = get_command_writer_prompt()
        self.assertIn('program', prompt.input_variables)
        self.assertIn('valid_keywords', prompt.input_variables)


# Run tests
if __name__ == '__main__':
    # Run with verbosity
    unittest.main(verbosity=2)
