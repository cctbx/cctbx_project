"""
Unit tests for the langchain_tools library using the built-in unittest framework.
"""
from __future__ import division
import unittest
from libtbx.langchain import langchain_tools as lct
from langchain_core.documents import Document
from langchain_core.retrievers import BaseRetriever
from unittest.mock import MagicMock, patch

class TestLangchainTools(unittest.TestCase):
    """A test suite for the langchain_tools library."""

    def test_get_log_map_prompt(self):
        prompt = lct.get_log_map_prompt()
        self.assertIn("You are a highly skilled data extraction bot.", prompt.template)
        # --- ADDED CHECK for new text ---
        self.assertIn("overall_best", prompt.template)

    def test_get_log_combine_prompt(self):
        prompt = lct.get_log_combine_prompt()
        self.assertIn("You are an expert research scientist.", prompt.template)
        # --- ADDED CHECK for new text ---
        self.assertIn("overall_best", prompt.template)

    def test_get_docs_query_prompt(self):
        prompt = lct.get_docs_query_prompt()
        self.assertIn("You are an expert assistant for the Phenix software suite.", prompt.template)
        # --- ADDED CHECK for new text ---
        self.assertIn("Do not discuss limitations of the sources.", prompt.template)
        self.assertIn("Name the tools that are to be used", prompt.template)


    def test_get_log_analysis_prompt(self):
        prompt = lct.get_log_analysis_prompt()
        # --- UPDATED CHECKS for new text ---
        self.assertIn("You are expert in crystallography and cryo-EM.", prompt.template)
        self.assertIn("You are a Phenix power-user.", prompt.template)
        self.assertIn("Do not suggest depositing the model.", prompt.template)
        self.assertIn("Name the tools that are to be used", prompt.template)


    def test_custom_chunker(self):
        docs = [
            Document(page_content="Some text before. List of all available keywords is here. Some text after."),
            Document(page_content="Another document without the keyword.")
        ]
        chunks = lct._custom_chunker(docs)
        self.assertEqual(len(chunks), 3)
        self.assertIn("Some text before.", chunks[0].page_content)
        self.assertIn("List of all available keywords is here.", chunks[1].page_content)
        self.assertIn("Another document", chunks[2].page_content)

    def test_find_text_block(self):
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
        extracted = lct.find_text_block(log_text, "**3", end_text="**")
        self.assertIn("This is the text to extract. Line 2 of the text.", extracted)

    def test_get_processed_log_dict(self):
        log_text = """
**3. The Phenix Program**
phenix.real_space_refine
**8. Key Output Files and Evaluation:**
This is the summary.
        """
        processed_dict = lct.get_processed_log_dict(log_text)
        self.assertIn("phenix.real_space_refine", processed_dict['phenix_program'])
        self.assertIn("This is the summary.", processed_dict['summary'])

    @patch('langchain_tools.ChatGoogleGenerativeAI')
    def test_create_reranking_retriever(self, mock_llm):
        # Use a context manager to patch the environment variable
        with patch.dict('os.environ', {'COHERE_API_KEY': 'test_key'}):
            mock_vectorstore = MagicMock()
            mock_retriever = MagicMock(spec=BaseRetriever)
            mock_vectorstore.as_retriever.return_value = mock_retriever

            retriever = lct.create_reranking_retriever(mock_vectorstore, mock_llm)
            self.assertIsNotNone(retriever)

    @patch('langchain_tools.ChatGoogleGenerativeAI')
    def test_create_reranking_rag_chain(self, mock_llm):
        mock_retriever = MagicMock()
        prompt = lct.get_docs_query_prompt()
        chain = lct.create_reranking_rag_chain(mock_retriever, mock_llm, prompt)
        self.assertIsNotNone(chain)

    @patch('langchain_tools.ChatGoogleGenerativeAI')
    def test_create_log_analysis_chain(self, mock_llm):
        mock_retriever = MagicMock()
        prompt = lct.get_log_analysis_prompt()
        chain = lct.create_log_analysis_chain(mock_retriever, mock_llm, prompt)
        self.assertIsNotNone(chain)

# This block allows the script to be run directly
if __name__ == '__main__':
    unittest.main()
