"""
Unit tests for word_index_generator (generate an index from html directories).
This module tests visible text extraction, word tokenization, filtering, and index building.
"""

from __future__ import division
import unittest
import os
import re

OK = True
try:
  from word_index_generator import get_visible_text, tokenize, build_word_index
  from collections import defaultdict
  from bs4 import BeautifulSoup
  from nltk.corpus import stopwords
  import nltk
  test = (BeautifulSoup, stopwords, nltk)
except Exception as e:
  OK = False


class TestWordIndexGenerator(unittest.TestCase):
    """
    Unit tests for the word index generator script.

    These tests verify the behavior of:
    - get_visible_text(): Ensuring hidden content is removed and visible content retained.
    - tokenize(): Splitting input strings into normalized words.
    - stopword and pattern filtering: Ignoring trivial or unwanted words.
    - Integration: Index generation from actual file content.
    """

    def test_visible_text_extraction(self):
        """
        Ensure get_visible_text() extracts only visible content,
        excluding scripts, styles, comments, and hidden elements.
        """
        html = '''<!DOCTYPE html>
        <html><head><style>h1{display:none;}</style></head>
        <body>
        <h1 style="display:none;">Hidden Title</h1>
        <h2>Visible Heading</h2>
        <script>console.log('hi')</script>
        <p>This is a <b>test</b> page.</p>
        </body></html>'''
        visible = get_visible_text(html)
        self.assertIn("Visible Heading", visible)
        self.assertIn("This is a test page.", visible)
        self.assertNotIn("Hidden Title", visible)
        self.assertNotIn("console.log", visible)

    def test_tokenize_words(self):
        """
        Ensure tokenize() splits text into lowercase words,
        stripping out punctuation and handling alphanumerics.
        """
        text = "This is a test: only_words-123 _should be tokenized."
        tokens = tokenize(text)
        self.assertIn("this", tokens)
        self.assertIn("only_words", tokens)
        self.assertIn("123", tokens)
        self.assertNotIn(":", tokens)

    def test_exclusion_filtering(self):
        """
        Ensure filtering logic excludes stopwords, underscores,
        digits, and excluded patterns.
        """
        stop_words = {"this", "is", "a"}
        exclude_pattern = re.compile(r'^[a-zA-Z][0-9_]')
        sample_files = {"/tmp/sample.html"}
        word_index = defaultdict(set)
        words = ["this", "word", "x1", "_private", "z_test", "keepme"]
        for word in words:
            if (
                word not in stop_words and
                len(word) > 2 and
                not word[0].isdigit() and
                not word.startswith("_") and
                not exclude_pattern.match(word)
            ):
                word_index[word].update(sample_files)

        self.assertIn("word", word_index)
        self.assertIn("keepme", word_index)
        self.assertNotIn("this", word_index)
        self.assertNotIn("x1", word_index)
        self.assertNotIn("_private", word_index)
        self.assertNotIn("z_test", word_index)

    def test_build_word_index_from_sample_html(self):
        """
        Full integration test: index a small HTML file and validate words and filtering.
        """
        sample_dir = "test_html"
        os.makedirs(sample_dir, exist_ok=True)
        sample_path = os.path.join(sample_dir, "sample.html")
        with open(sample_path, "w", encoding="utf-8") as f:
            f.write("<html><body><p>Hello world from unittest!</p></body></html>")

        stop_words = {"from"}
        exclude_pattern = re.compile(r'^[a-zA-Z][0-9_]')
        word_index = build_word_index(sample_dir, stop_words, exclude_pattern, index_dir="index_files")

        self.assertIn("hello", word_index)
        self.assertIn("world", word_index)
        self.assertNotIn("from", word_index)

        os.remove(sample_path)
        os.rmdir(sample_dir)

if __name__ == '__main__':
  if OK:
    unittest.main()
  else:
    print("Skipping unit tests (requires nltk beautifulsoup4)")
