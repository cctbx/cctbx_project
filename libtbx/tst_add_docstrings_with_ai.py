from __future__ import division
import unittest
import os
import shutil
from libtbx.command_line.add_docstrings_with_ai import (
    strip_markdown_fences,
    file_needs_processing,
    get_line_ranges_with_doc_strings,
    remove_doc_string_line_ranges,
    check_for_only_changes_in_doc_strings,
)

class TestAddDocstringsWithAI(unittest.TestCase):

    def setUp(self):
        """Set up a temporary directory for test files."""
        self.test_dir = "test_temp_dir"
        os.makedirs(self.test_dir, exist_ok=True)

    def tearDown(self):
        """Clean up the temporary directory."""
        shutil.rmtree(self.test_dir)

    def test_strip_markdown_fences(self):
        """Test the removal of markdown code fences."""
        text_with_fences = "```python\nprint('hello')\n```"
        expected_text = "print('hello')"
        self.assertEqual(strip_markdown_fences(text_with_fences), expected_text)

        text_without_fences = "print('hello')"
        self.assertEqual(strip_markdown_fences(text_without_fences), text_without_fences)
        print("OK test_strip_markdown_fences")

    def test_file_needs_processing(self):
        """Test the logic for determining if a file needs processing."""
        file_with_defs = os.path.join(self.test_dir, "file_with_defs.py")
        with open(file_with_defs, "w") as f:
            f.write("def my_func():\n    pass\n\nclass MyClass:\n    pass")

        file_without_defs = os.path.join(self.test_dir, "file_without_defs.py")
        with open(file_without_defs, "w") as f:
            f.write("a = 10\nb = 20")

        empty_file = os.path.join(self.test_dir, "empty.py")
        open(empty_file, "w").close()

        self.assertTrue(file_needs_processing(file_with_defs))
        self.assertFalse(file_needs_processing(file_without_defs))
        self.assertFalse(file_needs_processing(empty_file))
        print("OK test_file_needs_processing")

    def test_get_line_ranges_with_doc_strings(self):
        """Test the identification of docstring line ranges."""
        lines = [
            '"""Module docstring."""',
            'import os',
            '',
            'class MyClass:',
            "    '''Class docstring.'''",
            '    def my_method(self):',
            '        """Method docstring."""',
            '        pass'
        ]
        expected_ranges = [[0, 0], [4, 4], [6, 6]]
        self.assertEqual(get_line_ranges_with_doc_strings(lines), expected_ranges)
        print("OK test_get_line_ranges_with_doc_strings")

    def test_remove_doc_string_line_ranges(self):
        """Test the removal of docstring lines."""
        lines = [
            '"""Module docstring."""',
            'import os',
            '',
            'class MyClass:',
            "    '''Class docstring.'''",
            '    def my_method(self):',
            '        """Method docstring."""',
            '        pass'
        ]
        expected_lines = [
            'import os',
            '',
            'class MyClass:',
            '    def my_method(self):',
            '        pass'
        ]
        self.assertEqual(remove_doc_string_line_ranges(lines, get_comment_and_empty_lines_too=False), expected_lines)
        print("OK test_remove_doc_string_line_ranges")

    def test_check_for_only_changes_in_doc_strings(self):
        """Test the function that verifies only docstrings have changed."""
        original_lines = [
            'def my_func():',
            '    pass'
        ]
        modified_lines = [
            'def my_func():',
            '    """This is a docstring."""',
            '    pass'
        ]
        self.assertIsNotNone(check_for_only_changes_in_doc_strings(original_lines, modified_lines))

        modified_lines_with_code_change = [
            'def my_func():',
            '    """This is a docstring."""',
            '    x = 1'
        ]
        self.assertIsNone(check_for_only_changes_in_doc_strings(original_lines, modified_lines_with_code_change))
        print("OK test_check_for_only_changes_in_doc_strings")

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
