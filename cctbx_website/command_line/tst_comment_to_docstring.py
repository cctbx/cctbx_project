from __future__ import absolute_import, division, print_function
import unittest
import sys
import io
import textwrap

try:
    from comment_to_docstring import convert_comments_to_docstrings
except ImportError:
    print("Error: Make sure your script 'comment_to_docstring.py' is in the same directory.")
    sys.exit(1)

class TestCommentToDocstring(unittest.TestCase):
    """
    Final, corrected test suite for the comment_to_docstring.py script.
    """

    def assertCodeEqual(self, actual, expected):
        """
        Helper method to compare code blocks. It normalizes by dedenting the
        expected string and comparing the split lines of both strings.
        """
        self.assertEqual(
            actual.splitlines(),
            textwrap.dedent(expected).splitlines()
        )

    def test_simple_function_conversion(self):
        """Tests converting a single comment in a simple function."""
        source_code = """
            def my_func():
                # This is a test comment.
                a = 1
            """
        expected_output = '''
            def my_func():
                """ This is a test comment. """
                a = 1
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_multi_line_comment_in_function(self):
        """Tests converting a multi-line comment, ensuring line count is preserved."""
        source_code = """
            def my_func():
                # This is line one.
                # This is line two.
                a = 1
            """
        expected_output = '''
            def my_func():
                """ This is line one. This is line two.
                """
                a = 1
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_comment_with_blank_lines(self):
        """Tests that blank lines within a comment block are preserved."""
        source_code = """
            def my_func():
                # This is the first part.

                # This is the second part.
                b = 2
            """
        expected_output = '''
            def my_func():
                """ This is the first part. This is the second part.

                """
                b = 2
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_class_docstring_conversion(self):
        """Tests converting a comment into a class docstring."""
        source_code = """
            class MyClass:
              # This is a class comment.
              def __init__(self):
                pass
            """
        expected_output = '''
            class MyClass:
              """ This is a class comment. """
              def __init__(self):
                pass
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_no_change_if_docstring_exists(self):
        """Ensures functions with existing docstrings are not modified."""
        source_code = '''
            def my_func():
                """This is an existing docstring."""
                # This comment should be ignored.
                a = 1
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), source_code)

    def test_ignore_comment_after_code(self):
        """Ensures comments after the first line of code in a body are ignored."""
        source_code = """
            def my_func():
                a = 1
                # This comment describes the line above and should be ignored.
                return a
            """
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), textwrap.dedent(source_code))

    def test_promote_floating_module_comment(self):
        """Tests promoting a comment from between import statements to a module docstring."""
        # FIX: The source code now explicitly includes the leading blank line.
        source_code = """

            import a
            # This should be the module docstring.
            import b
            """
        # The expected output also has the leading blank line.
        expected_output = '''

            """ This should be the module docstring. """
            import a
            import b
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_shebang_and_coding_declaration_order(self):
        """Ensures the module docstring is placed correctly after special headers."""
        source_code = """#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This is the real module docstring.
import sys
"""
        expected_output = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" This is the real module docstring. """
import sys
'''
        self.assertCodeEqual(convert_comments_to_docstrings(source_code), expected_output)

    def test_multi_line_function_signature(self):
        """Tests the script with a multi-line function definition."""
        source_code = """
            class MyClass:
                def set_bounds(self, lower_bounds=None, upper_bounds=None,
                               map_data=None, box_size=None):
                    # This comment should be converted.
                    lower_bounds = list(lower_bounds)
            """
        expected_output = '''
            class MyClass:
                def set_bounds(self, lower_bounds=None, upper_bounds=None,
                               map_data=None, box_size=None):
                    """ This comment should be converted. """
                    lower_bounds = list(lower_bounds)
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_comment_ending_in_quote(self):
        """Tests the fix for comments ending in a quote to prevent invalid syntax."""
        source_code = """
            def my_func():
                # print "hello world"
                pass
            """
        expected_output = '''
            def my_func():
                """ print "hello world" """
                pass
            '''
        self.assertCodeEqual(convert_comments_to_docstrings(textwrap.dedent(source_code)), expected_output)

    def test_syntax_error_input(self):
        """Tests that the script handles a syntax error gracefully."""
        source_code = "def my_func()\n a = 1"

        captured_stderr = io.StringIO()
        original_stderr = sys.stderr
        sys.stderr = captured_stderr
        try:
            result = convert_comments_to_docstrings(source_code)
            self.assertEqual(result, source_code)
            self.assertIn("Error: Could not parse source.", captured_stderr.getvalue())
        finally:
            sys.stderr = original_stderr

if __name__ == '__main__':
    unittest.main(verbosity=2)
