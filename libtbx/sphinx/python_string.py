from __future__ import absolute_import, division, print_function

import docutils.statemachine
import libtbx.utils
from docutils.parsers.rst import Directive
from six import string_types

def setup(app):
  app.add_directive('python_string', PythonStringDirective)
  return {"parallel_read_safe": True}

class PythonStringDirective(Directive):
  # this disables content in the directive
  has_content = False
  required_arguments = 1

  def run(self):
    python_path = self.arguments[0]
    python_string = libtbx.utils.import_python_object(
      import_path=python_path,
      error_prefix="",
      target_must_be="",
      where_str="").object

    assert isinstance(python_string, string_types)

    include_lines = docutils.statemachine.string2lines(
        python_string, tab_width=2, convert_whitespace=True)
    self.state_machine.insert_input(include_lines, python_path)
    return []
