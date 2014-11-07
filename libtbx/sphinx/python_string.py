from __future__ import division

from docutils.parsers.rst import Directive
from docutils.parsers.rst import directives
from docutils import nodes


def setup(app):
  app.add_directive('python_string', PythonStringDirective)


class PythonStringDirective(Directive):

  # this disables content in the directive
  has_content = False
  required_arguments = 1

  def run(self):
    import libtbx.utils

    python_path = self.arguments[0]
    python_string = libtbx.utils.import_python_object(
      import_path=python_path,
      error_prefix="",
      target_must_be="",
      where_str="").object

    assert isinstance(python_string, basestring)


    from docutils import statemachine
    include_lines = statemachine.string2lines(python_string, 2,
                                              convert_whitespace=True)
    self.state_machine.insert_input(include_lines, python_path)
    return []

    text = python_string
    text_nodes, messages = self.state.inline_text(text, self.lineno)
    node = nodes.literal_block(text, '', *text_nodes)
    self.add_name(node)
    return [node] + messages
