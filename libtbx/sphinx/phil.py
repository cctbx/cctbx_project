from __future__ import division

from docutils.parsers.rst import Directive
from docutils.parsers.rst import directives
from docutils import nodes


def setup(app):
  app.add_directive('phil', PhilDirective)


class PhilDirective(Directive):

  # this disables content in the directive
  has_content = False
  required_arguments = 1
  option_spec = {'expert-level': directives.nonnegative_int,
                 'attributes-level': directives.nonnegative_int}

  def run(self):
    from cStringIO import StringIO

    phil_include = self.arguments[0]
    expert_level = self.options.get('expert-level', None)
    attributes_level = self.options.get('attributes-level', 0)

    self.master_params = self._find_phil_scope(phil_include)
    s = StringIO()
    self.master_params.show(
      s, expert_level=expert_level, attributes_level=attributes_level)

    text = s.getvalue()
    text_nodes, messages = self.state.inline_text(text, 0)
    node = nodes.literal_block(text, '', *text_nodes)
    self.add_name(node)
    return [node] + messages

  def _find_phil_scope(self, phil_include):
    from types import ModuleType
    import libtbx.utils
    import libtbx.phil
    import iotbx.phil

    # common phil variable names
    search = [
      "", "master_phil_scope", "master_params", "master_phil",
      "master_params_str", "master_phil_str", "get_master_phil"]

    master_params = None
    for i in search:
      try:
        if i == "":
          import_path = phil_include
        else:
          import_path = "%s.%s" %(phil_include, i)
        print import_path
        master_params = libtbx.utils.import_python_object(
          import_path=import_path,
          error_prefix="",
          target_must_be="",
          where_str="").object
        # check that we haven't just imported a module
        if isinstance(master_params, ModuleType):
          continue
        break
      except Exception, e:
        print e
        pass

    # Check if the module attribute is a string, a scope, ...
    if isinstance(master_params, libtbx.phil.scope):
      pass
    elif isinstance(master_params, (str, unicode)):
      master_params = iotbx.phil.parse(master_params, process_includes=True)
    elif hasattr(master_params, '__call__'):
      master_params = master_params()
    else:
      pass

    if not master_params:
      raise Exception("No PHIL command found for %s." % phil_include)
    return master_params
