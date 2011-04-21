import libtbx.phil
from libtbx.utils import Sorry, format_exception

class argument_interpreter(object):

  def __init__(self,
        master_phil=None,
        home_scope=None,
        argument_description=None,
        master_params=None):
    if (argument_description is None):
      argument_description = "command line "
    assert [master_params, master_phil].count(None) == 1
    if (master_phil is None):
      import warnings
      warnings.warn(
        message='The "master_params" keyword argument name is deprecated.'
                ' Please use "master_phil" instead.',
        category=DeprecationWarning,
        stacklevel=2)
      master_phil = master_params
    self.master_phil = master_phil
    self.home_scope = home_scope
    self.argument_description = argument_description
    self.target_paths = None

  def get_path_score(self, source_path, target_path):
    i = target_path.find(source_path)
    if (i < 0): return 0
    if (i == 0 and len(source_path) == len(target_path)): return 8
    target_path_start_with_home_scope = False
    if (self.home_scope is not None):
      if (self.home_scope+"."+source_path == target_path): return 7
      if (target_path.startswith(self.home_scope+".")):
        if (target_path.endswith("."+source_path)): return 6
        if (target_path.endswith(source_path)): return 5
        target_path_start_with_home_scope = True
    if (target_path_start_with_home_scope): return 2
    if (target_path.endswith("."+source_path)): return 4
    if (target_path.endswith(source_path)): return 3
    return 1

  def process(self, arg):
    try:
      params = libtbx.phil.parse(
        input_string=arg.replace(r"\n", "\n"),
        source_info=self.argument_description+"argument")
    except RuntimeError:
      raise Sorry((
        'Error interpreting %sargument as parameter definition:\n'
        '  "%s"\n  %s') % (self.argument_description, arg, format_exception()))
    if (self.target_paths is None):
      self.target_paths = [object_locator.path
        for object_locator in self.master_phil.all_definitions()]
    source_definitions = params.all_definitions()
    complete_definitions = ""
    for object_locator in source_definitions:
      object = object_locator.object
      scores = [self.get_path_score(object_locator.path, target_path)
        for target_path in self.target_paths]
      max_score = max(scores)
      if (max_score == 0):
        raise Sorry("Unknown %sparameter definition: %s" % (
          self.argument_description, object.as_str().strip()))
      if (scores.count(max_score) > 1):
        error = ["Ambiguous parameter definition: %s" %
          object.as_str().strip()]
        error.append("Best matches:")
        for target_path,score in zip(self.target_paths, scores):
          if (score == max_score):
            error.append("  " + target_path)
        raise Sorry("\n".join(error))
      complete_definitions += object.customized_copy(
        name=self.target_paths[scores.index(max_score)]).as_str()
    if (complete_definitions == ""):
      raise Sorry(('%sparameter definition has no effect: "%s"' % (
        self.argument_description, arg)).capitalize())
    return libtbx.phil.parse(
      input_string=complete_definitions,
      source_info=self.argument_description+"argument")
