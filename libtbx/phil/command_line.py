import libtbx.phil
from libtbx.utils import Sorry, format_exception
import os
op = os.path

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

  def process_arg(self, arg):
    try:
      params = libtbx.phil.parse(
        input_string=arg,  # XXX was: arg.replace(r"\n", "\n")
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

  def process_args(self, args, custom_processor=None):
    user_phils = []
    for arg in args:
      if (len(arg.strip()) == 0):
        continue
      if (arg.startswith("--")):
        arg_work = arg[2:]
        if (arg_work.find("=") < 0): arg_work += " = True"
        user_phils.append(self.process_arg(arg=arg_work))
        continue
      if (op.isfile(arg) and op.getsize(arg) > 0):
        try: user_phils.append(libtbx.phil.parse(file_name=arg))
        except Exception: pass
        else: continue
      if (arg.find("=") >= 0):
        try: user_phils.append(self.process_arg(arg=arg))
        except Exception: pass
        else: continue
      if (custom_processor is not None) :
        result = custom_processor(arg=arg)
        if (isinstance(result, libtbx.phil.scope)) :
          user_phils.append(result)
          continue
        elif (result is not None) and (result != False) :
          continue
      if (op.isfile(arg)):
        libtbx.phil.parse(file_name=arg) # exception expected
        from libtbx.str_utils import show_string
        raise RuntimeError(
          'Programming error or highly unusual situation'
          ' (while processing %sargument %s).' % (
            self.argument_description, show_string(arg)))
      from libtbx.str_utils import show_string
      raise Sorry(
        "Uninterpretable %sargument: %s" % (
          self.argument_description, show_string(arg)))
    return user_phils

  def process(self, arg=None, args=None, custom_processor=None):
    assert [arg, args].count(None) == 1
    if (arg is not None):
      assert custom_processor is None
      return self.process_arg(arg=arg)
    return self.process_args(args=args, custom_processor=custom_processor)

  def process_and_fetch(self, args, custom_processor=None):
    if (isinstance(custom_processor, str)):
      assert custom_processor == "collect_remaining"
      remaining_args = []
      def custom_processor(arg):
        remaining_args.append(arg)
        return True
    else:
      remaining_args = None
    result = self.master_phil.fetch(
      sources=self.process(args=args, custom_processor=custom_processor))
    if (remaining_args is None):
      return result
    return result, remaining_args

class process(object):

  def __init__(self, args, master_string, parse=None):
    if (parse is None): parse = libtbx.phil.parse
    self.parse = parse
    self.master = self.parse(input_string=master_string, process_includes=True)
    self.work, self.remaining_args = argument_interpreter(
      master_phil=self.master) \
        .process_and_fetch(
          args=args,
          custom_processor="collect_remaining")

  def show(self, out=None):
    self.work.show(out=out)
    return self
