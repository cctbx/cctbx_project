from __future__ import absolute_import, division, print_function
import os
from rstbx.phil import phil_preferences
from rstbx.phil.scope import scope
from libtbx.utils import Sorry

class Extract(object):
  def __init__(self,other):
    self.persist = other
  def __getattr__(self,item):
    if item=="persist":
      return super(Extract,self).__getattr__(item)
    return self.persist.commands.__getattribute__(item)
  def __setattr__(self,item,value):
    if item=="persist":
      super(Extract,self).__setattr__(item,value)
    else:
      self.persist.commands.__setattr__(item,value)

class RunTimePreferences(object):
  def __init__(self,scope = scope.value):

    "Parameters governed by PHIL module"
    self.phil_scope = phil_preferences.effective_param_generator.default(scope)
    self.new_scope_extract()
    self.command_extractor = Extract(self)

  def show(self):
    phil_preferences.effective_param_generator.show(self.commands)

  def rollback_dataset_preferences(self):
    self.phil_scope = phil_preferences.effective_param_generator.default(scope.value)
    self.new_scope_extract()

  def try_any_preferences_file(self,filename):
    if os.path.isfile(filename):
      import libtbx
      user_phil = libtbx.phil.parse(open(filename).read())
      self.phil_scope = self.phil_scope.fetch(source=user_phil)
      self.new_scope_extract()

  def try_dataset_preferences(self):
    filename = "dataset_preferences.py"
    self.try_any_preferences_file(filename)

  def merge_command_line(self,args):

    from libtbx.phil.command_line import argument_interpreter

    argument_interpreter = argument_interpreter(
      master_phil=phil_preferences.effective_param_generator.master(),
    )
    consume = []
    for arg in args:

      try:
        command_line_params = argument_interpreter.process(
          arg=arg
        )
        self.phil_scope = self.phil_scope.fetch(sources=[command_line_params,])
        consume.append(arg)

      except Sorry as e:
        pass

    for item in consume:
      args.remove(item)

    self.new_scope_extract()

  def new_scope_extract(self):
    self.commands = self.phil_scope.extract()
    if scope.value=="iotbx": # not libtbx
      phil_preferences.effective_param_generator.validation(self.commands)
