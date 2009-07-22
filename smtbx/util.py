import os, sys
import libtbx.option_parser
from libtbx import group_args
from cctbx.development import debug_utils
from cctbx import sgtbx


class space_group_processed_options(libtbx.option_parser.processed_options):

  def loop_over_space_groups(self, f):
    for sgi in self.space_group_info_list:
      f(sgi, **self.options.__dict__)


class space_group_option_parser(libtbx.option_parser.option_parser):

  space_group_sets = ['all', 'chiral', 'all_settings', 'unusual_settings']
  processed_options_type = space_group_processed_options

  def __init__(self, cmd=None, usage_insertion='', description_extension=''):
    if cmd is None:
      frame = sys._getframe(1)
      cmd = frame.f_code.co_filename
      root, ext = os.path.splitext(os.path.basename(cmd))
      cmd = root
    libtbx.option_parser.option_parser.__init__(self,
      usage="\n\t%s %s --space_group_set=SET"
            "\n\t%s %s SPACE_GROUP_SYMBOL"
            "\n\t%s %s" % ((cmd, usage_insertion, )*3),
      description="The 3rd form selects a predefined set of representative "
                  "spacegroups (c.f. module cctbx.development.debug_utils)"
                  + description_extension
    )
    self.option('-v', '--verbose', default=0)
    self.option('--space_group_set',
                metavar=' ' + ' | '.join(self.space_group_sets))

  def process(self, args):
    command_line = libtbx.option_parser.option_parser.process(self, args)
    opts = command_line.options
    sg_set_flags = [ opts.space_group_set == x
                     for x in self.space_group_sets ]
    if sg_set_flags.count(True) == 0 and command_line.args:
      sg_symbols = command_line.args
    else:
      sg_symbols = debug_utils.get_test_space_group_symbols(*sg_set_flags)
    command_line.space_group_info_list = [
      sgtbx.space_group_info(symbol) for symbol in sg_symbols ]
    del command_line.options.space_group_set
    for attr in command_line.options.__dict__.keys():
      if getattr(command_line.options, attr) in ('', None):
        delattr(command_line.options, attr)
    command_line.args = ()
    return command_line


class test_case(object):

  class __metaclass__(type):
    def __init__(cls, classname, bases, classdict):
      exercises = []
      for base in bases:
        try:
          exercises.extend(base.exercises)
        except AttributeError:
          pass
      for name, attr in classdict.items():
        if callable(attr) and name.startswith('exercise'):
          exercises.append(attr)
      dsu = [ (ex.__name__, ex) for ex in exercises ]
      dsu.sort()
      cls.exercises = [ ex for foo, ex in dsu ]

  def run(cls, verbose=False, *args, **kwds):
    if verbose: print cls.__name__
    for exercise in cls.exercises:
      if verbose: print "\t%s ... " % exercise.__name__,
      o = cls(*args, **kwds)
      exercise(o)
      if verbose: print "OK"
  run = classmethod(run)
