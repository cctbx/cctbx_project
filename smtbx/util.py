import os, sys
import libtbx.option_parser
from libtbx import group_args
from cctbx.development import debug_utils
from cctbx import sgtbx

class space_group_option_parser(libtbx.option_parser.option_parser):

  space_group_sets = ['all', 'chiral', 'all_settings', 'unusual_settings']

  def __init__(self, cmd, usage_insertion, description_extension):
    cmd = os.path.basename(cmd)
    libtbx.option_parser.option_parser.__init__(self,
      usage="   %s %s --space_group_set=SET"
            "or %s %s SPACE_GROUP_SYMBOL"
            "or %s %s" % ((cmd, usage_insertion, )*3),
      description="The 3rd form selects a predefined set of representative "
                  "spacegroups (c.f. module cctbx.development.debug_utils)"
                  + description_extension
    )
    self.option('-v', '--verbose')
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

  def exercise(self, exercise_type, args):
    command_line = self.process(args)
    for sgi in command_line.space_group_info_list:
      e = exercise_type(space_group_info=sgi, options=command_line.options)
      e.run()


class test_case(object):

  def run(cls, verbose=False):
    dsu = []
    for name, attr in cls.__dict__.iteritems():
      if not callable(attr) or not name.startswith('exercise'): continue
      dsu.append((name, attr))
    dsu.sort()
    exercises = [ attr for foo, attr in dsu ]
    if verbose: print cls.__name__
    for exercise in exercises:
      if verbose: print "\t%s ... " % exercise.__name__,
      o = cls()
      exercise(o)
      if verbose: print "OK"
  run = classmethod(run)
