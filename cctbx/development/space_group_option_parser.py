import os, sys
import libtbx.option_parser
from cctbx.development import debug_utils
from cctbx import sgtbx

class space_group_processed_options(libtbx.option_parser.processed_options):

  def loop_over_space_groups(self, f, **kwds):
    kwds.update(self.options.__dict__)
    for sgi in self.space_group_info_list:
      f(space_group_info=sgi, **kwds)


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
    self.option('-v', '--verbose', type="int", default=0)
    self.option('--space_group_set',
                metavar=' ' + ' | '.join(self.space_group_sets))
    self.option('-s', '--n_scatterers', default=None)
    self.option(None, '--general_positions_only',
                action="store_true", default=False)

  def process(self, args):
    adjusted_args = []
    for arg in args:
      if (arg == "--Verbose"): arg = "--verbose=1"
      adjusted_args.append(arg)
    command_line = libtbx.option_parser.option_parser.process(
      self, adjusted_args)
    opts = command_line.options
    if opts.n_scatterers is not None:
      opts.n_scatterers = [ int(x) for x in opts.n_scatterers.split(',') ]
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
