# optik141 is the unmodified distribution, except for:
# mv lib optik
import sys, os
if (sys.version_info[0] < 3 or sys.version_info[1] < 3):
  sys.path.insert(0, os.path.join(os.environ["LIBTBX_DIST"], "optik141"))
  from optik import *
else:
  from optparse import *

class option_parser(OptionParser):

  def option(self, *args, **kw):
    self.add_option(apply(make_option, args, kw))
    return self

  def process(self, args):
    try: (options, args) = self.parse_args(args)
    except OptionError, e:
      print >> sys.stderr, e
      sys.exit(1)
    return processed_options(self, options, args)

class processed_options:

  def __init__(self, parser, options, args):
    self.parser = parser
    self.options = options
    self.args = args
