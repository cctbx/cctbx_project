import sys, os

sys.path.insert(0, os.path.join(
  os.path.dirname(os.path.dirname(__file__)), "optik151"))
from optik import *
make_option = Option

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

class processed_options(object):

  def __init__(self, parser, options, args):
    self.parser = parser
    self.options = options
    self.args = args
