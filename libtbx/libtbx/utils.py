import sys

class UserError(Exception):

  def __init__(self, *args, **keyword_args):
    sys.tracebacklimit = 0
    Exception.__init__(self, *args, **keyword_args)
