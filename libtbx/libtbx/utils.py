import sys

class Keep: pass

class UserError(Exception):

  # trick to get just "UserError" instead of "libtbx.utils.UserError"
  __module__ = "exceptions"

  def __init__(self, *args, **keyword_args):
    self.previous_tracebacklimit = getattr(sys, "tracebacklimit", None)
    sys.tracebacklimit = 0
    Exception.__init__(self, *args, **keyword_args)

  def __str__(self):
    self.reset_tracebacklimit()
    return Exception.__str__(self)

  def __del__(self):
    self.reset_tracebacklimit()

  def reset_tracebacklimit(self):
    if (hasattr(sys, "tracebacklimit")):
      if (self.previous_tracebacklimit is None):
        del sys.tracebacklimit
      else:
        sys.tracebacklimit = self.previous_tracebacklimit
