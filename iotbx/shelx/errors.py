class error(RuntimeError):
  def __init__(self, msg, line, *args):
    msg = "ShelX: " + msg % args
    if (line is not None):
      msg += " at line %i" % (line+1)
    RuntimeError.__init__(self, msg)
