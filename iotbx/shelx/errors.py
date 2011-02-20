class error(RuntimeError):
  def __init__(self, msg, line, *args):
    RuntimeError.__init__(self, (('ShelX: ' + msg) % args)
                          + " at line %i" % (line+1))
