"""
Generic job scheduling
"""

from __future__ import absolute_import, division, print_function


class SchedulingError(Exception):
  """
  Package exception
  """


class SetupError(SchedulingError):
  """
  Error from the software environment
  """


class identifier(object):
  """
  Job identifier
  """

  def __init__(self):

    self.jobid = id( self )


  def __eq__(self, other):

    return self.jobid == other.jobid


  def __ne__(self, other):

    return not ( self == other )


  def __hash__(self):

    return hash( self.jobid )


  def __str__(self):

    return "identifier( jobid = %s )" % self.jobid


  def __repr__(self):

    return str( self )


class ignore(object):
  """
  Ignore extra stacktrace information

  Note this a singleton
  """

  @staticmethod
  def enter():

    pass


  @staticmethod
  def exit(exc_type, exc_val, exc_tb):

    if exc_type is None:
      from libtbx.scheduling import stacktrace
      stacktrace.cleanup()


class excepthook(object):
  """
  Use custom sys.excepthook to print error, replace when finished

  Note this is a singleton
  """

  @staticmethod
  def enter():

    pass


  @staticmethod
  def exit(exc_type, exc_val, exc_tb):

    from libtbx.scheduling import stacktrace

    if exc_type is not None:
      stacktrace.enable()

    else:
      stacktrace.cleanup()


class decorate(object):
  """
  Append stacktrace to exception message

  Note this is a singleton
  """

  @staticmethod
  def enter():

    pass


  @staticmethod
  def exit(exc_type, exc_val, exc_tb):

    from libtbx.scheduling import stacktrace

    if exc_type is not None:
      data = stacktrace.exc_info()

      if exc_val is data[0]:
        message = "%s. Embedded exception follows:\n%s" % (
          exc_val,
          "".join( data[1] ).rstrip(),
          )
        exc_val.args = ( message, ) + exc_val.args[1:]

    stacktrace.cleanup()


class holder(object):
  """
  Context manager for a manager
  """

  def __init__(self, creator, stacktrace = ignore):

    self.creator = creator
    self.manager = None
    self.stacktrace = stacktrace


  def __enter__(self):

    self.stacktrace.enter()
    self.manager = self.creator.create()
    return self.manager


  def __exit__(self, exc_type, exc_val, exc_tb):

    self.creator.destroy( manager = self.manager )
    self.stacktrace.exit( exc_type = exc_type, exc_val = exc_val, exc_tb = exc_tb )

    return False
