from __future__ import absolute_import, division, print_function

class success(object):
  """
  Valid result
  """

  def __init__(self, value):

    self.value = value


  def __call__(self):

    return self.value


  def __str__(self):

    return "success( value = %s )" % self.value

# Exception handling classes
class regular_exception(object):
  """
  Unexpected result
  """

  def __init__(self, exception):

    self.exception = exception


  def __call__(self):

    return self.exception


  def __str__(self):

    return "%s( message = '%s' )" % ( self.exception.__class__.__name__, self.exception )


class sorry_exception(object):
  """
  Special handling for Sorry
  """

  def __init__(self, exception):

    self.args = exception.args


  def __call__(self):

    return self.exception_class()( *self.args )


  def __str__(self):

    return "Sorry( message = '%s' )" % self.exception_class()( *self.args )


  @classmethod
  def exception_class(self):

    from libtbx.utils import Sorry
    return Sorry

# Extendable exception handling registry
class exception_to_result_registry(object):
  """
  Translates unpickleable exceptions
  """

  def __init__(self):

    self.handlers = []


  def add_handler(self, handler_class):

    self.handlers.append( handler_class )


  def convert(self, exception):

    for handler_class in self.handlers:
      if isinstance( exception, handler_class.exception_class() ):
        return handler_class( exception = exception )

    else:
      return regular_exception( exception = exception )


TRANSLATOR = exception_to_result_registry()


def register_unpickleable_exception(handler):

  TRANSLATOR.add_handler( handler_class = handler )

# Register Sorry
register_unpickleable_exception( handler = sorry_exception )


class failure(object):
  """
  Exception and stacktrace propagation
  """

  def __init__(self, exception, traceback):

    self.exception = exception
    self.traceback = traceback


  def __call__(self):

    self.traceback( exception = self.exception() )


  def __str__(self):

    return str( self.exception )


def get_exception(process, exit_code):

  return getattr( process, "err", RuntimeError( "exit code = %s" % exit_code ) )


def get_crash_info(process):

  from libtbx.scheduling import stacktrace
  printout = getattr( process, "stacktrace", None )

  if printout is None:
    return stacktrace.no_crash_info

  else:
    return stacktrace.stacktrace_info.from_stderr( message = printout )


def get_traceback_info():

  from libtbx.scheduling import stacktrace
  import sys
  return stacktrace.traceback_info( traceback = sys.exc_info()[2] )


def error(exception, traceback = None):

  if traceback is None:
    from libtbx.scheduling import stacktrace
    traceback = stacktrace.no_crash_info

  return failure(
    exception = TRANSLATOR.convert( exception = exception ),
    traceback = traceback,
    )
