from __future__ import division

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


class regular_exception(object):
  """
  Unexpected result
  """

  def __init__(self, exception):
    self.exception = exception

  def __call__(self):
    if hasattr(self.exception, 'exc_trace') and self.exception.exc_trace():
      # If exception contains full trace information, re-use it:
      raise self.exception.exc_trace()[0], self.exception.exc_trace()[1], self.exception.exc_trace()[2]
    if hasattr(self.exception, 'trace'):
      # If exception contains extended information, re-use it:
      raise Exception("%s. Embedded exception follows:\n%s" % (str(self.exception), self.exception.trace))
    # Otherwise just re-raise
    raise self.exception

  def __str__(self):
    return "%s( message = '%s' )" % ( self.exception.__class__.__name__, self.exception )


class sorry_exception(object):
  """
  Special handling for Sorry
  """

  def __init__(self, exception):
    self.args = exception.args

  def __call__(self):
    raise self.exception_class()( *self.args )

  def __str__(self):
    return "Sorry( message = '%s' )" % self.exception_class()( *self.args )

  @classmethod
  def exception_class(self):
    from libtbx.utils import Sorry
    return Sorry


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

def error(exception):
  return TRANSLATOR.convert( exception = exception )

register_unpickleable_exception( handler = sorry_exception )
