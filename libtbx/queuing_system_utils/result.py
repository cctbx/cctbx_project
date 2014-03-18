from __future__ import division

class Success(object):
  """
  Valid result
  """

  def __init__(self, value):

    self.value = value


  def __call__(self):

    return self.value


class Error(object):
  """
  Unexpected result
  """

  def __init__(self, exception):

    self.exception = exception


  def __call__(self):

    raise self.exception


class Sorry(object):
  """
  Unpickleable exception
  """

  def __init__(self, exception):

    self.args = exception.args


  def __call__(self):

    from libtbx.utils import Sorry
    raise Sorry( *self.args )


def AnyException(exception):

  import libtbx.utils

  if isinstance( exception, libtbx.utils.Sorry ):
    return Sorry( exception = exception )

  else:
    return Error( exception = exception )
