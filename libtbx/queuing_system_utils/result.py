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

