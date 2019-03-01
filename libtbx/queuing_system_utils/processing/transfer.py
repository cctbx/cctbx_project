from __future__ import absolute_import, division, print_function, with_statement

from six.moves import cPickle as pickle

class TemporaryFile(object):
  """
  Sends data by writing out a temporary file
  """

  SCRIPT= "( target, args, kwargs ) = pickle.load( open( \"%s\" ) )"

  def __init__(self, name, target, args, kwargs):

    self.target = "%s.target" % name

    with open( self.target, "wb" ) as ifile:
      pickle.dump( ( target, args, kwargs ), ifile )


  def script(self):

    return self.SCRIPT % self.target


  def files(self):

    return [ self.target ]


class Stdin(object):
  """
  Sends data by sending along the command file as a pickled string
  """

  SCRIPT = "( target, args, kwargs ) = pickle.loads( %r )"

  def __init__(self, name, target, args, kwargs):

    self.target = target
    self.args = args
    self.kwargs = kwargs


  def script(self):

    data = pickle.dumps( ( self.target, self.args, self.kwargs ) )
    return self.SCRIPT % data


  def files(self):

    return []
