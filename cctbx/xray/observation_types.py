from __future__ import division
from builtins import object
class any(object):

  def __repr__(self):
    return "xray." + self.__class__.__name__

class amplitude(any): pass
class intensity(any): pass
class reconstructed_amplitude(amplitude): pass
