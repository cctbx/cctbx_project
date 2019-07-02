""" Constraints on scatterer sites """
from __future__ import absolute_import, division, print_function
import six

class any(object):
  """ Base class for any constraints of scatterer sites """

  # default value of the attributes of an instance of this class
  staggered = False
  pivot = None
  bond_length = None
  rotating = False
  stretching = False
  flapping = False
  angle = None
  constrained_site_indices = ()
  need_pivot_neighbour_substituent = False

  def __init__(self, **kwds):
    cls = self.__class__
    for k,v in six.iteritems(kwds):
      if k == 'rotating':
        if v and not self.staggered: self.rotating = True
      elif v != getattr(cls, k): setattr(self, k, v)

  def finalise(self, first, last):
    """ finalise the construction with the scatterer index range [first, last)
    """

  @property
  def constrained_parameters(self):
    return tuple((idx, 'site') for idx in self.constrained_site_indices)

  def __repr__(self):
    return "%s(\n%s)" % (self.__class__.__name__,
                       '\n'.join(["  %s=%r" % (a,v)
                                 for (a,v) in six.iteritems(self.__dict__)]))

  def __eq__(self, other):
    """ For debugging purposes mostly as it is not needed by the
        constraint framework.
    """
    if type(self) != type(other): return False
    for attr, val in six.iteritems(self.__dict__):
      if getattr(other, attr) != val: return False
    return True
