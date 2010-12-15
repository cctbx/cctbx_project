""" Constraints on scatterer sites """

class any(object):
  """ Base class for any constraints of scatterer sites """

  # default value of the attributes of an instance of this class
  staggered = False
  pivot = None
  bond_length = None
  rotating = False
  stretching = False
  constrained_site_indices = ()
  need_pivot_neighbour_substituent = False

  def __init__(self, **kwds):
    cls = self.__class__
    for k,v in kwds.iteritems():
      if k == 'rotating':
        if v and not self.staggered: self.rotating = True
      elif v != getattr(cls, k): setattr(self, k, v)

  def finalise(self, first, last):
    """ finalise the construction with the scatterer index range [first, last)
    """

  def __eq__(self, other):
    return self.__dict__ == other.__dict__

  def __repr__(self):
    return "%s(\n%s)" % (self.__class__.__name__,
                       '\n'.join(["  %s=%r" % (a,v)
                                 for (a,v) in self.__dict__.iteritems()]))
