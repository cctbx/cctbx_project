""" Constraints on scatterer sites """
from __future__ import division

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
    for k,v in kwds.iteritems():
      if k == 'rotating':
        if v and not self.staggered: self.rotating = True
      elif v != getattr(cls, k): setattr(self, k, v)

  def finalise(self, first, last):
    """ finalise the construction with the scatterer index range [first, last)
    """

  def get_parameter_set(self, reparametrisation):
    rv_l = []
    for s in self.constrained_site_indices: rv_l.append("%s_xyz" %s)
    rv = set(rv_l)
    if len(rv_l) != len(rv) or len(reparametrisation.constrained_parameters&rv) != 0:
      print("Redundant atoms in %s - '%s' skipping" %(
        self.__class__.__name__,
        reparametrisation.format_scatter_list(self.constrained_site_indices)))
      return None
    return rv

  def __repr__(self):
    return "%s(\n%s)" % (self.__class__.__name__,
                       '\n'.join(["  %s=%r" % (a,v)
                                 for (a,v) in self.__dict__.iteritems()]))
