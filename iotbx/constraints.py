""" Constraint types supported by the iotbx """

class constraint(object):

  staggered = False
  pivot = None
  bond_length = None
  rotating = False
  constrained_site_indices = ()

  def __init__(self, **kwds):
    cls = self.__class__
    for k,v in kwds.iteritems():
      if k == 'rotating':
        if v and not self.staggered: self.rotating = True
      elif v != getattr(cls, k): setattr(self, k, v)

  def __eq__(self, other):
    return self.__dict__ == other.__dict__

  def __repr__(self):
    return "%s(\n%s)" % (self.__class__.__name__,
                       '\n'.join("  %s=%r" % (a,v)
                                 for (a,v) in self.__dict__.iteritems()))

class tertiary_ch_site(constraint): pass
class secondary_ch2_sites(constraint): pass
class staggered_terminal_tetrahedral_xh3_sites(constraint):
  staggered = True
class secondary_planar_xh_site(constraint): pass
class staggered_terminal_tetrahedral_xh_site(constraint):
  staggered = True
class terminal_planar_xh2_sites(constraint): pass
class terminal_tetrahedral_xh3_sites(constraint): pass
class terminal_tetrahedral_xh_site(constraint): pass
class polyhedral_bh_site(constraint): pass
class terminal_linear_ch_site(constraint): pass
