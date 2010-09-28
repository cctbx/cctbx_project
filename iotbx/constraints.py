""" Constraint types supported by the iotbx """

class constraint(object):

  staggered = False
  pivot = None
  bond_length = None
  rotating = False
  stretching = False
  constrained_site_indices = ()

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

class geometrical_hydrogens(constraint):

  def finalise(self, first, last):
    self.constrained_site_indices = tuple(xrange(first, last))

class tertiary_ch_site(geometrical_hydrogens):
  n_constrained_sites = 1

class secondary_ch2_sites(geometrical_hydrogens):
  n_constrained_sites = 2

class staggered_terminal_tetrahedral_xh3_sites(geometrical_hydrogens):
  n_constrained_sites = 3
  staggered = True

class secondary_planar_xh_site(geometrical_hydrogens):
  n_constrained_sites = 1

class staggered_terminal_tetrahedral_xh_site(geometrical_hydrogens):
  n_constrained_sites = 1
  staggered = True

class terminal_planar_xh2_sites(geometrical_hydrogens):
  n_constrained_sites = 2

class terminal_tetrahedral_xh3_sites(geometrical_hydrogens):
  n_constrained_sites = 3

class terminal_tetrahedral_xh_site(geometrical_hydrogens):
  n_constrained_sites = 1

class polyhedral_bh_site(geometrical_hydrogens):
  n_constrained_sites = 1

class terminal_linear_ch_site(geometrical_hydrogens):
  n_constrained_sites = 1
