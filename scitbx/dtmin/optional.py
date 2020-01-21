from __future__ import print_function
from __future__ import division
from scitbx.array_family import flex

class Optional:
  def target_gradient(self):
    """Return a tuple of (target, gradient)"""
    frac_large = 0.01
    return self.finite_difference_gradient(frac_large)

  def target_gradient_hessian(self):
    """
    Return a tuple of (target, gradient, hessian, is_diagonal)
    target"""
    (target, gradient) = self.target_gradient()

    is_diagonal = True
    hessian = flex.double(self.nmp * self.nmp, 0)
    hessian.reshape(flex.grid(self.nmp, self.nmp))
    ls = self.macrocycle_large_shifts()
    for i in range(self.nmp):
      hessian[i,i] = 1/(ls[i]**2)

    return (target, gradient, hessian, is_diagonal)

  def bounds(self):
    """Return a list of Bounds objects, of length nmp"""
    return []

  def reparameterize(self):
    """Return a list of Reparams objects, of length nmp"""
    return []

  def reject_outliers(self):
    """
    A place to put any code for rejecting outliers in your data before
    doing any target calculations
    """
    pass

  def setup(self): pass

  def cleanup(self): pass

  def finalize(self): pass

  def maximum_distance_special(self, x, p, bounded, dist, start_distance):
    return 0.
