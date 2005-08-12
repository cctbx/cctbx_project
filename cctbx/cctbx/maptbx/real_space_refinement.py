import scitbx.lbfgs
from cctbx.array_family import flex
from cctbx import maptbx

# Done by Erik McKee
class lbfgs(object):
  '''
     If you run the lbfgs with a non-crystallographic interpolator, and it
     throws an out-of-bounds error, then change the call to the intepolator to:
        get_non_crystallographic_interpolator(map,grid_mtx,True,0.0)
  '''
  def __init__(self,interpolator,sites_cart,delta_h=1.0):
    self.sites_cart = sites_cart
    self.interpolator = interpolator
    self.x = flex.double(sites_cart.size()*3, 0)
    self.weights = flex.double(sites_cart.size(),1.0)
    self.delta_h = delta_h
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
                      termination_params=
                        scitbx.lbfgs.termination_parameters(max_iterations=1000))
    self.apply_shifts()
    self.compute_target(compute_gradients=False)
    sites_cart.clear()
    sites_cart.extend(self.sites_shifted)
    del self.sites_shifted
    del self.x

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    self.compute_target(compute_gradients=True)
    return self.residual, self.gradients.as_double()

  def apply_shifts(self):
    self.sites_shifted = self.sites_cart + flex.vec3_double(self.x)

  def compute_target(self, compute_gradients):
    self.residual = maptbx.real_space_refinement_residual(
      interpolator=self.interpolator,
      sites=self.sites_shifted,
      weights=self.weights)
    if compute_gradients:
      self.gradients = maptbx.real_space_refinement_gradients(
        interpolator=self.interpolator,
        sites=self.sites_shifted,
        delta_h=self.delta_h)
    else:
      self.gradients = None
