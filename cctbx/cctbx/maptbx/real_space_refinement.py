import scitbx.lbfgs
from cctbx.array_family import flex
from cctbx import maptbx

# Done by Erik McKee
class lbfgs:
  def __init__(self,data,grid_mat,sites_cart,delta_h=1.0):
    self.sites_cart = sites_cart
    self.grid_mat = grid_mat
    self.data = data
    self.n = sites_cart.size()*3
    self.x = flex.double(self.n, 0)
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

  def __call__(self):
    self.apply_shifts()
    self.compute_target(compute_gradients=True)
    return self.x, self.residual, self.gradients.as_double()

  def apply_shifts(self):
    self.sites_shifted = self.sites_cart + flex.vec3_double(self.x)

  def compute_target(self, compute_gradients):
    self.residual = maptbx.real_space_refinement_residual(map=self.data,
      gridding_matrix=self.grid_mat,
      sites_cart=self.sites_shifted, weights=self.weights)
    if compute_gradients:
      self.gradients = maptbx.real_space_refinement_gradients(
        map=self.data,
        gridding_matrix=self.grid_mat,
        sites_cart=self.sites_shifted,
        delta_h=self.delta_h)
    else:
      self.gradients = None
