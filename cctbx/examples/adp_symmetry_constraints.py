from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import adptbx
import math

def run():
  #
  # try these Hall symbols:
  #   "P 1", "P 2", "P 3", "P 3*", "P 4", "P 6", "P 2 2 3"
  #
  space_group = sgtbx.space_group("P 3*") # Hall symbol

  #
  # initialization of space group symmetry constraints
  #
  adp_constraints = sgtbx.tensor_rank_2_constraints(
    space_group=space_group,
    reciprocal_space=True)

  #
  # number of independent u_star parameters
  #
  n_indep = adp_constraints.n_independent_params()

  #
  # arbitrary Miller index and u_star tensor
  #
  h = (3,1,2)
  u_star=(0.000004, 0.000004, 0.000007, 0.000002, 0.0000000, 0.0000000)
  # optional: enforce symmetry at the beginning
  u_star = space_group.average_u_star(u_star)

  #
  # pass u_indep to the minimizer
  #
  u_indep = adp_constraints.independent_params(all_params=u_star)
  assert len(u_indep) == n_indep

  #
  # "expand" the independent parameters modified by the minimizer
  #
  u_star = adp_constraints.all_params(independent_params=u_indep)
  assert len(u_star) == 6

  #
  # these calculations are completely independent of the symmetry
  #
  dwf = adptbx.debye_waller_factor_u_star(h, u_star)
  gc = adptbx.debye_waller_factor_u_star_gradient_coefficients(h)
  # all_gradients is an array of six values
  all_gradients = [-2*math.pi**2 * dwf * c for c in gc]
  assert len(all_gradients) == 6
  cc = adptbx.debye_waller_factor_u_star_curvature_coefficients(h)
  # all_curvatures is an array of 21 values (upper triangle of 6x6 matrix)
  all_curvatures = (-2*math.pi**2)**2 * dwf * cc
  assert len(all_curvatures) == 6*(6+1)//2

  #
  # here we apply the symmetry constraints to the gradients and curvatures
  #
  # g_indep is an array of n_indep values
  g_indep = adp_constraints.independent_gradients(
    all_gradients=all_gradients)
  assert len(g_indep) == n_indep
  # c_indep is an array of n_indep*(n_indep+1)/2 values (upper triangle)
  c_indep = adp_constraints.independent_curvatures(
    all_curvatures=all_curvatures)
  assert len(c_indep) == n_indep*(n_indep+1)//2
  # feed g_indep and c_indep to the minimizer

  #
  # initialization of site symmetry constraints
  # (for sites on special positions)
  #
  unit_cell = uctbx.unit_cell((12,12,15,90,90,120))
  space_group = sgtbx.space_group_info("P 6").group()
  site_symmetry = sgtbx.site_symmetry(
    unit_cell=unit_cell,
    space_group=space_group,
    original_site=(1/3.,2/3.,0), # site on 3-fold axis
    min_distance_sym_equiv=0.5)
  assert len(site_symmetry.matrices()) == 3
  adp_constraints = site_symmetry.adp_constraints()
  # use adp_constraints as before

  print("OK")

if (__name__ == "__main__"):
  run()
