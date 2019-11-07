from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx.array_family import flex

def run():
  #
  # Initialization of a crystal symmetry.
  #
  crystal_symmetry = crystal.symmetry(
    unit_cell=(12,12,15,90,90,120),
    space_group_symbol="P 6")
  crystal_symmetry.show_summary()

  #
  # Definition of tolerance for determination of special positions.
  # min_distance_sym_equiv is the minimum distance between
  # symmetry-equivalent sites. If the distance is smaller than
  # min_distance_sym_equiv, the site is moved to the exact
  # location of the closest special position.
  #
  special_position_settings = crystal_symmetry.special_position_settings(
    min_distance_sym_equiv=0.5)

  #
  # Computation of the site symmetry for a site located
  # (almost) on 3-fold axis.
  #
  site_symmetry = special_position_settings.site_symmetry(site=(0.33,0.67,0))
  print("special position operator:", site_symmetry.special_op())
  print("exact location of special position:", site_symmetry.exact_site())

  #
  # The site_symmetry object caches the site constraints, which
  # are generated on demand when the site_constraints() method
  # is called the first time.
  #
  site_constraints = site_symmetry.site_constraints()

  #
  # Number of independent coordinates.
  #
  n_indep = site_constraints.n_independent_params()
  print("n_indep:", n_indep)

  #
  # For refinement we need only the independent parameters.
  #
  site_indep = site_constraints.independent_params(
    all_params=site_symmetry.exact_site())
  assert len(site_indep) == n_indep

  #
  # During refinement the site is shifted.
  #
  site_indep_shifted = list(site_indep) # copy site_indep
  site_indep_shifted[0] += 0.1

  #
  # "Expand" the independent parameters modified by the refinement.
  # site_shifted is certain to move only along the 3-fold axis.
  #
  site_shifted = site_constraints.all_params(
    independent_params=site_indep_shifted)
  print("site_shifted:", site_shifted)

  #
  # During refinement gradients are calculated.
  # These calculations can be performed without considering the
  # site symmetry, which simplifies the algorithms. The full set
  # of gradients can easily be transformed to the smaller set of
  # gradients w.r.t. the independent parameters.
  #
  independent_gradients = site_constraints.independent_gradients(
    all_gradients=flex.double([-0.01, 0.03, 0.02]))
  print("independent_gradients:", independent_gradients)

  #
  # Refinement with second derivatives (curvatures) is supported
  # in a similar way.
  # all_curvatures is the upper triangle of the symmetric 3x3 matrix
  # of second derivatives, i.e. an array with 3*(3+1)/2 elements.
  # independent_curvatures is the upper triangle of the constraint
  # curvature matrix with n_indep*(n_indep+1)/2 values.
  #
  independent_curvatures = site_constraints.independent_curvatures(
    all_curvatures=flex.double([-1, 2, -3, 4, -5, 6]))
  print("independent_curvatures:", independent_curvatures)

  #
  # See also the comprehensive unit test exercising the
  # site_constraints class:
  #   cctbx/regression/tst_sgtbx_site_constraints.py
  #

  print("OK")

if (__name__ == "__main__"):
  run()
