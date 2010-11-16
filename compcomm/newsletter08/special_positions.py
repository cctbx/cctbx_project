from cctbx import crystal
from cctbx.array_family import flex

crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,10,90,90,90),
    space_group_symbol="Pm3m")
crystal_symmetry.show_summary()

special_position_settings = crystal_symmetry.special_position_settings(
  min_distance_sym_equiv=0.5)

site_symmetry = special_position_settings.site_symmetry(site=(0.3, 0.31, 0.111))
print "special position operator:", site_symmetry.special_op()
print "exact location of special position:", site_symmetry.exact_site()

site_constraints = site_symmetry.site_constraints()

n_indep = site_constraints.n_independent_params()
print "n_indep:", n_indep

site_indep = site_constraints.independent_params(
  all_params=site_symmetry.exact_site())
assert len(site_indep) == n_indep

site_indep_shifted = list(site_indep) # copy site_indep
site_indep_shifted[0] += 0.1

site_shifted = site_constraints.all_params(
  independent_params=site_indep_shifted)
print "site_shifted:", site_shifted

site_constraints_2 = special_position_settings.site_symmetry(
  site=(0.5, 0.111, 0.222)).site_constraints()


def f(u,v):
  x, y, z = site_constraints.all_params((u,v))
  return -x + 2*y + 3*z
u,v = 0.111, 0.888
h = 1e-12
df_du = (f(u+h,v) - f(u-h,v))/(2*h)
df_dv = (f(u,v+h) - f(u,v-h))/(2*h)
print df_du, df_dv
independent_gradients = site_constraints.independent_gradients(
  all_gradients=flex.double((-1,2,3)))
print "independent_gradients:", independent_gradients

frac_adp_constraints = site_symmetry.adp_constraints()
cart_adp_constraints = site_symmetry.cartesian_adp_constraints(
  crystal_symmetry.unit_cell())

