from cctbx import adp_restraints, uctbx, adptbx, sgtbx
import scitbx.random
from scitbx import linalg, matrix
from libtbx.test_utils import approx_equal
import libtbx.utils

site_coord = scitbx.random.variate(
  scitbx.random.uniform_distribution(0, 1))
u_eigenval = scitbx.random.variate(
  scitbx.random.uniform_distribution(0.0005, 0.003))
g_eigenval = scitbx.random.variate(
  scitbx.random.uniform_distribution(1,10))
direction = scitbx.random.variate(
  scitbx.random.uniform_distribution(0, 1))
variance_eigenval = scitbx.random.variate(
  scitbx.random.uniform_distribution(0.1, 10))
symm_mat = linalg.random_normal_matrix_generator(3,3)\
         .symmetric_matrix_with_eigenvalues
def as_sym_mat3(packed_u):
  return matrix.sqr(packed_u.matrix_packed_u_as_symmetric()).as_sym_mat3()


def exercise_mean_square_displacement(options, n_trials):
  if options.fix_random_seeds:
    scitbx.random.set_random_seed(0)

  # check adptbx.mean_square_displacement_difference
  # against adp_restraints.rigid_bond_pair
  for i in xrange(n_trials):
    x1 = matrix.col(site_coord(3))
    x2 = matrix.col(site_coord(3))
    u1 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
    u2 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
    g  = as_sym_mat3(symm_mat(g_eigenval(3)))
    uc = uctbx.unit_cell(metrical_matrix=g)
    hirshfeld = adptbx.mean_square_displacement(uc, x1-x2)
    rigid = adp_restraints.rigid_bond_pair(x1, x2, u1, u2, uc)
    del uc
    assert hirshfeld.well_defined
    h1 = hirshfeld(u1).value
    h2 = hirshfeld(u2).value
    assert approx_equal(abs(h1 - h2), rigid.delta_z(), eps=1e-12)

  # check gradients with finite difference
  finite_difference_computation = scitbx.math.finite_difference_computation()
  best_delta = finite_difference_computation.best_delta
  for i in xrange(n_trials):
    z = site_coord(3)
    dz = best_delta(z, site_coord(3))
    z = matrix.col(z)
    dz = matrix.col(dz)
    u_eigen = u_eigenval(3)
    du_eigen = best_delta(u_eigen, direction(3))
    u = matrix.col(as_sym_mat3(symm_mat(u_eigen)))
    du = matrix.col(as_sym_mat3(symm_mat(du_eigen)))
    g_eigen = g_eigenval(3)
    dg_eigen = best_delta(g_eigen, direction(3))
    g  = matrix.col(as_sym_mat3(symm_mat(g_eigen)))
    dg = matrix.col(as_sym_mat3(symm_mat(dg_eigen)))
    uc = uctbx.unit_cell(metrical_matrix=g)
    h = adptbx.mean_square_displacement(uc, z)(u)
    uc_p = uctbx.unit_cell(metrical_matrix=g+dg)
    h_p = adptbx.mean_square_displacement(uc_p, z+dz)(u+du).value
    uc_m = uctbx.unit_cell(metrical_matrix=g-dg)
    h_m = adptbx.mean_square_displacement(uc_m, z-dz)(u-du).value
    finite_diff = (h_p - h_m)/2
    taylor_diff = (  matrix.col(h.grad_u).dot(du)
                   + matrix.col(h.grad_z).dot(dz)
                   + matrix.col(h.grad_g).dot(dg) )
    assert approx_equal(taylor_diff, finite_diff,
                        eps=5*finite_difference_computation.precision),\
           (taylor_diff, finite_diff)

def exercise_hirshfeld_relative_difference(options, n_trials):
  if options.fix_random_seeds:
    scitbx.random.set_random_seed(0)

  operators = [ sgtbx.rt_mx(),
                sgtbx.rt_mx('-x, y, -z'),
                sgtbx.rt_mx('y, z, x') ]

  # check against adptbx.mean_square_displacement_difference
  for op in operators:
    for i in xrange(n_trials):
      x1 = matrix.col(site_coord(3))
      x2 = matrix.col(site_coord(3))
      u1 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
      u2 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
      g  = as_sym_mat3(symm_mat(g_eigenval(3)))
      uc = uctbx.unit_cell(metrical_matrix=g)
      hirshfeld = adptbx.mean_square_displacement(uc, x1 - matrix.col(op(x2)))
      h1 = hirshfeld(u1).value
      h2 = hirshfeld(op(u2)).value
      r = adptbx.relative_hirshfeld_difference(uc, x1, u1, x2, u2, op).value
      assert approx_equal(r, 2*(h1 - h2)/(h1 + h2), eps=1e-12)

  # check gradients with finite difference
  finite_difference_computation = scitbx.math.finite_difference_computation()
  best_delta = finite_difference_computation.best_delta
  for op in operators:
    for i in xrange(n_trials):
      x1 = site_coord(3)
      dx1 = best_delta(x1, site_coord(3))
      x1 = matrix.col(x1)
      dx1 = matrix.col(dx1)
      x2 = site_coord(3)
      dx2 = best_delta(x2, site_coord(3))
      x2 = matrix.col(x2)
      dx2 = matrix.col(dx2)
      u1_eigen = u_eigenval(3)
      du1_eigen = best_delta(u1_eigen, direction(3))
      u1 = matrix.col(as_sym_mat3(symm_mat(u1_eigen)))
      du1 = matrix.col(as_sym_mat3(symm_mat(du1_eigen)))
      u2_eigen = u_eigenval(3)
      du2_eigen = best_delta(u2_eigen, direction(3))
      u2 = matrix.col(as_sym_mat3(symm_mat(u2_eigen)))
      du2 = matrix.col(as_sym_mat3(symm_mat(du2_eigen)))
      g_eigen = g_eigenval(3)
      dg_eigen = best_delta(g_eigen, direction(3))
      g  = matrix.col(as_sym_mat3(symm_mat(g_eigen)))
      dg = matrix.col(as_sym_mat3(symm_mat(dg_eigen)))
      uc = uctbx.unit_cell(metrical_matrix=g)
      r = adptbx.relative_hirshfeld_difference(uc, x1, u1, x2, u2, op)
      uc_p = uctbx.unit_cell(metrical_matrix=g+dg)
      r_p = adptbx.relative_hirshfeld_difference(uc_p,
                                                 x1+dx1, u1+du1,
                                                 x2+dx2, u2+du2, op).value
      duc_p_params = matrix.col(uc_p.parameters()) - matrix.col(uc.parameters())
      uc_m = uctbx.unit_cell(metrical_matrix=g-dg)
      r_m = adptbx.relative_hirshfeld_difference(uc_m,
                                                 x1-dx1, u1-du1,
                                                 x2-dx2, u2-du2, op).value
      duc_m_params = matrix.col(uc_m.parameters()) - matrix.col(uc.parameters())
      finite_diff = (r_p - r_m)/2
      taylor_diff = (  matrix.col(r.grad_x1).dot(dx1)
                     + matrix.col(r.grad_x2).dot(dx2)
                     + matrix.col(r.grad_u1).dot(du1)
                     + matrix.col(r.grad_u2).dot(du2)
                     + matrix.col(r.grad_unit_cell_params).dot(
                       (duc_p_params - duc_m_params)/2))
    assert approx_equal(taylor_diff, finite_diff,
                        eps=5*finite_difference_computation.precision),\
           (taylor_diff, finite_diff)

  # check esd computation
  for op in operators:
    for i in xrange(n_trials):
      x1 = matrix.col(site_coord(3))
      x2 = matrix.col(site_coord(3))
      u1 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
      u2 = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
      g  = matrix.col(as_sym_mat3(symm_mat(g_eigenval(3))))
      r = adptbx.relative_hirshfeld_difference(uc, x1, u1, x2, u2, op)
      n = 2*(3+6) + 5
      v = linalg.random_normal_matrix_generator(n,n)\
         .symmetric_matrix_with_eigenvalues(variance_eigenval(n))
      sigma = r.esd(crystallographic_variance_matrix_packed_u=v,
                    index_x1=1, index_u1=5,
                    index_x2=13, index_u2=17,
                    a_b_c_alpha_beta_gamma_sigmas=(0,0,0,0,0,0))
      sigma_sq = sigma**2
      g = matrix.col(  (0,)  + r.grad_x1
                     + (0,)  + r.grad_u1
                     + (0,0) + r.grad_x2
                     + (0,)  + r.grad_u2)
      vv = matrix.rec(v.matrix_packed_u_as_symmetric(), (n,n))
      sigma_sq_0 = g.dot(vv*g)
      assert approx_equal(sigma_sq, sigma_sq_0, eps=1e-12)


def run():
  import optparse
  cmd_line = optparse.OptionParser()
  cmd_line.add_option("--fix-random-seeds",
                      action='store_true')
  cmd_line.add_option("--n-trials",
                      type='int',
                      default=None)
  cmd_line.add_option("--Verbose",
                      action='store_true')
  opts, args = cmd_line.parse_args()
  libtbx.utils.show_times_at_exit()
  exercise_mean_square_displacement(opts,
                                    n_trials=opts.n_trials or 50)
  exercise_hirshfeld_relative_difference(opts,
                                         n_trials=opts.n_trials or 50)

if __name__ == '__main__':
  run()
