from cctbx import adp_restraints, uctbx, adptbx
from scitbx import random, linalg, matrix
from libtbx.test_utils import approx_equal

def exercise(n_trials):
  site_coord = random.variate(random.uniform_distribution(0, 1))
  u_eigenval = random.variate(random.uniform_distribution(0.0005, 0.003))
  g_eigenval = random.variate(random.uniform_distribution(1,10))
  symm_mat = linalg.random_normal_matrix_generator(3,3)\
           .symmetric_matrix_with_eigenvalues
  def as_sym_mat3(packed_u):
    return matrix.sqr(packed_u.matrix_packed_u_as_symmetric()).as_sym_mat3()

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
  s = 1e-3
  for i in xrange(n_trials):
    z = matrix.col(site_coord(3))
    dz = matrix.col(site_coord(3))*s
    u = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))
    du = matrix.col(as_sym_mat3(symm_mat(u_eigenval(3))))*s
    g  = matrix.col(as_sym_mat3(symm_mat(g_eigenval(3))))
    dg = matrix.col(as_sym_mat3(symm_mat(g_eigenval(3))))*s
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
    r = (taylor_diff - finite_diff)/finite_diff
    assert abs(r) < 0.001, r

def run():
  exercise(n_trials=50)
  print 'OK'

if __name__ == '__main__':
  run()
