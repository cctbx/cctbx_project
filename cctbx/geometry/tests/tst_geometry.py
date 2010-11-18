from cctbx.array_family import flex
from cctbx import covariance, crystal, geometry, sgtbx, xray
from libtbx.test_utils import approx_equal
from scitbx import matrix

def exercise_geometry():
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))
  uc = xs.unit_cell()
  flags = xs.scatterer_flags()
  for f in flags:
    f.set_grad_site(True)
  xs.set_scatterer_flags(flags)
  cov = flex.double((1e-8,1e-9,2e-9,3e-9,4e-9,5e-9,
                          2e-8,1e-9,2e-9,3e-9,4e-9,
                               3e-8,1e-9,2e-9,3e-9,
                                    2e-8,1e-9,2e-9,
                                         3e-8,1e-9,
                                              4e-8))
  param_map = xs.parameter_map()
  cov_cart = covariance.orthogonalize_covariance_matrix(cov, uc, param_map)
  O = matrix.sqr(uc.orthogonalization_matrix())
  F = matrix.sqr(uc.fractionalization_matrix())
  sites_cart = xs.sites_cart()
  sites_frac = xs.sites_frac()
  # distances
  rt_mx_ji = sgtbx.rt_mx('-y,x-y,z-1/3')
  sites = (sites_cart[0], uc.orthogonalize(rt_mx_ji*sites_frac[1]))
  d = geometry.distance(sites)
  assert approx_equal(d.distance_model, 1.6159860469110217)
  v = matrix.col(sites[1]) - matrix.col(sites[0])
  r_inv_cart = (O * matrix.sqr(rt_mx_ji.r().inverse().as_double()) * F)
  g = d.d_distance_d_sites()
  g = matrix.row(g[0] + tuple(r_inv_cart*matrix.col(g[1])))
  f = g * matrix.sqr(cov_cart.matrix_packed_u_as_symmetric()) * g.transpose()
  assert approx_equal(d.variance(cov_cart, uc, rt_mx_ji), f[0], eps=1e-15)
  # angles
  rt_mx_ji = sgtbx.rt_mx('x-y,x,z-2/3')
  rt_mx_ki = sgtbx.rt_mx('-y,x-y,z-1/3')
  r_inv_cart_ji = (O * matrix.sqr(rt_mx_ji.r().inverse().as_double()) * F)
  r_inv_cart_ki = (O * matrix.sqr(rt_mx_ki.r().inverse().as_double()) * F)
  cov_a = covariance.extract_covariance_matrix_for_sites(flex.size_t([1,0,1]), cov_cart, param_map)
  sites = (uc.orthogonalize(rt_mx_ji*sites_frac[1]),
           sites_cart[0],
           uc.orthogonalize(rt_mx_ki*sites_frac[1]))
  a = geometry.angle(sites)
  assert approx_equal(a.angle_model, 101.30738566828551)
  g = a.d_angle_d_sites()
  g = matrix.row(tuple(r_inv_cart_ji*matrix.col(g[0])) + g[1] +
                 tuple(r_inv_cart_ki*matrix.col(g[2])))
  f = g * matrix.sqr(cov_a.matrix_packed_u_as_symmetric()) * g.transpose()
  assert approx_equal(
    a.variance(cov_a, uc, (rt_mx_ji, sgtbx.rt_mx(), rt_mx_ki)), f[0], eps=1e-15)

def run():
  exercise_geometry()
  print "OK"

if __name__ == '__main__':
  run()
