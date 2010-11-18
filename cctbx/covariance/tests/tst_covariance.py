from cctbx.array_family import flex
from cctbx import covariance, crystal, xray
from libtbx.test_utils import approx_equal
from scitbx import matrix

def exercise_covariance():
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
  assert approx_equal(cov,
    covariance.extract_covariance_matrix_for_sites(flex.size_t([0,1]), cov, param_map))
  cov_cart = covariance.orthogonalize_covariance_matrix(cov, uc, param_map)
  O = matrix.sqr(uc.orthogonalization_matrix())
  for i in range(param_map.n_scatterers):
    cov_i = covariance.extract_covariance_matrix_for_sites(flex.size_t([i]), cov, param_map)
    cov_i_cart = covariance.extract_covariance_matrix_for_sites(flex.size_t([i]), cov_cart, param_map)
    assert approx_equal(
      O * matrix.sym(sym_mat3=cov_i) * O.transpose(),
      matrix.sym(sym_mat3=cov_i_cart).as_mat3())

def run():
  exercise_covariance()
  print "OK"

if __name__ == '__main__':
  run()
