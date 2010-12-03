from cctbx.array_family import flex
from cctbx import covariance, crystal, geometry, sgtbx, uctbx, xray
from libtbx.test_utils import approx_equal, show_diff
from scitbx import matrix
import math
from cStringIO import StringIO

def quartz():
  return xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))

def quartz_p1(metrical_matrix=None):
  if metrical_matrix is not None:
    unit_cell = uctbx.unit_cell(metrical_matrix=metrical_matrix)
  else:
    unit_cell = uctbx.unit_cell(parameters=(5.01,5.01,5.47,90,90,120))
  return xray.structure(
    crystal_symmetry=crystal.symmetry(
    unit_cell=unit_cell,
    space_group_symbol='hall:  P 1'),
    scatterers=flex.xray_scatterer((
      xray.scatterer( #0
        label='Si',
        site=(0.500000, 0.500000, 0.333333),
        u=0.000000),
      xray.scatterer( #1
        label='Si',
        site=(0.000000, 0.500000, 0.666667),
        u=0.000000),
      xray.scatterer( #2
        label='Si',
        site=(0.500000, 0.000000, 0.000000),
        u=0.000000),
      xray.scatterer( #3
        label='O',
        site=(0.197000, 0.803000, 0.833333),
        u=0.000000),
      xray.scatterer( #4
        label='O',
        site=(0.394000, 0.197000, 0.166667),
        u=0.000000),
      xray.scatterer( #5
        label='O',
        site=(0.803000, 0.606000, 0.500000),
        u=0.000000),
      xray.scatterer( #6
        label='O',
        site=(0.197000, 0.394000, 0.500000),
        u=0.000000),
      xray.scatterer( #7
        label='O',
        site=(0.606000, 0.803000, 0.166667),
        u=0.000000),
      xray.scatterer( #8
        label='O',
        site=(0.803000, 0.197000, 0.833333),
        u=0.000000)
        )))

def exercise_geometry():
  xs = quartz()
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
  cell_vcv = flex.double((3e-2,3e-2,0,0,0,0,
                               3e-2,0,0,0,0,
                                 4e-2,0,0,0,
                                      0,0,0,
                                        0,0,
                                          0))
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
  assert approx_equal(
    0.0018054494791580823, d.variance(cov_cart, cell_vcv, uc, rt_mx_ji))
  rt_mx_ji = sgtbx.rt_mx('x+1,y,z')
  sites = (sites_cart[0], uc.orthogonalize(rt_mx_ji*sites_frac[0]))
  d = geometry.distance(sites)
  assert approx_equal(d.distance_model, uc.parameters()[0])
  assert approx_equal(cell_vcv.matrix_packed_u_diagonal()[0],
                      d.variance(cov_cart, cell_vcv, uc, rt_mx_ji))
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
  assert approx_equal(0.0042632511984529199,
    a.variance(cov_a, cell_vcv, uc, (rt_mx_ji, sgtbx.rt_mx(), rt_mx_ki)))

def exercise_grad_metrical_matrix():
  def calc_distance(metrical_matrix=None):
    xs = quartz_p1(metrical_matrix=metrical_matrix)
    uc = xs.unit_cell()
    sites_cart = xs.sites_cart()
    sites = (sites_cart[0], sites_cart[5])
    d = geometry.distance(sites)
    return d.distance_model
  def calc_angle(metrical_matrix=None):
    xs = quartz_p1(metrical_matrix=metrical_matrix)
    uc = xs.unit_cell()
    sites_cart = xs.sites_cart()
    sites = (sites_cart[4], sites_cart[0], sites_cart[5])
    a = geometry.angle(sites)
    return a.angle_model*(math.pi/180)
  xs = quartz_p1()
  uc = xs.unit_cell()
  sites_cart = xs.sites_cart()
  # distances
  sites = (sites_cart[0], sites_cart[5])
  d = geometry.distance(sites)
  grads = d.d_distance_d_metrical_matrix(uc)
  fd_grads = finite_differences(calc_distance, xs.unit_cell())
  assert approx_equal(grads, fd_grads)
  # angles
  sites = (sites_cart[4], sites_cart[0], sites_cart[5])
  a = geometry.angle(sites)
  grads = a.d_angle_d_metrical_matrix(uc)
  fd_grads = finite_differences(calc_angle, xs.unit_cell())
  assert approx_equal(grads, fd_grads)

def exercise_high_symmetry_case():
  """
  BF6 sitting on inversion centre in high symmetry space group.

  The angles formed between F20, B17 and either F18 or F19 are necessarily
  exactly equal to 90 degrees, and hence their variance should be exactly zero.
  The angles formed by F18, B17 and F19 are free to refine, and hence should
  have an associated variance.
  """
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
    unit_cell=(22.543, 22.543, 35.5167, 90, 90, 90),
    space_group_symbol='hall: -I 4bd 2'),
    scatterers=flex.xray_scatterer((
      xray.scatterer( #24
        label='F18',
        site=(0.500000, 0.970584, -0.041027),
        u=(0.000258, 0.000675, 0.000138, -0.000000, -0.000000, -0.000043),
        fp=0.017939,
        fdp=0.010330),
                     xray.scatterer( #25
        label='F19',
        site=(0.500000, 0.933678, 0.021766),
        u=(0.000178, 0.001046, 0.000158, -0.000000, -0.000000, 0.000144),
        fp=0.017939,
        fdp=0.010330),
                     xray.scatterer( #26
        label='F20',
        site=(0.438015, 1.000000, 0.000000),
        u=(0.000231, 0.000500, 0.000106, -0.000000, -0.000000, 0.000078),
        fp=0.017939,
        fdp=0.010330),
                     xray.scatterer( #27
        label='B17',
        site=(0.500000, 1.000000, 0.000000),
        u=(0.000133, 0.000624, 0.000051, -0.000000, -0.000000, 0.000048),
        fp=0.001413,
        fdp=0.000674)
        )))

  flags = xs.scatterer_flags()
  for f in flags: f.set_grad_site(True)
  xs.set_scatterer_flags(flags)
  sites_frac = xs.sites_frac()
  unit_cell = xs.unit_cell()
  cov = flex.double([
    0,0,0,0,0,0,0,0,0,0,0,0,5.5895145222321514e-07,-1.7223269587621895e-08,0,
    -2.0367609296594096e-07,-6.6231210988833185e-08,-2.6525243183929821e-09,0,0,
    0,0,0,1.1503438451957222e-07,0,-3.5665932804823073e-08,
    7.4165082206213702e-09,-6.8444182449209374e-09,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    9.0344019852871238e-07,1.2109549448751337e-07,-7.4794454124745421e-09,0,0,
    0,0,0,1.5840179985755122e-07,9.1193835484941833e-09,0,0,0,0,0,
    1.447679091961411e-07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  cell_cov = flex.double([
    2.4999999999999999e-07,2.4999999999999999e-07,0,0,0,0,
    2.4999999999999999e-07,0,0,0,0,1.4399999999999998e-06,0,0,0,0,0,0,0,0,0])
  param_map = xs.parameter_map()
  cov_cart = covariance.orthogonalize_covariance_matrix(
    cov, unit_cell, param_map)
  i_seqs = flex.size_t((2,3,1))
  cov_i_seqs = covariance.extract_covariance_matrix_for_sites(
    i_seqs, cov_cart, param_map)
  for sym_ops in (
      (sgtbx.rt_mx(),sgtbx.rt_mx(),sgtbx.rt_mx()),
      (sgtbx.rt_mx("1-x,2-y,-z"), sgtbx.rt_mx(), sgtbx.rt_mx()),
      (sgtbx.rt_mx(), sgtbx.rt_mx(), sgtbx.rt_mx("1-x,2-y,-z")),
      (sgtbx.rt_mx("1-x,2-y,-z"), sgtbx.rt_mx(), sgtbx.rt_mx("1-x,2-y,-z"))):
    site_frac_ji = sym_ops[0] * sites_frac[i_seqs[0]]
    site_frac_i = sym_ops[1] * sites_frac[i_seqs[1]]
    site_frac_ki = sym_ops[2] * sites_frac[i_seqs[2]]
    a = geometry.angle((unit_cell.orthogonalize(site_frac_ji),
                        unit_cell.orthogonalize(site_frac_i),
                        unit_cell.orthogonalize(site_frac_ki)))
    var = a.variance(cov_i_seqs, cell_cov, unit_cell, sym_ops)
    assert approx_equal(var, 0, eps=2e-16)
  pair_asu_table = xs.pair_asu_table(distance_cutoff=2)
  import libtbx.load_env
  if libtbx.env.has_module('iotbx'):
    import iotbx.cif
    s = StringIO()
    print >> s, iotbx.cif.distances_as_cif_loop(
      pair_asu_table,
      xs.scatterers().extract_labels(),
      sites_frac=xs.sites_frac(),
      covariance_matrix=cov,
      cell_covariance_matrix=cell_cov,
      parameter_map=param_map).loop
    assert not show_diff(s.getvalue(), """\
loop_
  _geom_bond_atom_site_label_1
  _geom_bond_atom_site_label_2
  _geom_bond_distance
  _geom_bond_site_symmetry_2
   F18 B17 1.601(13) 4_575
   F19 B17 1.683(18) .
   F20 B17 1.397(9) .

""")
    s = StringIO()
    print >> s, iotbx.cif.angles_as_cif_loop(
      pair_asu_table,
      xs.scatterers().extract_labels(),
      sites_frac=xs.sites_frac(),
      covariance_matrix=cov,
      cell_covariance_matrix=cell_cov,
      parameter_map=param_map).loop
    assert not show_diff(s.getvalue(), """\
loop_
  _geom_angle_atom_site_label_1
  _geom_angle_atom_site_label_2
  _geom_angle_atom_site_label_3
  _geom_angle
  _geom_angle_site_symmetry_1
  _geom_angle_site_symmetry_3
   F18 B17 F18 180.0 . 4_575
   F19 B17 F18 87.1(7) . 4_575
   F19 B17 F18 92.9(7) . .
   F19 B17 F18 92.9(7) 4_575 4_575
   F19 B17 F18 87.1(7) 4_575 .
   F19 B17 F19 180.0 4_575 .
   F20 B17 F18 90.0 . 4_575
   F20 B17 F18 90.0 . .
   F20 B17 F19 90.0 . .
   F20 B17 F19 90.0 . 4_575
   F20 B17 F18 90.0 9_675 4_575
   F20 B17 F18 90.0 9_675 .
   F20 B17 F19 90.0 9_675 .
   F20 B17 F19 90.0 9_675 4_575
   F20 B17 F20 180.0 9_675 .

""")


def finite_differences(func, unit_cell, eps=1.e-6):
  gradients = []
  for i in range(6):
    mm = list(unit_cell.metrical_matrix())
    mm[i] += 2*eps
    qm = func(metrical_matrix=mm)
    mm[i] -= 2*eps
    qp = func(metrical_matrix=mm)
    dq = (qm-qp)/(2*eps)
    gradients.append(dq)
  return gradients

def run():
  exercise_high_symmetry_case()
  exercise_grad_metrical_matrix()
  exercise_geometry()
  print "OK"

if __name__ == '__main__':
  run()
