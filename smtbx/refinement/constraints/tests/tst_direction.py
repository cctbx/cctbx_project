from __future__ import absolute_import, division, print_function

def exercise_direction():
  from cctbx.array_family import flex
  from cctbx import uctbx, xray, crystal
  from smtbx.refinement import constraints
  from scitbx.matrix import col, row
  from libtbx.test_utils import approx_equal

  uc = uctbx.unit_cell((1, 2, 3))
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=uc,
      space_group_symbol='hall: P 2x 2y'),
    scatterers=flex.xray_scatterer((
      xray.scatterer('C0', site=(0,0,0)),
      xray.scatterer('C1', site=(0,2,0)),
      xray.scatterer('C2', site=(1,1,0)),
      xray.scatterer('C3', site=(3,1,0)),
      )))
  r = constraints.ext.reparametrisation(xs.unit_cell())
  sc = xs.scatterers()
  site_0 = r.add(constraints.independent_site_parameter, sc[0])
  site_1 = r.add(constraints.independent_site_parameter, sc[1])
  site_2 = r.add(constraints.independent_site_parameter, sc[2])
  site_3 = r.add(constraints.independent_site_parameter, sc[3])
  d = constraints.vector_direction((site_0, site_1, site_2)).direction(uc)
  sd = constraints.static_direction.calc_best_line(uc, (site_0, site_1, site_2))
  assert approx_equal(d, sd, eps=1e-15)

  d = constraints.vector_direction((site_0, site_1)).direction(uc)
  assert approx_equal(d,
    row(uc.orthogonalize(col(sc[1].site)-col(sc[0].site))).normalize(),
    eps=1e-15)

  n = constraints.static_direction.calc_best_plane_normal(
    uc, (site_0, site_1, site_2))

  n1 = constraints.static_direction.calc_best_plane_normal(
    uc, (site_0, site_1, site_2, site_3))

  v01 = uc.orthogonalize(col(sc[0].site)-col(sc[1].site))
  v21 = uc.orthogonalize(col(sc[2].site)-col(sc[1].site))
  nc = row(v01).cross(row(v21)).normalize()
  assert approx_equal(n, n1, eps=1e-15)
  assert approx_equal(n, nc, eps=1e-15)

if __name__ == '__main__':
  exercise_direction()
  print('OK')
