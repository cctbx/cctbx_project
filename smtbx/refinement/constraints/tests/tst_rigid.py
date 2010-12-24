from __future__ import division
from scitbx.array_family import flex
from scitbx import sparse
from cctbx import uctbx, xray
from smtbx.refinement import constraints

def exercise_rigid_site_proxy(n=5):
  uc = uctbx.unit_cell((1, 2, 3))
  reparam = constraints.ext.reparametrisation(uc)
  independents = [ ]
  for name in ('C#', 'C##'):
    sc = xray.scatterer(name, site=tuple(flex.random_double(3)))
    sc.flags.set_grad_site(True)
    p = reparam.add(constraints.independent_site_parameter, sc)
    independents.append(p)
  pivot, pivot_neighbour = independents
  rigid_group_scatterers = [ ]
  for i in xrange(n):
    sc = xray.scatterer('C%i' %i,
                        site=tuple(flex.random_double(3)))
    sc.flags.set_grad_site(True)
    rigid_group_scatterers.append(sc)
  phi = reparam.add(constraints.independent_scalar_parameter,
                    value=0.1, variable=True)
  rigid_group = reparam.add(constraints.rigid_pivoted_rotable_group,
                            pivot, pivot_neighbour,
                            azimuth=phi,
                            scatterers=rigid_group_scatterers)
  proxies = [ ]
  for i in xrange(n):
    proxies.append(reparam.add(constraints.rigid_site_proxy,
                               parent=rigid_group,
                               index=i))
  reparam.finalise()

  assert str(reparam) == """\
digraph dependencies {
7 -> 0;
7 -> 3;
7 -> 6;
22 -> 7;
25 -> 7;
28 -> 7;
31 -> 7;
34 -> 7;
0 [label="independent_site_parameter (C#) #0"];
3 [label="independent_site_parameter (C##) #3"];
6 [label="independent_scalar_parameter #6"];
7 [label="rigid_pivoted_rotable_group (C0, C1, C2, C3, C4) #7"];
22 [label="rigid_site_proxy #22"];
25 [label="rigid_site_proxy #25"];
28 [label="rigid_site_proxy #28"];
31 [label="rigid_site_proxy #31"];
34 [label="rigid_site_proxy #34"]
}"""

  reparam.linearise()
  jt = reparam.jacobian_transpose

  q = 2*3 + 1 # pivot, its neighbour, azimuthal angle
  jt0 = sparse.matrix(q, q + 2*3*n) # + rigid_group + constrained site proxies
  assert jt.n_rows == jt0.n_rows
  assert jt.n_cols == jt0.n_cols
  for i,j in zip(xrange(q, q+3*n), xrange(q+3*n, jt0.n_cols)):
    assert jt.col(i) == jt.col(j)

def run():
  exercise_rigid_site_proxy()
  print 'OK'

if __name__ == '__main__':
  run()
