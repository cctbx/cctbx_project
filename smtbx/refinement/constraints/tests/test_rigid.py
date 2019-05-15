from __future__ import absolute_import, division, print_function

from cctbx import uctbx, xray, crystal
from cctbx.array_family import flex
import math
import pytest
from scitbx import matrix as mat
from scitbx import sparse
from scitbx.matrix import col
from smtbx.refinement import constraints
from smtbx.refinement.constraints import rigid
from six.moves import range

def test_rigid_site_proxy(n=5):
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
  for i in range(n):
    sc = xray.scatterer('C%i' %i,
                        site=tuple(flex.random_double(3)))
    sc.flags.set_grad_site(True)
    rigid_group_scatterers.append(sc)
  phi = reparam.add(constraints.independent_scalar_parameter,
                    value=0.1, variable=True)
  size = reparam.add(constraints.independent_scalar_parameter,
                    value=1, variable=True)
  rigid_group = reparam.add(constraints.rigid_pivoted_rotatable_group,
                            pivot, pivot_neighbour,
                            azimuth=phi,
                            size=size,
                            scatterers=rigid_group_scatterers)
  proxies = [ ]
  for i in range(n):
    proxies.append(reparam.add(constraints.rigid_site_proxy,
                               parent=rigid_group,
                               index=i))
  reparam.finalise()

  assert str(reparam) == """\
digraph dependencies {
8 -> 0;
8 -> 3;
8 -> 6;
8 -> 7;
23 -> 8;
26 -> 8;
29 -> 8;
32 -> 8;
35 -> 8;
0 [label="independent_site_parameter (C#) #0"];
3 [label="independent_site_parameter (C##) #3"];
6 [label="independent_scalar_parameter #6"];
7 [label="independent_scalar_parameter #7"];
8 [label="rigid_pivoted_rotatable_group (C0, C1, C2, C3, C4) #8"];
23 [label="rigid_site_proxy #23"];
26 [label="rigid_site_proxy #26"];
29 [label="rigid_site_proxy #29"];
32 [label="rigid_site_proxy #32"];
35 [label="rigid_site_proxy #35"]
}"""

  reparam.linearise()
  jt = reparam.jacobian_transpose

  q = 2*3 + 1 + 1 # pivot, its neighbour, azimuthal angle, size
  jt0 = sparse.matrix(q, q + 2*3*n) # + rigid_group + constrained site proxies
  assert jt.n_rows == jt0.n_rows
  assert jt.n_cols == jt0.n_cols
  for i, j in zip(range(q, q+3*n), range(q+3*n, jt0.n_cols)):
    assert jt.col(i) == jt.col(j)

def test_rigid_pivoted_rotatable():
  uc = uctbx.unit_cell((1, 1, 1))
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=uc,
      space_group_symbol='hall: P 2x 2y'),
    scatterers=flex.xray_scatterer(( #triangle
      xray.scatterer('C0', site=(0,0,0)),
      xray.scatterer('C1', site=(0,2,0)),
      xray.scatterer('C2', site=(1,1,0)),
      )))
  r = constraints.ext.reparametrisation(xs.unit_cell())
  sc = xs.scatterers()
  pivot = r.add(constraints.independent_site_parameter, sc[0])
  pivot_neighbour = r.add(constraints.independent_site_parameter, sc[1])
  azimuth = r.add(constraints.independent_scalar_parameter,
                  value=math.pi/2, variable=True)
  size = r.add(constraints.independent_scalar_parameter,
                  value=1, variable=False)
  rg = r.add(constraints.rigid_pivoted_rotatable_group,
                pivot=pivot,
                pivot_neighbour=pivot_neighbour,
                azimuth=azimuth,
                size=size,
                scatterers=(sc[1], sc[2]))
  site_proxy = r.add(constraints.rigid_site_proxy, rg, 1)
  r.finalise()
  r.linearise()
  r.store()
  #check that proxy and the final results are the same...
  assert uc.distance(col(site_proxy.value), col(sc[2].site)) == pytest.approx(0, abs=1e-15)
  #rotation happens around the center of gravity
  assert uc.distance(col((0,1,1)), col(sc[2].site)) == pytest.approx(0, abs=1e-15)

@pytest.fixture
def rr():
  class rigid_rotatable(object):
    def __init__(self):
      self.size_value = 9
      self.rx = math.pi
      self.ry = math.pi/2
      self.rz = math.pi/3
      self.sites = ((0,0,0), (1,0,0), (0,1,0), (0,0,1))
      self.uc = uctbx.unit_cell((1, 1, 1))
      self.xs = xray.structure(
        crystal_symmetry=crystal.symmetry(
          unit_cell=self.uc,
          space_group_symbol='hall: P 2x 2y'),
        scatterers=flex.xray_scatterer(( #triangle
          xray.scatterer('C0'),
          xray.scatterer('C1'),
          xray.scatterer('C2'),
          xray.scatterer('C3'),
          )))
      self.center = col((0,0,0))
      for s in self.sites:
        self.center = self.center + col(s)
      self.center = self.center / len(self.sites)
      sc = self.xs.scatterers()
      for i, s in enumerate(self.sites):
        sc[i].site = s
  return rigid_rotatable()

def test_rotatable_expansion(rr):
  r = constraints.ext.reparametrisation(rr.uc)
  sc = rr.xs.scatterers()
  pivot = r.add(constraints.independent_site_parameter, sc[0])
  size = r.add(constraints.independent_scalar_parameter,
                  value=rr.size_value, variable=True)
  r_x = r.add(constraints.independent_scalar_parameter,
                  value=0, variable=False)
  r_y = r.add(constraints.independent_scalar_parameter,
                  value=0, variable=False)
  r_z = r.add(constraints.independent_scalar_parameter,
                  value=0, variable=False)
  rg = r.add(constraints.rigid_rotatable_expandable_group,
                pivot=pivot,
                size = size,
                alpha = r_x,
                beta = r_y,
                gamma = r_z,
                scatterers=(sc[1], sc[2], sc[3]))
  r.finalise()
  r.linearise()
  r.store()
  shift = col(rr.sites[0]) - (col(rr.sites[0])-rr.center)*rr.size_value
  for i in range(1,4):
    calc_site = (col(rr.sites[i])-rr.center)*rr.size_value + shift
    assert rr.uc.distance(calc_site, col(sc[i].site)) == pytest.approx(0, abs=1e-14)

def test_rotatable_rotation(rr):
  r = constraints.ext.reparametrisation(rr.uc)
  sc = rr.xs.scatterers()
  pivot = r.add(constraints.independent_site_parameter, sc[0])
  size = r.add(constraints.independent_scalar_parameter,
                  value=1, variable=False)
  r_x = r.add(constraints.independent_scalar_parameter,
                  value=math.pi, variable=True)
  r_y = r.add(constraints.independent_scalar_parameter,
                  value=math.pi/2, variable=True)
  r_z = r.add(constraints.independent_scalar_parameter,
                  value=math.pi/3, variable=True)
  rg = r.add(constraints.rigid_rotatable_expandable_group,
                pivot=pivot,
                size = size,
                alpha = r_x,
                beta = r_y,
                gamma = r_z,
                scatterers=(sc[1], sc[2], sc[3]))
  r.finalise()
  r.linearise()
  r.store()
  rx_m = mat.sqr((1, 0, 0,
              0, math.cos(rr.rx), -math.sin(rr.rx),
              0, math.sin(rr.rx), math.cos(rr.rx)))
  ry_m = mat.sqr((math.cos(rr.ry), 0, math.sin(rr.ry),
                  0, 1, 0,
                  -math.sin(rr.ry), 0, math.cos(rr.ry)))
  rz_m = mat.sqr((math.cos(rr.rz), -math.sin(rr.rz), 0,
                  math.sin(rr.rz), math.cos(rr.rz), 0,
                  0, 0, 1))
  R = rx_m*ry_m*rz_m #comulative rotation matrix
  shift = col(rr.sites[0])-col(mat.row(col(rr.sites[0])-rr.center)*R)
  for i in range(1,4):
    calc_site = col(mat.row(col(rr.sites[i])-rr.center)*R) + shift
    assert rr.uc.distance(calc_site, col(sc[i].site)) == pytest.approx(0, abs=1e-14)

idealised_def_len_ref = {
      "Naphthalene" : ((0.695,-1.203775),
                       (-0.695,-1.203775),
                       (-1.39,-0),
                       (-0.695,1.203775),
                       (0.695,1.203775),
                       (1.39,0),
                       (2.78,0),
                       (3.475,1.203775),
                       (2.78,2.407551),
                       (1.39,2.407551)),
      "Cp*" :         ((0.373269,-1.148804),
                       (-0.977231,-0.71),
                       (-0.977231,0.71),
                       (0.373269,1.148804),
                       (1.207924,0),
                       (0.701754,-2.159777),
                       (-1.837216,-1.334816),
                       (-1.837216,1.334816),
                       (0.701754,2.159777),
                       (2.270924,0)),
      "Cp" :          ((0.373269,-1.148804),
                       (-0.977231,-0.71),
                       (-0.977231,0.71),
                       (0.373269,1.148804),
                       (1.207924,0)),
      "Ph" :          ((0.695,-1.203775),
                       (-0.695,-1.203775),
                       (-1.39,-0),
                       (-0.695,1.203775),
                       (0.695,1.203775),
                       (1.39,0))
}

@pytest.mark.parametrize('fragment', ("Cp", "Ph", "Cp*", "Naphthalene"), ids=("Cp", "Ph", "CpStar", "Naphthalene"))
def test_idealised_generation(fragment):
  generated = rigid.idealised_fragment().generate_fragment(fragment)
  reference = idealised_def_len_ref[fragment]
  for i, gen in enumerate(generated):
    assert gen.x == pytest.approx(reference[i][0], abs=1e-6)
    assert gen.y == pytest.approx(reference[i][1], abs=1e-6)

def test_idealised_fitting():
  tested = rigid.idealised_fragment()
  source_pts_indices = [1, 2, 5]
  ref_sites = ((-0.695,-1.203775,0), (-1.39,-0,0), (1.39,0,0))
  control_pts_indices = [0, 1, 4]
  fragment = tested.generate_fragment("Ph")
  crds = tested.fit(fragment, ref_sites, control_pts_indices)
  uc = uctbx.unit_cell((1, 1, 1))
  for i, p in enumerate(source_pts_indices):
    assert uc.distance(col((fragment[p].x,fragment[p].y,0)),
                       crds[control_pts_indices[i]]) == pytest.approx(0, abs=1e-6)
