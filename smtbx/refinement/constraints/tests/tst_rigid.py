from __future__ import division
from cctbx.array_family import flex
from scitbx import sparse
from cctbx import uctbx, xray, crystal
from smtbx.refinement import constraints
from math import pi
import math
from scitbx.matrix import col
from scitbx import matrix as mat
from libtbx.test_utils import approx_equal
from smtbx.refinement.constraints import rigid

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
  size = reparam.add(constraints.independent_scalar_parameter,
                    value=1, variable=True)
  rigid_group = reparam.add(constraints.rigid_pivoted_rotable_group,
                            pivot, pivot_neighbour,
                            azimuth=phi,
                            size=size,
                            scatterers=rigid_group_scatterers)
  proxies = [ ]
  for i in xrange(n):
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
8 [label="rigid_pivoted_rotable_group (C0, C1, C2, C3, C4) #8"];
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
  for i,j in zip(xrange(q, q+3*n), xrange(q+3*n, jt0.n_cols)):
    assert jt.col(i) == jt.col(j)

def exercise_rigid_pivoted_rotable():
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
                  value=pi/2, variable=True)
  size = r.add(constraints.independent_scalar_parameter,
                  value=1, variable=False)
  rg = r.add(constraints.rigid_pivoted_rotable_group,
                pivot=pivot,
                pivot_neighbour=pivot_neighbour,
                azimuth=azimuth,
                size=size,
                scatterers=(sc[1], sc[2]))
  site_proxy = r.add(constraints.rigid_site_proxy, rg, 1)
  r.finalise()
  r.linearise()
  r.store()
  #check that proxy an the final results are the same...
  assert approx_equal(
    uc.distance(col(site_proxy.value), col(sc[2].site)), 0, eps=1e-15)
  #rotation happens around the center of gravity
  assert approx_equal(
    uc.distance(col((0,1,1)), col(sc[2].site)), 0, eps=1e-15)

class rigid_rotable(object):
  def __init__(self):
    self.size_value = 9
    self.rx = pi
    self.ry = pi/2
    self.rz = pi/3
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
    self.reset_sites()
  def reset_sites(self):
    sc = self.xs.scatterers()
    for i, s in enumerate(self.sites):
      sc[i].site = s
  def exercise_expansion(self):
    self.reset_sites()
    r = constraints.ext.reparametrisation(self.uc)
    sc = self.xs.scatterers()
    pivot = r.add(constraints.independent_site_parameter, sc[0])
    size = r.add(constraints.independent_scalar_parameter,
                    value=self.size_value, variable=True)
    r_x = r.add(constraints.independent_scalar_parameter,
                    value=0, variable=False)
    r_y = r.add(constraints.independent_scalar_parameter,
                    value=0, variable=False)
    r_z = r.add(constraints.independent_scalar_parameter,
                    value=0, variable=False)
    rg = r.add(constraints.rigid_rotable_expandable_group,
                  pivot=pivot,
                  size = size,
                  alpha = r_x,
                  beta = r_y,
                  gamma = r_z,
                  scatterers=(sc[1], sc[2], sc[3]))
    r.finalise()
    r.linearise()
    r.store()
    shift = col(self.sites[0]) - (col(self.sites[0])-self.center)*self.size_value
    for i in xrange(1,4):
      calc_site = (col(self.sites[i])-self.center)*self.size_value + shift
      assert approx_equal(
        self.uc.distance(
          calc_site, col(sc[i].site)), 0, eps=1e-14)
  def exercise_rotation(self):
    self.reset_sites()
    r = constraints.ext.reparametrisation(self.uc)
    sc = self.xs.scatterers()
    pivot = r.add(constraints.independent_site_parameter, sc[0])
    size = r.add(constraints.independent_scalar_parameter,
                    value=1, variable=False)
    r_x = r.add(constraints.independent_scalar_parameter,
                    value=pi, variable=True)
    r_y = r.add(constraints.independent_scalar_parameter,
                    value=pi/2, variable=True)
    r_z = r.add(constraints.independent_scalar_parameter,
                    value=pi/3, variable=True)
    rg = r.add(constraints.rigid_rotable_expandable_group,
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
                0, math.cos(self.rx), -math.sin(self.rx),
                0, math.sin(self.rx), math.cos(self.rx)))
    ry_m = mat.sqr((math.cos(self.ry), 0, math.sin(self.ry),
                    0, 1, 0,
                    -math.sin(self.ry), 0, math.cos(self.ry)))
    rz_m = mat.sqr((math.cos(self.rz), -math.sin(self.rz), 0,
                    math.sin(self.rz), math.cos(self.rz), 0,
                    0, 0, 1))
    R = rx_m*ry_m*rz_m #comulative rotation matrix
    shift = col(self.sites[0])-col(mat.row(col(self.sites[0])-self.center)*R)
    for i in xrange(1,4):
      calc_site = col(mat.row(col(self.sites[i])-self.center)*R) + shift
      assert approx_equal(
        self.uc.distance(
          calc_site, col(sc[i].site)), 0, eps=1e-14)
  def excercise(self):
    self.exercise_expansion()
    self.exercise_rotation()

class idealised(object):
  def __init__(self):
    self.def_len_ref = {
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
    self.tested = rigid.idealised_fragment()

  def compare_results(self, _a , _b):
    for i, a in enumerate(_a):
      assert approx_equal(a.x, _b[i][0], 1e-6)
      assert approx_equal(a.y, _b[i][1], 1e-6)

  def excersise_generation(self):
    frags = ("Cp", "Ph", "Cp*", "Naphthalene")
    for i in frags:
      self.compare_results(
        self.tested.generate_fragment(i), self.def_len_ref[i])

  def excersise_fitting(self):
    source_pts_indices = [1, 2, 5]
    ref_sites = ((-0.695,-1.203775,0), (-1.39,-0,0), (1.39,0,0))
    control_pts_indices = [0, 1, 4]
    fragment = self.tested.generate_fragment("Ph")
    crds = self.tested.fit(fragment, ref_sites, control_pts_indices)
    uc = uctbx.unit_cell((1, 1, 1))
    for i, p in enumerate(source_pts_indices):
      approx_equal(
        uc.distance(
          col((fragment[p].x,fragment[p].y,0)),
          crds[control_pts_indices[i]]), 1e-14)

  def excersise(self):
    self.excersise_generation()
    self.excersise_fitting()

def run():
  exercise_rigid_site_proxy()
  exercise_rigid_pivoted_rotable()
  rigid_rotable().excercise()
  idealised().excersise()
  print 'OK'

if __name__ == '__main__':
  run()
