from libtbx.test_utils import approx_equal
from cctbx import uctbx
from cctbx.array_family import flex
import pickle

def exercise_flex_miller_index():
  from scitbx.array_family.flex import exercise_triple
  exercise_triple(flex_triple=flex.miller_index, flex_order=flex.order)
  a = flex.miller_index([(1,2,3), (2,3,4)])
  assert approx_equal(a.as_vec3_double(), [(1,2,3), (2,3,4)])

def exercise_flex_sym_mat3_double():
  a = flex.sym_mat3_double()
  a = flex.sym_mat3_double(((1,2,3,4,5,6), (2,3,4,5,6,7)))
  assert a.size() == 2
  assert tuple(a) == ((1,2,3,4,5,6), (2,3,4,5,6,7))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)
  assert approx_equal(tuple(a.as_double()), (1,2,3,4,5,6,2,3,4,5,6,7))
  b = flex.sym_mat3_double(a.as_double())
  assert tuple(a) == tuple(b)

def exercise_flex_hendrickson_lattman():
  a = flex.hendrickson_lattman()
  a = flex.hendrickson_lattman(((1,2,3,4), (2,3,4,5), (3,4,5,6)))
  assert a.size() == 3
  assert tuple(a) == ((1,2,3,4), (2,3,4,5), (3,4,5,6))
  assert tuple(a+a) == ((2,4,6,8), (4,6,8,10), (6,8,10,12))
  a += a
  assert tuple(a) == ((2,4,6,8), (4,6,8,10), (6,8,10,12))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)
  centric_flags = flex.bool([False, True])
  phase_integrals = flex.complex_double([complex(0.5,-0.7), complex(-0.3,0.4)])
  a = flex.hendrickson_lattman(
    centric_flags=centric_flags,
    phase_integrals=phase_integrals,
    max_figure_of_merit=1-1.e-6)
  assert approx_equal(a, [(2.2684820912654264, -3.1758749277715967, 0, 0),
                          (-0.3295836866004328, 0.43944491546724396, 0, 0)])

def exercise_flex_tiny_size_t_2():
  a = flex.tiny_size_t_2()
  a = flex.tiny_size_t_2(((1,2), (2,3), (3,4)))
  assert a.size() == 3
  assert tuple(a) == ((1,2), (2,3), (3,4))
  assert tuple(a.column(0)) == (1,2,3)
  assert tuple(a.column(1)) == (2,3,4)

def exercise_flex_xray_scatterer():
  from cctbx import uctbx, sgtbx, xray
  uc = uctbx.unit_cell((10,11,12))
  sg = sgtbx.space_group_info("P 2")
  a = flex.xray_scatterer()
  a = flex.xray_scatterer((
    xray.scatterer("Si1", (0.1,0.2,0.3)),
    xray.scatterer("O1", (0.2,0.3,0.4), (1,2,3,-0.1,0.2,-0.3), 0.9),
    xray.scatterer("K1", (0.3,0.4,0.5), (3,1,2,-0.2,0.3,-0.1), 0.8,
                         fp=5, fdp=7)))
  assert a.size() == 3
  assert a[1].multiplicity() == 0
  a[1].apply_symmetry(uc, sg.group())
  assert a[1].multiplicity() == 2
  assert approx_equal(a[1].weight(), 0.9)
  a.front().occupancy = 0.8
  assert approx_equal(a[0].occupancy, 0.8)
  a.back().occupancy = 0.7
  assert approx_equal(a[-1].occupancy, 0.7)
  p = pickle.dumps(a)
  b = pickle.loads(p)
  for i,ai in enumerate(a):
    bi = b[i]
    assert ai.label == bi.label
    assert ai.scattering_type == bi.scattering_type
    assert approx_equal(ai.fp, bi.fp)
    assert approx_equal(ai.fdp, bi.fdp)
    assert approx_equal(ai.site, bi.site)
    assert ai.anisotropic_flag == bi.anisotropic_flag
    assert ai.u_iso == bi.u_iso
    assert ai.u_star == bi.u_star
    assert ai.multiplicity() == bi.multiplicity()
    assert approx_equal(ai.weight(), bi.weight())
  assert approx_equal(a.extract_sites(),
                      ((0.1,0.2,0.3),(0.2,0.3,0.4),(0.3,0.4,0.5)))
  a.set_sites(sites=flex.vec3_double(
    ((-0.1,-0.2,-0.3),(-0.2,-0.3,-0.4),(-0.3,-0.4,-0.5))))
  assert approx_equal(a.extract_sites(),
    ((-0.1,-0.2,-0.3),(-0.2,-0.3,-0.4),(-0.3,-0.4,-0.5)))
  assert approx_equal(a[1].site, (-0.2,-0.3,-0.4))
  assert approx_equal(a.extract_occupancies(), (0.8,0.9,0.7))
  a.set_occupancies(occupancies=flex.double((0.1,0.2,0.3)))
  assert approx_equal(a.extract_occupancies(), (0.1,0.2,0.3))
  assert approx_equal(a[1].occupancy, 0.2)
  assert approx_equal(a.extract_u_iso(), (0.0, -1.0, -1.0))
  a.set_u_iso(u_iso=flex.double((3,4,5)))
  assert approx_equal(a.extract_u_iso(), (3,4,5))
  assert approx_equal(a[1].u_iso, 4)
  assert approx_equal(a.extract_u_star(),
     [(-1,-1,-1,-1,-1,-1),
      (1,2,3,-0.1,0.2,-0.3),
      (3,1,2,-0.2,0.3,-0.1)])
  a.set_u_star(u_star=flex.sym_mat3_double(
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)]))
  assert approx_equal(a.extract_u_star(),
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)])
  assert approx_equal(a[1].u_star, (1,2,3,-0.6,0.2,-0.3))
  unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  a.set_u_cart(
    unit_cell=unit_cell,
    u_cart=flex.sym_mat3_double(
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)]))
  assert approx_equal(a.extract_u_cart(unit_cell=unit_cell),
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)])
  unit_cell = uctbx.unit_cell((10,10,10,90,90,90))
  a.set_u_cart(
    unit_cell=unit_cell,
    u_cart=flex.sym_mat3_double(
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)]))
  assert approx_equal(a.extract_u_star(),
    [(-0.01, -0.02, -0.01, -0.01, -0.01, -0.01),
     (0.01, 0.02, 0.03, -0.006, 0.002, -0.003),
     (0.03, 0.01, 0.02, -0.002, 0.005, -0.001)])
  assert approx_equal(a.extract_u_iso(), [3,4,5])
  a.set_u_iso_from_u_star(unit_cell=unit_cell)
  assert approx_equal(a.extract_u_iso(), [3,2,2])
  assert a.count_anisotropic() == 2
  assert a.count_anomalous() == 1

def run():
  exercise_flex_miller_index()
  exercise_flex_sym_mat3_double()
  exercise_flex_hendrickson_lattman()
  exercise_flex_tiny_size_t_2()
  exercise_flex_xray_scatterer()
  print "OK"

if (__name__ == "__main__"):
  run()
