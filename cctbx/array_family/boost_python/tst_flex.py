from libtbx.test_utils import approx_equal
from cctbx import uctbx
from cctbx.array_family import flex
import pickle

def exercise_flex_miller_index():
  from scitbx.array_family.flex import exercise_triple
  exercise_triple(flex_triple=flex.miller_index, flex_order=flex.order)
  a = flex.miller_index([(1,2,3), (2,3,4)])
  assert approx_equal(a.as_vec3_double(), [(1,2,3), (2,3,4)])
  h, k, l = [flex.int((0,1,2,3)),
             flex.int((1,2,3,4)),
             flex.int((2,3,4,5))]
  b = flex.miller_index(h, k, l)
  assert approx_equal(b, ((0,1,2),(1,2,3),(2,3,4),(3,4,5)))

def exercise_flex_hendrickson_lattman():
  a = flex.hendrickson_lattman()
  assert a.size() == 0
  a = flex.hendrickson_lattman(132)
  for x in a:
    assert x == (0,0,0,0)
  a = flex.hendrickson_lattman(((1,2,3,4), (2,3,4,5), (3,4,5,6)))
  assert a.size() == 3
  assert a.count((1,2,3,4)) == 1
  assert a.count((0,0,0,0)) == 0
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
  assert approx_equal(
    [a.slice(i) for i in xrange(4)],
    [[2.2684820912654264, -0.3295836866004328],
     [-3.1758749277715967, 0.43944491546724396],
     [0.0, 0.0],
     [0.0, 0.0]])
  a = flex.hendrickson_lattman(3, (1,2,3,4))
  assert a.all_eq((1,2,3,4))
  assert not a.all_eq((1,2,0,4))
  assert approx_equal(a.conj(), [(1,-2,3,-4), (1,-2,3,-4), (1,-2,3,-4)])
  assert approx_equal(a.conj().conj(), a)
  #
  a = flex.double()
  h = flex.hendrickson_lattman(a=a, b=a)
  assert h.size() == 0
  a = flex.double([1,2,3])
  b = flex.double([-3,4,5])
  h = flex.hendrickson_lattman(a=a, b=b)
  assert approx_equal(h, [(1,-3,0,0), (2,4,0,0), (3,5,0,0)])
  assert approx_equal(h == (1,-3,0,0), (True,False,False))
  assert approx_equal(h != (1,-3,0,0), (False,True,True))
  assert approx_equal(h != (0,0,0,0), (True,True,True))
  assert approx_equal(h == h.deep_copy(), (True, True, True))
  assert approx_equal(
    h == flex.hendrickson_lattman(a=b, b=a), (False, False, False))
  assert approx_equal(
    h != flex.hendrickson_lattman(a=b, b=a), (True, True, True))
  assert approx_equal(
    h != flex.hendrickson_lattman(a=b, b=a), (True, True, True))
  assert approx_equal(
    h != h.deep_copy(), (False, False, False))


def exercise_flex_xray_scatterer():
  from cctbx import uctbx, sgtbx, xray
  uc = uctbx.unit_cell((10,11,12))
  sg = sgtbx.space_group_info("P 2")
  a = flex.xray_scatterer()
  assert a.size() == 0
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
  a[0].flags.set_grad_site(state=True)
  a[1].flags.set_grad_fp(state=True)
  a[2].flags.param = -234
  p = pickle.dumps(a)
  b = pickle.loads(p)
  a_ = a.deep_copy()
  assert a_.n_grad_u_iso() == a_.n_grad_u_aniso() == 0
  a_[0].flags.set_grad_u_iso(state=True)
  a_[1].flags.set_grad_u_aniso(state=True)
  a_[2].flags.set_grad_u_aniso(state=True)
  assert a_.n_grad_u_iso() == 1
  assert a_.n_grad_u_aniso() == 2
  for i,ai in enumerate(a):
    bi = b[i]
    assert ai.label == bi.label
    assert ai.scattering_type == bi.scattering_type
    assert approx_equal(ai.fp, bi.fp)
    assert approx_equal(ai.fdp, bi.fdp)
    assert approx_equal(ai.site, bi.site)
    assert ai.flags.use_u_aniso() == bi.flags.use_u_aniso()
    assert ai.u_iso == bi.u_iso
    assert ai.u_star == bi.u_star
    assert ai.multiplicity() == bi.multiplicity()
    assert approx_equal(ai.weight(), bi.weight())
    assert ai.flags.bits == bi.flags.bits
    assert ai.flags.param == bi.flags.param
  assert b[0].flags.grad_site()
  assert not b[0].flags.grad_fp()
  assert not b[1].flags.grad_site()
  assert b[1].flags.grad_fp()
  assert b[2].flags.param == -234
  assert list(a.extract_labels()) == ["Si1", "O1", "K1"]
  assert list(a.extract_scattering_types()) == ["Si", "O", "K"]
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
  assert approx_equal(a.extract_u_iso_or_u_equiv(unit_cell=uc),
    (0.0, 258, 236+1/3.))
  a.set_u_iso(u_iso=flex.double((3,4,5)), selection=flex.bool(a.size(), True),
    unit_cell = uc)
  assert approx_equal(a.extract_u_iso(), (3,-1,-1))
  assert approx_equal(a.extract_u_iso_or_u_equiv(unit_cell=uc),
    (3, 4, 5))
  assert approx_equal(a[1].u_iso, -1)
  u_cart_answer = [(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
                   (4, 4, 4, 0, 0, 0), (5, 5, 5, 0, 0, 0)]
  assert approx_equal(a.extract_u_cart(uc), u_cart_answer)
  a.set_u_star(u_star=flex.sym_mat3_double(
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)]))
  assert approx_equal(a.extract_u_star(),
     [(-1,-1,-1,-1,-1,-1),
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
     [(-1,-1,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)])
  #
  a.set_u_cart(unit_cell = unit_cell,
               u_cart    = flex.sym_mat3_double([(1,2,3,4,5,6),
                                                 (0,0,0,1,2,3),
                                                 (1,2,3,0,0,0)]),
               selection = flex.size_t([1,2]))
  assert approx_equal(a.extract_u_cart(unit_cell=unit_cell),
     [(-1,-1,-1,-1,-1,-1),
      (0,0,0,1,2,3),
      (1,2,3,0,0,0)])
  #
  unit_cell = uctbx.unit_cell((10,10,10,90,90,90))
  a.set_u_cart(
    unit_cell=unit_cell,
    u_cart=flex.sym_mat3_double(
     [(-1,-2,-1,-1,-1,-1),
      (1,2,3,-0.6,0.2,-0.3),
      (3,1,2,-0.2,0.5,-0.1)]))
  assert approx_equal(a.extract_u_star(),
    [(-1,-1,-1,-1,-1,-1),
     (0.01, 0.02, 0.03, -0.006, 0.002, -0.003),
     (0.03, 0.01, 0.02, -0.002, 0.005, -0.001)])
  assert approx_equal(a.extract_u_iso(), [3,-1,-1])
  a.scale_adps(2.0)
  assert approx_equal(a.extract_u_star(),
    [(-1,-1,-1,-1,-1,-1),
     (0.02, 0.04, 0.06, -0.012, 0.004, -0.006),
     (0.06, 0.02, 0.04, -0.004, 0.01,  -0.002)])
  assert approx_equal(a.extract_u_iso(), [6,-1,-1])
  assert a.count_anisotropic() == 2
  assert a.count_anomalous() == 1
  a.convert_to_isotropic(unit_cell=unit_cell)
  assert a.count_anisotropic() == 0
  a.convert_to_anisotropic(unit_cell=unit_cell)
  assert a.count_anisotropic() == 3
  m = a.sites_mod_positive()
  assert approx_equal(m.extract_sites(), [
    (0.9,0.8,0.7),
    (0.8,0.7,0.6),
    (0.7,0.6,0.5)])
  m[2].site = (0.7,0.6,1.4) # to avoid +-0.5 ambiguity
  m = m.sites_mod_short()
  assert approx_equal(m.extract_sites(), [
    (-0.1,-0.2,-0.3),
    (-0.2,-0.3,-0.4),
    (-0.3,-0.4,0.4)])
  #
  assert a.extract_grad_u_iso().all_eq(False)
  a[1].flags.set_grad_u_iso(state=True)
  assert list(a.extract_grad_u_iso()) == [False, True, False]

def exercise_extract_u_cart_plus_u_iso():
  from cctbx import uctbx, sgtbx, xray
  uc = uctbx.unit_cell((1,1,1))
  sg = sgtbx.space_group_info("P 1")
  a = flex.xray_scatterer()
  assert a.size() == 0
  s1 = xray.scatterer(label = "C", u = 0.1)
  s2 = xray.scatterer(label = "C", u = 0.1)
  s2.flags.set_use_u_iso(False)
  s3 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s4 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s4.flags.set_use_u_aniso(False)
  s5 = xray.scatterer(label = "C", u = 0.1)
  s5.u_star=(1,1,1,1,1,1)
  s5.flags.set_use_u_aniso(True)
  s6 = xray.scatterer(label = "C", u = 0.1)
  s6.u_star=(1,1,1,1,1,1)
  s7 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s7.u_iso=0.1
  s8 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s8.u_iso=0.1
  s8.flags.set_use_u_iso(True)
  s9 = xray.scatterer(label = "C")
  s10 = xray.scatterer(label = "C")
  s10.flags.set_use_u_iso(False)
  a = flex.xray_scatterer((s1,s2,s3,s4,s5,s6,s7,s8,s9,s10))
  u_cart_total = a.extract_u_cart_plus_u_iso(uc)
  assert approx_equal(u_cart_total,
    [(0.1,0.1,0.1,0,0,0),
    (0,0,0,0,0,0),
    (1,1,1,1,1,1),
    (0,0,0,0,0,0),
    (1.1,1.1,1.1,1,1,1),
    (0.1,0.1,0.1,0,0,0),
    (1,1,1,1,1,1),
    (1.1,1.1,1.1,1,1,1),
    (0,0,0,0,0,0),
    (0,0,0,0,0,0)])

def run():
  exercise_flex_miller_index()
  exercise_flex_hendrickson_lattman()
  exercise_flex_xray_scatterer()
  exercise_extract_u_cart_plus_u_iso()
  print "OK"

if (__name__ == "__main__"):
  run()
