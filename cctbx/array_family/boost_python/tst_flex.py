import pickle
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_flex_triple(flex_triple, ordered, as_double=False):
  a = flex_triple()
  a = flex_triple(((1,2,3), (2,3,4), (3,4,5)))
  assert a.size() == 3
  assert tuple(a) == ((1,2,3), (2,3,4), (3,4,5))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)
  if (ordered):
    assert flex.order(a, b) == 0
  if (as_double):
    assert approx_equal(tuple(a.as_double()), (1,2,3,2,3,4,3,4,5))
    b = flex_triple().from_double(a.as_double())
    assert tuple(a) == tuple(b)

def exercise_flex_sym_mat3_double():
  a = flex.sym_mat3_double()
  a = flex.sym_mat3_double(((1,2,3,4,5,6), (2,3,4,5,6,7)))
  assert a.size() == 2
  assert tuple(a) == ((1,2,3,4,5,6), (2,3,4,5,6,7))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)

def exercise_flex_hendrickson_lattman():
  a = flex.hendrickson_lattman()
  a = flex.hendrickson_lattman(((1,2,3,4), (2,3,4,5), (3,4,5,6)))
  assert a.size() == 3
  assert tuple(a) == ((1,2,3,4), (2,3,4,5), (3,4,5,6))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)

def exercise_flex_tiny_size_t_2():
  a = flex.tiny_size_t_2()
  a = flex.tiny_size_t_2(((1,2), (2,3), (3,4)))
  assert a.size() == 3
  assert tuple(a) == ((1,2), (2,3), (3,4))

def exercise_flex_xray_scatterer():
  from cctbx import uctbx, sgtbx, xray
  uc = uctbx.unit_cell((10,11,12))
  sg = sgtbx.space_group_info("P 2")
  a = flex.xray_scatterer()
  a = flex.xray_scatterer((
    xray.scatterer("Si1", (0.1,0.2,0.3)),
    xray.scatterer("O1", (0.2,0.3,0.4), (1,2,3,-0.1,0.2,-0.3), 0.9),
    xray.scatterer("K1", (0.3,0.4,0.5), (3,1,2,-0.2,0.3,-0.1), 0.8,
                         fp_fdp=5+7j)))
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
  for i,ai in a.items():
    bi = b[i]
    assert ai.label == bi.label
    assert ai.caasf.label() == bi.caasf.label()
    assert approx_equal(ai.fp_fdp, bi.fp_fdp)
    assert approx_equal(ai.site, bi.site)
    assert ai.anisotropic_flag == bi.anisotropic_flag
    assert ai.u_iso == bi.u_iso
    assert ai.u_star == bi.u_star
    assert ai.multiplicity() == bi.multiplicity()
    assert approx_equal(ai.weight(), bi.weight())

def run():
  exercise_flex_triple(flex.miller_index, ordered=True)
  exercise_flex_triple(flex.vec3_double, ordered=False, as_double=True)
  exercise_flex_sym_mat3_double()
  exercise_flex_hendrickson_lattman()
  exercise_flex_tiny_size_t_2()
  exercise_flex_xray_scatterer()
  print "OK"

if (__name__ == "__main__"):
  run()
