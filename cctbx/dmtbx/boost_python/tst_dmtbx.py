from cctbx import dmtbx
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_triplet_generator():
  sg = sgtbx.space_group_info("P 41").group()
  i = flex.miller_index(((1,2,3),(2,3,4)))
  a = flex.double((1,2))
  t = dmtbx.triplet_generator(sg, i)
  t = dmtbx.triplet_generator(sg, i, 00000)
  assert t.sigma_2_only() == 00000
  assert t.t_den() == sg.t_den()
  assert tuple(t.n_relations(00000, 00000)) == (0,0)
  assert t.relations_for(0, 00000) == ()
  assert approx_equal(tuple(t.sum_of_amplitude_products(i, a, 00000, 00000)),
                      (0,0))
  s = flex.bool()
  o = flex.size_t()
  r = t.raw_apply_tangent_formula(a, a, s, o, 00000, 00000, 00000, 1.e-15)
  assert approx_equal(tuple(r), (1,2))
  i = flex.miller_index(((4,6,0),(5,2,5),(6,1,5)))
  a = flex.double((1,2,3))
  t = dmtbx.triplet_generator(sg, i, 0001)
  assert t.sigma_2_only() == 0001
  assert tuple(t.n_relations(00000, 00000)) == (4,2,2)
  assert tuple(t.n_relations(00000, 0001)) == (2,2,2)
  assert tuple(t.n_relations(0001, 00000)) == (2,1,1)
  assert tuple(t.n_relations(0001, 0001)) == (1,1,1)
  assert [r.format(i, 0) for r in t.relations_for(0, 00000)] \
      == ["(4,6,0) (5,2,5) 1 (6,1,5) 0 3 2",
          "(4,6,0) (5,2,5) 0 (6,1,5) 1 9 2"]
  assert [r.format(i, 0) for r in t.relations_for(0, 0001)] \
      == ["(4,6,0) (5,2,5) 1 (6,1,5) 0 3 2"]
  assert approx_equal(tuple(t.sum_of_amplitude_products(i, a, 00000, 00000)),
                      (24,6,4))
  assert approx_equal(tuple(t.sum_of_amplitude_products(i, a, 00000, 0001)),
                      (12,6,4))
  assert approx_equal(tuple(t.sum_of_amplitude_products(i, a, 0001, 00000)),
                      (12,3,2))
  assert approx_equal(tuple(t.sum_of_amplitude_products(i, a, 0001, 0001)),
                      (6,3,2))

def run():
  exercise_triplet_generator()
  print "OK"

if (__name__ == "__main__"):
  run()
