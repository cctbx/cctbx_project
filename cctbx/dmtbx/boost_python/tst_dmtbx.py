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
  t = dmtbx.triplet_generator(sg, i, 00000, 00000)
  assert t.t_den() == sg.t_den()
  assert t.sigma_2_only() == 00000
  assert t.discard_weights() == 00000
  t = dmtbx.triplet_generator(sg, i, 0001, 00000)
  assert t.sigma_2_only() == 0001
  assert t.discard_weights() == 00000
  t = dmtbx.triplet_generator(sg, i, 00000, 0001)
  assert t.sigma_2_only() == 00000
  assert t.discard_weights() == 0001
  assert tuple(t.n_relations()) == (0,0)
  assert t.relations_for(0) == ()
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (0,0))
  s = flex.bool()
  o = flex.size_t()
  r = t.raw_apply_tangent_formula(a, a, s, o, 00000, 1.e-15)
  assert approx_equal(tuple(r), (1,2))
  i = flex.miller_index(((4,6,0),(5,2,5),(6,1,5)))
  a = flex.double((1,2,3))
  t = dmtbx.triplet_generator(sg, i, 00000, 00000)
  assert tuple(t.n_relations()) == (4,2,2)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,2,5) 1 (6,1,5) 0 3 2",
          "(4,6,0) (5,2,5) 0 (6,1,5) 1 9 2"]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,2,5) (4,6,0) 0 (6,1,5) 0 3 2"]
  assert [r.format(i, 2) for r in t.relations_for(2)] \
      == ["(6,1,5) (4,6,0) 0 (5,2,5) 0 9 2"]
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (24,6,4))
  t = dmtbx.triplet_generator(sg, i, 00000, 0001)
  assert tuple(t.n_relations()) == (1,1,1)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,2,5) 1 (6,1,5) 0 3 1"]
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (6,3,2))
  t = dmtbx.triplet_generator(sg, i, 00000, 00000)
  r0 = t.relations_for(0)
  r1 = t.relations_for(1)
  assert r0[0].is_sigma_2(0)
  assert r0[0].is_similar_to(r0[1])
  assert not r0[0].is_similar_to(r1[0])
  i = flex.miller_index(((4,6,0),(5,1,2)))
  t = dmtbx.triplet_generator(sg, i, 00000, 00000)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,1,2) 0 (5,1,2) 1 6 4"]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,1,2) (4,6,0) 0 (5,1,2) 0 6 4"]
  assert not t.relations_for(0)[0].is_sigma_2(0)
  t = dmtbx.triplet_generator(sg, i, 00000, 0001)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,1,2) 0 (5,1,2) 1 6 1"]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,1,2) (4,6,0) 0 (5,1,2) 0 6 1"]
  t = dmtbx.triplet_generator(sg, i, 0001, 00000)
  assert tuple(t.n_relations()) == (0,0)

def run():
  exercise_triplet_generator()
  print "OK"

if (__name__ == "__main__"):
  run()
