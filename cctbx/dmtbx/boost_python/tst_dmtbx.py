from __future__ import absolute_import, division, print_function
from cctbx import dmtbx
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise_triplet_generator():
  sg = sgtbx.space_group_info("P 41").group()
  i = flex.miller_index(((1,2,3),(2,3,4)))
  a = flex.double((1,2))
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, False)
  assert t.t_den() == sg.t_den()
  assert t.max_relations_per_reflection() == 0
  assert t.sigma_2_only() == False
  assert t.discard_weights() == False
  t = dmtbx.ext.triplet_generator(sg, i, a, 3, True, False)
  assert t.max_relations_per_reflection() == 3
  assert t.sigma_2_only() == True
  assert t.discard_weights() == False
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, True)
  assert t.sigma_2_only() == False
  assert t.discard_weights() == True
  assert tuple(t.n_relations()) == (0,0)
  assert t.relations_for(0) == ()
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (0,0))
  s = flex.bool()
  r = t.raw_apply_tangent_formula(a, a, s, False, False, 1.e-15)
  assert approx_equal(tuple(r), (1,2))
  i = flex.miller_index(((4,6,0),(5,2,5),(6,1,5)))
  a = flex.double((1,2,3))
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, False)
  assert tuple(t.n_relations()) == (4,2,2)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,2,5) %s (6,1,5) %s 3 2" % (True, False),
          "(4,6,0) (5,2,5) %s (6,1,5) %s 9 2" % (False, True)]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,2,5) (4,6,0) %s (6,1,5) %s 3 2" % (False, False)]
  assert [r.format(i, 2) for r in t.relations_for(2)] \
      == ["(6,1,5) (4,6,0) %s (5,2,5) %s 9 2" % (False, False)]
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (24,6,4))
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, True)
  assert tuple(t.n_relations()) == (1,1,1)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,2,5) %s (6,1,5) %s 3 1" % (True, False)]
  assert approx_equal(tuple(t.sums_of_amplitude_products(a)), (6,3,2))
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, False)
  r0 = t.relations_for(0)
  r1 = t.relations_for(1)
  assert r0[0].is_sigma_2(0)
  assert r0[0].is_similar_to(r0[1])
  assert not r0[0].is_similar_to(r1[0])
  i = flex.miller_index(((4,6,0),(5,1,2)))
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, False)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,1,2) %s (5,1,2) %s 6 4" % (False, True)]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,1,2) (4,6,0) %s (5,1,2) %s 6 4" % (False, False)]
  assert not t.relations_for(0)[0].is_sigma_2(0)
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, False, True)
  assert [r.format(i, 0) for r in t.relations_for(0)] \
      == ["(4,6,0) (5,1,2) %s (5,1,2) %s 6 1" % (False, True)]
  assert [r.format(i, 1) for r in t.relations_for(1)] \
      == ["(5,1,2) (4,6,0) %s (5,1,2) %s 6 1" % (False, False)]
  t = dmtbx.ext.triplet_generator(sg, i, None, 0, True, False)
  assert tuple(t.n_relations()) == (0,0)

def run():
  exercise_triplet_generator()
  print("OK")

if (__name__ == "__main__"):
  run()
