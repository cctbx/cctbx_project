from cctbx import dmtbx
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_triplet():
  sg = sgtbx.space_group_info("P 31")
  i = flex.miller_index(((1,2,3),(2,3,4)))
  e = flex.double((1,2))
  t = dmtbx.triplet_invariants(sg.type(), i, 00000, 00000)
  assert t.number_of_weighted_triplets() == 0
  assert t.total_number_of_triplets() == 0
  assert approx_equal(t.average_number_of_triplets_per_reflection(), 0)
  t.dump_triplets(i)
  assert t.n_relations(00000, 00000) == 0
  assert approx_equal(tuple(t.sum_of_e_products(i, e, 00000, 00000)), (0,0))
  s = flex.bool(2)
  o = flex.size_t()
  r = t.apply_tangent_formula(e, e, s, o, 00000, 00000, 00000)
  assert approx_equal(tuple(r), (1,2))

def run():
  exercise_triplet()
  print "OK"

if (__name__ == "__main__"):
  run()
