from cctbx import dmtbx
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_triplet():
  sg = sgtbx.space_group_info("P 31")
  i = flex.miller_index(((1,2,3),(2,3,4)))
  e = flex.double((1,2))
  t = dmtbx.triplet_invariants(sg.type(), i, False, False)
  assert t.number_of_weighted_triplets() == 0
  assert t.total_number_of_triplets() == 0
  assert approx_equal(t.average_number_of_triplets_per_reflection(), 0)
  assert t.n_relations(0) == 0
  t.dump_triplets(i)
  r = t.apply_tangent_formula(e, e)
  assert approx_equal(tuple(r), (1,2))
  t.weights_and_epsilon(sg.type(), i)

def run():
  exercise_triplet()
  print "OK"

if (__name__ == "__main__"):
  run()
