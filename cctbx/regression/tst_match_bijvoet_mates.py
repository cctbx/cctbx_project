from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected
import random

def reference_match(sg_type, indices, assert_is_unique_set_under_symmetry):
  """
  Transcription of the std::map algorithm that match_bijvoet_mates::match_
  used up to cctbx ad7fdc65b9. The C++ implementation must reproduce its
  pairs and singles exactly, in the same order, including the behaviour on
  duplicate indices (the last occurrence wins the lookup).
  """
  asu = sgtbx.reciprocal_space_asu(sg_type)
  lookup = {}
  for i, h in enumerate(indices):
    if assert_is_unique_set_under_symmetry and h in lookup:
      raise RuntimeError("miller array is not a unique set under symmetry")
    lookup[h] = i
  paired_already = [False] * len(indices)
  pairs = []
  singles_plus = []
  singles_minus = []
  for i, h in enumerate(indices):
    if paired_already[i]:
      continue
    if h == (0, 0, 0):
      singles_plus.append(i)
      continue
    asu_which = asu.which(h)
    assert asu_which != 0
    j = lookup.get((-h[0], -h[1], -h[2]))
    if j is None:
      if asu_which > 0:
        singles_plus.append(i)
      else:
        singles_minus.append(i)
    else:
      if asu_which > 0:
        pairs.append((i, j))
      else:
        pairs.append((j, i))
      paired_already[j] = True
  return pairs, singles_plus, singles_minus

def check(sg_type, indices, assert_is_unique_set_under_symmetry=True,
          compare_singles=True):
  pairs, singles_plus, singles_minus = reference_match(
    sg_type, list(indices), assert_is_unique_set_under_symmetry)
  matches = miller.match_bijvoet_mates(
    sg_type, indices, assert_is_unique_set_under_symmetry)
  assert list(matches.pairs()) == pairs, \
    "%d pairs from C++, %d from the reference" % (
      matches.pairs().size(), len(pairs))
  assert matches.n_singles() == len(singles_plus) + len(singles_minus)
  if not compare_singles:
    return len(pairs)
  assert list(matches.singles("+")) == singles_plus
  assert list(matches.singles("-")) == singles_minus
  assert list(matches.pairs_hemisphere_selection("+")) == [p[0] for p in pairs]
  assert list(matches.pairs_hemisphere_selection("-")) == [p[1] for p in pairs]
  return len(pairs)

def asu_indices(space_group_symbol, anomalous_flag, d_min):
  space_group_info = sgtbx.space_group_info(space_group_symbol)
  crystal_symmetry = space_group_info.any_compatible_crystal_symmetry(
    asu_volume=1500)
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=anomalous_flag,
    d_min=d_min).map_to_asu()
  return space_group_info.type(), miller_set.indices()

def shuffled(indices, rng):
  perm = list(range(indices.size()))
  rng.shuffle(perm)
  return indices.select(flex.size_t(perm))

def exercise_generated_sets():
  rng = random.Random(0)
  for symbol in ["P 1", "P -1", "P 21 21 21", "P 63", "I 41 3 2"]:
    for anomalous_flag in [False, True]:
      sg_type, indices = asu_indices(symbol, anomalous_flag, d_min=2.0)
      assert indices.size() > 100, (symbol, indices.size())
      n_pairs = check(sg_type, indices)
      if anomalous_flag and symbol != "P -1":
        assert n_pairs > 0, symbol
      if not anomalous_flag:
        assert n_pairs == 0, (symbol, n_pairs)
      # input order must not matter for correctness and must be preserved
      check(sg_type, shuffled(indices, rng))
      # a random half of the minus hemisphere removed: pairs and minus
      # singles mixed
      keep = flex.bool(indices.size(), True)
      asu = sgtbx.reciprocal_space_asu(sg_type)
      for i, h in enumerate(indices):
        if asu.which(h) < 0 and rng.random() < 0.5:
          keep[i] = False
      check(sg_type, indices.select(keep))

def exercise_special_cases():
  sg_type = sgtbx.space_group_info("P 21 21 21").type()
  # empty set
  check(sg_type, flex.miller_index())
  # the zero index is a plus single, never paired with itself
  check(sg_type, flex.miller_index([(0, 0, 0)]))
  check(sg_type, flex.miller_index([(1, 2, 3), (0, 0, 0), (-1, -2, -3)]))
  # a single pair, both orders of appearance
  check(sg_type, flex.miller_index([(1, 2, 3), (-1, -2, -3)]))
  check(sg_type, flex.miller_index([(-1, -2, -3), (1, 2, 3)]))
  # P1: every index is acentric and its own asu is a hemisphere
  check(sgtbx.space_group_info("P 1").type(),
        flex.miller_index([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, 0, -1)]))

def exercise_duplicates():
  rng = random.Random(1)
  sg_type, indices = asu_indices("P 21 21 21", anomalous_flag=True, d_min=3.0)
  extra = shuffled(indices, rng)[:indices.size() // 3]
  with_duplicates = indices.concatenate(extra)
  # non-unique sets are accepted when the assertion is off (the code path
  # of miller.set.n_bijvoet_pairs), with the old lookup semantics: a
  # duplicate pairs again with the last occurrence of its mate, so
  # size_processed() exceeds the input size and singles() fails its
  # intrinsic size assertion; only pairs() and n_singles() are defined
  check(sg_type, with_duplicates, assert_is_unique_set_under_symmetry=False,
        compare_singles=False)
  check(sg_type, shuffled(with_duplicates, rng),
        assert_is_unique_set_under_symmetry=False, compare_singles=False)
  # and rejected with the assertion on
  try:
    miller.match_bijvoet_mates(sg_type, with_duplicates, True)
  except RuntimeError as e:
    assert str(e).find("miller array is not a unique set under symmetry") >= 0, \
      str(e)
  else:
    raise Exception_expected

def exercise_miller_set_interface():
  # n_bijvoet_pairs (the reflection readers' auto_anomalous) and the
  # set-level match go through the same C++ class
  space_group_info = sgtbx.space_group_info("P 21 21 21")
  crystal_symmetry = space_group_info.any_compatible_crystal_symmetry(
    asu_volume=1500)
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry, anomalous_flag=True, d_min=2.0)
  asu = miller_set.map_to_asu()
  pairs, singles_plus, singles_minus = reference_match(
    space_group_info.type(), list(asu.indices()), False)
  assert miller_set.n_bijvoet_pairs() == len(pairs)
  assert miller_set.auto_anomalous().anomalous_flag()
  asu_set, matches = miller_set.match_bijvoet_mates()
  assert list(matches.pairs()) == pairs
  assert list(matches.singles("-")) == singles_minus

def run():
  exercise_generated_sets()
  exercise_special_cases()
  exercise_duplicates()
  exercise_miller_set_interface()
  print("OK")

if __name__ == "__main__":
  run()
