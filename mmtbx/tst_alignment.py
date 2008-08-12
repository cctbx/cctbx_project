from mmtbx.alignment import align
import sys
from scitbx.array_family import flex

def exercise():
  #
  i_seqs, j_seqs = align("EASYA",
                         "AETSYT").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 2, 3]) and j_seqs == flex.size_t([1, 3, 4])
  #
  i_seqs, j_seqs = align("AAAGGTT",
                         "AAATT").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 1, 2, 5, 6]) and \
         j_seqs == flex.size_t([0, 1, 2, 3, 4])
  #
  i_seqs, j_seqs = align("AESD",
                         "AEDK").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 1]) and j_seqs == flex.size_t([0, 1])
  #
  i_seqs, j_seqs = align("EASY",
                         "YSAE").extract_alignment().exact_match_selections()
  assert i_seqs.size()==0 and j_seqs.size()==0
  #
  i_seqs, j_seqs = align("EASY",
                         "KMT").extract_alignment().exact_match_selections()
  assert i_seqs.size()==0 and j_seqs.size()==0
  #
  i_seqs, j_seqs = align("EASY",
                         "KMST").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([2]) and j_seqs == flex.size_t([2])
  #
  i_seqs, j_seqs = align("EASY",
                         "KMMST").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([2]) and j_seqs == flex.size_t([3])
  #
  i_seqs, j_seqs = align("EASY",
                         "EATY").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 1, 3]) and j_seqs == flex.size_t([0, 1, 3])
  #
  i_seqs, j_seqs = align("EASIEST",
                         "EATY").extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 1]) and j_seqs == flex.size_t([0, 1])
  #
  i_seqs, j_seqs = align("EEEEEASIEST",
                         "EEEATYRRIESQQEIES"
                        ).extract_alignment().exact_match_selections()
  assert i_seqs == flex.size_t([0, 1, 2, 7, 8, 9]) and \
         j_seqs == flex.size_t([0, 1, 2, 8, 9, 10])


if (__name__ == "__main__"):
  exercise()
