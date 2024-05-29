from __future__ import absolute_import, division, print_function
from mmtbx.alignment import align
from mmtbx.alignment import amino_acid_codes, blosum62, dayhoff_mdm78_similarity_scores, \
    blosum50_similarity_scores, pairwise_global, dayhoff, blosum50, blosum62, \
    identity
import sys
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

import operator
import unittest
from functools import reduce
from six.moves import range

def exercise_align_mask():
  B="GCGAGATAAAGGGACCCATAAA" +"TGTCG"+ "TAGCATCGGGCTAATAGATAAGACACA"
  A=                          "TGTCG"+  "AGCATCGGGCTAATAGATAAGACACA"
  scores=[]
  for i in range(len(A)):
    masking_a=len(A)*[10]
    masking_a[i]=1
    obj = align(A,B,masking_a=masking_a)
    print("score=%.1f" % obj.score())
    alignment = obj.extract_alignment()
    print(alignment.match_codes)
    scores.append(obj.score())
  print(scores)
  assert scores==[2.0, 2.0, 3.0, 4.0, 5.0, 6.0, 5.0, 4.0, 3.0, 2.0, 2.0,
                  2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                  2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
  print("OK")
def exercise_align():
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

def exercise_exact_mismatch_selection():
  """ examples from exercise_aligh function
  """
  def get_mismatch_sel(s1, s2):
    i, j = align(s1, s2).extract_alignment().exact_mismatch_selection()
    return list(i), list(j)
  i_seqs, j_seqs = get_mismatch_sel("EASYA",
                                    "AETSYT")
  assert i_seqs == [1,4] and j_seqs == [2,5]

  i_seqs, j_seqs = get_mismatch_sel("EASY",
                                    "KMST")
  assert i_seqs == [0,1,3] and j_seqs == [0,1,3]

  i_seqs, j_seqs = get_mismatch_sel("EASY",
                                    "KMT")
  assert i_seqs == [1,2,3] and j_seqs == [0,1,2]
  i_seqs, j_seqs = get_mismatch_sel("EASY",
                                    "EASY")
  assert i_seqs == [] and j_seqs == []

def exercise_2():
  A = "AAAGGTT"
  B = "AAATT"
  obj = align(A,B)
  obj.show_matrices()

  print("score=%.1f" % obj.score())
  alignment = obj.extract_alignment()
  print(alignment.match_codes)
  print(alignment.a)
  print(alignment.identity_matches())
  print(alignment.b)

  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP".lower()
  obj = align(A,B,gap_opening_penalty=150,gap_extension_penalty=20,similarity_function="dayhoff",style="global")

  print("\n1rra vs. 1bli; GLOBAL allignment; mdm78")
  print("score=%.1f" % obj.score())
  alignment = obj.extract_alignment()

  print(alignment.match_codes)
  print(alignment.a)
  print(alignment.dayhoff_matches())
  print(alignment.b)
  assert approx_equal(alignment.calculate_sequence_identity(), 0.330645)


  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=150,gap_extension_penalty=20,similarity_function="dayhoff",style="local")

  print("\n1rra vs. 1bli; LOCAL allignment; mdm78")
  print("score=%.1f" % obj.score())
  alignment = obj.extract_alignment()

  print(alignment.match_codes)
  print(alignment.a)
  print(alignment.dayhoff_matches())
  print(alignment.b)
  assert approx_equal(alignment.calculate_sequence_identity(), 0.341880)



  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=10,gap_extension_penalty=2,similarity_function="blosum50",style="global")

  print("\n1rra vs. 1bli; GLOBAL allignment; blosum50")
  print("score=%.1f" % obj.score())
  alignment = obj.extract_alignment()

  print(alignment.match_codes)
  print(alignment.a)
  print(alignment.matches())
  print(alignment.b)
  assert approx_equal(alignment.calculate_sequence_identity(), 0.362903)

  # 1rra vs. 1bli
  A = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQAICSQGQVTCKNGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
  B = "DNSRYTHFLTQHYDAKPQGRDDRYCESIMRRRGLTSPCKDINTFIHGNKRSIKAICENKNGNPHRENLRISKSSFQVTTCKLHGGSPWPPCQYRATAGFRNVVVACENGLPVHLDQSIFRRP"
  obj = align(A,B,gap_opening_penalty=10,gap_extension_penalty=2,similarity_function="blosum50",style="local")

  print("\n1rra vs. 1bli; LOCAL allignment; blosum50")
  print("score=%.1f" % obj.score())
  alignment = obj.extract_alignment()

  print(alignment.match_codes)
  print(alignment.a)
  print(alignment.matches(similarity_function=blosum50, is_similar_threshold=0))
  print(alignment.b)
  assert approx_equal(alignment.calculate_sequence_identity(), 0.368852)
  print()
  alignment.pretty_print(
    matches = None,
    out = None,
    block_size = 50,
    n_block = 1,
    top_name = "1rra",
    bottom_name = "1bli",
    comment = """pretty_print is pretty pretty""")

  # example from PDB ID 2dex
  A = "GTLIRVTPEQPTHAVCVLGTLTQLDICSSAPXXXTSFSINASPGVVVDI"
  B = "GPLGSPEFMAQGTLIRVTPEQPTHAVCVLGTLTQLDICSSAPEDCTSFSINASPGVVVDI"
  obj = align(A, B, similarity_function="identity")
  alignment = obj.extract_alignment()
  assert alignment.match_codes == 'iiiiiiiiiiimmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm'
  assert alignment.a == '-----------GTLIRVTPEQPTHAVCVLGTLTQLDICSSAPXXXTSFSINASPGVVVDI'
  assert alignment.b == 'GPLGSPEFMAQGTLIRVTPEQPTHAVCVLGTLTQLDICSSAPEDCTSFSINASPGVVVDI'

  print("OK") # necessary for auto_build checking

def exercise_similarity_scores():
  from scitbx.array_family import flex
  for m in [dayhoff_mdm78_similarity_scores, blosum50_similarity_scores]:
    assert flex.double(m).matrix_is_symmetric(relative_epsilon=1e-15)

def exercise_ext():
  seq1="THEQUICKBOWNFOXJUMPSOVETHELAZY"
  seq2="QUICKBRWNFXJUMPSVERTH3LAZYDOG"
  pg = pairwise_global(seq1,seq2.lower())
  assert pg.result1 == "THEQUICKBOWNFOXJUMPSOVE-THELAZY---"
  assert pg.result2 == "---QUICKBRWNF-XJUMPS-VERTH3LAZYDOG"
  pg = pairwise_global_wrapper(seq1,seq2)
  assert pg.range_matches_from_aligned_sequences() == (
  [[3, 13], [14, 20], [21, 23], [23, 30]], [[0, 10], [10, 16], [16, 18], [19, 26]])
  assert ("%.2f" % pg.calculate_sequence_identity()) == "0.92"

class test_blosum62(unittest.TestCase):

  def testSymmetric(self):

    tests = reduce(
      operator.add,
      [
        [ ( amino_acid_codes[ l ], amino_acid_codes[ r ] )
          for r in range( l + 1 ) ]
        for l in range( len( amino_acid_codes ) )
        ]
      )

    for ( left, right ) in tests:
      self.assertEqual(
        blosum62( left, right ),
        blosum62( right, left )
        )


  def testSelected(self):

      for ( left, right, value ) in [
        ( "A", "C", 0 ), ( "E", "H", 0 ), ( "W", "W", 11 ), ( "F", "P", -4 )
        ]:
        self.assertEqual( blosum62( left, right ), value )


  def testUnknown(self):

    self.assertEqual( blosum62( "Q", "B" ), 0 )

suite_blosum62 = unittest.TestLoader().loadTestsFromTestCase(
  test_blosum62
  )

alltests = unittest.TestSuite(
  [
    suite_blosum62,
    ]
  )

def exercise():
  exercise_align_mask()
  exercise_align()
  exercise_exact_mismatch_selection()
  exercise_2()
  exercise_similarity_scores()
  exercise_ext
  sys.stdout.flush()
  unittest.TextTestRunner(stream=sys.stdout, verbosity = 2 ).run( alltests )

if (__name__ == "__main__"):
  exercise()
