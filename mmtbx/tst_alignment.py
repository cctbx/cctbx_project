from mmtbx.alignment import align
from mmtbx.alignment import amino_acid_codes, blosum62
import sys
from scitbx.array_family import flex

import operator
import unittest

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
  exercise_align()
  sys.stdout.flush()
  unittest.TextTestRunner(stream=sys.stdout, verbosity = 2 ).run( alltests )

if (__name__ == "__main__"):
  exercise()
