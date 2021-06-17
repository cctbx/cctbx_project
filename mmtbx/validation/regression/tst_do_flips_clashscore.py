from __future__ import absolute_import, division, print_function
from mmtbx.validation.clashscore import clashscore
import libtbx.load_env
import iotbx.pdb
import unittest
import os
from six.moves import map

__author__ = 'Youval'

test_pdb_str = '''\
CRYST1   77.977   77.977   66.800  90.00  90.00 120.00 P 31 2 1      6
SCALE1      0.012824  0.007404  0.000000        0.00000
SCALE2      0.000000  0.014808  0.000000        0.00000
SCALE3      0.000000  0.000000  0.014970        0.00000
ATOM    489  N   TRP A  96      27.616 -26.119   6.863  1.00 11.38           N
ATOM    490  CA  TRP A  96      28.969 -26.184   6.322  1.00 10.73           C
ATOM    491  C   TRP A  96      29.155 -25.427   5.023  1.00 10.54           C
ATOM    492  O   TRP A  96      28.250 -25.321   4.198  1.00 10.74           O
ATOM    493  CB  TRP A  96      29.451 -27.607   6.139  1.00 10.00           C
ATOM    494  CG  TRP A  96      29.781 -28.324   7.405  1.00  9.99           C
ATOM    495  CD1 TRP A  96      29.155 -28.212   8.613  1.00 11.22           C
ATOM    496  CD2 TRP A  96      30.773 -29.342   7.558  1.00 10.57           C
ATOM    497  NE1 TRP A  96      29.720 -29.094   9.512  1.00 11.57           N
ATOM    498  CE2 TRP A  96      30.713 -29.795   8.890  1.00 10.98           C
ATOM    499  CE3 TRP A  96      31.720 -29.904   6.696  1.00 10.60           C
ATOM    500  CZ2 TRP A  96      31.547 -30.793   9.373  1.00 12.52           C
ATOM    501  CZ3 TRP A  96      32.560 -30.881   7.186  1.00 12.06           C
ATOM    502  CH2 TRP A  96      32.479 -31.313   8.512  1.00 11.72           C
ATOM    503  N   THR A  97      30.383 -24.936   4.864  1.00 10.88           N
ATOM    504  CA  THR A  97      30.827 -24.204   3.693  1.00 10.65           C
ATOM    505  C   THR A  97      31.766 -25.066   2.858  1.00 10.84           C
ATOM    506  O   THR A  97      32.223 -26.130   3.296  1.00 10.69           O
ATOM    507  CB  THR A  97      31.602 -22.960   4.143  1.00 11.07           C
ATOM    508  OG1 THR A  97      32.717 -23.404   4.923  1.00 12.40           O
ATOM    509  CG2 THR A  97      30.701 -22.026   4.977  1.00 12.69           C
ATOM    510  N   ILE A  98      32.069 -24.587   1.655  1.00 10.75           N
ATOM    511  CA  ILE A  98      33.014 -25.241   0.774  1.00 10.93           C
ATOM    512  C   ILE A  98      34.349 -25.612   1.442  1.00 11.69           C
ATOM    513  O   ILE A  98      34.775 -26.756   1.341  1.00 11.92           O
ATOM    514  CB  ILE A  98      33.229 -24.393  -0.515  1.00 10.59           C
ATOM    515  CG1 ILE A  98      31.994 -24.515  -1.415  1.00 10.44           C
ATOM    516  CG2 ILE A  98      34.467 -24.832  -1.258  1.00 11.99           C
ATOM    517  CD1 ILE A  98      31.865 -23.451  -2.472  1.00 11.19           C
ATOM    518  N   PRO A  99      35.028 -24.654   2.103  1.00 12.86           N
ATOM    519  CA  PRO A  99      36.312 -25.036   2.697  1.00 12.90           C
ATOM    520  C   PRO A  99      36.208 -26.072   3.817  1.00 12.70           C
ATOM    521  O   PRO A  99      37.131 -26.870   4.015  1.00 13.09           O
ATOM    522  CB  PRO A  99      36.890 -23.697   3.193  1.00 13.38           C
ATOM    523  CG  PRO A  99      35.777 -22.758   3.263  1.00 13.70           C
ATOM    524  CD  PRO A  99      34.770 -23.207   2.227  1.00 13.18           C
ATOM    525  N   GLN A 100      35.081 -26.100   4.520  1.00 12.56           N
ATOM    526  CA  GLN A 100      34.892 -27.115   5.545  1.00 12.83           C
ATOM    527  C   GLN A 100      34.722 -28.505   4.925  1.00 12.30           C
ATOM    528  O   GLN A 100      35.279 -29.478   5.404  1.00 12.28           O
ATOM    529  CB  GLN A 100      33.705 -26.765   6.421  1.00 13.53           C
ATOM    530  CG  GLN A 100      33.883 -25.509   7.218  1.00 16.43           C
ATOM    531  CD  GLN A 100      32.648 -25.287   8.021  1.00 19.16           C
ATOM    532  NE2 GLN A 100      31.646 -24.831   7.508  1.00 17.90           O
ATOM    533  OE1 GLN A 100      32.675 -25.714   9.280  1.00 24.40           N
'''


class MyTestCase(unittest.TestCase):

  def setUp(self):
    self.file_to_delete = []
    # import files used in tests
    self.file_name = 'test_do_flips_clashscore.pdb'
    with open(self.file_name,'w') as f:
      f.write(test_pdb_str)
    self.file_to_delete.append(self.file_name)

  def test_identifying_and_addition_of_hydrogen(self):
    """ test identifying and addition of hydrogen """
    has_reduce = libtbx.env.has_module(name="reduce")
    if has_reduce:
      pdb_inp = iotbx.pdb.input(file_name=self.file_name)
      pdb_hierarchy = pdb_inp.construct_hierarchy()

      # don't do flip
      result = clashscore(
        pdb_hierarchy=pdb_hierarchy,
        keep_hydrogens=False)#verbose=True)

      self.assertAlmostEqual(result.clashscore,22.9885057471,places=4)

      # do flip
      result = clashscore(
        pdb_hierarchy=pdb_hierarchy,
        keep_hydrogens=False,#verbose=True)
        do_flips=True)

      self.assertEqual(result.clashscore,0)
    else:
      # Skip test if reduce is not present
      pass

  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_something']
  suite = unittest.TestSuite(list(map(MyTestCase, tests)))
  return suite


if __name__ == '__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
