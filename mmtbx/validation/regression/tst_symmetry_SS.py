from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import unittest
import os

__author__ = 'bhintze'

test_pdb_str = '''\
CRYST1   89.760   89.760  150.890  90.00  90.00  90.00 P 43 21 2
SCALE1      0.011141  0.000000  0.000000        0.00000
SCALE2     -0.000000  0.011141  0.000000        0.00000
SCALE3      0.000000 -0.000000  0.006627        0.00000
ATOM    376  N   ASP A 532       2.906  76.529 186.658  1.00134.73           N
ATOM    377  CA  ASP A 532       4.198  77.174 186.751  1.00132.33           C
ATOM    378  CB  ASP A 532       4.325  77.971 188.064  1.00131.20           C
ATOM    379  CG  ASP A 532       5.669  78.695 188.203  1.00117.21           C
ATOM    380  OD1 ASP A 532       6.714  78.045 188.046  1.00109.08           O
ATOM    381  OD2 ASP A 532       5.685  79.910 188.487  1.00112.03           O
ATOM    382  C   ASP A 532       5.202  76.022 186.702  1.00132.62           C
ATOM    383  O   ASP A 532       5.352  75.246 187.645  1.00149.81           O
ATOM    384  N   ILE A 533       5.929  75.960 185.604  1.00147.69           N
ATOM    385  CA  ILE A 533       6.846  74.859 185.312  1.00146.43           C
ATOM    386  CB  ILE A 533       7.535  75.096 183.938  1.00154.91           C
ATOM    387  CG1 ILE A 533       6.486  75.039 182.797  1.00159.66           C
ATOM    388  CD1 ILE A 533       5.871  76.379 182.402  1.00149.38           C
ATOM    389  CG2 ILE A 533       8.693  74.126 183.701  1.00154.51           C
ATOM    390  C   ILE A 533       7.882  74.683 186.433  1.00137.78           C
ATOM    391  O   ILE A 533       8.304  73.551 186.741  1.00123.42           O
ATOM    392  N   CYS A 534       8.263  75.797 187.056  1.00131.42           N
ATOM    393  CA  CYS A 534       9.333  75.808 188.057  1.00131.90           C
ATOM    394  CB  CYS A 534       9.810  77.251 188.285  1.00130.90           C
ATOM    395  SG  CYS A 534      10.891  77.555 189.708  1.00138.71           S
ATOM    396  C   CYS A 534       8.938  75.168 189.395  1.00139.83           C
ATOM    397  O   CYS A 534       9.692  74.343 189.945  1.00123.05           O
ATOM    398  N   LYS A 535       7.769  75.567 189.904  1.00132.70           N
ATOM    399  CA  LYS A 535       7.276  75.128 191.195  1.00113.57           C
ATOM    400  CB  LYS A 535       6.321  76.171 191.766  1.00104.76           C
ATOM    401  CG  LYS A 535       7.063  77.459 192.082  1.00103.62           C
ATOM    402  CD  LYS A 535       6.239  78.489 192.831  1.00111.95           C
ATOM    403  CE  LYS A 535       7.178  79.531 193.454  1.00118.75           C
ATOM    404  NZ  LYS A 535       6.522  80.730 194.049  1.00109.53           N
ATOM    405  C   LYS A 535       6.621  73.770 191.083  1.00112.61           C
ATOM    406  O   LYS A 535       6.876  72.888 191.904  1.00127.98           O
ATOM    407  N   GLN A 536       5.823  73.580 190.039  1.00115.26           N
ATOM    408  CA  GLN A 536       5.038  72.353 189.884  1.00113.55           C
ATOM    409  CB  GLN A 536       3.702  72.682 189.204  1.00 93.32           C
ATOM    410  C   GLN A 536       5.808  71.230 189.153  1.00116.48           C
ATOM    411  O   GLN A 536       5.180  70.277 188.684  1.00123.87           O
'''

class TestMPGeo(unittest.TestCase):

  def setUp(self):
    self.file_to_delete = ['try.geo']
    # import files used in tests
    self.file_name = 'symmetry_SS.pdb'
    with open(self.file_name,'w')as f:
      f.write(test_pdb_str)
    self.file_to_delete.append(self.file_name)

  def test_mpgeo(self):
    cmd = "mmtbx.mp_geo pdb=symmetry_SS.pdb "
    cmd+= "out_file=try.geo kinemage=True"
    er = easy_run.call(command=cmd)
    #er.show_stderr()
    assert er == 0, 'command "%s" failed' % cmd

  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)

if __name__ == '__main__' :
    unittest.main()
