from __future__ import absolute_import, division, print_function
import iotbx.ncs as ncs
import iotbx.pdb
from iotbx.ncs import ncs_group_master_phil
import iotbx.phil

pdb_str_1 = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     23  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     24  CA AGLN A   4       0.384   1.888   3.199  0.50 10.53           C
ATOM     25  CB AGLN A   4       0.656   2.148   1.711  0.50  9.80           C
ATOM     26  CG AGLN A   4       1.944   1.458   1.213  0.50 10.25           C
ATOM     27  CD AGLN A   4       2.504   2.044  -0.089  0.50 12.43           C
ATOM     28  OE1AGLN A   4       2.744   3.268  -0.190  0.50 14.62           O
ATOM     29  NE2AGLN A   4       2.750   1.161  -1.091  0.50  9.05           N
ATOM     30  CA BGLN A   4       0.484   1.888   3.199  0.50 10.53           C
ATOM     31  CB BGLN A   4       0.756   2.148   1.711  0.50  9.80           C
ATOM     32  CG BGLN A   4      -0.021   1.186   0.788  0.50 10.25           C
ATOM     33  CD BGLN A   4       0.170   1.463  -0.708  0.50 12.43           C
ATOM     34  OE1BGLN A   4       0.485   2.604  -1.113  0.50 14.62           O
ATOM     35  NE2BGLN A   4       0.013   0.404  -1.543  0.50  9.05           N
ATOM     36  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     37  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     38  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     39  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     40  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     41  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     42  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     43  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     44  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     45  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     46  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     47  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     48  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     49  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     50  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     51  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     52  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     53  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     54  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     55  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     56  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     57  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     58  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     59  CD1 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     60  CD2 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     61  CE1 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     62  CE2 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     63  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     64  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     65  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER
HETATM   66  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   67  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   68  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   69  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   70  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   71  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   72  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
ATOM      1  N   GLY B   1       9.009   7.045  -6.102  1.00 16.77           N
ATOM      2  CA  GLY B   1       9.052   6.640  -4.651  1.00 16.57           C
ATOM      3  C   GLY B   1       8.015   5.573  -4.419  1.00 16.16           C
ATOM      4  O   GLY B   1       7.523   4.954  -5.381  1.00 16.78           O
ATOM      5  N   ASN B   2       7.656   5.356  -3.155  1.00 15.02           N
ATOM      6  CA  ASN B   2       6.522   4.471  -2.831  1.00 14.10           C
ATOM      7  C   ASN B   2       5.241   4.970  -3.427  1.00 13.13           C
ATOM      8  O   ASN B   2       4.978   6.175  -3.426  1.00 11.91           O
ATOM      9  CB  ASN B   2       6.346   4.314  -1.341  1.00 15.38           C
ATOM     10  CG  ASN B   2       7.584   3.775  -0.692  1.00 14.08           C
ATOM     11  OD1 ASN B   2       8.025   2.660  -1.016  1.00 17.46           O
ATOM     12  ND2 ASN B   2       8.204   4.588   0.169  1.00 11.72           N
ATOM     13  N   ASN B   3       4.438   4.023  -3.905  1.00 12.26           N
ATOM     14  CA  ASN B   3       3.193   4.337  -4.589  1.00 11.74           C
ATOM     15  C   ASN B   3       1.955   3.765  -3.895  1.00 11.10           C
ATOM     16  O   ASN B   3       1.872   2.552  -3.648  1.00 10.42           O
ATOM     17  CB  ASN B   3       3.259   3.811  -6.042  1.00 12.15           C
ATOM     18  CG  ASN B   3       2.006   4.172  -6.861  1.00 12.82           C
ATOM     19  OD1 ASN B   3       1.702   5.358  -7.072  1.00 15.05           O
ATOM     20  ND2 ASN B   3       1.271   3.148  -7.306  1.00 13.48           N
ATOM     21  N   GLN B   4       1.005   4.661  -3.598  1.00 10.29           N
ATOM     22  C   GLN B   4      -1.435   5.039  -4.088  1.00 10.24           C
ATOM     23  O   GLN B   4      -1.547   6.276  -4.115  1.00  8.86           O
ATOM     24  CA AGLN B   4      -0.384   4.321  -3.199  0.50 10.53           C
ATOM     25  CB AGLN B   4      -0.656   4.581  -1.711  0.50  9.80           C
ATOM     26  CG AGLN B   4      -1.944   3.891  -1.213  0.50 10.25           C
ATOM     27  CD AGLN B   4      -2.504   4.477   0.089  0.50 12.43           C
ATOM     28  OE1AGLN B   4      -2.744   5.701   0.190  0.50 14.62           O
ATOM     29  NE2AGLN B   4      -2.750   3.594   1.091  0.50  9.05           N
ATOM     30  CA BGLN B   4      -0.484   4.321  -3.199  0.50 10.53           C
ATOM     31  CB BGLN B   4      -0.756   4.581  -1.711  0.50  9.80           C
ATOM     32  CG BGLN B   4       0.021   3.619  -0.788  0.50 10.25           C
ATOM     33  CD BGLN B   4      -0.170   3.896   0.708  0.50 12.43           C
ATOM     34  OE1BGLN B   4      -0.485   5.037   1.113  0.50 14.62           O
ATOM     35  NE2BGLN B   4      -0.013   2.837   1.543  0.50  9.05           N
ATOM     36  N   GLN B   5      -2.154   4.254  -4.871  1.00 10.38           N
ATOM     37  CA  GLN B   5      -3.270   4.794  -5.640  1.00 11.39           C
ATOM     38  C   GLN B   5      -4.594   4.201  -5.172  1.00 11.52           C
ATOM     39  O   GLN B   5      -4.768   2.979  -5.054  1.00 12.05           O
ATOM     40  CB  GLN B   5      -3.056   4.616  -7.147  1.00 11.96           C
ATOM     41  CG  GLN B   5      -1.829   5.383  -7.647  1.00 10.81           C
ATOM     42  CD  GLN B   5      -1.344   4.847  -8.954  1.00 13.10           C
ATOM     43  OE1 GLN B   5      -0.774   3.758  -9.002  1.00 10.65           O
ATOM     44  NE2 GLN B   5      -1.549   5.620 -10.039  1.00 12.30           N
ATOM     45  N   ASN B   6      -5.514   5.097  -4.856  1.00 11.99           N
ATOM     46  CA  ASN B   6      -6.831   4.743  -4.318  1.00 12.30           C
ATOM     47  C   ASN B   6      -7.854   5.194  -5.324  1.00 13.40           C
ATOM     48  O   ASN B   6      -8.219   6.376  -5.374  1.00 13.92           O
ATOM     49  CB  ASN B   6      -7.065   5.449  -2.993  1.00 12.13           C
ATOM     50  CG  ASN B   6      -5.961   5.168  -2.003  1.00 12.77           C
ATOM     51  OD1 ASN B   6      -5.798   4.037  -1.551  1.00 14.27           O
ATOM     52  ND2 ASN B   6      -5.195   6.180  -1.679  1.00 10.07           N
ATOM     53  N   TYR B   7      -8.292   4.250  -6.147  1.00 14.70           N
ATOM     54  CA  TYR B   7      -9.159   4.577  -7.299  1.00 15.18           C
ATOM     55  C   TYR B   7     -10.603   4.764  -6.885  1.00 15.91           C
ATOM     56  O   TYR B   7     -11.041   4.244  -5.855  1.00 15.76           O
ATOM     57  CB  TYR B   7      -9.061   3.498  -8.369  1.00 15.35           C
ATOM     58  CG  TYR B   7      -7.665   3.362  -8.902  1.00 14.45           C
ATOM     59  CD1 TYR B   7      -7.210   4.189  -9.920  1.00 14.80           C
ATOM     60  CD2 TYR B   7      -6.771   2.454  -8.327  1.00 15.68           C
ATOM     61  CE1 TYR B   7      -5.904   4.082 -10.416  1.00 14.33           C
ATOM     62  CE2 TYR B   7      -5.480   2.339  -8.796  1.00 13.46           C
ATOM     63  CZ  TYR B   7      -5.047   3.162  -9.831  1.00 15.09           C
ATOM     64  OH  TYR B   7      -3.766   3.022 -10.291  1.00 14.39           O
ATOM     65  OXT TYR B   7     -11.358   5.432  -7.612  1.00 17.49           O
TER
HETATM   66  O   HOH B   8       6.471   7.660  -7.124  1.00 22.62           O
HETATM   67  O   HOH B   9     -10.431   4.291  -3.216  1.00 19.71           O
HETATM   68  O   HOH B  10      11.286   4.189   1.468  1.00 17.08           O
HETATM   69  O   HOH B  11     -11.808   6.612  -9.970  1.00 23.99           O
HETATM   70  O   HOH B  12     -13.605   3.760  -9.198  1.00 26.17           O
HETATM   71  O   HOH B  13       2.749   5.862 -10.024  1.00 39.15           O
HETATM   72  O   HOH B  14       1.500   3.115 -10.967  1.00 43.49           O
END
"""


def exercise_1():
  """
  Pure NCS search with AC and water that should be dropped
  """
  h = iotbx.pdb.input(lines=pdb_str_1, source_info=None).construct_hierarchy()
  ncs_inp = ncs.input(hierarchy = h)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  assert ncs_groups[0].master_iselection.size() == 50, ncs_groups[0].master_iselection.size()
  assert ncs_groups[0].master_iselection.size() == ncs_groups[0].copies[0].iselection.size()

def exercise_2():
  """
  Same as 1, but include water
  """
  search_params = ncs.input.get_default_params()
  search_params.ncs_search.exclude_selection = "element H or element D"
  h = iotbx.pdb.input(lines=pdb_str_1, source_info=None).construct_hierarchy()
  ncs_inp = ncs.input(
      hierarchy = h,
      params=search_params.ncs_search)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  assert ncs_groups[0].master_iselection.size() == 57, ncs_groups[0].master_iselection.size()
  assert ncs_groups[0].master_iselection.size() == ncs_groups[0].copies[0].iselection.size()

def exercise_3():
  """
  Same as 2, but also provide groups and don't check them allowing AC to be
  included
  """

  phil_str="""
ncs_group {
  reference = chain A
  selection = chain B
}
"""

  search_params = ncs.input.get_default_params()
  search_params.ncs_search.exclude_selection = "element H or element D"
  search_params.ncs_search.validate_user_supplied_groups = False
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  h = iotbx.pdb.input(lines=pdb_str_1, source_info=None).construct_hierarchy()
  ncs_inp = ncs.input(
      hierarchy = h,
      params=search_params.ncs_search,
      ncs_phil_groups=phil_groups.ncs_group)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  assert ncs_groups[0].master_iselection.size() == 72, ncs_groups[0].master_iselection.size()
  assert ncs_groups[0].master_iselection.size() == ncs_groups[0].copies[0].iselection.size()
  assert h.atoms_size() == 72*2, h.atoms_size() # covers whole model

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  exercise_3()
  print("OK")
