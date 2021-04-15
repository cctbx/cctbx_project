from __future__ import absolute_import, division, print_function
from iotbx import pdb
from mmtbx.conformation_dependent_library.testing_utils import get_geometry_restraints_manager

#This is a very basic test of generate_residue_tuples for RNA/DNA
#sample structure lines are taken from 1ehz.pdb and have been modified at
#  The P atom of residue 3 to introduce a chain break

pdblines = """ATOM      1  OP3   G A   1      50.193  51.190  50.534  1.00 99.85           O
ATOM      2  P     G A   1      50.626  49.730  50.573  1.00100.19           P
ATOM      3  OP1   G A   1      49.854  48.893  49.562  1.00100.19           O
ATOM      4  OP2   G A   1      52.137  49.542  50.511  1.00 99.21           O
ATOM      5  O5'   G A   1      50.161  49.136  52.023  1.00 99.82           O
ATOM      6  C5'   G A   1      50.216  49.948  53.210  1.00 98.63           C
ATOM      7  C4'   G A   1      50.968  49.231  54.309  1.00 97.84           C
ATOM      8  O4'   G A   1      50.450  47.888  54.472  1.00 97.10           O
ATOM      9  C3'   G A   1      52.454  49.030  54.074  1.00 98.07           C
ATOM     10  O3'   G A   1      53.203  50.177  54.425  1.00 99.39           O
ATOM     11  C2'   G A   1      52.781  47.831  54.957  1.00 96.96           C
ATOM     12  O2'   G A   1      53.018  48.156  56.313  1.00 96.77           O
ATOM     13  C1'   G A   1      51.502  47.007  54.836  1.00 95.70           C
ATOM     14  N9    G A   1      51.628  45.992  53.798  1.00 93.67           N
ATOM     15  C8    G A   1      51.064  46.007  52.547  1.00 92.60           C
ATOM     16  N7    G A   1      51.379  44.966  51.831  1.00 91.19           N
ATOM     17  C5    G A   1      52.197  44.218  52.658  1.00 91.47           C
ATOM     18  C6    G A   1      52.848  42.992  52.425  1.00 90.68           C
ATOM     19  O6    G A   1      52.826  42.291  51.404  1.00 90.38           O
ATOM     20  N1    G A   1      53.588  42.588  53.534  1.00 90.71           N
ATOM     21  C2    G A   1      53.685  43.282  54.716  1.00 91.21           C
ATOM     22  N2    G A   1      54.452  42.733  55.671  1.00 91.23           N
ATOM     23  N3    G A   1      53.077  44.429  54.946  1.00 91.92           N
ATOM     24  C4    G A   1      52.356  44.836  53.879  1.00 92.62           C
ATOM     25  P     C A   2      54.635  50.420  53.741  1.00100.19           P
ATOM     26  OP1   C A   2      55.145  51.726  54.238  1.00100.19           O
ATOM     27  OP2   C A   2      54.465  50.204  52.269  1.00100.19           O
ATOM     28  O5'   C A   2      55.563  49.261  54.342  1.00 98.27           O
ATOM     29  C5'   C A   2      55.925  49.246  55.742  1.00 95.40           C
ATOM     30  C4'   C A   2      56.836  48.075  56.049  1.00 93.33           C
ATOM     31  O4'   C A   2      56.122  46.828  55.830  1.00 92.18           O
ATOM     32  C3'   C A   2      58.090  47.947  55.197  1.00 92.75           C
ATOM     33  O3'   C A   2      59.174  48.753  55.651  1.00 92.89           O
ATOM     34  C2'   C A   2      58.416  46.463  55.298  1.00 91.81           C
ATOM     35  O2'   C A   2      59.140  46.136  56.466  1.00 91.36           O
ATOM     36  C1'   C A   2      57.022  45.836  55.356  1.00 90.59           C
ATOM     37  N1    C A   2      56.570  45.364  54.029  1.00 88.84           N
ATOM     38  C2    C A   2      57.094  44.157  53.520  1.00 88.64           C
ATOM     39  O2    C A   2      57.921  43.516  54.198  1.00 88.97           O
ATOM     40  N3    C A   2      56.686  43.721  52.301  1.00 87.36           N
ATOM     41  C4    C A   2      55.802  44.437  51.597  1.00 87.11           C
ATOM     42  N4    C A   2      55.430  43.972  50.397  1.00 86.30           N
ATOM     43  C5    C A   2      55.259  45.660  52.089  1.00 86.87           C
ATOM     44  C6    C A   2      55.663  46.080  53.296  1.00 88.01           C
ATOM     45  P     G A   3     600.184  49.419  54.574  1.00 92.31           P
ATOM     46  OP1   G A   3      61.015  50.422  55.295  1.00 92.97           O
ATOM     47  OP2   G A   3      59.371  49.857  53.404  1.00 91.56           O
ATOM     48  O5'   G A   3      61.137  48.219  54.105  1.00 88.57           O
ATOM     49  C5'   G A   3      62.175  47.724  54.969  1.00 83.44           C
ATOM     50  C4'   G A   3      62.769  46.443  54.422  1.00 79.87           C
ATOM     51  O4'   G A   3      61.734  45.427  54.299  1.00 78.36           O
ATOM     52  C3'   G A   3      63.405  46.499  53.040  1.00 78.97           C
ATOM     53  O3'   G A   3      64.741  47.029  53.060  1.00 79.76           O
ATOM     54  C2'   G A   3      63.359  45.032  52.608  1.00 77.19           C
ATOM     55  O2'   G A   3      64.411  44.256  53.155  1.00 77.80           O
ATOM     56  C1'   G A   3      62.018  44.572  53.194  1.00 73.98           C
ATOM     57  N9    G A   3      60.934  44.675  52.202  1.00 68.20           N
ATOM     58  C8    G A   3      60.024  45.702  52.050  1.00 65.03           C
ATOM     59  N7    G A   3      59.252  45.556  51.003  1.00 62.99           N
ATOM     60  C5    G A   3      59.655  44.348  50.447  1.00 59.95           C
ATOM     61  C6    G A   3      59.189  43.675  49.292  1.00 55.65           C
ATOM     62  O6    G A   3      58.287  44.013  48.522  1.00 53.32           O
ATOM     63  N1    G A   3      59.893  42.491  49.072  1.00 54.00           N
ATOM     64  C2    G A   3      60.906  42.006  49.876  1.00 55.46           C
ATOM     65  N2    G A   3      61.512  40.873  49.479  1.00 48.16           N
ATOM     66  N3    G A   3      61.312  42.605  50.983  1.00 56.69           N
ATOM     67  C4    G A   3      60.666  43.774  51.193  1.00 61.76           C"""

def main():
  pdb_inp = pdb.input(lines=pdblines, source_info='rna_test_1ehz')
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(raw_records=pdblines)
  pdb_hierarchy.reset_i_seq_if_necessary()
  from mmtbx.conformation_dependent_library import generate_residue_tuples

  #P of G 3 has been changed from x=60.184 to x=600.184
  #Check that this breaks connectivity
  rna_residue_pairs = []
  for pair in generate_residue_tuples(hierarchy=pdb_hierarchy,
                                    geometry=geometry_restraints_manager,
                                    length=2,
                                    dna_rna_residues=True,
                                    include_non_linked=False, #default behavior
                                    backbone_only=False,
                                    include_non_standard_residues=False,
                                    verbose=False):
    residue1 = pair[0]
    residue2 = pair[1]

    #sample output so you can see the pairs
    chain1 = residue1.parent().parent().id
    resname1 = residue1.resname
    resseq1 = residue1.resseq+residue1.icode
    resid1 = ":".join([chain1,resname1,resseq1])
    chain2 = residue2.parent().parent().id
    resname2 = residue2.resname
    resseq2 = residue2.resseq+residue2.icode
    resid2 = ":".join([chain2,resname2,resseq2])
    #print(residue1.id_str(), "---", residue2.id_str())
    #print(resid1,"---",resid2)
    rna_residue_pairs.append(resid1+"---"+resid2)
  assert len(rna_residue_pairs)==1

  #Check that generate residue_tuples can find two total pairs
  rna_residue_pairs = []
  for pair in generate_residue_tuples(hierarchy=pdb_hierarchy,
                                    geometry=geometry_restraints_manager,
                                    length=2,
                                    dna_rna_residues=True,
                                    include_non_linked=True, #chnged for this test
                                    backbone_only=False,
                                    include_non_standard_residues=False,
                                    verbose=False):
    residue1 = pair[0]
    residue2 = pair[1]

    #sample output so you can see the pairs
    chain1 = residue1.parent().parent().id
    resname1 = residue1.resname
    resseq1 = residue1.resseq+residue1.icode
    resid1 = ":".join([chain1,resname1,resseq1])
    chain2 = residue2.parent().parent().id
    resname2 = residue2.resname
    resseq2 = residue2.resseq+residue2.icode
    resid2 = ":".join([chain2,resname2,resseq2])
    #print(residue1.id_str(), "---", residue2.id_str())
    #print(resid1,"---",resid2)
    rna_residue_pairs.append(resid1+"---"+resid2)
  assert len(rna_residue_pairs)==2
  print("OK")

if __name__ == '__main__':
  main()
