from __future__ import absolute_import, division, print_function
from mmtbx.geometry_restraints.afitt import finite_difference_test

pdb_lines = """
CRYST1   94.100   94.100  131.400  90.00  90.00 120.00 P 61 2 2
ATOM    890  N   ALA E 113      32.525  39.841  -1.660  1.00 15.64           N
ATOM    891  CA  ALA E 113      31.553  40.901  -1.907  1.00 15.53           C
ATOM    893  C   ALA E 113      31.742  41.398  -3.334  1.00 15.47           C
ATOM    894  O   ALA E 113      32.873  41.486  -3.811  1.00 15.46           O
ATOM    892  CB  ALA E 113      31.740  42.036  -0.918  1.00 15.53           C
ATOM   2602 ZN    ZN F   5      36.912  44.685  -7.131  1.00 13.76          ZN
ATOM   2606  P1  LEP X   1      35.374  42.336  -6.280  1.00 14.03           P
ATOM   2607  O2  LEP X   1      34.629  43.616  -5.992  1.00 14.13           O
ATOM   2608  O3  LEP X   1      36.737  42.562  -6.884  1.00 13.72           O
ATOM   2609  O4  LEP X   1      34.549  41.267  -6.954  1.00 14.26           O
ATOM   2610  N5  LEP X   1      35.736  41.643  -4.701  1.00 14.42           N
ATOM   2611  C6  LEP X   1      36.701  42.323  -3.815  1.00 14.69           C
ATOM   2612  C7  LEP X   1      37.996  41.614  -4.069  1.00 14.79           C
ATOM   2613  O8  LEP X   1      39.062  42.207  -4.117  1.00 15.05           O
ATOM   2614  C9  LEP X   1      36.278  42.128  -2.373  1.00 14.85           C
ATOM   2615  C10 LEP X   1      37.233  42.793  -1.403  1.00 15.12           C
ATOM   2616  C11 LEP X   1      37.112  44.297  -1.488  1.00 15.24           C
ATOM   2617  C12 LEP X   1      36.866  42.344  -0.010  1.00 15.23           C
ATOM   2618  N13 LEP X   1      37.835  40.262  -4.227  1.00 14.73           N
ATOM   2632  H1  LEP X   1      35.374  42.336  -6.280  1.00 20.00           H
ATOM   2633  H2  LEP X   1      35.163  44.358  -6.248  1.00 20.00           H
ATOM   2619  H5  LEP X   1      34.776  41.783  -4.384  1.00 14.42           H
ATOM   2620  H6  LEP X   1      36.777  43.375  -4.099  1.00 14.69           H
ATOM   2623  H10 LEP X   1      38.271  42.498  -1.590  1.00 15.12           H
ATOM   2624 H111 LEP X   1      37.719  44.781  -0.714  1.00 15.24           H
ATOM   2625 H112 LEP X   1      37.468  44.680  -2.449  1.00 15.24           H
ATOM   2626 H113 LEP X   1      36.078  44.634  -1.353  1.00 15.24           H
ATOM   2627 H121 LEP X   1      36.971  41.259   0.097  1.00 15.23           H
ATOM   2628 H122 LEP X   1      37.521  42.806   0.737  1.00 15.23           H
ATOM   2629 H123 LEP X   1      35.834  42.607   0.247  1.00 15.23           H
ATOM   2630 H131 LEP X   1      38.608  39.742  -4.600  1.00 14.73           H
ATOM   2631 H132 LEP X   1      36.892  40.088  -4.602  1.00 14.73           H
ATOM   2621  H92 LEP X   1      35.260  42.517  -2.233  1.00 14.85           H
ATOM   2622  H93 LEP X   1      36.213  41.054  -2.148  1.00 14.85           H
"""
cif_lines = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
LEP LEP '(2S)-2-[(hydroxy-oxido-oxo-$l^{6}-phosphanyl)amino]-4-methyl-pentanamide'  non-polymer   28   13
data_comp_LEP
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
LEP P1     P 'P1  '  0.000
LEP O2     O 'OP  '  0.000
LEP O3     O 'OP  '  0.000
LEP O4     O 'OP  '  0.000
LEP N5     N 'NT1 '  0.000
LEP C6     C 'CH1 '  0.000
LEP C7     C 'C   '  0.000
LEP O8     O 'O   '  0.000
LEP C9     C 'CH2 '  0.000
LEP C10    C 'CH1 '  0.000
LEP C11    C 'CH3 '  0.000
LEP C12    C 'CH3 '  0.000
LEP N13    N 'NH2 '  0.000
LEP H1     H 'HCH '  0.000
LEP H2     H 'HOH1'  0.000
LEP H5     H 'HNH1'  0.000
LEP H6     H 'HCH1'  0.000
LEP H10    H 'HCH1'  0.000
LEP H111   H 'HCH3'  0.000
LEP H112   H 'HCH3'  0.000
LEP H113   H 'HCH3'  0.000
LEP H121   H 'HCH3'  0.000
LEP H122   H 'HCH3'  0.000
LEP H123   H 'HCH3'  0.000
LEP H131   H 'HNH2'  0.000
LEP H132   H 'HNH2'  0.000
LEP H92    H 'HCH2'  0.000
LEP H93    H 'HCH2'  0.000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LEP     P1      O2 single 1.630 0.020
LEP     P1      O3 single 1.510 0.020
LEP     P1      O4 double 1.510 0.015
LEP     P1      N5 single 1.762 0.020
LEP     N5      C6 single 1.472 0.020
LEP     N5      H5 single 1.028 0.020
LEP     C6      C7 single 1.492 0.020
LEP     C6      C9 single 1.508 0.020
LEP     C6      H6 single 1.093 0.020
LEP     C7      O8 double 1.222 0.015
LEP     C7     N13 single 1.369 0.020
LEP     C9     C10 single 1.508 0.020
LEP     C9     H92 single 1.093 0.020
LEP     C9     H93 single 1.093 0.020
LEP    C10     C11 single 1.508 0.020
LEP    C10     C12 single 1.508 0.020
LEP    C10     H10 single 1.093 0.020
LEP    C11    H111 single 1.093 0.020
LEP    C11    H112 single 1.093 0.020
LEP    C11    H113 single 1.093 0.020
LEP    C12    H121 single 1.093 0.020
LEP    C12    H122 single 1.093 0.020
LEP    C12    H123 single 1.093 0.020
LEP    N13    H131 single 1.015 0.020
LEP    N13    H132 single 1.015 0.020
LEP     P1      H1 single 1.411 0.020
LEP     O2      H2 single 0.981 0.020
"""

def run(pdb_file, cif_file, ligand_names,atom,scale,verbose):
  finite_difference_test(pdb_file, cif_file, ligand_names,atom,scale,verbose)

if (__name__ == "__main__"):
  import sys
  if len(sys.argv)==1:
    f=file("afitt_fd.pdb", "wb")
    f.write(pdb_lines)
    f.close()
    f=file("afitt_fd.cif", "wb")
    f.write(cif_lines)
    f.close()
    run("afitt_fd.pdb", 'afitt_fd.cif', ["LEP"], 7, 1, True)

  else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="pdb file")
    parser.add_argument("cif_file", help="cif file", default=0)
    parser.add_argument("ligand_names", help="3-letter ligand names separated by commas")
    parser.add_argument("atom", help="i_seq of atom whose x-coord will be "
                                     "modified (to test AFITT this should "
                                     "be a ligand atom)")
    parser.add_argument( "-s", "--scale", default=1, help="weight applied" \
                        "to AFITT gradient/traget", type=float)
    parser.add_argument('-v', dest='verbose', action='store_true', help="verbose output")
    parser.epilog='Example: phenix.python afitt_fd.py vAla3.pdb vAla3.cif NME 37'
    args = parser.parse_args()
    ligand_names=args.ligand_names.split(',')
    run(args.pdb_file, args.cif_file, ligand_names, int(args.atom),
      args.scale, args.verbose)

