from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file
import mmtbx.monomer_library.server
from mmtbx import monomer_library


# from here: https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
# LINK
# Overview

# The LINK records specify connectivity between residues that is not implied by the primary structure. Connectivity is expressed in terms of the atom names. They also include the distance associated with the each linkage following the symmetry operations at the end of each record.

# This record supplements information given in CONECT records and is provided here for convenience in searching.

# Record Format

# COLUMNS         DATA TYPE      FIELD           DEFINITION
# ------------------------------------------------------------------------------------
#  1 -  6         Record name    "LINK  "
# 13 - 16         Atom           name1           Atom name.
# 17              Character      altLoc1         Alternate location indicator.
# 18 - 20         Residue name   resName1        Residue  name.
# 22              Character      chainID1        Chain identifier.
# 23 - 26         Integer        resSeq1         Residue sequence number.
# 27              AChar          iCode1          Insertion code.
# 43 - 46         Atom           name2           Atom name.
# 47              Character      altLoc2         Alternate location indicator.
# 48 - 50         Residue name   resName2        Residue name.
# 52              Character      chainID2        Chain identifier.
# 53 - 56         Integer        resSeq2         Residue sequence number.
# 57              AChar          iCode2          Insertion code.
# 60 - 65         SymOP          sym1            Symmetry operator atom 1.
# 67 - 72         SymOP          sym2            Symmetry operator atom 2.
# 74 - 78         Real(5.2)      Length          Link distance
# Details

# The atoms involved in bonds between HET groups or between a HET group and standard residue are listed.
# Inter-residue linkages not implied by the primary structure are listed (e.g., reduced peptide bond).
# Non-standard linkages between residues, e.g., side-chain to side-chain, are listed.
# Each LINK record specifies one linkage.
# These records do not specify connectivity within a HET group (see CONECT) or disulfide bridges (see SSBOND).
# sym1 and sym2 are right justified and are given as blank when the identity operator (and no cell translation) is to be applied to the atom.
# - For NMR entries, only one set (or model) of LINK records will be supplied.
# - Coordinate bonds are also listed as LINKs
# Verification/Validation/Value Authority Control

# The distance between the pair of atoms listed must be consistent with the bonding.

# Relationships to Other Record Types

# CONECT records are generated from LINKs when both atoms are present in the entry. If symmetry operators are given to generate one of the residues involved in the bond, REMARK 290 defines the symmetry transformation.

# Example

#          1         2         3         4         5         6         7         8
# 12345678901234567890123456789012345678901234567890123456789012345678901234567890
# LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98
# LINK         OG1 THR A  51                NA    NA A6001     1555   1555  2.72
# LINK         OD2 ASP A  66                NA    NA A6001     1555   1555  2.72
# LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93

# LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93
# LINK         C21 2EG A   7                 C22 2EG B  19     1555   1555  1.56


restr_cif_txt = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
MCM        MCM '7-amino-4-methyl-2H-chromen-2-one' ligand 22 13 .
#
data_comp_MCM
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
MCM         N      N   NH2   .         -3.8313    0.0126   -0.1203
MCM         CA     C   CR6   .         -2.4124    0.0126   -0.1203
MCM         C2     C   CR16  .         -1.7183    0.0126    1.0827
MCM         C3     C   CR66  .         -0.3102    0.0126    1.0833
MCM         C4     C   CR6   .          1.7583    0.0126    2.2891
MCM         C5     C   CR16  .          2.4557    0.0126    1.0847
MCM         C6     C   CR6   .          1.7703    0.0126   -0.1062
MCM         C7     C   CR66  .          0.3770    0.0126   -0.1073
MCM         C8     C   CR16  .         -0.3279    0.0126   -1.3263
MCM         C9     C   CR16  .         -1.7167    0.0126   -1.3257
MCM         C10    C   CH3   .          2.5277   -0.0477   -1.4146
MCM         O1     O   O     .          0.3860    0.0126    2.2882
MCM         O2     O   OC    .          2.3647    0.0146    3.3417
MCM         HN1    H   HNH2  .         -4.2873    0.8024   -0.1159
MCM         HN2    H   HNH2  .         -4.2873   -0.7772   -0.1160
MCM         H2     H   HCR6  .         -2.2066    0.0126    1.9279
MCM         H5     H   HCR6  .          3.4319    0.0126    1.0860
MCM         H8     H   HCR6  .          0.1597    0.0126   -2.1720
MCM         H9     H   HCR6  .         -2.2050    0.0131   -2.1709
MCM        H101    H   HCH3  .          2.3139    0.7394   -1.9509
MCM        H102    H   HCH3  .          3.4869   -0.0680   -1.2348
MCM        H103    H   HCH3  .          2.2709   -0.8533   -1.9024
"""

pdb_1c_id = """
HETATM 2888  OH  ALY C   5     -15.286  24.739   8.076  1.00 45.73      C    O
HETATM 2889  CH  ALY C   5     -16.455  25.116   7.877  1.00 41.87      C    C
HETATM 2890  CH3 ALY C   5     -16.828  25.980   6.695  1.00 36.75      C    C
HETATM 2891  NZ  ALY C   5     -17.470  24.749   8.740  1.00 48.84      C    N
HETATM 2892  CE  ALY C   5     -17.263  23.930   9.901  1.00 39.21      C    C
HETATM 2893  CD  ALY C   5     -18.562  23.629  10.640  1.00 41.37      C    C
HETATM 2894  CG  ALY C   5     -18.329  22.814  11.892  1.00 41.99      C    C
HETATM 2895  CB  ALY C   5     -19.662  22.656  12.590  1.00 45.24      C    C
HETATM 2896  CA  ALY C   5     -19.524  22.011  13.974  1.00 51.20      C    C
HETATM 2897  N   ALY C   5     -20.755  21.670  14.665  1.00 49.87      C    N
HETATM 2898  C   ALY C   5     -18.619  22.936  14.781  1.00 47.11      C    C
HETATM 2899  O   ALY C   5     -17.427  22.728  15.024  1.00 49.30      C    O
HETATM 2900  N   MCM C   6     -19.324  23.360  15.960  1.00 49.19      C    N
HETATM 2901  CA  MCM C   6     -18.698  24.321  16.695  1.00 44.46      C    C
HETATM 2902  C2  MCM C   6     -19.284  25.579  16.752  1.00 49.13      C    C
HETATM 2903  C3  MCM C   6     -18.708  26.580  17.525  1.00 48.36      C    C
HETATM 2904  C4  MCM C   6     -18.700  28.830  18.338  1.00 45.36      C    C
HETATM 2905  C5  MCM C   6     -17.544  28.627  19.099  1.00 42.82      C    C
HETATM 2906  C6  MCM C   6     -16.902  27.390  19.114  1.00 45.52      C    C
HETATM 2907  C7  MCM C   6     -17.484  26.294  18.302  1.00 46.42      C    C
HETATM 2908  C8  MCM C   6     -16.929  25.030  18.246  1.00 40.92      C    C
HETATM 2909  C9  MCM C   6     -17.552  24.060  17.457  1.00 43.70      C    C
HETATM 2910  C10 MCM C   6     -15.654  27.180  19.937  1.00 42.32      C    C
HETATM 2911  O1  MCM C   6     -19.260  27.832  17.578  1.00 46.61      C    O
HETATM 2912  O2  MCM C   6     -19.248  29.950  18.344  1.00 44.55      C    O
"""

pdb_2c_id = """
HETATM 2888  OH  ALYaC   5     -15.286  24.739   8.076  1.00 45.73      C    O
HETATM 2889  CH  ALYaC   5     -16.455  25.116   7.877  1.00 41.87      C    C
HETATM 2890  CH3 ALYaC   5     -16.828  25.980   6.695  1.00 36.75      C    C
HETATM 2891  NZ  ALYaC   5     -17.470  24.749   8.740  1.00 48.84      C    N
HETATM 2892  CE  ALYaC   5     -17.263  23.930   9.901  1.00 39.21      C    C
HETATM 2893  CD  ALYaC   5     -18.562  23.629  10.640  1.00 41.37      C    C
HETATM 2894  CG  ALYaC   5     -18.329  22.814  11.892  1.00 41.99      C    C
HETATM 2895  CB  ALYaC   5     -19.662  22.656  12.590  1.00 45.24      C    C
HETATM 2896  CA  ALYaC   5     -19.524  22.011  13.974  1.00 51.20      C    C
HETATM 2897  N   ALYaC   5     -20.755  21.670  14.665  1.00 49.87      C    N
HETATM 2898  C   ALYaC   5     -18.619  22.936  14.781  1.00 47.11      C    C
HETATM 2899  O   ALYaC   5     -17.427  22.728  15.024  1.00 49.30      C    O
HETATM 2900  N   MCMaC   6     -19.324  23.360  15.960  1.00 49.19      C    N
HETATM 2901  CA  MCMaC   6     -18.698  24.321  16.695  1.00 44.46      C    C
HETATM 2902  C2  MCMaC   6     -19.284  25.579  16.752  1.00 49.13      C    C
HETATM 2903  C3  MCMaC   6     -18.708  26.580  17.525  1.00 48.36      C    C
HETATM 2904  C4  MCMaC   6     -18.700  28.830  18.338  1.00 45.36      C    C
HETATM 2905  C5  MCMaC   6     -17.544  28.627  19.099  1.00 42.82      C    C
HETATM 2906  C6  MCMaC   6     -16.902  27.390  19.114  1.00 45.52      C    C
HETATM 2907  C7  MCMaC   6     -17.484  26.294  18.302  1.00 46.42      C    C
HETATM 2908  C8  MCMaC   6     -16.929  25.030  18.246  1.00 40.92      C    C
HETATM 2909  C9  MCMaC   6     -17.552  24.060  17.457  1.00 43.70      C    C
HETATM 2910  C10 MCMaC   6     -15.654  27.180  19.937  1.00 42.32      C    C
HETATM 2911  O1  MCMaC   6     -19.260  27.832  17.578  1.00 46.61      C    O
HETATM 2912  O2  MCMaC   6     -19.248  29.950  18.344  1.00 44.55      C    O
"""
def exercise_1(prefix="tst_pdb_link_records_1"):
  with open('%s_mcm.cif' % prefix, 'w') as f:
    f.write(restr_cif_txt)
  with open('%s.pdb' % prefix, 'w') as f:
    f.write(pdb_1c_id)
  cmd = "phenix.geometry_minimization %s.pdb %s_mcm.cif write_geo=False" % (prefix, prefix)
  assert not easy_run.call(cmd)
  assert_lines_in_file('%s_minimized.pdb' % prefix,
      "LINK         C   ALY C   5                 N   MCM C   6 ",
      remove_white_spaces=False)
#LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98
#LINK         C   ALY C   5                 N   MCM C   6
#LINK         C   ALYaC   5                 N   MCMaC   6
def exercise_2(prefix="tst_pdb_link_records_2"):
  with open('%s_mcm.cif' % prefix, 'w') as f:
    f.write(restr_cif_txt)
  with open('%s.pdb' % prefix, 'w') as f:
    f.write(pdb_2c_id)
  cmd = "phenix.geometry_minimization %s.pdb %s_mcm.cif write_geo=False" % (prefix, prefix)
  assert not easy_run.call(cmd)
  assert_lines_in_file('%s_minimized.pdb' % prefix,
      "LINK         C   ALYaC   5                 N   MCMaC   6 ",
      remove_white_spaces=False)

if __name__ == '__main__':
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_1()
    exercise_2()
    print("OK")
