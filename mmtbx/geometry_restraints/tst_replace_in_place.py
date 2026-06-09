from __future__ import division
from libtbx import easy_run

pdb_lines = '''
ATOM  11248  N   THR T  12     -68.780 -34.422  27.326  1.00 32.02           N
ATOM  11249  CA  THR T  12     -68.470 -33.115  27.894  1.00 33.26           C
ATOM  11250  C   THR T  12     -66.974 -32.798  27.988  1.00 33.30           C
ATOM  11251  O   THR T  12     -66.556 -31.730  27.571  1.00 42.35           O
ATOM  11252  CB  THR T  12     -69.119 -32.902  29.279  1.00 33.46           C
ATOM  11253  OG1 THR T  12     -68.685 -33.924  30.174  1.00 51.63           O
ATOM  11254  CG2 THR T  12     -70.618 -32.958  29.153  1.00 34.58           C
HETATM11255  N   FVU T  13     -66.196 -33.776  28.548  1.00 32.54           N
HETATM11256  CA AFVU T  13     -64.701 -33.455  28.454  0.50 31.52           C
HETATM11257  C   FVU T  13     -64.006 -34.747  28.463  1.00 22.94           C
HETATM11258  O   FVU T  13     -64.527 -35.765  28.860  1.00 21.83           O
HETATM11259  CB AFVU T  13     -64.210 -32.523  29.561  0.50 26.71           C
HETATM11260  C7 AFVU T  13     -64.157 -33.164  30.947  0.50 41.04           C
HETATM11261  C8 AFVU T  13     -63.580 -32.238  32.015  0.50 31.57           C
HETATM11262  C9 AFVU T  13     -63.648 -32.829  33.404  0.50 46.08           C
HETATM11263  N11AFVU T  13     -62.678 -32.479  34.239  0.50 33.37           N
HETATM11264  O10AFVU T  13     -64.561 -33.598  33.721  0.50 39.62           O
HETATM11265  CA BFVU T  13     -64.650 -33.480  28.449  0.50 31.52           C
HETATM11266  CB BFVU T  13     -64.157 -32.542  29.550  0.50 26.71           C
HETATM11267  C7 BFVU T  13     -64.023 -31.192  29.722  0.50 41.04           C
HETATM11268  C8 BFVU T  13     -63.336 -30.544  30.557  0.50 31.57           C
HETATM11269  C9 BFVU T  13     -62.571 -29.525  30.528  0.50 46.08           C
HETATM11270  N11BFVU T  13     -62.543 -28.357  31.548  0.50 33.37           N
HETATM11271  O10BFVU T  13     -62.108 -29.471  29.614  0.50 39.62           O
HETATM11272  N   A7U T  14     -60.604 -35.202  26.006  1.00 24.23           N
HETATM11273  CA  A7U T  14     -58.632 -36.413  25.079  1.00 19.95           C
HETATM11274  C   A7U T  14     -58.856 -37.861  24.744  1.00 31.05           C
HETATM11275  O   A7U T  14     -59.920 -38.241  24.303  1.00 25.22           O
HETATM11276  CB  A7U T  14     -59.853 -35.553  24.792  1.00 21.25           C
HETATM11277  C10 A7U T  14     -61.972 -35.956  27.900  1.00 22.51           C
HETATM11278  C12 A7U T  14     -60.995 -35.703  29.047  1.00 24.94           C
HETATM11279  C7  A7U T  14     -62.257 -33.705  27.010  1.00 22.98           C
HETATM11280  C8  A7U T  14     -61.599 -34.166  25.738  1.00 19.62           C
HETATM11281  C9  A7U T  14     -61.261 -36.355  26.624  1.00 20.11           C
HETATM11282  N11 A7U T  14     -62.882 -34.813  27.747  1.00 24.82           N
HETATM11283  N14 A7U T  14     -60.344 -36.757  29.510  1.00 18.68           N
HETATM11284  O13 A7U T  14     -60.859 -34.566  29.508  1.00 21.49           O'''

data_link_cif_lines = '''
data_link_pep_link_12_THR_13_FVU

loop_
_chem_link_bond.link_id
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
 pep_link_12_THR_13_FVU 1 C      2 N      single   1.338 0.010

loop_
_chem_link_angle.link_id
_chem_link_angle.atom_1_comp_id
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
 pep_link_12_THR_13_FVU 1 CA     1 C      2 N      116.30 1.50
 pep_link_12_THR_13_FVU 1 O      1 C      2 N      123.20 1.10
 pep_link_12_THR_13_FVU 1 C      2 N      2 CA     121.60 1.70

loop_
_chem_link_plane.link_id
_chem_link_plane.plane_id
_chem_link_plane.atom_comp_id
_chem_link_plane.atom_id
_chem_link_plane.dist_esd
 pep_link_12_THR_13_FVU plan1    1 CA     0.020
 pep_link_12_THR_13_FVU plan1    1 C      0.020
 pep_link_12_THR_13_FVU plan1    1 O      0.020
 pep_link_12_THR_13_FVU plan1    2 N      0.020
 pep_link_12_THR_13_FVU plan2    1 C      0.020
 pep_link_12_THR_13_FVU plan2    2 N      0.020
 pep_link_12_THR_13_FVU plan2    2 CA     0.020

data_link_cross_link_13_FVU_14_A7U

loop_
_chem_link_bond.link_id
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
 cross_link_13_FVU_14_A7U 1 C      2 N11    single   1.343 0.010

loop_
_chem_link_angle.link_id
_chem_link_angle.atom_1_comp_id
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
 cross_link_13_FVU_14_A7U 1 CA     1 C      2 N11    118.60 1.90
 cross_link_13_FVU_14_A7U 1 O      1 C      2 N11    121.50 0.90
 cross_link_13_FVU_14_A7U 1 C      2 N11    2 C10    121.30 3.00
 cross_link_13_FVU_14_A7U 1 C      2 N11    2 C7     123.70 3.90

loop_
_chem_link_plane.link_id
_chem_link_plane.plane_id
_chem_link_plane.atom_comp_id
_chem_link_plane.atom_id
_chem_link_plane.dist_esd
 cross_link_13_FVU_14_A7U plan1    1 CA     0.020
 cross_link_13_FVU_14_A7U plan1    1 C      0.020
 cross_link_13_FVU_14_A7U plan1    1 O      0.020
'''

apply_link_lines = '''
refinement.pdb_interpretation.apply_cif_link {
  data_link = pep_link_12_THR_13_FVU
  residue_selection_1 = chain T and resname THR and resseq 12
  residue_selection_2 = chain T and resname FVU and resseq 13 and altloc A
}

refinement.pdb_interpretation.apply_cif_link {
  data_link = pep_link_12_THR_13_FVU
  residue_selection_1 = chain T and resname THR and resseq 12
  residue_selection_2 = chain T and resname FVU and resseq 13 and altloc B
}

refinement.pdb_interpretation.apply_cif_link {
  data_link = cross_link_13_FVU_14_A7U
  residue_selection_1 = chain T and resname FVU and resseq 13 and altloc A
  residue_selection_2 = chain T and resname A7U and resseq 14
}

refinement.pdb_interpretation.apply_cif_link {
  data_link = cross_link_13_FVU_14_A7U
  residue_selection_1 = chain T and resname FVU and resseq 13 and altloc B
  residue_selection_2 = chain T and resname A7U and resseq 14
}
'''

restraints_one = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
FVU        FVU 'Unknown                  ' ligand 22 10 .
#
data_comp_FVU
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
FVU         N      N   NH2    0    .     -66.1960  -33.7760   28.5480
FVU         C      C   C1     0    .     -64.0060  -34.7470   28.4630
FVU         O      O   O      0    .     -64.5270  -35.7650   28.8600
FVU         CA     C   CH1    0    .     -64.7010  -33.4550   28.4540
FVU         CB     C   CH2    0    .     -64.2100  -32.5230   29.5610
FVU         C7     C   CH2    0    .     -64.1570  -33.1640   30.9470
FVU         C8     C   CH2    0    .     -63.5800  -32.2380   32.0150
FVU         C9     C   C      0    .     -63.6480  -32.8290   33.4040
FVU         N11    N   NH2    0    .     -62.6780  -32.4790   34.2390
FVU         O10    O   O      0    .     -64.5610  -33.5980   33.7210
FVU         H      H   HNH2   0    .     -66.7286  -32.9313   28.4277
FVU         H2     H   HNH2   0    .     -66.4436  -34.4338   27.8285
FVU         HC1    H   H      0    .     -62.9822  -34.7874   28.0962
FVU         HA     H   HCH1   0    .     -64.5144  -32.9725   27.4949
FVU         HB2    H   HCH2   0    .     -64.8739  -31.6602   29.6052
FVU         HB3    H   HCH2   0    .     -63.2095  -32.1772   29.3023
FVU         H71    H   HCH2   0    .     -63.5424  -34.0615   30.8937
FVU         H72    H   HCH2   0    .     -65.1671  -33.4482   31.2395
FVU         H81    H   HCH2   0    .     -64.1365  -31.3016   32.0049
FVU         H82    H   HCH2   0    .     -62.5387  -32.0302   31.7732
FVU        H111    H   HNH2   0    .     -62.6668  -32.8415   35.1768
FVU        H112    H   HNH2   0    .     -61.9531  -31.8536   33.9333
#
'''

restraints_two = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
A7U        A7U 'Unknown                  ' ligand 30 13 .
#
data_comp_A7U
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
A7U         N      N   NT     0    .     -60.6040  -35.2020   26.0060
A7U         CA     C   CH2    0    .     -58.6320  -36.4130   25.0790
A7U         C      C   CH2    0    .     -58.8560  -37.8610   24.7440
A7U         O      O   OH1    0    .     -59.9200  -38.2410   24.3030
A7U         CB     C   CH2    0    .     -59.8530  -35.5530   24.7920
A7U         C10    C   CH1    0    .     -61.9720  -35.9560   27.9000
A7U         C12    C   C      0    .     -60.9950  -35.7030   29.0470
A7U         C7     C   CH2    0    .     -62.2570  -33.7050   27.0100
A7U         C8     C   CH2    0    .     -61.5990  -34.1660   25.7380
A7U         C9     C   CH2    0    .     -61.2610  -36.3550   26.6240
A7U         N11    N   NH1    0    .     -62.8820  -34.8130   27.7470
A7U         N14    N   NH2    0    .     -60.3440  -36.7570   29.5100
A7U         O13    O   O      0    .     -60.8590  -34.5660   29.5080
A7U         HA     H   HCH2   0    .     -57.7933  -36.0412   24.4905
A7U         HA3    H   HCH2   0    .     -58.3844  -36.3308   26.1358
A7U         HC2    H   HCH2   0    .     -58.6669  -38.4309   25.6534
A7U         HC3    H   HCH2   0    .     -58.0921  -38.1420   24.0187
A7U         HO1    H   HOH1   0    .     -59.8343  -39.1340   23.9962
A7U         HB2    H   HCH2   0    .     -59.5278  -34.6346   24.3031
A7U         HB3    H   HCH2   0    .     -60.5127  -36.0955   24.1164
A7U        H101    H   HCH1   0    .     -62.5890  -36.8073   28.1896
A7U         H71    H   HCH2   0    .     -61.5054  -33.2395   27.6459
A7U         H72    H   HCH2   0    .     -63.0233  -32.9699   26.7634
A7U         H81    H   HCH2   0    .     -62.3599  -34.5676   25.0689
A7U         H82    H   HCH2   0    .     -61.1093  -33.3163   25.2621
A7U         H91    H   HCH2   0    .     -60.5120  -37.1096   26.8553
A7U         H92    H   HCH2   0    .     -61.9887  -36.7680   25.9246
A7U        H111    H   HNH1   0    .     -63.1624  -34.5028   28.6212
A7U        H141    H   HNH2   0    .     -60.5769  -37.6780   29.1779
A7U        H142    H   HNH2   0    .     -59.6373  -36.6411   30.2201
#
'''

def main():
  preamble='tst_replace_in_place'
  f=open(f'{preamble}.pdb', 'w')
  f.write(pdb_lines)
  del f
  f=open(f'{preamble}_data_link.cif', 'w')
  f.write(data_link_cif_lines)
  del f
  f=open(f'{preamble}_apply_link.def', 'w')
  f.write(apply_link_lines)
  del f
  f=open(f'{preamble}_one.cif', 'w')
  f.write(restraints_one)
  del f
  f=open(f'{preamble}_two.cif', 'w')
  f.write(restraints_two)
  del f
  cmd = f'phenix.pdb_interpretation {preamble}.pdb'
  cmd+= f' {preamble}_data_link.cif {preamble}_apply_link.def'
  cmd+= f' {preamble}_one.cif {preamble}_two.cif'
  rc=easy_run.go(cmd)
  print(cmd)
  assert rc.stdout_lines[-2].find('scattering factor at diffraction angle 0')>-1

if __name__ == '__main__':
  main()
