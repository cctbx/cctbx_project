from __future__ import absolute_import, division, print_function
'''
Tests for mmtbx.ligands.rdkit_utils.residue_molecule: the RDKit molecule of one
residue conformer of a processed model (bond orders and formal charges from
rdkit.Chem.rdDetermineBonds.DetermineBondOrders on the restraints' graph with the
total charge from the restraint file).

Restraint files inline (coordinates from RDKit embedding; nonbonded types by
element): ZAC acetate ('deloc' C-O, -1 on O2), ZZW taurine zwitterion ('deloc'
S-O, +1 on N1, -1 on O3), ZTP methyl triphosphate (total -4), ZPC acetate with
partial charges only and 'coval' bonds, ZQM the same with a charge column of
only '?', ZSH ethanethiol (thiol H), ZNC acetate with -1 on both O (inconsistent),
ZAH acetate with -1 on O2 and an H (HO21) on O2 in the file, ZGU methylguanidinium
with partial charges only (sum 1.548, rounds to 2), ZC5 a carbon with five H (no
valid structure at any total), ZNM nitromethane with partial charges only (sum
0; -2 valid for DetermineBondOrders, but its N1-O1 single contradicts the file);
derived in the exercises: ZAA acetic acid (from ZAH), ZZS ZSH with its Zn in the
residue, ZAC with an explicit C=C.
OFO (Fe-O-Fe-OH) inline from the CCD entry (ideal coordinates; Fe-O 'metal').
Library restraint files: ACT, GOL (GeoStd coordinates), GLY and SEP (chain), SER
(ester link), HEM (CCD ideal coordinates, no H).
'''
import os
import sys
import tempfile
import iotbx.cif
import iotbx.pdb
import iotbx.phil
import mmtbx.model
from libtbx.utils import null_out
from rdkit import Chem
from mmtbx.ligands import rdkit_utils

zac_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZAC  ZAC  'ZAC' ligand 7 4 .

data_comp_ZAC
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
 ZAC  C2   C  C     0   0.000  -0.6392  -0.0414   0.0138
 ZAC  C1   C  C     0   0.000   0.8764   0.0503  -0.0201
 ZAC  O1   O  O     0   0.000   1.4777  -1.0385  -0.2345
 ZAC  O2   O  O    -1   0.000   1.3412   1.2089   0.1715
 ZAC  HC21 H  H     0   0.000  -1.0086   0.2872   0.9898
 ZAC  HC22 H  H     0   0.000  -0.9807  -1.0667  -0.1580
 ZAC  HC23 H  H     0   0.000  -1.0668   0.6002  -0.7625

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZAC  C2   C1   single   1.519  0.020
 ZAC  C1   O1   deloc    1.262  0.020
 ZAC  C1   O2   deloc    1.263  0.020
 ZAC  C2   HC21 single   1.094  0.020
 ZAC  C2   HC22 single   1.094  0.020
 ZAC  C2   HC23 single   1.094  0.020
'''

zzw_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZZW  ZZW  'ZZW' ligand 14 7 .

data_comp_ZZW
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
 ZZW  N1   N  N     1   0.000  -1.5027   0.1813   0.4595
 ZZW  C1   C  C     0   0.000  -0.7130  -0.4747  -0.6253
 ZZW  C2   C  C     0   0.000   0.6175   0.2167  -0.8372
 ZZW  S1   S  S     0   0.000   1.5968   0.0032   0.6133
 ZZW  O1   O  O     0   0.000   1.8021  -1.4372   0.6304
 ZZW  O2   O  O     0   0.000   0.6590   0.4920   1.6344
 ZZW  O3   O  O    -1   0.000   2.7547   0.8408   0.3596
 ZZW  HN11 H  H     0   0.000  -2.2777  -0.3868   0.8070
 ZZW  HN12 H  H     0   0.000  -1.8052   1.1329   0.2383
 ZZW  HN13 H  H     0   0.000  -0.8405   0.3016   1.2769
 ZZW  HC11 H  H     0   0.000  -0.6067  -1.5254  -0.3391
 ZZW  HC12 H  H     0   0.000  -1.3323  -0.4180  -1.5251
 ZZW  HC21 H  H     0   0.000   0.4874   1.2898  -1.0104
 ZZW  HC22 H  H     0   0.000   1.1605  -0.2162  -1.6824

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZZW  N1   C1   single   1.494  0.020
 ZZW  C1   C2   single   1.514  0.020
 ZZW  C2   S1   single   1.763  0.020
 ZZW  S1   O1   deloc    1.455  0.020
 ZZW  S1   O2   deloc    1.470  0.020
 ZZW  S1   O3   deloc    1.451  0.020
 ZZW  N1   HN11 single   1.022  0.020
 ZZW  N1   HN12 single   1.023  0.020
 ZZW  N1   HN13 single   1.059  0.020
 ZZW  C1   HC11 single   1.094  0.020
 ZZW  C1   HC12 single   1.094  0.020
 ZZW  C2   HC21 single   1.095  0.020
 ZZW  C2   HC22 single   1.094  0.020
'''

ztp_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZTP  ZTP  'ZTP' ligand 17 14 .

data_comp_ZTP
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
 ZTP  C1   C  C     0   0.000   2.1660   0.7983   0.1225
 ZTP  O1   O  O     0   0.000   2.4232  -0.5322  -0.2750
 ZTP  P1   P  P     0   0.000   1.5478  -1.0917  -1.5320
 ZTP  O11  O  O     0   0.000   1.9457  -0.2167  -2.7224
 ZTP  O12  O  O    -1   0.000   1.8414  -2.5754  -1.6882
 ZTP  O2   O  O     0   0.000   0.0429  -0.7555  -1.1116
 ZTP  P2   P  P     0   0.000  -0.9076  -1.1054   0.1538
 ZTP  O21  O  O     0   0.000  -0.0307  -1.2273   1.3919
 ZTP  O22  O  O    -1   0.000  -1.7498  -2.3090  -0.2678
 ZTP  O3   O  O     0   0.000  -1.7918   0.2423   0.1382
 ZTP  P3   P  P     0   0.000  -2.7538   1.1195   1.0957
 ZTP  O31  O  O     0   0.000  -3.5193   2.0399   0.1272
 ZTP  O32  O  O    -1   0.000  -1.8528   1.9457   2.0181
 ZTP  O33  O  O    -1   0.000  -3.7069   0.1869   1.8421
 ZTP  HC11 H  H     0   0.000   1.1037   0.9758   0.3098
 ZTP  HC12 H  H     0   0.000   2.5221   1.4913  -0.6487
 ZTP  HC13 H  H     0   0.000   2.7199   0.9961   1.0464

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZTP  C1   O1   single   1.412  0.020
 ZTP  O1   P1   single   1.631  0.020
 ZTP  P1   O11  deloc    1.530  0.020
 ZTP  P1   O12  deloc    1.521  0.020
 ZTP  P1   O2   single   1.598  0.020
 ZTP  O2   P2   single   1.621  0.020
 ZTP  P2   O21  deloc    1.522  0.020
 ZTP  P2   O22  deloc    1.528  0.020
 ZTP  P2   O3   single   1.612  0.020
 ZTP  O3   P3   single   1.616  0.020
 ZTP  P3   O31  deloc    1.540  0.020
 ZTP  P3   O32  deloc    1.531  0.020
 ZTP  P3   O33  deloc    1.528  0.020
 ZTP  C1   HC11 single   1.093  0.020
 ZTP  C1   HC12 single   1.096  0.020
 ZTP  C1   HC13 single   1.095  0.020
'''

zpc_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZPC  ZPC  'ZPC' ligand 7 4 .

data_comp_ZPC
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZPC  C2   C  C   -0.300  -0.6392  -0.0414   0.0138
 ZPC  C1   C  C    0.400   0.8764   0.0503  -0.0201
 ZPC  O1   O  O   -0.700   1.4777  -1.0385  -0.2345
 ZPC  O2   O  O   -0.700   1.3412   1.2089   0.1715
 ZPC  HC21 H  H    0.100  -1.0086   0.2872   0.9898
 ZPC  HC22 H  H    0.100  -0.9807  -1.0667  -0.1580
 ZPC  HC23 H  H    0.100  -1.0668   0.6002  -0.7625

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZPC  C2   C1   coval    1.519  0.020
 ZPC  C1   O1   coval    1.262  0.020
 ZPC  C1   O2   coval    1.263  0.020
 ZPC  C2   HC21 coval    1.094  0.020
 ZPC  C2   HC22 coval    1.094  0.020
 ZPC  C2   HC23 coval    1.094  0.020
'''

zsh_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZSH  ZSH  'ZSH' ligand 9 3 .

data_comp_ZSH
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
 ZSH  C2   C  C     0   0.000  -0.8750   0.3315  -0.0739
 ZSH  C1   C  C     0   0.000   0.2659  -0.6636   0.0255
 ZSH  S1   S  S     0   0.000   1.7093   0.0923   0.8229
 ZSH  HC21 H  H     0   0.000  -0.5973   1.2033  -0.6760
 ZSH  HC22 H  H     0   0.000  -1.7394  -0.1397  -0.5534
 ZSH  HC23 H  H     0   0.000  -1.1934   0.6782   0.9151
 ZSH  HC11 H  H     0   0.000   0.5547  -1.0098  -0.9714
 ZSH  HC12 H  H     0   0.000  -0.0429  -1.5341   0.6120
 ZSH  HS11 H  H     0   0.000   1.9180   1.0419  -0.1009

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZSH  C2   C1   single   1.517  0.020
 ZSH  C1   S1   single   1.814  0.020
 ZSH  C2   HC21 single   1.095  0.020
 ZSH  C2   HC22 single   1.095  0.020
 ZSH  C2   HC23 single   1.095  0.020
 ZSH  C1   HC11 single   1.094  0.020
 ZSH  C1   HC12 single   1.094  0.020
 ZSH  S1   HS11 single   1.341  0.020
'''

zah_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZAH  ZAH  'ZAH' ligand 8 4 .

data_comp_ZAH
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
 ZAH  C2   C  C     0   0.000  -0.9565  -0.0956  -0.0643
 ZAH  C1   C  C     0   0.000   0.4861   0.2884  -0.0648
 ZAH  O1   O  O     0   0.000   0.9554   1.3471  -0.4445
 ZAH  O2   O  O    -1   0.000   1.2772  -0.6877   0.4173
 ZAH  HC21 H  H     0   0.000  -1.5507   0.7309  -0.4644
 ZAH  HC22 H  H     0   0.000  -1.2854  -0.2993   0.9577
 ZAH  HC23 H  H     0   0.000  -1.1067  -0.9726  -0.6988
 ZAH  HO21 H  H     0   0.000   2.1805  -0.3112   0.3617

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZAH  C2   C1   single   1.493  0.020
 ZAH  C1   O1   deloc    1.219  0.020
 ZAH  C1   O2   deloc    1.346  0.020
 ZAH  C2   HC21 single   1.094  0.020
 ZAH  C2   HC22 single   1.093  0.020
 ZAH  C2   HC23 single   1.093  0.020
 ZAH  O2   HO21 single   0.980  0.020
'''

zgu_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZGU  ZGU  'ZGU' ligand 13 5 .

data_comp_ZGU
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZGU  C1   C  C    0.000   1.7918  -0.2410  -0.1613
 ZGU  N1   N  N   -0.200   0.4158  -0.7338  -0.1517
 ZGU  C2   C  C    0.500  -0.6648   0.0343   0.0506
 ZGU  N2   N  N   -0.300  -0.5665   1.3482   0.2653
 ZGU  N3   N  N   -0.300  -1.8793  -0.5233   0.0384
 ZGU  HC11 H  H    0.231   2.0363   0.1874   0.8138
 ZGU  HC12 H  H    0.231   1.9163   0.4873  -0.9663
 ZGU  HC13 H  H    0.231   2.4579  -1.0877  -0.3489
 ZGU  HN11 H  H    0.231   0.2668  -1.7229  -0.3083
 ZGU  HN21 H  H    0.231   0.3325   1.8112   0.2828
 ZGU  HN22 H  H    0.231  -1.3824   1.9260   0.4177
 ZGU  HN31 H  H    0.231  -2.0100  -1.5142  -0.1198
 ZGU  HN32 H  H    0.231  -2.7144   0.0286   0.1877

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZGU  C1   N1   single   1.462  0.020
 ZGU  N1   C2   single   1.341  0.020
 ZGU  C2   N2   double   1.335  0.020
 ZGU  C2   N3   single   1.336  0.020
 ZGU  C1   HC11 single   1.093  0.020
 ZGU  C1   HC12 single   1.093  0.020
 ZGU  C1   HC13 single   1.094  0.020
 ZGU  N1   HN11 single   1.012  0.020
 ZGU  N2   HN21 single   1.011  0.020
 ZGU  N2   HN22 single   1.011  0.020
 ZGU  N3   HN31 single   1.012  0.020
 ZGU  N3   HN32 single   1.012  0.020
'''

znm_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZNM  ZNM  'ZNM' ligand 7 4 .

data_comp_ZNM
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZNM  C1   C  C   -0.300  -0.6519  -0.0458   0.0124
 ZNM  N1   N  N    0.700   0.8315   0.0577  -0.0234
 ZNM  O1   O  O   -0.350   1.4709  -0.9962   0.0763
 ZNM  O2   O  O   -0.350   1.3133   1.1917  -0.1300
 ZNM  HC11 H  H    0.100  -0.9554   0.0310   1.0581
 ZNM  HC12 H  H    0.100  -0.9400  -1.0097  -0.4127
 ZNM  HC13 H  H    0.100  -1.0683   0.7713  -0.5807

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZNM  C1   N1   single   1.487  0.020
 ZNM  N1   O1   double   1.237  0.020
 ZNM  N1   O2   single   1.237  0.020
 ZNM  C1   HC11 single   1.092  0.020
 ZNM  C1   HC12 single   1.092  0.020
 ZNM  C1   HC13 single   1.092  0.020
'''

zc5_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZC5  ZC5  'ZC5' ligand 6 1 .

data_comp_ZC5
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
 ZC5  C1   C  C     0   0.000   0.0000   0.0000   0.0000
 ZC5  H1   H  H     0   0.000   1.0900   0.0000   0.0000
 ZC5  H2   H  H     0   0.000  -1.0900   0.0000   0.0000
 ZC5  H3   H  H     0   0.000   0.0000   1.0900   0.0000
 ZC5  H4   H  H     0   0.000   0.0000  -1.0900   0.0000
 ZC5  H5   H  H     0   0.000   0.0000   0.0000   1.0900

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZC5  C1   H1   single   1.090  0.020
 ZC5  C1   H2   single   1.090  0.020
 ZC5  C1   H3   single   1.090  0.020
 ZC5  C1   H4   single   1.090  0.020
 ZC5  C1   H5   single   1.090  0.020
'''

znc_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZNC  ZNC  'ZNC' ligand 7 4 .

data_comp_ZNC
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
 ZNC  C2   C  C     0   0.000  -0.6392  -0.0414   0.0138
 ZNC  C1   C  C     0   0.000   0.8764   0.0503  -0.0201
 ZNC  O1   O  O    -1   0.000   1.4777  -1.0385  -0.2345
 ZNC  O2   O  O    -1   0.000   1.3412   1.2089   0.1715
 ZNC  HC21 H  H     0   0.000  -1.0086   0.2872   0.9898
 ZNC  HC22 H  H     0   0.000  -0.9807  -1.0667  -0.1580
 ZNC  HC23 H  H     0   0.000  -1.0668   0.6002  -0.7625

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZNC  C2   C1   single   1.519  0.020
 ZNC  C1   O1   deloc    1.262  0.020
 ZNC  C1   O2   deloc    1.263  0.020
 ZNC  C2   HC21 single   1.094  0.020
 ZNC  C2   HC22 single   1.094  0.020
 ZNC  C2   HC23 single   1.094  0.020
'''

ofo_cif = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 OFO  OFO  'OFO' ligand 5 4 .

data_comp_OFO
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
 OFO  FE1  FE FE    ?   0.000   1.9080   0.6530  -0.7540
 OFO  O    O  O    -2   0.000   0.8030  -0.5330   0.0040
 OFO  FE2  FE FE    ?   0.000  -0.2920   0.9040   0.8130
 OFO  OH   O  O    -1   0.000   0.3370   0.3450   2.6600
 OFO  HO   H  H     0   0.000  -0.3500   0.1930   3.1650

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 OFO  FE1  O    metal    1.900  0.050
 OFO  O    FE2  metal    1.900  0.050
 OFO  FE2  OH   metal    1.950  0.050
 OFO  OH   HO   single   0.960  0.020
'''

# GeoStd ACT (acetate; OXT -1) and GOL, at their dictionary coordinates
act_pdb = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C   ACT A   1      20.858  19.841  20.175  1.00 20.00           C
HETATM    2  O   ACT A   1      21.125  19.609  21.373  1.00 20.00           O
HETATM    3  OXT ACT A   1      21.670  19.900  19.230  1.00 20.00           O
HETATM    4  CH3 ACT A   1      19.375  20.104  19.853  1.00 20.00           C
HETATM    5  H1  ACT A   1      19.167  20.122  18.783  1.00 20.00           H
HETATM    6  H2  ACT A   1      18.733  19.354  20.319  1.00 20.00           H
HETATM    7  H3  ACT A   1      19.071  21.069  20.267  1.00 20.00           H
END
'''

gol_pdb = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  GOL A   1      19.128  19.265  20.755  1.00 20.00           C
HETATM    2  O1  GOL A   1      17.768  19.028  20.478  1.00 20.00           O
HETATM    3  C2  GOL A   1      19.737  20.066  19.624  1.00 20.00           C
HETATM    4  O2  GOL A   1      19.115  21.330  19.582  1.00 20.00           O
HETATM    5  C3  GOL A   1      21.239  20.197  19.824  1.00 20.00           C
HETATM    6  O3  GOL A   1      21.752  20.949  18.746  1.00 20.00           O
HETATM    7  H11 GOL A   1      19.263  19.814  21.697  1.00 20.00           H
HETATM    8  H12 GOL A   1      19.690  18.326  20.848  1.00 20.00           H
HETATM    9  HO1 GOL A   1      17.403  18.511  21.201  1.00 20.00           H
HETATM   10  H2  GOL A   1      19.567  19.516  18.684  1.00 20.00           H
HETATM   11  HO2 GOL A   1      19.581  21.862  18.928  1.00 20.00           H
HETATM   12  H31 GOL A   1      21.436  20.686  20.787  1.00 20.00           H
HETATM   13  H32 GOL A   1      21.687  19.197  19.868  1.00 20.00           H
HETATM   14  HO3 GOL A   1      22.634  21.254  18.977  1.00 20.00           H
END
'''

# HEM heavy atoms at the CCD ideal coordinates (monomer library names)
hem_pdb = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  CHA HEM A   1      22.451  20.827  20.221  1.00 20.00           C
HETATM    2  CHB HEM A   1      18.190  23.025  20.485  1.00 20.00           C
HETATM    3  CHC HEM A   1      15.972  18.774  20.043  1.00 20.00           C
HETATM    4  CHD HEM A   1      20.248  16.558  20.152  1.00 20.00           C
HETATM    5  C1A HEM A   1      21.505  21.846  20.321  1.00 20.00           C
HETATM    6  C2A HEM A   1      21.753  23.203  20.408  1.00 20.00           C
HETATM    7  C3A HEM A   1      20.543  23.826  20.482  1.00 20.00           C
HETATM    8  C4A HEM A   1      19.574  22.851  20.430  1.00 20.00           C
HETATM    9  CMA HEM A   1      20.337  25.316  20.587  1.00 20.00           C
HETATM   10  CAA HEM A   1      23.106  23.870  20.422  1.00 20.00           C
HETATM   11  CBA HEM A   1      23.658  24.190  19.036  1.00 20.00           C
HETATM   12  CGA HEM A   1      25.033  24.851  19.038  1.00 20.00           C
HETATM   13  O1A HEM A   1      26.038  24.112  19.090  1.00 20.00           O
HETATM   14  O2A HEM A   1      25.083  26.098  18.986  1.00 20.00           O
HETATM   15  C1B HEM A   1      17.164  22.091  20.335  1.00 20.00           C
HETATM   16  C2B HEM A   1      15.800  22.363  20.274  1.00 20.00           C
HETATM   17  C3B HEM A   1      15.119  21.141  20.109  1.00 20.00           C
HETATM   18  C4B HEM A   1      16.120  20.169  20.131  1.00 20.00           C
HETATM   19  CMB HEM A   1      15.147  23.720  20.299  1.00 20.00           C
HETATM   20  CAB HEM A   1      13.635  21.025  20.019  1.00 20.00           C
HETATM   21  CBB HEM A   1      12.835  20.079  19.586  1.00 20.00           C
HETATM   22  C1C HEM A   1      16.921  17.751  20.009  1.00 20.00           C
HETATM   23  C2C HEM A   1      16.664  16.387  19.886  1.00 20.00           C
HETATM   24  C3C HEM A   1      17.898  15.709  19.888  1.00 20.00           C
HETATM   25  C4C HEM A   1      18.856  16.711  20.045  1.00 20.00           C
HETATM   26  CMC HEM A   1      15.317  15.732  19.722  1.00 20.00           C
HETATM   27  CAC HEM A   1      18.035  14.227  19.788  1.00 20.00           C
HETATM   28  CBC HEM A   1      19.042  13.439  19.486  1.00 20.00           C
HETATM   29  C1D HEM A   1      21.274  17.504  20.145  1.00 20.00           C
HETATM   30  C2D HEM A   1      22.629  17.270  20.120  1.00 20.00           C
HETATM   31  C3D HEM A   1      23.256  18.479  20.132  1.00 20.00           C
HETATM   32  C4D HEM A   1      22.268  19.444  20.186  1.00 20.00           C
HETATM   33  CMD HEM A   1      23.335  15.939  20.068  1.00 20.00           C
HETATM   34  CAD HEM A   1      24.746  18.711  20.116  1.00 20.00           C
HETATM   35  CBD HEM A   1      25.371  18.845  21.502  1.00 20.00           C
HETATM   36  CGD HEM A   1      26.848  19.231  21.497  1.00 20.00           C
HETATM   37  O1D HEM A   1      27.135  20.446  21.449  1.00 20.00           O
HETATM   38  O2D HEM A   1      27.694  18.313  21.540  1.00 20.00           O
HETATM   39  NA  HEM A   1      20.161  21.624  20.337  1.00 20.00           N
HETATM   40  NB  HEM A   1      17.380  20.750  20.265  1.00 20.00           N
HETATM   41  NC  HEM A   1      18.259  17.968  20.123  1.00 20.00           N
HETATM   42  ND  HEM A   1      21.042  18.847  20.193  1.00 20.00           N
HETATM   43  FE  HEM A   1      19.209  19.793  20.272  1.00 20.00          FE
END
'''

# Gly-SEP-Gly (RDKit embedding): H on SEP only (H, HA, HB2, HB3, HOP2, HOP3);
# peptide bonds from pdb_interpretation
chain_pdb = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1      17.870  21.775  17.752  1.00 20.00           N
ATOM      2  CA  GLY A   1      18.413  22.821  18.643  1.00 20.00           C
ATOM      3  C   GLY A   1      19.403  22.190  19.625  1.00 20.00           C
ATOM      4  O   GLY A   1      20.441  22.738  19.988  1.00 20.00           O
HETATM    5  N   SEP A   2      19.028  20.929  20.039  1.00 20.00           N
HETATM    6  CA  SEP A   2      19.963  20.053  20.747  1.00 20.00           C
HETATM    7  CB  SEP A   2      19.219  18.983  21.544  1.00 20.00           C
HETATM    8  OG  SEP A   2      20.108  18.358  22.460  1.00 20.00           O
HETATM    9  P   SEP A   2      19.454  17.349  23.513  1.00 20.00           P
HETATM   10  O1P SEP A   2      18.547  16.296  22.986  1.00 20.00           O
HETATM   11  O2P SEP A   2      18.793  18.302  24.609  1.00 20.00           O
HETATM   12  O3P SEP A   2      20.721  16.782  24.300  1.00 20.00           O
HETATM   13  C   SEP A   2      20.866  19.417  19.670  1.00 20.00           C
HETATM   14  O   SEP A   2      20.700  18.287  19.218  1.00 20.00           O
ATOM     15  N   GLY A   3      21.822  20.285  19.169  1.00 20.00           N
ATOM     16  CA  GLY A   3      22.511  19.972  17.932  1.00 20.00           C
ATOM     17  C   GLY A   3      21.662  20.349  16.723  1.00 20.00           C
ATOM     18  O   GLY A   3      20.590  20.939  16.718  1.00 20.00           O
ATOM     19  OXT GLY A   3      22.242  19.981  15.560  1.00 20.00           O
HETATM   20  H   SEP A   2      18.353  20.501  19.402  1.00 20.00           H
HETATM   21  HA  SEP A   2      20.594  20.650  21.417  1.00 20.00           H
HETATM   22  HB2 SEP A   2      18.405  19.447  22.113  1.00 20.00           H
HETATM   23  HB3 SEP A   2      18.777  18.221  20.892  1.00 20.00           H
HETATM   24 HOP2 SEP A   2      18.360  17.778  25.305  1.00 20.00           H
HETATM   25 HOP3 SEP A   2      21.189  16.114  23.766  1.00 20.00           H
END
'''

# Ser A 1 O-acetyl ester: ACT B 2 without OXT, its C bonded to Ser OG (restraint
# edit below)
ester_pdb = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   1      18.708  21.872  19.624  1.00 20.00           N
ATOM      2  CA  SER A   1      18.478  20.406  19.606  1.00 20.00           C
ATOM      3  CB  SER A   1      19.740  19.647  19.173  1.00 20.00           C
ATOM      4  OG  SER A   1      20.821  19.902  20.087  1.00 20.00           O
HETATM    5  C   ACT B   2      22.015  19.380  19.703  1.00 20.00           C
HETATM    6  CH3 ACT B   2      23.072  19.723  20.707  1.00 20.00           C
HETATM    7  O   ACT B   2      22.209  18.721  18.691  1.00 20.00           O
ATOM      8  C   SER A   1      17.970  19.937  20.973  1.00 20.00           C
ATOM      9  O   SER A   1      17.488  20.643  21.848  1.00 20.00           O
ATOM     10  OXT SER A   1      17.978  18.598  21.128  1.00 20.00           O
ATOM     11  HA  SER A   1      17.670  20.202  18.893  1.00 20.00           H
ATOM     12  HB2 SER A   1      19.551  18.568  19.143  1.00 20.00           H
ATOM     13  HB3 SER A   1      20.030  19.976  18.167  1.00 20.00           H
HETATM   14  H1  ACT B   2      23.172  20.809  20.783  1.00 20.00           H
HETATM   15  H2  ACT B   2      22.815  19.291  21.678  1.00 20.00           H
HETATM   16  H3  ACT B   2      24.029  19.307  20.382  1.00 20.00           H
END
'''

ester_edits = '''pdb_interpretation.geometry_restraints.edits {
  bond {
    atom_selection_1 = chain B and resseq 2 and name C
    atom_selection_2 = chain A and resseq 1 and name OG
    distance_ideal = 1.34
    sigma = 0.02
  }
}'''

zinc_edits = '''pdb_interpretation.geometry_restraints.edits {
  bond {
    atom_selection_1 = chain A and resseq 1 and name S1
    atom_selection_2 = chain Z and resseq 1 and name ZN
    distance_ideal = 2.30
    sigma = 0.05
  }
}'''

# ------------------------------------------------------------------------------

def get_model(pdb_str, cifs=(), edits=None):
  '''cifs: (code, restraint cif text); edits: pdb_interpretation PHIL.'''
  ro = [(code + '.cif', iotbx.cif.reader(input_string=t).model()) for code, t in cifs]
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None),
    restraint_objects=ro or None, log=null_out())
  if edits is None:
    model.process(make_restraints=True)
  else:
    p = mmtbx.model.manager.get_default_pdb_interpretation_scope().fetch(
      iotbx.phil.parse(edits)).extract()
    model.process(make_restraints=True, pdb_interpretation_params=p)
  return model

def pdb_from_cif(code, cif_text, drop=(), extra=()):
  '''One residue (chain A, resseq 1) at the restraint file's coordinates.'''
  b = iotbx.cif.reader(input_string=cif_text).model()['comp_%s' % code]
  ids, el = b['_chem_comp_atom.atom_id'], b['_chem_comp_atom.type_symbol']
  x, y, z = [b['_chem_comp_atom.%s' % c] for c in 'xyz']
  out = ['CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1']
  for k in range(len(ids)):
    if ids[k] in drop:
      continue
    n = ids[k] if len(ids[k]) == 4 else ' ' + ids[k]
    out.append('HETATM%5d %-4s %3s A   1    %8.3f%8.3f%8.3f  1.00 20.00          %2s' % (
      k + 1, n, code, float(x[k]) + 20, float(y[k]) + 20, float(z[k]) + 20, el[k]))
  return '\n'.join(out + list(extra) + ['END'])

def build(model, chain='A', resseq=1, altloc=''):
  for rg in model.get_hierarchy().residue_groups():
    if rg.parent().id == chain and rg.resseq_as_int() == resseq:
      return rdkit_utils.residue_molecule(model, rg, altloc=altloc)
  raise AssertionError('no residue %s %d' % (chain, resseq))

def smiles(r):
  '''Canonical SMILES of the fragment molecule, without stereo.'''
  return Chem.MolToSmiles(Chem.RemoveHs(r.fragment_mol), isomericSmiles=False)

def names(r):
  return sorted([a.GetProp('_Name') for a in r.mol.GetAtoms()])

def charged(r):
  return sorted([(a.GetProp('_Name'), a.GetFormalCharge()) for a in r.mol.GetAtoms()
    if a.GetFormalCharge()])

def zzs_cif():
  '''ZSH with a Zn in the residue bonded to S1 ('metal').'''
  t = zsh_cif.replace('ZSH', 'ZZS').replace("'ZZS' ligand 9 3", "'ZZS' ligand 10 4")
  t = t.replace(' ZZS  HS11 H  H     0   0.000   1.9180   1.0419  -0.1009\n',
    ' ZZS  HS11 H  H     0   0.000   1.9180   1.0419  -0.1009\n'
    ' ZZS  ZN1  ZN ZN    0   0.000   3.5394   1.0507   1.8340\n')
  t = t.rstrip('\n') + '\n ZZS  S1   ZN1  metal    2.300  0.050\n'
  assert t.count('ZN1') == 2
  return t

class captured(object):
  '''Everything written to the process's stdout and stderr (file descriptors).'''
  def __enter__(self):
    sys.stdout.flush()
    sys.stderr.flush()
    self.tmp = tempfile.TemporaryFile(mode='w+')
    self.saved = [os.dup(1), os.dup(2)]
    os.dup2(self.tmp.fileno(), 1)
    os.dup2(self.tmp.fileno(), 2)
    return self
  def __exit__(self, *args):
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(self.saved[0], 1)
    os.dup2(self.saved[1], 2)
    for f in self.saved:
      os.close(f)
    self.tmp.seek(0)
    self.text = self.tmp.read()
    self.tmp.close()

# ------------------------------------------------------------------------------

def exercise_carboxylate():
  ''''deloc' C-O and -1 on one O (G39/ACT style): carboxylate, total -1.'''
  for code, text in (('ZAC', zac_cif), ('ACT', None)):
    if text is None:
      r = build(get_model(act_pdb))
    else:
      r = build(get_model(pdb_from_cif(code, text), cifs=((code, text),)))
    assert r.ok, r.reason
    assert (r.total_charge, r.total_charge_source) == (-1, 'formal charges')
    assert r.charge_certain is True and r.hydrogens == 'model'
    assert r.h_differences == [] and r.uncertain_atoms == []
    assert smiles(r) == 'CC(=O)[O-]', smiles(r)
    assert r.differences == dict(bonds=[], charges=[]), r.differences
    assert r.caps == [] and r.added_h == []
    assert sorted(r.rdkit_to_iseq.values()) == list(range(r.mol.GetNumAtoms()))
    assert r.seconds < 1
    # sanitized: hybridization assigned
    c = [a for a in r.mol.GetAtoms() if a.GetSymbol() == 'C' and
      [n for n in a.GetNeighbors() if n.GetSymbol() == 'O']][0]
    assert c.GetHybridization() == Chem.HybridizationType.SP2
  # the -1 on the other O: resonance-equivalent, no difference reported either way
  other = zac_cif.replace(' ZAC  O1   O  O     0 ', ' ZAC  O1   O  O    -1 ').replace(
    ' ZAC  O2   O  O    -1 ', ' ZAC  O2   O  O     0 ')
  assert other != zac_cif
  r = build(get_model(pdb_from_cif('ZAC', other), cifs=(('ZAC', other),)))
  assert r.ok and r.total_charge == -1 and r.differences['charges'] == [], r.differences

def exercise_zwitterion():
  '''MES-style zwitterion: sulfonate ('deloc'), ammonium; total 0.'''
  r = build(get_model(pdb_from_cif('ZZW', zzw_cif), cifs=(('ZZW', zzw_cif),)))
  assert r.ok, r.reason
  assert r.total_charge == 0
  assert smiles(r) == '[NH3+]CCS(=O)(=O)[O-]', smiles(r)
  assert [(n[0], q) for n, q in charged(r)] == [('N', 1), ('O', -1)], charged(r)
  assert r.differences['charges'] == [], r.differences

def exercise_triphosphate():
  '''Methyl triphosphate: total -4, P neutral (no P+), the charges on terminal O.'''
  r = build(get_model(pdb_from_cif('ZTP', ztp_cif), cifs=(('ZTP', ztp_cif),)))
  assert r.ok, r.reason
  assert r.total_charge == -4
  q = charged(r)
  assert len(q) == 4 and set([c for n, c in q]) == set([-1]), q
  assert [n for n, c in q if n.startswith('P')] == []
  assert r.differences['charges'] == [], r.differences

def exercise_partial_charges():
  '''
  Partial charges only, 'coval' bonds: no formal charges, so the totals -4..+4 are
  searched; -1 is the only plausible one (-3 charges a carbon), uncertain; a charge
  column of only '?' counts as missing.
  '''
  r = build(get_model(pdb_from_cif('ZPC', zpc_cif), cifs=(('ZPC', zpc_cif),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source) == (-1, 'search'), r.total_charge_source
  assert r.charge_certain is False
  assert r.search['valid'] == [-3, -1] and r.search['set_aside'] == [-3], r.search
  assert r.search['calls'] == 9
  assert smiles(r) == 'CC(=O)[O-]'
  assert [n for n in r.charge_notes if n.startswith('sum of partial charges -1.000')]
  zqm = zpc_cif.replace('ZPC', 'ZQM').replace(
    '_chem_comp_atom.partial_charge', '_chem_comp_atom.charge\n_chem_comp_atom.partial_charge')
  lines = []
  for l in zqm.split('\n'):
    if l.startswith(' ZQM ') and len(l.split()) == 8:
      p = l.split()
      l = ' '.join(p[:4] + ['?'] + p[4:])
    lines.append(l)
  zqm = '\n'.join(lines)
  r = build(get_model(pdb_from_cif('ZQM', zqm), cifs=(('ZQM', zqm),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source, r.charge_certain) == (-1, 'search', False)

def exercise_chain_residue():
  '''
  SEP between two Gly: caps on N and C (polymer bonds), OXT absent and not capped,
  the phosphate as modelled (HOP2 and HOP3: neutral; GeoStd SEP has no charge
  column, so the total comes from the search: 0 is the only valid one; -2, with a
  single C-O, contradicts the file's C=O).
  '''
  model = get_model(chain_pdb)
  r = build(model, resseq=2)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert sorted([(c['kind'], atoms[c['on']].name.strip(), atoms[c['partner']].name.strip(),
    atoms[c['partner']].parent().parent().resseq_as_int()) for c in r.caps]) == [
    ('linked', 'C', 'N', 3), ('linked', 'N', 'C', 1)], r.caps
  assert sorted([atoms[i].name.strip() for i in r.linked]) == ['C', 'N']
  assert r.missing_neighbour == []
  assert 'OXT' not in names(r)
  assert (r.total_charge, r.total_charge_source) == (0, 'search'), r.search
  assert r.search['disagree'] == {-2: ['C-O: restraint file 2, RDKit 1']}, r.search
  assert r.hydrogens == 'model' and r.h_differences == []
  assert charged(r) == []
  p = [b for b in r.mol.GetBonds() if 'P' in (b.GetBeginAtom().GetSymbol(),
    b.GetEndAtom().GetSymbol()) and b.GetBondType() == Chem.BondType.DOUBLE]
  assert len(p) == 1
  caps = [a for a in r.mol.GetAtoms() if a.HasProp('cap')]
  assert len(caps) == 2 and r.fragment_mol.GetNumAtoms() == r.mol.GetNumAtoms() - 2
  assert [r.rdkit_to_iseq.get(a.GetIdx()) for a in caps] == [None, None]
  # without any H: the restraint file's H; the peptide N keeps its one H (the
  # polymer bond does not replace it)
  no_h = '\n'.join([l for l in chain_pdb.split('\n') if not (l.startswith('HETATM') and
    l[76:78].strip() == 'H')])
  model = get_model(no_h)
  r = build(model, resseq=2)
  assert r.ok, r.reason
  assert len(r.added_h) == 6 and len(r.caps) == 2
  assert r.hydrogens == 'restraint file (no H in the model)' and r.h_differences == []
  assert sorted(r.uncertain_atoms) == ['N', 'O2P', 'O3P'], r.uncertain_atoms
  n = [a for a in r.mol.GetAtoms() if a.HasProp('_Name') and a.GetProp('_Name') == 'N'][0]
  assert sorted([x.GetSymbol() for x in n.GetNeighbors()]) == ['C', 'H', 'H']
  assert charged(r) == []

def exercise_covalent_link():
  '''
  An acyl ester to Ser OG: a cap on the ligand's C, OXT (leaving) not capped: no
  carboxylate. An ether GOL C3-O1 GOL: the link replaces the leaving O3, C3 keeps
  both H.
  '''
  model = get_model(ester_pdb, edits=ester_edits)
  r = build(model, chain='B', resseq=2)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert [(c['kind'], atoms[c['on']].name.strip(), atoms[c['partner']].name.strip())
    for c in r.caps] == [('linked', 'C', 'OG')]
  assert r.missing_neighbour == []
  assert r.total_charge == 0 and charged(r) == []
  assert smiles(r) == 'CC=O', smiles(r)
  # GOL A without O3, its C3 bonded to O1 of GOL B (O3 the leaving atom): the link
  # replaces O3, not an H, so C3 keeps H31 and H32 (model H, and without H)
  a = [l for l in gol_pdb.split('\n') if l[12:16].strip() not in ('O3', 'HO3', 'END')]
  b = [l[:21] + 'B' + l[22:30] + '%8.3f' % (float(l[30:38]) + 6.0) + l[38:]
    for l in gol_pdb.split('\n') if l.startswith('HETATM')]
  edits = zinc_edits.replace('chain A and resseq 1 and name S1',
    'chain A and resseq 1 and name C3').replace('chain Z and resseq 1 and name ZN',
    'chain B and resseq 1 and name O1').replace('2.30', '1.43')
  model = get_model('\n'.join(a + b + ['END']), edits=edits)
  r = build(model)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert [(c['kind'], atoms[c['on']].name.strip(), atoms[c['partner']].name.strip())
    for c in r.caps] == [('linked', 'C3', 'O1')]
  assert r.hydrogens == 'model' and r.h_differences == []
  assert smiles(r) == 'CC(O)CO', smiles(r)
  no_h = [l for l in a if not (l.startswith('HETATM') and l[12:16].strip().startswith('H'))]
  r = build(get_model('\n'.join(no_h + b + ['END']), edits=edits))
  assert r.ok, r.reason
  assert sorted([r.mol.GetAtomWithIdx(k).GetProp('_Name') for k in r.added_h]) == [
    'H11', 'H12', 'H2', 'H31', 'H32', 'HO1', 'HO2']

def exercise_metal():
  '''A thiolate on Zn: the thiol H missing on a metal-bound S counts as a deprotonation.'''
  zn = 'HETATM   99 ZN    ZN Z   1    %8.3f%8.3f%8.3f  1.00 20.00          ZN'
  b = iotbx.cif.reader(input_string=zsh_cif).model()['comp_ZSH']
  xyz = dict([(n, [float(b['_chem_comp_atom.%s' % c][k]) + 20 for c in 'xyz'])
    for k, n in enumerate(b['_chem_comp_atom.atom_id'])])
  s, c = xyz['S1'], xyz['C1']
  d = [s[k] - c[k] for k in range(3)]
  n = sum([v * v for v in d]) ** 0.5
  site = [s[k] + 2.3 * d[k] / n for k in range(3)]
  pdb = pdb_from_cif('ZSH', zsh_cif, drop=('HS11',), extra=(zn % tuple(site),))
  model = get_model(pdb, cifs=(('ZSH', zsh_cif),), edits=zinc_edits)
  r = build(model)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert [atoms[i].name.strip() for i in r.metal_bound] == ['S1']
  assert r.total_charge == -1
  assert [x for x in r.charge_notes if 'deprotonation' in x], r.charge_notes
  assert r.h_differences == [dict(atom='S1', model=0, restraint_file=1,
    kind='deprotonation', added=[])], r.h_differences
  assert r.hydrogens == 'model' and r.uncertain_atoms == []
  assert smiles(r) == 'CC[S-]', smiles(r)
  assert 'ZN' not in [a.GetSymbol().upper() for a in r.mol.GetAtoms()]
  # the Zn belongs to another residue: not in the fragments
  rg = [g for g in model.get_hierarchy().residue_groups() if g.parent().id == 'A'][0]
  rc = rdkit_utils.residue_rigid_components(model, rg)
  assert rc.approximate is None
  assert [sorted([atoms[i].name.strip() for i in c]) for c in rc.components] == [
    ['C1', 'C2', 'HC11', 'HC12', 'HC21', 'HC22', 'HC23', 'S1']]

def exercise_no_h():
  '''A residue without H: the restraint file's H added (acetate: three).'''
  pdb = '\n'.join([l for l in act_pdb.split('\n') if not l[12:16].strip().startswith('H')])
  r = build(get_model(pdb))
  assert r.ok, r.reason
  assert len(r.added_h) == 3 and r.total_charge == -1
  assert r.hydrogens == 'restraint file (no H in the model)'
  assert r.h_differences == [] and r.uncertain_atoms == []
  assert smiles(r) == 'CC(=O)[O-]'
  assert [r.rdkit_to_iseq.get(k) for k in r.added_h] == [None] * 3

def exercise_missing_heavy_atoms():
  '''
  GOL without O3: a cap on C3 (C3 keeps its two H), HO3 left out (noted); acetate
  without OXT: a cap on C, total 0 (OXT's -1 drops out).
  '''
  pdb = '\n'.join([l for l in gol_pdb.split('\n') if l[12:16].strip() != 'O3'])
  model = get_model(pdb)
  r = build(model)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert [(c['kind'], atoms[c['on']].name.strip(), c['partner']) for c in r.caps] == [
    ('missing', 'C3', 'O3')]
  assert [atoms[i].name.strip() for i in r.missing_neighbour] == ['C3']
  assert r.total_charge == 0 and smiles(r) == 'CC(O)CO', smiles(r)
  assert 'HO3' not in names(r)
  assert 'H without its heavy atom, left out: HO3' in r.charge_notes, r.charge_notes
  pdb = '\n'.join([l for l in act_pdb.split('\n') if l[12:16].strip() != 'OXT'])
  r = build(get_model(pdb))
  assert r.ok, r.reason
  assert [(c['kind'], c['partner']) for c in r.caps] == [('missing', 'OXT')]
  assert r.total_charge == 0 and smiles(r) == 'CC=O', smiles(r)

def exercise_split_oxygen():
  '''ACT with OXT in two conformers: one molecule per conformer, each with its own OXT.'''
  lines = []
  for l in act_pdb.split('\n'):
    if l[12:16].strip() == 'OXT':
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      x = float(l[30:38]) + 0.3
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % x + l[38:54] + '  0.50' + l[60:])
    else:
      lines.append(l)
  model = get_model('\n'.join(lines))
  atoms = model.get_hierarchy().atoms()
  for alt in ('A', 'B'):
    r = build(model, altloc=alt)
    assert r.ok, r.reason
    assert r.total_charge == -1 and smiles(r) == 'CC(=O)[O-]'
    oxt = [i for k, i in r.rdkit_to_iseq.items() if atoms[i].name.strip() == 'OXT']
    assert [atoms[i].parent().altloc for i in oxt] == [alt]
    assert r.mol.GetNumAtoms() == 7

def exercise_formal_total():
  '''
  With formal charges the total is fixed (no search): ZAH (-1 on O2, the model has
  the file's HO21) and ZNC (-1 on both O) fail with the reason, after the input,
  canonical and 10 random atom orders (12 calls).
  '''
  r = build(get_model(pdb_from_cif('ZAH', zah_cif), cifs=(('ZAH', zah_cif),)))
  assert not r.ok and r.mol is None and r.fragment_mol is None
  assert r.reason.startswith('DetermineBondOrders fails for ZAH with the formal total -1:'), \
    r.reason
  assert r.search['calls'] == 12 and r.total_charge is None and r.charge_certain is None
  assert r.reason.endswith(' (also in the canonical and 10 random atom orders)'), r.reason
  r = build(get_model(pdb_from_cif('ZNC', znc_cif), cifs=(('ZNC', znc_cif),)))
  assert not r.ok
  assert r.reason.startswith('DetermineBondOrders fails for ZNC with the formal total -2:'), \
    r.reason
  assert r.search['calls'] == 12 and r.search['order'] is None

def exercise_file_first():
  '''
  Every bond between the heavy atoms present has an explicit order and the file has
  formal charges: the molecule from the file (ZAA acetic acid; ZSH thiolate on a
  Zn of the residue, the deprotonation applied to S1), no DetermineBondOrders call.
  Otherwise DetermineBondOrders, retried in other atom orders: GeoStd AQS with
  C1A-C2A made 'coval' fails in the input order at +2 and succeeds in a random one.
  '''
  zaa = zah_cif.replace('ZAH', 'ZAA').replace(' ZAA  O2   O  O    -1 ',
    ' ZAA  O2   O  O     0 ').replace(' ZAA  C1   O1   deloc ', ' ZAA  C1   O1   double'
    ).replace(' ZAA  C1   O2   deloc ', ' ZAA  C1   O2   single')
  r = build(get_model(pdb_from_cif('ZAA', zaa), cifs=(('ZAA', zaa),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source, r.charge_certain) == (0,
    'restraint file', True)
  assert r.search['calls'] == 0 and smiles(r) == 'CC(=O)O'
  model = get_model(pdb_from_cif('ZZS', zzs_cif(), drop=('HS11',)), cifs=(('ZZS', zzs_cif()),))
  r = build(model)
  assert r.ok and (r.total_charge, r.total_charge_source) == (-1, 'restraint file'), r.reason
  assert charged(r) == [('S1', -1)]
  # ZAC ('deloc' C-O): DetermineBondOrders
  r = build(get_model(pdb_from_cif('ZAC', zac_cif), cifs=(('ZAC', zac_cif),)))
  assert r.total_charge_source == 'formal charges' and r.search['order'] == 'input'
  import re
  from mmtbx.monomer_library import server
  path = os.path.join(server.server().geostd_path, 'a', 'data_AQS.cif')
  text = open(path).read()
  coval = re.sub(r'( AQS\s+C1A\s+C2A\s+)single', r'\1coval ', text)
  assert coval != text
  r = build(get_model(pdb_from_cif('AQS', text), cifs=(('AQS', text),)))
  assert r.ok and r.total_charge_source == 'restraint file', r.reason
  r = build(get_model(pdb_from_cif('AQS', coval), cifs=(('AQS', coval),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source) == (2, 'formal charges')
  assert r.search['order'].startswith('random (seed ') and r.search['calls'] > 2, r.search
  assert 'DetermineBondOrders succeeded in the %s atom order' % r.search['order'] in \
    r.charge_notes
  assert r.differences['bonds'] == [], r.differences

def exercise_peptide_ends():
  '''
  A peptide restraint file whose unlinked N or C has a valence open gets a cap for
  the absent neighbour residue: GeoStd 4GJ (N-CA and N-H single, no H2) alone:
  capped N, ok; GeoStd MH6 (imine N=CA, no H2): no cap, ok.
  '''
  from mmtbx.monomer_library import server
  geostd = server.server().geostd_path
  for code, capped in (('4GJ', True), ('MH6', False)):
    text = open(os.path.join(geostd, code[0].lower(), 'data_%s.cif' % code)).read()
    model = get_model(pdb_from_cif(code, text), cifs=((code, text),))
    r = build(model)
    assert r.ok, (code, r.reason)
    atoms = model.get_hierarchy().atoms()
    caps = [(atoms[c['on']].name.strip(), c['partner']) for c in r.caps]
    assert caps == ([('N', '(no preceding residue)')] if capped else []), (code, caps)

def exercise_search():
  '''
  No formal charges: the totals -4..+4, one plausible total taken, uncertain. ZGU
  (partial sum 1.548): +1 (-1 contradicts the file's C2=N2); ZNM: 0 (-2 as
  CN([O-])[O-] contradicts the file's N1=O1); ZNM with 'coval' bonds: -2 and 0 both
  valid, ambiguous.
  '''
  r = build(get_model(pdb_from_cif('ZGU', zgu_cif), cifs=(('ZGU', zgu_cif),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source, r.charge_certain) == (1, 'search', False)
  assert 'sum of partial charges 1.548' in r.charge_notes, r.charge_notes
  assert smiles(r) in ('CNC(=[NH2+])N', 'CN=C([NH3+])N', 'C[NH+]=C(N)N', 'CNC(N)=[NH2+]'), smiles(r)
  assert r.search['calls'] == 9 and r.search['valid'] == [1], r.search
  assert r.search['disagree'] == {-1: ['C2-N2: restraint file 2, RDKit 1']}, r.search
  assert r.search['seconds'] >= 0
  r = build(get_model(pdb_from_cif('ZNM', znm_cif), cifs=(('ZNM', znm_cif),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source, r.charge_certain) == (0, 'search', False)
  assert r.search['valid'] == [0], r.search
  assert r.search['disagree'] == {-2: ['N1-O1: restraint file 2, RDKit 1']}, r.search
  assert smiles(r) == 'C[N+](=O)[O-]'
  coval = znm_cif.replace(' double ', ' coval  ').replace(' single ', ' coval  ')
  assert coval != znm_cif
  r = build(get_model(pdb_from_cif('ZNM', coval), cifs=(('ZNM', coval),)))
  assert not r.ok and r.reason == 'ambiguous total charge for ZNM (no formal charges): ' \
    '-2 +0', r.reason

def exercise_hydrogen_contract():
  '''
  The restraint file is the reference protonation. ZAC without HC21: HC21 added,
  total -1, C-C single, certain. ZAA (acetic acid) without the O-H: HO21 added, O2
  uncertain. A ZSH disulfide with the thiol H kept: the link cap replaces HS11, so
  the model's H is extra: failure, approximate fragments flagged.
  '''
  r = build(get_model(pdb_from_cif('ZAC', zac_cif, drop=('HC21',)), cifs=(('ZAC', zac_cif),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.total_charge_source, r.charge_certain) == (-1, 'formal charges',
    True)
  assert smiles(r) == 'CC(=O)[O-]', smiles(r)
  c = [b for b in r.mol.GetBonds() if b.GetBeginAtom().GetSymbol() == 'C' and
    b.GetEndAtom().GetSymbol() == 'C']
  assert len(c) == 1 and c[0].GetBondType() == Chem.BondType.SINGLE
  assert r.hydrogens == 'completed from the restraint file'
  assert r.h_differences == [dict(atom='C2', model=2, restraint_file=3, kind='added',
    added=['HC21'])], r.h_differences
  assert r.uncertain_atoms == []
  assert [r.mol.GetAtomWithIdx(k).GetProp('_Name') for k in r.added_h] == ['HC21']
  assert [r.rdkit_to_iseq.get(k) for k in r.added_h] == [None]
  assert 'H added from the restraint file: HC21' in r.charge_notes, r.charge_notes
  zaa = zah_cif.replace('ZAH', 'ZAA').replace(' ZAA  O2   O  O    -1 ',
    ' ZAA  O2   O  O     0 ').replace(' ZAA  C1   O1   deloc ', ' ZAA  C1   O1   double'
    ).replace(' ZAA  C1   O2   deloc ', ' ZAA  C1   O2   single')
  assert zaa.count('double') == 1 and zaa.count('-1 ') == 0
  r = build(get_model(pdb_from_cif('ZAA', zaa, drop=('HO21',)), cifs=(('ZAA', zaa),)))
  assert r.ok, r.reason
  assert (r.total_charge, r.charge_certain) == (0, True)
  assert smiles(r) == 'CC(=O)O', smiles(r)
  assert r.h_differences == [dict(atom='O2', model=0, restraint_file=1, kind='added',
    added=['HO21'])], r.h_differences
  assert r.uncertain_atoms == ['O2']
  # disulfide: two ZSH, S1-S1 bonded by an edit
  a = pdb_from_cif('ZSH', zsh_cif).split('\n')
  b = [l[:21] + 'B' + l[22:30] + '%8.3f' % (-float(l[30:38]) + 45.4686) + l[38:]
    for l in a if l.startswith('HETATM')]
  edits = zinc_edits.replace('chain Z and resseq 1 and name ZN',
    'chain B and resseq 1 and name S1').replace('2.30', '2.05')
  model = get_model('\n'.join(a[:-1] + b + ['END']), cifs=(('ZSH', zsh_cif),), edits=edits)
  rg = [g for g in model.get_hierarchy().residue_groups() if g.parent().id == 'A'][0]
  with captured() as c:
    rc = rdkit_utils.residue_rigid_components(model, rg)
  r = rc.molecule
  assert not r.ok and r.reason == 'H not in the restraint file: S1: 1 H in the model, ' \
    '0 in the restraint file', r.reason
  assert r.h_differences == [dict(atom='S1', model=1, restraint_file=0, kind='extra',
    added=[])], r.h_differences
  assert rc.approximate == 'approximate: ' + r.reason
  assert sorted([i for comp in rc.components for i in comp]) == list(range(9))
  assert c.text == '', repr(c.text)

def exercise_bond_order_agreement():
  '''
  A structure is accepted only if its bonds agree with the file's explicit orders:
  ZAC with C2-C1 'double' (the model's CH3 allows only a single bond): failure
  listing the bond. Resonance partners apart (ZGU's C2=N3 against the file's C2=N2,
  exercise_search).
  '''
  cc = zac_cif.replace(' ZAC  C2   C1   single ', ' ZAC  C2   C1   double ')
  assert cc != zac_cif
  r = build(get_model(pdb_from_cif('ZAC', cc), cifs=(('ZAC', cc),)))
  assert not r.ok and r.reason == 'bond orders disagree with the restraint file for ZAC ' \
    'at the formal total -1: C2-C1: restraint file 2, RDKit 1 (also in the canonical ' \
    'and 10 random atom orders)', r.reason

def exercise_metal_fragments():
  '''
  The residue's metals are in the fragments (dative bonds from the residue's atoms,
  not cut). HEM (no H): all 43 heavy atoms, Fe with the porphyrin, as with the
  pre-B1 route (22f62b3fc6). ZZS (ZSH with its Zn in the residue, thiol H missing):
  deprotonation, S1 with the Zn (as pre-B1). OFO: from the file (oxide -2, hydroxide
  -1), one component (as pre-B1); with the oxide's charge '?' (counted 0) no
  structure, and the approximate molecule keeps the metals (dative bonds): one
  component.
  '''
  model = get_model(hem_pdb)
  atoms = model.get_hierarchy().atoms()
  rc = rdkit_utils.residue_rigid_components(model, model.get_hierarchy().only_residue_group())
  assert rc.approximate is None, rc.approximate
  part = sorted([sorted([atoms[i].name.strip() for i in c]) for c in rc.components])
  assert part == [['C1A', 'C1B', 'C1C', 'C1D', 'C2A', 'C2B', 'C2C', 'C2D', 'C3A', 'C3B',
    'C3C', 'C3D', 'C4A', 'C4B', 'C4C', 'C4D', 'CAA', 'CAD', 'CHA', 'CHB', 'CHC', 'CHD',
    'CMA', 'CMB', 'CMC', 'CMD', 'FE', 'NA', 'NB', 'NC', 'ND'], ['CAB', 'CBB'],
    ['CAC', 'CBC'], ['CBA', 'CGA', 'O1A', 'O2A'], ['CBD', 'CGD', 'O1D', 'O2D']], part
  r = rc.molecule
  assert sorted([atoms[i].name.strip() for i in r.metal_bound]) == ['NA', 'NB', 'NC', 'ND']
  fe = [a for a in r.fragment_mol.GetAtoms() if a.GetSymbol() == 'Fe']
  assert len(fe) == 1 and atoms[r.fragment_to_iseq[fe[0].GetIdx()]].name.strip() == 'FE'
  assert sorted([(b.GetBeginAtom().GetSymbol(), str(b.GetBondType()))
    for b in fe[0].GetBonds()]) == [('N', 'DATIVE')] * 4
  assert 'Fe' not in [a.GetSymbol() for a in r.mol.GetAtoms()]
  assert r.charge_certain is False
  # Zn in the residue
  zzs = zzs_cif()
  model = get_model(pdb_from_cif('ZZS', zzs, drop=('HS11',)), cifs=(('ZZS', zzs),))
  atoms = model.get_hierarchy().atoms()
  rc = rdkit_utils.residue_rigid_components(model, model.get_hierarchy().only_residue_group())
  r = rc.molecule
  assert r.ok and rc.approximate is None, r.reason
  assert r.total_charge == -1 and [d['kind'] for d in r.h_differences] == ['deprotonation']
  assert Chem.MolToSmiles(Chem.RemoveHs(r.fragment_mol)) == 'CC[S-]->[Zn]', \
    Chem.MolToSmiles(r.fragment_mol)
  part = sorted([sorted([atoms[i].name.strip() for i in c if atoms[i].element.strip() != 'H'])
    for c in rc.components])
  assert part == [['C1', 'C2'], ['S1', 'ZN1']], part
  assert sorted([i for c in rc.components for i in c]) == list(range(model.get_number_of_atoms()))
  # OFO: from the file
  model = get_model(pdb_from_cif('OFO', ofo_cif), cifs=(('OFO', ofo_cif),))
  atoms = model.get_hierarchy().atoms()
  rc = rdkit_utils.residue_rigid_components(model, model.get_hierarchy().only_residue_group())
  assert rc.approximate is None and rc.molecule.total_charge_source == 'restraint file'
  assert Chem.MolToSmiles(rc.molecule.fragment_mol) == '[H][O-]->[Fe]<-[O-2]->[Fe]', \
    Chem.MolToSmiles(rc.molecule.fragment_mol)
  assert [sorted([atoms[i].name.strip() for i in comp]) for comp in rc.components] == [
    ['FE1', 'FE2', 'HO', 'O', 'OH']], rc.components
  # the oxide's charge '?': approximate, with the metals
  q = ofo_cif.replace(' OFO  O    O  O    -2 ', ' OFO  O    O  O     ? ')
  assert q != ofo_cif
  model = get_model(pdb_from_cif('OFO', q), cifs=(('OFO', q),))
  atoms = model.get_hierarchy().atoms()
  rg = model.get_hierarchy().only_residue_group()
  with captured() as c:
    rc = rdkit_utils.residue_rigid_components(model, rg)
  assert rc.approximate.startswith('approximate: DetermineBondOrders fails for OFO with '
    'the formal total -1:'), rc.approximate
  assert [sorted([atoms[i].name.strip() for i in comp]) for comp in rc.components] == [
    ['FE1', 'FE2', 'HO', 'O', 'OH']], rc.components
  mol, rdkit_to_iseq = rdkit_utils.approximate_residue_molecule(model, rg)
  assert sorted([atoms[i].name.strip() for i in rdkit_to_iseq.values()]) == [
    'FE1', 'FE2', 'HO', 'O', 'OH']
  bonds = sorted([(b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol(),
    str(b.GetBondType())) for b in mol.GetBonds()])
  assert bonds == [('O', 'Fe', 'DATIVE')] * 3 + [('O', 'H', 'SINGLE')], bonds
  assert c.text == '', repr(c.text)

def exercise_split_neighbour():
  '''The Gly after SEP split into A and B: SEP's C still gets one cap (one partner N).'''
  lines = []
  for l in chain_pdb.split('\n'):
    if l.startswith('ATOM') and l[22:26] == '   3':
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.3) + l[38:54] +
        '  0.50' + l[60:])
    else:
      lines.append(l)
  model = get_model('\n'.join(lines))
  r = build(model, resseq=2)
  assert r.ok, r.reason
  atoms = model.get_hierarchy().atoms()
  assert sorted([(atoms[c['on']].name.strip(), atoms[c['partner']].name.strip())
    for c in r.caps]) == [('C', 'N'), ('N', 'C')], r.caps

def exercise_failure():
  '''
  A carbon with five bonds (ZC5): DetermineBondOrders fails at the formal total
  (no search); failure with the reason, no molecule, nothing printed.
  '''
  model = get_model(pdb_from_cif('ZC5', zc5_cif), cifs=(('ZC5', zc5_cif),))
  with captured() as c:
    r = build(model)
  assert not r.ok and r.mol is None and r.fragment_mol is None
  assert r.reason.startswith('DetermineBondOrders fails for ZC5 with the formal total 0:'), \
    r.reason
  assert r.search['calls'] == 12 and r.search['valid'] == []
  assert c.text == '', repr(c.text)

def exercise_rigid_components():
  '''
  residue_rigid_components: fragments from the builder's fragment_mol (ZTP: the
  phosphates split off; added H have no i_seq), and approximate fragments with the
  reason when the builder fails (ZC5), also written below the PNG.
  '''
  model = get_model(pdb_from_cif('ZTP', ztp_cif), cifs=(('ZTP', ztp_cif),))
  rg = model.get_hierarchy().only_residue_group()
  rc = rdkit_utils.residue_rigid_components(model, rg)
  assert rc.approximate is None and rc.molecule.ok
  n = sum([c.size() for c in rc.components])
  assert n == model.get_number_of_atoms(), (n, model.get_number_of_atoms())
  assert len(rc.components) > 1
  assert rdkit_utils.get_cctbx_isel_for_rigid_components(model, rg)[0].size() == \
    rc.components[0].size()
  # without H: added H are in the molecule but not in the components (ZTP is cut)
  b = iotbx.cif.reader(input_string=ztp_cif).model()['comp_ZTP']
  hs = [n for n, e in zip(b['_chem_comp_atom.atom_id'], b['_chem_comp_atom.type_symbol'])
    if e == 'H']
  model = get_model(pdb_from_cif('ZTP', ztp_cif, drop=hs), cifs=(('ZTP', ztp_cif),))
  rc = rdkit_utils.residue_rigid_components(model, model.get_hierarchy().only_residue_group())
  assert len(rc.molecule.added_h) == 3 and len(rc.components) > 1
  assert sorted([i for c in rc.components for i in c]) == list(range(model.get_number_of_atoms()))
  # failure: approximate fragments, the reason kept
  model = get_model(pdb_from_cif('ZC5', zc5_cif), cifs=(('ZC5', zc5_cif),))
  rg = model.get_hierarchy().only_residue_group()
  with captured() as c:
    rc = rdkit_utils.residue_rigid_components(model, rg)
  assert rc.approximate.startswith('approximate: DetermineBondOrders fails for ZC5'), \
    rc.approximate
  assert sorted([i for comp in rc.components for i in comp]) == list(range(6))
  assert c.text == '', repr(c.text)
  png = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
  png.close()
  try:
    model = get_model(pdb_from_cif('ZTP', ztp_cif), cifs=(('ZTP', ztp_cif),))
    rc = rdkit_utils.residue_rigid_components(model, model.get_hierarchy().only_residue_group())
    from PIL import Image
    rdkit_utils.draw_colored_fragments(rc.mol, rc.frags, png.name)
    h_plain = Image.open(png.name).size[1]
    rdkit_utils.draw_colored_fragments(rc.mol, rc.frags, png.name,
      note='approximate: a reason written below the figure')
    assert Image.open(png.name).size[1] > h_plain
  finally:
    os.unlink(png.name)

def run():
  exercise_carboxylate()
  exercise_zwitterion()
  exercise_triphosphate()
  exercise_partial_charges()
  exercise_chain_residue()
  exercise_covalent_link()
  exercise_metal()
  exercise_no_h()
  exercise_missing_heavy_atoms()
  exercise_split_oxygen()
  exercise_formal_total()
  exercise_file_first()
  exercise_peptide_ends()
  exercise_search()
  exercise_hydrogen_contract()
  exercise_bond_order_agreement()
  exercise_metal_fragments()
  exercise_split_neighbour()
  exercise_failure()
  exercise_rigid_components()
  print('OK')

if __name__ == '__main__':
  run()
