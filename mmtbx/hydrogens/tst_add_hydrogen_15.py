from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from libtbx import group_args
from scitbx.array_family import flex

str_1S6_B_401_4jxg = """
data_default
_cell.length_a                    26.871
_cell.length_b                    26.430
_cell.length_c                    30.452
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      21627.027
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                21.55
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   1S6
   SER

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   ATOM 6 CA . SER B 64 ? 8.81200 8.39800 16.90000 1.000 15.06000 C ? A ? 2 CA 1
   ATOM 7 C . SER B 64 ? 7.91100 7.42700 16.14500 1.000 14.26000 C ? A ? 2 C 1
   ATOM 8 O . SER B 64 ? 6.73200 7.73500 16.00400 1.000 14.43000 O ? A ? 2 O 1
   ATOM 9 CB . SER B 64 ? 8.99800 9.60700 16.06900 1.000 14.19000 C ? A ? 2 CB 1
   ATOM 10 OG . SER B 64 ? 9.58500 10.62900 16.88900 1.000 18.69000 O ? A ? 2 OG 1
   HETATM 120 N . 1S6 B 401 ? 12.79500 11.93900 17.59000 1.000 18.86000 N ? B ? . N 1
   HETATM 121 CA . 1S6 B 401 ? 11.42400 12.12000 17.09200 1.000 19.06000 C ? B ? . CA 1
   HETATM 122 C . 1S6 B 401 ? 10.94500 10.77700 16.75200 1.000 20.57000 C ? B ? . C 1
   HETATM 123 CB . 1S6 B 401 ? 11.11600 12.95200 15.88600 1.000 20.89000 C ? B ? . CB 1
   HETATM 124 OXT . 1S6 B 401 ? 11.61000 9.86800 16.33900 1.000 18.23000 O ? B ? . OXT 1
   HETATM 125 CAA . 1S6 B 401 ? 14.03300 11.13000 21.26800 1.000 32.84000 C ? B ? . CAA 1
   HETATM 126 CAB . 1S6 B 401 ? 14.15300 12.81400 13.28200 1.000 25.99000 C ? B ? . CAB 1
   HETATM 127 CAC . 1S6 B 401 ? 12.14600 13.93700 12.44400 1.000 25.23000 C ? B ? . CAC 1
   HETATM 128 CAI . 1S6 B 401 ? 16.74200 14.85600 14.79300 1.000 40.47000 C ? B ? . CAI 1
   HETATM 129 CAJ . 1S6 B 401 ? 17.25700 13.55800 14.88000 1.000 42.42000 C ? B ? . CAJ 1
   HETATM 130 CAK . 1S6 B 401 ? 15.94700 15.38000 15.83700 1.000 37.93000 C ? B ? . CAK 1
   HETATM 131 CAL . 1S6 B 401 ? 16.99500 12.74200 16.00800 1.000 40.13000 C ? B ? . CAL 1
   HETATM 132 CAM . 1S6 B 401 ? 15.73300 14.56600 16.96500 1.000 33.61000 C ? B ? . CAM 1
   HETATM 133 CAT . 1S6 B 401 ? 11.40700 11.19200 12.64000 1.000 27.95000 C ? B ? . CAT 1
   HETATM 134 CAU . 1S6 B 401 ? 13.30600 12.76400 18.43000 1.000 25.25000 C ? B ? . CAU 1
   HETATM 135 CAV . 1S6 B 401 ? 14.85300 11.69600 20.19200 1.000 29.93000 C ? B ? . CAV 1
   HETATM 136 CAW . 1S6 B 401 ? 16.22500 13.25900 17.04000 1.000 34.14000 C ? B ? . CAW 1
   HETATM 137 CAX . 1S6 B 401 ? 15.94500 12.52600 18.25000 1.000 30.72000 C ? B ? . CAX 1
   HETATM 138 CAY . 1S6 B 401 ? 14.64300 12.41800 18.95200 1.000 26.52000 C ? B ? . CAY 1
   HETATM 139 CBB . 1S6 B 401 ? 11.93000 11.82800 13.87300 1.000 23.23000 C ? B ? . CBB 1
   HETATM 140 CBC . 1S6 B 401 ? 12.69300 13.18100 13.63400 1.000 25.84000 C ? B ? . CBC 1
   HETATM 141 NAN . 1S6 B 401 ? 16.89500 11.99800 19.11900 1.000 28.67000 N ? B ? . NAN 1
   HETATM 142 NAP . 1S6 B 401 ? 10.81600 12.16100 14.72600 1.000 22.89000 N ? B ? . NAP 1
   HETATM 143 OAE . 1S6 B 401 ? 12.15000 10.48600 11.98300 1.000 29.05000 O ? B ? . OAE 1
   HETATM 144 OAF . 1S6 B 401 ? 12.77300 13.75600 18.76700 1.000 24.58000 O ? B ? . OAF 1
   HETATM 145 OAH . 1S6 B 401 ? 10.23400 11.45800 12.25800 1.000 26.16000 O ? B ? . OAH 1
   HETATM 146 OAQ . 1S6 B 401 ? 16.23400 11.51900 20.23100 1.000 35.16000 O ? B ? . OAQ 1
   HETATM 147 SAR . 1S6 B 401 ? 12.56500 13.90400 15.25800 1.000 25.33000 S ? B ? . SAR 1
"""

str_3W9_A_601_4x0u = """
data_default
_cell.length_a                    25.096
_cell.length_b                    22.062
_cell.length_c                    30.107
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      16669.281
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                26.21
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   3W9
   CYS

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   ATOM 66 N . CYS A 302 ? 8.55600 8.63000 20.64900 1.000 24.10815 N ? A ? 8 N 1
   ATOM 67 CA . CYS A 302 ? 8.21100 7.20500 20.49100 1.000 23.66863 C ? A ? 8 CA 1
   ATOM 68 C . CYS A 302 ? 9.07800 6.44400 21.47600 1.000 21.50521 C ? A ? 8 C 1
   ATOM 69 O . CYS A 302 ? 8.74900 5.33300 21.81800 1.000 25.28987 O ? A ? 8 O 1
   ATOM 70 CB . CYS A 302 ? 8.34800 6.68500 19.04800 1.000 27.98757 C ? A ? 8 CB 1
   ATOM 71 SG . CYS A 302 ? 10.02600 6.75500 18.49100 1.000 38.77044 S ? A ? 8 SG 1
   HETATM 72 C03 . 3W9 A 601 ? 14.40300 10.25200 13.03700 1.000 46.87404 C ? B ? . C03 1
   HETATM 73 C04 . 3W9 A 601 ? 14.32800 7.79100 12.50800 1.000 41.86555 C ? B ? . C04 1
   HETATM 74 C05 . 3W9 A 601 ? 12.96900 8.72600 14.39700 1.000 44.66325 C ? B ? . C05 1
   HETATM 75 C06 . 3W9 A 601 ? 15.64800 10.53800 13.85800 1.000 47.41621 C ? B ? . C06 1
   HETATM 76 C07 . 3W9 A 601 ? 13.32700 7.46400 11.41600 1.000 43.70261 C ? B ? . C07 1
   HETATM 77 C08 . 3W9 A 601 ? 12.61400 9.79600 15.21900 1.000 42.65775 C ? B ? . C08 1
   HETATM 78 C09 . 3W9 A 601 ? 12.41100 7.47100 14.65100 1.000 45.77128 C ? B ? . C09 1
   HETATM 79 C10 . 3W9 A 601 ? 11.71900 9.61800 16.27100 1.000 40.92070 C ? B ? . C10 1
   HETATM 80 C11 . 3W9 A 601 ? 11.51600 7.28800 15.70400 1.000 42.76039 C ? B ? . C11 1
   HETATM 81 C12 . 3W9 A 601 ? 11.17000 8.36300 16.51700 1.000 41.86291 C ? B ? . C12 1
   HETATM 82 C13 . 3W9 A 601 ? 10.20300 8.20300 17.65700 1.000 43.81052 C ? B ? . C13 1
   HETATM 83 N02 . 3W9 A 601 ? 13.88000 8.91600 13.33000 1.000 43.96053 N ? B ? . N02 1
   HETATM 84 O01 . 3W9 A 601 ? 9.53600 9.17400 17.98000 1.000 36.78336 O ? B ? . O01 1
"""

str_ARA_A_601_4zeb = """
data_default
_cell.length_a                    28.645
_cell.length_b                    29.382
_cell.length_c                    24.874
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      20935.137
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                26.26
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   ARA
   TT7

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   HETATM 1 C1 . ARA A 601 ? 16.14500 14.34200 13.30000 1.000 27.53000 C ? A ? . C1 1
   HETATM 2 C2 . ARA A 601 ? 15.36700 15.64100 13.15300 1.000 27.34000 C ? A ? . C2 1
   HETATM 3 C3 . ARA A 601 ? 13.94200 15.39600 13.20100 1.000 26.79000 C ? A ? . C3 1
   HETATM 4 C4 . ARA A 601 ? 13.45300 14.34100 12.37000 1.000 26.91000 C ? A ? . C4 1
   HETATM 5 C5 . ARA A 601 ? 14.24000 13.06700 12.61600 1.000 28.63000 C ? A ? . C5 1
   HETATM 6 O1 . ARA A 601 ? 17.48600 14.50100 13.04200 1.000 28.53000 O ? A ? . O1 1
   HETATM 7 O2 . ARA A 601 ? 15.74300 16.54200 14.23000 1.000 26.19000 O ? A ? . O2 1
   HETATM 8 O3 . ARA A 601 ? 13.16100 16.61900 12.96600 1.000 23.08000 O ? A ? . O3 1
   HETATM 9 O4 . ARA A 601 ? 13.62800 14.63900 10.96400 1.000 25.91000 O ? A ? . O4 1
   HETATM 10 O5 . ARA A 601 ? 15.63900 13.29400 12.48100 1.000 27.40000 O ? A ? . O5 1
   HETATM 119 C1 . TT7 C 1 ? 13.62400 23.30300 14.55900 1.000 38.99000 C ? G ? . C1 1
   HETATM 120 C2 . TT7 C 1 ? 13.84200 21.95500 15.22700 1.000 36.76000 C ? G ? . C2 1
   HETATM 121 C3 . TT7 C 1 ? 14.99800 21.23800 14.56100 1.000 35.09000 C ? G ? . C3 1
   HETATM 122 C4 . TT7 C 1 ? 14.76800 19.89600 14.87200 1.000 31.23000 C ? G ? . C4 1
   HETATM 123 C5 . TT7 C 1 ? 13.29300 19.78000 14.79100 1.000 35.07000 C ? G ? . C5 1
   HETATM 124 C6 . TT7 C 1 ? 12.78700 18.80500 15.79500 1.000 36.08000 C ? G ? . C6 1
   HETATM 125 O1 . TT7 C 1 ? 14.37100 24.38200 15.04400 1.000 40.98000 O ? G ? . O1 1
   HETATM 126 O2 . TT7 C 1 ? 14.23300 22.14800 16.62800 1.000 35.46000 O ? G ? . O2 1
   HETATM 127 O3 . TT7 C 1 ? 16.27200 21.68300 14.97000 1.000 34.75000 O ? G ? . O3 1
   HETATM 128 O4 . TT7 C 1 ? 15.32900 19.06600 13.85200 1.000 24.92000 O ? G ? . O4 1
   HETATM 129 O5 . TT7 C 1 ? 12.77400 21.06000 15.04500 1.000 35.93000 O ? G ? . O5 1
   HETATM 130 O6 . TT7 C 1 ? 11.40000 18.62300 15.76600 1.000 36.67000 O ? G ? . O6 1
   HETATM 131 OP1 . TT7 C 1 ? 17.38700 17.97400 12.87300 1.000 23.50000 O ? G ? . OP1 1
   HETATM 132 OP2 . TT7 C 1 ? 17.25900 18.19300 15.31900 1.000 23.88000 O ? G ? . OP2 1
   HETATM 133 P4 . TT7 C 1 ? 16.47800 17.94500 14.07000 1.000 24.90000 P ? G ? . P4 1
"""

str_BGS_A_961_2b5z = """
data_default
_cell.length_a                    23.935
_cell.length_b                    25.452
_cell.length_c                    27.733
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      16894.767
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                19.46
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   BGS
   LYS

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   HETATM 1 C1 . BGS A 961 ? 12.79100 13.02600 14.85400 1.000 19.40000 C ? A ? . C1 1
   HETATM 2 C2 . BGS A 961 ? 11.84600 12.24700 15.79600 1.000 19.59000 C ? A ? . C2 1
   HETATM 3 C2' . BGS A 961 ? 14.53600 11.26800 13.42800 1.000 20.00000 C ? A ? . C2' 1
   HETATM 4 C3 . BGS A 961 ? 10.45400 12.88300 15.73600 1.000 19.93000 C ? A ? . C3 1
   HETATM 5 C4 . BGS A 961 ? 10.56600 14.36000 16.15400 1.000 20.24000 C ? A ? . C4 1
   HETATM 6 C5 . BGS A 961 ? 11.55100 15.08500 15.20400 1.000 21.11000 C ? A ? . C5 1
   HETATM 7 C6 . BGS A 961 ? 11.74500 16.54600 15.58200 1.000 23.72000 C ? A ? . C6 1
   HETATM 8 CS . BGS A 961 ? 14.74700 12.12600 12.18900 1.000 20.47000 C ? A ? . CS 1
   HETATM 9 O1' . BGS A 961 ? 15.53500 13.39400 14.89600 1.000 21.70000 O ? A ? . O1' 1
   HETATM 10 O2 . BGS A 961 ? 11.78200 10.89900 15.40700 1.000 19.92000 O ? A ? . O2 1
   HETATM 11 O2' . BGS A 961 ? 14.68200 11.52300 16.13400 1.000 20.30000 O ? A ? . O2' 1
   HETATM 12 O3 . BGS A 961 ? 9.58200 12.19900 16.58600 1.000 17.47000 O ? A ? . O3 1
   HETATM 13 O4 . BGS A 961 ? 9.26700 14.99700 16.10500 1.000 21.50000 O ? A ? . O4 1
   HETATM 14 O5 . BGS A 961 ? 12.86000 14.42800 15.23600 1.000 19.23000 O ? A ? . O5 1
   HETATM 15 O6 . BGS A 961 ? 12.69600 16.81300 16.59200 1.000 26.79000 O ? A ? . O6 1
   HETATM 16 S1 . BGS A 961 ? 14.47700 12.35600 14.91800 1.000 19.43000 S ? A ? . S1 1
   ATOM 25 N . LYS A 96 ? 15.93800 18.91900 11.46200 1.000 14.14000 N ? B ? 3 N 1
   ATOM 26 CA . LYS A 96 ? 14.50700 18.82600 11.74300 1.000 14.75000 C ? B ? 3 CA 1
   ATOM 27 C . LYS A 96 ? 14.12400 19.82400 12.82900 1.000 17.03000 C ? B ? 3 C 1
   ATOM 28 O . LYS A 96 ? 13.06700 20.45200 12.77500 1.000 16.73000 O ? B ? 3 O 1
   ATOM 29 CB . LYS A 96 ? 14.13900 17.43400 12.23400 1.000 14.91000 C ? B ? 3 CB 1
   ATOM 30 CG . LYS A 96 ? 14.14800 16.35800 11.19400 1.000 16.40000 C ? B ? 3 CG 1
   ATOM 31 CD . LYS A 96 ? 13.77000 15.05000 11.85500 1.000 17.35000 C ? B ? 3 CD 1
   ATOM 32 CE . LYS A 96 ? 13.69900 13.93600 10.84900 1.000 19.02000 C ? B ? 3 CE 1
   ATOM 33 NZ . LYS A 96 ? 13.54200 12.61000 11.50900 1.000 20.40000 N ? B ? 3 NZ 1
   HETATM 36 C1 . BGS A 962 ? 8.82900 11.73300 10.97700 0.980 25.71000 C ? D ? . C1 1
   HETATM 37 C2 . BGS A 962 ? 7.76400 12.16200 12.01700 0.980 28.35000 C ? D ? . C2 1
   HETATM 38 C2' . BGS A 962 ? 11.57200 11.78500 10.18000 0.980 22.10000 C ? D ? . C2' 1
   HETATM 39 C3 . BGS A 962 ? 6.37300 12.01100 11.39200 0.980 30.04000 C ? D ? . C3 1
   HETATM 40 C4 . BGS A 962 ? 6.30400 12.88600 10.13100 0.980 31.03000 C ? D ? . C4 1
   HETATM 41 C5 . BGS A 962 ? 7.41100 12.44200 9.13600 0.980 31.62000 C ? D ? . C5 1
   HETATM 42 C6 . BGS A 962 ? 7.43300 13.29800 7.88000 0.980 32.48000 C ? D ? . C6 1
   HETATM 43 CS . BGS A 962 ? 12.99600 11.56800 10.64200 0.980 22.06000 C ? D ? . CS 1
   HETATM 44 O1' . BGS A 962 ? 10.69800 13.24800 12.34900 0.980 21.48000 O ? D ? . O1' 1
   HETATM 45 O2 . BGS A 962 ? 7.86300 11.37200 13.17300 0.980 28.61000 O ? D ? . O2 1
   HETATM 46 O2' . BGS A 962 ? 10.79800 10.85700 12.63200 0.980 25.60000 O ? D ? . O2' 1
   HETATM 47 O3 . BGS A 962 ? 5.38900 12.39700 12.30800 0.980 30.87000 O ? D ? . O3 1
   HETATM 48 O4 . BGS A 962 ? 5.00000 12.76900 9.52300 0.980 33.18000 O ? D ? . O4 1
   HETATM 49 O5 . BGS A 962 ? 8.72300 12.54700 9.77100 0.980 28.70000 O ? D ? . O5 1
   HETATM 50 O6 . BGS A 962 ? 8.42300 12.96500 6.92000 0.980 36.56000 O ? D ? . O6 1
   HETATM 51 S1 . BGS A 962 ? 10.48900 11.95400 11.66700 0.980 24.10000 S ? D ? . S1 1
"""

str_VPH_A_304_8oji = """
data_default
_cell.length_a                    21.124
_cell.length_b                    26.009
_cell.length_c                    30.339
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      16668.675
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                19.89
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   VPH
   YIO

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   HETATM 1 C1' . VPH A 304 ? 12.21900 12.90000 17.87400 1.000 27.25000 C ? A ? . C1' 1
   HETATM 2 C2' . VPH A 304 ? 11.48400 12.44800 16.61100 1.000 22.73000 C ? A ? . C2' 1
   HETATM 3 C3' . VPH A 304 ? 12.27400 12.81800 15.37300 1.000 22.83000 C ? A ? . C3' 1
   HETATM 4 C4' . VPH A 304 ? 12.60500 14.30200 15.38400 1.000 24.04000 C ? A ? . C4' 1
   HETATM 5 C5' . VPH A 304 ? 13.29800 14.64900 16.69300 1.000 26.60000 C ? A ? . C5' 1
   HETATM 6 C6' . VPH A 304 ? 13.65900 16.10500 16.82200 1.000 32.60000 C ? A ? . C6' 1
   HETATM 7 C7 . VPH A 304 ? 11.45900 12.53700 19.14200 1.000 31.13000 C ? A ? . C7 1
   HETATM 8 C9 . VPH A 304 ? 9.75600 13.50900 20.43200 1.000 40.81000 C ? A ? . C9 1
   HETATM 9 O3' . VPH A 304 ? 11.52400 12.47400 14.21400 1.000 20.61000 O ? A ? . O3' 1
   HETATM 10 O4' . VPH A 304 ? 11.41300 15.05500 15.26900 1.000 20.95000 O ? A ? . O4' 1
   HETATM 11 O5' . VPH A 304 ? 12.45200 14.30800 17.80100 1.000 29.20000 O ? A ? . O5' 1
   HETATM 12 O6' . VPH A 304 ? 14.68000 16.25400 17.79500 1.000 37.26000 O ? A ? . O6' 1
   HETATM 13 O7 . VPH A 304 ? 11.27100 11.41300 19.51400 1.000 31.75000 O ? A ? . O7 1
   HETATM 14 O8 . VPH A 304 ? 11.03600 13.62400 19.77000 1.000 34.46000 O ? A ? . O8 1
   HETATM 26 C1 . YIO A 303 ? 8.71100 11.76400 17.01800 1.000 17.50000 C ? E ? . C1 1
   HETATM 27 C2 . YIO A 303 ? 7.23600 12.11600 16.87500 1.000 17.74000 C ? E ? . C2 1
   HETATM 28 C3 . YIO A 303 ? 6.37200 10.88200 17.11700 1.000 15.97000 C ? E ? . C3 1
   HETATM 29 C4 . YIO A 303 ? 6.81500 9.74100 16.23400 1.000 16.04000 C ? E ? . C4 1
   HETATM 30 C5 . YIO A 303 ? 8.28600 9.47600 16.50700 1.000 15.71000 C ? E ? . C5 1
   HETATM 31 C6 . YIO A 303 ? 8.87100 8.34000 15.69300 1.000 16.06000 C ? E ? . C6 1
   HETATM 32 O2 . YIO A 303 ? 6.86900 13.10900 17.82500 1.000 19.24000 O ? E ? . O2 1
   HETATM 33 O3 . YIO A 303 ? 5.00000 11.18300 16.89300 1.000 16.51000 O ? E ? . O3 1
   HETATM 34 O4 . YIO A 303 ? 6.63900 10.13500 14.88000 1.000 13.62000 O ? E ? . O4 1
   HETATM 35 O5 . YIO A 303 ? 9.04100 10.65600 16.19700 1.000 16.94000 O ? E ? . O5 1
   HETATM 36 O6 . YIO A 303 ? 10.27400 8.25200 15.87800 1.000 15.91000 O ? E ? . O6 1
   HETATM 37 S1 . YIO A 303 ? 9.79300 13.13600 16.52800 1.000 19.82000 S ? E ? . S1 1
"""

str_Z5L_L_1_2ww3 = """
data_default
_cell.length_a                    26.702
_cell.length_b                    22.647
_cell.length_c                    22.939
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      13871.677
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
_refine.pdbx_stereochemistry_target_values 'GeoStd + Monomer Library + CDL v1.2'
_refine.B_iso_mean                31.07
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _chem_comp.id
   MAN
   Z5L

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   HETATM 1 C1 . Z5L L 1 ? 15.96600 11.15300 12.05900 1.000 35.55000 C ? A ? . C1 1
   HETATM 2 C1M . Z5L L 1 ? 17.33700 9.44000 11.13900 1.000 37.34000 C ? A ? . C1M 1
   HETATM 3 C2 . Z5L L 1 ? 15.44800 11.64300 13.41200 1.000 34.32000 C ? A ? . C2 1
   HETATM 4 C3 . Z5L L 1 ? 14.17400 10.91400 13.84900 1.000 33.45000 C ? A ? . C3 1
   HETATM 5 C4 . Z5L L 1 ? 13.16500 10.85100 12.72300 1.000 33.94000 C ? A ? . C4 1
   HETATM 6 C5 . Z5L L 1 ? 13.80100 10.39700 11.41000 1.000 35.07000 C ? A ? . C5 1
   HETATM 7 C6 . Z5L L 1 ? 12.77700 10.51700 10.28300 1.000 36.69000 C ? A ? . C6 1
   HETATM 8 O1 . Z5L L 1 ? 16.41800 9.79500 12.17900 1.000 37.89000 O ? A ? . O1 1
   HETATM 9 O3 . Z5L L 1 ? 13.55200 11.59600 14.94600 1.000 34.39000 O ? A ? . O3 1
   HETATM 10 O4 . Z5L L 1 ? 12.08900 9.99900 13.11700 1.000 32.39000 O ? A ? . O4 1
   HETATM 11 O5 . Z5L L 1 ? 14.92500 11.22900 11.08700 1.000 35.01000 O ? A ? . O5 1
   HETATM 12 O6 . Z5L L 1 ? 13.39400 10.11000 9.05400 1.000 38.00000 O ? A ? . O6 1
   HETATM 13 S2 . Z5L L 1 ? 15.21200 13.40100 13.36900 1.000 34.73000 S ? A ? . S2 1
   HETATM 14 C1 . MAN L 2 ? 15.46900 13.84900 15.05200 1.000 34.60000 C ? B ? . C1 1
   HETATM 15 C2 . MAN L 2 ? 14.80400 15.20100 15.36300 1.000 33.70000 C ? B ? . C2 1
   HETATM 16 C3 . MAN L 2 ? 15.61100 16.42200 14.93500 1.000 32.27000 C ? B ? . C3 1
   HETATM 17 C4 . MAN L 2 ? 17.11000 16.20900 15.17100 1.000 33.89000 C ? B ? . C4 1
   HETATM 18 C5 . MAN L 2 ? 17.58000 14.85100 14.64400 1.000 34.23000 C ? B ? . C5 1
   HETATM 19 C6 . MAN L 2 ? 19.08800 14.66000 14.80600 1.000 33.66000 C ? B ? . C6 1
   HETATM 20 O2 . MAN L 2 ? 14.55600 15.23000 16.77900 1.000 32.97000 O ? B ? . O2 1
   HETATM 21 O3 . MAN L 2 ? 15.19600 17.59400 15.67300 1.000 31.20000 O ? B ? . O3 1
   HETATM 22 O4 . MAN L 2 ? 17.86300 17.26900 14.54100 1.000 33.04000 O ? B ? . O4 1
   HETATM 23 O5 . MAN L 2 ? 16.88200 13.81000 15.32000 1.000 34.76000 O ? B ? . O5 1
   HETATM 24 O6 . MAN L 2 ? 19.44100 14.53800 16.18600 1.000 33.01000 O ? B ? . O6 1
"""

args = [

  group_args(
    prefix  = "1S6_B_401_4jxg",
    pdb_str = str_1S6_B_401_4jxg,
    sel_A   = "resname SER and name OG",
    sel_B   = "resname 1S6 and name C",
    d_AB    = 1.37487,
    badh    = ["chain B resseq  401 resname 1S6 name  H1"]
  ),

  group_args(
    prefix  = "3W9_A_601_4x0u",
    pdb_str = str_3W9_A_601_4x0u,
    sel_A   = "resname CYS and name SG",
    sel_B   = "resname 3W9 and name C13",
    d_AB    = 1.68035383178,
    badh    = ["chain A resseq  601 resname 3W9 name  H1 "]
  ),

  group_args(
    prefix  = "ARA_A_601_4zeb",
    pdb_str = str_ARA_A_601_4zeb,
    sel_A   = "resname ARA and name O2",
    sel_B   = "resname TT7 and name P4",
    d_AB    = 1.59192776,
    badh    = ["chain A resseq  601 resname ARA name  HO2"]
  ),

  group_args(
    prefix  = "BGS_A_961_2b5z",
    pdb_str = str_BGS_A_961_2b5z,
    sel_A   = "resname LYS and name NZ",
    sel_B   = "resname BGS and resseq 961 and name CS",
    d_AB    = 1.465837985,
    badh    = ["chain A resseq   96 resname LYS name  HZ1"]
  ),

  group_args(
    prefix  = "VPH_A_304_8oji",
    pdb_str = str_VPH_A_304_8oji,
    sel_A   = "resname YIO and name S1",
    sel_B   = "resname VPH and name C2' ",
    d_AB    = 1.82748844,
    badh    = ["chain A resseq  303 resname YIO name  HS1",
               "chain A resseq  304 resname VPH name  H4"]
  ),

  group_args(
    prefix  = "Z5L_L_1_2ww3",
    pdb_str = str_Z5L_L_1_2ww3,
    sel_A   = "resname Z5L and name S2",
    sel_B   = "resname MAN and name C1",
    d_AB    = 1.7604664,
    badh    = ["chain L resseq    1 resname Z5L name  HS2"]
  ),

]

def add_h(model):
  o = reduce_hydrogen.place_hydrogens(
    model                 = model,
    use_neutron_distances = True,
    keep_existing_H       = False,
    validate_e            = False # crworking with altloc
    )
  o.run()
  return o.get_model()

def workaround_003(model):
  pairs = get_ligand_interactions(
    model=model, dist_min=1.1, cutoff_cno=1.6, cutoff_sp=1.9)
  remove_selection = []
  for pair in pairs:
    badH_i_seqs = get_incorrect_hydrogens_for_bond(
      model   = model,
      i_seq_A = pair[0],
      i_seq_B = pair[1])
    if len(badH_i_seqs)>0:
      remove_selection.extend(badH_i_seqs)
  removed = 0
  if len(remove_selection)>0:
    badH_i_seqs = flex.size_t(badH_i_seqs)
    removed = badH_i_seqs.size()
    keep_selection = ~flex.bool(model.size(), badH_i_seqs)
    model = model.select(keep_selection)
  return model, removed

def run_one(prefix, pdb_str, badh):
  #print(prefix)
  fin = "%s.cif"%prefix
  with open(fin, "w") as fo:
    fo.write(pdb_str)
  pdb_inp = iotbx.pdb.input(file_name = fin)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model = add_h(model = model)
  for ss in badh:
    sel = model.selection(string = ss)
    assert sel.count(True)==0

if (__name__ == "__main__"):
  t0 = time.time()
  for g in args:
    run_one(
      prefix  = g.prefix,
      pdb_str = g.pdb_str,
      badh    = g.badh)
  print("OK. Time: %8.3f"%(time.time()-t0))
