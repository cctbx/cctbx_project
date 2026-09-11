from __future__ import division
from __future__ import print_function
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

# Answer strings: reduce2 (placement + optimization) output, checked for
# protonation and geometry. Removing H and re-running reduce2 must reproduce them.

ARA_A_601_4zeb_H = '''
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
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A
   B

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
   HETATM 11 H1 . ARA A 601 ? 16.01224 14.08818 14.22674 1.000 27.53000 H ? A ? . H1 1
   HETATM 12 H2 . ARA A 601 ? 15.60621 15.98923 12.27984 1.000 27.34000 H ? A ? . H2 1
   HETATM 13 H3 . ARA A 601 ? 13.89016 15.09270 14.12090 1.000 26.79000 H ? A ? . H3 1
   HETATM 14 H4 . ARA A 601 ? 12.52884 14.25771 12.65265 1.000 26.91000 H ? A ? . H4 1
   HETATM 15 H51 . ARA A 601 ? 14.03242 12.74092 13.50565 1.000 28.63000 H ? A ? . H51 1
   HETATM 16 H52 . ARA A 601 ? 13.94242 12.39193 11.98621 1.000 28.63000 H ? A ? . H52 1
   HETATM 17 HO1 . ARA A 601 ? 17.68776 15.32443 13.10334 1.000 28.53000 H ? A ? . HO1 1
   HETATM 18 HO3 . ARA A 601 ? 13.61326 17.28780 13.23184 1.000 23.08000 H ? A ? . HO3 1
   HETATM 19 HO4 . ARA A 601 ? 14.44925 14.53969 10.76859 1.000 25.91000 H ? A ? . HO4 1
   HETATM 20 C1 . TT7 C 1 ? 13.62400 23.30300 14.55900 1.000 38.99000 C ? B ? . C1 1
   HETATM 21 C2 . TT7 C 1 ? 13.84200 21.95500 15.22700 1.000 36.76000 C ? B ? . C2 1
   HETATM 22 C3 . TT7 C 1 ? 14.99800 21.23800 14.56100 1.000 35.09000 C ? B ? . C3 1
   HETATM 23 C4 . TT7 C 1 ? 14.76800 19.89600 14.87200 1.000 31.23000 C ? B ? . C4 1
   HETATM 24 C5 . TT7 C 1 ? 13.29300 19.78000 14.79100 1.000 35.07000 C ? B ? . C5 1
   HETATM 25 C6 . TT7 C 1 ? 12.78700 18.80500 15.79500 1.000 36.08000 C ? B ? . C6 1
   HETATM 26 O1 . TT7 C 1 ? 14.37100 24.38200 15.04400 1.000 40.98000 O ? B ? . O1 1
   HETATM 27 O2 . TT7 C 1 ? 14.23300 22.14800 16.62800 1.000 35.46000 O ? B ? . O2 1
   HETATM 28 O3 . TT7 C 1 ? 16.27200 21.68300 14.97000 1.000 34.75000 O ? B ? . O3 1
   HETATM 29 O4 . TT7 C 1 ? 15.32900 19.06600 13.85200 1.000 24.92000 O ? B ? . O4 1
   HETATM 30 O5 . TT7 C 1 ? 12.77400 21.06000 15.04500 1.000 35.93000 O ? B ? . O5 1
   HETATM 31 O6 . TT7 C 1 ? 11.40000 18.62300 15.76600 1.000 36.67000 O ? B ? . O6 1
   HETATM 32 OP1 . TT7 C 1 ? 17.38700 17.97400 12.87300 1.000 23.50000 O ? B ? . OP1 1
   HETATM 33 OP2 . TT7 C 1 ? 17.25900 18.19300 15.31900 1.000 23.88000 O ? B ? . OP2 1
   HETATM 34 P4 . TT7 C 1 ? 16.47800 17.94500 14.07000 1.000 24.90000 P ? B ? . P4 1
   HETATM 35 H3 . TT7 C 1 ? 14.85088 21.37304 13.61178 1.000 35.09000 H ? B ? . H3 1
   HETATM 36 H11 . TT7 C 1 ? 13.80834 23.16591 13.61660 1.000 38.99000 H ? B ? . H11 1
   HETATM 37 H12 . TT7 C 1 ? 12.67566 23.48997 14.64020 1.000 38.99000 H ? B ? . H12 1
   HETATM 38 H4 . TT7 C 1 ? 15.04177 19.71616 15.78502 1.000 31.23000 H ? B ? . H4 1
   HETATM 39 H5 . TT7 C 1 ? 13.03438 19.52574 13.89135 1.000 35.07000 H ? B ? . H5 1
   HETATM 40 H61 . TT7 C 1 ? 13.24650 17.96649 15.63172 1.000 36.08000 H ? B ? . H61 1
   HETATM 41 H62 . TT7 C 1 ? 13.07332 19.12148 16.66607 1.000 36.08000 H ? B ? . H62 1
   HETATM 42 HO1 . TT7 C 1 ? 14.90218 24.08799 15.63890 1.000 40.98000 H ? B ? . HO1 1
   HETATM 43 HO2 . TT7 C 1 ? 15.01912 21.83055 16.68919 1.000 35.46000 H ? B ? . HO2 1
   HETATM 44 HO3 . TT7 C 1 ? 16.17632 22.51958 15.08612 1.000 34.75000 H ? B ? . HO3 1
   HETATM 45 HO6 . TT7 C 1 ? 11.11355 18.88441 15.00962 1.000 36.67000 H ? B ? . HO6 1
'''

Z5L_L_1_2ww3_H = """
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
   HETATM 14 H1 . Z5L L 1 ? 16.69231 11.73360 11.78281 1.000 35.55000 H ? A ? . H1 1
   HETATM 15 H2 . Z5L L 1 ? 16.17705 11.46593 14.02685 1.000 34.32000 H ? A ? . H2 1
   HETATM 16 H3 . Z5L L 1 ? 14.48489 10.02775 14.09149 1.000 33.45000 H ? A ? . H3 1
   HETATM 17 H31 . Z5L L 1 ? 17.59548 8.51669 11.28592 1.000 37.34000 H ? A ? . H31 1
   HETATM 18 H32 . Z5L L 1 ? 18.09689 10.03875 11.20944 1.000 37.34000 H ? A ? . H32 1
   HETATM 19 H33 . Z5L L 1 ? 16.87283 9.55048 10.29446 1.000 37.34000 H ? A ? . H33 1
   HETATM 20 H4 . Z5L L 1 ? 12.84454 11.75776 12.59654 1.000 33.94000 H ? A ? . H4 1
   HETATM 21 H5 . Z5L L 1 ? 14.10149 9.47721 11.47781 1.000 35.07000 H ? A ? . H5 1
   HETATM 22 H61 . Z5L L 1 ? 12.00977 9.96474 10.50041 1.000 36.69000 H ? A ? . H61 1
   HETATM 23 H62 . Z5L L 1 ? 12.46586 11.43485 10.24254 1.000 36.69000 H ? A ? . H62 1
   HETATM 24 HO3 . Z5L L 1 ? 13.17427 12.29153 14.63606 1.000 34.39000 H ? A ? . HO3 1
   HETATM 25 HO4 . Z5L L 1 ? 11.97091 10.07168 13.95561 1.000 32.39000 H ? A ? . HO4 1
   HETATM 26 HO6 . Z5L L 1 ? 13.88400 9.43493 9.21737 1.000 38.00000 H ? A ? . HO6 1
   HETATM 27 C1 . MAN L 2 ? 15.46900 13.84900 15.05200 1.000 34.60000 C ? B ? . C1 1
   HETATM 28 C2 . MAN L 2 ? 14.80400 15.20100 15.36300 1.000 33.70000 C ? B ? . C2 1
   HETATM 29 C3 . MAN L 2 ? 15.61100 16.42200 14.93500 1.000 32.27000 C ? B ? . C3 1
   HETATM 30 C4 . MAN L 2 ? 17.11000 16.20900 15.17100 1.000 33.89000 C ? B ? . C4 1
   HETATM 31 C5 . MAN L 2 ? 17.58000 14.85100 14.64400 1.000 34.23000 C ? B ? . C5 1
   HETATM 32 C6 . MAN L 2 ? 19.08800 14.66000 14.80600 1.000 33.66000 C ? B ? . C6 1
   HETATM 33 O2 . MAN L 2 ? 14.55600 15.23000 16.77900 1.000 32.97000 O ? B ? . O2 1
   HETATM 34 O3 . MAN L 2 ? 15.19600 17.59400 15.67300 1.000 31.20000 O ? B ? . O3 1
   HETATM 35 O4 . MAN L 2 ? 17.86300 17.26900 14.54100 1.000 33.04000 O ? B ? . O4 1
   HETATM 36 O5 . MAN L 2 ? 16.88200 13.81000 15.32000 1.000 34.76000 O ? B ? . O5 1
   HETATM 37 O6 . MAN L 2 ? 19.44100 14.53800 16.18600 1.000 33.01000 O ? B ? . O6 1
   HETATM 38 H1 . MAN L 2 ? 15.02283 13.16624 15.57704 1.000 34.60000 H ? B ? . H1 1
   HETATM 39 H2 . MAN L 2 ? 13.99464 15.19794 14.82837 1.000 33.70000 H ? B ? . H2 1
   HETATM 40 H3 . MAN L 2 ? 15.42002 16.51134 13.98819 1.000 32.27000 H ? B ? . H3 1
   HETATM 41 H4 . MAN L 2 ? 17.20847 16.23661 16.13559 1.000 33.89000 H ? B ? . H4 1
   HETATM 42 H5 . MAN L 2 ? 17.40102 14.79465 13.69232 1.000 34.23000 H ? B ? . H5 1
   HETATM 43 H61 . MAN L 2 ? 19.53392 15.41688 14.39468 1.000 33.66000 H ? B ? . H61 1
   HETATM 44 H62 . MAN L 2 ? 19.34764 13.87062 14.30563 1.000 33.66000 H ? B ? . H62 1
   HETATM 45 HO2 . MAN L 2 ? 15.29353 15.39124 17.16958 1.000 32.97000 H ? B ? . HO2 1
   HETATM 46 HO3 . MAN L 2 ? 15.12384 17.37728 16.49174 1.000 31.20000 H ? B ? . HO3 1
   HETATM 47 HO4 . MAN L 2 ? 17.99491 17.07411 13.72423 1.000 33.04000 H ? B ? . HO4 1
   HETATM 48 HO6 . MAN L 2 ? 18.76323 14.24548 16.60738 1.000 33.01000 H ? B ? . HO6 1
"""

VPH_A_304_8oji_H = """
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
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A
   B
   C
   D
   E
   F
   G
   H

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
   HETATM 15 H1 . VPH A 304 ? 14.12396 14.14184 16.73148 1.000 26.60000 H ? A ? . H1 1
   HETATM 16 H2 . VPH A 304 ? 13.20474 14.46537 14.63934 1.000 24.04000 H ? A ? . H2 1
   HETATM 17 H3 . VPH A 304 ? 13.13736 12.37610 15.35751 1.000 22.83000 H ? A ? . H3 1
   HETATM 18 H10 . VPH A 304 ? 10.76430 14.58367 15.55103 1.000 20.95000 H ? A ? . H10 1
   HETATM 19 H11 . VPH A 304 ? 10.90934 13.05297 14.11662 1.000 20.61000 H ? A ? . H11 1
   HETATM 20 H12 . VPH A 304 ? 9.79007 12.81044 21.10413 1.000 40.81000 H ? A ? . H12 1
   HETATM 21 H13 . VPH A 304 ? 9.06571 13.28963 19.78680 1.000 40.81000 H ? A ? . H13 1
   HETATM 22 H14 . VPH A 304 ? 9.52717 14.34742 20.86281 1.000 40.81000 H ? A ? . H14 1
   HETATM 23 H5 . VPH A 304 ? 11.37411 11.48469 16.64017 1.000 22.73000 H ? A ? . H5 1
   HETATM 24 H6 . VPH A 304 ? 13.04509 12.39826 17.95602 1.000 27.25000 H ? A ? . H6 1
   HETATM 25 H7 . VPH A 304 ? 13.95739 16.42710 15.95706 1.000 32.60000 H ? A ? . H7 1
   HETATM 26 H8 . VPH A 304 ? 12.86427 16.59883 17.07781 1.000 32.60000 H ? A ? . H8 1
   HETATM 27 H9 . VPH A 304 ? 15.06223 15.50139 17.89490 1.000 37.26000 H ? A ? . H9 1
   HETATM 28 C1 . YIO A 303 ? 8.71100 11.76400 17.01800 1.000 17.50000 C ? E ? . C1 1
   HETATM 29 C2 . YIO A 303 ? 7.23600 12.11600 16.87500 1.000 17.74000 C ? E ? . C2 1
   HETATM 30 C3 . YIO A 303 ? 6.37200 10.88200 17.11700 1.000 15.97000 C ? E ? . C3 1
   HETATM 31 C4 . YIO A 303 ? 6.81500 9.74100 16.23400 1.000 16.04000 C ? E ? . C4 1
   HETATM 32 C5 . YIO A 303 ? 8.28600 9.47600 16.50700 1.000 15.71000 C ? E ? . C5 1
   HETATM 33 C6 . YIO A 303 ? 8.87100 8.34000 15.69300 1.000 16.06000 C ? E ? . C6 1
   HETATM 34 O2 . YIO A 303 ? 6.86900 13.10900 17.82500 1.000 19.24000 O ? E ? . O2 1
   HETATM 35 O3 . YIO A 303 ? 5.00000 11.18300 16.89300 1.000 16.51000 O ? E ? . O3 1
   HETATM 36 O4 . YIO A 303 ? 6.63900 10.13500 14.88000 1.000 13.62000 O ? E ? . O4 1
   HETATM 37 O5 . YIO A 303 ? 9.04100 10.65600 16.19700 1.000 16.94000 O ? E ? . O5 1
   HETATM 38 O6 . YIO A 303 ? 10.27400 8.25200 15.87800 1.000 15.91000 O ? E ? . O6 1
   HETATM 39 S1 . YIO A 303 ? 9.79300 13.13600 16.52800 1.000 19.82000 S ? E ? . S1 1
   HETATM 40 H1 . YIO A 303 ? 8.87346 11.61738 17.96299 1.000 17.50000 H ? E ? . H1 1
   HETATM 41 H2 . YIO A 303 ? 7.06728 12.40088 15.96326 1.000 17.74000 H ? E ? . H2 1
   HETATM 42 H3 . YIO A 303 ? 6.50194 10.56132 18.02319 1.000 15.97000 H ? E ? . H3 1
   HETATM 43 HA . YIO A 303 ? 7.53275 13.62455 17.95204 1.000 19.24000 H ? E ? . HA 1
   HETATM 44 HB . YIO A 303 ? 4.86472 11.23625 16.05553 1.000 16.51000 H ? E ? . HB 1
   HETATM 45 H4 . YIO A 303 ? 6.34716 8.91935 16.45061 1.000 16.04000 H ? E ? . H4 1
   HETATM 46 H5 . YIO A 303 ? 8.35945 9.22873 17.44207 1.000 15.71000 H ? E ? . H5 1
   HETATM 47 HC . YIO A 303 ? 7.32214 10.58109 14.64163 1.000 13.62000 H ? E ? . HC 1
   HETATM 48 HO6 . YIO A 303 ? 10.56476 9.02802 16.06706 1.000 15.91000 H ? E ? . HO6 1
   HETATM 49 H61C . YIO A 303 ? 8.42067 7.52493 15.96460 1.000 16.06000 H ? E ? . H61C 1
   HETATM 50 H62C . YIO A 303 ? 8.63999 8.49797 14.76425 1.000 16.06000 H ? E ? . H62C 1
"""

BGS_A_961_2b5z_H = """
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
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A
   B
   C
   D
   E
   F
   G
   H
   I
   J

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
   HETATM 17 H1 . BGS A 961 ? 12.46274 12.92207 13.94717 1.000 19.40000 H ? A ? . H1 1
   HETATM 18 H2 . BGS A 961 ? 12.18514 12.25723 16.70472 1.000 19.59000 H ? A ? . H2 1
   HETATM 19 H3 . BGS A 961 ? 10.07672 12.81137 14.84525 1.000 19.93000 H ? A ? . H3 1
   HETATM 20 H4 . BGS A 961 ? 10.87230 14.42813 17.07185 1.000 20.24000 H ? A ? . H4 1
   HETATM 21 H5 . BGS A 961 ? 11.16199 15.07631 14.31546 1.000 21.11000 H ? A ? . H5 1
   HETATM 22 H61 . BGS A 961 ? 12.02002 17.02113 14.78231 1.000 23.72000 H ? A ? . H61 1
   HETATM 23 H62 . BGS A 961 ? 10.88825 16.88649 15.88356 1.000 23.72000 H ? A ? . H62 1
   HETATM 24 HO2 . BGS A 961 ? 11.46124 10.84553 14.62166 1.000 19.92000 H ? A ? . HO2 1
   HETATM 25 HO3 . BGS A 961 ? 8.99498 12.73806 16.88150 1.000 17.47000 H ? A ? . HO3 1
   HETATM 26 HO4 . BGS A 961 ? 8.71364 14.48189 15.71646 1.000 21.50000 H ? A ? . HO4 1
   HETATM 27 HO6 . BGS A 961 ? 12.65513 16.20649 17.18612 1.000 26.79000 H ? A ? . HO6 1
   HETATM 28 HS1 . BGS A 961 ? 15.25348 11.60016 11.55036 1.000 20.47000 H ? A ? . HS1 1
   HETATM 29 HS2 . BGS A 961 ? 15.26639 12.90132 12.45360 1.000 20.47000 H ? A ? . HS2 1
   HETATM 30 H2'1 . BGS A 961 ? 13.69991 10.78297 13.34670 1.000 20.00000 H ? A ? . H2'1 1
   HETATM 31 H2'2 . BGS A 961 ? 15.26921 10.63880 13.51411 1.000 20.00000 H ? A ? . H2'2 1
   ATOM 32 N . LYS A 96 ? 15.93800 18.91900 11.46200 1.000 14.14000 N ? B ? 3 N 1
   ATOM 33 CA . LYS A 96 ? 14.50700 18.82600 11.74300 1.000 14.75000 C ? B ? 3 CA 1
   ATOM 34 C . LYS A 96 ? 14.12400 19.82400 12.82900 1.000 17.03000 C ? B ? 3 C 1
   ATOM 35 O . LYS A 96 ? 13.06700 20.45200 12.77500 1.000 16.73000 O ? B ? 3 O 1
   ATOM 36 CB . LYS A 96 ? 14.13900 17.43400 12.23400 1.000 14.91000 C ? B ? 3 CB 1
   ATOM 37 CG . LYS A 96 ? 14.14800 16.35800 11.19400 1.000 16.40000 C ? B ? 3 CG 1
   ATOM 38 CD . LYS A 96 ? 13.77000 15.05000 11.85500 1.000 17.35000 C ? B ? 3 CD 1
   ATOM 39 CE . LYS A 96 ? 13.69900 13.93600 10.84900 1.000 19.02000 C ? B ? 3 CE 1
   ATOM 40 NZ . LYS A 96 ? 13.54200 12.61000 11.50900 1.000 20.40000 N ? B ? 3 NZ 1
   ATOM 41 H . LYS A 96 ? 16.22384 18.33047 10.90385 1.000 14.14000 H ? B ? 3 H 1
   ATOM 42 HA . LYS A 96 ? 14.01705 19.01093 10.92651 1.000 14.75000 H ? B ? 3 HA 1
   ATOM 43 HB2 . LYS A 96 ? 14.77303 17.17591 12.92123 1.000 14.91000 H ? B ? 3 HB2 1
   ATOM 44 HB3 . LYS A 96 ? 13.24344 17.46875 12.60502 1.000 14.91000 H ? B ? 3 HB3 1
   ATOM 45 HG2 . LYS A 96 ? 13.50217 16.56109 10.49933 1.000 16.40000 H ? B ? 3 HG2 1
   ATOM 46 HG3 . LYS A 96 ? 15.03436 16.27472 10.80887 1.000 16.40000 H ? B ? 3 HG3 1
   ATOM 47 HD2 . LYS A 96 ? 14.43718 14.81997 12.52047 1.000 17.35000 H ? B ? 3 HD2 1
   ATOM 48 HD3 . LYS A 96 ? 12.90011 15.13994 12.27464 1.000 17.35000 H ? B ? 3 HD3 1
   ATOM 49 HE2 . LYS A 96 ? 12.93676 14.07828 10.26620 1.000 19.02000 H ? B ? 3 HE2 1
   ATOM 50 HE3 . LYS A 96 ? 14.51700 13.92188 10.32787 1.000 19.02000 H ? B ? 3 HE3 1
   HETATM 51 C1 . BGS A 962 ? 8.82900 11.73300 10.97700 0.980 25.71000 C ? D ? . C1 1
   HETATM 52 C2 . BGS A 962 ? 7.76400 12.16200 12.01700 0.980 28.35000 C ? D ? . C2 1
   HETATM 53 C2' . BGS A 962 ? 11.57200 11.78500 10.18000 0.980 22.10000 C ? D ? . C2' 1
   HETATM 54 C3 . BGS A 962 ? 6.37300 12.01100 11.39200 0.980 30.04000 C ? D ? . C3 1
   HETATM 55 C4 . BGS A 962 ? 6.30400 12.88600 10.13100 0.980 31.03000 C ? D ? . C4 1
   HETATM 56 C5 . BGS A 962 ? 7.41100 12.44200 9.13600 0.980 31.62000 C ? D ? . C5 1
   HETATM 57 C6 . BGS A 962 ? 7.43300 13.29800 7.88000 0.980 32.48000 C ? D ? . C6 1
   HETATM 58 CS . BGS A 962 ? 12.99600 11.56800 10.64200 0.980 22.06000 C ? D ? . CS 1
   HETATM 59 O1' . BGS A 962 ? 10.69800 13.24800 12.34900 0.980 21.48000 O ? D ? . O1' 1
   HETATM 60 O2 . BGS A 962 ? 7.86300 11.37200 13.17300 0.980 28.61000 O ? D ? . O2 1
   HETATM 61 O2' . BGS A 962 ? 10.79800 10.85700 12.63200 0.980 25.60000 O ? D ? . O2' 1
   HETATM 62 O3 . BGS A 962 ? 5.38900 12.39700 12.30800 0.980 30.87000 O ? D ? . O3 1
   HETATM 63 O4 . BGS A 962 ? 5.00000 12.76900 9.52300 0.980 33.18000 O ? D ? . O4 1
   HETATM 64 O5 . BGS A 962 ? 8.72300 12.54700 9.77100 0.980 28.70000 O ? D ? . O5 1
   HETATM 65 O6 . BGS A 962 ? 8.42300 12.96500 6.92000 0.980 36.56000 O ? D ? . O6 1
   HETATM 66 S1 . BGS A 962 ? 10.48900 11.95400 11.66700 0.980 24.10000 S ? D ? . S1 1
   HETATM 67 H1 . BGS A 962 ? 8.71012 10.79328 10.76795 0.980 25.71000 H ? D ? . H1 1
   HETATM 68 H2 . BGS A 962 ? 7.91604 13.07819 12.29697 0.980 28.35000 H ? D ? . H2 1
   HETATM 69 H3 . BGS A 962 ? 6.18938 11.08358 11.17501 0.980 30.04000 H ? D ? . H3 1
   HETATM 70 H4 . BGS A 962 ? 6.42171 13.82207 10.35642 0.980 31.03000 H ? D ? . H4 1
   HETATM 71 H5 . BGS A 962 ? 7.20826 11.53590 8.85534 0.980 31.62000 H ? D ? . H5 1
   HETATM 72 H61 . BGS A 962 ? 6.56511 13.21879 7.45410 0.980 32.48000 H ? D ? . H61 1
   HETATM 73 H62 . BGS A 962 ? 7.57830 14.21752 8.15251 0.980 32.48000 H ? D ? . H62 1
   HETATM 74 HO2 . BGS A 962 ? 8.66499 11.10658 13.26711 0.980 28.61000 H ? D ? . HO2 1
   HETATM 75 HO3 . BGS A 962 ? 4.74585 12.76622 11.89263 0.980 30.87000 H ? D ? . HO3 1
   HETATM 76 HO4 . BGS A 962 ? 4.98819 12.09426 9.00619 0.980 33.18000 H ? D ? . HO4 1
   HETATM 77 HO6 . BGS A 962 ? 8.68825 12.16835 7.05224 0.980 36.56000 H ? D ? . HO6 1
   HETATM 78 HS1 . BGS A 962 ? 13.02817 10.72855 11.12697 0.980 22.06000 H ? D ? . HS1 1
   HETATM 79 HS2 . BGS A 962 ? 13.55973 11.50846 9.85488 0.980 22.06000 H ? D ? . HS2 1
   HETATM 80 H2'1 . BGS A 962 ? 11.52256 12.59304 9.64565 0.980 22.10000 H ? D ? . H2'1 1
   HETATM 81 H2'2 . BGS A 962 ? 11.28275 11.02652 9.64902 0.980 22.10000 H ? D ? . H2'2 1
"""

SER_1S6_B_401_4jxg_H = """
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
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z


loop_
  _struct_asym.id
   A
   B


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
   ATOM 1 CA . SER B 64 ? 8.81200 8.39800 16.90000 1.000 15.06000 C ? A ? 1 CA 1
   ATOM 2 C . SER B 64 ? 7.91100 7.42700 16.14500 1.000 14.26000 C ? A ? 1 C 1
   ATOM 3 O . SER B 64 ? 6.73200 7.73500 16.00400 1.000 14.43000 O ? A ? 1 O 1
   ATOM 4 CB . SER B 64 ? 8.99800 9.60700 16.06900 1.000 14.19000 C ? A ? 1 CB 1
   ATOM 5 OG . SER B 64 ? 9.58500 10.62900 16.88900 1.000 18.69000 O ? A ? 1 OG 1
   ATOM 6 HA . SER B 64 ? 8.39678 8.64128 17.74220 1.000 15.06000 H ? A ? 1 HA 1
   ATOM 7 HB2 . SER B 64 ? 9.58884 9.40492 15.32672 1.000 14.19000 H ? A ? 1 HB2 1
   ATOM 8 HB3 . SER B 64 ? 8.13857 9.90943 15.73611 1.000 14.19000 H ? A ? 1 HB3 1
   HETATM 9 N . 1S6 B 401 ? 12.79500 11.93900 17.59000 1.000 18.86000 N ? B ? . N 1
   HETATM 10 CA . 1S6 B 401 ? 11.42400 12.12000 17.09200 1.000 19.06000 C ? B ? . CA 1
   HETATM 11 C . 1S6 B 401 ? 10.94500 10.77700 16.75200 1.000 20.57000 C ? B ? . C 1
   HETATM 12 CB . 1S6 B 401 ? 11.11600 12.95200 15.88600 1.000 20.89000 C ? B ? . CB 1
   HETATM 13 OXT . 1S6 B 401 ? 11.61000 9.86800 16.33900 1.000 18.23000 O ? B ? . OXT 1
   HETATM 14 CAA . 1S6 B 401 ? 14.03300 11.13000 21.26800 1.000 32.84000 C ? B ? . CAA 1
   HETATM 15 CAB . 1S6 B 401 ? 14.15300 12.81400 13.28200 1.000 25.99000 C ? B ? . CAB 1
   HETATM 16 CAC . 1S6 B 401 ? 12.14600 13.93700 12.44400 1.000 25.23000 C ? B ? . CAC 1
   HETATM 17 CAI . 1S6 B 401 ? 16.74200 14.85600 14.79300 1.000 40.47000 C ? B ? . CAI 1
   HETATM 18 CAJ . 1S6 B 401 ? 17.25700 13.55800 14.88000 1.000 42.42000 C ? B ? . CAJ 1
   HETATM 19 CAK . 1S6 B 401 ? 15.94700 15.38000 15.83700 1.000 37.93000 C ? B ? . CAK 1
   HETATM 20 CAL . 1S6 B 401 ? 16.99500 12.74200 16.00800 1.000 40.13000 C ? B ? . CAL 1
   HETATM 21 CAM . 1S6 B 401 ? 15.73300 14.56600 16.96500 1.000 33.61000 C ? B ? . CAM 1
   HETATM 22 CAT . 1S6 B 401 ? 11.40700 11.19200 12.64000 1.000 27.95000 C ? B ? . CAT 1
   HETATM 23 CAU . 1S6 B 401 ? 13.30600 12.76400 18.43000 1.000 25.25000 C ? B ? . CAU 1
   HETATM 24 CAV . 1S6 B 401 ? 14.85300 11.69600 20.19200 1.000 29.93000 C ? B ? . CAV 1
   HETATM 25 CAW . 1S6 B 401 ? 16.22500 13.25900 17.04000 1.000 34.14000 C ? B ? . CAW 1
   HETATM 26 CAX . 1S6 B 401 ? 15.94500 12.52600 18.25000 1.000 30.72000 C ? B ? . CAX 1
   HETATM 27 CAY . 1S6 B 401 ? 14.64300 12.41800 18.95200 1.000 26.52000 C ? B ? . CAY 1
   HETATM 28 CBB . 1S6 B 401 ? 11.93000 11.82800 13.87300 1.000 23.23000 C ? B ? . CBB 1
   HETATM 29 CBC . 1S6 B 401 ? 12.69300 13.18100 13.63400 1.000 25.84000 C ? B ? . CBC 1
   HETATM 30 NAN . 1S6 B 401 ? 16.89500 11.99800 19.11900 1.000 28.67000 N ? B ? . NAN 1
   HETATM 31 NAP . 1S6 B 401 ? 10.81600 12.16100 14.72600 1.000 22.89000 N ? B ? . NAP 1
   HETATM 32 OAE . 1S6 B 401 ? 12.15000 10.48600 11.98300 1.000 29.05000 O ? B ? . OAE 1
   HETATM 33 OAF . 1S6 B 401 ? 12.77300 13.75600 18.76700 1.000 24.58000 O ? B ? . OAF 1
   HETATM 34 OAH . 1S6 B 401 ? 10.23400 11.45800 12.25800 1.000 26.16000 O ? B ? . OAH 1
   HETATM 35 OAQ . 1S6 B 401 ? 16.23400 11.51900 20.23100 1.000 35.16000 O ? B ? . OAQ 1
   HETATM 36 SAR . 1S6 B 401 ? 12.56500 13.90400 15.25800 1.000 25.33000 S ? B ? . SAR 1
   HETATM 37 H2 . 1S6 B 401 ? 13.26111 11.27083 17.31450 1.000 18.86000 H ? B ? . H2 1
   HETATM 38 H3 . 1S6 B 401 ? 10.99886 12.62089 17.80562 1.000 19.06000 H ? B ? . H3 1
   HETATM 39 H10 . 1S6 B 401 ? 14.73583 13.58930 13.29238 1.000 25.99000 H ? B ? . H10 1
   HETATM 40 H11 . 1S6 B 401 ? 12.28065 13.43222 11.62671 1.000 25.23000 H ? B ? . H11 1
   HETATM 41 H12 . 1S6 B 401 ? 12.59092 14.79370 12.34903 1.000 25.23000 H ? B ? . H12 1
   HETATM 42 H13 . 1S6 B 401 ? 11.19522 14.09813 12.54869 1.000 25.23000 H ? B ? . H13 1
   HETATM 43 H15 . 1S6 B 401 ? 16.92349 15.37715 14.04442 1.000 40.47000 H ? B ? . H15 1
   HETATM 44 H16 . 1S6 B 401 ? 17.78009 13.22191 14.18839 1.000 42.42000 H ? B ? . H16 1
   HETATM 45 H17 . 1S6 B 401 ? 15.57676 16.23151 15.78447 1.000 37.93000 H ? B ? . H17 1
   HETATM 46 H18 . 1S6 B 401 ? 17.33651 11.87787 16.04748 1.000 40.13000 H ? B ? . H18 1
   HETATM 47 H19 . 1S6 B 401 ? 15.24962 14.90823 17.68203 1.000 33.61000 H ? B ? . H19 1
   HETATM 48 H20 . 1S6 B 401 ? 10.46662 11.38623 14.99014 1.000 22.89000 H ? B ? . H20 1
   HETATM 49 H22 . 1S6 B 401 ? 12.55077 11.21208 14.29274 1.000 23.23000 H ? B ? . H22 1
   HETATM 50 H4 . 1S6 B 401 ? 10.41422 13.55999 16.16664 1.000 20.89000 H ? B ? . H4 1
   HETATM 51 H5 . 1S6 B 401 ? 13.92351 11.76473 21.99328 1.000 32.84000 H ? B ? . H5 1
   HETATM 52 H6 . 1S6 B 401 ? 14.44127 10.32929 21.63279 1.000 32.84000 H ? B ? . H6 1
   HETATM 53 H7 . 1S6 B 401 ? 13.14977 10.89257 20.94485 1.000 32.84000 H ? B ? . H7 1
   HETATM 54 H8 . 1S6 B 401 ? 14.22646 12.42015 12.39860 1.000 25.99000 H ? B ? . H8 1
   HETATM 55 H9 . 1S6 B 401 ? 14.52793 12.17332 13.90639 1.000 25.99000 H ? B ? . H9 1
"""

CYS_1E8_A_701_5p9j_H = """
data_default
_cell.length_a                    20.405
_cell.length_b                    18.624
_cell.length_c                    24.173
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A
   B

loop_
  _chem_comp.id
   1E8
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
   ATOM 1 N . CYS A 481 ? 20.64500 10.93100 -2.08000 1.000 12.93000 N ? A ? 1 N 1
   ATOM 2 CA . CYS A 481 ? 20.49100 12.18700 -1.34600 1.000 13.38000 C ? A ? 1 CA 1
   ATOM 3 C . CYS A 481 ? 20.15500 13.39300 -2.28300 1.000 12.24000 C ? A ? 1 C 1
   ATOM 4 O . CYS A 481 ? 20.70800 13.48400 -3.37800 1.000 12.28000 O ? A ? 1 O 1
   ATOM 5 CB . CYS A 481 ? 21.73800 12.53100 -0.60400 1.000 16.61000 C ? A ? 1 CB 1
   ATOM 6 SG . CYS A 481 ? 23.18400 12.84300 -1.46700 1.000 21.18000 S ? A ? 1 SG 1
   ATOM 7 H . CYS A 481 ? 21.38034 10.51695 -1.91435 1.000 12.93000 H ? A ? 1 H 1
   ATOM 8 HA . CYS A 481 ? 19.75846 12.05695 -0.72361 1.000 13.38000 H ? A ? 1 HA 1
   ATOM 9 HB2 . CYS A 481 ? 21.55411 13.33106 -0.08729 1.000 16.61000 H ? A ? 1 HB2 1
   ATOM 10 HB3 . CYS A 481 ? 21.93307 11.78971 -0.00958 1.000 16.61000 H ? A ? 1 HB3 1
   HETATM 11 C2 . 1E8 A 701 ? 19.02900 5.21000 -0.50200 1.000 13.21000 C ? B ? . C2 1
   HETATM 12 C4 . 1E8 A 701 ? 20.28700 6.08700 1.10600 1.000 13.59000 C ? B ? . C4 1
   HETATM 13 C5 . 1E8 A 701 ? 19.24400 5.77000 2.04500 1.000 13.49000 C ? B ? . C5 1
   HETATM 14 C6 . 1E8 A 701 ? 18.03500 5.20200 1.58300 1.000 12.52000 C ? B ? . C6 1
   HETATM 15 CAA . 1E8 A 701 ? 24.63600 11.91400 -0.79600 1.000 20.77000 C ? B ? . CAA 1
   HETATM 16 CAD . 1E8 A 701 ? 24.92600 10.40100 -1.22500 1.000 21.68000 C ? B ? . CAD 1
   HETATM 17 CAE . 1E8 A 701 ? 14.68100 8.43700 10.02900 1.000 18.37000 C ? B ? . CAE 1
   HETATM 18 CAF . 1E8 A 701 ? 14.53900 8.01300 8.70900 1.000 17.58000 C ? B ? . CAF 1
   HETATM 19 CAG . 1E8 A 701 ? 15.73100 7.99400 10.79500 1.000 20.38000 C ? B ? . CAG 1
   HETATM 20 CAI . 1E8 A 701 ? 15.45800 7.20100 8.14300 1.000 17.54000 C ? B ? . CAI 1
   HETATM 21 CAJ . 1E8 A 701 ? 16.63400 7.08200 10.22600 1.000 18.02000 C ? B ? . CAJ 1
   HETATM 22 CAK . 1E8 A 701 ? 18.22700 4.86000 6.45400 1.000 17.33000 C ? B ? . CAK 1
   HETATM 23 CAL . 1E8 A 701 ? 18.40300 7.30100 6.65500 1.000 17.38000 C ? B ? . CAL 1
   HETATM 24 CAM . 1E8 A 701 ? 18.76600 4.95500 5.19700 1.000 16.86000 C ? B ? . CAM 1
   HETATM 25 CAN . 1E8 A 701 ? 19.00400 7.37700 5.43000 1.000 17.48000 C ? B ? . CAN 1
   HETATM 26 CAO . 1E8 A 701 ? 24.94400 7.93300 1.76700 1.000 23.32000 C ? B ? . CAO 1
   HETATM 27 CAP . 1E8 A 701 ? 23.79100 6.99400 2.16100 1.000 22.36000 C ? B ? . CAP 1
   HETATM 28 CAQ . 1E8 A 701 ? 24.37500 9.39500 1.60100 1.000 21.53000 C ? B ? . CAQ 1
   HETATM 29 CAR . 1E8 A 701 ? 22.23400 8.54000 0.82600 1.000 21.89000 C ? B ? . CAR 1
   HETATM 30 CAW . 1E8 A 701 ? 23.70200 9.56100 -0.81100 1.000 21.67000 C ? B ? . CAW 1
   HETATM 31 CAY . 1E8 A 701 ? 16.53900 6.77700 8.91700 1.000 16.54000 C ? B ? . CAY 1
   HETATM 32 CAZ . 1E8 A 701 ? 18.00500 6.03200 7.18600 1.000 17.48000 C ? B ? . CAZ 1
   HETATM 33 CBA . 1E8 A 701 ? 19.10900 6.22100 4.67500 1.000 15.69000 C ? B ? . CBA 1
   HETATM 34 CBB . 1E8 A 701 ? 19.75700 6.27700 3.36000 1.000 14.48000 C ? B ? . CBB 1
   HETATM 35 CBE . 1E8 A 701 ? 22.53200 7.12700 1.26100 1.000 17.72000 C ? B ? . CBE 1
   HETATM 36 N1 . 1E8 A 701 ? 17.99600 4.90600 0.25100 1.000 12.40000 N ? B ? . N1 1
   HETATM 37 N3 . 1E8 A 701 ? 20.16500 5.77400 -0.20900 1.000 13.78000 N ? B ? . N3 1
   HETATM 38 NAB . 1E8 A 701 ? 16.99300 4.96200 2.34400 1.000 12.25000 N ? B ? . NAB 1
   HETATM 39 NAU . 1E8 A 701 ? 20.95600 6.78500 3.15200 1.000 15.16000 N ? B ? . NAU 1
   HETATM 40 NBF . 1E8 A 701 ? 21.24000 6.66000 1.80100 1.000 15.50000 N ? B ? . NBF 1
   HETATM 41 NBG . 1E8 A 701 ? 23.47700 9.23700 0.49700 1.000 21.51000 N ? B ? . NBG 1
   HETATM 42 OAC . 1E8 A 701 ? 22.89600 9.27100 -1.71900 1.000 19.14000 O ? B ? . OAC 1
   HETATM 43 OAV . 1E8 A 701 ? 17.47300 5.89900 8.44600 1.000 18.42000 O ? B ? . OAV 1
   HETATM 44 H2 . 1E8 A 701 ? 18.92474 4.98035 -1.39715 1.000 13.21000 H ? B ? . H2 1
   HETATM 45 HAA . 1E8 A 701 ? 24.54193 11.91159 0.12923 1.000 20.77000 H ? B ? . HAA 1
   HETATM 46 HAD . 1E8 A 701 ? 25.69512 10.09215 -1.64686 1.000 21.68000 H ? B ? . HAD 1
   HETATM 47 HAE . 1E8 A 701 ? 14.05938 9.02476 10.39373 1.000 18.37000 H ? B ? . HAE 1
   HETATM 48 HAF . 1E8 A 701 ? 13.80547 8.28939 8.20856 1.000 17.58000 H ? B ? . HAF 1
   HETATM 49 HAG . 1E8 A 701 ? 15.84523 8.28814 11.66983 1.000 20.38000 H ? B ? . HAG 1
   HETATM 50 HAI . 1E8 A 701 ? 15.37420 6.92935 7.25752 1.000 17.54000 H ? B ? . HAI 1
   HETATM 51 HAJ . 1E8 A 701 ? 17.29672 6.68712 10.74540 1.000 18.02000 H ? B ? . HAJ 1
   HETATM 52 HAK . 1E8 A 701 ? 18.01806 4.02225 6.79957 1.000 17.33000 H ? B ? . HAK 1
   HETATM 53 HAL . 1E8 A 701 ? 18.25037 8.07442 7.14839 1.000 17.38000 H ? B ? . HAL 1
   HETATM 54 HAM . 1E8 A 701 ? 18.90046 4.18084 4.69951 1.000 16.86000 H ? B ? . HAM 1
   HETATM 55 HAN . 1E8 A 701 ? 19.33341 8.18735 5.11421 1.000 17.48000 H ? B ? . HAN 1
   HETATM 56 HAO . 1E8 A 701 ? 25.63758 7.92235 2.44504 1.000 23.32000 H ? B ? . HAO 1
   HETATM 57 HAP . 1E8 A 701 ? 24.10072 6.07559 2.12228 1.000 22.36000 H ? B ? . HAP 1
   HETATM 58 HAQ . 1E8 A 701 ? 23.90860 9.70947 2.39124 1.000 21.53000 H ? B ? . HAQ 1
   HETATM 59 HAR . 1E8 A 701 ? 21.65458 8.53129 0.04812 1.000 21.89000 H ? B ? . HAR 1
   HETATM 60 HBE . 1E8 A 701 ? 22.77149 6.53823 0.52827 1.000 17.72000 H ? B ? . HBE 1
   HETATM 61 HAAA . 1E8 A 701 ? 25.39181 12.40303 -1.02946 1.000 20.77000 H ? B ? . HAAA 1
   HETATM 62 HAOA . 1E8 A 701 ? 25.35221 7.64009 0.93726 1.000 23.32000 H ? B ? . HAOA 1
   HETATM 63 HAPA . 1E8 A 701 ? 23.53309 7.17385 3.07863 1.000 22.36000 H ? B ? . HAPA 1
   HETATM 64 HAQA . 1E8 A 701 ? 25.06494 10.04723 1.40233 1.000 21.53000 H ? B ? . HAQA 1
   HETATM 65 HARA . 1E8 A 701 ? 21.77717 9.01561 1.53734 1.000 21.89000 H ? B ? . HARA 1
   HETATM 66 HNAA . 1E8 A 701 ? 16.99677 5.14091 3.18518 1.000 12.25000 H ? B ? . HNAA 1
   HETATM 67 HNAB . 1E8 A 701 ? 16.27145 4.61777 2.02704 1.000 12.25000 H ? B ? . HNAB 1
"""

def tst_00():
  '''
    4zeb: ARA O2 - TT7 P4 phosphodiester (1.59 A, no LINK record) -> no HO2.
  '''
  compare_models(pdb_str  = ARA_A_601_4zeb_H,
                 sel_str  = "resname ARA and resseq 601 and name HO2",
                 optimize = True)

def tst_01():
  '''
    2ww3: Z5L S2 - MAN C1 thioglycoside (1.76 A) -> no HS2; anomeric H1 stays.
  '''
  compare_models(pdb_str  = Z5L_L_1_2ww3_H,
                 sel_str  = "resname Z5L and name HS2",
                 optimize = True)

def tst_02():
  '''
    8oji: YIO S1 - VPH C2' thioglycoside (1.83 A) -> no HS1; C2' keeps H5,
    H4 is the position S1 took (geostd ideal chirality).
  '''
  compare_models(pdb_str  = VPH_A_304_8oji_H,
                 sel_str  = "(resname YIO and name HS1) or (resname VPH and name H4)",
                 optimize = True)

def tst_03():
  '''
    2b5z: LYS 96 NZ alkylated by two BGS (both NZ-CS LINKs deposited) ->
    tertiary amine without HZ; each CS a CH2. Needs NZ's second link past
    linking_setup.maximum_per_atom_links (covalent range, free valence).
  '''
  compare_models(pdb_str  = BGS_A_961_2b5z_H,
                 sel_str  = "resname LYS and resseq 96 and name HZ*",
                 optimize = True)

def tst_04():
  '''
    4jxg: SER 64 OG - 1S6 C acyl-enzyme ester (1.38 A, LINK deposited). C is
    then a carbonyl (CA, =OXT, OG) -> no H1; SER OG -> no HG.
    Known issue, kept on purpose: H20 (ring N NAP) is on the face 1.92 A from C.
  '''
  compare_models(pdb_str  = SER_1S6_B_401_4jxg_H,
                 sel_str  = "(resname 1S6 and name H1) or (resname SER and name HG)",
                 optimize = True)

def tst_05():
  '''
    5p9j: CYS 481 SG - CAA ibrutinib (Michael adduct), labelled with the free
    drug code 1E8. Its dictionary still has CAA=CAD, which the reaction consumed
    (CAA-CAD 1.60 A): CAA keeps both H (HAA, HAAA); SG -> no HG.
    Limitation, kept on purpose: 1E8 gives CAD one trigonal H (HAD); the adduct
    has two. The PDB uses the reacted form 8E8.
  '''
  compare_models(pdb_str  = CYS_1E8_A_701_5p9j_H,
                 sel_str  = "resname CYS and name HG",
                 optimize = True)

def run():
  tst_00()
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  tst_05()

if __name__ == '__main__':
  run()
