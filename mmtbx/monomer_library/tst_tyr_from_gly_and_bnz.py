from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
from itertools import count
from six.moves import cStringIO as StringIO
import sys
from six.moves import zip

pdb_records_gly_bnz = """\
CRYST1   56.470   56.470   80.520  90.00  90.00  90.00 P 43 21 2     8
ATOM    141  N   GLY     1      19.851  38.691  38.859  1.00 14.08           N
ATOM    142  CA  GLY     1      19.095  38.743  40.099  1.00 15.88           C
ATOM    143  C   GLY     1      17.786  37.966  39.974  1.00 18.62           C
ATOM    144  O   GLY     1      17.440  37.187  40.865  1.00 18.84           O
ATOM    145  CB  GLY     1      18.798  40.190  40.494  1.00 15.26           C
ATOM    146  CG  BNZ     2      17.817  40.296  41.640  1.00 16.71           C
ATOM    147  CD1 BNZ     2      18.228  40.074  42.953  1.00 19.01           C
ATOM    148  CD2 BNZ     2      16.473  40.584  41.411  1.00 16.09           C
ATOM    149  CE1 BNZ     2      17.322  40.135  44.011  1.00 17.90           C
ATOM    150  CE2 BNZ     2      15.560  40.647  42.461  1.00 16.31           C
ATOM    151  CZ  BNZ     2      15.994  40.423  43.758  1.00 19.02           C
ATOM    152  OH  BNZ     2      15.111  40.507  44.807  1.00 21.58           O
END
""".splitlines()

pdb_records_tyr = [line
  .replace(" GLY ", " TYR ")
  .replace(" BNZ     2  ", " TYR     1 ") for line in pdb_records_gly_bnz]

cif_string = """\
#
# Copy of monomer library entries 2006_01_03, to ensure the tests below
# are independent of changes to the BNZ, GLY, TYR definitions.
# Changes with respect to the original monomer library definitions
# are marked with "# ADJUSTMENT".
#

data_comp_BNZ
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 BNZ           C6     C    CR16      0.000
 BNZ           H6     H    HCR6      0.000
 BNZ           C5     C    CR16      0.000
 BNZ           H5     H    HCR6      0.000
 BNZ           C4     C    CR16      0.000
 BNZ           H4     H    HCR6      0.000
 BNZ           C3     C    CR16      0.000
 BNZ           H3     H    HCR6      0.000
 BNZ           C2     C    CR16      0.000
 BNZ           H2     H    HCR6      0.000
 BNZ           C1     C    CR16      0.000
 BNZ           H1     H    HCR6      0.000
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 BNZ      C6     n/a    C5     START
 BNZ      H6     C6     .      .
 BNZ      C5     C6     C4     .
 BNZ      H5     C5     .      .
 BNZ      C4     C5     C3     .
 BNZ      H4     C4     .      .
 BNZ      C3     C4     C2     .
 BNZ      H3     C3     .      .
 BNZ      C2     C3     C1     .
 BNZ      H2     C2     .      .
 BNZ      C1     C2     H1     .
 BNZ      H1     C1     .      END
 BNZ      C6     C1     .    ADD
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 BNZ      H6     C6        coval       0.960    0.020
 BNZ      C5     C6        coval       1.390    0.020
 BNZ      H5     C5        coval       0.960    0.020
 BNZ      C4     C5        coval       1.390    0.020
 BNZ      H4     C4        coval       0.960    0.020
 BNZ      C3     C4        coval       1.390    0.020
 BNZ      H3     C3        coval       0.960    0.020
 BNZ      C2     C3        coval       1.390    0.020
 BNZ      H2     C2        coval       0.960    0.020
 BNZ      C1     C2        coval       1.390    0.020
 BNZ      C1     C6        coval       1.390    0.020
 BNZ      H1     C1        coval       0.960    0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 BNZ      H6     C6     C5      120.000    3.000
 BNZ      H6     C6     C1      120.000    3.000
 BNZ      C5     C6     C1      120.000    3.000
 BNZ      C6     C5     H5      120.000    3.000
 BNZ      C6     C5     C4      120.000    3.000
 BNZ      H5     C5     C4      120.000    3.000
 BNZ      C5     C4     H4      120.000    3.000
 BNZ      C5     C4     C3      120.000    3.000
 BNZ      H4     C4     C3      120.000    3.000
 BNZ      C4     C3     H3      120.000    3.000
 BNZ      C4     C3     C2      120.000    3.000
 BNZ      H3     C3     C2      120.000    3.000
 BNZ      C3     C2     H2      120.000    3.000
 BNZ      C3     C2     C1      120.000    3.000
 BNZ      H2     C2     C1      120.000    3.000
 BNZ      C2     C1     H1      120.000    3.000
 BNZ      C2     C1     C6      120.000    3.000
 BNZ      H1     C1     C6      120.000    3.000
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 BNZ      CONST_1  C2     C1     C6     C5         0.000    0.000   0
 BNZ      CONST_2  C1     C6     C5     C4         0.000    0.000   0
 BNZ      CONST_4  C6     C5     C4     C3         0.000    0.000   0
 BNZ      CONST_5  C5     C4     C3     C2         0.000    0.000   0
 BNZ      CONST_6  C4     C3     C2     C1         0.000    0.000   0
 BNZ      CONST_7  C3     C2     C1     C6         0.000    0.000   0
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 BNZ      plan    C1        0.020 # ADJUSTMENT: all "plan-1" renamed to "plan"
 BNZ      plan    C2        0.020
 BNZ      plan    C3        0.020
 BNZ      plan    C4        0.020
 BNZ      plan    C5        0.020
 BNZ      plan    C6        0.020
 BNZ      plan    H1        0.020
 BNZ      plan    H2        0.020
 BNZ      plan    H3        0.020
 BNZ      plan    H4        0.020
 BNZ      plan    H5        0.020
 BNZ      plan    H6        0.020

data_comp_GLY
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 GLY           N      N    NH1      -0.204
 GLY           H      H    HNH1      0.204
 GLY           CA     C    CH2       0.002
 GLY           HA1    H    HCH2      0.051
 GLY           HA2    H    HCH2      0.051
 GLY           C      C    C         0.318
 GLY           O      O    O        -0.422
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 GLY      N      n/a    CA     START
 GLY      H      N      .      .
 GLY      CA     N      C      .
 GLY      HA1    CA     .      .
 GLY      HA2    CA     .      .
 GLY      C      CA     .      END
 GLY      O      C      .      .
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 GLY      N      H         coval       0.860    0.020
 GLY      N      CA        coval       1.451    0.016
 GLY      CA     HA1       coval       0.970    0.020
 GLY      CA     HA2       coval       0.970    0.020
 GLY      CA     C         coval       1.516    0.018
 GLY      C      O         coval       1.231    0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 GLY      H      N      CA      114.000    3.000
 GLY      HA1    CA     HA2     109.000    3.000
 GLY      HA2    CA     C       109.000    3.000
 GLY      HA1    CA     C       109.000    3.000
 GLY      N      CA     HA1     110.000    3.000
 GLY      N      CA     HA2     110.000    3.000
 GLY      N      CA     C       112.500    2.900
 GLY      CA     C      O       120.800    2.100

data_comp_TYR
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 TYR           N      N    NH1      -0.204
 TYR           H      H    HNH1      0.204
 TYR           CA     C    CH1       0.058
 TYR           HA     H    HCH1      0.046
 TYR           CB     C    CH2      -0.054
 TYR           HB1    H    HCH2      0.049
 TYR           HB2    H    HCH2      0.049
 TYR           CG     C    CR6      -0.044
 TYR           CD1    C    CR16     -0.053
 TYR           HD1    H    HCR6      0.053
 TYR           CE1    C    CR16     -0.099
 TYR           HE1    H    HCR6      0.054
 TYR           CZ     C    CR6       0.176
 TYR           OH     O    OH1      -0.391
 TYR           HH     H    HOH1      0.305
 TYR           CE2    C    CR16     -0.099
 TYR           HE2    H    HCR6      0.054
 TYR           CD2    C    CR16     -0.053
 TYR           HD2    H    HCR6      0.053
 TYR           C      C    C         0.318
 TYR           O      O    O        -0.422
loop_
_chem_comp_tree.comp_id
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
 TYR      N      n/a    CA     START
 TYR      H      N      .      .
 TYR      CA     N      C      .
 TYR      HA     CA     .      .
 TYR      CB     CA     CG     .
 TYR      HB1    CB     .      .
 TYR      HB2    CB     .      .
 TYR      CG     CB     CD1    .
 TYR      CD1    CG     CE1    .
 TYR      HD1    CD1    .      .
 TYR      CE1    CD1    CZ     .
 TYR      HE1    CE1    .      .
 TYR      CZ     CE1    CE2    .
 TYR      OH     CZ     HH     .
 TYR      HH     OH     .      .
 TYR      CE2    CZ     CD2    .
 TYR      HE2    CE2    .      .
 TYR      CD2    CE2    HD2    .
 TYR      HD2    CD2    .      .
 TYR      C      CA     .      END
 TYR      O      C      .      .
 TYR      CD2    CG     .    ADD
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 TYR      N      H         coval       0.860    0.020
 TYR      N      CA        coval       1.458    0.019
 TYR      CA     HA        coval       0.980    0.020
 TYR      CA     CB        coval       1.530    0.020
 TYR      CB     HB1       coval       0.970    0.020
 TYR      CB     HB2       coval       0.970    0.020
 TYR      CB     CG        coval       1.512    0.022
 TYR      CG     CD2       coval       1.389    0.021
 TYR      CG     CD1       coval       1.389    0.021
 TYR      CD1    HD1       coval       0.930    0.020
 TYR      CD1    CE1       coval       1.382    0.030
 TYR      CE1    HE1       coval       0.930    0.020
 TYR      CE1    CZ        coval       1.378    0.024
 TYR      CZ     OH        coval       1.376    0.021
 TYR      OH     HH        coval       0.820    0.020
 TYR      CZ     CE2       coval       1.378    0.024
 TYR      CE2    HE2       coval       0.930    0.020
 TYR      CE2    CD2       coval       1.382    0.030
 TYR      CD2    HD2       coval       0.930    0.020
 TYR      CA     C         coval       1.525    0.021
 TYR      C      O         coval       1.231    0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 TYR      H      N      CA      114.000    3.000
 TYR      HA     CA     CB      109.000    3.000
 TYR      CB     CA     C       110.100    1.900
 TYR      HA     CA     C       109.000    3.000
 TYR      N      CA     HA      110.000    3.000
 TYR      N      CA     CB      110.500    1.700
 TYR      HB1    CB     HB2     110.000    3.000
 TYR      HB2    CB     CG      108.000    3.000
 TYR      HB1    CB     CG      108.000    3.000
 TYR      CA     CB     HB1     109.000    3.000
 TYR      CA     CB     HB2     109.000    3.000
 TYR      CA     CB     CG      113.900    1.800
 TYR      CB     CG     CD2     120.800    1.500
 TYR      CD1    CG     CD2     118.100    1.500
 TYR      CB     CG     CD1     120.800    1.500
 TYR      HD1    CD1    CE1     119.400    3.000
 TYR      CG     CD1    HD1     119.400    3.000
 TYR      CG     CD1    CE1     121.200    1.500
 TYR      HD2    CD2    CE2     119.400    3.000
 TYR      CG     CD2    HD2     119.400    3.000
 TYR      CG     CD2    CE2     121.200    1.500
 TYR      HE1    CE1    CZ      120.200    3.000
 TYR      CD1    CE1    HE1     120.200    3.000
 TYR      CD1    CE1    CZ      119.600    1.800
 TYR      OH     CZ     CE2     119.900    3.000
 TYR      CE1    CZ     OH      119.900    3.000
 TYR      CZ     OH     HH      110.000    3.000
 TYR      CE1    CZ     CE2     120.300    2.000
 TYR      HE2    CE2    CD2     120.200    3.000
 TYR      CZ     CE2    HE2     120.200    3.000
 TYR      CZ     CE2    CD2     119.600    1.800
 TYR      N      CA     C       111.200    2.800
 TYR      CA     C      O       120.800    1.700
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 TYR      chi1     N      CA     CB     CG       180.000   15.000   3
 TYR      chi2     CA     CB     CG     CD1       90.000   20.000   2
 TYR      CONST_01 CB     CG     CD1    CE1      180.000    0.000   0
 TYR      CONST_02 CG     CD1    CE1    CZ         0.000    0.000   0
 TYR      CONST_03 CD1    CE1    CZ     CE2        0.000    0.000   0
 TYR      hh1      CE1    CZ     OH     HH        60.000   30.000   2
 TYR      CONST_04 CE1    CZ     CE2    CD2        0.000    0.000   0
 TYR      CONST_05 CZ     CE2    CD2    CG         0.000    0.000   0
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 TYR      chir_01  CA     N      CB     C         negativ
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 TYR      plan      CB        0.020
 TYR      plan      CG        0.020
 TYR      plan      CD1       0.020
 TYR      plan      CE1       0.020
 TYR      plan      CZ        0.020
 TYR      plan      CE2       0.020
 TYR      plan      CD2       0.020
 TYR      plan      OH        0.020
 TYR      plan      HD1       0.020
 TYR      plan      HE1       0.020
 TYR      plan      HE2       0.020
 TYR      plan      HD2       0.020

#
# end copy of monomer library entries
#

#
# Modification of GLY, to add a CB
#
data_mod_gly_plus_c_beta
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_partial_charge
 gly_plus_c_beta   change   CA     .      .    CH1     0.058
 gly_plus_c_beta   change   HA1    HA     .    HCH1    0.046
 gly_plus_c_beta   delete   HA2    .      .    .        .
 gly_plus_c_beta   add      .      CB     C    CH2    -0.054
 gly_plus_c_beta   add      .      HB1    H    HCH2    0.049
 gly_plus_c_beta   add      .      HB2    H    HCH2    0.049
loop_
_chem_mod_bond.mod_id
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist
_chem_mod_bond.new_value_dist_esd
 gly_plus_c_beta   change   CA       C         .            1.525    0.021
 gly_plus_c_beta   change   CA       N         .            1.458    0.019
 gly_plus_c_beta   change   CA       HA        .            0.980     .
 gly_plus_c_beta   add      CA       CB        coval        1.530    0.020
 gly_plus_c_beta   add      CB       HB1       coval        0.970    0.020
 gly_plus_c_beta   add      CB       HB2       coval        0.970    0.020
loop_
_chem_mod_angle.mod_id
_chem_mod_angle.function
_chem_mod_angle.atom_id_1
_chem_mod_angle.atom_id_2
_chem_mod_angle.atom_id_3
_chem_mod_angle.new_value_angle
_chem_mod_angle.new_value_angle_esd
 gly_plus_c_beta   add      C       CA      CB      110.1      1.900
 gly_plus_c_beta   change   C       CA      N       111.2      2.800
 gly_plus_c_beta   change   CA      C       O          .       1.700
 gly_plus_c_beta   add      CA      CB      HB1     109.0      3.000
 gly_plus_c_beta   add      CA      CB      HB2     109.0      3.000
 gly_plus_c_beta   add      CB      CA      HA      109.0      3.000
 gly_plus_c_beta   add      CB      CA      N       110.5      1.700
 gly_plus_c_beta   add      HB1     CB      HB2     110.0      3.000
loop_
_chem_mod_chir.mod_id
_chem_mod_chir.function
_chem_mod_chir.atom_id_centre
_chem_mod_chir.atom_id_1
_chem_mod_chir.atom_id_2
_chem_mod_chir.atom_id_3
_chem_mod_chir.new_volume_sign
 gly_plus_c_beta  add  CA     N      CB     C         negativ

#
# Modification of BNZ (benzene), to make it into a TYR sidechain,
# but without the CB
#
data_mod_bnz_to_tyr_sidechain
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_partial_charge
 bnz_to_tyr_sidechain   change   C1     CZ     .    CR6     0.176
 bnz_to_tyr_sidechain   delete   H1     .      .    .        .
 bnz_to_tyr_sidechain   add      .      OH     O    OH1    -0.391
 bnz_to_tyr_sidechain   add      .      HH     H    HOH1    0.305
 bnz_to_tyr_sidechain   change   C2     CE1    .    CR16   -0.099
 bnz_to_tyr_sidechain   change   H2     HE1    .    HCR6    0.054
 bnz_to_tyr_sidechain   change   C3     CD1    .    CR16   -0.053
 bnz_to_tyr_sidechain   change   H3     HD1    .    HCR6    0.053
 bnz_to_tyr_sidechain   change   C4     CG     .    CR6    -0.044
 bnz_to_tyr_sidechain   delete   H4     .      .    .        .
 bnz_to_tyr_sidechain   change   C5     CD2    .    CR16   -0.053
 bnz_to_tyr_sidechain   change   H5     HD2    .    HCR6    0.053
 bnz_to_tyr_sidechain   change   C6     CE2    .    CR16   -0.099
 bnz_to_tyr_sidechain   change   H6     HE2    .    HCR6    0.054
loop_
_chem_mod_bond.mod_id
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist
_chem_mod_bond.new_value_dist_esd
 bnz_to_tyr_sidechain   change   CD1      CE1       .            1.382    0.030
 bnz_to_tyr_sidechain   change   CD1      CG        .            1.389    0.021
 bnz_to_tyr_sidechain   change   CD1      HD1       .            0.930    0.020
 bnz_to_tyr_sidechain   change   CD2      CE2       .            1.382    0.030
 bnz_to_tyr_sidechain   change   CD2      CG        .            1.389    0.021
 bnz_to_tyr_sidechain   change   CD2      HD2       .            0.930    0.020
 bnz_to_tyr_sidechain   change   CE1      CZ        .            1.378    0.024
 bnz_to_tyr_sidechain   change   CE1      HE1       .            0.930    0.020
 bnz_to_tyr_sidechain   change   CE2      CZ        .            1.378    0.024
 bnz_to_tyr_sidechain   change   CE2      HE2       .            0.930    0.020
 bnz_to_tyr_sidechain   add      CZ       OH        coval        1.376    0.021
 bnz_to_tyr_sidechain   add      HH       OH        coval        0.820    0.020
loop_
_chem_mod_angle.mod_id
_chem_mod_angle.function
_chem_mod_angle.atom_id_1
_chem_mod_angle.atom_id_2
_chem_mod_angle.atom_id_3
_chem_mod_angle.new_value_angle
_chem_mod_angle.new_value_angle_esd
 bnz_to_tyr_sidechain   change   CD1     CE1     CZ      119.6      1.800
 bnz_to_tyr_sidechain   change   CD1     CE1     HE1     120.2      3.000
 bnz_to_tyr_sidechain   change   CD1     CG      CD2     118.1      1.500
 bnz_to_tyr_sidechain   change   CD2     CE2     CZ      119.6      1.800
 bnz_to_tyr_sidechain   change   CD2     CE2     HE2     120.2       .
 bnz_to_tyr_sidechain   change   CE1     CD1     CG      121.2      1.500
 bnz_to_tyr_sidechain   change   CE1     CD1     CG      121.2      1.500
 bnz_to_tyr_sidechain   change   CE1     CD1     HD1     119.4       .
 bnz_to_tyr_sidechain   change   CE1     CZ      CE2     120.3      2.000
 bnz_to_tyr_sidechain   add      CE1     CZ      OH      119.9      3.000
 bnz_to_tyr_sidechain   change   CE2     CD2     CG      121.2      1.500
 bnz_to_tyr_sidechain   change   CE2     CD2     HD2     119.4       .
 bnz_to_tyr_sidechain   add      CE2     CZ      OH      119.9      3.000
 bnz_to_tyr_sidechain   change   CG      CD1     HD1     119.4       .
 bnz_to_tyr_sidechain   change   CG      CD2     HD2     119.4       .
 bnz_to_tyr_sidechain   change   CZ      CE1     HE1     120.2       .
 bnz_to_tyr_sidechain   change   CZ      CE2     HE2     120.2       .
 bnz_to_tyr_sidechain   add      CZ      OH      HH      110.0      3.000
loop_
_chem_mod_tor.mod_id
_chem_mod_tor.function
_chem_mod_tor.id
_chem_mod_tor.atom_id_1
_chem_mod_tor.atom_id_2
_chem_mod_tor.atom_id_3
_chem_mod_tor.atom_id_4
_chem_mod_tor.new_value_angle
_chem_mod_tor.new_value_angle_esd
_chem_mod_tor.new_period
 bnz_to_tyr_sidechain  delete  CONST_1   CE1   CZ    CE2   CD2      .     .   .
 bnz_to_tyr_sidechain  add     CONST_04  CE1   CZ    CE2   CD2     0.0   0.0  0
 bnz_to_tyr_sidechain  delete  CONST_2   CZ    CE2   CD2   CG       .     .   .
 bnz_to_tyr_sidechain  add     CONST_05  CZ    CE2   CD2   CG      0.0   0.0  0
 bnz_to_tyr_sidechain  delete  CONST_4   CE2   CD2   CG    CD1      .     .   .
 bnz_to_tyr_sidechain  delete  CONST_5   CD2   CG    CD1   CE1      .     .   .
 bnz_to_tyr_sidechain  add     hh1       CE1   CZ    OH    HH     60.0  30.0  2
 bnz_to_tyr_sidechain  delete  CONST_6   CG    CD1   CE1   CZ       .     .   .
 bnz_to_tyr_sidechain  add     CONST_02  CG    CD1   CE1   CZ      0.0   0.0  0
 bnz_to_tyr_sidechain  delete  CONST_7   CD1   CE1   CZ    CE2      .     .   .
 bnz_to_tyr_sidechain  add     CONST_03  CD1   CE1   CZ    CE2     0.0   0.0  0
loop_
_chem_mod_plane_atom.mod_id
_chem_mod_plane_atom.function
_chem_mod_plane_atom.plane_id
_chem_mod_plane_atom.atom_id
_chem_mod_plane_atom.new_dist_esd
 bnz_to_tyr_sidechain  add       plan   OH        0.020

data_link_list
loop_
_chem_link.id
_chem_link.comp_id_1
_chem_link.mod_id_1
_chem_link.group_comp_1
_chem_link.comp_id_2
_chem_link.mod_id_2
_chem_link.group_comp_2
_chem_link.name
gly_bnz_to_tyr  GLY  gly_plus_c_beta  .  BNZ  bnz_to_tyr_sidechain  . .

data_link_gly_bnz_to_tyr
loop_
_chem_link_bond.link_id
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
 gly_bnz_to_tyr  1 CB      2 CG        .           1.512    0.022
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
 gly_bnz_to_tyr  1 CA      1 CB      2 CG      113.900    1.800
 gly_bnz_to_tyr  1 CB      2 CG      2 CD1     120.800    1.500
 gly_bnz_to_tyr  1 CB      2 CG      2 CD2     120.800    1.500
 gly_bnz_to_tyr  2 CG      1 CB      1 HB1     108.000    3.000
 gly_bnz_to_tyr  2 CG      1 CB      1 HB2     108.000    3.000
loop_
_chem_link_tor.link_id
_chem_link_tor.id
_chem_link_tor.atom_1_comp_id
_chem_link_tor.atom_id_1
_chem_link_tor.atom_2_comp_id
_chem_link_tor.atom_id_2
_chem_link_tor.atom_3_comp_id
_chem_link_tor.atom_id_3
_chem_link_tor.atom_4_comp_id
_chem_link_tor.atom_id_4
_chem_link_tor.value_angle
_chem_link_tor.value_angle_esd
_chem_link_tor.period
 gly_bnz_to_tyr  CONST_01  1 CB     2 CG     2 CD1    2 CE1     180.00   0.0 0
 gly_bnz_to_tyr  chi1      1 N      1 CA     1 CB     2 CG      180.00  15.0 3
 gly_bnz_to_tyr  chi2      1 CA     1 CB     2 CG     2 CD1      90.00  20.0 2
loop_
_chem_link_plane.link_id
_chem_link_plane.plane_id
_chem_link_plane.atom_comp_id
_chem_link_plane.atom_id
_chem_link_plane.dist_esd
 gly_bnz_to_tyr  plan   1 CB    0.020
"""

def exercise(args):
  if ("--verbose" in args):
    out = sys.stdout
  else:
    out = StringIO()
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  with open("tmp.cif", "w") as f:
    f.write(cif_string)
  mon_lib_srv.process_cif(file_name="tmp.cif")
  mod_gly_plus_c_beta = mon_lib_srv.mod_mod_id_dict[
    "gly_plus_c_beta"]
  mod_bnz_to_tyr_sidechain = mon_lib_srv.mod_mod_id_dict[
    "bnz_to_tyr_sidechain"]
  #
  pdb_inp = pdb.input(
    source_info=None, lines=flex.std_string(pdb_records_gly_bnz))
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  residues = pdb_hierarchy.only_conformer().residues()
  assert len(residues) == 2
  monomer_definitions = []
  for residue in residues:
    md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
      residue_name=residue.resname,
      atom_names=residue.atoms().extract_name())
    if (ani is not None):
      assert residue.resname == "GLY"
      assert ani.mon_lib_names() == ["N", "CA", "C", "O", None]
    monomer_definitions.append(md)
  gly_plus_c_beta = monomer_definitions[0].apply_mod(mod=mod_gly_plus_c_beta)
  gly_plus_c_beta.show(f=out)
  sidechain = monomer_definitions[1].apply_mod(mod=mod_bnz_to_tyr_sidechain)
  sidechain.show(f=out)
  #
  tyr_definition = mon_lib_srv.get_comp_comp_id_direct("TYR")
  tyr_atoms = []
  for atom in tyr_definition.atom_list:
    s = StringIO()
    atom.show(f=s)
    tyr_atoms.append(s.getvalue())
  tyr_atoms.sort()
  #
  mod_atoms = []
  for atom in gly_plus_c_beta.atom_list + sidechain.atom_list:
    s = StringIO()
    atom.show(f=s)
    mod_atoms.append(s.getvalue())
  mod_atoms.sort()
  assert not show_diff("".join(mod_atoms), "".join(tyr_atoms))
  #
  tyr_bonds = []
  for bond in tyr_definition.normalized_bond_list():
    s = StringIO()
    bond.show(f=s)
    tyr_bonds.append(s.getvalue())
  tyr_bonds.sort()
  #
  mod_bonds = []
  for bond in   gly_plus_c_beta.normalized_bond_list() \
              + sidechain.normalized_bond_list():
    s = StringIO()
    bond.show(f=s)
    mod_bonds.append(s.getvalue())
  mod_bonds.sort()
  mod_bonds.insert(5, None)
  mod_bonds += [None] * (len(tyr_bonds) - len(mod_bonds))
  for i,t,m in zip(count(), tyr_bonds, mod_bonds):
    print("bond index:", i, file=out)
    if (m is None):
      print(t, file=out)
    else:
      t = t.splitlines()
      m = m.splitlines()
      assert len(m) == len(t)
      for ti,mi in zip(t, m):
        print("%-39s %-39s" % (ti,mi), file=out)
    print(file=out)
  assert not show_diff("".join(mod_bonds[:5]), "".join(tyr_bonds[:5]))
  assert not show_diff("".join(mod_bonds[6]), "".join(tyr_bonds[6]))
  assert not show_diff("".join(mod_bonds[8:]), "".join(tyr_bonds[8:]))
  #
  tyr_angles = []
  for angle in tyr_definition.normalized_angle_list():
    s = StringIO()
    angle.show(f=s)
    tyr_angles.append(s.getvalue())
  tyr_angles.sort()
  #
  mod_angles = []
  for angle in  gly_plus_c_beta.normalized_angle_list() \
              + sidechain.normalized_angle_list():
    s = StringIO()
    angle.show(f=s)
    mod_angles.append(s.getvalue())
  mod_angles.sort()
  mod_angles.insert(4, None)
  mod_angles.insert(10, None)
  mod_angles.insert(11, None)
  mod_angles.insert(24, None)
  mod_angles.insert(25, None)
  mod_angles += [None] * (len(tyr_angles) - len(mod_angles))
  for i,t,m in zip(count(), tyr_angles, mod_angles):
    print("angle index:", i, file=out)
    if (m is None):
      print(t, file=out)
    else:
      t = t.splitlines()
      m = m.splitlines()
      assert len(m) == len(t)
      for ti,mi in zip(t, m):
        print("%-39s %-39s" % (ti,mi), file=out)
      print(file=out)
  assert not show_diff("".join(mod_angles[:4]), "".join(tyr_angles[:4]))
  assert not show_diff("".join(mod_angles[5:10]), "".join(tyr_angles[5:10]))
  assert not show_diff("".join(mod_angles[12:24]), "".join(tyr_angles[12:24]))
  assert not show_diff("".join(mod_angles[26:]), "".join(tyr_angles[26:]))
  #
  tyr_tors = []
  for tor in tyr_definition.tor_list:
    s = StringIO()
    tor.show(f=s)
    tyr_tors.append(s.getvalue())
  tyr_tors.sort()
  #
  mod_tors = []
  for tor in  gly_plus_c_beta.tor_list \
              + sidechain.tor_list:
    s = StringIO()
    tor.show(f=s)
    mod_tors.append(s.getvalue())
  mod_tors.sort()
  mod_tors.insert(0, None)
  mod_tors.insert(5, None)
  mod_tors.insert(6, None)
  for i,t,m in zip(count(), tyr_tors, mod_tors):
    print("tor index:", i, file=out)
    if (m is None):
      print(t, file=out)
    else:
      t = t.splitlines()
      m = m.splitlines()
      assert len(m) == len(t)
      for ti,mi in zip(t, m):
        print("%-39s %-39s" % (ti,mi), file=out)
      print(file=out)
  assert not show_diff("".join(mod_tors[1:5]), "".join(tyr_tors[1:5]))
  assert not show_diff("".join(mod_tors[7]), "".join(tyr_tors[7]))
  #
  assert len(tyr_definition.chir_list) == 1
  assert len(gly_plus_c_beta.chir_list) == 1
  assert len(sidechain.chir_list) == 0
  gly_plus_c_beta.chir_list[0].id = tyr_definition.chir_list[0].id
  s = StringIO()
  tyr_definition.chir_list[0].show(f=s)
  tyr_chir = s.getvalue()
  s = StringIO()
  gly_plus_c_beta.chir_list[0].show(f=s)
  mod_chir = s.getvalue()
  assert not show_diff(mod_chir, tyr_chir)
  #
  tyr_plane_atoms = []
  for plane_atom in tyr_definition.plane_atom_list:
    s = StringIO()
    plane_atom.show(f=s)
    tyr_plane_atoms.append(s.getvalue())
  tyr_plane_atoms.sort()
  #
  mod_plane_atoms = []
  for plane_atom in  gly_plus_c_beta.plane_atom_list \
                   + sidechain.plane_atom_list:
    s = StringIO()
    plane_atom.show(f=s)
    mod_plane_atoms.append(s.getvalue())
  mod_plane_atoms.sort()
  mod_plane_atoms.insert(0, None)
  mod_plane_atoms += [None] * (len(tyr_plane_atoms) - len(mod_plane_atoms))
  for i,t,m in zip(count(), tyr_plane_atoms, mod_plane_atoms):
    print("plane_atom index:", i, file=out)
    if (m is None):
      print(t, file=out)
    else:
      t = t.splitlines()
      m = m.splitlines()
      assert len(m) == len(t)
      for ti,mi in zip(t, m):
        print("%-39s %-39s" % (ti,mi), file=out)
      print(file=out)
  assert not show_diff(
    "".join(mod_plane_atoms[1:]), "".join(tyr_plane_atoms[1:]))
  #
  processed_tyr = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=None,
    raw_records=pdb_records_tyr,
    log=sys.stdout)
  print()
  processed_gly_bnz = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=None,
    raw_records=pdb_records_gly_bnz,
    log=sys.stdout)
  print()
  print("geo_tyr")
  geo_tyr = processed_tyr.geometry_restraints_manager()
  print()
  print("geo_gly_bnz")
  geo_gly_bnz = processed_gly_bnz.geometry_restraints_manager()
  print()
  #
  print("TODO: compare geometry restraints")
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
