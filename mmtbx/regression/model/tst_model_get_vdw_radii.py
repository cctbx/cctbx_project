from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.utils
import mmtbx.model
from libtbx.utils import null_out

pdb = """
CRYST1   31.000   18.000   43.000  90.00  90.00  90.00 P 2 21 21
HETATM    1  C24 8ZB C   1     -11.000  13.211  10.894  1.00 30.00           C
HETATM    2  C28 8ZB C   1     -11.000  11.722  10.295  1.00 30.00           C
HETATM    3  C31 8ZB C   1     -10.000  11.074  10.706  1.00 30.00           C
HETATM    4  C36 8ZB C   1     -10.000  19.626  10.001  1.00 30.00           C
HETATM    5  C38 8ZB C   1      -8.000  18.259  10.474  1.00 30.00           C
HETATM    6  C39 8ZB C   1      -8.000  18.457  10.949  1.00 30.00           C
HETATM    7  C40 8ZB C   1      -9.000  19.327  10.279  1.00 30.00           C
HETATM    8  C49 8ZB C   1      -6.000  19.546  10.455  1.00 30.00           C
HETATM    9  C50 8ZB C   1      -6.000  18.069  10.155  1.00 30.00           C
HETATM   10  C51 8ZB C   1      -7.000  17.489  10.934  1.00 30.00           C
HETATM   11  C53 8ZB C   1      -6.000  15.539  10.000  1.00 30.00           C
HETATM   12  C55 8ZB C   1      -5.000  17.242  10.321  1.00 30.00           C
HETATM   13  N47 8ZB C   1      -7.000  18.437  10.617  1.00 30.00           N
HETATM   14  N48 8ZB C   1      -5.000  19.513  10.551  1.00 30.00           N
HETATM   15  N52 8ZB C   1      -7.000  16.248  10.879  1.00 30.00           N
HETATM   16  N54 8ZB C   1      -5.000  15.993  10.237  1.00 30.00           N
HETATM   17  N57 8ZB C   1      -4.000  17.770  10.517  1.00 30.00           N
HETATM   18  O37 8ZB C   1      -9.000  19.010  10.951  1.00 30.00           O
HETATM   19  O41 8ZB C   1      -7.000  18.355  10.894  1.00 30.00           O
HETATM   20  O42 8ZB C   1     -10.000  18.384  10.923  1.00 30.00           O
TER
END
"""

cif = """
data_comp_list
#
#
#
#
#
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
8ZB 8ZB
'_2__R_,3__R_,4__S_,5__R__-2-_6-aminopurin-9-yl_-5-propyl-oxolane-3,4-diol'
non-polymer 37 20 .
data_comp_8ZB
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
8ZB C24 C CH3  -0.065 5.193  -2.668 -0.514
8ZB C28 C CH2  -0.054 4.269  -1.547 -0.941
8ZB C31 C CH2  -0.024 3.677  -0.800 0.239
8ZB C36 C CH1  0.090  2.849  0.415  -0.131
8ZB C38 C CH1  0.167  0.586  0.822  -0.595
8ZB C39 C CH1  0.128  1.070  1.901  0.383
8ZB C40 C CH1  0.111  2.276  1.217  1.029
8ZB C49 C CR16 0.100  -0.400 -1.049 0.802
8ZB C50 C CR6  0.147  -2.470 -0.881 0.363
8ZB C51 C CR6  0.167  -1.837 0.107  -0.374
8ZB C53 C CR16 0.120  -3.743 0.809  -1.235
8ZB C55 C CR6  0.147  -3.871 -0.970 0.239
8ZB N47 N NR6  -0.285 -0.499 -0.007 -0.082
8ZB N48 N NRD6 -0.231 -1.549 -1.607 1.103
8ZB N52 N NRD6 -0.216 -2.419 0.988  -1.196
8ZB N54 N NRD6 -0.219 -4.486 -0.092 -0.584
8ZB N57 N NH2  -0.341 -4.618 -1.871 0.885
8ZB O37 O O2   -0.348 1.695  -0.016 -0.883
8ZB O41 O OH1  -0.385 1.438  3.076  -0.331
8ZB O42 O OH1  -0.387 3.232  2.149  1.524
8ZB H2  H HCH3 0.023  6.163  -2.277 -0.210
8ZB H1  H HCH3 0.023  4.799  -3.242 0.324
8ZB H29 H HCH3 0.023  5.383  -3.375 -1.320
8ZB H28 H HCH2 0.026  4.777  -0.858 -1.614
8ZB H33 H HCH2 0.026  3.458  -1.960 -1.539
8ZB H35 H HCH2 0.029  3.085  -1.495 0.832
8ZB H31 H HCH2 0.029  4.465  -0.458 0.907
8ZB H61 H HCH1 0.062  3.420  1.059  -0.799
8ZB H62 H HCH1 0.087  0.266  1.231  -1.553
8ZB H44 H HCH1 0.067  0.317  2.156  1.128
8ZB H43 H HCH1 0.065  1.967  0.592  1.865
8ZB H58 H HCR  0.103  0.549  -1.387 1.216
8ZB H56 H HCR  0.105  -4.305 1.478  -1.885
8ZB H60 H HNH2 0.144  -4.185 -2.546 1.507
8ZB H59 H HNH2 0.144  -5.627 -1.913 0.779
8ZB H46 H HOH1 0.210  0.856  3.798  0.026
8ZB H45 H HOH1 0.210  2.872  2.596  2.335
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
8ZB C24 C28 single   1.513 0.038
8ZB C24 H1  single   1.089 0.010
8ZB C24 H2  single   1.089 0.010
8ZB C24 H29 single   1.089 0.010
8ZB C28 C31 single   1.516 0.024
8ZB C28 H28 single   1.089 0.010
8ZB C28 H33 single   1.089 0.010
8ZB C31 C36 single   1.516 0.013
8ZB C31 H31 single   1.089 0.010
8ZB C31 H35 single   1.089 0.010
8ZB C36 C40 single   1.526 0.013
8ZB C36 H61 single   1.089 0.010
8ZB C36 O37 single   1.444 0.011
8ZB C38 C39 single   1.531 0.014
8ZB C38 H62 single   1.089 0.010
8ZB C38 N47 single   1.459 0.012
8ZB C38 O37 single   1.416 0.011
8ZB C39 C40 single   1.531 0.013
8ZB C39 H44 single   1.089 0.010
8ZB C39 O41 single   1.423 0.012
8ZB C40 H43 single   1.089 0.010
8ZB C40 O42 single   1.423 0.012
8ZB C49 H58 single   1.089 0.010
8ZB C49 N47 aromatic 1.370 0.008
8ZB C49 N48 aromatic 1.311 0.008
8ZB C50 C51 aromatic 1.388 0.011
8ZB C50 C55 aromatic 1.409 0.010
8ZB C50 N48 aromatic 1.387 0.007
8ZB C51 N47 aromatic 1.374 0.008
8ZB C51 N52 aromatic 1.338 0.012
8ZB C53 H56 single   1.089 0.010
8ZB C53 N52 aromatic 1.334 0.011
8ZB C53 N54 aromatic 1.334 0.011
8ZB C55 N54 aromatic 1.350 0.010
8ZB C55 N57 single   1.337 0.014
8ZB N57 H59 single   1.015 0.010
8ZB N57 H60 single   1.015 0.010
8ZB O41 H46 single   0.993 0.010
8ZB O42 H45 single   0.993 0.010
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
8ZB C28 C24 H1  112.7 3.0
8ZB C28 C24 H2  110.9 3.0
8ZB H1  C24 H2  107.3 3.0
8ZB C28 C24 H29 112.1 3.0
8ZB H1  C24 H29 106.8 3.0
8ZB H2  C24 H29 106.6 3.0
8ZB C24 C28 C31 112.7 1.9
8ZB C24 C28 H28 110.8 3.0
8ZB C31 C28 H28 110.5 3.0
8ZB C24 C28 H33 109.3 3.0
8ZB C31 C28 H33 109.0 3.0
8ZB H28 C28 H33 104.4 3.0
8ZB C28 C31 C36 114.6 1.8
8ZB C28 C31 H31 110.5 3.0
8ZB C36 C31 H31 107.1 3.0
8ZB C28 C31 H35 108.9 3.0
8ZB C36 C31 H35 110.4 3.0
8ZB H31 C31 H35 105.2 3.0
8ZB C31 C36 C40 116.2 1.2
8ZB C31 C36 H61 109.6 3.0
8ZB C40 C36 H61 110.7 3.0
8ZB C31 C36 O37 108.8 1.1
8ZB C40 C36 O37 105.2 1.3
8ZB H61 C36 O37 106.3 3.0
8ZB C39 C38 H62 113.0 3.0
8ZB C39 C38 N47 114.2 1.5
8ZB H62 C38 N47 107.6 3.0
8ZB C39 C38 O37 106.5 1.2
8ZB H62 C38 O37 105.9 3.0
8ZB N47 C38 O37 108.3 1.3
8ZB C38 C39 C40 101.5 1.2
8ZB C38 C39 H44 112.9 3.0
8ZB C40 C39 H44 111.5 3.0
8ZB C38 C39 O41 110.5 2.8
8ZB C40 C39 O41 111.9 2.7
8ZB H44 C39 O41 109.4 3.0
8ZB C36 C40 C39 102.5 1.0
8ZB C36 C40 H43 113.0 3.0
8ZB C39 C40 H43 111.0 3.0
8ZB C36 C40 O42 111.0 2.5
8ZB C39 C40 O42 111.9 2.7
8ZB H43 C40 O42 107.5 3.0
8ZB H58 C49 N47 122.9 3.0
8ZB H58 C49 N48 123.1 3.0
8ZB N47 C49 N48 114.1 0.7
8ZB C51 C50 C55 117.0 0.6
8ZB C51 C50 N48 110.7 0.5
8ZB C55 C50 N48 132.1 1.2
8ZB C50 C51 N47 105.8 0.5
8ZB C50 C51 N52 126.8 0.7
8ZB N47 C51 N52 127.0 1.1
8ZB H56 C53 N52 116.3 3.0
8ZB H56 C53 N54 114.9 3.0
8ZB N52 C53 N54 128.7 0.9
8ZB C50 C55 N54 117.6 1.0
8ZB C50 C55 N57 123.5 1.1
8ZB N54 C55 N57 118.2 1.2
8ZB C38 N47 C49 126.8 1.9
8ZB C38 N47 C51 126.8 1.8
8ZB C49 N47 C51 105.7 0.5
8ZB C49 N48 C50 103.7 0.5
8ZB C51 N52 C53 111.9 2.2
8ZB C53 N54 C55 118.6 1.0
8ZB C55 N57 H59 122.2 3.0
8ZB C55 N57 H60 120.4 3.0
8ZB H59 N57 H60 117.3 3.0
8ZB C36 O37 C38 109.4 1.6
8ZB C39 O41 H46 105.4 3.0
8ZB C40 O42 H45 109.5 3.0
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
8ZB chir_04 C36 O37 C40 C31 negative
8ZB chir_05 C38 O37 N47 C39 negative
8ZB chir_06 C39 O41 C40 C38 positive
8ZB chir_07 C40 O42 C39 C36 negative
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
8ZB csd-C50  C50 0.020
8ZB csd-C50  C51 0.020
8ZB csd-C50  C55 0.020
8ZB csd-C50  N48 0.020
8ZB csd-C51  C51 0.020
8ZB csd-C51  C50 0.020
8ZB csd-C51  N47 0.020
8ZB csd-C51  N52 0.020
8ZB csd-C55  C55 0.020
8ZB csd-C55  C50 0.020
8ZB csd-C55  N54 0.020
8ZB csd-C55  N57 0.020
8ZB csd-N47  N47 0.020
8ZB csd-N47  C38 0.020
8ZB csd-N47  C49 0.020
8ZB csd-N47  C51 0.020
8ZB qm-C49   C49 0.020
8ZB qm-C49   N47 0.020
8ZB qm-C49   N48 0.020
8ZB qm-C49   H58 0.020
8ZB qm-C53   C53 0.020
8ZB qm-C53   N52 0.020
8ZB qm-C53   N54 0.020
8ZB qm-C53   H56 0.020
8ZB qm-N57   N57 0.020
8ZB qm-N57   C55 0.020
8ZB qm-N57   H59 0.020
8ZB qm-N57   H60 0.020
8ZB qmf-01   C50 0.020
8ZB qmf-01   C55 0.020
8ZB qmf-01   N57 0.020
8ZB qmf-01   H59 0.020
8ZB ring5A-1 C49 0.020
8ZB ring5A-1 N47 0.020
8ZB ring5A-1 C51 0.020
8ZB ring5A-1 C50 0.020
8ZB ring5A-2 N47 0.020
8ZB ring5A-2 C51 0.020
8ZB ring5A-2 C50 0.020
8ZB ring5A-2 N48 0.020
8ZB ring5A-3 C51 0.020
8ZB ring5A-3 C50 0.020
8ZB ring5A-3 N48 0.020
8ZB ring5A-3 C49 0.020
8ZB ring5A-4 C50 0.020
8ZB ring5A-4 N48 0.020
8ZB ring5A-4 C49 0.020
8ZB ring5A-4 N47 0.020
8ZB ring5A-5 N48 0.020
8ZB ring5A-5 C49 0.020
8ZB ring5A-5 N47 0.020
8ZB ring5A-5 C51 0.020
8ZB ring6A-1 C50 0.020
8ZB ring6A-1 C51 0.020
8ZB ring6A-1 N52 0.020
8ZB ring6A-1 C53 0.020
8ZB ring6A-2 C51 0.020
8ZB ring6A-2 N52 0.020
8ZB ring6A-2 C53 0.020
8ZB ring6A-2 N54 0.020
8ZB ring6A-3 N52 0.020
8ZB ring6A-3 C53 0.020
8ZB ring6A-3 N54 0.020
8ZB ring6A-3 C55 0.020
8ZB ring6A-4 C53 0.020
8ZB ring6A-4 N54 0.020
8ZB ring6A-4 C55 0.020
8ZB ring6A-4 C50 0.020
8ZB ring6A-5 N54 0.020
8ZB ring6A-5 C55 0.020
8ZB ring6A-5 C50 0.020
8ZB ring6A-5 C51 0.020
8ZB ring6A-6 C55 0.020
8ZB ring6A-6 C50 0.020
8ZB ring6A-6 C51 0.020
8ZB ring6A-6 N52 0.020
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
8ZB other-001      H1  C24 C28 C31 0.0   1000000.0 10
8ZB csd-sp3sp3-01  C24 C28 C31 C36 180.0 12.3      3
8ZB csd-sp3sp3-02  C28 C31 C36 C40 180.0 12.3      3
8ZB other-002      C31 C36 C40 C39 0.0   1000000.0 10
8ZB other-003      C31 C36 O37 C38 0.0   1000000.0 10
8ZB other-004      N47 C38 C39 C40 0.0   1000000.0 10
8ZB other-005      C39 C38 N47 C49 0.0   1000000.0 10
8ZB other-006      C39 C38 O37 C36 0.0   1000000.0 10
8ZB other-007      C38 C39 C40 C36 0.0   1000000.0 10
8ZB other-008      C38 C39 O41 H46 0.0   1000000.0 10
8ZB other-009      C36 C40 O42 H45 0.0   1000000.0 10
8ZB CONST_ring5A-5 N48 C49 N47 C51 0.0   1000000.0 0
8ZB CONST_ring5A-4 N47 C49 N48 C50 0.0   1000000.0 0
8ZB CONST_ring6A-6 C55 C50 C51 N52 0.0   1000000.0 0
8ZB CONST_ring6A-5 C51 C50 C55 N54 0.0   1000000.0 0
8ZB CONST_ring5A-3 C51 C50 N48 C49 0.0   1000000.0 0
8ZB CONST_ring5A-1 C50 C51 N47 C49 0.0   1000000.0 0
8ZB CONST_ring6A-1 C50 C51 N52 C53 0.0   1000000.0 0
8ZB CONST_ring6A-2 N54 C53 N52 C51 0.0   1000000.0 0
8ZB CONST_ring6A-3 N52 C53 N54 C55 0.0   1000000.0 0
8ZB CONST_ring6A-4 C50 C55 N54 C53 0.0   1000000.0 0
8ZB qmf-01         C50 C55 N57 H59 180.0 1000000.0 2

#
#
"""

def run(prefix = "tst_model_get_vdw_radii"):
  fo = open("%s.pdb"%prefix,"w")
  print(pdb, file=fo)
  fo.close()
  fo = open("%s.cif"%prefix,"w")
  print(cif, file=fo)
  fo.close()
  #
  i = mmtbx.utils.process_command_line_args(args = ["%s.cif"%prefix])
  pdb_inp = iotbx.pdb.input(file_name = "%s.pdb"%prefix)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    restraint_objects = i.cif_objects,
    log               = null_out())
  model.process()
  derived_keys = model.get_vdw_radii()
  # compare with answer
  for atom_name in pdb_inp.atoms().extract_name():
    assert atom_name.strip() in derived_keys

if (__name__ == "__main__"):
  run()
