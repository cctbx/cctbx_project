"""
mmtbx.monomer_library.terminal_nh3 through pdb_interpretation, without
reduce2: an N-terminal NH3 of a peptide monomer whose dictionary does not
fit mod_NH3 (no atom H, or H1-H3 on another atom than N) is restrained with
the dictionary's N-H names plus computed names (generated modification
NH3-<comp>); mod_NH3 is not applied where it would take a carbon H's name;
SCY (H2 on N) and standard residues keep mod_NH3.
"""
from __future__ import absolute_import, division, print_function
import iotbx.cif
import iotbx.pdb
import mmtbx.monomer_library.server
from libtbx.utils import null_out
from mmtbx.monomer_library import pdb_interpretation, terminal_nh3

# 0A0-GLY: HN1, HN2 (dictionary) and HN3 (computed) on N
pdb_0a0 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   0A0 A   1      24.101   8.113 -14.105  1.00 20.00           N
HETATM    2  CA  0A0 A   1      23.417   9.374 -13.816  1.00 20.00           C
HETATM    3  C   0A0 A   1      24.034  10.606 -14.545  1.00 20.00           C
HETATM    4  O   0A0 A   1      23.726  11.736 -14.108  1.00 20.00           O
HETATM    5  CB  0A0 A   1      23.377   9.577 -12.304  1.00 20.00           C
HETATM    6  CG  0A0 A   1      24.721   9.728 -11.566  1.00 20.00           C
HETATM    7  OD1 0A0 A   1      25.780   9.685 -12.231  1.00 20.00           O
HETATM    8  OD2 0A0 A   1      24.641   9.879 -10.329  1.00 20.00           O
HETATM    9  CM  0A0 A   1      21.979   9.248 -14.335  1.00 20.00           C
HETATM   10  HN1 0A0 A   1      24.064   7.949 -14.979  1.00 20.00           H
HETATM   11  HN2 0A0 A   1      23.704   7.452 -13.662  1.00 20.00           H
HETATM   12  HN3 0A0 A   1      24.952   8.172 -13.850  1.00 20.00           H
HETATM   13 HBC1 0A0 A   1      22.843  10.362 -12.106  1.00 20.00           H
HETATM   14 HBC2 0A0 A   1      22.901   8.835 -11.899  1.00 20.00           H
HETATM   15 HMC1 0A0 A   1      21.513   8.517 -13.900  1.00 20.00           H
HETATM   16 HMC2 0A0 A   1      21.475  10.060 -14.170  1.00 20.00           H
HETATM   17 HMC3 0A0 A   1      21.965   9.081 -15.290  1.00 20.00           H
ATOM     18  N   GLY A   2      23.223  11.325 -15.249  1.00  5.27           N
ATOM     19  CA  GLY A   2      23.645  12.533 -15.945  1.00  5.73           C
ATOM     20  C   GLY A   2      23.253  13.729 -15.098  1.00  5.60           C
ATOM     21  O   GLY A   2      22.111  13.831 -14.628  1.00  5.67           O
ATOM     22  H   GLY A   2      22.389  11.143 -15.355  1.00  5.27           H
ATOM     23  HA2 GLY A   2      24.606  12.529 -16.074  1.00  5.73           H
ATOM     24  HA3 GLY A   2      23.210  12.593 -16.810  1.00  5.73           H
END
"""

# 2KZ-GLY: H1, H3 on N (names 2KZ lacks), H2 on CB (the dictionary's)
pdb_2kz = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   2KZ A   1      21.584  21.923  19.934  1.00 20.00           N
HETATM    2  CA  2KZ A   1      20.688  20.877  19.924  1.00 20.00           C
HETATM    3  C   2KZ A   1      21.181  19.479  19.862  1.00 20.00           C
HETATM    4  O   2KZ A   1      20.664  18.687  19.102  1.00 20.00           O
HETATM    5  CB  2KZ A   1      19.356  21.125  19.970  1.00 20.00           C
HETATM    6  CG  2KZ A   1      18.411  20.018  20.124  1.00 20.00           C
HETATM    7  OD1 2KZ A   1      17.085  20.247  20.036  1.00 20.00           O
HETATM    8  OD2 2KZ A   1      18.822  18.892  20.330  1.00 20.00           O
HETATM    9  H1  2KZ A   1      21.225  22.604  20.317  1.00 20.00           H
HETATM   10  H2  2KZ A   1      19.023  21.991  19.905  1.00 20.00           H
HETATM   12  H3  2KZ A   1      22.291  21.691  20.366  1.00 20.00           H
ATOM     13  N   GLY A   2      20.825  18.621  20.776  1.00  5.27           N
ATOM     14  CA  GLY A   2      21.259  17.231  20.813  1.00  5.73           C
ATOM     15  C   GLY A   2      20.158  16.376  20.213  1.00  5.60           C
ATOM     16  O   GLY A   2      18.979  16.512  20.569  1.00  5.67           O
ATOM     17  H   GLY A   2      20.302  18.818  21.430  1.00  5.27           H
ATOM     18  HA2 GLY A   2      22.073  17.118  20.298  1.00  5.73           H
ATOM     19  HA3 GLY A   2      21.423  16.953  21.728  1.00  5.73           H
END
"""

# SCY (GeoStd ideal coordinates): H1 (the dictionary's H), H2 (on N in the
# dictionary), H3
pdb_scy = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   SCY A   1      23.056 -11.316   8.334  1.00 20.00           N
HETATM    2  CA  SCY A   1      24.398 -11.863   8.288  1.00 20.00           C
HETATM    3  CB  SCY A   1      24.407 -13.293   8.816  1.00 20.00           C
HETATM    4  SG  SCY A   1      23.858 -13.382  10.547  1.00 20.00           S
HETATM    5  CD  SCY A   1      22.196 -13.986  10.444  1.00 20.00           C
HETATM    6  OCD SCY A   1      21.678 -14.373   9.429  1.00 20.00           O
HETATM    7  CE  SCY A   1      21.503 -13.943  11.778  1.00 20.00           C
HETATM    8  C   SCY A   1      24.985 -11.809   6.848  1.00 20.00           C
HETATM    9  O   SCY A   1      25.985 -12.521   6.638  1.00 20.00           O
HETATM   10  OXT SCY A   1      24.438 -11.027   6.048  1.00 20.00           O
HETATM   11  H1  SCY A   1      23.014 -10.536   7.691  1.00 20.00           H
HETATM   12  H2  SCY A   1      22.858 -10.953   9.257  1.00 20.00           H
HETATM   13  HA  SCY A   1      25.105 -11.286   8.907  1.00 20.00           H
HETATM   14  HB2 SCY A   1      23.756 -13.940   8.227  1.00 20.00           H
HETATM   15  HB3 SCY A   1      25.413 -13.707   8.795  1.00 20.00           H
HETATM   16  HE1 SCY A   1      20.730 -14.708  11.828  1.00 20.00           H
HETATM   17  HE2 SCY A   1      22.194 -14.069  12.610  1.00 20.00           H
HETATM   18  HE3 SCY A   1      21.024 -12.968  11.890  1.00 20.00           H
HETATM   19  H3  SCY A   1      22.472 -11.954   8.125  1.00 20.00           H
END
"""

# ALA-GLY: H1, H2, H3
pdb_ala = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -1.405   9.032  -6.919  1.00  5.08           N
ATOM      2  CA  ALA A   1      -2.526   8.411  -7.691  1.00  4.71           C
ATOM      3  C   ALA A   1      -2.739   9.202  -8.960  1.00  5.06           C
ATOM      4  O   ALA A   1      -2.863  10.449  -8.918  1.00  5.50           O
ATOM      5  CB  ALA A   1      -3.787   8.434  -6.835  1.00  4.85           C
ATOM      6  H1  ALA A   1      -0.664   9.023  -7.412  1.00  5.08           H
ATOM      7  H2  ALA A   1      -1.615   9.873  -6.715  1.00  5.08           H
ATOM      8  H3  ALA A   1      -1.269   8.570  -6.170  1.00  5.08           H
ATOM      9  HA  ALA A   1      -2.328   7.490  -7.920  1.00  4.71           H
ATOM     10  HB1 ALA A   1      -4.004   9.354  -6.617  1.00  4.85           H
ATOM     11  HB2 ALA A   1      -4.515   8.032  -7.334  1.00  4.85           H
ATOM     12  HB3 ALA A   1      -3.626   7.931  -6.022  1.00  4.85           H
ATOM     13  N   GLY A   2      -2.773   8.494 -10.089  1.00  5.27           N
ATOM     14  CA  GLY A   2      -3.068   9.061 -11.398  1.00  5.73           C
ATOM     15  C   GLY A   2      -4.514   8.748 -11.734  1.00  5.60           C
ATOM     16  O   GLY A   2      -4.966   7.602 -11.604  1.00  5.67           O
ATOM     17  OXT GLY A   2      -5.200   9.700 -12.100  1.00  5.60           O
ATOM     18  H   GLY A   2      -2.622   7.648 -10.122  1.00  5.27           H
ATOM     19  HA2 GLY A   2      -2.939  10.022 -11.386  1.00  5.73           H
ATOM     20  HA3 GLY A   2      -2.488   8.673 -12.072  1.00  5.73           H
END
"""

# NPH-VAL as ready_set (old reduce) writes it: H and H2 on N; NPH's H3
# is on C3
pdb_nph_old_reduce = """\
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N   NPH A   1      23.870  15.268 -50.490  1.00 25.06           N
HETATM    2  CA  NPH A   1      23.515  14.664 -49.210  1.00 22.94           C
HETATM    3  C   NPH A   1      22.061  14.206 -49.245  1.00 23.82           C
HETATM    4  O   NPH A   1      21.158  15.023 -49.423  1.00 24.62           O
HETATM    5  CB  NPH A   1      23.658  15.702 -48.113  1.00 21.44           C
HETATM    6  SG  NPH A   1      25.281  16.410 -47.839  1.00 31.29           S
HETATM    7  CD  NPH A   1      25.498  16.310 -46.059  0.50 33.70           C
HETATM    8  CE  NPH A   1      26.971  16.492 -45.820  0.50 39.92           C
HETATM    9  NZ  NPH A   1      27.411  17.709 -45.649  0.50 43.87           N
HETATM   10  OZ  NPH A   1      27.815  15.348 -45.806  0.50 43.26           O
HETATM   11  C10 NPH A   1      30.543  19.684 -45.064  0.50 54.23           C
HETATM   12  C1A NPH A   1      31.522  18.669 -45.046  0.50 55.19           C
HETATM   13  C2  NPH A   1      33.832  17.924 -44.852  0.50 55.79           C
HETATM   14  C3  NPH A   1      33.439  16.599 -45.020  0.50 53.96           C
HETATM   15  C4  NPH A   1      32.087  16.314 -45.202  0.50 53.40           C
HETATM   16  C4A NPH A   1      31.137  17.332 -45.219  0.50 52.42           C
HETATM   17  C5  NPH A   1      29.800  17.009 -45.411  0.50 51.26           C
HETATM   18  C6  NPH A   1      28.830  18.015 -45.436  0.50 51.14           C
HETATM   19  C6A NPH A   1      29.195  19.344 -45.262  0.50 52.54           C
HETATM   20  C7  NPH A   1      28.201  20.331 -45.288  0.50 53.28           C
HETATM   21  C8  NPH A   1      28.549  21.670 -45.111  0.50 54.83           C
HETATM   22  C9  NPH A   1      29.889  22.011 -44.911  0.50 55.29           C
HETATM   23  N1  NPH A   1      32.879  18.963 -44.861  0.50 57.78           N
HETATM   24  N10 NPH A   1      30.892  21.027 -44.889  0.50 53.87           N
HETATM   25  H   NPH A   1      24.112  16.093 -50.524  1.00 25.06           H
HETATM   26  H2  NPH A   1      23.837  14.799 -51.210  1.00 25.06           H
HETATM   27  H3  NPH A   1      34.108  15.897 -45.009  0.50 53.96           H
HETATM   28  HA  NPH A   1      24.099  13.908 -49.042  1.00 22.94           H
HETATM   29  HB2 NPH A   1      23.355  15.287 -47.290  1.00 21.44           H
HETATM   30  HB3 NPH A   1      23.062  16.441 -48.312  1.00 21.44           H
HETATM   31  HD2 NPH A   1      25.226  15.451 -45.699  0.50 33.70           H
HETATM   32  HD3 NPH A   1      25.032  17.028 -45.603  0.50 33.70           H
HETATM   33  HZ  NPH A   1      26.833  18.345 -45.684  0.50 43.87           H
HETATM   34  H4  NPH A   1      31.747  15.412 -45.311  0.50 53.40           H
HETATM   35  H5  NPH A   1      29.643  16.059 -45.528  0.50 51.26           H
HETATM   36  H7  NPH A   1      27.261  20.125 -45.413  0.50 53.28           H
HETATM   37  H8  NPH A   1      27.874  22.367 -45.116  0.50 54.83           H
HETATM   38  H9  NPH A   1      30.089  22.949 -44.766  0.50 55.29           H
HETATM   39  HC2 NPH A   1      34.778  18.099 -44.728  0.50 55.79           H
ATOM     40  N   VAL A   2      21.831  12.913 -49.037  1.00 28.67           N
ATOM     41  CA  VAL A   2      20.474  12.371 -49.062  1.00 33.96           C
ATOM     42  C   VAL A   2      20.069  11.729 -47.737  1.00 36.17           C
ATOM     43  O   VAL A   2      20.854  11.018 -47.109  1.00 40.12           O
ATOM     44  CB  VAL A   2      20.297  11.336 -50.209  1.00 39.91           C
ATOM     45  CG1 VAL A   2      18.831  10.915 -50.336  1.00 38.68           C
ATOM     46  CG2 VAL A   2      20.791  11.918 -51.524  1.00 43.67           C
ATOM     47  H   VAL A   2      22.442  12.328 -48.879  1.00 28.67           H
ATOM     48  HA  VAL A   2      19.848  13.087 -49.232  1.00 33.96           H
ATOM     49  HB  VAL A   2      20.823  10.546 -50.010  1.00 39.91           H
ATOM     50 HG11 VAL A   2      18.723  10.392 -51.146  1.00 38.68           H
ATOM     51 HG12 VAL A   2      18.585  10.381 -49.564  1.00 38.68           H
ATOM     52 HG13 VAL A   2      18.277  11.710 -50.381  1.00 38.68           H
ATOM     53 HG21 VAL A   2      20.576  11.302 -52.242  1.00 43.67           H
ATOM     54 HG22 VAL A   2      20.352  12.770 -51.676  1.00 43.67           H
ATOM     55 HG23 VAL A   2      21.751  12.045 -51.474  1.00 43.67           H
END
"""

# GeoStd 0A0, 01W and 004 as they were until 2026-11-24 (free forms: 0A0
# with HN1, HN2 on N, 01W with H1N, H2N, H3N, 004 with HN, HNA); GeoStd now
# has one N-H named H (0A0, 01W) or H, HNA (004)
cif_0a0 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 0A0  0A0  '2-methyl-L-aspartic acid '  L-peptide  17  10  .  2012-08-28  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_0A0
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
 0A0  N     N  NH2    0  -0.814  24.1008   8.1126  -14.1049
 0A0  CA    C  CT     0  -0.033  23.4170   9.3744  -13.8160
 0A0  C     C  C      0   0.409  24.0336  10.6057  -14.5448
 0A0  O     O  O      0  -0.620  23.7262  11.7356  -14.1079
 0A0  CB    C  CH2    0  -0.541  23.3767   9.5772  -12.3036
 0A0  CG    C  C      0   0.407  24.7209   9.7277  -11.5659
 0A0  OD1   O  O      0  -0.601  25.7799   9.6850  -12.2307
 0A0  OD2   O  OC    -1  -0.623  24.6411   9.8793  -10.3289
 0A0  OXT   O  OC    -1  -0.617  24.7248  10.3756  -15.5605
 0A0  CM    C  CH3    0  -0.640  21.9789   9.2475  -14.3349
 0A0  HN1   H  HNH2   0   0.311  24.2673   8.0856  -15.1048
 0A0  HN2   H  HNH2   0   0.324  25.0263   8.1971  -13.6983
 0A0  HBC1  H  HCH2   0   0.211  22.7740  10.4605  -12.0793
 0A0  HBC2  H  HCH2   0   0.201  22.8378   8.7366  -11.8484
 0A0  HMC1  H  HCH3   0   0.207  21.4774   8.3908  -13.8786
 0A0  HMC2  H  HCH3   0   0.208  21.3925  10.1397  -14.1125
 0A0  HMC3  H  HCH3   0   0.210  21.9621   9.1017  -15.4186

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 0A0  N   CA    single  1.458  0.019  1.458
 0A0  N   HN1   single  0.890  0.020  1.040
 0A0  N   HN2   single  0.890  0.020  1.040
 0A0  CA  C     single  1.525  0.021  1.525
 0A0  CA  CB    single  1.530  0.020  1.530
 0A0  CA  CM    single  1.534  0.020  1.534
 0A0  C   O     deloc   1.231  0.020  1.231
 0A0  C   OXT   deloc   1.250  0.020  1.250
 0A0  CB  CG    single  1.541  0.020  1.541
 0A0  CB  HBC1  single  0.970  0.020  1.090
 0A0  CB  HBC2  single  0.970  0.020  1.090
 0A0  CG  OD1   deloc   1.251  0.020  1.251
 0A0  CG  OD2   deloc   1.249  0.020  1.249
 0A0  CM  HMC1  single  0.970  0.020  1.090
 0A0  CM  HMC2  single  0.970  0.020  1.090
 0A0  CM  HMC3  single  0.970  0.020  1.090

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 0A0  HN2   N   HN1   104.33  3.000
 0A0  HN2   N   CA    105.98  3.000
 0A0  HN1   N   CA    107.11  3.000
 0A0  CM    CA  CB    108.75  3.000
 0A0  CM    CA  C     106.15  3.000
 0A0  CB    CA  C     110.10  1.900
 0A0  CM    CA  N     107.44  3.000
 0A0  CB    CA  N     110.50  1.700
 0A0  C     CA  N     111.00  2.800
 0A0  OXT   C   O     125.91  3.000
 0A0  OXT   C   CA    116.97  3.000
 0A0  O     C   CA    120.80  1.700
 0A0  HBC2  CB  HBC1  105.26  3.000
 0A0  HBC2  CB  CG    107.74  3.000
 0A0  HBC1  CB  CG    107.70  3.000
 0A0  HBC2  CB  CA    108.79  3.000
 0A0  HBC1  CB  CA    108.99  3.000
 0A0  CG    CB  CA    117.67  3.000
 0A0  OD2   CG  OD1   125.77  3.000
 0A0  OD2   CG  CB    115.49  3.000
 0A0  OD1   CG  CB    118.74  3.000
 0A0  HMC3  CM  HMC2  107.63  3.000
 0A0  HMC3  CM  HMC1  107.59  3.000
 0A0  HMC2  CM  HMC1  108.03  3.000
 0A0  HMC3  CM  CA    111.14  3.000
 0A0  HMC2  CM  CA    111.55  3.000
 0A0  HMC1  CM  CA    110.73  3.000

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
 0A0  Var_01  OD1   CG  CB  CA     0.36  30.0  3
 0A0  Var_02  HBC1  CB  CA  N   -173.92  30.0  3
 0A0  Var_03  HBC2  CB  CA  N    -59.65  30.0  3
 0A0  Var_04  HMC1  CM  CA  N     57.06  30.0  3
 0A0  Var_05  HMC2  CM  CA  N    177.41  30.0  3
 0A0  Var_06  HMC3  CM  CA  N    -62.46  30.0  3
 0A0  Var_07  HN1   N   CA  C    -44.45  30.0  3
 0A0  Var_08  HN2   N   CA  C     66.50  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 0A0  chir_01  CA  C  CB  CM  negativ

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 0A0  plan-1  CA   0.020
 0A0  plan-1  C    0.020
 0A0  plan-1  O    0.020
 0A0  plan-1  OXT  0.020
 0A0  plan-2  CB   0.020
 0A0  plan-2  CG   0.020
 0A0  plan-2  OD1  0.020
 0A0  plan-2  OD2  0.020
"""

cif_01w = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 01W  01W  (2S)-2-ammonio-4-[(2,4-dinitrophenyl)amino]butanoate  L-peptide  32  20  .  2012-08-28  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_01W
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
 01W  O52   O  O     -1  -0.387  -13.9087  38.8704  -21.2752
 01W  N5    N  N      1   0.320  -13.5582  38.9413  -20.1180
 01W  O5    O  O      0  -0.384  -13.9641  39.7780  -19.3418
 01W  C5    C  CR6    0   0.243  -12.6074  37.9788  -19.6293
 01W  C6    C  CR16   0  -0.246  -12.1797  38.0269  -18.2967
 01W  C7    C  CR16   0  -0.356  -11.2891  37.1047  -17.8342
 01W  C4    C  CR16   0  -0.342  -12.1239  37.0105  -20.4753
 01W  C3    C  CR6    0   0.208  -11.2119  36.0724  -20.0136
 01W  N3    N  N      1   0.348  -10.7606  35.1008  -20.9746
 01W  O32   O  O     -1  -0.410   -9.9391  34.2727  -20.6256
 01W  O3    O  O      0  -0.376  -11.2090  35.1382  -22.0965
 01W  C2    C  CR6    0   0.432  -10.7559  36.0747  -18.6583
 01W  N13   N  NH1    0  -0.658   -9.8885  35.1995  -18.1551
 01W  C1    C  CH2    0  -0.325   -9.4693  35.1359  -16.7712
 01W  CB    C  CH2    0  -0.466  -10.5069  34.4330  -15.8965
 01W  CA    C  CH1    0  -0.228  -10.0428  34.2662  -14.4516
 01W  N     N  NT3    1  -0.790  -11.0306  33.4248  -13.7175
 01W  C     C  C      0   0.453   -9.9654  35.6304  -13.7108
 01W  O     O  O      0  -0.580  -10.9605  35.9098  -13.0157
 01W  OXT   O  OC    -1  -0.564   -8.9480  36.2995  -13.9254
 01W  H6    H  HCR6   0   0.307  -12.5465  38.7870  -17.6221
 01W  H7    H  HCR6   0   0.302  -10.9821  37.1722  -16.8013
 01W  H4    H  HCR6   0   0.328  -12.4474  36.9737  -21.5047
 01W  HN13  H  HNH1   0   0.422   -9.5540  34.4766  -18.7714
 01W  H11   H  HCH2   0   0.283   -8.5357  34.5747  -16.7558
 01W  H21   H  HCH2   0   0.269   -9.2370  36.1311  -16.3986
 01W  H1B   H  HCH2   0   0.260  -11.4543  34.9795  -15.9134
 01W  H2B   H  HCH2   0   0.268  -10.7001  33.4477  -16.3305
 01W  HA    H  HCH1   0   0.294   -9.0846  33.7474  -14.4219
 01W  H1N   H  HNT3   0   0.466  -10.6007  32.9349  -12.9334
 01W  H2N   H  HNT3   0   0.447  -11.7269  34.0709  -13.3202
 01W  H3N   H  HNT3   0   0.461  -11.4947  32.7226  -14.2933

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 01W  O3   N3    double    1.209  0.020  1.209
 01W  N3   O32   single    1.218  0.020  1.218
 01W  N3   C3    single    1.439  0.020  1.439
 01W  O52  N5    single    1.211  0.020  1.211
 01W  C4   C3    aromatic  1.387  0.020  1.387
 01W  C4   C5    aromatic  1.374  0.020  1.374
 01W  C3   C2    aromatic  1.430  0.020  1.430
 01W  N5   C5    single    1.439  0.020  1.439
 01W  N5   O5    double    1.211  0.020  1.211
 01W  C5   C6    aromatic  1.400  0.020  1.400
 01W  C2   N13   single    1.331  0.020  1.331
 01W  C2   C7    aromatic  1.423  0.020  1.423
 01W  N13  C1    single    1.447  0.020  1.447
 01W  C6   C7    aromatic  1.363  0.020  1.363
 01W  C1   CB    single    1.528  0.020  1.528
 01W  CB   CA    single    1.521  0.033  1.521
 01W  CA   C     single    1.525  0.021  1.525
 01W  CA   N     single    1.458  0.019  1.458
 01W  C    O     deloc     1.231  0.020  1.231
 01W  C    OXT   deloc     1.237  0.020  1.237
 01W  C6   H6    single    0.930  0.020  1.080
 01W  C7   H7    single    0.930  0.020  1.080
 01W  C4   H4    single    0.930  0.020  1.080
 01W  N13  HN13  single    0.860  0.020  1.020
 01W  C1   H11   single    0.970  0.020  1.090
 01W  C1   H21   single    0.970  0.020  1.090
 01W  CB   H1B   single    0.970  0.020  1.090
 01W  CB   H2B   single    0.970  0.020  1.090
 01W  CA   HA    single    0.970  0.020  0.970
 01W  N    H1N   single    0.890  0.020  1.040
 01W  N    H2N   single    0.890  0.020  1.040
 01W  N    H3N   single    0.890  0.020  1.040

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 01W  C5    N5   O5   117.77  3.000
 01W  O5    N5   O52  123.76  3.000
 01W  C5    N5   O52  118.46  3.000
 01W  C4    C5   C6   120.18  3.000
 01W  C6    C5   N5   120.14  3.000
 01W  C4    C5   N5   119.67  3.000
 01W  H6    C6   C7   119.07  3.000
 01W  H6    C6   C5   120.98  3.000
 01W  C7    C6   C5   119.95  3.000
 01W  H7    C7   C2   119.52  3.000
 01W  H7    C7   C6   117.91  3.000
 01W  C2    C7   C6   122.56  3.000
 01W  H4    C4   C3   119.42  3.000
 01W  H4    C4   C5   120.38  3.000
 01W  C3    C4   C5   120.20  3.000
 01W  C2    C3   N3   122.27  3.000
 01W  C2    C3   C4   121.60  3.000
 01W  N3    C3   C4   116.13  3.000
 01W  O3    N3   O32  122.50  3.000
 01W  O32   N3   C3   118.65  3.000
 01W  O3    N3   C3   118.85  3.000
 01W  N13   C2   C3   124.41  3.000
 01W  N13   C2   C7   120.09  3.000
 01W  C3    C2   C7   115.50  3.000
 01W  HN13  N13  C1   117.22  3.000
 01W  HN13  N13  C2   117.20  3.000
 01W  C1    N13  C2   125.39  3.000
 01W  H21   C1   H11  106.47  3.000
 01W  H21   C1   CB   111.70  3.000
 01W  H11   C1   CB   109.69  3.000
 01W  H21   C1   N13  110.43  3.000
 01W  H11   C1   N13  106.52  3.000
 01W  CB    C1   N13  111.77  3.000
 01W  H2B   CB   H1B  106.91  3.000
 01W  H2B   CB   CA   109.31  3.000
 01W  H1B   CB   CA   109.42  3.000
 01W  H2B   CB   C1   107.88  3.000
 01W  H1B   CB   C1   110.45  3.000
 01W  CA    CB   C1   112.68  3.000
 01W  HA    CA   C    110.00  3.000
 01W  HA    CA   N    110.00  3.000
 01W  C     CA   N    111.00  2.800
 01W  HA    CA   CB   109.00  3.000
 01W  N     CA   CB   110.40  1.500
 01W  C     CA   CB   110.50  1.500
 01W  H3N   N    H2N  110.01  3.000
 01W  H3N   N    H1N  107.18  3.000
 01W  H2N   N    H1N  106.85  3.000
 01W  H3N   N    CA   114.33  3.000
 01W  H2N   N    CA   106.50  3.000
 01W  H1N   N    CA   111.74  3.000
 01W  OXT   C    O    129.26  3.000
 01W  OXT   C    CA   120.80  1.700
 01W  O     C    CA   120.80  1.700

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
 01W  CONST_01  C2    C3   C4   C5     0.00   0.0  0
 01W  CONST_02  C2    C7   C6   C5     0.00   0.0  0
 01W  CONST_03  C3    C4   C5   C6     0.00   0.0  0
 01W  CONST_04  C3    C2   C7   C6     0.00   0.0  0
 01W  CONST_05  C4    C3   C2   C7     0.00   0.0  0
 01W  CONST_06  C4    C5   C6   C7     0.00   0.0  0
 01W  CONST_07  C3    C4   C5   N5   180.00   0.0  0
 01W  CONST_08  C7    C6   C5   N5   180.00   0.0  0
 01W  CONST_09  N3    C3   C4   C5   180.00   0.0  0
 01W  CONST_10  N13   C2   C7   C6   180.00   0.0  0
 01W  CONST_11  N3    C3   C2   C7   180.00   0.0  0
 01W  CONST_12  N13   C2   C3   C4   180.00   0.0  0
 01W  CONST_13  H7    C7   C6   C5   180.00   0.0  0
 01W  CONST_14  H4    C4   C5   C6   180.00   0.0  0
 01W  CONST_15  H6    C6   C5   C4   180.00   0.0  0
 01W  CONST_16  C4    C5   N5   O52    0.00   0.0  0
 01W  CONST_17  C6    C5   N5   O52  180.00   0.0  0
 01W  CONST_18  C4    C5   N5   O5   180.00   0.0  0
 01W  CONST_19  C6    C5   N5   O5     0.00   0.0  0
 01W  CONST_20  C1    N13  C2   C7     0.00   0.0  0
 01W  CONST_21  O3    N3   C3   C4     0.00   0.0  0
 01W  CONST_22  O32   N3   C3   C4   180.00   0.0  0
 01W  CONST_23  C1    N13  C2   C3   180.00   0.0  0
 01W  CONST_24  C2    C3   N3   O32    0.00   0.0  0
 01W  CONST_25  C2    C3   N3   O3   180.00   0.0  0
 01W  CONST_26  HN13  N13  C2   C7   180.00   0.0  0
 01W  Var_01    CB    C1   N13  C2    79.46  30.0  3
 01W  Var_02    CA    CB   C1   N13  175.93  30.0  3
 01W  Var_03    HA    CA   CB   C1   -55.12  30.0  3
 01W  Var_04    H1N   N    CA   CB   153.77  30.0  3
 01W  Var_05    H2N   N    CA   CB   -89.90  30.0  3
 01W  Var_06    H3N   N    CA   CB    31.81  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 01W  chir_01  CA  CB  C  N  negativ

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 01W  plan-1  N5    0.020
 01W  plan-1  C5    0.020
 01W  plan-1  C6    0.020
 01W  plan-1  C7    0.020
 01W  plan-1  C4    0.020
 01W  plan-1  C3    0.020
 01W  plan-1  N3    0.020
 01W  plan-1  C2    0.020
 01W  plan-1  N13   0.020
 01W  plan-1  H6    0.020
 01W  plan-1  H7    0.020
 01W  plan-1  H4    0.020
 01W  plan-2  O52   0.020
 01W  plan-2  N5    0.020
 01W  plan-2  O5    0.020
 01W  plan-2  C5    0.020
 01W  plan-3  C3    0.020
 01W  plan-3  N3    0.020
 01W  plan-3  O32   0.020
 01W  plan-3  O3    0.020
 01W  plan-4  C2    0.020
 01W  plan-4  N13   0.020
 01W  plan-4  C1    0.020
 01W  plan-4  HN13  0.020
 01W  plan-5  CA    0.020
 01W  plan-5  C     0.020
 01W  plan-5  O     0.020
 01W  plan-5  OXT   0.020
"""

cif_004 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 004  004  '(2S)-amino(phenyl)ethanoic acid'  L-peptide  19  11  .  2012-08-28  2025-11-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as PERFECT
Added dihedrals for adding hydrogens : 2023-10-30
Idealise E
;

data_comp_004
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
 004  C    C  C      0   0.408  30.7797   2.4184  52.6942
 004  N    N  NH2    0  -0.804  32.2461   0.6198  53.5997
 004  O    O  O      0  -0.591  31.4760   3.1699  53.4008
 004  CA   C  CH1    0  -0.297  31.2046   0.9146  52.6343
 004  CB   C  CR6    0   0.172  29.9962   0.0351  52.8450
 004  CE   C  CR16   0  -0.277  27.6949  -1.5239  53.2362
 004  CD1  C  CR16   0  -0.242  28.3991  -1.0237  54.3249
 004  CD2  C  CR16   0  -0.243  28.1412  -1.2423  51.9521
 004  CG1  C  CR16   0  -0.293  29.5376  -0.2537  54.1295
 004  CG2  C  CR16   0  -0.326  29.2810  -0.4712  51.7636
 004  OXT  O  OC    -1  -0.593  29.8013   2.7283  51.9926
 004  HN   H  HNH2   0   0.333  32.8126  -0.1535  53.2782
 004  HNA  H  HNH2   0   0.338  32.8455   1.4299  53.6796
 004  HA   H  HCH1   0   0.217  31.5389   0.7611  51.5951
 004  HE   H  HCR6   0   0.244  26.8119  -2.1314  53.3879
 004  HD1  H  HCR6   0   0.242  28.0608  -1.2377  55.3314
 004  HD2  H  HCR6   0   0.243  27.6032  -1.6266  51.0943
 004  HG1  H  HCR6   0   0.235  30.0898   0.1152  54.9841
 004  HG2  H  HCR6   0   0.234  29.6262  -0.2668  50.7561

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 004  C    O    deloc     1.231  0.020  1.231
 004  C    OXT  deloc     1.243  0.019  1.243
 004  N    HN   single    0.890  0.010  1.040
 004  N    HNA  single    0.890  0.010  1.040
 004  CA   C    single    1.525  0.022  1.525
 004  CA   N    single    1.458  0.022  1.458
 004  CA   CB   single    1.521  0.012  1.521
 004  CA   HA   single    0.970  0.010  0.970
 004  CB   CG1  aromatic  1.394  0.010  1.394
 004  CE   CD1  aromatic  1.390  0.010  1.390
 004  CE   HE   single    0.930  0.014  1.080
 004  CD1  HD1  single    0.930  0.014  1.080
 004  CD2  CE   aromatic  1.388  0.014  1.388
 004  CD2  HD2  single    0.930  0.014  1.080
 004  CG1  CD1  aromatic  1.388  0.014  1.388
 004  CG1  HG1  single    0.930  0.010  1.080
 004  CG2  CB   aromatic  1.392  0.010  1.392
 004  CG2  CD2  aromatic  1.389  0.010  1.389
 004  CG2  HG2  single    0.930  0.010  1.080

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 004  OXT  C    CA   120.80  3.000
 004  CA   C    O    120.80  3.000
 004  OXT  C    O    127.60  1.432
 004  HNA  N    HN   107.80  1.000
 004  HNA  N    CA   108.39  1.000
 004  HN   N    CA   110.25  1.000
 004  HA   CA   CB   109.00  3.000
 004  HA   CA   N    110.00  3.000
 004  CB   CA   N    110.40  2.166
 004  HA   CA   C    110.00  3.000
 004  CB   CA   C    110.50  3.000
 004  N    CA   C    111.00  3.000
 004  CG2  CB   CG1  118.13  1.000
 004  CG2  CB   CA   120.99  1.785
 004  CG1  CB   CA   120.84  1.785
 004  HE   CE   CD2  120.38  3.000
 004  CD2  CE   CD1  119.26  1.000
 004  HE   CE   CD1  120.37  3.000
 004  HD1  CD1  CG1  119.77  3.000
 004  CG1  CD1  CE   120.33  1.000
 004  HD1  CD1  CE   119.91  3.000
 004  HD2  CD2  CG2  119.81  3.000
 004  CG2  CD2  CE   120.13  1.000
 004  HD2  CD2  CE   120.06  3.000
 004  HG1  CG1  CD1  119.76  3.000
 004  HG1  CG1  CB   119.28  3.000
 004  CD1  CG1  CB   120.95  1.000
 004  HG2  CG2  CD2  119.47  3.000
 004  HG2  CG2  CB   119.32  3.000
 004  CD2  CG2  CB   121.21  1.000

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
 004  CONST_01  CE   CD1  CG1  CB      0.00   0.0  0
 004  CONST_02  CE   CD2  CG2  CB      0.00   0.0  0
 004  CONST_03  CG2  CD2  CE   CD1     0.00   0.0  0
 004  CONST_04  CG2  CB   CG1  CD1     0.00   0.0  0
 004  CONST_05  CG1  CD1  CE   CD2     0.00   0.0  0
 004  CONST_06  CG1  CB   CG2  CD2     0.00   0.0  0
 004  CONST_07  CD1  CG1  CB   CA    180.00   0.0  0
 004  CONST_08  CD2  CG2  CB   CA    180.00   0.0  0
 004  CONST_09  HD1  CD1  CG1  CB    180.00   0.0  0
 004  CONST_10  HD2  CD2  CG2  CB    180.00   0.0  0
 004  CONST_11  HG1  CG1  CD1  CE    180.00   0.0  0
 004  CONST_12  HG2  CG2  CD2  CE    180.00   0.0  0
 004  CONST_13  HE   CE   CD1  CG1   180.00   0.0  0
 004  Var_01    HA   CA   CB   CG1   164.02  30.0  3
 004  Var_02    HN   N    CA   CB     84.84  30.0  3
 004  Var_03    HNA  N    CA   CB   -157.37  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 004  chir_01  CA  C  N  CB  negativ

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 004  plan-1  CA   0.020
 004  plan-1  CB   0.020
 004  plan-1  CE   0.020
 004  plan-1  CD1  0.020
 004  plan-1  CD2  0.020
 004  plan-1  CG1  0.020
 004  plan-1  CG2  0.020
 004  plan-1  HE   0.020
 004  plan-1  HD1  0.020
 004  plan-1  HD2  0.020
 004  plan-1  HG1  0.020
 004  plan-1  HG2  0.020
 004  plan-2  C    0.020
 004  plan-2  O    0.020
 004  plan-2  CA   0.020
 004  plan-2  OXT  0.020
"""

# GeoStd CTE as it was until 2026-11-25 (HN, HN1 on N; now H, HN1) and 0AA
# (HN1, HA2 on N), whose computed third name HN2 is also a name of
# pdb_interpretation's HN1/HN2 -> mod_NH2 rule
cif_cte = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 CTE  CTE  '7-chloro-L-tryptophan    '  L-peptide  26  16  .  2012-08-28  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND (side chain)
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_CTE
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
 CTE  O    O   O      0  -0.607  2.7665  5.8434  14.0246
 CTE  C    C   C      0   0.397  3.9059  5.5785  14.4579
 CTE  OXT  O   OC    -1  -0.609  4.2013  4.7160  15.3089
 CTE  CA   C   CH1    0  -0.154  5.0781  6.4201  13.8978
 CTE  N    N   NH2    0  -0.808  5.0364  6.4704  12.4391
 CTE  CB   C   CH2    0  -0.529  5.0024  7.8596  14.4237
 CTE  CG   C   CR5    0   0.015  5.0761  7.9706  15.9113
 CTE  CD1  C   CR15   0  -0.085  6.2014  7.9460  16.6860
 CTE  NE1  N   NR15   0  -0.684  5.8720  8.0544  18.0111
 CTE  CE2  C   CR56   0   0.458  4.5228  8.1609  18.1222
 CTE  CD2  C   CR56   0   0.019  3.9764  8.1070  16.8204
 CTE  CE3  C   CR16   0  -0.353  2.5900  8.2011  16.6550
 CTE  CZ3  C   CR16   0  -0.278  1.7917  8.3520  17.7705
 CTE  CH2  C   CR16   0  -0.267  2.3405  8.4058  19.0623
 CTE  CZ2  C   CR6    0  -0.122  3.7036  8.3090  19.2378
 CTE  CL   CL  CL     0  -0.131  4.4161  8.3700  20.8165
 CTE  HA   H   HCH1   0   0.173  6.0140  5.9813  14.2740
 CTE  HN1  H   HNH2   0   0.332  5.8823  6.9179  12.1054
 CTE  HN   H   HNH2   0   0.335  5.0660  5.5299  12.0646
 CTE  HB1  H   HCH2   0   0.231  5.8179  8.4358  13.9742
 CTE  HB2  H   HCH2   0   0.226  4.0741  8.3090  14.0623
 CTE  HD1  H   HCR5   0   0.278  7.2352  7.8596  16.3890
 CTE  HE1  H   HNR5   0   0.405  6.5300  8.0711  18.7725
 CTE  HE3  H   HCR6   0   0.240  2.1531  8.1438  15.6669
 CTE  HZ3  H   HCR6   0   0.251  0.7178  8.4264  17.6605
 CTE  HH2  H   HCR6   0   0.267  1.6933  8.5226  19.9212

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 CTE  O    C    deloc     1.231  0.020  1.231
 CTE  C    OXT  deloc     1.247  0.020  1.247
 CTE  C    CA   single    1.525  0.021  1.525
 CTE  CA   N    single    1.458  0.019  1.458
 CTE  CA   CB   single    1.530  0.020  1.530
 CTE  CA   HA   single    0.970  0.020  0.970
 CTE  N    HN1  single    0.890  0.020  1.040
 CTE  N    HN   single    0.890  0.020  1.040
 CTE  CB   CG   single    1.494  0.020  1.494
 CTE  CB   HB1  single    0.970  0.020  1.090
 CTE  CB   HB2  single    0.970  0.020  1.090
 CTE  CG   CD1  aromatic  1.366  0.020  1.366
 CTE  CG   CD2  aromatic  1.433  0.020  1.433
 CTE  CD1  NE1  aromatic  1.370  0.020  1.370
 CTE  CD1  HD1  single    0.930  0.020  1.080
 CTE  NE1  CE2  aromatic  1.358  0.020  1.358
 CTE  NE1  HE1  single    0.860  0.020  1.020
 CTE  CE2  CD2  aromatic  1.413  0.020  1.413
 CTE  CE2  CZ2  aromatic  1.392  0.020  1.392
 CTE  CD2  CE3  aromatic  1.399  0.020  1.399
 CTE  CE3  CZ3  aromatic  1.380  0.020  1.380
 CTE  CE3  HE3  single    0.930  0.020  1.080
 CTE  CZ3  CH2  aromatic  1.405  0.020  1.405
 CTE  CZ3  HZ3  single    0.930  0.020  1.080
 CTE  CH2  CZ2  aromatic  1.378  0.020  1.378
 CTE  CH2  HH2  single    0.930  0.020  1.080
 CTE  CZ2  CL   single    1.733  0.020  1.733

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 CTE  CA   C    OXT  116.33  3.000
 CTE  CA   C    O    120.80  1.700
 CTE  OXT  C    O    126.88  3.000
 CTE  HA   CA   CB   109.00  3.000
 CTE  HA   CA   N    110.00  3.000
 CTE  CB   CA   N    110.50  1.700
 CTE  HA   CA   C    109.00  3.000
 CTE  CB   CA   C    110.10  1.900
 CTE  N    CA   C    111.00  2.800
 CTE  HN   N    HN1  105.30  3.000
 CTE  HN   N    CA   109.67  3.000
 CTE  HN1  N    CA   108.68  3.000
 CTE  HB2  CB   HB1  106.29  3.000
 CTE  HB2  CB   CG   109.93  3.000
 CTE  HB1  CB   CG   109.44  3.000
 CTE  HB2  CB   CA   108.32  3.000
 CTE  HB1  CB   CA   108.43  3.000
 CTE  CG   CB   CA   114.12  3.000
 CTE  CD2  CG   CD1  105.90  3.000
 CTE  CD2  CG   CB   126.93  3.000
 CTE  CD1  CG   CB   127.16  3.000
 CTE  HD1  CD1  NE1  120.20  3.000
 CTE  HD1  CD1  CG   129.37  3.000
 CTE  NE1  CD1  CG   110.43  3.000
 CTE  HE1  NE1  CE2  125.90  3.000
 CTE  HE1  NE1  CD1  125.17  3.000
 CTE  CE2  NE1  CD1  108.92  3.000
 CTE  CZ2  CE2  CD2  121.00  3.000
 CTE  CZ2  CE2  NE1  131.20  3.000
 CTE  CD2  CE2  NE1  107.81  3.000
 CTE  CE3  CD2  CE2  119.31  3.000
 CTE  CE3  CD2  CG   133.75  3.000
 CTE  CE2  CD2  CG   106.94  3.000
 CTE  HE3  CE3  CZ3  120.69  3.000
 CTE  HE3  CE3  CD2  120.30  3.000
 CTE  CZ3  CE3  CD2  119.00  3.000
 CTE  HZ3  CZ3  CH2  118.60  3.000
 CTE  HZ3  CZ3  CE3  119.96  3.000
 CTE  CH2  CZ3  CE3  121.44  3.000
 CTE  HH2  CH2  CZ2  119.89  3.000
 CTE  HH2  CH2  CZ3  120.04  3.000
 CTE  CZ2  CH2  CZ3  120.07  3.000
 CTE  CL   CZ2  CH2  121.35  3.000
 CTE  CL   CZ2  CE2  119.46  3.000
 CTE  CH2  CZ2  CE2  119.18  3.000

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
 CTE  CONST_01  CE2  NE1  CD1  CG      0.00   0.0  0
 CTE  CONST_02  NE1  CE2  CD2  CG      0.00   0.0  0
 CTE  CONST_03  CZ2  CE2  CD2  CG    180.00   0.0  0
 CTE  CONST_04  CZ3  CE3  CD2  CG    180.00   0.0  0
 CTE  CONST_05  CE2  CD2  CG   CD1     0.00   0.0  0
 CTE  CONST_06  CE3  CD2  CG   CD1   180.00   0.0  0
 CTE  CONST_07  CD2  CE2  NE1  CD1     0.00   0.0  0
 CTE  CONST_08  CZ2  CE2  NE1  CD1   180.00   0.0  0
 CTE  CONST_09  CD2  CG   CD1  NE1     0.00   0.0  0
 CTE  CONST_10  CE3  CD2  CE2  NE1   180.00   0.0  0
 CTE  CONST_11  CH2  CZ2  CE2  NE1   180.00   0.0  0
 CTE  CONST_12  CZ3  CE3  CD2  CE2     0.00   0.0  0
 CTE  CONST_13  CZ3  CH2  CZ2  CE2     0.00   0.0  0
 CTE  CONST_14  CH2  CZ2  CE2  CD2     0.00   0.0  0
 CTE  CONST_15  CH2  CZ3  CE3  CD2     0.00   0.0  0
 CTE  CONST_16  CZ2  CE2  CD2  CE3     0.00   0.0  0
 CTE  CONST_17  CZ2  CH2  CZ3  CE3     0.00   0.0  0
 CTE  CONST_18  NE1  CD1  CG   CB    180.00   0.0  0
 CTE  CONST_19  CE2  CD2  CG   CB    180.00   0.0  0
 CTE  CONST_20  CE3  CD2  CG   CB      0.00   0.0  0
 CTE  CONST_21  CL   CZ2  CE2  NE1     0.00   0.0  0
 CTE  CONST_22  CL   CZ2  CE2  CD2   180.00   0.0  0
 CTE  CONST_23  CL   CZ2  CH2  CZ3   180.00   0.0  0
 CTE  CONST_24  HE1  NE1  CD1  CG    180.00   0.0  0
 CTE  CONST_25  HE3  CE3  CD2  CG      0.00   0.0  0
 CTE  CONST_26  HD1  CD1  NE1  CE2   180.00   0.0  0
 CTE  CONST_27  HH2  CH2  CZ2  CE2   180.00   0.0  0
 CTE  CONST_28  HZ3  CZ3  CE3  CD2   180.00   0.0  0
 CTE  Var_01    CD1  CG   CB   CA     81.30  30.0  2
 CTE  Var_02    HA   CA   CB   CG    -56.94  30.0  2
 CTE  Var_03    HN1  N    CA   C    -172.33  30.0  3
 CTE  Var_04    HN   N    CA   C     -57.70  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 CTE  chir_01  CA  C  N  CB  negativ

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 CTE  plan-1  CB   0.020
 CTE  plan-1  CG   0.020
 CTE  plan-1  CD1  0.020
 CTE  plan-1  NE1  0.020
 CTE  plan-1  CE2  0.020
 CTE  plan-1  CD2  0.020
 CTE  plan-1  CE3  0.020
 CTE  plan-1  CZ3  0.020
 CTE  plan-1  CH2  0.020
 CTE  plan-1  CZ2  0.020
 CTE  plan-1  CL   0.020
 CTE  plan-1  HD1  0.020
 CTE  plan-1  HE1  0.020
 CTE  plan-1  HE3  0.020
 CTE  plan-1  HZ3  0.020
 CTE  plan-1  HH2  0.020
 CTE  plan-2  O    0.020
 CTE  plan-2  C    0.020
 CTE  plan-2  OXT  0.020
 CTE  plan-2  CA   0.020
"""

cif_0aa = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 0AA  0AA  'methyl L-valinate        '  L-peptide  22  9  .  2022-03-10  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_0AA
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
 0AA  N     N  NH2   0  -0.794  -0.8712   1.8113   0.3810
 0AA  CA    C  CH1   0  -0.176  -0.5496   0.3925   0.3672
 0AA  C     C  C     0   0.597   0.8569   0.2736  -0.1892
 0AA  O     O  O     0  -0.456   1.1762   0.6793  -1.2802
 0AA  CB    C  CH1   0  -0.257  -1.5063  -0.5084  -0.4384
 0AA  CG1   C  CH3   0  -0.664  -1.0730  -1.9708  -0.3788
 0AA  CG2   C  CH3   0  -0.661  -2.9371  -0.3542   0.0633
 0AA  O1    O  O2    0  -0.447   1.6975  -0.3166   0.6461
 0AA  C1    C  CH3   0  -0.428   3.0499  -0.4791   0.2138
 0AA  HN1   H  HNH2  0   0.347  -1.1120   2.1256  -0.5518
 0AA  HA2   H  HNH2  0   0.361  -1.6867   1.9704   0.9574
 0AA  HB    H  HCH1  0   0.253  -0.5467   0.0461   1.4044
 0AA  HA    H  HCH1  0   0.226  -1.4719  -0.1775  -1.4830
 0AA  HG11  H  HCH3  0   0.230  -0.0838  -2.1367  -0.8096
 0AA  HG12  H  HCH3  0   0.222  -1.0512  -2.3315   0.6521
 0AA  HG13  H  HCH3  0   0.236  -1.7688  -2.6002  -0.9343
 0AA  HG21  H  HCH3  0   0.226  -3.3004   0.6686  -0.0401
 0AA  HG22  H  HCH3  0   0.235  -3.6130  -0.9951  -0.5040
 0AA  HG23  H  HCH3  0   0.219  -3.0200  -0.6372   1.1157
 0AA  H1C1  H  HCH3  0   0.239   3.1111  -1.1141  -0.6700
 0AA  H1C2  H  HCH3  0   0.240   3.5177   0.4811  -0.0020
 0AA  H1C3  H  HCH3  0   0.252   3.5717  -0.9568   1.0374

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 0AA  N    CA    single  1.458  0.019  1.458
 0AA  N    HN1   single  0.890  0.020  1.040
 0AA  N    HA2   single  0.890  0.020  1.040
 0AA  CA   C     single  1.525  0.021  1.525
 0AA  CA   CB    single  1.540  0.027  1.540
 0AA  CA   HB    single  0.970  0.020  1.090
 0AA  C    O     double  1.231  0.020  1.231
 0AA  C    O1    single  1.324  0.020  1.324
 0AA  CB   CG1   single  1.526  0.020  1.526
 0AA  CB   CG2   single  1.524  0.020  1.524
 0AA  CB   HA    single  0.970  0.020  1.090
 0AA  CG1  HG11  single  0.970  0.020  1.090
 0AA  CG1  HG12  single  0.970  0.020  1.090
 0AA  CG1  HG13  single  0.970  0.020  1.090
 0AA  CG2  HG21  single  0.970  0.020  1.090
 0AA  CG2  HG22  single  0.970  0.020  1.090
 0AA  CG2  HG23  single  0.970  0.020  1.090
 0AA  O1   C1    single  1.429  0.020  1.429
 0AA  C1   H1C1  single  0.970  0.020  1.090
 0AA  C1   H1C2  single  0.970  0.020  1.090
 0AA  C1   H1C3  single  0.970  0.020  1.090

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 0AA  HA2   N    HN1   106.51  3.000
 0AA  HA2   N    CA    109.71  3.000
 0AA  HN1   N    CA    110.26  3.000
 0AA  HB    CA   CB    108.20  3.000
 0AA  HB    CA   C     108.69  3.000
 0AA  CB    CA   C     111.40  1.900
 0AA  HB    CA   N     107.48  3.000
 0AA  CB    CA   N     111.50  1.700
 0AA  C     CA   N     111.00  2.800
 0AA  O1    C    O     123.52  3.000
 0AA  O1    C    CA    113.09  3.000
 0AA  O     C    CA    120.80  1.700
 0AA  HA    CB   CG2   108.22  3.000
 0AA  HA    CB   CG1   108.51  3.000
 0AA  CG2   CB   CG1   110.52  3.000
 0AA  HA    CB   CA    109.00  3.000
 0AA  CG2   CB   CA    110.58  3.000
 0AA  CG1   CB   CA    111.31  3.000
 0AA  HG13  CG1  HG12  107.63  3.000
 0AA  HG13  CG1  HG11  106.83  3.000
 0AA  HG12  CG1  HG11  107.71  3.000
 0AA  HG13  CG1  CB    110.61  3.000
 0AA  HG12  CG1  CB    111.03  3.000
 0AA  HG11  CG1  CB    112.79  3.000
 0AA  HG23  CG2  HG22  107.56  3.000
 0AA  HG23  CG2  HG21  107.99  3.000
 0AA  HG22  CG2  HG21  107.19  3.000
 0AA  HG23  CG2  CB    111.22  3.000
 0AA  HG22  CG2  CB    110.56  3.000
 0AA  HG21  CG2  CB    112.12  3.000
 0AA  C1    O1   C     117.43  3.000
 0AA  H1C3  C1   H1C2  109.37  3.000
 0AA  H1C3  C1   H1C1  109.38  3.000
 0AA  H1C2  C1   H1C1  109.19  3.000
 0AA  H1C3  C1   O1    105.98  3.000
 0AA  H1C2  C1   O1    111.46  3.000
 0AA  H1C1  C1   O1    111.39  3.000

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
 0AA  Var_01  HG11  CG1  CB  CA   -61.63  30.0  3
 0AA  Var_02  HG21  CG2  CB  CA    60.71  30.0  3
 0AA  Var_03  H1C1  C1   O1  C    -62.29  30.0  3
 0AA  Var_04  HA    CB   CA  N     61.37  30.0  3
 0AA  Var_05  HN1   N    CA  C     77.70  30.0  3
 0AA  Var_06  HA2   N    CA  C   -165.30  30.0  3
 0AA  Var_07  HB    CA   C   O   -174.12  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 0AA  chir_01  CA  N  C  CB  positiv

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 0AA  plan-1  CA  0.020
 0AA  plan-1  C   0.020
 0AA  plan-1  O   0.020
 0AA  plan-1  O1  0.020
"""

# GeoStd BB8, MH8 and TY2 as they were until 2026-11-25 (HN1, HA2 on N; now
# H, HA2)
cif_bb8 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 BB8  BB8  (betaS)-beta-hydroxy-L-phenylalanine  L-peptide  23  13  .  2012-08-28  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as PERFECT
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_BB8
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
 BB8  N    N  NH2    0  -0.781  10.0155  -1.6053  33.9253
 BB8  CA   C  CH1    0  -0.281  10.6147  -0.7121  32.9533
 BB8  C    C  C      0   0.440  11.3961  -1.4163  31.8102
 BB8  CB   C  CH1    0  -0.017  11.4935   0.3358  33.6760
 BB8  OB   O  OH1    0  -0.707  12.0423   1.2495  32.7658
 BB8  CG   C  CR6    0   0.159  10.7011   1.0861  34.7172
 BB8  CD1  C  CR16   0  -0.304   9.8146   2.0929  34.3408
 BB8  CE1  C  CR16   0  -0.240   9.0705   2.7747  35.2934
 BB8  CZ   C  CR16   0  -0.273   9.2049   2.4648  36.6416
 BB8  CE2  C  CR16   0  -0.242  10.0917   1.4688  37.0275
 BB8  CD2  C  CR16   0  -0.323  10.8323   0.7883  36.0692
 BB8  O    O  O      0  -0.561  11.2321  -2.6305  31.6499
 BB8  OXT  O  OC    -1  -0.617  12.1297  -0.6705  31.1002
 BB8  HN1  H  HNH2   0   0.339   9.4927  -2.3135  33.4252
 BB8  HA2  H  HNH2   0   0.327  10.7517  -2.1201  34.3964
 BB8  HA   H  HCH1   0   0.234   9.8007  -0.1533  32.4691
 BB8  HB   H  HCH1   0   0.190  12.2889  -0.2233  34.2034
 BB8  HC   H  HOH1   0   0.435  12.2083   0.6926  31.9589
 BB8  HD1  H  HCR6   0   0.241   9.7030   2.3491  33.2945
 BB8  HD2  H  HCR6   0   0.240  11.5266   0.0157  36.3791
 BB8  HE1  H  HCR6   0   0.246   8.3850   3.5534  34.9831
 BB8  HZ   H  HCR6   0   0.248   8.6263   2.9991  37.3844
 BB8  HE2  H  HCR6   0   0.247  10.2092   1.2218  38.0751

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 BB8  N    CA   single    1.458  0.019  1.458
 BB8  CA   C    single    1.525  0.021  1.525
 BB8  CA   CB   single    1.530  0.020  1.530
 BB8  C    OXT  deloc     1.264  0.020  1.264
 BB8  CB   OB   single    1.402  0.020  1.402
 BB8  CB   CG   single    1.508  0.020  1.508
 BB8  CG   CD1  aromatic  1.393  0.020  1.393
 BB8  CG   CD2  aromatic  1.391  0.020  1.391
 BB8  CD1  CE1  aromatic  1.388  0.020  1.388
 BB8  CE1  CZ   aromatic  1.390  0.020  1.390
 BB8  CZ   CE2  aromatic  1.388  0.020  1.388
 BB8  CE2  CD2  aromatic  1.389  0.020  1.389
 BB8  N    HN1  single    0.890  0.020  1.040
 BB8  N    HA2  single    0.890  0.020  1.040
 BB8  CA   HA   single    0.970  0.020  0.970
 BB8  C    O    deloc     1.231  0.020  1.231
 BB8  CB   HB   single    0.970  0.020  1.090
 BB8  OB   HC   single    0.850  0.020  0.980
 BB8  CD1  HD1  single    0.930  0.020  1.080
 BB8  CD2  HD2  single    0.930  0.020  1.080
 BB8  CE1  HE1  single    0.930  0.020  1.080
 BB8  CZ   HZ   single    0.930  0.020  1.080
 BB8  CE2  HE2  single    0.930  0.020  1.080

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 BB8  HA2  N    HN1  104.43  3.000
 BB8  HA2  N    CA   108.91  3.000
 BB8  HN1  N    CA   108.25  3.000
 BB8  HA   CA   CB   109.00  3.000
 BB8  HA   CA   C    109.00  3.000
 BB8  CB   CA   C    110.10  1.900
 BB8  HA   CA   N    110.00  3.000
 BB8  CB   CA   N    110.50  1.700
 BB8  C    CA   N    111.00  2.800
 BB8  OXT  C    O    125.72  3.000
 BB8  O    C    CA   120.80  1.700
 BB8  OXT  C    CA   115.96  3.000
 BB8  HB   CB   CG   107.46  3.000
 BB8  HB   CB   OB   110.95  3.000
 BB8  CG   CB   OB   109.25  3.000
 BB8  HB   CB   CA   106.79  3.000
 BB8  CG   CB   CA   111.16  3.000
 BB8  OB   CB   CA   111.14  3.000
 BB8  HC   OB   CB   103.14  3.000
 BB8  CD2  CG   CD1  118.52  3.000
 BB8  CD2  CG   CB   121.00  3.000
 BB8  CD1  CG   CB   120.48  3.000
 BB8  HD1  CD1  CE1  119.45  3.000
 BB8  HD1  CD1  CG   119.84  3.000
 BB8  CE1  CD1  CG   120.71  3.000
 BB8  HE1  CE1  CZ   119.97  3.000
 BB8  HE1  CE1  CD1  119.74  3.000
 BB8  CZ   CE1  CD1  120.29  3.000
 BB8  HZ   CZ   CE2  120.31  3.000
 BB8  HZ   CZ   CE1  120.26  3.000
 BB8  CE2  CZ   CE1  119.43  3.000
 BB8  HE2  CE2  CD2  119.86  3.000
 BB8  HE2  CE2  CZ   120.13  3.000
 BB8  CD2  CE2  CZ   120.02  3.000
 BB8  HD2  CD2  CE2  119.56  3.000
 BB8  HD2  CD2  CG   119.41  3.000
 BB8  CE2  CD2  CG   121.03  3.000

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
 BB8  CONST_01  CZ   CE1  CD1  CG      0.00   0.0  0
 BB8  CONST_02  CZ   CE2  CD2  CG      0.00   0.0  0
 BB8  CONST_03  CE2  CD2  CG   CD1     0.00   0.0  0
 BB8  CONST_04  CE2  CZ   CE1  CD1     0.00   0.0  0
 BB8  CONST_05  CD2  CG   CD1  CE1     0.00   0.0  0
 BB8  CONST_06  CD2  CE2  CZ   CE1     0.00   0.0  0
 BB8  CONST_07  CE1  CD1  CG   CB    180.00   0.0  0
 BB8  CONST_08  CE2  CD2  CG   CB    180.00   0.0  0
 BB8  CONST_09  HE1  CE1  CD1  CG    180.00   0.0  0
 BB8  CONST_10  HE2  CE2  CD2  CG    180.00   0.0  0
 BB8  CONST_11  HD2  CD2  CG   CD1   180.00   0.0  0
 BB8  CONST_12  HZ   CZ   CE1  CD1   180.00   0.0  0
 BB8  CONST_13  HD1  CD1  CE1  CZ    180.00   0.0  0
 BB8  Var_01    HB   CB   CG   CD1   166.52  30.0  2
 BB8  Var_02    HA   CA   CB   CG     60.58  30.0  2
 BB8  Var_03    HC   OB   CB   CG   -157.85  30.0  3
 BB8  Var_04    HN1  N    CA   C     -50.04  30.0  3
 BB8  Var_05    HA2  N    CA   C      62.94  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 BB8  chir_01  CA  N   C   CB  positiv
 BB8  chir_02  CB  OB  CA  CG  positiv

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 BB8  plan-1  CB   0.020
 BB8  plan-1  CG   0.020
 BB8  plan-1  CD1  0.020
 BB8  plan-1  CE1  0.020
 BB8  plan-1  CZ   0.020
 BB8  plan-1  CE2  0.020
 BB8  plan-1  CD2  0.020
 BB8  plan-1  HD1  0.020
 BB8  plan-1  HD2  0.020
 BB8  plan-1  HE1  0.020
 BB8  plan-1  HZ   0.020
 BB8  plan-1  HE2  0.020
 BB8  plan-2  CA   0.020
 BB8  plan-2  C    0.020
 BB8  plan-2  O    0.020
 BB8  plan-2  OXT  0.020
"""

cif_mh8 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 MH8  MH8  '(2S)-2-amino-2-methylhept-6-enoic acid'  L-peptide  25  11  .  2016-04-26  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as PERFECT
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_MH8
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
 MH8  N     N  NH2    0  -0.799  92.3434  85.9697   8.1019
 MH8  CA    C  CT     0  -0.028  91.5256  86.6361   7.0836
 MH8  CB1   C  CH3    0  -0.635  90.7790  85.5528   6.3144
 MH8  CB2   C  CH2    0  -0.408  90.5228  87.6078   7.7406
 MH8  CG    C  CH2    0  -0.445  89.6482  88.4327   6.8024
 MH8  CD    C  CH2    0  -0.485  88.6511  89.3214   7.5530
 MH8  CE    C  C1     0  -0.148  87.6193  88.5418   8.3046
 MH8  CZ    C  C2     0  -0.537  87.4038  88.6133   9.6136
 MH8  C     C  C      0   0.427  92.4875  87.4530   6.1746
 MH8  O     O  O      0  -0.597  92.3256  87.3843   4.9411
 MH8  OXT   O  OC    -1  -0.612  93.3400  88.1437   6.7729
 MH8  HN1   H  HNH2   0   0.325  92.9195  86.6713   8.5515
 MH8  HA2   H  HNH2   0   0.321  91.7359  85.6002   8.8259
 MH8  HB11  H  HCH3   0   0.213  90.1508  85.9625   5.5262
 MH8  HB12  H  HCH3   0   0.212  91.4740  84.8493   5.8542
 MH8  HB13  H  HCH3   0   0.198  90.1320  84.9859   6.9901
 MH8  HB21  H  HCH2   0   0.204  89.8777  87.0211   8.4059
 MH8  HB22  H  HCH2   0   0.205  91.0792  88.2941   8.3890
 MH8  HG1   H  HCH2   0   0.225  90.2726  89.0714   6.1732
 MH8  HG2   H  HCH2   0   0.215  89.0972  87.7776   6.1214
 MH8  HD1   H  HCH2   0   0.223  89.1853  89.9852   8.2393
 MH8  HD2   H  HCH2   0   0.237  88.1449  89.9652   6.8259
 MH8  HE    H  H      0   0.233  87.0105  87.8658   7.7081
 MH8  HZ1   H  H      0   0.230  86.6368  88.0227  10.0980
 MH8  HZ2   H  H      0   0.226  87.9884  89.2681  10.2500

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 MH8  N    CA    single  1.458  0.019  1.458
 MH8  CA   CB1   single  1.524  0.020  1.524
 MH8  CA   CB2   single  1.543  0.020  1.543
 MH8  CA   C     single  1.525  0.021  1.525
 MH8  CB2  CG    single  1.525  0.020  1.525
 MH8  CG   CD    single  1.532  0.020  1.532
 MH8  CD   CE    single  1.496  0.020  1.496
 MH8  CE   CZ    double  1.329  0.020  1.329
 MH8  C    O     deloc   1.231  0.020  1.231
 MH8  C    OXT   deloc   1.250  0.020  1.250
 MH8  N    HN1   single  0.890  0.020  1.040
 MH8  N    HA2   single  0.890  0.020  1.040
 MH8  CB1  HB11  single  0.970  0.020  1.090
 MH8  CB1  HB12  single  0.970  0.020  1.090
 MH8  CB1  HB13  single  0.970  0.020  1.090
 MH8  CB2  HB21  single  0.970  0.020  1.090
 MH8  CB2  HB22  single  0.970  0.020  1.090
 MH8  CG   HG1   single  0.970  0.020  1.090
 MH8  CG   HG2   single  0.970  0.020  1.090
 MH8  CD   HD1   single  0.970  0.020  1.090
 MH8  CD   HD2   single  0.970  0.020  1.090
 MH8  CE   HE    single  0.930  0.020  1.080
 MH8  CZ   HZ1   single  0.930  0.020  1.080
 MH8  CZ   HZ2   single  0.930  0.020  1.080

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 MH8  HA2   N    HN1   106.02  3.000
 MH8  HA2   N    CA    109.09  3.000
 MH8  HN1   N    CA    108.10  3.000
 MH8  C     CA   CB2   108.66  3.000
 MH8  C     CA   CB1   112.42  3.000
 MH8  CB2   CA   CB1   110.13  3.000
 MH8  C     CA   N     111.00  2.800
 MH8  CB2   CA   N     110.67  3.000
 MH8  CB1   CA   N     107.50  3.000
 MH8  HB13  CB1  HB12  107.65  3.000
 MH8  HB13  CB1  HB11  107.53  3.000
 MH8  HB12  CB1  HB11  107.77  3.000
 MH8  HB13  CB1  CA    110.27  3.000
 MH8  HB12  CB1  CA    111.05  3.000
 MH8  HB11  CB1  CA    112.38  3.000
 MH8  HB22  CB2  HB21  105.94  3.000
 MH8  HB22  CB2  CG    108.45  3.000
 MH8  HB21  CB2  CG    108.97  3.000
 MH8  HB22  CB2  CA    108.44  3.000
 MH8  HB21  CB2  CA    107.67  3.000
 MH8  CG    CB2  CA    116.83  3.000
 MH8  HG2   CG   HG1   106.22  3.000
 MH8  HG2   CG   CD    108.94  3.000
 MH8  HG1   CG   CD    108.36  3.000
 MH8  HG2   CG   CB2   110.37  3.000
 MH8  HG1   CG   CB2   110.04  3.000
 MH8  CD    CG   CB2   112.68  3.000
 MH8  HD2   CD   HD1   106.59  3.000
 MH8  HD2   CD   CE    108.73  3.000
 MH8  HD1   CD   CE    109.75  3.000
 MH8  HD2   CD   CG    108.45  3.000
 MH8  HD1   CD   CG    109.97  3.000
 MH8  CE    CD   CG    113.12  3.000
 MH8  HE    CE   CZ    118.88  3.000
 MH8  HE    CE   CD    115.75  3.000
 MH8  CZ    CE   CD    125.37  3.000
 MH8  HZ2   CZ   HZ1   116.67  3.000
 MH8  HZ2   CZ   CE    121.56  3.000
 MH8  HZ1   CZ   CE    121.77  3.000
 MH8  OXT   C    O     126.37  3.000
 MH8  OXT   C    CA    120.80  1.700
 MH8  O     C    CA    120.80  1.700

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
 MH8  Var_01  CD    CG   CB2  CA    177.69  30.0  3
 MH8  Var_02  CE    CD   CG   CB2   -65.87  30.0  3
 MH8  Var_03  CZ    CE   CD   CG    121.53  30.0  3
 MH8  Var_04  HZ1   CZ   CE   CD    179.81   0.0  2
 MH8  Var_05  HB11  CB1  CA   N    -179.42  30.0  3
 MH8  Var_06  HB12  CB1  CA   N     -58.61  30.0  3
 MH8  Var_07  HB13  CB1  CA   N      60.63  30.0  3
 MH8  Var_08  HB21  CB2  CA   N     -60.17  30.0  3
 MH8  Var_09  HB22  CB2  CA   N      54.05  30.0  3
 MH8  Var_10  HN1   N    CA   CB1   173.94  30.0  3
 MH8  Var_11  HA2   N    CA   CB1   -71.21  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 MH8  chir_01  CA  CB1  CB2  C  positiv

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 MH8  plan-1  CA   0.020
 MH8  plan-1  C    0.020
 MH8  plan-1  O    0.020
 MH8  plan-1  OXT  0.020
 MH8  plan-2  CD   0.020
 MH8  plan-2  CE   0.020
 MH8  plan-2  CZ   0.020
 MH8  plan-2  HE   0.020
 MH8  plan-2  HZ1  0.020
 MH8  plan-2  HZ2  0.020
"""

cif_ty2 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 TY2  TY2  '3-amino-L-tyrosine       '  L-peptide  25  14  .  2012-08-28  2023-12-07
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as PERFECT (side chain)
Added dihedrals for adding hydrogens : 2023-12-07
;

data_comp_TY2
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
 TY2  N     N  NH2    0  -0.798  1.5024  27.8639  41.2758
 TY2  CA    C  CH1    0  -0.173  2.0607  29.2050  41.4084
 TY2  C     C  C      0   0.412  1.3654  30.1315  40.3857
 TY2  O     O  O      0  -0.614  1.4469  29.7658  39.1926
 TY2  CB    C  CH2    0  -0.527  3.5883  29.2342  41.2098
 TY2  CG    C  CR6    0   0.125  4.2085  30.5845  41.4318
 TY2  CD1   C  CR16   0  -0.340  4.5815  30.9959  42.7050
 TY2  CD2   C  CR16   0  -0.372  4.4023  31.4666  40.3701
 TY2  CE1   C  CR16   0  -0.364  5.1339  32.2571  42.9049
 TY2  CE2   C  CR6    0   0.234  4.9495  32.7353  40.5447
 TY2  CZ    C  CR6    0   0.280  5.3183  33.1259  41.8429
 TY2  OH    O  OH1    0  -0.641  5.8553  34.3665  41.9734
 TY2  NE2   N  NH2    0  -0.852  5.0812  33.6258  39.4878
 TY2  OXT   O  OC    -1  -0.603  0.7943  31.1498  40.8175
 TY2  HN1   H  HNH2   0   0.325  1.6250  27.5651  40.3151
 TY2  HA2   H  HNH2   0   0.331  2.0362  27.2134  41.8411
 TY2  HB    H  HCH1   0   0.194  1.8441  29.5639  42.4196
 TY2  HA1   H  HCH2   0   0.225  4.0349  28.5130  41.9020
 TY2  HB2   H  HCH2   0   0.221  3.8072  28.8770  40.1998
 TY2  HA    H  HCR6   0   0.234  4.4523  30.3316  43.5512
 TY2  HC    H  HCR6   0   0.236  4.1145  31.1609  39.3694
 TY2  HD    H  HCR6   0   0.256  5.4275  32.5657  43.9031
 TY2  HE21  H  HNH2   0   0.370  5.8529  34.2698  39.5677
 TY2  HE22  H  HNH2   0   0.371  5.0827  33.1841  38.5814
 TY2  HH    H  HOH1   0   0.470  6.1164  34.5229  42.8874

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 TY2  N    CA    single    1.458  0.019  1.458
 TY2  CA   C     single    1.525  0.021  1.525
 TY2  CA   CB    single    1.530  0.020  1.530
 TY2  C    O     deloc     1.231  0.020  1.231
 TY2  C    OXT   deloc     1.245  0.020  1.245
 TY2  CB   CG    single    1.502  0.020  1.502
 TY2  CG   CD1   aromatic  1.389  0.020  1.389
 TY2  CG   CD2   aromatic  1.394  0.020  1.394
 TY2  CD1  CE1   aromatic  1.391  0.020  1.391
 TY2  CD2  CE2   aromatic  1.393  0.020  1.393
 TY2  CE1  CZ    aromatic  1.384  0.020  1.384
 TY2  CE2  CZ    aromatic  1.405  0.020  1.405
 TY2  CE2  NE2   single    1.388  0.020  1.388
 TY2  CZ   OH    single    1.358  0.020  1.358
 TY2  N    HN1   single    0.890  0.020  1.040
 TY2  N    HA2   single    0.890  0.020  1.040
 TY2  CA   HB    single    0.970  0.020  1.090
 TY2  CB   HA1   single    0.970  0.020  1.090
 TY2  CB   HB2   single    0.970  0.020  1.090
 TY2  CD1  HA    single    0.930  0.020  1.080
 TY2  CD2  HC    single    0.930  0.020  1.080
 TY2  CE1  HD    single    0.930  0.020  1.080
 TY2  NE2  HE21  single    0.890  0.020  1.040
 TY2  NE2  HE22  single    0.890  0.020  1.040
 TY2  OH   HH    single    0.850  0.020  0.980

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 TY2  HA2   N    HN1   106.01  3.000
 TY2  HA2   N    CA    109.74  3.000
 TY2  HN1   N    CA    108.11  3.000
 TY2  HB    CA   CB    108.00  3.000
 TY2  HB    CA   C     109.01  3.000
 TY2  CB    CA   C     110.10  1.900
 TY2  HB    CA   N     108.04  3.000
 TY2  CB    CA   N     110.50  1.700
 TY2  C     CA   N     111.00  2.800
 TY2  OXT   C    O     126.88  3.000
 TY2  OXT   C    CA    117.86  3.000
 TY2  O     C    CA    120.80  1.700
 TY2  HB2   CB   HA1   106.69  3.000
 TY2  HB2   CB   CG    110.33  3.000
 TY2  HA1   CB   CG    109.28  3.000
 TY2  HB2   CB   CA    108.14  3.000
 TY2  HA1   CB   CA    108.09  3.000
 TY2  CG    CB   CA    114.03  3.000
 TY2  CD2   CG   CD1   118.25  3.000
 TY2  CD2   CG   CB    120.91  3.000
 TY2  CD1   CG   CB    120.83  3.000
 TY2  HA    CD1  CE1   119.39  3.000
 TY2  HA    CD1  CG    120.15  3.000
 TY2  CE1   CD1  CG    120.45  3.000
 TY2  HC    CD2  CE2   118.45  3.000
 TY2  HC    CD2  CG    119.16  3.000
 TY2  CE2   CD2  CG    122.39  3.000
 TY2  HD    CE1  CZ    119.41  3.000
 TY2  HD    CE1  CD1   119.83  3.000
 TY2  CZ    CE1  CD1   120.77  3.000
 TY2  NE2   CE2  CZ    120.01  3.000
 TY2  NE2   CE2  CD2   121.75  3.000
 TY2  CZ    CE2  CD2   118.18  3.000
 TY2  OH    CZ   CE2   116.52  3.000
 TY2  OH    CZ   CE1   123.52  3.000
 TY2  CE2   CZ   CE1   119.96  3.000
 TY2  HH    OH   CZ    110.27  3.000
 TY2  HE22  NE2  HE21  110.48  3.000
 TY2  HE22  NE2  CE2   113.80  3.000
 TY2  HE21  NE2  CE2   114.96  3.000

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
 TY2  CONST_01  CZ    CE1  CD1  CG     0.00   0.0  0
 TY2  CONST_02  CZ    CE2  CD2  CG     0.00   0.0  0
 TY2  CONST_03  CE2   CD2  CG   CD1    0.00   0.0  0
 TY2  CONST_04  CE2   CZ   CE1  CD1    0.00   0.0  0
 TY2  CONST_05  CE1   CD1  CG   CD2    0.00   0.0  0
 TY2  CONST_06  CE1   CZ   CE2  CD2    0.00   0.0  0
 TY2  CONST_07  CE1   CD1  CG   CB   180.00   0.0  0
 TY2  CONST_08  CE2   CD2  CG   CB   180.00   0.0  0
 TY2  CONST_09  NE2   CE2  CD2  CG   180.00   0.0  0
 TY2  CONST_10  OH    CZ   CE1  CD1  180.00   0.0  0
 TY2  CONST_11  OH    CZ   CE2  CD2  180.00   0.0  0
 TY2  CONST_12  NE2   CE2  CZ   CE1  180.00   0.0  0
 TY2  CONST_13  HD    CE1  CD1  CG   180.00   0.0  0
 TY2  CONST_14  HC    CD2  CG   CD1  180.00   0.0  0
 TY2  CONST_15  HA    CD1  CG   CD2  180.00   0.0  0
 TY2  CONST_16  HE21  NE2  CE2  CD2  180.00   0.0  0
 TY2  CONST_17  HE22  NE2  CE2  CD2    0.00   0.0  0
 TY2  Var_01    CD1   CG   CB   CA    85.98  30.0  2
 TY2  Var_02    HH    OH   CZ   CE1    2.56  30.0  2
 TY2  Var_03    HB    CA   CB   CG   -56.82  30.0  2
 TY2  Var_04    HN1   N    CA   C     56.65  30.0  3
 TY2  Var_05    HA2   N    CA   C    171.83  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 TY2  chir_01  CA  N  C  CB  positiv

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 TY2  plan-1  CB    0.020
 TY2  plan-1  CG    0.020
 TY2  plan-1  CD1   0.020
 TY2  plan-1  CD2   0.020
 TY2  plan-1  CE1   0.020
 TY2  plan-1  CE2   0.020
 TY2  plan-1  CZ    0.020
 TY2  plan-1  OH    0.020
 TY2  plan-1  NE2   0.020
 TY2  plan-1  HA    0.020
 TY2  plan-1  HC    0.020
 TY2  plan-1  HD    0.020
 TY2  plan-2  CE2   0.020
 TY2  plan-2  NE2   0.020
 TY2  plan-2  HE21  0.020
 TY2  plan-2  HE22  0.020
 TY2  plan-3  CA    0.020
 TY2  plan-3  C     0.020
 TY2  plan-3  O     0.020
 TY2  plan-3  OXT   0.020
"""

OLD_DICTIONARIES = (('0A0', cif_0a0), ('01W', cif_01w), ('004', cif_004),
  ('CTE', cif_cte), ('0AA', cif_0aa), ('BB8', cif_bb8), ('MH8', cif_mh8),
  ('TY2', cif_ty2))

def old_residue_names(pdb_str):
  """the residue names in pdb_str that have an old dictionary here"""
  found = set(l[17:20] for l in pdb_str.splitlines()
    if l.startswith(('ATOM', 'HETATM')))
  return [comp for comp, text in OLD_DICTIONARIES if comp in found]

def old_dictionaries(pdb_str):
  """restraint objects of the old dictionaries of the residues in pdb_str;
  None if it has none of them"""
  names = old_residue_names(pdb_str)
  result = [(comp + '_old.cif', iotbx.cif.reader(input_string=text).model())
    for comp, text in OLD_DICTIONARIES if comp in names]
  return result or None

def server():
  """a monomer library server with the old dictionaries"""
  srv = mmtbx.monomer_library.server.server()
  for comp, text in OLD_DICTIONARIES:
    srv.process_cif_object(cif_object=iotbx.cif.reader(
      input_string=text).model(), file_name=comp + '_old.cif')
  return srv

_srv = []
def interpret(pdb_str):
  """residue names with their modifications; {atom name of residue 1:
  sorted bonded atom names}; {(name, name): N-CA ideal}"""
  if not _srv:
    _srv.append(server())
    _srv.append(mmtbx.monomer_library.server.ener_lib())
  p = pdb_interpretation.process(mon_lib_srv=_srv[0], ener_lib=_srv[1],
    pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_str),
    log=null_out(),
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  acp = p.all_chain_proxies
  names = set(mm.residue_name for mm in acp.all_monomer_mappings)
  atoms = acp.pdb_hierarchy.atoms()
  grm = p.geometry_restraints_manager(show_energies=False)
  bonded, n_ca = {}, []
  for prx in grm.pair_proxies(sites_cart=atoms.extract_xyz()
      ).bond_proxies.simple:
    a, b = [atoms[i] for i in prx.i_seqs]
    if a.parent().parent().resseq.strip() != '1': continue
    if b.parent().parent().resseq.strip() != '1': continue
    bonded.setdefault(a.name.strip(), []).append(b.name.strip())
    bonded.setdefault(b.name.strip(), []).append(a.name.strip())
    if sorted([a.name.strip(), b.name.strip()]) == ['CA', 'N']:
      n_ca.append(prx.distance_ideal)
  first = [a for a in atoms if a.parent().parent().resseq.strip() == '1']
  h = [a.name.strip() for a in first if a.element.strip() == 'H']
  return names, dict((k, sorted(v)) for k, v in bonded.items()), h, n_ca

def test_0a0():
  # the generated modification: every H on N restrained, N-CA as mod_NH3
  # (CA sp3)
  names, bonded, h, n_ca = interpret(pdb_0a0)
  assert '0A0%NH3-0A0' in names, names
  for x in h:
    assert x in bonded, (x, bonded)
  for x in ('HN1', 'HN2', 'HN3'):
    assert bonded[x] == ['N'], (x, bonded[x])
  assert len(n_ca) == 1 and abs(n_ca[0] - 1.491) < 1.e-6, n_ca

def test_2kz_carbon_h():
  # H1 and H3 on N: mod_NH3 is not applied (it would add H2 to N and take
  # CB's H2); H2 stays restrained to CB; H1 and H3 are unknown
  names, bonded, h, _ = interpret(pdb_2kz)
  assert '2KZ' in names and not [x for x in names if '%' in x and
    x.startswith('2KZ')], names
  assert bonded['H2'] == ['CB'], bonded['H2']
  assert 'H1' not in bonded and 'H3' not in bonded, bonded

def test_scy():
  # H2 on N in the dictionary: mod_NH3 still applies and replaces it. Without
  # OXT in GeoStd's SCY (phenix-project/geostd#1) the model's OXT also takes
  # mod_COO, as for a standard residue
  names, bonded, h, _ = interpret(pdb_scy)
  assert 'SCY%NH3' in names or 'SCY%COO%NH3' in names, names
  for x in ('H1', 'H2', 'H3'):
    assert bonded.get(x) == ['N'], (x, bonded.get(x))

def test_standard_terminus():
  names, bonded, h, n_ca = interpret(pdb_ala)
  assert 'ALA%NH3' in names, names
  for x in h:
    assert x in bonded, (x, bonded)
  assert len(n_ca) == 1 and abs(n_ca[0] - 1.491) < 1.e-6, n_ca

def test_partial_and_sp3():
  # NPH needs two computed names (H2, HN3): one of them present is enough,
  # as for mod_NH3's H1-H3; CA sp3 by neighbour count, also for
  # mon_lib-format entries whose bonds are all "coval" (MSE)
  from mmtbx.monomer_library import terminal_nh3
  interpret(pdb_ala)
  srv = _srv[0]
  nph = srv.get_comp_comp_id_direct('NPH')
  assert terminal_nh3.extra_h_names(nph) == ['H2', 'HN3']
  for unexpected in (['H2'], ['HN3'], ['H2', 'HN3']):
    assert terminal_nh3.terminal_mod(nph, unexpected) is not None, \
      unexpected
  assert terminal_nh3.terminal_mod(nph, ['HXX']) is None
  for comp, sp3 in (('0A0', True), ('MSE', True), ('2KZ', False)):
    assert terminal_nh3.ca_sp3(srv.get_comp_comp_id_direct(comp)) == sp3, \
      comp

def n_h_values(srv, ener_lib):
  """{H name: (neutron ideal, esd)} of the N-H bonds of residue 1"""
  params = pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  p = pdb_interpretation.process(mon_lib_srv=srv, ener_lib=ener_lib,
    params=params, pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_0a0),
    log=null_out(),
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  atoms = p.all_chain_proxies.pdb_hierarchy.atoms()
  grm = p.geometry_restraints_manager(show_energies=False)
  result = {}
  for prx in grm.pair_proxies(sites_cart=atoms.extract_xyz()
      ).bond_proxies.simple:
    a, b = sorted([atoms[i] for i in prx.i_seqs],
      key=lambda x: x.name.strip() != 'N')
    if (a.parent().parent().resseq.strip() == '1' and a.name.strip() == 'N'
        and b.name.strip().startswith('HN')):
      result[b.name.strip()] = (round(prx.distance_ideal, 3),
        round(prx.weight**-0.5, 3))
  return result

def test_replacement_cif():
  # a CIF for 0A0 loaded into the same server after a first interpretation
  # (N-H esd 0.03, neutron 1.08): the generated HN3 takes the new values
  # with the dictionary's HN1 and HN2, as in a fresh server
  srv = server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  assert n_h_values(srv, ener_lib) == dict((h, (1.04, 0.02))
    for h in ('HN1', 'HN2', 'HN3')), n_h_values(srv, ener_lib)
  text = cif_0a0
  old = 'single  0.890  0.020  1.040'
  assert text.count(old) == 2, text.count(old)
  cif = iotbx.cif.reader(input_string=text.replace(old,
    'single  0.890  0.030  1.080')).model()
  srv.process_cif_object(cif_object=cif, file_name='0A0_replacement.cif')
  values = n_h_values(srv, ener_lib)
  assert values == dict((h, (1.08, 0.03)) for h in ('HN1', 'HN2', 'HN3')), \
    values

def mod_n_h(srv, comp):
  """{N-H name: (x-ray, esd, neutron)} of comp with its generated
  modification applied"""
  mon = srv.get_comp_comp_id_direct(comp)
  mod = terminal_nh3.nh3_mod(mon, terminal_nh3.extra_h_names(mon),
    srv.mod_mod_id_dict['NH3'])
  new = mon.apply_mod(mod)
  result = {}
  for b in new.bond_list:
    ids = set([b.atom_id_1, b.atom_id_2])
    if 'N' in ids and (ids - set(['N'])) & set(
        terminal_nh3.n_h_names(mon) + terminal_nh3.extra_h_names(mon)):
      result[list(ids - set(['N']))[0]] = (b.value_dist, b.value_dist_esd,
        b.value_dist_neutron)
  return result

def test_n_h_values():
  # one value set for the three N-H: the dictionary's own where its N-H
  # agree (0A0), else mod_NH3's (TSY and 2KZ: 0.86 and 0.89; NPH: no
  # neutron value)
  srv = server()
  own = mod_n_h(srv, '0A0')
  assert sorted(own) == ['HN1', 'HN2', 'HN3'], own
  assert set(own.values()) == set([(0.89, 0.02, 1.04)]), own
  assert terminal_nh3.n_h_values(srv.get_comp_comp_id_direct('0A0'))[1]
  for comp in ('TSY', 'NPH', '2KZ'):
    values = mod_n_h(srv, comp)
    assert len(values) == 3, (comp, values)
    assert set(values.values()) == set([(0.89, 0.02, 1.04)]), (comp, values)
  mon = srv.get_comp_comp_id_direct('TSY')
  assert terminal_nh3.n_h_values(mon)[1] is False
  # one nuclear length in both columns is the dictionary's own; an x-ray
  # length repeated (NPH 0.86/0.86) is not
  mon = srv.get_comp_comp_id_direct('0A0')
  for b in mon.bond_list:
    if 'N' in (b.atom_id_1, b.atom_id_2) and set([b.atom_id_1,
        b.atom_id_2]) & set(terminal_nh3.n_h_names(mon)):
      b.value_dist = b.value_dist_neutron = 1.03
  assert terminal_nh3.n_h_values(mon) == ((1.03, 0.02, 1.03), True)
  for b in mon.bond_list:
    if 'N' in (b.atom_id_1, b.atom_id_2) and set([b.atom_id_1,
        b.atom_id_2]) & set(terminal_nh3.n_h_names(mon)):
      b.value_dist = b.value_dist_neutron = 0.86
  assert terminal_nh3.n_h_values(mon) == ((0.89, 0.02, 1.04), False)

def test_imine_and_bond_order():
  # an N=CA imine is not an amine terminus (HHI, ITR); 2KZ's N-CA (single,
  # 1.402, sp2 CA) is. Without a bond order (coval) the ideal length
  # decides: N_CA_SINGLE_MIN
  srv = mmtbx.monomer_library.server.server()
  for comp in ('HHI', 'ITR'):
    mon = srv.get_comp_comp_id_direct(comp)
    assert not terminal_nh3.n_ca_single(mon), comp
    assert not terminal_nh3.applies(mon), comp
  mon = srv.get_comp_comp_id_direct('2KZ')
  assert terminal_nh3.n_ca_single(mon) and terminal_nh3.applies(mon)
  mon = srv.get_comp_comp_id_direct('0A0')
  n_ca = [b for b in mon.bond_list if set([b.atom_id_1, b.atom_id_2]) ==
    set(['N', 'CA'])][0]
  n_ca.type = 'coval'
  for d, expected in ((1.30, False), (terminal_nh3.N_CA_SINGLE_MIN, True),
                      (1.458, True)):
    n_ca.value_dist = d
    assert terminal_nh3.n_ca_single(mon) == expected, (d, expected)
  n_ca.value_dist = None
  assert not terminal_nh3.n_ca_single(mon)

def test_computed_rule_name():
  # 0AA: the computed HN2 is also a name of the HN1/HN2 -> mod_NH2 rule;
  # NH3-0AA restrains it, never mod_NH2
  srv = server()
  mon = srv.get_comp_comp_id_direct('0AA')
  assert terminal_nh3.extra_h_names(mon) == ['HN2']
  nh3 = srv.mod_mod_id_dict['NH3']
  assert terminal_nh3.computed_names_mod(mon, ['HN2'], nh3) is not None
  assert terminal_nh3.computed_names_mod(mon, ['DN2'], nh3) is not None
  # a supplied rule name that is not computed: the rules as before
  assert terminal_nh3.computed_names_mod(mon, ['HN2', 'H1'], nh3) is None
  assert terminal_nh3.computed_names_mod(mon, ['OXT'], nh3) is None

def test_h_series():
  # dictionary N-H named H: the missing H continue the series H, H2, H3
  # (CCD, old reduce) where free, then the fallback names. NPH H3 is on
  # C3: H, H2, HN3. Old reduce's H2 is restrained to N, never by mod_NH3
  # (it would take C3's H3); as D (all_d), D and D2. 0G6 (H1 on another
  # atom): H2, H3
  interpret(pdb_ala)
  srv = _srv[0]
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('NPH')) == ['H2', 'HN3']
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('0G6')) == ['H2', 'H3']
  # 2KZ (H and H4 on N) is not a lone H: H5, as before
  assert terminal_nh3.extra_h_names(
    srv.get_comp_comp_id_direct('2KZ')) == ['H5']
  d = pdb_nph_old_reduce.replace(' H   NPH', ' D   NPH').replace(
    ' H2  NPH', ' D2  NPH')
  d = '\n'.join(l[:76] + ' D' + l[78:] if l[12:16] in (' D  ', ' D2 ')
    else l for l in d.split('\n'))
  for pdb_str, h, h2 in ((pdb_nph_old_reduce, 'H', 'H2'), (d, 'D', 'D2')):
    names, bonded, hs, _ = interpret(pdb_str)
    # every atom typed (old reduce's H2 stopped model.process())
    p = pdb_interpretation.process(mon_lib_srv=_srv[0], ener_lib=_srv[1],
      pdb_inp=iotbx.pdb.input(source_info=None, lines=pdb_str),
      log=null_out(),
      substitute_non_crystallographic_unit_cell_if_necessary=True)
    assert p.all_chain_proxies.fatal_problems_message() is None, \
      p.all_chain_proxies.fatal_problems_message()
    assert 'NPH%NH3-NPH' in names, names
    assert not [x for x in names if x.startswith('NPH%NH3') and
      x != 'NPH%NH3-NPH'], names
    for x in hs:
      assert x in bonded, (x, bonded)
    assert bonded[h] == ['N'] and bonded[h2] == ['N'], bonded
    assert bonded['H3'] == ['C3'], bonded['H3']

if __name__ == '__main__':
  test_0a0()
  test_2kz_carbon_h()
  test_scy()
  test_standard_terminus()
  test_partial_and_sp3()
  test_replacement_cif()
  test_n_h_values()
  test_imine_and_bond_order()
  test_computed_rule_name()
  test_h_series()
  print('OK')
