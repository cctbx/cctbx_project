from __future__ import absolute_import, division, print_function
import sys
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file
from iotbx.pdb import pdb_input

from mmtbx.hydrogens.specialised_hydrogen_atoms import add_side_chain_acid_hydrogens
from mmtbx.hydrogens.specialised_hydrogen_atoms import add_disulfur_hydrogen_atoms

pdb_strings = {'ASP' :'''
CRYST1  112.354   50.667   75.061  90.00 127.09  90.00 C 1 2 1
SCALE1      0.008900  0.000000  0.006729        0.00000
SCALE2      0.000000  0.019737  0.000000        0.00000
SCALE3      0.000000  0.000000  0.016701        0.00000
ATOM      1  OE1 GLU A 249      23.084 -14.146  55.991  1.00 28.06           O
ATOM      2  OD2 ASP A 258      20.555 -12.607  61.969  1.00 37.35           O
ATOM      3  CA  ASN A 269      25.385 -16.651  56.738  1.00 15.81           C
ATOM      4  C   ASN A 269      24.892 -16.548  58.188  1.00 15.44           C
ATOM      5  O   ASN A 269      25.355 -17.281  59.072  1.00 15.27           O
ATOM      6  ND2 ASN A 269      23.311 -19.122  57.209  1.00 18.88           N
ATOM      7  N   ASP A 270      23.952 -15.610  58.412  1.00 18.53           N
ATOM      8  CA  ASP A 270      23.251 -15.460  59.662  1.00 19.81           C
ATOM      9  C   ASP A 270      22.320 -16.673  59.865  1.00 17.50           C
ATOM     10  O   ASP A 270      21.839 -17.297  58.893  1.00 18.62           O
ATOM     11  CB  ASP A 270      22.308 -14.226  59.624  1.00 25.44           C
ATOM     12  CG  ASP A 270      22.950 -13.023  59.939  1.00 30.80           C
ATOM     13  OD1 ASP A 270      24.170 -12.930  59.713  1.00 37.22           O
ATOM     14  OD2 ASP A 270      22.211 -12.121  60.434  1.00 33.61           O
ATOM     15  N   ALA A 271      22.030 -16.956  61.120  1.00 15.85           N
ATOM     16  CA  ALA A 271      21.096 -17.989  61.489  1.00 16.38           C
''',
  'GLU' : '''
CRYST1   39.183   52.539  207.533  90.00  90.00  90.00 P 21 21 21
SCALE1      0.025521  0.000000  0.000000        0.00000
SCALE2      0.000000  0.019033  0.000000        0.00000
SCALE3      0.000000  0.000000  0.004819        0.00000
ATOM      3  N   TRP A 194      -7.810  10.259  62.633  1.00  9.35           N
ATOM      4  CA  TRP A 194      -8.280  11.585  62.268  1.00  8.84           C
ATOM      5  C   TRP A 194      -9.778  11.593  61.945  1.00  8.77           C
ATOM      6  O   TRP A 194     -10.269  10.714  61.243  1.00  8.79           O
ATOM      7  CB  TRP A 194      -7.513  12.112  61.045  1.00  9.05           C
ATOM      8  CG  TRP A 194      -6.007  12.290  61.215  1.00  8.60           C
ATOM      9  N   GLU A 195     -10.509  12.592  62.436  1.00  8.58           N
ATOM     10  CA  GLU A 195     -11.870  12.816  61.951  1.00  8.31           C
ATOM     11  C   GLU A 195     -11.832  13.095  60.469  1.00  7.49           C
ATOM     12  O   GLU A 195     -12.710  12.683  59.735  1.00  7.68           O
ATOM     13  CB  GLU A 195     -12.550  13.986  62.659  1.00  8.68           C
ATOM     14  CG  GLU A 195     -12.799  13.775  64.152  1.00 11.46           C
ATOM     15  CD  GLU A 195     -14.045  14.499  64.657  1.00 12.48           C
ATOM     16  OE1 GLU A 195     -14.622  15.293  63.883  1.00 11.73           O
ATOM     17  OE2 GLU A 195     -14.441  14.257  65.826  1.00 12.64           O
ATOM     18  N   ARG A 196     -10.811  13.811  60.021  1.00  7.15           N
ATOM     19  CA  ARG A 196     -10.703  14.129  58.605  1.00  6.38           C
ATOM     20  C   ARG A 196      -9.350  13.639  58.107  1.00  6.27           C
ATOM     21  O   ARG A 196      -8.343  14.350  58.179  1.00  6.00           O
ATOM     22  CB  ARG A 196     -10.981  15.619  58.367  1.00  6.18           C
ATOM     23  CG  ARG A 196     -12.144  16.088  59.223  1.00  4.12           C
ATOM     24  CD  ARG A 196     -12.782  17.328  58.735  1.00  3.90           C
ATOM     25  N   PRO A 197      -9.326  12.382  57.633  1.00  6.31           N
ATOM     26  CG  PRO A 197     -10.015  10.265  56.861  1.00  5.41           C
ATOM     27  CD  PRO A 197     -10.537  11.559  57.405  1.00  5.91           C
TER
HETATM   49  O5  CIT A 500     -11.052  16.238  62.773  0.75 26.13           O
HETATM   52  O   HOH A 512     -14.767  12.104  68.190  1.00 29.69           O
''',
  'CYS_1' : '''
CRYST1   61.314   37.470   44.915  90.00  94.55  90.00 P 1 21 1
SCALE1      0.016309  0.000000  0.001298        0.00000
SCALE2      0.000000  0.026688  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022335        0.00000
ATOM      8  N   CYS A  26      21.129  -6.109  34.429  1.00  4.21           N
ATOM      9  CA  CYS A  26      20.061  -5.459  33.689  1.00  4.17           C
ATOM     10  C   CYS A  26      18.827  -6.331  33.559  1.00  4.48           C
ATOM     11  O   CYS A  26      17.694  -5.822  33.620  1.00  4.91           O
ATOM     12  CB  CYS A  26      20.554  -5.076  32.290  1.00  4.37           C
ATOM     13  SG  CYS A  26      21.811  -3.780  32.291  1.00  4.75           S
ATOM     14  H   CYS A  26      21.889  -6.315  34.036  1.00  5.05           H
ATOM     15  HA  CYS A  26      19.806  -4.622  34.173  1.00  5.01           H
ATOM     16  HB2 CYS A  26      20.927  -5.881  31.852  1.00  5.24           H
ATOM     17  HB3 CYS A  26      19.782  -4.771  31.752  1.00  5.24           H
ATOM     18  N   ASN A  27      19.008  -7.633  33.365  1.00  4.35           N
ATOM     19  CA  ASN A  27      17.845  -8.529  33.284  1.00  4.69           C
ATOM     20  C   ASN A  27      16.993  -8.399  34.537  1.00  5.07           C
ATOM     28  CB  CYS A  84      20.151  -1.930  30.554  1.00  5.45           C
ATOM     29  SG  CYS A  84      20.664  -2.084  32.297  1.00  5.40           S
''',
  'CYS_2' : '''
CRYST1   61.314   37.470   44.915  90.00  94.55  90.00 P 1 21 1
SCALE1      0.016309  0.000000  0.001298        0.00000
SCALE2      0.000000  0.026688  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022335        0.00000
ATOM      8  N   CYS A  26      21.129  -6.109  34.429  1.00  4.21           N
ATOM      9  CA  CYS A  26      20.061  -5.459  33.689  1.00  4.17           C
ATOM     10  C   CYS A  26      18.827  -6.331  33.559  1.00  4.48           C
ATOM     11  O   CYS A  26      17.694  -5.822  33.620  1.00  4.91           O
ATOM     12  CB  CYS A  26      20.554  -5.076  32.290  1.00  4.37           C
ATOM     13  SG  CYS A  26      21.811  -3.780  32.291  1.00  4.75           S
ATOM     14  H   CYS A  26      21.889  -6.315  34.036  1.00  5.05           H
ATOM     15  HA  CYS A  26      19.806  -4.622  34.173  1.00  5.01           H
ATOM     16  HB2 CYS A  26      20.927  -5.881  31.852  1.00  5.24           H
ATOM     17  HB3 CYS A  26      19.782  -4.771  31.752  1.00  5.24           H
''',
}
for key, item in list(pdb_strings.items()):
  pdb_strings[key+'_D'] = item
geo_strings = {'ASP' :
  {0 : '''nonbonded pdb=" OD2 ASP A 258 "
          pdb=" HD2 ASP A 270 "
   model   vdw
   1.612 1.970''',
   1 : '''nonbonded pdb=" OD2 ASP A 258 "
          pdb=" HD2 ASP A 270 "
   model   vdw
   2.750 1.970''',
   2 : '''nonbonded pdb=" N   ASP A 270 "
          pdb=" HD2 ASP A 270 "
   model   vdw
   2.262 2.650''',
   3 : '''nonbonded pdb=" OD1 ASP A 270 "
          pdb=" HD2 ASP A 270 "
   model   vdw
   2.475 2.056''',
  },
               'ASP_D' :
  {
  0 : '''nonbonded pdb=" OD2 ASP A 258 "
          pdb=" DD2 ASP A 270 "
   model   vdw
   1.583 1.970''',
  1 : '''nonbonded pdb=" OD2 ASP A 258 "
          pdb=" OD2 ASP A 270 "
   model   vdw
   2.310 3.040''',
  2 : '''nonbonded pdb=" N   ASP A 270 "
          pdb=" DD2 ASP A 270 "
   model   vdw
   2.228 2.650''',
  3 : '''nonbonded pdb=" N   ASP A 270 "
          pdb=" DD2 ASP A 270 "
   model   vdw
   3.904 2.650''',
  },
               'GLU' :
  {
  0 : '''nonbonded pdb=" HE2 GLU A 195 "
          pdb=" O   HOH A 512 "
   model   vdw
   2.511 1.970''',
  1 : '''nonbonded pdb=" HE2 GLU A 195 "
          pdb=" O   HOH A 512 "
   model   vdw
   3.316 1.970''',
  2 : '''nonbonded pdb=" HE2 GLU A 195 "
          pdb=" O5  CIT A 500 "
   model   vdw
   3.338 1.970''',
  3 : '''nonbonded pdb=" HE2 GLU A 195 "
          pdb=" O5  CIT A 500 "
   model   vdw
   4.592 1.970''',
  },
               'CYS_2' : '''bond pdb=" SG  CYS A  26 "
     pdb=" HG  CYS A  26 "
  ideal  model  delta    sigma   weight residual
  1.200  1.200  0.000 2.00e-02 2.50e+03 3.33e-04''',
               'CYS_2_D' : '''bond pdb=" SG  CYS A  26 "
     pdb=" DG  CYS A  26 "
  ideal  model  delta    sigma   weight residual
  1.200  1.200  0.000 2.00e-02 2.50e+03 3.33e-04''',
}

def tst_adding_side_chain_acid_hydrogen_atoms(switch):
  element='H'
  if switch.endswith('_D'): element='D'
  for i in range(4):
    fn = 'tst_ready_hydrogens_%s_%d.pdb' % (switch, i)
    with open(fn, 'w') as f:
      f.write(pdb_strings[switch])
    pdb_inp = pdb_input(fn)
    hierarchy = pdb_inp.construct_hierarchy()
    add_side_chain_acid_hydrogens(hierarchy, configuration_index=i, element=element)
    hierarchy.write_pdb_file(fn.replace('.pdb', '_updated.pdb'))
    cmd='phenix.pdb_interpretation const_shrink_donor_acceptor=0.6 %s flip_sym=0 write_geo=True' % fn.replace('.pdb', '_updated.pdb')
    print(cmd)
    easy_run.go(cmd)
    assert_lines_in_file(file_name='%s.geo' % fn.replace('.pdb', '_updated.pdb'),
                         lines = geo_strings[switch][i])

def tst_adding_disulfur_hydrogen_atoms(switch):
  element='H'
  if switch.endswith('_D'): element='D'
  fn = 'tst_ready_hydrogens_%s.pdb' % (switch)
  with open(fn, 'w') as f:
    f.write(pdb_strings[switch])
  pdb_inp = pdb_input(fn)
  hierarchy = pdb_inp.construct_hierarchy()
  from mmtbx.conformation_dependent_library.testing_utils import get_geometry_restraints_manager
  grm = get_geometry_restraints_manager(pdb_filename=fn)
  add_disulfur_hydrogen_atoms(grm, hierarchy, element=element)
  hierarchy.write_pdb_file(fn.replace('.pdb', '_updated.pdb'))
  if switch=='CYS_1':
    with open(fn.replace('.pdb', '_updated.pdb'), 'r') as f:
      lines = f.read()
    assert lines.find('HG')==-1
  else:
    cmd='phenix.pdb_interpretation const_shrink_donor_acceptor=0.6 %s write_geo=True' % fn.replace('.pdb', '_updated.pdb')
    print(cmd)
    easy_run.go(cmd)
    assert_lines_in_file(file_name='%s.geo' % fn.replace('.pdb', '_updated.pdb'),
                         lines = geo_strings[switch])

def switcher(switch):
  if switch in ['ASP', 'GLU', 'ASP_D']:
    tst_adding_side_chain_acid_hydrogen_atoms(switch)
  elif switch in ['CYS_1', 'CYS_2', 'CYS_2_D']:
    tst_adding_disulfur_hydrogen_atoms(switch)

def main(switch=None):
  if switch is None:
    for switch in ['ASP', 'GLU', 'CYS_1', 'CYS_2', 'ASP_D', 'CYS_2_D']:
      switcher(switch)
  else:
    switcher(switch)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
