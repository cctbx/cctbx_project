from __future__ import division
import sys
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

from mmtbx.hydrogens.specialised_hydrogen_atoms import add_side_chain_acid_hydrogens

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
}
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
  }
}

def tst_adding_side_chain_acid_hydrogen_atoms(switch):
  from iotbx.pdb import pdb_input
  for i in range(4):
    fn = 'tst_ready_hydrogens_%s_%d.pdb' % (switch, i)
    f=file(fn, 'wb')
    f.write(pdb_strings[switch])
    del f
    pdb_inp = pdb_input(fn)
    hierarchy = pdb_inp.construct_hierarchy()
    add_side_chain_acid_hydrogens(hierarchy, configuration_index=i)
    hierarchy.write_pdb_file(fn.replace('.pdb', '_updated.pdb'))
    cmd='phenix.pdb_interpretation %s write_geo=True' % fn.replace('.pdb', '_updated.pdb')
    print(cmd)
    easy_run.go(cmd)
    # print('%s.geo' % fn.replace('.pdb', '_updated.pdb'))
    # print(geo_strings[switch][i])
    assert_lines_in_file(file_name='%s.geo' % fn.replace('.pdb', '_updated.pdb'),
                         lines = geo_strings[switch][i])

def switcher(switch):
  if switch in ['ASP', 'GLU']:
    tst_adding_side_chain_acid_hydrogen_atoms(switch)

def main(switch=None):
  if switch is None:
    for switch in ['ASP', 'GLU']:
      tst_adding_side_chain_acid_hydrogen_atoms(switch)
  else:
    tst_adding_side_chain_acid_hydrogen_atoms(switch)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
