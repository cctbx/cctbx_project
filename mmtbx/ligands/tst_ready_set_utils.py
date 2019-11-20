from __future__ import division
import sys

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

def main(switch=None):
  if switch=='ASP':
    tst_adding_side_chain_acid_hydrogen_atoms(switch)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
