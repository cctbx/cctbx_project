from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons

pdb_str = '''
CRYST1   62.350   62.350   62.350  90.00  90.00  90.00 I 41 3 2
SCALE1      0.016038  0.000000  0.000000        0.00000
SCALE2      0.000000  0.016038  0.000000        0.00000
SCALE3      0.000000  0.000000  0.016038        0.00000
HETATM   42  O1  SO4 A  13      31.243  38.502  17.238  1.00 29.92           O
HETATM   43  O2  SO4 A  13      30.616  40.133  15.527  1.00 29.92           O
HETATM   44  O3  SO4 A  13      31.158  37.816  14.905  1.00 29.92           O
HETATM   45  O4  SO4 A  13      32.916  39.343  15.640  1.00 29.92           O
HETATM   46  S   SO4 A  13      31.443  38.969  15.855  1.00 29.92           S
HETATM   47  C   ACE B   0      27.108  39.352  15.639  1.00 29.92           C
HETATM   48  O   ACE B   0      26.447  40.060  14.845  1.00 29.92           O
HETATM   49  CH3 ACE B   0      27.538  37.960  15.306  1.00 29.92           C
HETATM   50  H1  ACE B   0      28.506  37.905  15.332  1.00 29.92           H
HETATM   51  H2  ACE B   0      27.163  37.342  15.952  1.00 29.92           H
HETATM   52  H3  ACE B   0      27.226  37.728  14.417  1.00 29.92           H
ATOM     53  N   GLU C   3      28.216  35.663  23.151  1.00 29.92           N
ATOM     54  CA  GLU C   3      27.676  36.900  22.573  1.00 29.92           C
ATOM     55  C   GLU C   3      26.228  37.153  22.978  1.00 29.92           C
ATOM     56  O   GLU C   3      25.741  38.287  22.876  1.00 29.92           O
ATOM     57  CB  GLU C   3      27.655  36.779  21.044  1.00 29.92           C
ATOM     58  CG  GLU C   3      28.212  37.900  20.176  1.00 29.92           C
ATOM     59  CD  GLU C   3      29.424  37.457  19.404  1.00 29.92           C
ATOM     60  OE1 GLU C   3      30.236  36.902  20.182  1.00 29.92           O
ATOM     61  OE2 GLU C   3      29.571  37.577  18.208  1.00 29.92           O
ATOM     62  OXT GLU C   3      25.539  36.228  23.408  1.00 29.92           O
ATOM     63  H   GLU C   3      28.143  34.981  22.632  1.00 29.92           H
ATOM     64  H2  GLU C   3      29.091  35.772  23.330  1.00 29.92           H
ATOM     65  H3  GLU C   3      27.776  35.470  23.912  1.00 29.92           H
ATOM     66  HA  GLU C   3      28.243  37.619  22.894  1.00 29.92           H
ATOM     67  HB2 GLU C   3      28.138  35.970  20.813  1.00 29.92           H
ATOM     68  HB3 GLU C   3      26.732  36.643  20.779  1.00 29.92           H
ATOM     69  HG2 GLU C   3      27.527  38.202  19.559  1.00 29.92           H
ATOM     70  HG3 GLU C   3      28.444  38.658  20.735  1.00 29.92           H
'''

def main():
  model_filename = 'tst_electron_00.pdb'
  f=open(model_filename, 'w')
  f.write(pdb_str)
  del f
  rc = run_program(program_class=electrons.Program,
                   args=[model_filename],
                   logger=null_out(),
                  )
  print(rc)
  assert rc.total_charge
  print('OK')

if __name__ == '__main__':
  main()
