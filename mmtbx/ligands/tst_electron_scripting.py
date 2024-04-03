from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons

pdb_1x24_modified = '''
CRYST1  146.891  146.891  146.891  90.00  90.00  90.00 I 21 3
ATOM      1  N   PRO A   9      39.761-121.668 105.645  1.00 65.50           N
ATOM      2  CA  PRO A   9      38.280-121.485 105.595  1.00 65.50           C
ATOM      3  C   PRO A   9      37.853-120.009 105.641  1.00 65.50           C
ATOM      4  O   PRO A   9      37.565-119.454 106.709  1.00 65.50           O
ATOM      5  CB  PRO A   9      37.692-122.267 106.767  1.00 82.07           C
ATOM      6  CG  PRO A   9      38.901-122.331 107.732  1.00 82.07           C
ATOM      7  CD  PRO A   9      40.132-122.487 106.816  1.00 82.07           C
ATOM      8  N   VAL A  10      37.821-119.386 104.464  1.00 67.31           N
ATOM      9  CA  VAL A  10      37.429-117.984 104.324  1.00 67.31           C
ATOM     10  C   VAL A  10      35.942-117.869 103.987  1.00 67.31           C
ATOM     11  O   VAL A  10      35.428-118.537 103.086  1.00 67.31           O
ATOM     12  CB  VAL A  10      38.256-117.261 103.215  1.00 55.89           C
ATOM     13  CG1 VAL A  10      39.729-117.340 103.529  1.00 55.89           C
ATOM     14  CG2 VAL A  10      37.996-117.887 101.864  1.00 55.89           C
ATOM     15  N   GLU A  11      35.257-117.017 104.732  1.00 38.26           N
ATOM     16  CA  GLU A  11      33.832-116.787 104.541  1.00 38.26           C
ATOM     17  C   GLU A  11      33.661-115.357 104.025  1.00 38.26           C
ATOM     18  O   GLU A  11      34.303-114.430 104.515  1.00 38.26           O
ATOM     19  CB  GLU A  11      33.113-116.993 105.880  1.00 47.28           C
ATOM     20  CG  GLU A  11      31.647-116.714 105.867  1.00 47.28           C
ATOM     21  CD  GLU A  11      30.927-117.240 107.101  1.00 47.28           C
ATOM     22  OE1 GLU A  11      31.427-117.062 108.233  1.00 47.28           O
ATOM     23  OE2 GLU A  11      29.834-117.825 106.935  1.00 47.28           O
'''

pdb_1x24_modified_complete = '''
CRYST1  146.891  146.891  146.891  90.00  90.00  90.00 I 21 3
ATOM      1  N   PRO A   9      39.761  25.223 105.645  1.00 65.50           N
ATOM      2  CA  PRO A   9      38.280  25.406 105.595  1.00 65.50           C
ATOM      3  C   PRO A   9      37.853  26.882 105.641  1.00 65.50           C
ATOM      4  O   PRO A   9      37.565  27.437 106.709  1.00 65.50           O
ATOM      5  CB  PRO A   9      37.692  24.624 106.767  1.00 82.07           C
ATOM      6  CG  PRO A   9      38.901  24.560 107.732  1.00 82.07           C
ATOM      7  CD  PRO A   9      40.132  24.404 106.816  1.00 82.07           C
ATOM      8  H2  PRO A   9      40.220  26.153 105.717  0.00 65.50           H
ATOM      9  H3  PRO A   9      40.078  24.746 104.777  0.00 65.50           H
ATOM     10  HA  PRO A   9      37.908  24.935 104.685  0.00 65.50           H
ATOM     11  HB2 PRO A   9      36.851  25.160 107.207  0.00 82.07           H
ATOM     12  HB3 PRO A   9      37.375  23.631 106.449  0.00 82.07           H
ATOM     13  HG2 PRO A   9      38.962  25.480 108.314  0.00 82.07           H
ATOM     14  HG3 PRO A   9      38.800  23.704 108.399  0.00 82.07           H
ATOM     15  HD2 PRO A   9      41.030  24.792 107.297  0.00 82.07           H
ATOM     16  HD3 PRO A   9      40.281  23.362 106.534  0.00 82.07           H
ATOM     17  N   VAL A  10      37.821  27.505 104.464  1.00 67.31           N
ATOM     18  CA  VAL A  10      37.429  28.907 104.324  1.00 67.31           C
ATOM     19  C   VAL A  10      35.942  29.022 103.987  1.00 67.31           C
ATOM     20  O   VAL A  10      35.428  28.354 103.086  1.00 67.31           O
ATOM     21  CB  VAL A  10      38.256  29.630 103.215  1.00 55.89           C
ATOM     22  CG1 VAL A  10      39.729  29.551 103.529  1.00 55.89           C
ATOM     23  CG2 VAL A  10      37.996  29.004 101.864  1.00 55.89           C
ATOM     24  H   VAL A  10      38.064  27.061 103.578  0.00 67.31           H
ATOM     25  HA  VAL A  10      37.621  29.395 105.280  0.00 67.31           H
ATOM     26  HB  VAL A  10      37.952  30.676 103.182  0.00 55.89           H
ATOM     27 HG11 VAL A  10      40.287  30.061 102.744  0.00 55.89           H
ATOM     28 HG12 VAL A  10      39.913  30.033 104.489  0.00 55.89           H
ATOM     29 HG13 VAL A  10      40.026  28.503 103.575  0.00 55.89           H
ATOM     30 HG21 VAL A  10      38.585  29.528 101.111  0.00 55.89           H
ATOM     31 HG22 VAL A  10      38.286  27.954 101.898  0.00 55.89           H
ATOM     32 HG23 VAL A  10      36.934  29.089 101.632  0.00 55.89           H
ATOM     33  N   GLU A  11      35.257  29.874 104.732  1.00 38.26           N
ATOM     34  CA  GLU A  11      33.832  30.104 104.541  1.00 38.26           C
ATOM     35  C   GLU A  11      33.661  31.534 104.025  1.00 38.26           C
ATOM     36  O   GLU A  11      34.303  32.461 104.515  1.00 38.26           O
ATOM     37  CB  GLU A  11      33.113  29.898 105.880  1.00 47.28           C
ATOM     38  CG  GLU A  11      31.647  30.177 105.867  1.00 47.28           C
ATOM     39  CD  GLU A  11      30.927  29.651 107.101  1.00 47.28           C
ATOM     40  OE1 GLU A  11      31.427  29.829 108.233  1.00 47.28           O
ATOM     41  OE2 GLU A  11      29.834  29.066 106.935  1.00 47.28           O
ATOM     42  OXT GLU A  11      32.873  31.768 103.108  1.00 38.26           O
ATOM     43  H   GLU A  11      35.663  30.428 105.486  0.00 38.26           H
ATOM     44  HA  GLU A  11      33.384  29.414 103.826  0.00 38.26           H
ATOM     45  HB2 GLU A  11      33.244  28.859 106.184  0.00 47.28           H
ATOM     46  HB3 GLU A  11      33.564  30.562 106.618  0.00 47.28           H
ATOM     47  HG2 GLU A  11      31.493  31.255 105.823  0.00 47.28           H
ATOM     48  HG3 GLU A  11      31.204  29.701 104.992  0.00 47.28           H
'''

def main():
  f=open('pdb_1x24_modified.pdb', 'w')
  f.write(pdb_1x24_modified)
  del f
  rc = run_program(electrons.Program, args=['pdb_1x24_modified.pdb'])
  print('-'*80)
  print(rc.keys())
  print(rc)

  f=open('pdb_1x24_modified_complete.pdb', 'w')
  f.write(pdb_1x24_modified_complete)
  del f
  rc = run_program(electrons.Program, args=['pdb_1x24_modified_complete.pdb'])
  print('-'*80)
  print(rc.keys())
  print(rc)

if __name__ == '__main__':
  main()
