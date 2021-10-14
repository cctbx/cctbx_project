from __future__ import absolute_import, division, print_function
import os

pdb_lines = ['''
CRYST1   45.938   57.554   63.431  90.00  90.00  90.00 P 21 21 21    4
ATOM    865  N   THR A  72       9.071 -34.112 -15.835  1.00 14.66           N
ANISOU  865  N   THR A  72     2090   1604   1878   -469   -585    286       N
ATOM    866  CA  THR A  72       7.946 -34.093 -14.896  1.00 14.79           C
ANISOU  866  CA  THR A  72     2217   1589   1815   -685   -469    195       C
ATOM    867  C   THR A  72       8.327 -33.683 -13.479  1.00 15.00           C
ANISOU  867  C   THR A  72     2303   1582   1815   -693   -182     27       C
ATOM    868  O   THR A  72       7.465 -33.583 -12.605  1.00 19.09           O
ANISOU  868  O   THR A  72     2603   2288   2364   -983   -106   -162       O
ATOM    869  CB  THR A  72       7.251 -35.461 -14.825  1.00 17.03           C
ANISOU  869  CB  THR A  72     2487   1931   2053   -651   -443     79       C
ATOM    870  OG1 THR A  72       8.205 -36.457 -14.441  1.00 18.25           O
ANISOU  870  OG1 THR A  72     2882   1816   2235   -133    -90    157       O
ATOM    871  CG2 THR A  72       6.662 -35.825 -16.174  1.00 18.10           C
ANISOU  871  CG2 THR A  72     2424   2248   2204   -906   -567     37       C
ATOM    872  H  ATHR A  72       9.611 -34.774 -15.731  0.73 49.72           H
ATOM    873  H  BTHR A  72       9.561 -34.817 -15.782  0.27 14.64           H
ATOM    874  N   GLY A  73       9.615 -33.455 -13.258  1.00 11.31           N
ANISOU  874  N   GLY A  73     1642   1142   1514   -333   -417    118       N
ATOM    875  CA  GLY A  73      10.114 -33.087 -11.948  1.00 10.77           C
ANISOU  875  CA  GLY A  73     1613   1063   1416   -193   -278     31       C
ATOM    876  C   GLY A  73      11.312 -33.932 -11.572  1.00 11.19           C
ANISOU  876  C   GLY A  73     2044    715   1494     68   -328   -139       C
ATOM    877  O   GLY A  73      11.706 -34.822 -12.321  1.00 12.19           O
ANISOU  877  O   GLY A  73     2296    807   1529    149   -161   -267       O
ATOM    878  HA2 GLY A  73      10.382 -32.155 -11.960  0.78  6.73           H
ATOM    879  HA3 GLY A  73       9.429 -33.206 -11.273  1.00 12.28           H
ATOM    880  N   PRO A  74      11.893 -33.676 -10.395  1.00 10.97           N
ANISOU  880  N   PRO A  74     1928    761   1478    243   -483   -116       N
ATOM    881  CA  PRO A  74      11.404 -32.709  -9.413  1.00  9.71           C
ANISOU  881  CA  PRO A  74     1557    872   1261     79   -246    -42       C
ATOM    882  C   PRO A  74      11.673 -31.272  -9.813  1.00  8.45           C
ANISOU  882  C   PRO A  74     1278    767   1166     10   -162   -146       C
ATOM    883  O   PRO A  74      12.628 -30.968 -10.533  1.00  8.97           O
ANISOU  883  O   PRO A  74     1245    764   1398     96   -168   -189       O
ATOM    884  CB  PRO A  74      12.210 -33.057  -8.160  1.00 11.48           C
ANISOU  884  CB  PRO A  74     1644   1305   1413    129   -290    -11       C
ATOM    885  CG  PRO A  74      13.489 -33.597  -8.697  1.00 12.84           C
ANISOU  885  CG  PRO A  74     1823   1499   1557    270   -573   -138       C
ATOM    886  CD  PRO A  74      13.104 -34.371  -9.927  1.00 12.36           C
ANISOU  886  CD  PRO A  74     1792   1312   1592    310   -500   -269       C
ATOM    887  HA  PRO A  74      10.458 -32.837  -9.242  0.77 10.11           H
ATOM    888  HD2 PRO A  74      13.805 -34.272 -10.591  1.00 13.53           H
ATOM    889  N   PHE A  75      10.815 -30.390  -9.330  1.00  7.53           N
ANISOU  889  N   PHE A  75     1141    749    970    -68   -209   -160       N
ATOM    890  CA  PHE A  75      11.088 -28.971  -9.410  1.00  6.73           C
ANISOU  890  CA  PHE A  75      969    707    882     13   -225    -91       C
ATOM    891  C   PHE A  75      11.768 -28.530  -8.116  1.00  6.60           C
ANISOU  891  C   PHE A  75      971    680    858     -2   -212    -68       C
ATOM    892  O   PHE A  75      12.132 -29.369  -7.289  1.00  8.29           O
ANISOU  892  O   PHE A  75     1418    744    990    -25   -440     12       O
ATOM    893  CB  PHE A  75       9.804 -28.226  -9.778  1.00  8.14           C
ANISOU  893  CB  PHE A  75     1048    911   1134     21   -317   -127       C
ATOM    894  CG  PHE A  75       9.300 -28.635 -11.125  1.00  8.35           C
ANISOU  894  CG  PHE A  75     1213    753   1208     81   -564   -118       C
ATOM    895  CD1 PHE A  75       9.930 -28.159 -12.259  1.00  8.58           C
ANISOU  895  CD1 PHE A  75     1336    733   1191     68   -605    -52       C
ATOM    896  CD2 PHE A  75       8.286 -29.573 -11.261  1.00  9.92           C
ANISOU  896  CD2 PHE A  75     1416   1020   1335    -21   -673    -27       C
ATOM    897  CE1 PHE A  75       9.518 -28.569 -13.511  1.00  9.43           C
ANISOU  897  CE1 PHE A  75     1520    854   1208    103   -714    -17       C
ATOM    898  CE2 PHE A  75       7.870 -29.983 -12.509  1.00 10.44           C
ANISOU  898  CE2 PHE A  75     1544   1039   1384     29   -668   -155       C
ATOM    899  CZ  PHE A  75       8.488 -29.487 -13.634  1.00 10.58           C
ANISOU  899  CZ  PHE A  75     1631   1050   1339    141   -831    -86       C
ATOM    900  HB3 PHE A  75       9.989 -27.275  -9.808  1.00 13.04           H
ATOM    901  HD1 PHE A  75      10.624 -27.547 -12.177  1.00  8.63           H
ATOM    902  HE2 PHE A  75       7.176 -30.598 -12.589  0.67 10.49           H
ATOM    903  HZ  PHE A  75       8.202 -29.755 -14.477  1.00 25.56           H
ATOM    904  N   GLY A  76      12.000 -27.235  -7.966  1.00  6.65           N
ANISOU  904  N   GLY A  76      900    774    854     -7   -139   -133       N
ATOM    905  CA  GLY A  76      12.650 -26.734  -6.775  1.00  5.67           C
ANISOU  905  CA  GLY A  76      746    701    709      3   -133   -107       C
ATOM    906  C   GLY A  76      11.672 -26.515  -5.645  1.00  5.32           C
ANISOU  906  C   GLY A  76      704    619    697    -49   -170    -30       C
ATOM    907  O   GLY A  76      10.542 -27.013  -5.659  1.00  6.39           O
ANISOU  907  O   GLY A  76      758    892    780   -237   -106   -122       O
ATOM    908  H   GLY A  76      11.791 -26.629  -8.538  1.00 10.05           H
ATOM    909  HA3 GLY A  76      13.113 -25.903  -6.960  1.00  6.44           H
''',
'''
ATOM    729  N   PHE A 150     241.493 239.740 110.966  1.00  0.73           N
ATOM    730  CA  PHE A 150     241.495 240.006 112.396  1.00  0.73           C
ATOM    731  C   PHE A 150     241.654 238.693 113.139  1.00  0.73           C
ATOM    732  O   PHE A 150     241.102 237.680 112.690  1.00  0.73           O
ATOM    733  CB  PHE A 150     240.203 240.701 112.838  1.00  0.73           C
ATOM    734  N   PRO A 151     242.398 238.656 114.253  1.00  0.76           N
ATOM    735  CA  PRO A 151     243.254 239.710 114.805  1.00  0.76           C
ATOM    736  C   PRO A 151     244.678 239.604 114.277  1.00  0.76           C
ATOM    737  O   PRO A 151     245.073 238.585 113.720  1.00  0.76           O
ATOM    738  CB  PRO A 151     243.170 239.452 116.311  1.00  0.76           C
ATOM    739  N   GLU A 152     245.470 240.656 114.436  1.00  0.77           N
ATOM    740  CA  GLU A 152     246.888 240.568 114.161  1.00  0.77           C
ATOM    741  C   GLU A 152     247.540 239.683 115.222  1.00  0.77           C
ATOM    742  O   GLU A 152     246.934 239.400 116.257  1.00  0.77           O
ATOM    743  CB  GLU A 152     247.496 241.970 114.140  1.00  0.77           C
''',
]

from libtbx import easy_run

def main():
  tf = 'test_cis_127.pdb'
  with open(tf, 'w') as f:
    f.write(pdb_lines[0])
  cmd = 'phenix.pdb_interpretation write_geo=True cis_pro_eh99=False %s' % tf
  print(cmd)
  easy_run.go(cmd)
  with open('%s.geo' % tf, 'r') as f:
    lines = f.read()
  find = '''angle pdb=" C   GLY A  73 "
      pdb=" N   PRO A  74 "
      pdb=" CA  PRO A  74 "
    ideal   model   delta    sigma   weight residual
   122.60  124.91   -2.31 5.00e+00 4.00e-02 2.14e-01'''
  assert lines.find(find)>-1
  cmd = 'phenix.pdb_interpretation write_geo=True cis_pro_eh99=True %s' % tf
  print(cmd)
  easy_run.go(cmd)
  with open('%s.geo' % tf, 'r') as f:
    lines = f.read()
  find = '''angle pdb=" C   GLY A  73 "
      pdb=" N   PRO A  74 "
      pdb=" CA  PRO A  74 "
    ideal   model   delta    sigma   weight residual
   127.00  124.91    2.09 2.40e+00 1.74e-01 7.57e-01'''
  assert lines.find(find)>-1

  tf = 'test_cis_127_1.pdb'
  with open(tf, 'w') as f:
    f.write(pdb_lines[1])
  cmd = 'phenix.pdb_interpretation write_geo=True cis_pro_eh99=True %s' % tf
  print(cmd)
  easy_run.go(cmd)
  assert os.path.exists('%s.geo' % tf), 'command failed due to missing CD'


if __name__ == '__main__':
  main()
