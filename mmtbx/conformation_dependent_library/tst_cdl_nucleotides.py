from __future__ import division
from libtbx import easy_run

pdb = '''
ATOM      1  P     C A   1       0.226 -10.894   0.425  1.00 76.88           P
ATOM      2  OP1   C A   1       0.482 -11.761   1.602  1.00 78.59           O
ATOM      3  OP2   C A   1      -0.746  -9.788   0.472  1.00 75.20           O
ATOM      4  O5'   C A   1       1.619 -10.251  -0.015  1.00 78.76           O
ATOM      5  C5'   C A   1       2.889 -10.871   0.240  1.00 77.93           C
ATOM      6  C4'   C A   1       4.016  -9.909  -0.104  1.00 78.66           C
ATOM      7  O4'   C A   1       3.904  -9.426  -1.486  1.00 79.49           O
ATOM      8  C3'   C A   1       4.077  -8.641   0.739  1.00 78.23           C
ATOM      9  O3'   C A   1       4.714  -8.922   1.993  1.00 77.17           O
ATOM     10  C2'   C A   1       4.872  -7.722  -0.207  1.00 79.34           C
ATOM     11  O2'   C A   1       6.269  -7.992  -0.276  1.00 79.34           O
ATOM     12  C1'   C A   1       4.252  -8.053  -1.569  1.00 78.22           C
ATOM     13  N1    C A   1       3.007  -7.180  -1.992  1.00 20.00           N
ATOM     14  C2    C A   1       3.150  -5.813  -2.248  1.00 20.00           C
ATOM     15  N3    C A   1       2.064  -5.098  -2.618  1.00 20.00           N
ATOM     16  C4    C A   1       0.877  -5.697  -2.729  1.00 20.00           C
ATOM     17  C5    C A   1       0.707  -7.085  -2.458  1.00 20.00           C
ATOM     18  C6    C A   1       1.784  -7.778  -2.098  1.00 20.00           C
ATOM     19  O2    C A   1       4.275  -5.294  -2.125  1.00 20.00           O
ATOM     20  N4    C A   1      -0.163  -4.952  -3.107  1.00 20.00           N
ATOM     21  P     G A   2       4.132  -8.090   3.235  1.00 76.88           P
ATOM     22  OP1   G A   2       4.816  -8.681   4.412  1.00 78.59           O
ATOM     23  OP2   G A   2       2.717  -7.684   3.282  1.00 75.20           O
ATOM     24  O5'   G A   2       4.957  -6.796   2.795  1.00 78.76           O
ATOM     25  C5'   G A   2       6.360  -6.632   3.050  1.00 77.93           C
ATOM     26  C4'   G A   2       6.789  -5.214   2.706  1.00 78.66           C
ATOM     27  O4'   G A   2       6.434  -4.868   1.324  1.00 79.49           O
ATOM     28  C3'   G A   2       6.156  -4.114   3.549  1.00 78.23           C
ATOM     29  O3'   G A   2       6.843  -4.006   4.803  1.00 77.17           O
ATOM     30  C2'   G A   2       6.328  -2.911   2.603  1.00 79.34           C
ATOM     31  O2'   G A   2       7.650  -2.383   2.534  1.00 79.34           O
ATOM     32  C1'   G A   2       5.985  -3.524   1.241  1.00 78.22           C
ATOM     33  N9    G A   2       4.566  -3.455   0.848  1.00 20.00           N
ATOM     34  C8    G A   2       3.693  -4.506   0.681  1.00 20.00           C
ATOM     35  N7    G A   2       2.484  -4.139   0.318  1.00 20.00           N
ATOM     36  C5    G A   2       2.545  -2.754   0.231  1.00 20.00           C
ATOM     37  C4    G A   2       3.812  -2.311   0.556  1.00 20.00           C
ATOM     38  N1    G A   2       2.016  -0.497  -0.062  1.00 20.00           N
ATOM     39  C2    G A   2       3.298  -0.143   0.266  1.00 20.00           C
ATOM     40  N3    G A   2       4.242  -1.013   0.586  1.00 20.00           N
ATOM     41  C6    G A   2       1.539  -1.808  -0.102  1.00 20.00           C
ATOM     42  O6    G A   2       0.354  -2.021  -0.408  1.00 20.00           O
ATOM     43  N2    G A   2       3.585   1.155   0.252  1.00 20.00           N
ATOM     44  P     A A   3       5.904  -3.620   6.045  1.00 76.88           P
ATOM     45  OP1   A A   3       6.799  -3.748   7.222  1.00 78.59           O
ATOM     46  OP2   A A   3       4.494  -4.043   6.092  1.00 75.20           O
ATOM     47  O5'   A A   3       5.899  -2.086   5.605  1.00 78.76           O
ATOM     48  C5'   A A   3       6.991  -1.189   5.860  1.00 77.93           C
ATOM     49  C4'   A A   3       6.586   0.236   5.516  1.00 78.66           C
ATOM     50  O4'   A A   3       6.100   0.335   4.134  1.00 79.49           O
ATOM     51  C3'   A A   3       5.459   0.819   6.359  1.00 78.23           C
ATOM     52  O3'   A A   3       5.979   1.282   7.613  1.00 77.17           O
ATOM     53  C2'   A A   3       4.954   1.925   5.413  1.00 79.34           C
ATOM     54  O2'   A A   3       5.781   3.083   5.344  1.00 79.34           O
ATOM     55  C1'   A A   3       4.997   1.223   4.051  1.00 78.22           C
ATOM     56  N9    A A   3       3.745   0.511   3.653  1.00 20.00           N
ATOM     57  C8    A A   3       3.592  -0.836   3.491  1.00 20.00           C
ATOM     58  N7    A A   3       2.386  -1.187   3.127  1.00 20.00           N
ATOM     59  C5    A A   3       1.698   0.009   3.047  1.00 20.00           C
ATOM     60  C4    A A   3       2.524   1.076   3.368  1.00 20.00           C
ATOM     61  N1    A A   3      -0.009   1.601   2.733  1.00 20.00           N
ATOM     62  C2    A A   3       0.899   2.529   3.062  1.00 20.00           C
ATOM     63  N3    A A   3       2.181   2.373   3.392  1.00 20.00           N
ATOM     64  C6    A A   3       0.372   0.307   2.716  1.00 20.00           C
ATOM     65  N6    A A   3      -0.525  -0.629   2.388  1.00 20.00           N
ATOM     66  P     A B   1     -10.978  -1.171   1.417  1.00 76.88           P
ATOM     67  OP1   A B   1     -11.872  -1.039   0.240  1.00 78.59           O
ATOM     68  OP2   A B   1      -9.746  -1.977   1.370  1.00 75.20           O
ATOM     69  O5'   A B   1     -10.538   0.299   1.857  1.00 78.76           O
ATOM     70  C5'   A B   1     -11.330   1.469   1.602  1.00 77.93           C
ATOM     71  C4'   A B   1     -10.536   2.720   1.946  1.00 78.66           C
ATOM     72  O4'   A B   1     -10.043   2.677   3.328  1.00 79.49           O
ATOM     73  C3'   A B   1      -9.290   2.959   1.103  1.00 78.23           C
ATOM     74  O3'   A B   1      -9.658   3.550  -0.151  1.00 77.17           O
ATOM     75  C2'   A B   1      -8.492   3.875   2.049  1.00 79.34           C
ATOM     76  O2'   A B   1      -8.956   5.220   2.118  1.00 79.34           O
ATOM     77  C1'   A B   1      -8.732   3.215   3.411  1.00 78.22           C
ATOM     78  N9    A B   1      -7.735   2.176   3.810  1.00 20.00           N
ATOM     79  C8    A B   1      -7.971   0.841   3.972  1.00 20.00           C
ATOM     80  N7    A B   1      -6.914   0.162   4.336  1.00 20.00           N
ATOM     81  C5    A B   1      -5.915   1.113   4.416  1.00 20.00           C
ATOM     82  C4    A B   1      -6.403   2.371   4.095  1.00 20.00           C
ATOM     83  N1    A B   1      -3.825   2.154   4.730  1.00 20.00           N
ATOM     84  C2    A B   1      -4.432   3.302   4.401  1.00 20.00           C
ATOM     85  N3    A B   1      -5.706   3.517   4.071  1.00 20.00           N
ATOM     86  C6    A B   1      -4.559   1.022   4.747  1.00 20.00           C
ATOM     87  N6    A B   1      -3.964  -0.130   5.075  1.00 20.00           N
ATOM     88  P     U B   2      -8.752   3.091  -1.393  1.00 76.88           P
ATOM     89  OP1   U B   2      -9.433   3.685  -2.570  1.00 78.59           O
ATOM     90  OP2   U B   2      -8.151   1.747  -1.440  1.00 75.20           O
ATOM     91  O5'   U B   2      -7.587   4.090  -0.953  1.00 78.76           O
ATOM     92  C5'   U B   2      -7.622   5.503  -1.208  1.00 77.93           C
ATOM     93  C4'   U B   2      -6.279   6.127  -0.864  1.00 78.66           C
ATOM     94  O4'   U B   2      -5.886   5.824   0.518  1.00 79.49           O
ATOM     95  C3'   U B   2      -5.100   5.654  -1.707  1.00 78.23           C
ATOM     96  O3'   U B   2      -5.091   6.351  -2.961  1.00 77.17           O
ATOM     97  C2'   U B   2      -3.934   5.995  -0.761  1.00 79.34           C
ATOM     98  O2'   U B   2      -3.598   7.377  -0.692  1.00 79.34           O
ATOM     99  C1'   U B   2      -4.493   5.569   0.601  1.00 78.22           C
ATOM    100  N1    U B   2      -4.200   4.068   1.024  1.00 20.00           N
ATOM    101  C2    U B   2      -2.904   3.668   1.275  1.00 20.00           C
ATOM    102  N3    U B   2      -2.757   2.346   1.647  1.00 20.00           N
ATOM    103  C4    U B   2      -3.761   1.415   1.783  1.00 20.00           C
ATOM    104  C5    U B   2      -5.081   1.922   1.496  1.00 20.00           C
ATOM    105  C6    U B   2      -5.264   3.207   1.136  1.00 20.00           C
ATOM    106  O2    U B   2      -1.952   4.421   1.179  1.00 20.00           O
ATOM    107  O4    U B   2      -3.505   0.253   2.122  1.00 20.00           O
ATOM    108  P     C B   3      -4.576   5.475  -4.203  1.00 76.88           P
ATOM    109  OP1   C B   3      -4.829   6.343  -5.380  1.00 78.59           O
ATOM    110  OP2   C B   3      -4.796   4.019  -4.250  1.00 75.20           O
ATOM    111  O5'   C B   3      -3.056   5.687  -3.763  1.00 78.76           O
ATOM    112  C5'   C B   3      -2.323   6.894  -4.018  1.00 77.93           C
ATOM    113  C4'   C B   3      -0.855   6.694  -3.674  1.00 78.66           C
ATOM    114  O4'   C B   3      -0.688   6.227  -2.292  1.00 79.49           O
ATOM    115  C3'   C B   3      -0.118   5.660  -4.517  1.00 78.23           C
ATOM    116  O3'   C B   3       0.266   6.240  -5.771  1.00 77.17           O
ATOM    117  C2'   C B   3       1.047   5.316  -3.571  1.00 79.34           C
ATOM    118  O2'   C B   3       2.077   6.298  -3.502  1.00 79.34           O
ATOM    119  C1'   C B   3       0.347   5.259  -2.209  1.00 78.22           C
ATOM    120  N1    C B   3      -0.229   3.852  -1.785  1.00 20.00           N
ATOM    121  C2    C B   3       0.630   2.778  -1.530  1.00 20.00           C
ATOM    122  N3    C B   3       0.102   1.591  -1.159  1.00 20.00           N
ATOM    123  C4    C B   3      -1.221   1.453  -1.048  1.00 20.00           C
ATOM    124  C5    C B   3      -2.113   2.530  -1.319  1.00 20.00           C
ATOM    125  C6    C B   3      -1.581   3.694  -1.679  1.00 20.00           C
ATOM    126  O2    C B   3       1.856   2.950  -1.652  1.00 20.00           O
ATOM    127  N4    C B   3      -1.693   0.265  -0.670  1.00 20.00           N
ATOM    128  P     G B   4       0.226   5.226  -7.013  1.00 76.88           P
ATOM    129  OP1   G B   4       0.482   6.093  -8.190  1.00 78.59           O
ATOM    130  OP2   G B   4      -0.746   4.120  -7.060  1.00 75.20           O
ATOM    131  O5'   G B   4       1.619   4.583  -6.573  1.00 78.76           O
ATOM    132  C5'   G B   4       2.889   5.203  -6.828  1.00 77.93           C
ATOM    133  C4'   G B   4       4.016   4.241  -6.484  1.00 78.66           C
ATOM    134  O4'   G B   4       3.904   3.758  -5.102  1.00 79.49           O
ATOM    135  C3'   G B   4       4.077   2.973  -7.327  1.00 78.23           C
ATOM    136  O3'   G B   4       4.714   3.254  -8.581  1.00 77.17           O
ATOM    137  C2'   G B   4       4.872   2.054  -6.381  1.00 79.34           C
ATOM    138  O2'   G B   4       6.269   2.324  -6.312  1.00 79.34           O
ATOM    139  C1'   G B   4       4.252   2.385  -5.019  1.00 78.22           C
ATOM    140  N9    G B   4       3.095   1.560  -4.626  1.00 20.00           N
ATOM    141  C8    G B   4       1.793   1.973  -4.459  1.00 20.00           C
ATOM    142  N7    G B   4       0.973   1.011  -4.096  1.00 20.00           N
ATOM    143  C5    G B   4       1.773  -0.122  -4.009  1.00 20.00           C
ATOM    144  C4    G B   4       3.079   0.190  -4.334  1.00 20.00           C
ATOM    145  N1    G B   4       2.547  -2.307  -3.716  1.00 20.00           N
ATOM    146  C2    G B   4       3.817  -1.912  -4.044  1.00 20.00           C
ATOM    147  N3    G B   4       4.142  -0.670  -4.364  1.00 20.00           N
ATOM    148  C6    G B   4       1.437  -1.461  -3.676  1.00 20.00           C
ATOM    149  O6    G B   4       0.325  -1.922  -3.370  1.00 20.00           O
ATOM    150  N2    G B   4       4.760  -2.849  -4.030  1.00 20.00           N
'''

std_asserts = [
  [
  'enable = 0',
  'esd = *phenix csd',
  ],
  [
  'enable = 0',
  'esd = phenix *csd',
  ],
  [
  'enable = 1',
  'esd = *phenix csd',
  ],
  [
  'enable = 1',
  'esd = phenix *csd',
  ],
]
geo_asserts = []
geo_asserts.append(
  ['''bond pdb=" P     G A   2 "
     pdb=" OP1   G A   2 "
  ideal  model  delta    sigma   weight residual
  1.485  1.484  0.001 2.00e-02 2.50e+03 2.16e-03''',
  '''angle pdb=" O3'   C B   3 "
      pdb=" P     G B   4 "
      pdb=" O5'   G B   4 "
    ideal   model   delta    sigma   weight residual
   104.00   91.11   12.89 1.50e+00 4.44e-01 7.38e+01''',
  ])
geo_asserts.append(geo_asserts[0])
geo_asserts.append(
  ['''bond pdb=" P     G A   2 "
     pdb=" OP1   G A   2 "
  ideal  model  delta    sigma   weight residual
  1.484  1.484 -0.000 2.00e-02 2.50e+03 1.25e-05''',
  '''angle pdb=" O3'   C B   3 "
      pdb=" P     G B   4 "
      pdb=" O5'   G B   4 "
    ideal   model   delta    sigma   weight residual
   104.20   91.11   13.09 1.50e+00 4.44e-01 7.61e+01''',
  ])
geo_asserts.append(
  ['''bond pdb=" P     G A   2 "
     pdb=" OP1   G A   2 "
  ideal  model  delta    sigma   weight residual
  1.484  1.484 -0.000 1.20e-02 6.94e+03 3.48e-05''',
  '''angle pdb=" O3'   C B   3 "
      pdb=" P     G B   4 "
      pdb=" O5'   G B   4 "
    ideal   model   delta    sigma   weight residual
   104.20   91.11   13.09 1.50e+00 4.44e-01 7.61e+01''',
  ])

restraintlib_installed=True
try:
  from restraintlib import launcher  # special import
  from restraintlib.printer import TuplePrinter  # special import
  from restraintlib.restraints import analyze_pdb_hierarhy  # special import
except ImportError as e:
  restraintlib_installed = False

def main():
  i=0
  filename = 'cdl_nucleotide.pdb'
  with open(filename, 'w') as f:
    f.write(pdb)
  for param1 in range(2):
    if param1 and not restraintlib_installed: break
    for param2 in ['phenix', 'csd']:
      cmd = 'phenix.pdb_interpretation %s write_geo=True cdl_nucleotides.enable=%s cdl_nucleotides.esd=%s' % (
        filename,
        param1,
        param2,
        )
      cmd += ' cdl_nucleotides.factor=1.'
      print(i, param1, param2, cmd)
      rc = easy_run.go(cmd)
      lines=[]
      for line in rc.stdout_lines:
        lines.append(line.strip())
      for line in std_asserts[i]:
        assert line in lines, 'missing %s' % line

      with open('%s.geo' % filename, 'r') as f:
        lines=f.read()
      for line in geo_asserts[i]:
        assert lines.find(line)>-1, 'missing %s' % line

      i+=1
    print('OK')
  print('OK')

if __name__ == '__main__':
  main()
