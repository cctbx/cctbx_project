
from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from io import StringIO

def exercise () :
  open("tmp_interpolation_start.pdb", "w").write("""\
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
""")
  open("tmp_interpolation_end.pdb", "w").write("""\
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7      10.002   1.319   9.510  1.00 14.45           C
ATOM     53  CD1 TYR A   7       9.703   2.284  10.475  1.00 15.68           C
ATOM     54  CD2 TYR A   7      11.222   0.661   9.584  1.00 14.80           C
ATOM     55  CE1 TYR A   7      10.573   2.540  11.514  1.00 13.46           C
ATOM     56  CE2 TYR A   7      12.116   0.900  10.635  1.00 14.33           C
ATOM     57  CZ  TYR A   7      11.782   1.856  11.583  1.00 15.09           C
ATOM     58  OH  TYR A   7      12.625   2.114  12.630  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
""")
  from mmtbx.command_line import interpolate
  args = [
    "tmp_interpolation_start.pdb",
    "tmp_interpolation_end.pdb",
    "interpolate_dihedrals=True",
    "align_atoms=None",
  ]
  out = StringIO()
  result = interpolate.run(args=args, out=out)
  print("OK")

if (__name__ == "__main__") :
  exercise()
