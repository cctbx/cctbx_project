"""
Andrey A. Lebedev, Alexei A. Vagin & Garib N. Murshudov
Acta Cryst. (2006). D62, 83-95.
http://journals.iucr.org/d/issues/2006/01/00/ba5089/index.html
Appendix A1. Algorithms used in the determination of twinning operators
and their type of merohedry

The formula for the perturbation mentioned in the appendix and the
matrices_from_email below kindly provided by Andrey.

This script reproduces the perturbations ("scores") given by Andrey.
It also shows the Le Page (1982, J. Appl. Cryst. 15, 255-259) deltas
in radians for comparison.
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx
from cctbx import uctbx
from scitbx import matrix
from six.moves import range

matrices_from_email = """
CRYST1   82.053   66.612   84.904  99.00 111.29 105.00 P 1

S:
   6732.6948  -1414.6310  -2529.5033
  -1414.6310   4437.1585   -884.7347
  -2529.5033   -884.7347   7208.6892

M:
           0           0          -1
           0          -1           0
          -1           0           0

Score:
  0.06996537

------------------------------------------
CRYST1   85.053   36.612   22.904  91.00 101.29 110.00 P 1

S:
   6176.9399    260.7664    128.5781
    260.7664   1340.4385    -14.6349
    128.5781    -14.6349    524.5932

M:
          -1           0           0
           0          -1           0
           1           0           1

Score:
  0.06177181

------------------------------------------
CRYST1   87.053   86.612   84.904  92.00 101.29 105.00 P 1

S:
   7578.2248  -1951.4527  -1447.0019
  -1951.4527   7501.6385   -256.6406
  -1447.0019   -256.6406   7208.6892

M:
          -1           0           0
           0           0          -1
           0          -1           0

Score:
  0.03867617

------------------------------------------
CRYST1   89.053   96.612   84.904  93.00 101.29 110.00 P 1

S:
   7930.4368  -2942.6005  -1480.2461
  -2942.6005   9333.8785   -429.2985
  -1480.2461   -429.2985   7208.6892

M:
           0          -1           0
          -1           0           0
           0           0          -1

Score:
  0.11207609

------------------------------------------
CRYST1   92.053   66.612   84.904  94.00  91.29 105.00 P 1

S:
   8473.7548  -1587.0355   -175.9529
  -1587.0355   4437.1585   -394.5165
   -175.9529   -394.5165   7208.6892

M:
          -1           0           0
           0          -1           0
           0           0           1

Score:
  0.06714734

------------------------------------------
CRYST1  102.053   46.612   74.904  95.00  71.29  95.00 P 1

S:
  10414.8148   -414.5907   2452.0864
   -414.5907   2172.6785   -304.2978
   2452.0864   -304.2978   5610.6092

M:
           1           0           0
           0          -1           0
          -1           0          -1

Score:
  0.06543032

------------------------------------------
CRYST1  102.053   66.612   64.904  96.00 111.29  90.00 P 1

S:
   9817.4017   -451.9168   1807.5581
   -451.9168   4437.1585   -451.9168
   1807.5581   -451.9168   4212.5292

M:
           1           0           0
           0          -1           0
          -1           0          -1

Score:
  0.05202084

------------------------------------------
CRYST1   53.053   96.612   54.904  97.00  91.29  85.00 P 1

S:
   2814.6208    446.7217    -65.5759
    446.7217   9333.8785   -646.4419
    -65.5759   -646.4419   3014.4492

M:
           0           0           1
           0          -1           0
           1           0           0

Score:
  0.03361885

------------------------------------------
CRYST1  102.053   36.612   44.904  98.00 111.29  80.00 P 1

S:
   9103.4130    420.0088    352.4837
    420.0088   1340.4385   -228.8041
    352.4837   -228.8041   2016.3692

M:
           1           0           0
          -1          -1           0
           0           0          -1

Score:
  0.10323493

------------------------------------------
CRYST1  102.053   66.612   34.904 100.00 111.29  75.00 P 1

S:
   9046.4187   1355.7037    -75.0535
   1355.7037   4437.1585   -403.7364
    -75.0535   -403.7364   1218.2892

M:
          -1           0           0
           0          -1           0
           0          -1           1

Score:
  0.08193190
"""

def run():
  lines = iter(matrices_from_email.splitlines())
  while True:
    for line in lines:
      if (line.rstrip() == "S:"): break
    else:
      break
    s = []
    for i in range(3):
      s.extend([float(v) for v in next(lines).split()])
    for line in lines:
      if (line.rstrip() == "M:"): break
    else:
      raise RuntimeError("S: found but not M:")
    m = []
    for i in range(3):
      m.extend([int(v) for v in next(lines).split()])
    for line in lines:
      if (line.rstrip() == "Score:"): break
    else:
      raise RuntimeError("S: and M: found but not Score:")
    score = float(next(lines).strip())
    s = matrix.sqr(s)
    m = matrix.sqr(m)
    print(s.mathematica_form(label="s", one_row_per_line=True))
    print(m.mathematica_form(label="m", one_row_per_line=True))
    r = sgtbx.rot_mx(m, 1)
    print("rotation type:", r.info().type())
    print("axis direction:", r.info().ev())
    u = uctbx.unit_cell(metrical_matrix=s.as_sym_mat3())
    p = r.lebedev_2005_perturbation(reduced_cell=u)
    print("score given:     ", score)
    print("score reproduced:", p)
    assert abs(p-score) < 1.e-6
    delta = r.le_page_1982_delta(reduced_cell=u)
    print("Le Page delta:   ", delta)
    print()
  print("OK")

if (__name__ == "__main__"):
  run()
