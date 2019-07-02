from __future__ import absolute_import, division, print_function
from mmtbx.tls import tools
import math
import time

pdb_str_1 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
ATOM      1  CA  THR A   6       0.000   0.000   0.000  1.00  0.00           C
ATOM      1  CA  THR B   6       3.000   0.000   0.000  1.00  0.00           C
"""

pdb_str_2 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
ATOM      1  CA  THR A   6       0.000   0.000   0.000  1.00  0.00           C
ATOM      1  CA  THR B   6       0.000   3.000   0.000  1.00  0.00           C
"""

pdb_str_3 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
ATOM      1  CA  THR A   6       0.000   0.000   0.000  1.00  0.00           C
ATOM      1  CA  THR B   6       0.000   0.000   3.000  1.00  0.00           C
"""

pdb_str_4 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
ATOM      1  CA  THR A   6       0.000   0.000   0.000  1.00  0.00           C
ATOM      1  CA  THR B   6       1.000   2.000   3.000  1.00  0.00           C
"""

def exercise_03():
  sqrt = math.sqrt
  vs = []
  vs.append( [(sqrt(2)/2, sqrt(2)/2, 0), (-sqrt(2)/2, sqrt(2)/2, 0), (0,0,1)] )
  vs.append( [(1,0,0), (0, sqrt(2)/2, sqrt(2)/2), (0, -sqrt(2)/2, sqrt(2)/2)] )
  vs.append( [(sqrt(3)/2, 1/2, 0), (-1/2, sqrt(3)/2, 0), (0,0,1)] )
  vs.append( [(1,0,0), (0, sqrt(3)/2, 1/2), (0, -1/2, sqrt(3)/2)] )
  for pdb_str in [pdb_str_1, pdb_str_2, pdb_str_3, pdb_str_4]:
    for vs_ in vs:
      vx,vy,vz = vs_
      print(vx,vy,vz)
      tools.u_tls_vs_u_ens(pdb_str=pdb_str,
          tx=0.05,ty=0.07,tz=0.09,
          vx=vx, vy=vy, vz=vz,
          n_models=1000)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_03()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
