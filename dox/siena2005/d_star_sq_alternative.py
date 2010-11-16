from cctbx import uctbx

def d_star_sq(unit_cell, miller_index):
  """\
h = {h0,h1,h2}
g = {{aa,ab,ac},{ab,bb,bc},{ac,bc,cc}}
dss = FullSimplify[h.Inverse[g].h]
FortranForm[dss]
"""
  aa,bb,cc,ab,ac,bc = unit_cell.metrical_matrix()
  h0,h1,h2 = miller_index
  return (
    (bc**2*h0**2 + h1*(2*ab*cc*h0 + ac**2*h1 - aa*cc*h1) -
     2*ab*ac*h1*h2 + ab**2*h2**2 -
     2*bc*(ac*h0*h1 + ab*h0*h2 - aa*h1*h2) -
     bb*(cc*h0**2 - 2*ac*h0*h2 + aa*h2**2))/
    (ac**2*bb - 2*ab*ac*bc + ab**2*cc + aa*(bc**2 - bb*cc)))

def exercise():
  unit_cell = uctbx.unit_cell([11,13,7,83,98,105])
  for miller_index in [(1,2,3),(-3,4,2),(-1,-2,1)]:
    dss1 = unit_cell.d_star_sq(miller_index)
    dss2 = d_star_sq(unit_cell, miller_index)
    print dss1, dss2

if (__name__ == "__main__"):
  exercise()
