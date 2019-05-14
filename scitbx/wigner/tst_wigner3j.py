from __future__ import absolute_import, division, print_function
from scitbx import wigner
import math

def tst1():
  wigner_fast = wigner.wigner3j_fast(30)

  result=math.sqrt(1.0/5.0)
  w= wigner.wigner3j(1,1,2,1,1,-2)
  value = wigner_fast.compute(1,1,2,1,1,-2)
  assert abs(result-value)<1e-6
  assert abs(w.get_value()-value)<1e-6

  result=math.sqrt(1.0/14.0)
  w=wigner.wigner3j(2,2,3,1,2,-3)
  value = wigner_fast.compute(2,2,3,1,2,-3)
  assert abs(result-value)<1e-6
  assert abs(w.get_value()-value)<1e-6

  result=math.sqrt(15.0/52003.0)*(-5.0)
  w=wigner.wigner3j(5, 8, 10, 3, -3, 0)
  value = wigner_fast.compute(5, 8, 10, 3, -3, 0)
  assert abs(result-value)<1e-6
  assert abs(w.get_value()-value)<1e-6

  w=wigner.wigner3j(5,8,10,2,4,-6)
  result=-(32.0/3.0)*math.sqrt(26.0/2860165.0)
  value = wigner_fast.compute(5,8,10,2,4,-6)
  assert abs(result-value)<1e-6
  assert abs(w.get_value()-value)<1e-6

if __name__=="__main__":
  tst1()
  print("OK")
