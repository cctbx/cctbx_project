"""
Approximating Log(BesselK[1/4,x]) using a Pade form.
Coefficients were fitted with gnuplot on numbers obtained from mathematica.

f(x) = (a0+a1*x+a2*x*x+a3*x*x*x+a4*x**4)/(1+b1*x+b2*x*x+b3*x*x*x+b4*x**4+b5*x**5)

The total range has been split in 3 regions:
0<x<=0.1 : Pade approx, 5+6 terms
0.1<x<=10: Pade approx, 4+5 terms
10<x<infty: Linear approx, 2 terms

No efforst have been made to ensure continuity, nor a rigorous error analyses has been
carried out. The error should be lower then a 5% over the x range 1E-4 to 100.
"""
from __future__ import absolute_import, division, print_function

import libtbx # for enumerate forward compatibility
from six.moves import zip

def log_besselk_1_4(x):
  result = None
  def funct(a_arr, b_arr, x):
    """ This computes an empirical pade approximation given coefficients """
    top = 0
    bottom = 0
    for n,coef in enumerate(a_arr):
      top+=coef*(x**n)
    for n,coef in enumerate(b_arr):
      bottom+=coef*(x**n)
    result = top/bottom
    return result

  if x < 0:
    x = -x
  if x < 0.1:
    a=[ 3.87347,99859.3, 1.12333e+08,7.44947e+09,-6.36247e+10]
    b=[ 1, 33946.7, 5.05513e+07, 5.19137e+09, -8.88956e+09,-2.56258e+11]
    result = funct(a,b,x)
  if (x >= 0.1) & (x<10):
    a = [1.67362 ,4.93003,-13.8875,-7.82595 ]
    b = [1,9.7653,7.14523 ,0.0341337,-0.000926896]
    result = funct(a,b,x)
  if x>=10:
    """ For large values, a linear approximation is good enough"""
    result = -1.10517 + -1.0109*x

  return result



def test(debug=False):
  # test various values
  xx = [ 0.01,    0.05,    0.1,       0.5,        5.0   ,  50]
  yy = [ 1.81901, 1.27752, 0.987739, -0.0404925, -5.5961, -51.7321] # from mathematica.
  for x,y in zip(xx,yy):
    f = log_besselk_1_4(x)
    if debug: print(x,y,f,100*abs((f-y)/y))
    assert abs((f-y)/y)*100.0 < 2.0

if __name__ == "__main__":
  import sys
  test('-debug' in sys.argv[1:])
  print('OK')
