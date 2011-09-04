"""
Inspired by Press, Teukolsky, Vetterling, Flannery (1992).
Numerical Recipes in C.  Cambridge University Press.
Chapter 3.3.  Cubic Spline Interpolation.
"""

class cubic_spline:  #Should be re-implemented in C++/Boost Python
  def __init__(self,xarr,yarr,lo_deriv1=None,hi_deriv1=None):
    from scitbx.array_family import flex
    assert len(xarr)==len(yarr)
    tempu = flex.double(len(xarr))
    self.deriv2 = flex.double(len(xarr))

    if lo_deriv1==None:
      self.deriv2[0] = 0.0; tempu[0] = 0.0
    else:
      self.deriv2[0] = -0.5;
      tempu[0] = (3./(xarr[1]-xarr[0]))*((yarr[1]-yarr[0])/(xarr[1]-xarr[0])-lo_deriv1)

    for i in xrange(1,len(xarr)-1):
      sig = (xarr[i]-xarr[i-1])/(xarr[i+1]-xarr[i-1])
      p = sig*self.deriv2[i-1]+2.0
      self.deriv2[i] = (sig-1.0)/p
      tempu[i] = (yarr[i+1]-yarr[i])/(xarr[i+1]-xarr[i]) - \
                 (yarr[i]-yarr[i-1])/(xarr[i]-xarr[i-1])
      tempu[i] = (6.*tempu[i] / (xarr[i+1]-xarr[i-1]) - sig *tempu[i-1])/p

    NX = len(xarr)
    if hi_deriv1==None:
      qn = 0.0; un = 0.0
    else:
      qn = 0.5;
      un = (3./(xarr[NX-1]-xarr[NX-2]))* \
           (hi_deriv1-(yarr[NX-1]-yarr[NX-2])/(xarr[NX-1]-xarr[NX-2]))

    self.deriv2[NX-1] = (un-qn*tempu[NX-2])/(qn*self.deriv2[NX-2]+1.0)
    for k in xrange(NX-2,-1,-1):
      self.deriv2[k]=self.deriv2[k]*self.deriv2[k+1]+tempu[k]
    self.xarr = xarr
    self.yarr = yarr

  def get_y(self,x):
    klo = 0
    khi = len(self.xarr)-1
    while (khi-klo > 1):
      k = khi+klo >> 1
      if self.xarr[k] > x: khi = k
      else: klo = k
    h = self.xarr[khi]-self.xarr[klo]
    assert h !=0
    a = (self.xarr[khi]-x) / h
    b = (x - self.xarr[klo])/ h
    return a * self.yarr[klo] + b * self.yarr[khi] + \
           (( a*a*a -a ) * self.deriv2[klo] + ( b*b*b - b )* self.deriv2[khi]) * (h*h)/6.

def test_case_1(plt):
  import math
  xarr = xrange(0,11)
  yarr = [math.sin(x) for x in xarr]
  CS = cubic_spline(xarr,yarr)
  xsam = [(x+.5)/10. for x in xrange(-1,101)]
  ysam = [CS.get_y(x) for x in xsam]
  plt.plot(xarr,yarr,"r+")
  plt.plot(xsam,ysam,"g.")
  plt.show()

def test_case_2(plt):
  import math
  xarr = xrange(-4,5)
  yarr = [0., 0.15, 1.12, 2.36, 2.36, 1.46, .49, .06, 0.]
  CS = cubic_spline(xarr,yarr,lo_deriv1=0.,hi_deriv1=0.)
  xsam = [(x+.5)/10. for x in xrange(-41,41)]
  ysam = [CS.get_y(x) for x in xsam]
  plt.plot(xarr,yarr,"r+")
  plt.plot(xsam,ysam,"g.")
  plt.show()

if __name__=="__main__":
  from matplotlib import pyplot as plt
  test_case_1(plt)
  test_case_2(plt)
