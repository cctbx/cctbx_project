from __future__ import absolute_import, division, print_function
import math
import scitbx.math
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import range

def halton_x(n=100):
  halton_object = scitbx.math.halton(1)
  x = flex.double()
  for ii in range(n):
    x.append( halton_object.nth_given_base(3,ii) )
  return x

x_precomputed = halton_x(40)

def log_p(z,sigz,centric_flag,n=40,lim=5.0):
  max_e = None
  max_s = None
  def experimental_array(x,z,sigz):
    dd = (z-x*x)*(z-x*x)/(2.0*sigz*sigz)
    sel = flex.bool(dd > 100)
    dd.set_selected( sel, 100)
    result = (1.0/(math.sqrt(2.0*math.pi)*sigz))*flex.exp( -dd )
    return result
  def experimental_single(x,z,sigz):
    dd = (z-x*x)*(z-x*x)/(2.0*sigz*sigz)
    if dd >= 100.0:
      dd = 100.0
    result = math.exp( -dd )*(1.0/(math.sqrt(2.0*math.pi)*sigz))
    return result

  def centric(x):
    result = math.sqrt(2.0/math.pi)*flex.exp( -x*x/2.0 )
    return result
  def centric_single(x):
    result = math.sqrt(2.0/math.pi)*math.exp( -x*x/2.0 )
    return result

  def acentric(x):
    result = 2.0*x*flex.exp( -x*x )
    return result
  def acentric_single(x):
    result = 2.0*x*math.exp( -x*x )
    return result

  def f_single(x,z,sigz,centric_flag):
    if centric_flag:
      result = -x*x/2.0 - (z-x*x)**2/(2*sigz*sigz)
    if not centric_flag:
      result = -x*x + math.log(x+1e-9)-(z-x*x)**2/(2*sigz*sigz)
    return -result

  def f_rat(x,z,sigz,centric_flag,level=5):
    result = f_single(x,z,sigz,centric_flag)-max_s
    result = (result -level*level)**2.0
    return result


  def golden_section_max(z,sigz,centric_flag,function,a=0,b=6,eps=1E-5,max_count=100):
    c   = (-1+math.sqrt(5))/2.0
    x1  = c*a+(1-c)*b
    fx1 = function(x1,z,sigz,centric_flag)
    x2  = (1-c)*a + c*b
    fx2 = function(x2,z,sigz,centric_flag)
    delta = b-a
    count = 0
    while delta> eps:
      if fx1 < fx2:
        b   = x2
        x2  = x1
        fx2 = fx1
        x1  = c*a+(1-c)*b
        fx1 = function(x1,z,sigz,centric_flag)
        delta = b-a
      else:
        a   = x1
        x1  = x2
        fx1 = fx2
        x2  = (1-c)*a+c*b
        fx2 = function(x2,z,sigz,centric_flag)
        delta = b-a
      count+=1
      if count >= max_count:
        delta=eps
      #print a,b,function(a,z,sigz,centric_flag),function(b,z,sigz,centric_flag)

    return (a+b)/2.0

  result = 0
  x = None
  if n == 40:
    x = x_precomputed
  else:
    x = halton_x(n)
  jac = lim
  if centric_flag == None: #absent
    result = experimental_single(0,z,sigz)
  else:
    max_e = golden_section_max(z,sigz,centric_flag,f_single)
    max_s = f_single(max_e,z,sigz,centric_flag)
    low_e = golden_section_max(z,sigz,centric_flag,f_rat,0,max_e)
    high_e = golden_section_max(z,sigz,centric_flag,f_rat,max_e,6)#math.sqrt(abs(z))+5.0*sigz )
    if centric_flag==True:
      x = x*(high_e-low_e)+low_e
      jac = (high_e-low_e)
      result = experimental_array(x,z,sigz)*centric(x)
      result = flex.sum(result)*jac/n
    if centric_flag==False:
      x = x*(high_e-low_e)+low_e
      jac = (high_e-low_e)
      result = experimental_array(x,z,sigz)*acentric(x)
      result = flex.sum(result)*jac/n
  return -math.log(result+1E-12)



def test():
  sigz = 0.000050
  for ii in range(600):
    z = ii/50.0 +0.3
    pac = log_p(z=z,sigz=sigz,centric_flag=False,n=20)
    pcc = log_p(z=z,sigz=sigz,centric_flag=True,n=20)
    tac = 0
    tcc = 0
    if z >=0:
     tac  =  z
     tmp  =  math.sqrt( 1.0 / (2.0*math.pi*max(z,1E-6)) )
     if z <=0:
       tmp = 0.0
     tcc  =  -math.log(math.exp(-z/2.0)*tmp+1E-12)
    assert approx_equal(pac,tac,eps=1e-1)
    assert approx_equal(pcc,tcc,eps=1e-1)
  print("OK")


if __name__ == "__main__":
  test()
