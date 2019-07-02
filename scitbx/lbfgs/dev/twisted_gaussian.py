from __future__ import absolute_import, division, print_function
from scitbx import lbfgs as scitbx_lbfgs
from scitbx.array_family import flex
from libtbx import adopt_init_args
import random
import math
import sys
from six.moves import range

def gauss2d0(xy, s11, s12, s22):
  (x,y) = xy
  return -math.log(1/math.sqrt(4*math.pi**2*(s11*s22-s12**2))
           *math.exp(-(s22*x**2-2*s12*x*y+s11*y**2)/(2*(s11*s22-s12**2))))

def twisted_gauss2d0(xy, s11, s12, s22, twist):
  (x,y) = xy
  arg = twist*math.sqrt(x**2+y**2)
  c = math.cos(arg)
  s = math.sin(arg)
  xt = x*c - y*s
  yt = y*c + x*s
  return gauss2d0((xt,yt), s11, s12, s22)

Cos = math.cos
Sin = math.sin
Sqrt = math.sqrt
Pi = math.pi

def analytic_grad_x(xy, s11, s12, s22, twist):
  (x,y) = xy
  if (x == 0 and y == 0): return 0.
  return (
         (-2*s12*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (Cos(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (Cos(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)))/
         (2.*(-s12**2 + s11*s22)))

def analytic_grad_y(xy, s11, s12, s22, twist):
  (x,y) = xy
  if (x == 0 and y == 0): return 0.
  return (
         (-2*s12*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (Cos(twist*Sqrt(x**2 + y**2)) +
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (Cos(twist*Sqrt(x**2 + y**2)) +
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)))/
         (2.*(-s12**2 + s11*s22)))

def analytic_curv_xx(xy, s11, s12, s22, twist):
  (x,y) = xy
  if (x == 0 and y == 0): return None
  return (
         (-2*s12*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 +
              (twist**2*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 +
              (twist**2*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(Cos(twist*Sqrt(x**2 + y**2)) -
               (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
               (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))**2
             + 2*s11*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x**3*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (3*twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 -
              (twist**2*x**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x**3*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (3*twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 -
              (twist**2*x**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           4*s12*(Cos(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            ((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/
                Sqrt(x**2 + y**2) + Sin(twist*Sqrt(x**2 + y**2)) -
               (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))**2)
          /(2.*(-s12**2 + s11*s22)))

def analytic_curv_yy(xy, s11, s12, s22, twist):
  (x,y) = xy
  if (x == 0 and y == 0): return None
  return (
         (-2*s12*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*y**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 +
              (twist**2*y**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*y**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 +
              (twist**2*y**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*y**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*y**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*y**3*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*y**3*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (3*twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/
                  Sqrt(x**2 + y**2)) - Sin(twist*Sqrt(x**2 + y**2)) -
               (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))**2\
            - 4*s12*(-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/
                 Sqrt(x**2 + y**2)) - Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            (Cos(twist*Sqrt(x**2 + y**2)) +
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*(Cos(twist*Sqrt(x**2 + y**2)) +
               (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
               (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))**2
           )/(2.*(-s12**2 + s11*s22)))

def analytic_curv_xy(xy, s11, s12, s22, twist):
  (x,y) = xy
  if (x == 0 and y == 0): return None
  return (
         (2*s11*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 -
              (twist**2*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            (-((twist*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/
                 (x**2 + y**2)**1.5) -
              (twist**2*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) +
              (twist*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 -
              (twist**2*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(y*Cos(twist*Sqrt(x**2 + y**2)) +
              x*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 +
              (twist**2*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(x*Cos(twist*Sqrt(x**2 + y**2)) -
              y*Sin(twist*Sqrt(x**2 + y**2)))*
            ((twist*x*y**2*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2)**1.5 -
              (twist**2*x**2*y*Cos(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*x*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              (twist*x**2*y*Sin(twist*Sqrt(x**2 + y**2)))/
               (x**2 + y**2)**1.5 +
              (twist**2*x*y**2*Sin(twist*Sqrt(x**2 + y**2)))/(x**2 + y**2) -
              (twist*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s22*(Cos(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            (-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(-((twist*y**2*Cos(twist*Sqrt(x**2 + y**2)))/
                 Sqrt(x**2 + y**2)) - Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            ((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) +
              Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) -
           2*s12*(Cos(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*x**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            (Cos(twist*Sqrt(x**2 + y**2)) +
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)) +
           2*s11*((twist*x**2*Cos(twist*Sqrt(x**2 + y**2)))/
               Sqrt(x**2 + y**2) + Sin(twist*Sqrt(x**2 + y**2)) -
              (twist*x*y*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2))*
            (Cos(twist*Sqrt(x**2 + y**2)) +
              (twist*x*y*Cos(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2) -
              (twist*y**2*Sin(twist*Sqrt(x**2 + y**2)))/Sqrt(x**2 + y**2)))/
         (2.*(-s12**2 + s11*s22)))

def finite_grad_x(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = twisted_gauss2d0((x-eps,y), s11, s12, s22, twist)
  tp = twisted_gauss2d0((x+eps,y), s11, s12, s22, twist)
  return (tp-tm)/(2*eps)

def finite_grad_y(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = twisted_gauss2d0((x,y-eps), s11, s12, s22, twist)
  tp = twisted_gauss2d0((x,y+eps), s11, s12, s22, twist)
  return (tp-tm)/(2*eps)

def finite_curv_xx(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = finite_grad_x((x-eps,y), s11, s12, s22, twist, eps)
  tp = finite_grad_x((x+eps,y), s11, s12, s22, twist, eps)
  return (tp-tm)/(2*eps)

def finite_curv_yy(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = finite_grad_y((x,y-eps), s11, s12, s22, twist, eps)
  tp = finite_grad_y((x,y+eps), s11, s12, s22, twist, eps)
  return (tp-tm)/(2*eps)

def finite_curv_xy(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = finite_grad_x((x,y-eps), s11, s12, s22, twist, eps)
  tp = finite_grad_x((x,y+eps), s11, s12, s22, twist, eps)
  return (tp-tm)/(2*eps)

def finite_curv_yx(xy, s11, s12, s22, twist, eps=1.e-6):
  (x,y) = xy
  tm = finite_grad_y((x-eps,y), s11, s12, s22, twist, eps)
  tp = finite_grad_y((x+eps,y), s11, s12, s22, twist, eps)
  return (tp-tm)/(2*eps)

def verify_derivatives(n=5, s11=1, s12=1.2, s22=2, twist=0.5, verbose=0):
  for ix in range(-n,n+1):
    for iy in range(-n,n+1):
      xy = [i/float(n) for i in (ix,iy)]
      if (0 or verbose):
        print("value: %5.3f" % twisted_gauss2d0(xy, s11, s12, s22, twist))
        print()
      for f,a in ((finite_grad_x,analytic_grad_x),
                  (finite_grad_y,analytic_grad_y)):
        fg = f(xy, s11, s12, s22, twist)
        ag = a(xy, s11, s12, s22, twist)
        if (0 or verbose):
          print("fg:", fg)
          print("ag:", ag)
          print()
        assert abs(fg-ag) < 1.e-5
      for f,a in ((finite_curv_xx,analytic_curv_xx),
                  (finite_curv_yy,analytic_curv_yy),
                  (finite_curv_xy,analytic_curv_xy),
                  (finite_curv_yx,analytic_curv_xy)):
        fc = f(xy, s11, s12, s22, twist)
        ac = a(xy, s11, s12, s22, twist)
        if (0 or verbose):
          print("fc:", fc)
          print("ac:", ac)
          print()
        if (xy != [0,0]):
          assert abs(fc-ac)/max(1,min(abs(fc),abs(ac))) < 1.e-3

class fortran_minimizer:

  def __init__(self, scitbx_minimizer):
    import numpy
    self.n = scitbx_minimizer.n()
    self.m = scitbx_minimizer.m()
    self.x = numpy.array(numpy.arange(self.n), numpy.float64)
    self.g = numpy.array(numpy.arange(self.n), numpy.float64)
    self.diag = numpy.array(numpy.arange(self.n), numpy.float64)
    self.iprint = [1, 0]
    self.eps = 1.e-5 # convergence test
    self.xtol = scitbx_minimizer.xtol()
    size_w = self.n*(2*self.m+1)+2*self.m
    self.w = numpy.array(numpy.arange(size_w), numpy.float64)
    self.iflag = numpy.array([0], numpy.int32)

  def __call__(self, x, f, g, diag=None, diagco=False):
    for i,v in enumerate(x): self.x[i] = v
    for i,v in enumerate(g): self.g[i] = v
    if (diag is not None):
      for i,v in enumerate(diag): self.diag[i] = v
    from fortran_lbfgs import lbfgs as fortran_lbfgs
    fortran_lbfgs(
      self.n, self.m,
      self.x, f, self.g, diagco, self.diag,
      self.iprint, self.eps, self.xtol, self.w, self.iflag)
    for i,v in enumerate(self.x): x[i] = v
    for i,v in enumerate(self.g): g[i] = v
    if (diag is not None):
      for i,v in enumerate(self.diag): diag[i] = v

def fortran_lbfgs_run(target_evaluator,
                      max_calls=100000,
                      use_curvatures=False):
  ext = scitbx_lbfgs.ext
  scitbx_minimizer = ext.minimizer(target_evaluator.n)
  minimizer = fortran_minimizer(scitbx_minimizer)
  icall = 0
  requests_f_and_g = True
  requests_diag = use_curvatures
  while 1:
    if (not use_curvatures):
      assert not requests_diag
    x, f, g, d = target_evaluator(
      requests_f_and_g=requests_f_and_g,
      requests_diag=requests_diag)
    if (requests_diag):
      print("x,f,d:", tuple(x), f, tuple(d))
    else:
      print("x,f:", tuple(x), f)
    sys.stdout.flush()
    sys.stderr.flush()
    minimizer(x, f, g, diag=d, diagco=use_curvatures)
    print("iflag:", minimizer.iflag[0])
    if (minimizer.iflag[0] <= 0): break
    requests_f_and_g = minimizer.iflag[0] == 1
    requests_diag = minimizer.iflag[0] == 2
    if (requests_f_and_g):
      icall += 1
      if (icall > max_calls): break
  minimizer.n_calls = icall
  return minimizer

def lbfgs_run(target_evaluator,
              min_iterations=0,
              max_iterations=None,
              traditional_convergence_test=1,
              use_curvatures=False):
  ext = scitbx_lbfgs.ext
  minimizer = ext.minimizer(target_evaluator.n)
  minimizer.error = None
  if (traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(target_evaluator.n)
  else:
    raise RuntimeError
    is_converged = ext.drop_convergence_test(min_iterations)
  try:
    icall = 0
    requests_f_and_g = True
    requests_diag = use_curvatures
    while 1:
      if (requests_f_and_g):
        icall += 1
      x, f, g, d = target_evaluator(
        requests_f_and_g=requests_f_and_g,
        requests_diag=requests_diag)
      if (requests_diag):
        print("x,f,d:", tuple(x), f, tuple(d))
      else:
        print("x,f:", tuple(x), f)
      if (use_curvatures):
        if (d is None): d = flex.double(x.size())
        have_request = minimizer.run(x, f, g, d)
      else:
        have_request = minimizer.run(x, f, g)
      if (have_request):
        requests_f_and_g = minimizer.requests_f_and_g()
        requests_diag = minimizer.requests_diag()
        continue
      assert not minimizer.requests_f_and_g()
      assert not minimizer.requests_diag()
      if (traditional_convergence_test):
        if (minimizer.iter() >= min_iterations and is_converged(x, g)): break
      else:
        if (is_converged(f)): break
      if (max_iterations is not None and minimizer.iter() >= max_iterations):
        break
      if (use_curvatures):
        have_request = minimizer.run(x, f, g, d)
      else:
        have_request = minimizer.run(x, f, g)
      if (not have_request): break
      requests_f_and_g = minimizer.requests_f_and_g()
      requests_diag = minimizer.requests_diag()
  except RuntimeError as e:
    minimizer.error = str(e)
  minimizer.n_calls = icall
  return minimizer

class twisted_gaussian_minimizer:

  def __init__(self, x, s11=1, s12=1.2, s22=2, twist=0.5,
               min_iterations=0, max_iterations=10000):
    adopt_init_args(self, locals())
    self.x = flex.double(x)
    self.n = self.x.size()

  def run(self, use_fortran=0, use_curvatures=0):
    if (not use_fortran):
      self.minimizer = lbfgs_run(
        target_evaluator=self,
        min_iterations=self.min_iterations,
        max_iterations=self.max_iterations,
        use_curvatures=use_curvatures)
    else:
      self.minimizer = fortran_lbfgs_run(
        target_evaluator=self,
        use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=False)
    return self

  def __call__(self, requests_f_and_g, requests_diag):
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = twisted_gauss2d0(self.x, self.s11,self.s12,self.s22, self.twist)
      self.g = flex.double(
        (finite_grad_x(self.x, self.s11, self.s12, self.s22, self.twist),
         finite_grad_y(self.x, self.s11, self.s12, self.s22, self.twist)))
      self.d = None
    if (requests_diag):
      self.d = flex.double(
        (analytic_curv_xx(self.x, self.s11, self.s12, self.s22, self.twist),
         analytic_curv_yy(self.x, self.s11, self.s12, self.s22, self.twist)))
      self.df = flex.double(
        (finite_curv_xx(self.x, self.s11, self.s12, self.s22, self.twist),
         finite_curv_yy(self.x, self.s11, self.s12, self.s22, self.twist)))
      assert self.d.all_ne(0)
      print(tuple(self.df), "finite")
      print(tuple(self.d), "analytic")
      self.d = 1 / self.d
    return self.x, self.f, self.g, self.d

def run(scale=2, twist=0.5):
  random.seed(0)
  if ("--verify" in sys.argv[1:]):
    verify_derivatives()
  use_fortran = "--fortran" in sys.argv[1:]
  for iteration in range(100):
    x = [random.random()*scale for i in (0,1)]
    print(x, "start")
    for use_curvatures in (False, True):
      m = twisted_gaussian_minimizer(x=x, twist=twist).run(
        use_fortran=False,
        use_curvatures=use_curvatures)
      print(x)
      print(tuple(m.x), "final")
      if (use_fortran):
        mf = twisted_gaussian_minimizer(x=x, twist=twist).run(
          use_fortran=True,
          use_curvatures=use_curvatures)
        assert mf.x.all_eq(m.x)
        print(mf.minimizer.n_calls, m.minimizer.n_calls)
        assert mf.minimizer.n_calls+1 == m.minimizer.n_calls
      if (abs(m.x[0]) > 1.e-4 or abs(m.x[1]) > 1.e-4):
        print(tuple(m.x), "failure, use_curvatures="+str(use_curvatures))
      print("iter,exception:", m.minimizer.iter(), m.minimizer.error)
      print("n_calls:", m.minimizer.n_calls)
      assert m.minimizer.n_calls == m.minimizer.nfun()
      print()

if (__name__ == "__main__"):
  run()
