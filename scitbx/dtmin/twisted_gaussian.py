#This file is broadly comparable to scitbx/lbfgs/dev/twisted_gaussian.py
#so may be a useful as a way for someone familiar with scitbx.lbfgs to
#see how to do things in dtmin.

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from scitbx.array_family import flex
from scitbx.dtmin.minimizer import Minimizer
from scitbx.dtmin.reparams import Reparams
from scitbx.dtmin.bounds import Bounds
from scitbx.dtmin.refinebase import RefineBase

import math
import random

class RefineTG(RefineBase):
  def __init__(self, start_x, use_curvatures):
    RefineBase.__init__(self)
    self.start_x = start_x
    self.xy = start_x
    self.use_curvatures = use_curvatures
    self.s11 = 1.0
    self.s12 = 1.2
    self.s22 = 2.0
    self.twist = 0.5

  def target(self):
    return twisted_gauss2d0(self.xy, self.s11, self.s12, self.s22, self.twist)

  def get_macrocycle_parameters(self):
    return self.xy

  def set_macrocycle_parameters(self, newx):
    self.xy = newx

  def target_gradient(self):
    f = twisted_gauss2d0(self.xy, self.s11, self.s12, self.s22, self.twist)
    grad_x = analytic_grad_x(self.xy, self.s11, self.s12, self.s22, self.twist)
    grad_y = analytic_grad_y(self.xy, self.s11, self.s12, self.s22, self.twist)
    g = flex.double([grad_x, grad_y])
    return (f, g)

  def target_gradient_hessian(self):
    (f,g) = self.target_gradient()
    h = flex.double(self.nmp * self.nmp, 0)
    h.reshape(flex.grid(self.nmp, self.nmp))
    if self.use_curvatures == False:
      h[0,0] = h[1,1] = 0.
      h[0,1] = h[1,0] = 1.
    elif self.use_curvatures == True:
      h[0,0] = analytic_curv_xx(self.xy, self.s11, self.s12, self.s22, self.twist)
      h[1,1] = analytic_curv_yy(self.xy, self.s11, self.s12, self.s22, self.twist)
      h[0,1] = h[1,0] = 0.
    return (f,g,h,True)

  def macrocycle_large_shifts(self):
    return [5., 5.]

  def set_macrocycle_protocol(self, macrocycle_protocol):
    if macrocycle_protocol == ["all"]:
      self.nmp = 2
    else:
      self.nmp = 0

  def macrocycle_parameter_names(self):
    return ["parameter 1", "parameter 2"]

#  def reparameterize(self):
#    rep_x = Reparams(True, 5.)
#    rep_y = Reparams(True, 5.)
#    return [rep_x, rep_y]

#  def bounds(self):
#    bnd_x = Bounds()
#    bnd_y = Bounds()
#    bnd_x.on(-5,5)
#    bnd_y.on(-5,5)
#    return [bnd_x, bnd_y]

  def current_statistics(self):
    if self.use_curvatures == False:
      print("x,f: " + str(tuple(self.get_macrocycle_parameters())) + " " + str(self.target()))

  def final_statistics(self):
    print(self.start_x)
    print(tuple(self.get_macrocycle_parameters()), "final")
    if (abs(self.get_macrocycle_parameters()[0]) > 1.e-4 or abs(self.get_macrocycle_parameters()[1]) > 1.e-4):
      print(tuple(self.get_macrocycle_parameters()), "failure, use_curvatures="+str(self.use_curvatures))
    #print("iter,exception:", m.minimizer.iter(), m.minimizer.error)
    #print("n_calls:", m.minimizer.n_calls)
    print()

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

def run():
  # create inputs for the minimizer's run method
  macro1 = ["all"]            # protocol for the first macrocycle
  MACRO2 = ["all"]            # protocol for the second macrocycle
  protocol = [macro1, MACRO2] # overall minimization protocol
  ncyc = 50                   # maximum number of microcycles per macrocycle
  minimizer_type = "bfgs"     # minimizer, bfgs or newton
  study_params = False        # flag for calling studyparams procedure

  #run the minimization
  random.seed(0)
  for iteration in range(100):
    scale = 2
    x = [random.random()*scale for i in (0,1)]
    print(x, "start")
    for use_curvatures in (False, True):
      refineTG = RefineTG(start_x=x, use_curvatures=use_curvatures)
      minimizer = Minimizer(0) # 0 for MUTE output see Minimizer.py
      minimizer.run(refineTG, protocol, ncyc, minimizer_type, study_params)

if (__name__ == "__main__"):
  run()
