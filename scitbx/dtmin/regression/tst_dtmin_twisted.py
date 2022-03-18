from __future__ import print_function
from __future__ import division
from scitbx.array_family import flex
from scitbx.dtmin.minimizer import Minimizer
from scitbx.dtmin.refinebase import RefineBase
from scitbx.dtmin.reparams import Reparams
from scitbx.dtmin.bounds import Bounds
import math

class RefineTG(RefineBase):
  def __init__(self, start_x):
    RefineBase.__init__(self)
    self.start_x = start_x
    self.xy = start_x
    self.s11 = 1.0
    self.s12 = 1.2
    self.s22 = 2.0
    self.twist = 0.5

  def get_macrocycle_parameters(self):
    return self.xy

  def set_macrocycle_parameters(self, newx):
    self.xy = newx

  def set_macrocycle_protocol(self, macrocycle_protocol):
    if macrocycle_protocol == ["all"]:
      self.nmp = 2

  def macrocycle_parameter_names(self):
    return ["x", "y"]

  def macrocycle_large_shifts(self):
    return [1., 1.]

  def target(self):
    return twisted_gauss2d0(self.xy, self.s11, self.s12, self.s22, self.twist)

  def target_gradient(self):
    f = self.target()
    grad_x = analytic_grad_x(self.xy, self.s11, self.s12, self.s22, self.twist)
    grad_y = analytic_grad_y(self.xy, self.s11, self.s12, self.s22, self.twist)
    g = flex.double([grad_x, grad_y])
    return (f, g)

  def target_gradient_hessian(self):
    (f,g) = self.target_gradient()
    h = flex.double(self.nmp * self.nmp, 0)
    h.reshape(flex.grid(self.nmp, self.nmp))
    h[0,0] = analytic_curv_xx(self.xy, self.s11, self.s12, self.s22, self.twist)
    h[1,1] = analytic_curv_yy(self.xy, self.s11, self.s12, self.s22, self.twist)
    h[0,1] = analytic_curv_xy(self.xy, self.s11, self.s12, self.s22, self.twist)
    h[1,0] = h[0,1]
    return (f,g,h,False)

  def current_statistics(self):
    print("x,f: " + str(tuple(self.get_macrocycle_parameters())) + " " + str(self.target()))

#  def final_statistics(self):
#    print(self.start_x)
#    print(tuple(self.get_macrocycle_parameters()), "final")

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

def gauss2d0(xy, s11, s12, s22):
  (x,y) = xy
#  return -math.log(1/math.sqrt(4*math.pi**2*(s11*s22-s12**2))
#           *math.exp(-(s22*x**2-2*s12*x*y+s11*y**2)/(2*(s11*s22-s12**2))))
  return -math.log(1/math.sqrt(4*math.pi**2*(s11*s22-s12**2))) + \
           (s22*x**2-2*s12*x*y+s11*y**2)/(2*(s11*s22-s12**2))

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
  """"""

  x = [5,7]
  print(x, "start")

  # create inputs for the minimizer's run method
  refine_tg = RefineTG(start_x=x) # refine object
  macro1 = ["all"]                # protocol for the first macrocycle
  protocol = [macro1]             # overall minimization protocol
  ncyc = 50                       # maximum number of microcycles per macrocycle
  minimizer_type = "bfgs"         # minimizer, bfgs or newton
  study_params = True            # flag for calling studyparams procedure

  #create the minimizer object
  minimizer = Minimizer(output_level=0) # 0 for MUTE output see Minimizer.py

  #run the minimization
  minimizer.run(refine_tg, protocol, ncyc, minimizer_type, study_params)

if (__name__ == "__main__"):
  run()
