from __future__ import print_function
from __future__ import division
from scitbx.array_family import flex
from scitbx.dtmin.minimizer import Minimizer
from scitbx.dtmin.refinebase import RefineBase
from scitbx.dtmin.reparams import Reparams
from scitbx.dtmin.bounds import Bounds

class RefineTG(RefineBase):
  def __init__(self, start_x):
    RefineBase.__init__(self)
    self.start_x = start_x
    self.xy = start_x

  def target(self):
    """The simplest two variable quadratic function"""
    x = self.xy[0]
    y = self.xy[1]
    return x**2 + y**2

  def get_macrocycle_parameters(self):
    return self.xy

  def set_macrocycle_parameters(self, newx):
    self.xy = newx

#  def target_gradient(self):
#    f = self.target()
#    grad_x = 2*self.xy[0]
#    grad_y = 2*self.xy[1]
#    g = flex.double([grad_x, grad_y])
#    return (f, g)

#  def target_gradient_hessian(self):
#    (f,g) = self.target_gradient()
#    h = flex.double(self.nmp * self.nmp, 0)
#    h.reshape(flex.grid(self.nmp, self.nmp))
#    h[0,0] = h[1,1] = 2.
#    h[0,1] = h[1,0] = 0.
#    return (f,g,h,True)

  def macrocycle_large_shifts(self):
    return [3., 3.]

  def set_macrocycle_protocol(self, macrocycle_protocol):
    if macrocycle_protocol == ["all"]:
      self.nmp = 2

  def macrocycle_parameter_names(self):
    return ["x", "y"]

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

  def initial_statistics(self):
    #see min_logging.py for logging functions that roughly mirror those in
    #phasertng's minimiser.
    self.log_tab(1,0,"Demonstrate dtmin logging call")

  def current_statistics(self):
    print("x,f: " + str(tuple(self.get_macrocycle_parameters())) + " " + str(self.target()))

def run():
  """
  Run the dtmin minimizer with only the target() implemented.
  (no target_gradient or target_gradient_hessian).
  The target function is simply x**2 + y**2.
  Gradients are done by finite difference and the hessian is estimated
  by one over large_shifts squared. (This is handled by the minimizer)
  There are commented out target_gradient and target_gradient_hessian
  implementations should you wish play around.
  """

  x = [5,8]
  print(x, "start")

  # create inputs for the minimizer's run method
  refineTG = RefineTG(start_x=x) #refine object
  macro1 = ["all"]               # protocol for the first macrocycle
  protocol = [macro1, macro1, macro1] # overall minimization protocol
  ncyc = 50                      # maximum number of microcycles per macrocycle
  minimizer_type = "bfgs"        # minimizer, bfgs or newton
  study_params = False           # flag for calling studyparams procedure

  #create the minimizer object
  minimizer = Minimizer(output_level=0) # 0 for MUTE output see Minimizer.py

  #run the minimization
  minimizer.run(refineTG, protocol, ncyc, minimizer_type, study_params)

if (__name__ == "__main__"):
  run()
