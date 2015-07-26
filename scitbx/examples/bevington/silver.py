from __future__ import division
import scitbx.lbfgs
from scitbx.array_family import flex
from math import exp,pow,sqrt
from scitbx.matrix import sqr

data = """
Perform parameter fit three ways:
1) LBFGS without curvatures
2) Levenberg Marq. with dense matrix algebra
3) Levenberg Marq. with Eigen sparse matrix algebra

Take example from
Bevington, P.R. & Robinson, D.K., Data Reduction and Error Analysis for the Physical Sciences,
Third Edition.  New York: McGraw Hill, 2003.
Chapter 8: Least-Squares Fit to an Arbitrary Function.
The data describe beta-emission counts from the radioactive decay of silver.
Time Counts
15   775
30   479
45   380
60   302
75   185
90   157
105  137
120  119
135  110
150  89
165  74
180  61
195  66
210  68
225  48
240  54
255  51
270  46
285  55
300  29
315  28
330  37
345  49
360  26
375  35
390  29
405  31
420  24
435  25
450  35
465  24
480  30
495  26
510  28
525  21
540  18
555  20
570  27
585  17
600  17
615  14
630  17
645  24
660  11
675  22
690  17
705  12
720  10
735  13
750  16
765  9
780  9
795  14
810  21
825  17
840  13
855  12
870  18
885  10

Will fit these data to the following functional form (bi-exponential decay):
y(x) = a0 + a1 * exp( -x/ a3 ) + a2 * exp( -x/ a4)

Initial values for parameters, to be refined:
10 900 80 27 225
"""

raw_strings = data.split("\n")
data_index = raw_strings.index("Time Counts")
x_obs = flex.double([float(line.split()[0]) for line in raw_strings[data_index+1:data_index+60]])
y_obs = flex.double([float(line.split()[1]) for line in raw_strings[data_index+1:data_index+60]])
w_obs = 1./(y_obs)
initial = flex.double([10, 900, 80, 27, 225])
#initial = flex.double([10.4, 958.3, 131.4, 33.9, 205])

class lbfgs_biexponential_fit:
  def __init__(self, x_obs, y_obs, w_obs, initial):
    assert x_obs.size() == y_obs.size()
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.w_obs = w_obs
    self.n = len(initial)
    self.x = initial.deep_copy()
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x

  def print_step(pfh,message,target):
    print "%s %10.4f"%(message,target),
    print "["," ".join(["%10.4f"%a for a in pfh.x]),"]"

  def compute_functional_and_gradients(self):
    self.a = self.x

    y_calc = self.x[1] * flex.exp( -x_obs/self.x[3] ) + \
             self.x[2] * flex.exp( -x_obs/self.x[4] ) + \
             self.x[0]
    y_diff = self.y_obs - y_calc
    f = flex.sum(y_diff*y_diff*self.w_obs)
    self.print_step("LBFGS stp",f)
    g = flex.double(self.n,0)
    self.c = flex.double(self.n,0)
    for i in xrange(len(self.x_obs)):
      prefactor = -2. * self.w_obs[i] * y_diff[i]
      g[0] += prefactor
      g[1] += prefactor * exp(-self.x_obs[i]/self.x[3])
      g[2] += prefactor * exp(-self.x_obs[i]/self.x[4])
      g[3] += prefactor * self.x[1] * exp(-self.x_obs[i]/self.x[3]) * (self.x_obs[i]/(self.x[3]*self.x[3]))
      g[4] += prefactor * self.x[2] * exp(-self.x_obs[i]/self.x[4]) * (self.x_obs[i] * pow(self.x[4],-2))
      self.c[0] += 2 * self.w_obs[i]
      self.c[1] += 2 * self.w_obs[i] * (exp(-self.x_obs[i]/self.x[3]))**2
      self.c[2] += 2 * self.w_obs[i] * (exp(-self.x_obs[i]/self.x[4]))**2
      self.c[3] += 2 * self.w_obs[i] * (self.x[1] * exp(-self.x_obs[i]/self.x[3]) * (self.x_obs[i]/(self.x[3]*self.x[3])))**2
      self.c[4] += 2 * self.w_obs[i] * (self.x[2] * exp(-self.x_obs[i]/self.x[4]) * (self.x_obs[i] * pow(self.x[4],-2)))**2
    return f, g

def lbfgs_example(verbose):

  fit = lbfgs_biexponential_fit(x_obs=x_obs,y_obs=y_obs,w_obs=w_obs,initial=initial)
  print "------------------------------------------------------------------------- "
  print "       Initial and fitted coeffcients, and inverse-curvature e.s.d.'s"
  print "------------------------------------------------------------------------- "

  for i in range(initial.size()):

    print "%2d %10.4f %10.4f %10.4f"%(
           i, initial[i], fit.a[i], sqrt(2./fit.c[i]))

from scitbx.lstbx import normal_eqns
from scitbx.lstbx import normal_eqns_solving
class levenberg_helper(normal_eqns.non_linear_ls,normal_eqns.non_linear_ls_mixin):
  def __init__(pfh, initial_estimates):
    super(levenberg_helper, pfh).__init__(n_parameters=len(initial_estimates))
    pfh.x_0 = flex.double(initial_estimates)
    pfh.restart()
    pfh.counter = 0

  def set_data(pfh,x,y,w):
    pfh.xobs = x
    pfh.yobs = y
    pfh.wobs = w # set this equal to 1/variance

  def restart(pfh):
    pfh.x = pfh.x_0.deep_copy()
    pfh.old_x = None

  def step_forward(pfh):
    pfh.old_x = pfh.x.deep_copy()
    pfh.x += pfh.step()

  def step_backward(pfh):
    assert pfh.old_x is not None
    pfh.x, pfh.old_x = pfh.old_x, None

  def parameter_vector_norm(pfh):
    return pfh.x.norm()

  def fvec_callable(pfh,current):
    ycalc = flex.double(len(pfh.xobs))
    for iobs in xrange(len(pfh.xobs)):
      ycalc[iobs] = current[0] + \
                    current[1] * exp( -pfh.xobs[iobs]/ current[3]) + \
                    current[2] * exp( -pfh.xobs[iobs]/ current[4])
    y_diff = pfh.yobs - ycalc
    return y_diff

  def print_step(pfh,message,target):
    print "%s %10.4f"%(message,flex.sum(target)),
    print "["," ".join(["%10.4f"%a for a in pfh.x]),"]"

  def build_up(pfh, objective_only=False):
    if not objective_only: pfh.counter+=1
    residuals = pfh.fvec_callable(pfh.x)
    target = residuals*residuals*pfh.wobs
    pfh.reset()
    if objective_only:
      pfh.add_residuals(residuals, weights=pfh.wobs)
    else:
      pfh.print_step("LM  dense",target)
      for iobs in xrange(len(pfh.xobs)):
        jacobian_one_row = flex.double(len(pfh.x))
        prefactor = 1.
        jacobian_one_row [0] = prefactor
        jacobian_one_row [1] = prefactor * exp(-pfh.xobs[iobs]/pfh.x[3])
        jacobian_one_row [2] = prefactor * exp(-pfh.xobs[iobs]/pfh.x[4])
        jacobian_one_row [3] = prefactor * pfh.x[1] * exp(-pfh.xobs[iobs]/pfh.x[3]) * (pfh.xobs[iobs]/(pfh.x[3]*pfh.x[3]))
        jacobian_one_row [4] = prefactor * pfh.x[2] * exp(-pfh.xobs[iobs]/pfh.x[4]) * (pfh.xobs[iobs] * pow(pfh.x[4],-2))
        pfh.add_equation(-residuals[iobs], jacobian_one_row, pfh.wobs[iobs])

class dense_worker(object):

  def __init__(self,x_obs,y_obs,w_obs,initial):
    self.counter = 0
    self.x = initial.deep_copy()
    self.helper = levenberg_helper(initial_estimates = self.x)
    self.helper.set_data(x_obs,y_obs,w_obs)
    self.helper.restart()
    iterations = normal_eqns_solving.levenberg_marquardt_iterations(
               non_linear_ls = self.helper,
               n_max_iterations = 5000,
               track_all=True,
               step_threshold = 0.0001
    )
    ###### get esd's
    self.helper.build_up()
    upper = self.helper.step_equations().normal_matrix_packed_u()
    nm_elem = flex.double(25)
    self.c = flex.double(5)
    ctr = 0
    for x in xrange(5):
      x_0 = ctr
      for y in xrange(4,x-1,-1):
        nm_elem[ 5*x+y ] = upper[x_0+(y-x)]
        ctr += 1
        if x!= y:
          nm_elem[ 5*y+x ] = upper[x_0+(y-x)]
        else:
          self.c[x]=upper[x_0+(y-x)]
    NM = sqr(nm_elem)
    self.helper.solve()
    #print list(self.helper.step_equations().cholesky_factor_packed_u())
    error_matrix = NM.inverse()
    self.error_diagonal = [error_matrix(a,a) for a in xrange(5)]
    print "End of minimization: Converged", self.helper.counter,"cycles"

  def show_summary(self):
    print "%d cycles"%self.counter

def levenberg_example(verbose):

  fit = dense_worker(x_obs=x_obs,y_obs=y_obs,w_obs=w_obs,initial=initial)
  print "-------------------------------------------------------------------------------------- "
  print " Initial and fitted parameters, full-matrix e.s.d.'s, and inverse-curvature e.s.d.'s"
  print "-------------------------------------------------------------------------------------- "

  for i in range(initial.size()):

    print "%2d %10.4f %10.4f %10.4f %10.4f"%(
      i, initial[i], fit.helper.x[i], sqrt(fit.error_diagonal[i]), sqrt(1./fit.c[i]))

from scitbx.examples.bevington import bevington_base_class as base_class
class eigen_helper(base_class,normal_eqns.non_linear_ls_mixin):
  def __init__(pfh, initial_estimates):
    super(eigen_helper, pfh).__init__(n_parameters=len(initial_estimates))
    pfh.x_0 = flex.double(initial_estimates)
    pfh.restart()
    pfh.counter = 0

  def restart(pfh):
    pfh.x = pfh.x_0.deep_copy()
    pfh.old_x = None

  def step_forward(pfh):
    pfh.old_x = pfh.x.deep_copy()
    pfh.x += pfh.step()

  def step_backward(pfh):
    assert pfh.old_x is not None
    pfh.x, pfh.old_x = pfh.old_x, None

  def parameter_vector_norm(pfh):
    return pfh.x.norm()

  def print_step(pfh,message,target):
    print "%s %10.4f"%(message,flex.sum(target)),
    print "["," ".join(["%10.4f"%a for a in pfh.x]),"]"

  def build_up(pfh, objective_only=False):
    if not objective_only: pfh.counter+=1
    pfh.reset()
    #print list(pfh.x),objective_only
    if not objective_only:
      xresid = pfh.fvec_callable(pfh.x)
      target = xresid*xresid*pfh.w_obs
      pfh.print_step("LM sparse",target)
    pfh.access_cpp_build_up_directly_eigen_eqn(objective_only, current_values = pfh.x)

class eigen_worker(object):

  def __init__(self,x_obs,y_obs,w_obs,initial):
    self.counter = 0
    self.x = initial.deep_copy()
    self.helper = eigen_helper(initial_estimates = self.x)
    self.helper.set_cpp_data(x_obs,y_obs,w_obs)
    self.helper.w_obs = w_obs
    self.helper.restart()
    iterations = normal_eqns_solving.levenberg_marquardt_iterations_encapsulated_eqns(
               non_linear_ls = self.helper,
               n_max_iterations = 5000,
               track_all=True,
               step_threshold = 0.0001
    )
    ###### get esd's
    self.helper.build_up()
    self.error_diagonal = self.helper.solve_returning_error_mat_diagonal()
    print "End of minimization: Converged", self.helper.counter,"cycles"

  def show_summary(self):
    print "%d cycles"%self.counter

def eigen_example(verbose):

  fit = eigen_worker(x_obs=x_obs,y_obs=y_obs,w_obs=w_obs,initial=initial)
  print "----------------------------------------------------------------- "
  print "       Initial and fitted parameters & full-matrix e.s.d.'s"
  print "----------------------------------------------------------------- "

  for i in range(initial.size()):

    print "%2d %10.4f %10.4f %10.4f"%(
          i, initial[i], fit.helper.x[i], sqrt(fit.error_diagonal[i]) )

if (__name__ == "__main__"):
  verbose=True
  print "\n LBFGS:"
  lbfgs_example(verbose)
  print "\n DENSE MATRIX:"
  levenberg_example(verbose)
  print "\n SPARSE MATRIX:"
  eigen_example(verbose)
'''
  To do list:
  is there a quick calculation of the diag elements?  try to implement the algorithm & see
  implement an eigen version
  refactor this code, implement a common class for expressing
    1) residual terms
    2) gradient terms
    3) curvature approximations
    4) diagonal elements for full matrix
    5) diagaonal approximatinos in the eigen case
'''
