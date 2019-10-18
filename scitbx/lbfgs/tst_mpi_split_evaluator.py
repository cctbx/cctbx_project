from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.lbfgs import core_parameters, termination_parameters
from scitbx.lbfgs import exception_handling_parameters, ext
from scitbx.array_family import flex
import scitbx

"""mpi_split_evaluator_run(), supports an LBFGS parameter optimization scenario where
the target (functional and gradients) are significantly rate limiting, and moreover where
the requisite terms of f and g can be load balanced by distributing the data over parallel
evaluator instances, each of which can be handled by an MPI worker rank.  Rank 0 then
performs a simple MPI.reduce sum to obtain the full f and g.  There has been
no low-level redesign to support MPI.  In particular, the ext.minimizer is
run (wastefully) by every worker rank, using the same x parameters, f, and g. A simple
working example is given."""

# based on scitbx/lbfgs/__init__.py, run_c_plus_plus
def mpi_split_evaluator_run(target_evaluator,
                    termination_params=None,
                    core_params=None,
                    exception_handling_params=None,
                    log=None,
                    #---> Insertion starts
                    gradient_only=False,
                    line_search=True):
                    #<--- Insertion ends
  """The supported scenario is that each MPI worker rank has a target evaluator
  that has part of the data.  Each rank calculates a bit of the functional and
  gradients, but then mpi reduce is used to sum them all up.  There has been
  no low-level redesign to support MPI.  In particular, the ext.minimizer is
  run (wastefully) by every worker rank, using the same data.  It is assumed that
  the calculation of compute_functional_and_gradients() is overwhelmingly the rate
  limiting step, and that is what MPI parallelism is intended to distribute here."""
  from libtbx.mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  if (termination_params is None):
    termination_params = termination_parameters()
  if (core_params is None):
    core_params = core_parameters()
  if (exception_handling_params is None):
    exception_handling_params = exception_handling_parameters()
  x = target_evaluator.x
  if (log is not None):
    print("lbfgs minimizer():", file=log)
    print("  x.size():", x.size(), file=log)
    print("  m:", core_params.m, file=log)
    print("  maxfev:", core_params.maxfev, file=log)
    print("  gtol:", core_params.gtol, file=log)
    print("  xtol:", core_params.xtol, file=log)
    print("  stpmin:", core_params.stpmin, file=log)
    print("  stpmax:", core_params.stpmax, file=log)
    print("lbfgs traditional_convergence_test:", \
      termination_params.traditional_convergence_test, file=log)
  minimizer = ext.minimizer(
    x.size(),
    core_params.m,
    core_params.maxfev,
    core_params.gtol,
    core_params.xtol,
    core_params.stpmin,
    core_params.stpmax)
  if (termination_params.traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(
      x.size(),
      termination_params.traditional_convergence_test_eps)
  else:
    is_converged = ext.drop_convergence_test(
      n_test_points=termination_params.drop_convergence_test_n_test_points,
      max_drop_eps=termination_params.drop_convergence_test_max_drop_eps,
      iteration_coefficient
        =termination_params.drop_convergence_test_iteration_coefficient)
  callback_after_step = getattr(target_evaluator, "callback_after_step", None)
  diag_mode = getattr(target_evaluator, "diag_mode", None)
  if (diag_mode is not None): assert diag_mode in ["once", "always"]
  f_min, x_min = None, None
  f, g = None, None
  try:
    while 1:
      if (diag_mode is None):
        #XXX Only the diag_mode==None case is currently implemented, just as example
        f_term, g_term = target_evaluator.compute_functional_and_gradients()
        f_total = comm.reduce(f_term, MPI.SUM, 0)
        g_total = comm.reduce(g_term, MPI.SUM, 0)
        if rank==0: transmit = (f_total,g_total)
        else: transmit = None
        f, g = comm.bcast(transmit, root=0)
        if False and rank==0: # for debug
          print ("%s %10.4f"%("MPI stp",f),"["," ".join(["%10.4f"%a for a in x]),"]")
        d = None
      else:
        f_term, g_term, c_term = target_evaluator.compute_functional_gradients_curvatures()
        f_total = comm.reduce(f_term, MPI.SUM, 0)
        g_total = comm.reduce(g_term, MPI.SUM, 0)
        c_total = comm.reduce(c_term, MPI.SUM, 0)
        transmit = None
        if rank == 0:
          transmit = f_total, g_total, c_total
        f, g, curv = comm.bcast(transmit, root=0)

        sel = (g != 0)
        curv.set_selected(~sel, 1000)  # however, if we've decided to not refine a certain parameter, we
        # can indicate this to LBFGS by setting the gradient to zero.
        # Then we can set the curvature to an arbitrary positive value that
        # will be tested for positivity but otherwise ignored.

        assert curv.select(sel).all_gt(0)
        d = 1 / curv
        if (diag_mode == "once"):
          diag_mode = None
      if (f_min is None):
        if (not termination_params.traditional_convergence_test):
          is_converged(f)
        f_min, x_min = f, x.deep_copy()
      elif (f_min > f):
        f_min, x_min = f, x.deep_copy()
      if (log is not None):
        print("lbfgs minimizer.run():" \
          " f=%.6g, |g|=%.6g, x_min=%.6g, x_mean=%.6g, x_max=%.6g" % (
          f, g.norm(), flex.min(x), flex.mean(x), flex.max(x)), file=log)
      if (d is None):
        #---> Insertion starts
        if (minimizer.run(x, f, g, gradient_only,line_search)): continue
        #<--- Insertion ends
      else:
        #---> Insertion starts
        if (minimizer.run(x, f, g, d, gradient_only,line_search)): continue
        #<--- Insertion ends
      if (log is not None):
        print("lbfgs minimizer step", file=log)
      if (callback_after_step is not None):
        if (callback_after_step(minimizer) is True):
          if (log is not None):
            print("lbfgs minimizer stop: callback_after_step is True", file=log)
          break
      if (termination_params.traditional_convergence_test):
        if (    minimizer.iter() >= termination_params.min_iterations
            and is_converged(x, g)):
          if (log is not None):
            print("lbfgs minimizer stop: traditional_convergence_test", file=log)
          break
      else:
        if (is_converged(f)):
          if (log is not None):
            print("lbfgs minimizer stop: drop_convergence_test", file=log)
          break
      if (    termination_params.max_iterations is not None
          and minimizer.iter() >= termination_params.max_iterations):
        if (log is not None):
          print("lbfgs minimizer stop: max_iterations", file=log)
        break
      if (    termination_params.max_calls is not None
          and minimizer.nfun() > termination_params.max_calls):
        if (log is not None):
          print("lbfgs minimizer stop: max_calls", file=log)
        break
      if (d is None):
        #---> Insertion starts
        if (not minimizer.run(x, f, g, gradient_only,line_search)): break
        #<--- Insertion ends
      else:
        #---> Insertion starts
        if (not minimizer.run(x, f, g, d, gradient_only,line_search)): break
        #<--- Insertion ends
  except RuntimeError as e:
    minimizer.error = str(e)
    if (log is not None):
      print("lbfgs minimizer exception:", str(e), file=log)
    if (x_min is not None):
      x.clear()
      x.extend(x_min)
    error_classification = exception_handling_params.filter(
      minimizer.error, x.size(), x, g)
    if (error_classification > 0):
      raise
    elif (error_classification < 0):
      minimizer.is_unusual_error = True
    else:
      minimizer.is_unusual_error = False
  else:
    minimizer.error = None
    minimizer.is_unusual_error = None
  if (log is not None):
    print("lbfgs minimizer done.", file=log)
  return minimizer

class simple_quadratic(object):
  def __init__(self):
    self.datax = flex.double(range(-15,17))
    self.datay = flex.double([20,15,18,12,10, 10,5,5,1,2, -3,-1,-4,-5,-4,
                          -6,-4,-6,-4,-4,  -4,-5,-1,0,-1,  1,5,4,9,10, 13,15])
    abc = 0.1,-0.3,-5.0 # The expected parameters, y = a*x*x + b*x +  c
    self.n = 3
    self.x = flex.double([1,1,1])#lay out the parameter estimates.

  def run(self):
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-3,
        max_calls=1000)
    )
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),"["," ".join(["%10.4f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.a = self.x
    residuals = self.datay - self.a[0]*self.datax*self.datax - self.a[1]*self.datax - self.a[2]
    f = flex.sum( 0.5 * residuals * residuals )
    g = flex.double(self.n)
    dR_da = -self.datax*self.datax
    dR_db = -self.datax
    dR_dc = flex.double(len(self.datax),-1)
    g[0] = flex.sum( residuals * dR_da )
    g[1] = flex.sum( residuals * dR_db )
    g[2] = flex.sum( residuals * dR_dc )
    # self.print_step("LBFGS stp",f)
    return f,g

class mpi_quadratic(simple_quadratic):

  def reinitialize(self,idx,logical_size):
    if idx >= logical_size:
      self.skip_flag = True
    else:
      self.skip_flag = False
      self.datax = self.datax[idx]
      self.datay = self.datay[idx]

  def compute_functional_and_gradients(self):
    if self.skip_flag: return 0,flex.double(self.n)
    a = self.x
    residual = (self.datay - a[0]*self.datax*self.datax - a[1]*self.datax - a[2])
    f = 0.5 * residual * residual
    g = flex.double(self.n)
    dR_da = -self.datax*self.datax
    dR_db = -self.datax
    dR_dc = -1.
    g[0] = residual * dR_da
    g[1] = residual * dR_db
    g[2] = residual * dR_dc
    return f,g

def run_mpi():
  from libtbx.mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  #print ("hello from rank %d of %d"%(rank,size))

  W = simple_quadratic()
  if rank==0:
    W.run()
    print(list(W.a), "Single process final answer")
  else:
    pass
  comm.barrier()


  M = mpi_quadratic()
  M.reinitialize(idx=rank, logical_size=len(W.datax))
  minimizer = mpi_split_evaluator_run(target_evaluator=M,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-3,
        max_calls=1000)
      )
  if rank==0:
    print(list(M.x), "MPI final answer")
    try:
      from libtbx.test_utils import approx_equal
      assert approx_equal(M.x,W.a)
      assert approx_equal(M.x,[0.09601410216133123, -0.28424727078557327, -4.848332140888606])
      print ("OK")
    except Exception:
      print ("FAIL")

if __name__=="__main__":
  Usage = """
srun -n 32 -c 2 libtbx.python scitbx/lbfgs/tst_mpi_split_evaluator.py #small test case, 1 node
...only works when MPI is present, e.g., salloc -C haswell -N1 -q interactive -t 00:15:00
"""
  run_mpi()
