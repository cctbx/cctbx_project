from __future__ import absolute_import, division, print_function
from scitbx import lbfgsb
from scitbx.array_family import flex
import scitbx.math
from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
import sys
from six.moves import range

def exercise_minimizer_interface():
  n = 25
  m = 5
  l = flex.double(n, -1)
  u = flex.double(n, 1)
  nbd = flex.int(n)
  factr=1.0e+7
  pgtol=1.0e-5
  iprint = -1
  for enable_stp_init in [False, True]:
    minimizer = lbfgsb.ext.minimizer(
      n, m, l, u, nbd, enable_stp_init, factr, pgtol, iprint)
    assert minimizer.n() == n
    assert minimizer.m() == m
    assert minimizer.l().id() == l.id()
    assert minimizer.u().id() == u.id()
    assert minimizer.nbd().id() == nbd.id()
    assert minimizer.enable_stp_init() == enable_stp_init
    assert eps_eq(minimizer.factr(), factr)
    assert eps_eq(minimizer.pgtol(), pgtol)
    assert eps_eq(minimizer.iprint(), iprint)
    assert not minimizer.requests_f_and_g()
    assert not minimizer.is_terminated()
    assert minimizer.task() == "START"
    x = flex.double(n, 0)
    f = 1
    g = flex.double(n, -1)
    assert minimizer.process(x, f, g, False)
    assert minimizer.task() == "FG_START"
    assert minimizer.f() == 0
    if (not enable_stp_init):
      try:
        minimizer.requests_stp_init()
      except RuntimeError as e:
        assert str(e).endswith(": SCITBX_ASSERT(enable_stp_init()) failure.")
      else: raise Exception_expected
    else:
      assert not minimizer.process(x, f, g)
      assert minimizer.requests_stp_init()
      assert approx_equal(minimizer.relative_step_length_line_search(), 0.2)
      assert approx_equal(minimizer.current_search_direction(), [1]*n)
      minimizer.set_relative_step_length_line_search(value=0.3)
      assert approx_equal(minimizer.relative_step_length_line_search(), 0.3)
    assert minimizer.process(x, f, g)
    assert minimizer.f() == 1
    assert minimizer.task() == "FG_LNSRCH"
    assert not minimizer.is_terminated()
    if (not enable_stp_init):
      assert approx_equal(x, [0.2]*n)
    else:
      assert approx_equal(x, [0.3]*n)
    minimizer.request_stop()
    assert not minimizer.process(x, f, g)
    assert minimizer.task() == "STOP: NO RESTORE"
    minimizer.request_stop_with_restore()
    assert minimizer.task() == "STOP: CPU"
    minimizer.request_restart()
    assert not minimizer.requests_f_and_g()
    assert not minimizer.is_terminated()
    assert minimizer.task() == "START"
    minimizer = lbfgsb.minimizer(n=n)
    assert minimizer.l().size() == n
    assert minimizer.u().size() == n
    assert minimizer.nbd().size() == n
    assert minimizer.nbd().all_eq(0)

def driver1(use_fortran_library=False):
  n = 25
  nbd = flex.int(n)
  x = flex.double(n)
  l = flex.double(n)
  u = flex.double(n)
  g = flex.double(n)
  if ("--Verbose" in sys.argv[1:]):
    iprint = 1000
  else:
    iprint = -1
  for i in range(0,n,2):
    nbd[i] = 2
    l[i] = 1.0e0
    u[i] = 1.0e2
  for i in range(1,n,2):
    nbd[i] = 2
    l[i] = -1.0e2
    u[i] = 1.0e2
  for i in range(n):
    x[i] = 3.0e0
  minimizer = lbfgsb.minimizer(
    n=n,
    m=5,
    l=l,
    u=u,
    nbd=nbd,
    factr=1.0e+7,
    pgtol=1.0e-5,
    iprint=iprint)
  f = 0
  while True:
    if (minimizer.process(x, f, g, use_fortran_library)):
      f=.25e0*(x[0]-1.e0)**2
      for i in range(1,n):
        f=f+(x[i]-x[i-1]**2)**2
      f=4.e0*f
      t1=x[1]-x[0]**2
      g[0]=2.e0*(x[0]-1.e0)-1.6e1*x[0]*t1
      for i in range(1,n-1):
        t2=t1
        t1=x[i+1]-x[i]**2
        g[i]=8.e0*t2-1.6e1*x[i]*t1
      g[n-1]=8.e0*t1
    elif (minimizer.is_terminated()):
      break
  assert minimizer.task() == "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  assert minimizer.f_list().size() == minimizer.n_iteration() + 1
  assert minimizer.f_list()[-1] == minimizer.f()
  assert minimizer.f_list()[-2] == minimizer.f_previous_iteration()
  assert not minimizer.initial_x_replaced_by_projection()
  assert minimizer.is_constrained()
  assert minimizer.is_fully_constrained()
  if ("--soft" not in sys.argv[1:]):
    assert minimizer.n_fg_evaluations_total() == 27
    assert minimizer.n_fg_evaluations_iter() == 1
    assert minimizer.n_intervals_explored_cauchy_search_total() == 48
    assert minimizer.n_intervals_explored_cauchy_search_iter() == 1
    assert minimizer.n_skipped_bfgs_updates_total() == 0
    assert minimizer.n_bfgs_updates_total() == 22
    assert minimizer.subspace_argmin_is_within_box() == 1
    assert minimizer.n_free_variables() == 25
    assert minimizer.n_active_constraints() == 0
    assert minimizer.n_variables_leaving_active_set() == 0
    assert minimizer.n_variables_entering_active_set() == 0
    assert eps_eq(minimizer.theta_bfgs_matrix_current(), 23.2674689856)
    assert minimizer.f_previous_iteration() >= 0
    assert minimizer.floating_point_epsilon() \
        == scitbx.math.floating_point_epsilon_double_get()
    assert eps_eq(minimizer.factr_times_floating_point_epsilon(),
                  minimizer.factr() * minimizer.floating_point_epsilon())
    assert eps_eq(minimizer.two_norm_line_search_direction_vector(),
                  1.56259327735e-05)
    assert eps_eq(minimizer.two_norm_line_search_direction_vector_sq(),
      minimizer.two_norm_line_search_direction_vector()**2)
    assert minimizer.accumulated_time_cauchy_search() > -0.02
    assert minimizer.accumulated_time_subspace_minimization() > -0.02
    assert minimizer.accumulated_time_line_search() > -0.02
    assert eps_eq(minimizer.slope_line_search_function_current(),
                  9.31496613169e-10)
    assert eps_eq(minimizer.slope_line_search_function_start(),
                  -3.6762390377e-09)
    assert eps_eq(minimizer.maximum_relative_step_length(), 1.37484166517)
    assert eps_eq(minimizer.relative_step_length_line_search(), 1.0)
    assert eps_eq(minimizer.infinity_norm_projected_gradient(),
                  0.000143369780806)

def run():
  exercise_minimizer_interface()
  driver1(use_fortran_library=("--fortran" in sys.argv[1:]))
  print("OK")

if (__name__ == "__main__"):
  run()
