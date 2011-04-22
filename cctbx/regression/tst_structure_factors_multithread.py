from cctbx import sgtbx, xray, crystal
from cctbx.development import random_structure
try:
  from cctbx.xray.ext import structure_factors_raw_multithreaded_direct
except ImportError:
  structure_factors_raw_multithreaded_direct = None
from cctbx import math_module
import omptbx
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import wall_clock_time, show_times_at_exit
from libtbx.introspection import number_of_processors
from libtbx import group_args
import sys, os

def show_times_vs_cpu(times, header):
  cols = [ "# CPU" + " "*3,
           "run-time" + " "*3,
           "speed-up" ]
  print header
  print ''.join(cols)
  fmt = "%%-%ii%%-%i.2fx" % tuple([ len(c) for c in cols[:2] ])
  t0 = times[0][1]
  for i,t in times:
    if (t == 0): s = "Inf"
    else: s = "%.2f" % (t0/t)
    print fmt % (i,t), s

def exercise_openmp(space_group_info,
                    elements,
                    anomalous_flag=False,
                    verbose=0):
  if (omptbx.have_omp_h): print "omptbx.have_omp_h"
  if (omptbx.have_stubs_h): print "omptbx.have_stubs_h"
  xs = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    volume_per_atom=18.6,
    min_distance=1.2)
  if 0:
    cos_sin_table = math_module.cos_sin_table(2**10)
  if 1:
    cos_sin_table = False
  timer = wall_clock_time()
  times = []
  omptbx.env.num_threads = 1
  if (verbose): d_min = 0.5
  else:         d_min = 2.0
  single_threaded_calc = xs.structure_factors(d_min=d_min,
                                              algorithm="direct",
                                              cos_sin_table=cos_sin_table)
  times.append((1, timer.elapsed()))
  single_threaded_f = single_threaded_calc.f_calc()
  for n in xrange(2, max(2, number_of_processors()) + 1):
    omptbx.env.num_threads = n
    timer = wall_clock_time()
    multi_threaded_calc = xs.structure_factors(d_min=d_min,
                                               algorithm='direct',
                                               cos_sin_table=cos_sin_table)
    times.append((n, timer.elapsed()))
    multi_threaded_f = multi_threaded_calc.f_calc()
    assert approx_equal(single_threaded_f.data(), multi_threaded_f.data())
    if (not omptbx.have_omp_h):
      break
  if verbose:
    show_times_vs_cpu(times, header="OpenMP")

def exercise_openmp_resilience_to_adptbx_exception(space_group_info,
                                                   verbose=0):
  if not omptbx.have_omp_h:
    print ('OpenMP not available:'
           ' test of resilience to adptbx exception skipped.')
    return
  import signal
  def handler(signum, frame):
    print ('Error: '
           ' exception thrown in OpenMP parallel region was not caught.\n'
           'The program aborted as a result.')
    signal.signal(signal.SIGABRT, signal.SIG_DFL)
    os.abort()

  signal.signal(signal.SIGABRT, handler)
  isinstance(space_group_info, sgtbx.space_group_info)
  cs = space_group_info.any_compatible_crystal_symmetry(volume=1000)
  xs = xray.structure(crystal.special_position_settings(cs))
  xs.add_scatterer(xray.scatterer('C*', site=(0,0,0), u=-100))
  xs.add_scatterer(xray.scatterer('C1', site=(0.1, 0.2, 0.3), u=2))
  xs.add_scatterer(xray.scatterer('C2', site=(0, 0, 0.5), u=1.5))
  if verbose:
    print 'OpenMP using %s processors' % omptbx.env.num_procs
  try:
    xs.structure_factors(d_min=2,
                         algorithm='direct',
                         cos_sin_table=False)
  except RuntimeError, e:
    assert str(e).find('cctbx::adptbx::debye_waller_factor_exp:'
                       ' arg_limit exceeded') != -1
  else:
    raise Exception_expected


def exercise_raw(space_group_info,
                 elements,
                 anomalous_flag=False,
                 verbose=0):
  if structure_factors_raw_multithreaded_direct is None:
    print "Skipping raw multithreaded structure factor computation tests"
    return
  xs = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    volume_per_atom=18.6,
    min_distance=1.2)
  if 0:
    cos_sin_table = math_module.cos_sin_table(2**10)
  if 1:
    cos_sin_table = False
  timer = wall_clock_time()
  times = []
  if (verbose): d_min = 0.5
  else:         d_min = 2.0
  omptbx.env.num_threads = 1
  single_threaded_calc = xs.structure_factors(d_min=d_min,
                                              algorithm="direct",
                                              cos_sin_table=cos_sin_table)
  times.append((1, timer.elapsed()))
  single_threaded_f = single_threaded_calc.f_calc()
  args = group_args(unit_cell=xs.unit_cell(),
                    space_group=xs.space_group(),
                    miller_indices=single_threaded_f.indices(),
                    scatterers=xs.scatterers(),
                    scattering_type_registry=xs.scattering_type_registry())
  for i in xrange(1, max(2, number_of_processors()) + 1):
    arg1s = group_args(n_threads=i, **args.__dict__)
    if cos_sin_table:
      arg1s = group_args(cos_sin_table=cos_sin_table, **args.__dict__)
    timer = wall_clock_time()
    multithreaded_calc = structure_factors_raw_multithreaded_direct(
      **arg1s.__dict__)
    times.append((i, timer.elapsed()))
    multithreaded_f = multithreaded_calc.f_calc()
    assert approx_equal(single_threaded_f.data(), multithreaded_f)
  if verbose:
    show_times_vs_cpu(
      times,
      header=("Boost.Thread "
              "(first line: ext.structure_factors_direct running on 1 thread)"
              ))

def run(args):
  show_times_at_exit()
  verbose = '--verbose' in args
  sgi = sgtbx.space_group_info("P21/n")
  exercise_openmp_resilience_to_adptbx_exception(sgi, verbose=verbose)
  elements = ['O']*15 + ['N']*9 + ['C']*100
  exercise_openmp(sgi, elements, verbose=verbose)
  exercise_raw(sgi, elements, verbose=verbose)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
