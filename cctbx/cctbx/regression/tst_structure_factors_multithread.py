from cctbx import sgtbx
from cctbx.development import random_structure
try:
  from cctbx.xray.ext import structure_factors_raw_multithreaded_direct
except ImportError:
  structure_factors_raw_multithreaded_direct = None
from scitbx import openmp
from cctbx import math_module
from libtbx.test_utils import approx_equal
from libtbx.utils import wall_clock_time, show_times_at_exit
from libtbx.introspection import number_of_processors
from libtbx import group_args

def show_times_vs_cpu(times, header):
  cols = [ "# CPU" + " "*3,
           "run-time" + " "*3,
           "speed-up" ]
  print header
  print ''.join(cols)
  fmt = "%%-%ii%%-%i.2fx" % tuple([ len(c) for c in cols[:2] ])
  for i,t in enumerate(times):
    if (times[i] == 0): s = "Inf"
    else: s = "%.2f" % (times[0]/times[i])
    print fmt % (i+1, times[i]), s

def exercise_openmp(space_group_info,
                    elements,
                    anomalous_flag=False,
                    verbose=0):
  if not openmp.available:
    print "Skipping OpenMP structure factor computation tests"
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
  openmp.environment.num_threads = 1
  if (verbose): d_min = 0.5
  else:         d_min = 2.0
  single_threaded_calc = xs.structure_factors(d_min=d_min,
                                              algorithm="direct",
                                              cos_sin_table=cos_sin_table)
  times.append(timer.elapsed())
  single_threaded_f = single_threaded_calc.f_calc()
  for n in xrange(2, max(2, number_of_processors()) + 1):
    openmp.environment.num_threads = n
    timer = wall_clock_time()
    multi_threaded_calc = xs.structure_factors(d_min=d_min,
                                               algorithm='direct',
                                               cos_sin_table=cos_sin_table)
    times.append(timer.elapsed())
    multi_threaded_f = multi_threaded_calc.f_calc()
    assert approx_equal(single_threaded_f.data(), multi_threaded_f.data())
  if verbose:
    show_times_vs_cpu(times, header="OpenMP")

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
  if openmp.available: openmp.environment.num_threads = 1
  if (verbose): d_min = 0.5
  else:         d_min = 2.0
  single_threaded_calc = xs.structure_factors(d_min=d_min,
                                              algorithm="direct",
                                              cos_sin_table=cos_sin_table)
  times.append(timer.elapsed())
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
    times.append(timer.elapsed())
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
  elements = ['O']*15 + ['N']*9 + ['C']*100
  exercise_openmp(sgi, elements, verbose=verbose)
  exercise_raw(sgi, elements, verbose=verbose)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
