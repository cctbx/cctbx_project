from cctbx import sgtbx
from cctbx.development import random_structure
try:
  from cctbx.xray.ext import structure_factors_multithreaded_direct
except ImportError:
  structure_factors_multithreaded_direct = None
from cctbx import math_module
from libtbx.test_utils import approx_equal
from libtbx.utils import wall_clock_time, show_times_at_exit
from libtbx.introspection import number_of_processors
from libtbx import group_args

def exercise_structure_factors(space_group_info,
                               elements,
                               d_min=0.5,
                               anomalous_flag=False,
                               verbose=0):
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
  for i in xrange(1, number_of_processors() + 1):
    arg1s = group_args(n_threads=i, **args.__dict__)
    if cos_sin_table:
      arg1s = group_args(cos_sin_table=cos_sin_table, **args.__dict__)
    timer = wall_clock_time()
    multithreaded_calc = structure_factors_multithreaded_direct(
      **arg1s.__dict__)
    times.append(timer.elapsed())
    multithreaded_f = multithreaded_calc.f_calc()
    assert approx_equal(single_threaded_f.data(), multithreaded_f)
  if verbose:
    cols = [ "# CPU" + " "*3,
             "run-time" + " "*3,
             "speed-up" ]
    print ("First line is the non-parallel ext.structure_factors_direct"
           " and the speed-up is with respect to that")
    print ''.join(cols)
    fmt = "%%-%ii%%-%i.2fx %%-%i.2f" % tuple([ len(c) for c in cols ])
    for i,t in enumerate(times):
      print fmt % (i, times[i], times[0]/times[i])

def run(args):
  if structure_factors_multithreaded_direct is None:
    print "Skipping multithreaded structure factor computation tests"
    return
  verbose = '--verbose' in args
  show_times_at_exit()
  sgi = sgtbx.space_group_info("P21/n")
  elements = ['O']*15 + ['N']*9 + ['C']*100
  exercise_structure_factors(sgi, elements, verbose=verbose)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
