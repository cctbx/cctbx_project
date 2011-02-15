import random
random.seed(0)
from scitbx.array_family import flex
flex.set_random_seed(0)
from math import pow,pi
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.xray.structure_factors.from_scatterers_direct_parallel import direct_summation_simple
from cctbx.xray.structure_factors.from_scatterers_direct_parallel import direct_summation_cuda_platform
from libtbx.utils import wall_clock_time, show_times_at_exit
from libtbx.test_utils import approx_equal

def show_times_vs_complexity(times, header):
  cols = [ "# HKL" + " "*3,
           "cpu-time" + " "*3,
           "simple-tm" + " "*3,
           "fft-time" + " "*3,
           "gpu-time" + " "*3,
           "d-min(angstrom)" ]
  print header
  print ''.join(cols)
  fmt = "%%-%ii%%-%i.2f%%-%i.2f%%-%i.2f%%-%i.2f%%5.2f" % tuple(
   [ len(c) for c in cols[:5] ])
  for i,t,d,g,s,f in times:
    print fmt % (i,t,s,f,g,d)

def show_diagnostics(xs):
  #help(xs.scatterers())
  print list(xs.scatterers().extract_labels())
  print list(xs.scatterers().extract_occupancies())
  print list(xs.scatterers().extract_scattering_types())
  print list(xs.scatterers().extract_sites())
  print
  print list(xs.scatterers().extract_u_cart(xs.unit_cell()))
  print
  print list(xs.scatterers().extract_u_iso())
  print
  print list(xs.scatterers().extract_u_star())

def exercise_direct(space_group_info,
                    elements,
                    anomalous_flag=False,
                    verbose=0):
  xs = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    volume_per_atom=18.6,
    min_distance=1.2)
  #show_diagnostics(xs)

  reciprocal_volume = xs.unit_cell().reciprocal().volume()
  times = []

  cuda_platform = direct_summation_cuda_platform()

  for x in xrange(1,7):
    print "There are %d scatterers"%len(elements)
    number_of_reflections = pow(10.,x)
    Volume = number_of_reflections *  reciprocal_volume * 2. #2 P1 asymmetric units
    Volume *= space_group_info.group().order_z() # take space group into acct.
    recip_radius = pow(3.*Volume/(4.*pi),1./3.)
    d_min = 1./recip_radius

    if 0:
      cos_sin_table = math_module.cos_sin_table(2**10)
    if 1:
      cos_sin_table = False

    timer = wall_clock_time()
    cpu_direct = xs.structure_factors(d_min=d_min,algorithm="direct",
                                            cos_sin_table=cos_sin_table)
    cpu_time = timer.elapsed()
    cpu_direct_f = cpu_direct.f_calc()

    timer = wall_clock_time()
    gpu_direct = xs.structure_factors(d_min=d_min,algorithm=cuda_platform)
    gpu_time = timer.elapsed()
    gpu_direct_f = gpu_direct.f_calc()

    timer = wall_clock_time()
    fft_algorithm = xs.structure_factors(d_min=d_min,algorithm="fft")
    fft_time = timer.elapsed()
    fft_f = fft_algorithm.f_calc()

    timer = wall_clock_time()
    simple_direct = xs.structure_factors(d_min=d_min,algorithm=direct_summation_simple())
    simple_time = timer.elapsed()
    simple_direct_f = simple_direct.f_calc()

    times.append((number_of_reflections,cpu_time,d_min,gpu_time,simple_time,fft_time))
    assert approx_equal(cpu_direct_f.data(), gpu_direct_f.data(), eps=1e-4)
    # doesn't assert correctly assert approx_equal(cpu_direct_f.data(), fft_f.data())
    assert approx_equal(cpu_direct_f.data(), simple_direct_f.data(), eps=1e-6)

  show_times_vs_complexity(times, header="run time vs. # reflections")

def run_scattering_type_tests():
  for C in xrange(35):
    for N in xrange(35):
      for total in xrange(93,99):
        elements = ['C']*C + ['N']*N + ['O']*(total-N-C)
        print "".join(elements)
        exercise_direct(sgtbx.space_group_info("P1"),elements)

def run(args,multiplier):
  show_times_at_exit()
  verbose = '--verbose' in args
  #count from 1hmg.pdb, chain A: C, 1583; N, 445; O, 495, S, 13
  elements = ['O']*19 + ['N']*18 + ['C']*62 + ['S']*1
  allelements = elements*multiplier

  if 0:
    for sn in xrange(1,231):
      try:
        sgi = sgtbx.space_group_info(sn)
        print "Space group",sgi,"number",sn
        exercise_direct(sgi, allelements, verbose=verbose)
      except Exception, e:
        print e
    return

  if 0:
    for symbol in ["P1","P3","P41","P212121","I41","F432"]:
      sgi = sgtbx.space_group_info(symbol)
      print "Space group",sgi
      exercise_direct(sgi, allelements, verbose=verbose)

  if 1:
    sgi = sgtbx.space_group_info("P1")
    print "Space group",sgi
    exercise_direct(sgi, allelements, verbose=verbose)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:],10)
  run(sys.argv[1:],20)
  run(sys.argv[1:],40)
  run(sys.argv[1:],80)
  run(sys.argv[1:],160)
