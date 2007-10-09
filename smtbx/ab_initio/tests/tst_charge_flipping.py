from __future__ import division

import os.path
import sys
import random

from cctbx import sgtbx
from cctbx import miller
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import itertbx

from smtbx.ab_initio import charge_flipping

def random_structure_factors(space_group_info, elements,
                             anomalous_flag,
                             d_min, grid_resolution_factor
                             ):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=False,
    use_u_aniso=False,
  )
  for s in structure.scatterers():
    assert s.flags.use_u_iso() and s.u_iso == 0
    assert not s.flags.use_u_aniso()

  f_indices = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=anomalous_flag,
    d_min=d_min)
  f_target = f_indices.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc()

  f_000 = flex.miller_index(1)
  f_000[0] = (0,0,0)
  f_000 = miller.set(structure, f_000).structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc().data()[0]
  return f_target, f_000



def exercise_fixed_delta(space_group_info, elements,
                         anomalous_flag,
                         d_min=0.8, grid_resolution_factor=1./2,
                         debug=False,
                         verbose=False
                         ):
  f_target, f_000 = random_structure_factors(space_group_info, elements,
                                             anomalous_flag,
                                             d_min, grid_resolution_factor)

  f_target_in_p1 = f_target.expand_to_p1()\
                           .as_non_anomalous_array()\
                           .merge_equivalents().array()
  rho_target_map = f_target_in_p1.fft_map(f_000=f_000)
  rho_target = rho_target_map.real_map_unpadded()
  max_rho_target = flex.max(rho_target)
  delta = 0.01 * max_rho_target # got that by running the debug branch just
                                # below
  if debug:
    # Let the user figure out graphically the value of delta using
    # a graph like the ones on figure 1 of ref [1] given in module
    # charge_flipping.
    # This makes a Mathematica .m file which creates the plot.
    rho_target /= max_rho_target
    voxels_sorting = flex.sort_permutation(rho_target.as_1d())
    sorted_voxels = rho_target.as_1d().select(voxels_sorting)
    mf = open(os.path.expanduser('~/Desktop/foo.m'), 'w')
    print >> mf, 'rho = { %s 1};' % ",\n".join(
      [ "%.3f" % v for v in sorted_voxels])
    print >> mf, 'ListPlot[rho,'
    print >> mf, '         PlotRange->{Automatic,{-0.05,0.5}},'
    print >> mf, '         GridLines->{None,{-0.02,0.04,0.001}}]'
    mf.close()
    print "Enter delta:",
    delta = float(sys.stdin.readline().strip()) * max_rho_target

  if verbose:
    print "c_tot / c_flip of target = %.3f" % (
      rho_target_map.c_tot() / rho_target_map.c_flip(delta) )

  f_obs = abs(f_target)

  basic_iter = charge_flipping.basic_iterator
  advanced_iter = charge_flipping.weak_reflection_improved_iterator
  final_r1 = []
  for flipped in (basic_iter(f_obs, delta),
                  advanced_iter(f_obs, delta),
                  ):
    if verbose:
      print "Starting %s" % type(flipped).__name__
      print "\tdelta/max(rho_target) = %.4f" % (flipped.delta/max_rho_target)
      print "\tR_1 = %.3f" % f_target_in_p1.r1_factor(
        f_target_in_p1.structure_factors_from_map(flipped.rho_map.real_map()))
    r1s = flex.double()
    for i,state in enumerate(itertbx.islice(flipped, 10)):
      isinstance(state, charge_flipping.density_modification_iterator)
      r1 = f_target_in_p1.r1_factor(state.g)
      r1s.append(r1)
      if debug:
        rho = state.rho_map.real_map()
        pos = (rho > 0).count(True) / rho.size()
        print "\t%i: delta=%.5f" % (i,state.delta)
        indent = " "*(len("%i" % i) + 1)
        print indent, "\tR1=%.5f" % r1
        print indent, "\tfraction of positive pixels=%.3f" % pos
        print indent, "\ttotal charge=%.3f" % state.g_000

    # r1 decreasing function of the iteration number
    assert list(flex.sort_permutation(r1s, reverse=True))\
           == list(xrange(r1s.size()))

    assert r1s[-1] < 0.12
    diffs = r1s[:-1] - r1s[1:]
    assert diffs[-3:].all_lt(0.01)

    final_r1.append(r1s[-1])
    if verbose:
      print "Finished"
      print "\tR_1=%.3f" % r1s[-1]
      print

  # The both of the basic and advanced charge flipping are compatible
  assert (approx_equal(final_r1[0], final_r1[1], 0.01, out=None))


def exercise(flags, space_group_info):
  for i in xrange(2):
    n_C = random.randint(10,20)
    n_O = random.randint(1,5)
    n_N = random.randint(1,5)
    if flags.Verbose:
      print "C%i O%i N%i" % (n_C, n_O, n_N)

    exercise_fixed_delta(
      space_group_info=space_group_info,
      elements=["C"]*n_C + ["O"]*n_O + ["N"]*n_N,
      anomalous_flag=False,
      debug=flags.Debug,
      verbose=flags.Verbose
    )

    if flags.Verbose: print

def run():
  import sys
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], exercise)
  print 'OK'

if __name__ == '__main__':
  import sys
  if '--profile' in sys.argv:
    import profile
    import pstats
    sys.argv.remove('--profile')
    profile.run('run()', 'charge_flipping.prof')
    p = pstats.Stats('charge_flipping.prof')
    p.strip_dirs().sort_stats('time').print_stats(10)
  else:
    run()
