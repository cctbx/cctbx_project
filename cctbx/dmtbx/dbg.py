import sys, os, math
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx_boost import dmtbx

def ampl_phase_rad(f):
  amplidutes = shared.double()
  phases = shared.double()
  for i in xrange(f.size()):
    a, p = xutils.f_as_ampl_phase(f[i], deg=0)
    amplidutes.append(a)
    phases.append(p)
  return amplidutes, phases

def inplace_divide(array, arg):
  for i in xrange(array.size()):
    array[i] = array[i] / arg

def erase_small(miller_indices, data, cutoff):
  assert miller_indices.size() == data.size()
  for i in xrange(miller_indices.size()-1, -1, -1):
    if (data[i] < cutoff):
      miller_indices.erase(i)
      data.erase(i)

def compute_mean_weighted_phase_error(h, f, phi1, phi2, verbose=0):
  sum_w_phase_error = 0
  sum_w = 0
  for i in xrange(h.size()):
    phase_error = debug_utils.phase_error(
      phi1[i], phi2[i], deg=0) * 180/math.pi
    if (0 or verbose):
      print h[i], "%.3f %.3f %.3f" % (
        phi1[i], phi2[i], phase_error)
    sum_w_phase_error += f[i] * phase_error
    sum_w += f[i]
  return sum_w_phase_error / sum_w

def square_emap(xtal, e000, p1_miller_indices, miller_indices, e_values, phases):
  index_span = sftbx.index_span(p1_miller_indices)
  if (1):
    print "index_span:"
    print "  min:", index_span.min()
    print "  max:", index_span.max()
    print "  abs_range:", index_span.abs_range()
    print "  map_grid:", index_span.map_grid()
  grid_logical = list(index_span.map_grid())
  if (grid_logical[2] % 2): grid_logical[2] -= 1
  grid_logical = [n * 1 for n in grid_logical]
  print "grid_logical:", grid_logical
  rfft = fftbx.real_to_complex_3d(grid_logical)
  friedel_flag = 1
  old_e_complex = shared.polar_rad(e_values, phases)
  n_complex = rfft.Ncomplex()
  print "n_complex:", n_complex
  conjugate = 0
  map = sftbx.structure_factor_map(
    xtal.SgOps, friedel_flag, miller_indices,
    old_e_complex, n_complex, conjugate)
  map[0] = e000
  rfft.backward(map)
  rmap = shared.reinterpret_complex_as_real(map)
  shared.set_if_less_than(rmap, 0, 0)
  shared.square(rmap)
  rfft.forward(map)
  new_e_complex = sftbx.collect_structure_factors(
    friedel_flag, miller_indices, map, n_complex, conjugate)
  new_phases = shared.arg_rad(new_e_complex)
  if (0):
    for i in xrange(miller_indices.size()):
      print miller_indices[i], "%.2f %.2f" % (
        phases[i]*180/math.pi,
        new_phases[i]*180/math.pi)
  return new_phases

def test_triplet_invariants(sginfo, miller_indices_h, e_values, phases,
                            verbose):
  tprs = dmtbx.triplet_invariants(sginfo, miller_indices_h, e_values)
  utprs = tprs.unique_triplets()
  print "number_of_weighted_triplets:", \
       tprs.number_of_weighted_triplets(), \
       utprs.number_of_weighted_triplets()
  print "total_number_of_triplets:", \
       tprs.total_number_of_triplets(), \
       utprs.total_number_of_triplets()
  print "average_number_of_triplets_per_reflection: %.2f %.2f" % (
       tprs.average_number_of_triplets_per_reflection(),
       utprs.average_number_of_triplets_per_reflection())
  means = []
  for t in (tprs, utprs):
    mean_weighted_phase_error = []
    for ignore_weights in (0, 1):
      new_phases = t.refine_phases(
        miller_indices_h, e_values, phases, ignore_weights)
      mean_weighted_phase_error.append(
        compute_mean_weighted_phase_error(
          miller_indices_h, e_values, phases, new_phases, verbose))
    print "mean weighted phase error: %.2f %.2f delta: %.2f" % (
      mean_weighted_phase_error[0],
      mean_weighted_phase_error[1],
      mean_weighted_phase_error[1] - mean_weighted_phase_error[0])
    means.append(mean_weighted_phase_error)
  print "        delta phase error: %.2f %.2f" % (
    means[1][0] - means[0][0],
    means[1][1] - means[0][1])
  if (0 or verbose):
    tprs.dump_triplets(miller_indices_h)
  return tprs

def exercise(SgInfo,
             number_of_point_atoms = 10,
             d_min=1.,
             e_min=1.4,
             exercise_triplets=0,
             exercise_squaring=0,
             verbose=0):
  elements = ["const"] * number_of_point_atoms
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    no_random_u=1)
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  MillerIndices = xutils.build_miller_indices(xtal, friedel_flag=1,d_min=d_min)
  MillerIndices.H.append((0,0,0))
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal, abs_F=1)
  print "F000:", Fcalc.F[MillerIndices.H.size()-1]
  e_values = Fcalc.F
  inplace_divide(
    e_values, math.sqrt(xtal.SgOps.OrderZ() * number_of_point_atoms))
  e000 = e_values[MillerIndices.H.size()-1]
  print "E000:", e000
  MillerIndices.H.pop_back()
  e_values.pop_back()
  dmtbx.inplace_sort(MillerIndices.H, e_values, 1)
  s = shared.statistics(e_values)
  print "mean2:", s.mean2()
  print "number of structure factors:", e_values.size()
  erase_small(MillerIndices.H, e_values, e_min)
  print "number of structure factors:", e_values.size()
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal, abs_F=0)
  dummy, phases = ampl_phase_rad(Fcalc.F)
  Fcalc.F = e_values
  if (0 or verbose):
    debug_utils.print_structure_factors(Fcalc)
  p1_H = shared.Miller_Index()
  p1_e_values = shared.double()
  p1_phases = shared.double()
  sgtbx.expand_to_p1(
    xtal.SgOps,
    MillerIndices.H, e_values, phases,
    p1_H, p1_e_values, p1_phases,
    0, 0)
  print "number of structure factors p1:", p1_H.size()
  if (exercise_triplets):
    tprs_sg = test_triplet_invariants(
      xtal.SgInfo, MillerIndices.H, e_values, phases, verbose)
    tprs_p1 = test_triplet_invariants(
      sgtbx.SpaceGroup().Info(), p1_H, p1_e_values, p1_phases, verbose)
    js = shared.join_sets(MillerIndices.H, p1_H)
    assert js.pairs().size() == MillerIndices.H.size()
    sg_new_phases = tprs_sg.refine_phases(MillerIndices.H, e_values, phases)
    p1_new_phases = tprs_p1.refine_phases(p1_H, p1_e_values, p1_phases)
    if (0 or verbose):
      for i,j in js.pairs():
        h = str(MillerIndices.H[i])
        print "%-15s %.2f %.2f" % (
          h, phases[i]*180/math.pi, sg_new_phases[i]*180/math.pi)
        print "%-15s %.2f %.2f" % (
          "", p1_phases[j]*180/math.pi, p1_new_phases[j]*180/math.pi)
      print
  if (exercise_squaring):
    new_phases = square_emap(xtal, e000, p1_H, MillerIndices.H, e_values, phases)
    mwpe = compute_mean_weighted_phase_error(
      MillerIndices.H, e_values, phases, new_phases, verbose)
    print "squaring mean weighted phase error: %.2f" % (mwpe,)
  print

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "IncludeVeryHighSymmetry",
    "ShowSymbolOnly",
    "triplets",
    "squaring",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  symbols_to_stdout = 0
  if (len(sys.argv) > 1 + Flags.n):
    symbols = sys.argv[1:]
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
  for RawSgSymbol in symbols:
    if (RawSgSymbol.startswith("--")): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(RawSgSymbol)
    SgInfo = sgtbx.SpaceGroup(SgSymbols).Info()
    LookupSymbol = SgInfo.BuildLookupSymbol()
    sys.stdout.flush()
    print >> sys.stderr, LookupSymbol
    sys.stderr.flush()
    if (symbols_to_stdout):
      print LookupSymbol
      sys.stdout.flush()
    if (SgInfo.SgOps().OrderZ() > 48 and not Flags.IncludeVeryHighSymmetry):
      print "High symmetry space group skipped."
      print
      sys.stdout.flush()
      continue
    if (Flags.ShowSymbolOnly):
      print "Space group:", LookupSymbol
    else:
      exercise(
        SgInfo,
        exercise_triplets=Flags.triplets,
        exercise_squaring=Flags.squaring)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
