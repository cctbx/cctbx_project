import sys, os, math
from cctbx_boost.arraytbx import flex
from cctbx_boost.arraytbx import flex_utils
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx_boost import dmtbx
from cctbx import euclidean_model_matching as emma

def ampl_phase(f, deg=0):
  amplidutes = flex.double()
  phases = flex.double()
  for i in xrange(f.size()):
    a, p = xutils.f_as_ampl_phase(f[i], deg)
    amplidutes.append(a)
    phases.append(p)
  return amplidutes, phases

def show_ampl_phases(miller_indices, ampl, phases, deg=0):
  for i in xrange(miller_indices.size()):
    phi = phases[i]
    if (not deg): phi *= 180/math.pi
    print miller_indices[i], "%.2f %.2f" % (ampl[i], phi)

def inplace_divide(array, arg):
  for i in xrange(array.size()):
    array[i] = array[i] / arg

def erase_small(miller_indices, data, cutoff):
  assert miller_indices.size() == data.size()
  for i in xrange(miller_indices.size()-1, -1, -1):
    if (data[i] < cutoff):
      miller_indices.erase(i)
      data.erase(i)

def normalize_quasi_normalized(sgops, miller_indices, data):
  epsilons = sgops.epsilon(miller_indices)
  for i in xrange(miller_indices.size()):
    data[i] /= epsilons[i]

def compute_mean_weighted_phase_error(tprs, sgops, h, f, phi1, phi2, verbose=0):
  sum_w_phase_error = 0
  sum_w = 0
  for i in xrange(h.size()):
    if (tprs != None and tprs.n_relations(i) == 0): continue
    phase_error = debug_utils.phase_error(
      phi1[i], phi2[i], deg=0) * 180/math.pi
    if (0 or verbose):
      print h[i], "%.3f %.3f %.3f" % (
        phi1[i], phi2[i], phase_error)
    m = sgops.multiplicity(h[i], 1)
    sum_w_phase_error += m * f[i] * phase_error
    sum_w += m * f[i]
  return sum_w_phase_error / sum_w

def raise_emap(xtal, e000, p1_miller_indices,
               miller_indices, e_values, phases,
               use_e000, zero_out_negative,
               exponent=2,
               verbose=0):
  index_span = miller.index_span(p1_miller_indices)
  # XXX index_span = xtal.SgOps.get_index_span(miller_indices)
  if (0 or verbose):
    print "index_span:"
    print "  min:", index_span.min()
    print "  max:", index_span.max()
    print "  abs_range:", index_span.abs_range()
    print "  map_grid:", index_span.map_grid()
  grid_logical_min = list(index_span.map_grid())
  if (grid_logical_min[2] % 2): grid_logical_min[2] -= 1
  if (0 or verbose): print "grid_logical_min:", grid_logical_min
  grid_logical_n = [n * 2 + 0 for n in grid_logical_min]
  if (0 or verbose): print "grid_logical_n:", grid_logical_n
  if (0):
    mandatory_gridding_factors = xtal.SgOps.refine_gridding()
  else:
    mandatory_gridding_factors = (1,1,1)
  grid_logical = fftbx.adjust_gridding_triple(
    grid_logical_n, 1, mandatory_gridding_factors)
  if (0 or verbose): print "grid_logical:", grid_logical
  rfft = fftbx.real_to_complex_3d(grid_logical)
  friedel_flag = 1
  old_e_complex = flex.polar(e_values, phases)
  n_complex = rfft.Ncomplex()
  if (0 or verbose): print "n_complex:", n_complex
  conjugate = 0
  map = sftbx.structure_factor_map(
    xtal.SgOps, friedel_flag, miller_indices,
    old_e_complex, n_complex, conjugate)
  if (use_e000):
    map[0] = e000
  rmap = rfft.backward(map)
  if (zero_out_negative):
    flex.set_if_less_than(rmap, 0, 0)
  flex_utils.in_place_pow(rmap, exponent)
  sf_map = rfft.forward(map)
  new_e_complex = sftbx.collect_structure_factors(
    friedel_flag, miller_indices, sf_map, conjugate)
  new_phases = flex.arg(new_e_complex)
  if (0 or verbose):
    for i in xrange(miller_indices.size()):
      print miller_indices[i], "%.2f %.2f" % (
        phases[i]*180/math.pi,
        new_phases[i]*180/math.pi)
  return new_phases

def test_triplet_invariants(sginfo, miller_indices, e_values, phases,
                            other_than_sigma_2,
                            verbose):
  tprs = dmtbx.triplet_invariants(
    sginfo, miller_indices, 1, other_than_sigma_2)
  print "number_of_weighted_triplets:", \
       tprs.number_of_weighted_triplets()
  print "total_number_of_triplets:", \
       tprs.total_number_of_triplets()
  print "average_number_of_triplets_per_reflection: %.2f" % (
       tprs.average_number_of_triplets_per_reflection(),)
  new_phases = tprs.apply_tangent_formula(e_values, phases)
  mean_weighted_phase_error = compute_mean_weighted_phase_error(
    tprs,
    sginfo.SgOps(), miller_indices, e_values, phases, new_phases,
    verbose)
  print "mean weighted phase error: %.2f" % (mean_weighted_phase_error,)
  if (0 or verbose):
    tprs.dump_triplets(miller_indices)
  if (0 or verbose):
    tprs.weights_and_epsilon(sginfo, miller_indices)
  return tprs

class simulated_data:

  def __init__(self, SgInfo,
               number_of_point_atoms = 10,
               d_min=1.,
               e_min=1.8,
               verbose=0):
    elements = ["const"] * number_of_point_atoms
    print "random.getstate():", debug_utils.random.getstate()
    if (0):
      debug_utils.random.setstate((1, (14401, 20036, 16503), None))
    xtal = debug_utils.random_structure(
      SgInfo, elements,
      volume_per_atom=50.,
      min_distance=1.5,
      general_positions_only=0,
      no_random_u=1)
    print "Unit cell:", xtal.UnitCell
    print "Space group:", xtal.SgInfo.BuildLookupSymbol()
    debug_utils.print_sites(xtal)
    miller_set = xutils.build_miller_set(
      xtal, friedel_flag=1,d_min=d_min)
    if (0):
      miller_set.H = flex.miller_Index()
      miller_set.H.append((0,2,0))
      miller_set.H.append((1,1,15))
    miller_set.H.append((0,0,0))
    Fcalc = xutils.calculate_structure_factors_direct(
      miller_set, xtal, abs_F=1)
    print "F000:", Fcalc.F[miller_set.H.size()-1]
    e_values = Fcalc.F
    inplace_divide(
      e_values, math.sqrt(xtal.SgOps.OrderZ() * number_of_point_atoms))
    e000 = e_values[miller_set.H.size()-1]
    print "E000:", e000
    miller_set.H.pop_back()
    e_values.pop_back()
    dmtbx.inplace_sort(miller_set.H, e_values, 1)
    s = flex_utils.statistics(e_values)
    print "mean2:", s.mean2()
    print "number of structure factors:", e_values.size()
    erase_small(miller_set.H, e_values, e_min)
    print "number of structure factors:", e_values.size()
    if (0): # XXX
      normalize_quasi_normalized(xtal.SgOps, miller_set.H, e_values)
    Fcalc = xutils.calculate_structure_factors_direct(
      miller_set, xtal, abs_F=0)
    dummy, phases = ampl_phase(Fcalc.F) # XXX use flex.arg()
    Fcalc.F = e_values
    if (0 or verbose):
      debug_utils.print_structure_factors(Fcalc)
    if (0 or verbose):
      show_ampl_phases(Fcalc.H, e_values, phases)
    self.xtal = xtal
    self.miller_set = miller_set
    self.e_values = e_values
    self.phases = phases
    self.e000 = e000

def exercise(SgInfo,
             number_of_point_atoms = 10,
             d_min=1.,
             e_min=1.8,
             exercise_triplets=0,
             other_than_sigma_2=0,
             exercise_squaring=0,
             use_e000=0,
             zero_out_negative=0,
             verbose=0):
  sim = simulated_data(SgInfo, number_of_point_atoms, d_min, e_min, verbose)
  p1_H = flex.miller_Index()
  p1_e_values = flex.double()
  p1_phases = flex.double()
  miller.expand_to_p1(
    sim.xtal.SgOps, 1,
    sim.miller_set.H, sim.e_values, sim.phases,
    p1_H, p1_e_values, p1_phases)
  print "number of structure factors p1:", p1_H.size()
  if (0 or verbose):
    show_ampl_phases(p1_H, p1_e_values, p1_phases)
  tprs_sg = None
  if (exercise_triplets):
    tprs_sg = test_triplet_invariants(
      sim.xtal.SgInfo, sim.miller_set.H, sim.e_values, sim.phases,
      other_than_sigma_2, verbose)
    if (other_than_sigma_2):
      tprs_p1 = test_triplet_invariants(
        sgtbx.SpaceGroup().Info(), p1_H, p1_e_values, p1_phases, 1, verbose)
      sg_new_phases = tprs_sg.apply_tangent_formula(sim.e_values, sim.phases)
      p1_new_phases = tprs_p1.apply_tangent_formula(p1_e_values, p1_phases)
      ref_p1_H = flex.miller_Index()
      ref_p1_e_values = flex.double()
      ref_p1_phases = flex.double()
      miller.expand_to_p1(
        sim.xtal.SgOps, 1,
        sim.miller_set.H, sim.e_values, sg_new_phases,
        ref_p1_H, ref_p1_e_values, ref_p1_phases)
      js = miller.join_sets(ref_p1_H, p1_H)
      assert not js.have_singles()
      for i,j in js.pairs():
        phase_error = debug_utils.phase_error(
          ref_p1_phases[i], p1_new_phases[j], deg=0)
        if (phase_error >= 1.e-4):
          print "Error: phase mismatch" # XXX assert
        if (0 or verbose or phase_error >= 1.e-4):
          assert ref_p1_H[i] == p1_H[j]
          h = str(ref_p1_H[i])
          print "%-15s %.2f" % (h, p1_phases[i]*180/math.pi)
          print "%-15s %.2f sg" % ("", ref_p1_phases[i]*180/math.pi)
          print "%-15s %.2f p1" % ("", p1_new_phases[j]*180/math.pi)
          print "tprs_p1.n_relations():", tprs_p1.n_relations(j)
      if (0 or verbose):
        print
  if (exercise_squaring):
    new_phases = raise_emap(
      sim.xtal, sim.e000, p1_H, sim.miller_set.H, sim.e_values, sim.phases,
      use_e000, zero_out_negative)
    tprs_plus = [tprs_sg]
    if (tprs_sg != None): tprs_plus.append(None)
    for t in tprs_plus:
      mwpe = compute_mean_weighted_phase_error(
        t,
        sim.miller_set.SgOps, sim.miller_set.H,
        sim.e_values, sim.phases, new_phases,
        verbose)
      print "squaring mean weighted phase error: %.2f" % (mwpe,)
  print

# XXX consolidate with xutils.fft_map
class map_from_ampl_phases:

  def __init__(self, xtal, d_min,
               grid_resolution_factor = 1./3):
    self.xtal = xtal
    max_q = 1. / (d_min**2)
    max_prime = 5
    mandatory_grid_factors = xtal.SgOps.refine_gridding()
    grid_logical = sftbx.determine_grid(
      xtal.UnitCell,
      max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
    self.rfft = fftbx.real_to_complex_3d(grid_logical)
    self.tags = sftbx.grid_tags(self.rfft.Nreal())
    self.sym_flags = sftbx.map_symmetry_flags(1)
    self.tags.build(xtal.SgInfo, self.sym_flags)
    print "tags.n_independent:", self.tags.n_independent()

  def __call__(self, miller_indices, ampl, phases, deg=0):
    f = flex.polar(ampl, phases, deg)
    friedel_flag = 1
    n_complex = self.rfft.Ncomplex()
    conjugate = 1
    map = sftbx.structure_factor_map(
      self.xtal.SgOps, friedel_flag, miller_indices,
      f, n_complex, conjugate)
    return self.rfft.backward(map)

  def get_peak_list(self, map, peak_search_level, max_peaks):
    flex_utils.inplace_unpad(map)
    return sftbx.get_peak_list(map, self.tags, peak_search_level, max_peaks)

def peak_list_as_emma_model(xtal, n_real, peak_list):
  positions = []
  i_label = 0
  for peak in peak_list:
    i_label += 1
    label = "peak%d" % (i_label,)
    coordinates = [peak["index"][i] / float(n_real[i]) for i in xrange(3)]
    positions.append(emma.labeled_position(label, coordinates))
  return emma.model(xtal, positions)

def recycle(SgInfo,
            number_of_point_atoms = 10,
            d_min=1.,
            e_min=1.2,
            use_triplets=0,
            other_than_sigma_2=0,
            use_squaring=0,
            use_e000=0,
            zero_out_negative=0,
            n_trials=10,
            n_cycles_per_trial=10,
            verbose=0):
  sim = simulated_data(SgInfo, number_of_point_atoms, d_min, e_min, verbose)
  sim_emma_model = sim.xtal.as_emma_model()
  if (0 or verbose):
    sim_emma_model.show("Random model")
  p1_H = flex.miller_Index()
  miller.expand_to_p1(sim.xtal.SgOps, 1, sim.miller_set.H, p1_H)
  map_calculator = map_from_ampl_phases(sim.xtal, d_min)
  print "Nreal:", map_calculator.rfft.Nreal()
  print "Mreal:", map_calculator.rfft.Mreal()
  LookupSymbol = SgInfo.BuildLookupSymbol()
  for i_trial in xrange(n_trials):
    random_phi = debug_utils.random_phases(
      sim.xtal.SgOps, sim.miller_set.H, sim.e_values)
    if (0 or verbose):
      show_ampl_phases(sim.miller_set.H, sim.e_values, random_phi)
    for exponent in [2]:
      new_phases = random_phi
      if (0): # XXX
        new_phases = sim.phases
      n_matches = []
      for cycle in xrange(n_cycles_per_trial):
        print LookupSymbol,
        print "i_trial:", i_trial,
        print "exponent:", exponent,
        print "cycle:", cycle
        sys.stdout.flush()
        new_phases = raise_emap(
          sim.xtal, sim.e000,
          p1_H, sim.miller_set.H, sim.e_values, new_phases,
          use_e000, zero_out_negative, exponent=exponent)
        use_emma = (cycle == 0 or cycle == n_cycles_per_trial-1)
        if (use_emma):
          if (0): # XXX
            new_phases = sim.phases
          map = map_calculator(sim.miller_set.H, sim.e_values, new_phases)
          peak_list = map_calculator.get_peak_list(map, 3, 20)
          if (0 or verbose):
            for peak in peak_list:
              print "peak (%d,%d,%d)" % peak["index"], "%.6g" % (
                peak["value"],)
          peaks_emma_model = peak_list_as_emma_model(
            sim.xtal, map_calculator.rfft.Nreal(), peak_list)
          if (0 or verbose):
            peaks_emma_model.show("Peak list")
          refined_matches = emma.match_models(sim_emma_model, peaks_emma_model)
          if (len(refined_matches)):
            refined_matches[0].show()
            n_matches.append(len(refined_matches[0].pairs))
          else:
            print "No matches"
            n_matches.append(0)
          sys.stdout.flush()
      print LookupSymbol,
      print "i_trial:", i_trial,
      print "exponent:", exponent,
      print "n_matches:", n_matches
      print
      sys.stdout.flush()

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "IncludeVeryHighSymmetry",
    "ShowSymbolOnly",
    "triplets",
    "sigma_2",
    "squaring",
    "use_e000",
    "zero_out_negative",
    "recycle",
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
    if (SgInfo.SgOps().isCentric() and Flags.triplets):
      print "Centric space group skipped."
      print
      sys.stdout.flush()
      continue
    if (SgInfo.SgOps().OrderZ() > 24 and not Flags.IncludeVeryHighSymmetry):
      print "High symmetry space group skipped."
      print
      sys.stdout.flush()
      continue
    if (Flags.ShowSymbolOnly):
      print "Space group:", LookupSymbol
    elif (Flags.recycle):
      recycle(
        SgInfo,
        use_triplets=Flags.triplets,
        other_than_sigma_2=not Flags.sigma_2,
        use_squaring=Flags.squaring,
        use_e000=Flags.use_e000,
        zero_out_negative=Flags.zero_out_negative)
    else:
      exercise(
        SgInfo,
        exercise_triplets=Flags.triplets,
        other_than_sigma_2=not Flags.sigma_2,
        exercise_squaring=Flags.squaring,
        use_e000=Flags.use_e000,
        zero_out_negative=Flags.zero_out_negative)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
