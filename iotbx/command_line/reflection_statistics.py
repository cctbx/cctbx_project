# LIBTBX_SET_DISPATCHER_NAME phenix.reflection_statistics

from iotbx import reflection_file_reader
from iotbx.option_parser import option_parser
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.crystal import reindex
from cctbx.array_family import flex
from libtbx.utils import Sorry
from itertools import count
import sys

def show_average_of_binned_data(binned_data_list):
  l = len(binned_data_list[0].binner.bin_legend(0))
  print " "*(l-9), "average:",
  for binned_data in binned_data_list:
    data = flex.double()
    for d in binned_data.data[1:-1]:
      if (d is not None): data.append(d)
    if (data.size() == 0): print " "*7,
    else: print "%7.4f" % flex.mean(data),
  print

class array_cache(object):

  def __init__(self, input, n_bins, lattice_symmetry_max_delta):
    self.input = input.eliminate_sys_absent(integral_only=True, log=sys.stdout)
    self.lattice_symmetry_max_delta = lattice_symmetry_max_delta
    if (not self.input.is_unique_set_under_symmetry()):
      print "Merging symmetry-equivalent reflections:"
      merged = self.input.merge_equivalents()
      merged.show_summary(prefix="  ")
      print
      self.input = merged.array()
      del merged
      if (input.info() is not None):
        self.input.set_info(input.info().customized_copy(merged=True))
      else:
        self.input.set_info(miller.array_info(merged=True))
    self.input.show_comprehensive_summary()
    print
    self.input.setup_binner(n_bins=n_bins)
    self.resolution_range = self.input.resolution_range()
    self.change_of_basis_op_to_minimum_cell \
      = self.input.change_of_basis_op_to_minimum_cell()
    self.observations = self.input.change_basis(
      cb_op=self.change_of_basis_op_to_minimum_cell) \
        .expand_to_p1() \
        .map_to_asu()
    if (self.input.anomalous_flag()):
      self.anom_diffs = abs(self.input.anomalous_differences()).change_basis(
        cb_op=self.change_of_basis_op_to_minimum_cell) \
          .expand_to_p1() \
          .map_to_asu()
    else:
      self.anom_diffs = None
    self.minimum_cell_symmetry = crystal.symmetry.change_basis(
      self.input,
      cb_op=self.change_of_basis_op_to_minimum_cell)
    self.intensity_symmetry = \
      self.minimum_cell_symmetry.reflection_intensity_symmetry(
        anomalous_flag=self.input.anomalous_flag())
    self.lattice_group = sgtbx.lattice_symmetry.group(
      self.minimum_cell_symmetry.unit_cell(),
      max_delta=self.lattice_symmetry_max_delta)
    self.lattice_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    self.lattice_group.make_tidy()
    self.lattice_symmetry = crystal.symmetry(
      unit_cell=self.minimum_cell_symmetry.unit_cell(),
      space_group_info=sgtbx.space_group_info(group=self.lattice_group),
      assert_is_compatible_unit_cell=False)

  def show_completeness(self):
    print "Completeness of %s:" % str(self.input.info())
    no_sys_abs = self.input.eliminate_sys_absent()
    no_sys_abs.use_binning_of(self.input)
    no_sys_abs.completeness(use_binning=True).show()
    print

  def idealized_input_unit_cell(self):
    return self.lattice_symmetry.unit_cell().change_basis(
      cb_op=self.change_of_basis_op_to_minimum_cell.inverse())

  def possible_twin_laws(self):
    cosets = sgtbx.cosets.left_decomposition_point_groups_only(
      g=self.lattice_group,
      h=self.intensity_symmetry.space_group()
          .build_derived_acentric_group()
          .make_tidy())
    return cosets.best_partition_representatives(
      cb_op=self.change_of_basis_op_to_minimum_cell.inverse(),
      omit_first_partition=True,
      omit_negative_determinants=True)

  def show_possible_twin_laws(self):
    print "Space group of the intensities:", \
      self.intensity_symmetry.space_group_info() \
        .as_reference_setting()
    print "Space group of the metric:     ", \
      self.lattice_symmetry.space_group_info() \
        .as_reference_setting()
    d = "%.6g" % self.lattice_symmetry_max_delta
    if (d == "1"): d += " degree"
    else: d += " degrees"
    print "  Tolerance used in the determination of the"
    print "  lattice symmetry:", d
    twin_laws = self.possible_twin_laws()
    if (len(twin_laws) == 0):
      print "Possible twin laws: None"
    else:
      idealized_cell = self.idealized_input_unit_cell()
      s = str(idealized_cell)
      if (s != str(self.input.unit_cell())):
        print "Idealized unit cell:", s
      print "Possible twin laws:"
      for s in twin_laws:
        hkl_str = s.r().as_hkl()
        cb_op = sgtbx.change_of_basis_op(hkl_str)
        assert idealized_cell.change_basis(
          cb_op=cb_op).is_similar_to(idealized_cell)
        info = ""
        try:
          cb_sg = self.input.space_group_info().change_basis(cb_op=cb_op)
        except RuntimeError, e:
          if (str(e).find("Unsuitable value for rational") < 0): raise
          info = " # Info: this changes the setting of the input space group"
        else:
          if (cb_sg.group() != self.input.space_group()):
            info = " # Info: new setting: %s" % str(cb_sg)
        print " ", hkl_str + info
      print "  Note:"
      print "    phenix.xtriage provides comprehensive twinning analysis"
      print "    facilities for macromolecular structures."
      print "    For more information enter: phenix.xtriage --help"
    print

  def show_patterson_peaks(self,
        min_relative_peak_height=0.1,
        show_at_least=3):
    print "Patterson peaks for %s:" % str(self.input.info())
    reciprocal_map = self.input
    if (reciprocal_map.anomalous_flag()):
      reciprocal_map = reciprocal_map.average_bijvoet_mates()
    patterson_map = reciprocal_map.patterson_map(
      symmetry_flags=maptbx.use_space_group_symmetry)
    patterson_map.apply_sigma_scaling()
    peak_list = patterson_map.tags().peak_search(
      map=patterson_map.real_map(),
      parameters=maptbx.peak_search_parameters())
    max_height = peak_list.heights()[0]
    sym_equiv_origin = sgtbx.sym_equiv_sites(
      unit_cell=patterson_map.unit_cell(),
      space_group=patterson_map.space_group(),
      original_site=(0,0,0))
    print "      Fractional coordinates     Height  Distance from origin"
    for i_peak in xrange(peak_list.size()):
      height = peak_list.heights()[i_peak]
      if (height < max_height * min_relative_peak_height
          and i_peak > show_at_least): break
      site = peak_list.sites()[i_peak]
      dist_info = sgtbx.min_sym_equiv_distance_info(sym_equiv_origin, site)
      print "  %8.4f %8.4f %8.4f" % (dist_info.sym_op()*site),
      print "  %8.3f  %8.3f" % (height, dist_info.dist())
    print

  def show_perfect_merohedral_twinning_test(self, n_bins=None):
    assert not self.input.space_group().is_centric()
    print "Perfect merohedral twinning test for %s:"%str(self.input.info())
    acentric = self.input.select_acentric().as_intensity_array()
    acentric = acentric.array(
      data=acentric.data()/acentric.epsilons().data().as_double())
    acentric.set_observation_type_xray_intensity()
    if (n_bins is not None):
      acentric.setup_binner(n_bins=n_bins)
    else:
      acentric.setup_binner(auto_binning=True)
      if (acentric.binner().n_bins_used() > 30):
        acentric.setup_binner(n_bins=30)
    acentric.binner().counts_complete(include_centric=False)
    sm = acentric.second_moment_of_intensities(use_binning=True)
    wr = acentric.wilson_ratio(use_binning=True)
    print acentric.second_moment_of_intensities.__doc__
    print acentric.wilson_ratio.__doc__
    print "See also: http://www.doe-mbi.ucla.edu/Services/Twinning/intro.html"
    for i_bin,s,w in zip(count(), sm.data, wr.data):
      print sm.binner.bin_legend(i_bin),
      for v in s, w:
        if (v is None): print " "*7,
        else: print "%7.4f" % v,
      print
    show_average_of_binned_data([sm, wr])
    print

  def show_second_moments_of_intensities(self, n_bins=None):
    print "Second moments of intensities for %s:"%str(self.input.info())
    f = self.input.as_intensity_array()
    f = f.array(data=f.data()/f.epsilons().data().as_double())
    f.set_observation_type_xray_intensity()
    if (n_bins is not None):
      f.setup_binner(n_bins=n_bins)
    else:
      f.setup_binner(auto_binning=True)
      if (f.binner().n_bins_used() > 30):
        f.setup_binner(n_bins=30)
    print f.second_moment_of_intensities.__doc__.split()[0]
    sm = f.second_moment_of_intensities(use_binning=True)
    sm.show()
    show_average_of_binned_data([sm])
    print

  def show_measurability(self, cutoff=3):
    if (self.input.sigmas() is not None):
      work_array = self.input.select(self.input.sigmas() > 0)
      if (work_array.size() > 0):
        print "Observed measurabilities of %s:" % str(self.input.info())
        print self.input.measurability.__doc__.replace("cutoff", str(cutoff))
        work_array.use_binning_of(self.input)
        meas_obs = work_array.measurability(use_binning=True, cutoff=cutoff)
        meas_obs.show()
      print

  def unique_reindexing_operators(self,
        other,
        relative_length_tolerance,
        absolute_angle_tolerance):
    #make crystal symmetries
    self_xs = crystal.symmetry( unit_cell = self.input.unit_cell(),
                                space_group = self.input.space_group() )
    other_xs = crystal.symmetry( unit_cell = other.input.unit_cell(),
                                space_group = other.input.space_group() )
    # some downstream routines expect things to be in minimum cell
    self_xs = self_xs.change_basis(
      self.change_of_basis_op_to_minimum_cell )
    other_xs = other_xs.change_basis(
      other.change_of_basis_op_to_minimum_cell )

    double_cosets = reindex.reindexing_operators(
      self_xs,
      other_xs,
      relative_length_tolerance,
      absolute_angle_tolerance,
      self.input.anomalous_flag() )
    result = double_cosets.combined_cb_ops()
    return result

  def combined_cb_op(self, other, cb_op):
    sc = self.change_of_basis_op_to_minimum_cell
    oc = other.change_of_basis_op_to_minimum_cell

    cb_op = cb_op.new_denominators(sc)
    best_choice = None
    best_choice_as_hkl = None
    for s_symop in self.minimum_cell_symmetry.space_group():
      s_symop = sgtbx.change_of_basis_op(sgtbx.rt_mx(
        s_symop.r())).new_denominators(sc)
      for o_symop in other.minimum_cell_symmetry.space_group():
        o_symop = sgtbx.change_of_basis_op(sgtbx.rt_mx(
          o_symop.r())).new_denominators(sc)
        possible_choice = sc.inverse() * s_symop * cb_op * o_symop * oc
        possible_choice_as_hkl = possible_choice.as_hkl()
        if (best_choice_as_hkl is None
            or sgtbx.compare_cb_op_as_hkl(
                 best_choice_as_hkl, possible_choice_as_hkl) > 0):
          best_choice = possible_choice
          best_choice_as_hkl = possible_choice_as_hkl
    assert best_choice is not None
    return best_choice

  def setup_common_binner(self,
        other,
        auto_binning=False,
        reflections_per_bin=0,
        n_bins=0,
        d_tolerance=1.e-6):
    d_max = min(self.resolution_range[0], other.resolution_range[0])
    d_min = max(self.resolution_range[1], other.resolution_range[1])
    if (d_max == d_min):
      d_max += d_max*0.5
      d_min -= d_min*0.5
    else:
      d_max *= (1+d_tolerance)
      d_min *= (1-d_tolerance)
    self.observations.setup_binner(
      d_max=d_max,
      d_min=d_min,
      auto_binning=auto_binning,
      reflections_per_bin=reflections_per_bin,
      n_bins=n_bins)
    if (self.anom_diffs is not None and other.anom_diffs is not None):
      self.anom_diffs.setup_binner(
        d_max=d_max,
        d_min=d_min,
        auto_binning=auto_binning,
        reflections_per_bin=reflections_per_bin,
        n_bins=n_bins)

def _process_miller_arrays(
      command_line,
      input_miller_arrays,
      active_miller_arrays):
  n_f_sq_as_f = 0
  for miller_array in input_miller_arrays:
    info = miller_array.info()
    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))
    if (miller_array.indices().size() == 0): continue
    if (miller_array.is_xray_intensity_array()):
      miller_array = miller_array.f_sq_as_f()
      n_f_sq_as_f += 1
    elif (miller_array.is_complex_array()):
      miller_array = abs(miller_array)
    if (miller_array.is_real_array()):
      if (miller_array.unit_cell() is None):
        print
        print "*" * 79
        print "Unknown unit cell parameters:", miller_array.info()
        print "Use --symmetry or --unit_cell to define unit cell:"
        print "*" * 79
        print
        command_line.parser.show_help()
        return -1
      if (miller_array.space_group_info() is None):
        print
        print "*" * 79
        print "Unknown space group:", miller_array.info()
        print "Use --symmetry or --space_group to define space group:"
        print "*" * 79
        print
        command_line.parser.show_help()
        return -1
      if (   command_line.options.resolution is not None
          or command_line.options.low_resolution is not None):
        miller_array = miller_array.resolution_filter(
          d_max=command_line.options.low_resolution,
          d_min=command_line.options.resolution)
      miller_array = miller_array.map_to_asu()
      miller_array.set_info(info=info)
      active_miller_arrays.append(miller_array)
  return n_f_sq_as_f

def run(
      args,
      command_name="phenix.reflection_statistics",
      additional_miller_arrays=[]):
  print "Command line arguments:",
  for arg in args: print arg,
  print
  print
  command_line = (option_parser(
    usage=command_name+" [options] reflection_file [...]",
    description="Example: %s data1.mtz data2.sca" % command_name)
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--quick",
      action="store_true",
      help="Do not compute statistics between pairs of data arrays")
    .enable_resolutions()
    .option(None, "--bins",
      action="store",
      type="int",
      dest="n_bins",
      default=10,
      help="Number of bins",
      metavar="INT")
    .option(None, "--bins_twinning_test",
      action="store",
      type="int",
      dest="n_bins_twinning_test",
      default=None,
      help="Number of bins for twinning test",
      metavar="INT")
    .option(None, "--bins_second_moments",
      action="store",
      type="int",
      dest="n_bins_second_moments",
      default=None,
      help="Number of bins for second moments of intensities",
      metavar="INT")
    .option(None, "--lattice_symmetry_max_delta",
      action="store",
      type="float",
      default=3.,
      help="angular tolerance in degrees used in the determination"
           " of the lattice symmetry")
  ).process(args=args)
  if (len(command_line.args) == 0 and len(additional_miller_arrays) == 0):
    command_line.parser.show_help()
    return
  active_miller_arrays = []
  n_f_sq_as_f = 0
  for file_name in command_line.args:
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name=file_name)
    miller_arrays = None
    if (reflection_file.file_type() is not None):
      try:
        miller_arrays = reflection_file.as_miller_arrays(
          crystal_symmetry=command_line.symmetry,
          force_symmetry=not command_line.options.weak_symmetry,
          merge_equivalents=False)
      except Sorry, KeyboardInterrupt: raise
      except Exception: pass
    if (miller_arrays is None):
      print >> sys.stderr, "Warning: unknown file format:", file_name
      print >> sys.stderr
      sys.stderr.flush()
    else:
      n = _process_miller_arrays(
        command_line=command_line,
        input_miller_arrays=miller_arrays,
        active_miller_arrays=active_miller_arrays)
      if (n < 0): return
      n_f_sq_as_f += n
  if (additional_miller_arrays is not None):
    n = _process_miller_arrays(
      command_line=command_line,
      input_miller_arrays=additional_miller_arrays,
      active_miller_arrays=active_miller_arrays)
    if (n < 0): return
    n_f_sq_as_f += n
  if (n_f_sq_as_f > 0):
    if (n_f_sq_as_f == 1):
      print "Note: Intensity array has been converted to an amplitude array."
    else:
      print "Note: Intensity arrays have been converted to amplitude arrays."
    print
  if (len(active_miller_arrays) > 2 and not command_line.options.quick):
    print "Array indices (for quick searching):"
    for i_0,input_0 in enumerate(active_miller_arrays):
      print "  %2d:" % (i_0+1), input_0.info()
    print
    print "Useful search patterns are:"
    print "    Summary i"
    print "    CC Obs i j"
    print "    CC Ano i j"
    print "  i and j are the indices shown above."
    print
  n_bins = command_line.options.n_bins
  array_caches = []
  for i_0,input_0 in enumerate(active_miller_arrays):
    print "Summary", i_0+1
    print
    cache_0 = array_cache(
      input=input_0,
      n_bins=n_bins,
      lattice_symmetry_max_delta=command_line.options.lattice_symmetry_max_delta)
    cache_0.show_possible_twin_laws()
    cache_0.show_completeness()
    cache_0.show_patterson_peaks()
    if (not cache_0.input.space_group().is_centric()):
      cache_0.show_perfect_merohedral_twinning_test(
        n_bins=command_line.options.n_bins_twinning_test)
    else:
      cache_0.show_second_moments_of_intensities(
        n_bins=command_line.options.n_bins_second_moments)
    if (cache_0.input.anomalous_flag()):
      print "Anomalous signal of %s:" % str(cache_0.input.info())
      print cache_0.input.anomalous_signal.__doc__
      anom_signal = cache_0.input.anomalous_signal(use_binning=True)
      anom_signal.show()
      print
      cache_0.show_measurability()
    for i_1,cache_1 in enumerate(array_caches):
      unique_reindexing_operators = cache_1.unique_reindexing_operators(
        other=cache_0,
        relative_length_tolerance=0.05,
        absolute_angle_tolerance=5)
      if (len(unique_reindexing_operators) == 0):
        print "Incompatible unit cells:"
        print "  %2d:" % (i_1+1), cache_1.input.info()
        print "  %2d:" % (i_0+1), cache_0.input.info()
        print "No comparison."
        print
      else:
        ccs = flex.double()
        for cb_op in unique_reindexing_operators:
          similar_array_0 = cache_0.observations \
            .change_basis(cb_op) \
            .map_to_asu()
          ccs.append(cache_1.observations.correlation(
            other=similar_array_0,
            assert_is_similar_symmetry=False).coefficient())
        permutation = flex.sort_permutation(ccs, reverse=True)
        ccs = ccs.select(permutation)
        unique_reindexing_operators = flex.select(
          unique_reindexing_operators, permutation=permutation)
        for i_cb_op,cb_op,cc in zip(count(),
                                    unique_reindexing_operators,
                                    ccs):
          combined_cb_op = cache_1.combined_cb_op(other=cache_0, cb_op=cb_op)
          if (not combined_cb_op.c().is_unit_mx()):
            reindexing_note = "  after reindexing %d using %s" % (
              i_0+1, combined_cb_op.as_hkl())
          else:
            reindexing_note = ""
          print "CC Obs %d %d %6.3f  %s" % (
            i_1+1, i_0+1, cc, combined_cb_op.as_hkl())
          print "Correlation of:"
          print "  %2d:" % (i_1+1), cache_1.input.info()
          print "  %2d:" % (i_0+1), cache_0.input.info()
          print "Overall correlation: %6.3f%s" % (cc, reindexing_note)
          show_in_bins = False
          if (i_cb_op == 0 or (cc >= 0.3 and cc >= ccs[0]-0.2)):
            show_in_bins = True
            similar_array_0 = cache_0.observations \
              .change_basis(cb_op) \
              .map_to_asu()
            cache_1.setup_common_binner(cache_0, n_bins=n_bins)
            correlation = cache_1.observations.correlation(
              other=similar_array_0,
              use_binning=True,
              assert_is_similar_symmetry=False)
            correlation.show()
          print
          if (    cache_0.anom_diffs is not None
              and cache_1.anom_diffs is not None):
            similar_anom_diffs_0 = cache_0.anom_diffs \
              .change_basis(cb_op) \
              .map_to_asu()
            correlation = cache_1.anom_diffs.correlation(
              other=similar_anom_diffs_0,
              assert_is_similar_symmetry=False)
            print "CC Ano %d %d %6.3f  %s" % (
              i_1+1, i_0+1, correlation.coefficient(), combined_cb_op.as_hkl())
            print "Anomalous difference correlation of:"
            print "  %2d:" % (i_1+1), cache_1.input.info()
            print "  %2d:" % (i_0+1), cache_0.input.info()
            print "Overall anomalous difference correlation: %6.3f%s" % (
              correlation.coefficient(), reindexing_note)
            if (show_in_bins):
              correlation = cache_1.anom_diffs.correlation(
                other=similar_anom_diffs_0,
                use_binning=True,
                assert_is_similar_symmetry=False)
              correlation.show()
            print
    if (not command_line.options.quick):
      array_caches.append(cache_0)
    print "=" * 79
    print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
