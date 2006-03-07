from iotbx import reflection_file_reader
from iotbx.option_parser import iotbx_option_parser
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry
from libtbx.itertbx import count
import math
import sys

def show_average_of_binned_data(binned_data_list):
  l = len(binned_data_list[0].binner.bin_legend(0))
  print " "*(l-9), "average:",
  for binned_data in binned_data_list:
    data = flex.double()
    for d in binned_data.data[1:-1]:
      if (d is not None): data.append(d)
    if (data.size() == 0): print " "*7,
    print "%7.4f" % flex.mean(data),
  print

class array_cache(object):

  def __init__(self, input, n_bins, lattice_symmetry_max_delta):
    self.input = input.eliminate_sys_absent(integral_only=True, log=sys.stdout)
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
      max_delta=lattice_symmetry_max_delta)
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
    return self.change_of_basis_op_to_minimum_cell.inverse().apply(
      self.lattice_symmetry.unit_cell())

  def possible_twin_laws(self):
    result = []
    cb_op = self.change_of_basis_op_to_minimum_cell.inverse()
    for partition in sgtbx.cosets.left_decomposition(
      g=self.lattice_group,
      h=self.intensity_symmetry.space_group()
          .build_derived_acentric_group()
          .make_tidy()).partitions[1:]:
      if (partition[0].r().determinant() > 0):
        result.append(cb_op.apply(partition[0]))
    return result

  def show_possible_twin_laws(self):
    print "Space group of the intensities:", \
      self.intensity_symmetry.space_group_info() \
        .as_reference_setting()
    print "Space group of the metric:     ", \
      self.lattice_symmetry.space_group_info() \
        .as_reference_setting()
    twin_laws = self.possible_twin_laws()
    if (len(twin_laws) == 0):
      print "Possible twin laws: None"
    else:
      idealized_cell = self.idealized_input_unit_cell()
      s = str(idealized_cell)
      if (s != str(self.input.unit_cell())):
        print "Idealized unit cell:", s
      print "Possible twin laws:"
      n_changed_settings = 0
      for s in twin_laws:
        hkl_str = s.r().as_hkl()
        cb_op = sgtbx.change_of_basis_op(hkl_str)
        assert cb_op.apply(idealized_cell).is_similar_to(idealized_cell)
        info = ""
        try:
          cb_sg = self.input.space_group_info().change_basis(cb_op=cb_op)
        except RuntimeError, e:
          if (str(e).find("Unsuitable value for rational") < 0): raise
          info = " # Info: this changes the setting of the input space group"
          n_changed_settings += 1
        else:
          if (cb_sg.group() != self.input.space_group()):
            info = " # Info: new setting: %s" % str(cb_sg)
            n_changed_settings += 1
        print " ", hkl_str + info
      if (n_changed_settings > 0):
        print "  ***************************************************"
        print "  This is a very interesting combination of unit cell"
        print "  parameters and space group symmetry."
        print "  PLEASE send this output to:"
        print
        print "    cctbx@cci.lbl.gov"
        print
        print "  Thank you in advance!"
        print "  ***************************************************"
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
    c_inv_rs = self.minimum_cell_symmetry.unit_cell() \
      .similarity_transformations(
        other=other.minimum_cell_symmetry.unit_cell(),
        relative_length_tolerance=relative_length_tolerance,
        absolute_angle_tolerance=absolute_angle_tolerance)
    min_bases_msd = None
    similarity_cb_op = None
    for c_inv_r in c_inv_rs:
      c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r))
      cb_op = sgtbx.change_of_basis_op(c_inv).inverse()
      bases_msd = self.minimum_cell_symmetry.unit_cell() \
        .bases_mean_square_difference(
          other=cb_op.apply(other.minimum_cell_symmetry.unit_cell()))
      if (min_bases_msd is None
          or min_bases_msd > bases_msd):
        min_bases_msd = bases_msd
        similarity_cb_op = cb_op
    if (similarity_cb_op is None): return []
    common_lattice_group = sgtbx.space_group(self.lattice_group)
    for s in other.lattice_group.build_derived_acentric_group() \
               .change_basis(similarity_cb_op):
      try: common_lattice_group.expand_smx(s)
      except RuntimeError: return []
    common_lattice_group.make_tidy()
    result = []
    for s in sgtbx.cosets.double_unique(
               g=common_lattice_group,
               h1=self.intensity_symmetry.space_group()
                   .build_derived_acentric_group()
                   .make_tidy(),
               h2=other.intensity_symmetry.space_group()
                   .build_derived_acentric_group()
                   .change_basis(similarity_cb_op)
                   .make_tidy()):
      if (s.r().determinant() > 0):
        result.append(sgtbx.change_of_basis_op(s) * similarity_cb_op)
    return result

  def combined_cb_op(self, other, cb_op):
    s = self.change_of_basis_op_to_minimum_cell
    o = other.change_of_basis_op_to_minimum_cell
    return s.inverse() * cb_op.new_denominators(s) * o

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

def run(args):
  print "Command line arguments:",
  for arg in args: print arg,
  print
  print
  command_line = (iotbx_option_parser(
    usage="iotbx.reflection_statistics [options] reflection_file [...]",
    description="Example: iotbx.reflection_statistics data1.mtz data2.sca")
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      dest="weak_symmetry",
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--quick",
      action="store_true",
      dest="quick",
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
  ).process(args=args)
  if (len(command_line.args) == 0):
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
      except: pass
    if (miller_arrays is None):
      print >> sys.stderr, "Warning: unknown file format:", file_name
      print >> sys.stderr
      sys.stderr.flush()
    else:
      for miller_array in miller_arrays:
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
            return
          if (miller_array.space_group_info() is None):
            print
            print "*" * 79
            print "Unknown space group:", miller_array.info()
            print "Use --symmetry or --space_group to define space group:"
            print "*" * 79
            print
            command_line.parser.show_help()
            return
          if (   command_line.options.resolution is not None
              or command_line.options.low_resolution is not None):
            miller_array = miller_array.resolution_filter(
              d_max=command_line.options.low_resolution,
              d_min=command_line.options.resolution)
          miller_array = miller_array.map_to_asu()
          miller_array.set_info(info=info)
          active_miller_arrays.append(miller_array)
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
      lattice_symmetry_max_delta=3.0)
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
    if (not command_line.options.quick):
      for i_1,cache_1 in enumerate(array_caches):
        unique_reindexing_operators = cache_0.unique_reindexing_operators(
          other=cache_1,
          relative_length_tolerance=0.05,
          absolute_angle_tolerance=5)
        if (len(unique_reindexing_operators) == 0):
          print "Incompatible unit cells:"
          print " ", cache_0.input.info()
          print " ", cache_1.input.info()
          print "No comparison."
          print
        else:
          ccs = flex.double()
          for cb_op in unique_reindexing_operators:
            similar_array_1 = cache_1.observations \
              .change_basis(cb_op) \
              .map_to_asu()
            ccs.append(cache_0.observations.correlation(
              other=similar_array_1,
              assert_is_similar_symmetry=False).coefficient())
          permutation = flex.sort_permutation(ccs, reverse=True)
          ccs = ccs.select(permutation)
          unique_reindexing_operators = flex.select(
            unique_reindexing_operators, permutation=permutation)
          for i_cb_op,cb_op,cc in zip(count(),
                                      unique_reindexing_operators,
                                      ccs):
            combined_cb_op = cache_0.combined_cb_op(other=cache_1, cb_op=cb_op)
            if (not combined_cb_op.c().is_unit_mx()):
              reindexing_note = " with reindexing"
              hkl_str = " "+combined_cb_op.as_hkl()
            else:
              reindexing_note = ""
              hkl_str = ""
            print "CC Obs", i_0+1, i_1+1, "%6.3f"%cc, combined_cb_op.as_hkl()
            print "Correlation of:"
            print " ", cache_0.input.info()
            print " ", cache_1.input.info()
            print "Overall correlation%s: %6.3f%s" % (
              reindexing_note, cc, hkl_str)
            show_in_bins = False
            if (i_cb_op == 0 or (cc >= 0.3 and cc >= ccs[0]-0.2)):
              show_in_bins = True
              similar_array_1 = cache_1.observations \
                .change_basis(cb_op) \
                .map_to_asu()
              cache_0.setup_common_binner(cache_1, n_bins=n_bins)
              correlation = cache_0.observations.correlation(
                other=similar_array_1,
                use_binning=True,
                assert_is_similar_symmetry=False)
              correlation.show()
            print
            if (    cache_0.anom_diffs is not None
                and cache_1.anom_diffs is not None):
              similar_anom_diffs_1 = cache_1.anom_diffs \
                .change_basis(cb_op) \
                .map_to_asu()
              correlation = cache_0.anom_diffs.correlation(
                other=similar_anom_diffs_1,
                assert_is_similar_symmetry=False)
              print "CC Ano", i_0+1, i_1+1, \
                "%6.3f"%correlation.coefficient(), combined_cb_op.as_hkl()
              print "Anomalous difference correlation of:"
              print " ", cache_0.input.info()
              print " ", cache_1.input.info()
              print "Overall anomalous difference correlation%s: %6.3f%s" % (
                reindexing_note, correlation.coefficient(), hkl_str)
              if (show_in_bins):
                correlation = cache_0.anom_diffs.correlation(
                  other=similar_anom_diffs_1,
                  use_binning=True,
                  assert_is_similar_symmetry=False)
                correlation.show()
              print
      array_caches.append(cache_0)
    print "=" * 79
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
