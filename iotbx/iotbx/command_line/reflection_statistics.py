from iotbx import reflection_file_reader
from iotbx.option_parser import iotbx_option_parser
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.itertbx import count
import math
import sys

class array_cache:

  def __init__(self, input, lattice_symmetry_max_delta):
    self.input = input
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

  def idealized_input_unit_cell(self):
    return self.change_of_basis_op_to_minimum_cell.inverse().apply(
      self.lattice_symmetry.unit_cell())

  def possible_twin_laws(self):
    result = []
    for partition in sgtbx.cosets.left_decomposition(
      g=self.lattice_group,
      h=self.intensity_symmetry.space_group()
          .build_derived_acentric_group()
          .make_tidy()).partitions[1:]:
      if (partition[0].r().determinant() > 0):
        result.append(partition[0])
    return result

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
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  array_caches = []
  for file_name in command_line.args:
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name=file_name)
    miller_arrays = None
    if (reflection_file.file_type() is not None):
      try:
        miller_arrays = reflection_file.as_miller_arrays(
          crystal_symmetry=command_line.symmetry)
      except:
        pass
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
          array_caches.append(array_cache(
            input=miller_array,
            lattice_symmetry_max_delta=0.1))
  if (len(array_caches) > 2):
    print "Array indices (for quick searching):"
    for i_0,cache_0 in enumerate(array_caches):
      print "  %2d:" % (i_0+1), cache_0.input.info()
    print
    print "Useful search patterns are:"
    print "    Summary i"
    print "    CC Obs i j"
    print "    CC Ano i j"
    print "  i and j are the indices shown above."
    print
  n_bins = command_line.options.n_bins
  for i_0,cache_0 in enumerate(array_caches):
    print "Summary", i_0+1
    cache_0.input.show_comprehensive_summary()
    print
    print "Space group of the intensities:", \
      cache_0.intensity_symmetry.space_group_info() \
        .as_reference_setting()
    print "Space group of the metric:     ", \
      cache_0.lattice_symmetry.space_group_info() \
        .as_reference_setting()
    twin_laws = cache_0.possible_twin_laws()
    if (len(twin_laws) == 0):
      print "Possible twin laws: None"
      print
    else:
      s = str(cache_0.idealized_input_unit_cell())
      if (s != str(cache_0.input.unit_cell())):
        print "Idealized unit cell:", s
      print "Possible twin laws:"
      for s in twin_laws:
        print " ", s.r().as_hkl()
      print
    print "Completeness of %s:" % str(cache_0.input.info())
    cache_0.input.setup_binner(n_bins=n_bins)
    cache_0.input.completeness(use_binning=True).show()
    print
    print "Perfect merohedral twinning test for %s:"%str(cache_0.input.info())
    acentric = cache_0.input.select_acentric()
    assert acentric.observation_type() is cache_0.input.observation_type()
    acentric.setup_binner(auto_binning=True).counts_complete(
      include_centric=False)
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
    print
    del acentric
    if (cache_0.input.anomalous_flag()):
      print "Anomalous signal of %s:" % str(cache_0.input.info())
      print cache_0.input.anomalous_signal.__doc__
      anom_signal = cache_0.input.anomalous_signal(use_binning=True)
      anom_signal.show()
      print
    if (not command_line.options.quick):
      for j_1,cache_1 in enumerate(array_caches[i_0+1:]):
        i_1 = j_1+i_0+1
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
    print "=" * 79
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
