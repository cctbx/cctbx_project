from iotbx import reflection_file_reader
from iotbx.option_parser import iotbx_option_parser
from cctbx import crystal
from cctbx import sgtbx
from libtbx.itertbx import count
import sys

class array_cache:

  def __init__(self, input):
    self.input = input
    self.change_of_basis_op_to_minimum_cell \
      = self.input.change_of_basis_op_to_minimum_cell()
    self.observations = self.input.change_basis(
      cb_op=self.change_of_basis_op_to_minimum_cell) \
        .expand_to_p1() \
        .map_to_asu()
    if (self.input.anomalous_flag()):
      self.anom_diffs = self.input.anomalous_differences().change_basis(
        cb_op=self.change_of_basis_op_to_minimum_cell) \
          .expand_to_p1() \
          .map_to_asu()
    else:
      self.anom_diffs = None
    self.minimum_cell_symmetry = crystal.symmetry.change_basis(
      self.input,
      cb_op=self.change_of_basis_op_to_minimum_cell)
    self.patterson_group = self.minimum_cell_symmetry.space_group() \
      .build_derived_patterson_group()
    self.patterson_group.make_tidy()

  def similarity_transformations(self,
        other,
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=2):
    c_inv_rs = self.minimum_cell_symmetry.unit_cell() \
      .similarity_transformations(
        other=other.minimum_cell_symmetry.unit_cell(),
        relative_length_tolerance=relative_length_tolerance,
        absolute_angle_tolerance=absolute_angle_tolerance)
    expanded_groups = [self.patterson_group]
    if (other.patterson_group != self.patterson_group):
      expanded_groups.append(other.patterson_group)
    patterson_groups = tuple(expanded_groups)
    result = []
    for c_inv_r in c_inv_rs:
      c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r))
      if (c_inv.is_unit_mx()):
        result.append(sgtbx.change_of_basis_op(c_inv))
      else:
        for patterson_group in patterson_groups:
          expanded_group = sgtbx.space_group(patterson_group)
          expanded_group.expand_smx(c_inv)
          expanded_group.make_tidy()
          def is_in_expanded_groups():
            for g in expanded_groups:
              if (g == expanded_group): return True
            return False
          if (not is_in_expanded_groups()):
            expanded_groups.append(expanded_group)
            result.append(sgtbx.change_of_basis_op(c_inv).inverse())
    return result

  def combined_cb_op(self, other, cb_op):
    s = self.change_of_basis_op_to_minimum_cell
    o = other.change_of_basis_op_to_minimum_cell
    return s.inverse() * cb_op.new_denominators(s) * o

def binned_correlation_fmt(correlation):
  return "%6.3f" % correlation.coefficient()

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
    if (reflection_file.file_type() is not None):
      try:
        miller_arrays = reflection_file.as_miller_arrays(
          crystal_symmetry=command_line.symmetry)
      except:
        print >> sys.stderr, "Warning: unknown file format:", file_name
        print >> sys.stderr
        sys.stderr.flush()
        miller_array = None
      if (miller_arrays is not None):
        for miller_array in miller_arrays:
          info = miller_array.info()
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
            array_caches.append(array_cache(input=miller_array))
  n_bins = command_line.options.n_bins
  for i_0,cache_0 in enumerate(array_caches):
    cache_0.input.show_comprehensive_summary()
    print
    reindexing_matrices = cache_0.input.reindexing_matrices()
    print "Possible twin laws:",
    if (len(reindexing_matrices) == 0):
      print "None"
    else:
      print
      for c in reindexing_matrices:
        print " ", c.r().as_hkl()
      print
    print "Completeness of %s:" % str(cache_0.input.info())
    if (cache_0.input.binner() is None):
      cache_0.input.setup_binner(n_bins=n_bins)
    cache_0.input.show_completeness_in_bins()
    print
    if (cache_0.input.anomalous_flag()):
      print "Anomalous signal of %s:" % str(cache_0.input.info())
      print cache_0.input.anomalous_signal.__doc__
      anom_signal = cache_0.input.anomalous_signal(use_binning=True)
      anom_signal.show(data_fmt="%.4f")
      print
    if (not command_line.options.quick):
      for info_1 in array_caches[i_0+1:]:
        similarity_transformations = cache_0.similarity_transformations(
          other=info_1,
          relative_length_tolerance=0.05,
          absolute_angle_tolerance=5)
        if (len(similarity_transformations) == 0):
          print "Incompatible symmetries:"
          print " ", cache_0.input.info()
          print " ", info_1.input.info()
          print "No comparison."
          print
        else:
          for cb_op in similarity_transformations:
            print "Correlation of:"
            print " ", cache_0.input.info()
            print " ", info_1.input.info()
            similar_array_1 = info_1.observations \
              .change_basis(cb_op) \
              .map_to_asu()
            correlation = cache_0.observations.correlation(
              other=similar_array_1,
              assert_is_similar_symmetry=False)
            combined_cb_op = cache_0.combined_cb_op(other=info_1, cb_op=cb_op)
            if (not combined_cb_op.c().is_unit_mx()):
              reindexing_note = " with reindexing"
              hkl_str = " "+combined_cb_op.as_hkl()
            else:
              reindexing_note = ""
              hkl_str = ""
            print "Overall correlation%s: %6.3f%s" % (
              reindexing_note, correlation.coefficient(), hkl_str)
            if (cache_0.observations.binner() is None):
              cache_0.observations.setup_binner(n_bins=n_bins)
            correlation = cache_0.observations.correlation(
              other=similar_array_1,
              use_binning=True,
              assert_is_similar_symmetry=False)
            correlation.show(data_fmt=binned_correlation_fmt)
            print
            if (    cache_0.anom_diffs is not None
                and info_1.anom_diffs is not None):
              print "Anomalous difference correlation of:"
              print " ", cache_0.input.info()
              print " ", info_1.input.info()
              similar_anom_diffs_1 = info_1.anom_diffs \
                .change_basis(cb_op) \
                .map_to_asu()
              correlation = cache_0.anom_diffs.correlation(
                other=similar_anom_diffs_1,
                assert_is_similar_symmetry=False)
              print "Overall correlation%s: %6.3f%s" % (
                reindexing_note, correlation.coefficient(), hkl_str)
              if (cache_0.anom_diffs.binner() is None):
                cache_0.anom_diffs.setup_binner(n_bins=n_bins)
              correlation = cache_0.anom_diffs.correlation(
                other=similar_anom_diffs_1,
                use_binning=True,
                assert_is_similar_symmetry=False)
              correlation.show(data_fmt=binned_correlation_fmt)
              print
    print "=" * 79
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
