from iotbx import reflection_file_reader
from iotbx.option_parser import iotbx_option_parser
from libtbx.itertbx import count
import sys

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
    .option(None, "--resolution",
      action="store",
      type="float",
      dest="resolution",
      help="High resolution limit",
      metavar="FLOAT")
    .option(None, "--bins",
      action="store",
      type="int",
      dest="n_bins",
      default=10,
      help="Number of bins",
      metavar="INT")
  ).process()
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  all_miller_arrays = []
  all_anom_diffs = []
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
          elif (miller_array.is_complex()):
            miller_array = abs(miller_array)
          if (miller_array.is_real()):
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
            if (command_line.options.resolution is not None):
              miller_array = miller_array.resolution_filter(
                d_min=command_line.options.resolution)
            miller_array.set_info(info=info)
            all_miller_arrays.append(miller_array)
            if (miller_array.anomalous_flag()):
              anom_diffs = miller_array.anomalous_differences()
            else:
              anom_diffs = None
            all_anom_diffs.append(anom_diffs)
  assert len(all_miller_arrays) == len(all_anom_diffs)
  for i_0,array_0,anom_diffs_0 in zip(count(), all_miller_arrays,
                                               all_anom_diffs):
    array_0.show_comprehensive_summary()
    print
    array_0.setup_binner(n_bins=command_line.options.n_bins)
    if (anom_diffs_0 is not None):
      anom_diffs_0.setup_binner(n_bins=command_line.options.n_bins)
    print "Completeness of %s:" % str(array_0.info())
    array_0.show_completeness_in_bins()
    print
    if (array_0.anomalous_flag()):
      print "Anomalous signal of %s:" % str(array_0.info())
      print array_0.anomalous_signal.__doc__
      anom_signal = array_0.anomalous_signal(use_binning=0001)
      anom_signal.show(data_fmt="%.4f")
      print
    for array_1,anom_diffs_1 in zip(all_miller_arrays[i_0+1:],
                                    all_anom_diffs[i_0+1:]):
      if (not array_0.is_similar_symmetry(array_1)):
        print "Incompatible symmetries:"
        print " ", array_0.info()
        print " ", array_1.info()
        print "No comparison."
        print
      elif (not command_line.options.quick):
        print "Correlation of:"
        print " ", array_0.info()
        print " ", array_1.info()
        correlation = array_0.correlation(
          other=array_1)
        print "Overall correlation: %6.3f" % correlation.coefficient()
        correlation = array_0.correlation(
          other=array_1,
          use_binning=0001)
        correlation.show(data_fmt=binned_correlation_fmt)
        print
        if (anom_diffs_0 is not None and anom_diffs_1 is not None):
          print "Anomalous difference correlation of:"
          print " ", array_0.info()
          print " ", array_1.info()
          correlation = anom_diffs_0.correlation(
            other=anom_diffs_1)
          print "Overall correlation: %6.3f" % correlation.coefficient()
          correlation = anom_diffs_0.correlation(
            other=anom_diffs_1,
            use_binning=0001)
          correlation.show(data_fmt=binned_correlation_fmt)
          print
    print "=" * 79
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
