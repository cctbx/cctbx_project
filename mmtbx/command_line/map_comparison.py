from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_comparison

from iotbx import file_reader, phil

from cctbx import maptbx
import iotbx.ccp4_map
from cctbx import crystal
from libtbx.utils import Sorry
from scitbx.array_family import flex
import os, sys

master_phil = phil.parse("""
input
{
  map_1 = None
    .type = path
    .short_caption = Map 1
    .help = A CCP4-formatted map
    .style = file_type:ccp4_map bold input_file
  map_2 = None
    .type = path
    .short_caption = Map 2
    .help = A CCP4-formatted map
    .style = file_type:ccp4_map bold input_file
}
""")

master_params = master_phil

def show_overall_statistics(s, header):
  print header
  print "  min/max/mean: %6.4f %6.4f %6.4f"%(s.min(), s.max(), s.mean())
  print "  kurtosis    : %6.4f" % s.kurtosis()
  print "  skewness    : %6.4f" % s.skewness()
  print "  sigma       : %6.4f" % s.sigma()

def create_statistics_dict(s):
  statistics_dict = dict()
  statistics_dict['min'] = s.min()
  statistics_dict['max'] = s.max()
  statistics_dict['mean'] = s.mean()
  statistics_dict['kurtosis'] = s.kurtosis()
  statistics_dict['skewness'] = s.skewness()
  statistics_dict['sigma'] = s.sigma()
  return statistics_dict

def show_citation(out=sys.stdout):
  print >> out, "-"*79
  msg = """Map comparison and statistics. For details see:
  Acta Cryst. (2014). D70, 2593-2606
  Metrics for comparison of crystallographic maps
  A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams"""
  print >> out, msg
  print >> out, "-"*79

# =============================================================================
def run(args, validated=False):
  show_citation()
  if ( (len(args) == 0) or (len(args) > 2) ):
    print '\nUsage: phenix.map_comparison map_1=<first map> map_2=<second map>\n'
    sys.exit()

  # process arguments
  try: # automatic parsing
    params = phil.process_command_line_with_files(
      args=args, master_phil=master_phil).work.extract()
  except Exception: # map_file_def only handles one map phil
    from libtbx.phil.command_line import argument_interpreter
    arg_int = argument_interpreter(master_phil=master_phil)
    command_line_args = list()
    map_files = list()
    for arg in args:
      if (os.path.isfile(arg)):
        map_files.append(arg)
      else:
        command_line_args.append(arg_int.process(arg))
    params = master_phil.fetch(sources=command_line_args).extract()
    for map_file in map_files:
      if (params.input.map_1 is None):
        params.input.map_1 = map_file
      else:
        params.input.map_2 = map_file

  # validate arguments (GUI sets validated to true, no need to run again)
  if (not validated):
    validate_params(params)

  # ---------------------------------------------------------------------------
  # map 1
  ccp4_map_1 = iotbx.ccp4_map.map_reader(file_name=params.input.map_1)
  cs_1 = crystal.symmetry(ccp4_map_1.unit_cell().parameters(),
    ccp4_map_1.space_group_number)
  m1 = ccp4_map_1.map_data()

  # map 2
  ccp4_map_2 = iotbx.ccp4_map.map_reader(file_name=params.input.map_2)
  cs_2 = crystal.symmetry(ccp4_map_2.unit_cell().parameters(),
    ccp4_map_2.space_group_number)
  m2 = ccp4_map_2.map_data()

  # show general statistics
  s1 = maptbx.more_statistics(m1)
  s2 = maptbx.more_statistics(m2)
  show_overall_statistics(s=s1, header="Map 1 (%s):"%params.input.map_1)
  show_overall_statistics(s=s2, header="Map 2 (%s):"%params.input.map_2)
  cc_input_maps = flex.linear_correlation(x = m1.as_1d(),
                                          y = m2.as_1d()).coefficient()
  print "CC, input maps: %6.4f" % cc_input_maps

  # compute CCpeak
  cc_peaks = list()
  m1_he = maptbx.volume_scale(map = m1,  n_bins = 10000).map_data()
  m2_he = maptbx.volume_scale(map = m2,  n_bins = 10000).map_data()
  cc_quantile = flex.linear_correlation(x = m1_he.as_1d(),
                                        y = m2_he.as_1d()).coefficient()
  print "CC, quantile rank-scaled (histogram equalized) maps: %6.4f" % \
    cc_quantile
  print "Peak correlation:"
  print "  cutoff  CCpeak"
  for cutoff in [i/100. for i in range(0,100,5)]+[0.99, 1.0]:
    cc_peak = maptbx.cc_peak(map_1=m1_he, map_2=m2_he, cutoff=cutoff)
    print "  %3.2f   %7.4f" % (cutoff, cc_peak)
    cc_peaks.append((cutoff, cc_peak))

  # compute discrepancy function (D-function)
  discrepancies = list()
  cutoffs = flex.double([i/20. for i in range(1,20)])
  df = maptbx.discrepancy_function(map_1=m1_he, map_2=m2_he, cutoffs=cutoffs)
  print "Discrepancy function:"
  print "  cutoff  D"
  for c, d in zip(cutoffs, df):
    print "  %3.2f   %7.4f" % (c,d)
    discrepancies.append((c, d))

  # compute and output histograms
  h1 = maptbx.histogram(map=m1, n_bins=10000)
  h2 = maptbx.histogram(map=m2, n_bins=10000)
  print "Map histograms:"
  print "Map 1 (%s)     Map 2 (%s)"%(params.input.map_1,params.input.map_2)
  print "(map_value,cdf,frequency) <> (map_value,cdf,frequency)"
  for a1,c1,v1, a2,c2,v2 in zip(h1.arguments(), h1.c_values(), h1.values(),
                                h2.arguments(), h2.c_values(), h2.values()):
    print "(%9.5f %9.5f %9.5f) <> (%9.5f %9.5f %9.5f)"%(a1,c1,v1, a2,c2,v2)

  # store results
  s1_dict = create_statistics_dict(s1)
  s2_dict = create_statistics_dict(s2)
  results = dict()
  results['map_files'] = (params.input.map_1, params.input.map_2)
  results['map_statistics'] = (s1_dict, s2_dict)
  results['cc_input_maps'] = cc_input_maps
  results['cc_quantile'] = cc_quantile
  results['cc_peaks'] = cc_peaks
  results['discrepancies'] = discrepancies
  results['map_histograms'] = ( (h1.arguments(), h1.c_values(), h1.values()),
                                (h2.arguments(), h2.c_values(), h2.values()) )

  return results

# =============================================================================
# Parameter validation for CLI and GUI
def validate_params(params):

  if ( (params.input.map_1 is None) or (params.input.map_2 is None) ):
    raise Sorry('Two CCP4-formatted maps are required.')

  # check files
  p = [params.input.map_1, params.input.map_2]
  maps = [None, None]
  for i in xrange(2):
    maps[i] = file_reader.any_file(p[i])
    if (maps[i].file_type != 'ccp4_map'):
      raise Sorry('Please input a CCP4-formatted map for %s.' % p[i])

  # check symmetry
  m1 = maps[0].file_object
  m2 = maps[1].file_object
  cs1 = crystal.symmetry(m1.unit_cell().parameters(), m1.space_group_number)
  cs2 = crystal.symmetry(m2.unit_cell().parameters(), m2.space_group_number)
  if (cs1.is_similar_symmetry(cs2) is False):
    raise Sorry('The symmetry of the two maps is not similar.')

  # check maps
  m1 = m1.map_data()
  m2 = m2.map_data()
  if ( (m1.accessor().all() != m2.accessor().all()) or
       (m1.accessor().focus() != m2.accessor().focus()) or
       (m1.accessor().origin() != m2.accessor().origin()) ):
    raise Sorry('The two maps are not similar.')

  return True

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    # os.mkdir(self.output_dir)
    # os.chdir(self.output_dir)
    result = run(args=self.args, validated=True)
    return result

# =============================================================================
if (__name__ == "__main__"):
  run(sys.argv[1:])
