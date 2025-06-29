"""Map comparison and statistics"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_comparison

from cctbx import crystal
from cctbx import maptbx, miller
from cctbx.sgtbx import space_group_info
from iotbx import file_reader, phil
import iotbx.ccp4_map
from libtbx.utils import Sorry
from scitbx.array_family import flex
import os, sys
from six.moves import zip
from six.moves import range

master_phil = phil.parse("""
include scope libtbx.phil.interface.tracking_params
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
  mtz_1 = None
    .type = path
    .short_caption = Map 1
    .help = MTZ file containing map
    .style = file_type:hkl bold input_file process_hkl child:map_labels:mtz_label_1
  mtz_2 = None
    .type = path
    .short_caption = Map 2
    .help = MTZ file containing map
    .style = file_type:hkl bold input_file process_hkl child:map_labels:mtz_label_2
  mtz_label_1 = None
    .type = str
    .short_caption = Data label
    .help = Data label for complex map coefficients in MTZ file
    .style = renderer:draw_map_arrays_widget
  mtz_label_2 = None
    .type = str
    .short_caption = Data label
    .help = Data label for complex map coefficients in MTZ file
    .style = renderer:draw_map_arrays_widget
}
options
{
  resolution_factor = 0.25
    .type = float
    .short_caption = Resolution gridding factor
    .help = Determines grid spacing in map
  shift_origin = True
    .type = bool
    .short_caption = Shift origin(s) to (0,0,0)
    .help = Shift origin if necessary
  contour_to_match = None
    .type = float
    .short_caption = Contour to match
    .help = Contour level in map1 to match in map2 by volume equalization.

}
""", process_includes=True)

master_params = master_phil

def show_overall_statistics(out=sys.stdout, s=None, header=None):
  print(header, file=out)
  print("  min/max/mean: %6.4f %6.4f %6.4f"%(s.min(), s.max(), s.mean()), file=out)
  print("  kurtosis    : %6.4f" % s.kurtosis(), file=out)
  print("  skewness    : %6.4f" % s.skewness(), file=out)
  print("  sigma       : %6.4f" % s.sigma(), file=out)

def create_statistics_dict(out=sys.stdout, s=None):
  statistics_dict = dict()
  statistics_dict['min'] = s.min()
  statistics_dict['max'] = s.max()
  statistics_dict['mean'] = s.mean()
  statistics_dict['kurtosis'] = s.kurtosis()
  statistics_dict['skewness'] = s.skewness()
  statistics_dict['sigma'] = s.sigma()
  return statistics_dict

def show_citation(out=sys.stdout):
  print("-"*79, file=out)
  msg = """Map comparison and statistics. For details see:
  Acta Cryst. (2014). D70, 2593-2606
  Metrics for comparison of crystallographic maps
  A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams"""
  print(msg, file=out)
  print("-"*79, file=out)
def match_contour_level(m1=None,m2=None,
       contour_to_match=None,results=None):
  # just find map value in m2 that is bigger than same number of grid points
  #   as contour_to_match is for m1
  s=(m1>=contour_to_match)
  from cctbx.maptbx.segment_and_split_map import find_threshold_in_map
  contour_in_map_2=find_threshold_in_map(target_points=s.count(True),
         map_data=m2)
  s2=(m2>=contour_in_map_2)
  results['matching_contour']=contour_in_map_2
  results['v1']=s.count(True)/s.size()
  results['v2']=s2.count(True)/s2.size()
  return results

# =============================================================================
def run(args, out=sys.stdout, validated=False):
  show_citation(out=out)
  if (len(args) == 0):
    master_phil.show(out=out)
    print('\nUsage: phenix.map_comparison <CCP4> <CCP4>\n',\
      '       phenix.map_comparison <CCP4> <MTZ> mtz_label_1=<label>\n',\
      '       phenix.map_comparison <MTZ 1> mtz_label_1=<label 1> <MTZ 2> mtz_label_2=<label 2>\n', file=out)
    sys.exit()

  # process arguments
  params = None
  input_attributes = ['map_1', 'mtz_1', 'map_2', 'mtz_2']
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

    # check if more files are necessary
    n_defined = 0
    for attribute in input_attributes:
      if (getattr(params.input, attribute) is not None):
        n_defined += 1

    # matches files to phil scope, stops once there is sufficient data
    for map_file in map_files:
      if (n_defined < 2):
        current_map = file_reader.any_file(map_file)
        if (current_map.file_type == 'ccp4_map'):
          n_defined += 1
          if (params.input.map_1 is None):
            params.input.map_1 = map_file
          elif (params.input.map_2 is None):
            params.input.map_2 = map_file
        elif (current_map.file_type == 'hkl'):
          n_defined += 1
          if (params.input.mtz_1 is None):
            params.input.mtz_1 = map_file
          elif (params.input.mtz_2 is None):
            params.input.mtz_2 = map_file
      else:
        print('WARNING: only the first two files are used', file=out)
        break

  # validate arguments (GUI sets validated to true, no need to run again)
  assert (params is not None)
  if (not validated):
    validate_params(params)

  # ---------------------------------------------------------------------------
  # check if maps need to be generated from mtz
  n_maps = 0
  maps = list()
  map_names = list()
  for attribute in input_attributes:
    filename = getattr(params.input, attribute)
    if (filename is not None):
      map_names.append(filename)
      current_map = file_reader.any_file(filename)
      maps.append(current_map)
      if (current_map.file_type == 'ccp4_map'):
        n_maps += 1

  # construct maps, if necessary
  crystal_gridding = None
  m1 = None
  m2 = None
  assert params.options.shift_origin==True
  # 1 map, 1 mtz file
  if (n_maps == 1):
    for current_map in maps:
      if (current_map.file_type == 'ccp4_map'):
        uc = current_map.file_object.unit_cell()
        sg_info = space_group_info(current_map.file_object.space_group_number)
        n_real = current_map.file_object.unit_cell_grid
        crystal_gridding = maptbx.crystal_gridding(
          uc, space_group_info=sg_info, pre_determined_n_real=n_real)
        m1 = current_map.file_object.map_data()
        if params.options.shift_origin:
          m1=m1.shift_origin()
    if (crystal_gridding is not None):
      label = None
      for attribute in [('mtz_1', 'mtz_label_1'),
                        ('mtz_2', 'mtz_label_2')]:
        filename = getattr(params.input, attribute[0])
        label = getattr(params.input, attribute[1])
        if ( (filename is not None) and (label is not None) ):
          break
      # labels will match currently open mtz file
      for current_map in maps:
        if (current_map.file_type == 'hkl'):
          m2 = miller.fft_map(
            crystal_gridding=crystal_gridding,
            fourier_coefficients=current_map.file_server.get_miller_array(
              label)).apply_sigma_scaling().real_map_unpadded()
    else:
      raise Sorry('Gridding is not defined.')

  # 2 mtz files
  elif (n_maps == 0):
    crystal_symmetry = get_crystal_symmetry(maps[0])
    d_min = min(get_d_min(maps[0]), get_d_min(maps[1]))
    crystal_gridding = maptbx.crystal_gridding(
      crystal_symmetry.unit_cell(), d_min=d_min,
      resolution_factor=params.options.resolution_factor,
      space_group_info=crystal_symmetry.space_group_info())
    m1 = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=maps[0].file_server.get_miller_array(
        params.input.mtz_label_1)).apply_sigma_scaling().real_map_unpadded()
    m2 = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=maps[1].file_server.get_miller_array(
        params.input.mtz_label_2)).apply_sigma_scaling().real_map_unpadded()

  # 2 maps
  else:
    m1 = maps[0].file_object.map_data()
    m2 = maps[1].file_object.map_data()
    if params.options.shift_origin:
      m1=m1.shift_origin()
      m2=m2.shift_origin()

  # ---------------------------------------------------------------------------
  # analyze maps
  assert ( (m1 is not None) and (m2 is not None) )
  if (list(m1.origin()) != [0,0,0]) or  \
     (list(m2.origin()) != [0,0,0]):
    raise Sorry(
       "Shift_origin must be set if maps do not have origin at (0,0,0)")
  if m1.size() != m2.size():
    raise Sorry ("Maps must be the same size")
  results=dict()
  results['map_files'] = None
  results['map_statistics'] = None
  results['cc_input_maps'] = None
  results['cc_quantile'] = None
  results['cc_peaks'] = None
  results['discrepancies'] = None
  results['map_histograms'] = None

  if params.options.contour_to_match:
    match_contour_level(m1=m1,m2=m2,
       contour_to_match=params.options.contour_to_match,
        results=results)
    print ("Contour level map 1: %.4f (fractional volume of %.3f ) " %(
       params.options.contour_to_match,results['v1']),\
       "\nmatches enclosed volume of "+\
       "contour level map 2 of : %.4f (volume %.3f )" %(
       results['matching_contour'],results['v2']),file=out)
    return results
  # show general statistics
  s1 = maptbx.more_statistics(m1)
  s2 = maptbx.more_statistics(m2)
  show_overall_statistics(out=out, s=s1, header="Map 1 (%s):"%map_names[0])
  show_overall_statistics(out=out, s=s2, header="Map 2 (%s):"%map_names[1])
  cc_input_maps = flex.linear_correlation(x = m1.as_1d(),
                                          y = m2.as_1d()).coefficient()
  print("CC, input maps: %6.4f" % cc_input_maps, file=out)

  # compute CCpeak
  cc_peaks = list()
  m1_he = maptbx.volume_scale(map = m1,  n_bins = 10000).map_data()
  m2_he = maptbx.volume_scale(map = m2,  n_bins = 10000).map_data()
  cc_quantile = flex.linear_correlation(x = m1_he.as_1d(),
                                        y = m2_he.as_1d()).coefficient()
  print("CC, quantile rank-scaled (histogram equalized) maps: %6.4f" % \
    cc_quantile, file=out)
  print("Peak correlation:", file=out)
  print("  cutoff  CCpeak", file=out)
  cutoffs = [i/100.  for i in range(1,90)]+ [i/1000 for i in range(900,1000)]
  for cutoff in cutoffs:
    cc_peak = maptbx.cc_peak(map_1=m1_he, map_2=m2_he, cutoff=cutoff)
    print("  %3.2f   %7.4f" % (cutoff, cc_peak), file=out)
    cc_peaks.append((cutoff, cc_peak))

  # compute discrepancy function (D-function)
  discrepancies = list()
  cutoffs = flex.double(cutoffs)
  df = maptbx.discrepancy_function(map_1=m1_he, map_2=m2_he, cutoffs=cutoffs)
  print("Discrepancy function:", file=out)
  print("  cutoff  D", file=out)
  for c, d in zip(cutoffs, df):
    print("  %3.2f   %7.4f" % (c,d), file=out)
    discrepancies.append((c, d))

  # compute and output histograms
  h1 = maptbx.histogram(map=m1, n_bins=10000)
  h2 = maptbx.histogram(map=m2, n_bins=10000)
  print("Map histograms:", file=out)
  print("Map 1 (%s)     Map 2 (%s)"%\
    (params.input.map_1,params.input.map_2), file=out)
  print("(map_value,cdf,frequency) <> (map_value,cdf,frequency)", file=out)
  for a1,c1,v1, a2,c2,v2 in zip(h1.arguments(), h1.c_values(), h1.values(),
                                h2.arguments(), h2.c_values(), h2.values()):
    print("(%9.5f %9.5f %9.5f) <> (%9.5f %9.5f %9.5f)"%\
      (a1,c1,v1, a2,c2,v2), file=out)

  # store results
  s1_dict = create_statistics_dict(s=s1)
  s2_dict = create_statistics_dict(s=s2)
  results = dict()
  inputs = list()
  for attribute in input_attributes:
    filename = getattr(params.input,attribute)
    if (filename is not None):
      inputs.append(filename)
  assert (len(inputs) == 2)
  results['map_files'] = inputs
  results['map_statistics'] = (s1_dict, s2_dict)
  results['cc_input_maps'] = cc_input_maps
  results['cc_quantile'] = cc_quantile
  results['cc_peaks'] = cc_peaks
  results['discrepancies'] = discrepancies
  # TODO, verify h1,h2 are not dicts, e.g. .values is py2/3 compat. I assume it is here
  results['map_histograms'] = ( (h1.arguments(), h1.c_values(), h1.values()),
                                (h2.arguments(), h2.c_values(), h2.values()) )

  return results

# -----------------------------------------------------------------------------
def get_crystal_symmetry(file_handle):
  '''
  Helper function for get crystal symmetry from files
  '''
  file_object = file_handle.file_object
  cs = None
  if (hasattr(file_object, 'space_group_number')):     # CCP4 map
    cs = crystal.symmetry(file_object.unit_cell().parameters(),
                          file_object.space_group_number)
  elif (hasattr(file_object, 'as_miller_arrays')):     # MTZ file
    ma = file_object.as_miller_arrays()
    for a in ma:
      if (a.is_complex_array()):
        cs = a.crystal_symmetry()
        break
    if (cs is None):
      raise Sorry('No map coefficients found in %s.' %
                file_handle.file_name)
  if (cs is None):
    raise Sorry('Could not find crystal symmetry in %s.' %
                file_handle.file_name)
  return cs

def get_d_min(file_handle):
  '''
  Helper function for getting d_min from mtz file
  '''
  miller_arrays = file_handle.file_server.miller_arrays
  d_min = 10.0
  for miller_array in miller_arrays:
    d_min = min(d_min, miller_array.d_min())
  return d_min

def get_mtz_labels(file_handle):
  '''
  Helper function for getting data labels for complex miller arrays
  Returns list of labels for complex arrays
  '''
  miller_arrays = file_handle.file_server.miller_arrays
  labels = list()
  for miller_array in miller_arrays:
    if (miller_array.is_complex_array()):
      labels.append(miller_array.info().label_string())
  return labels

# =============================================================================
# Parameter validation for CLI and GUI
def validate_params(params):

  # check that only 2 files, in any combination, are provided
  input_attributes = ['map_1', 'mtz_1', 'map_2', 'mtz_2']
  n_defined = 0
  for attribute in input_attributes:
    if (getattr(params.input, attribute) is not None):
      n_defined += 1
  if (n_defined != 2):
    raise Sorry('Insufficient data, please provide 2 files' +
                ' (CCP4-formated map or MTZ)')

  # check file type
  maps = list()
  for attribute in input_attributes:
    filename = getattr(params.input, attribute)
    if (filename is not None):
      file_handle = file_reader.any_file(filename)
      if ( (file_handle.file_type != 'ccp4_map') and
           (file_handle.file_type != 'hkl') ):
        raise Sorry('Please input a CCP4-formatted map or MTZ file for %s.'\
                    % filename)
      else:
        maps.append(file_handle)

  # check symmetry
  cs1 = get_crystal_symmetry(maps[0])
  cs2 = get_crystal_symmetry(maps[1])
  if (cs1.is_similar_symmetry(cs2) is False):
    raise Sorry('The symmetry of the two files is not similar.')

  # check gridding if 2 map files are provided
  if ( (maps[0].file_type == 'ccp4_map') and
       (maps[1].file_type == 'ccp4_map') ):
    m1 = maps[0].file_object.map_data()
    m2 = maps[1].file_object.map_data()
    if ( (m1.accessor().all() != m2.accessor().all()) or
         (m1.accessor().focus() != m2.accessor().focus()) or
         (m1.accessor().origin() != m2.accessor().origin()) ):
      raise Sorry('The gridding of the two maps is not compatible.')
  else:
  # check if MTZ files have complex arrays and labels
    for i in range(len(maps)):
      if (maps[i].file_type == 'hkl'):
        labels = get_mtz_labels(maps[i])
        if (len(labels) == 0):
          raise Sorry('%s does not have complex map coefficients' %
                      maps[i].file_name)
        label_phil = getattr(params.input, 'mtz_label_' + str(i+1))
        if (label_phil is None):
          raise Sorry('No labels were specified for %s.' % maps[i].file_name)
        elif (label_phil not in labels):
          raise Sorry('%s does not exist in %s' %
                      (label_phil, maps[i].file_name))

  # check for valid resolution gridding
  if (params.options.resolution_factor < 0.0):
    raise Sorry(
      'Please use a positive value for the resolution gridding factor.')
  return True

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher(runtime_utils.target_with_save_result):
  def run(self):
    result = run(args=self.args, validated=True, out=sys.stdout)
    return result

# =============================================================================
if (__name__ == "__main__"):
  run(sys.argv[1:])

