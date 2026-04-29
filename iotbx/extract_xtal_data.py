"""
Extract crystallographic data from an array.
  Encapsulates logic for extracting experimental amplitudes and R-free flags
  from the given input file(s) via using reflection_file_server.
"""
from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx.array_family import flex
from libtbx.utils import Sorry
import iotbx.phil
from iotbx import reflection_file_utils
from libtbx import adopt_init_args
import os
from cctbx import french_wilson
import libtbx.callbacks # import dependency
from six.moves import zip
from six.moves import range
from libtbx import group_args
from six.moves import cStringIO as StringIO

def miller_array_symmetry_safety_check(miller_array, data_description,
                                       working_point_group):
  """
  Returns None or string.
  """
  return miller_array.crystal_symmetry_is_compatible_with_symmetry_from_file(
    working_point_group = working_point_group).format_error_message(
      data_description = data_description)

def explain_how_to_generate_array_of_r_free_flags(scope = '<parent scope>.r_free_flags.generate'):
  part1 = """\
If previously used R-free flags are available run this command again
with the name of the file containing the original flags as an
additional input. If the structure was never refined before, or if the
original R-free flags are unrecoverable, run this command again with
the additional definition:

"""
  part3 = """

If the structure was refined previously using different R-free flags,
the values for R-free will become meaningful only after many cycles of
refinement.
"""
  return part1 + """%s=True""" %(scope) + part3

data_and_flags_str_part1a = """\
  file_name = None
    .type=path
    .short_caption=Reflections file
    .style = bold input_file file_type:hkl noauto process_hkl \
      child:fobs:labels child:d_min:high_resolution \
      child:d_max:low_resolution child:rfree_file:r_free_flags.file_name
    .expert_level = 0
  labels = None
    .type=strings
    .input_size = 160
    .short_caption = Data labels
    .style = bold renderer:draw_fobs_label_widget noauto \
      OnChange:auto_update_label_choice child:d_min:high_resolution \
      child:d_max:low_resolution parent:file_name:file_name
    .expert_level = 0
"""

data_and_flags_str_part1b = """\
  high_resolution = None
    .type=float
    .input_size = 80
    .style = bold renderer:draw_resolution_widget noauto
    .expert_level = 0
  low_resolution = None
    .type=float
    .input_size = 80
    .style = bold renderer:draw_resolution_widget noauto
    .expert_level = 0
  twin_law = None
    .type=str
    .input_size = 80
    .style = bold noauto
  outliers_rejection = True
    .type=bool
    .short_caption = Reject outliers
    .help = Remove "basic wilson outliers", "extreme wilson outliers", and \
              "beamstop shadow outliers"
    .expert_level = 0
  french_wilson_scale = True
    .type=bool
    .short_caption = use French-Wilson method to handle negative intensities
  french_wilson
  {
     include scope cctbx.french_wilson.master_phil
  }
  sigma_fobs_rejection_criterion = None
    .type=float
    .short_caption = Sigma(Fobs) rejection criterion
    .expert_level = 0
  sigma_iobs_rejection_criterion = None
    .type=float
    .short_caption = Sigma(Iobs) rejection criterion
    .expert_level = 0
"""

data_and_flags_str_part1 = """\
  %s
  %s
"""%(data_and_flags_str_part1a, data_and_flags_str_part1b)

data_and_flags_str_part2a = """\
  file_name = None
    .type=path
    .short_caption=File with R(free) flags
    .help = This is normally the same as the file containing Fobs and is \
      usually selected automatically.
    .input_size = 200
    .style = noauto input_file file_type:hkl process_hkl child:rfree:label
    .expert_level = 0
  label = None
    .type=str
    .short_caption = R-free label
    .input_size = 160
    .style = bold renderer:draw_rfree_label_widget noauto \
             OnChange:update_rfree_flag_value
    .expert_level = 0
"""

data_and_flags_str_part2b = """\
  required = True
    .type = bool
    .help = Specify if free-r flags are must be present (or else generated)
  test_flag_value = None
    .type=int
    .help = This value is usually selected automatically - do not change \
      unless you really know what you're doing!
    .style = noauto
    .expert_level = 0
  ignore_r_free_flags = False
    .type=bool
    .short_caption = Ignore R-free flags
    .help = Use all reflections in refinement (work and test)
    .expert_level=0
"""

data_and_flags_str_part2 = """\
  %s
  %s
"""%(data_and_flags_str_part2a, data_and_flags_str_part2b)

misc1 = """\
  ignore_all_zeros = True
    .type=bool
    .short_caption = Ignore all-zero arrays
    .expert_level = 1
  force_anomalous_flag_to_be_equal_to = None
    .type=bool
    .short_caption = Use anomalous data
    .style = tribool
    .expert_level = 1
  convert_to_non_anomalous_if_ratio_pairs_lone_less_than_threshold=0.5
    .type = float
    .expert_level = 2
  r_free_flags
    .expert_level=0
    .style = box auto_align
    .caption = This information will be extracted automatically if possible. \
      If no test set is present in the reflections file, one can be generated \
      automatically, or you can use the reflection file editor to combine an \
      existing set with your X-ray or neutron data.
"""
misc2 = """\
    disable_suitability_test = False
      .type=bool
      .expert_level = 2
    ignore_pdb_hexdigest = False
      .type=bool
      .short_caption = Ignore PDB hexdigest sanity check
      .help=If True, disables safety check based on MD5 hexdigests stored in \
            PDB files produced by previous runs.
      .expert_level=2
    generate = False
      .type=bool
      .short_caption = Generate new R-free flags
      .help = Generate R-free flags (if not available in input files)
"""

data_and_flags_str = """\
  %s
  %s
  {
    %s
    %s
    %s
  }
""" % (data_and_flags_str_part1,
       misc1,
       data_and_flags_str_part2,
       misc2,
       miller.generate_r_free_params_str)

data_and_flags_str_no_filenames = """\
  %s
  %s
  {
    %s
    %s
    %s
  }
""" % (data_and_flags_str_part1b,
       misc1,
       data_and_flags_str_part2b,
       misc2,
       miller.generate_r_free_params_str)

xray_data_str = """\
xray_data
  .help=Scope of X-ray data and free-R flags
  .style = scrolled auto_align
{
  %s
}
"""%data_and_flags_str

xray_data_str_no_filenames = """\
xray_data
  .help=Scope of X-ray data and free-R flags
  .style = scrolled auto_align
{
  %s
}
"""%data_and_flags_str_no_filenames

neutron_data_str = """\
neutron_data
  .help=Scope of neutron data and neutron free-R flags
  .style = scrolled auto_align
{
  ignore_xn_free_r_mismatch = False
    .type = bool
    .expert_level=2
    .short_caption = Ignore Xray/neutron R-free flags set mismatch
  %s
}

"""%data_and_flags_str

neutron_data_str_no_filenames = """\
neutron_data
  .help=Scope of neutron data and neutron free-R flags
  .style = scrolled auto_align
{
  ignore_xn_free_r_mismatch = False
    .type = bool
    .expert_level=2
    .short_caption = Ignore Xray/neutron R-free flags set mismatch
  %s
}

"""%data_and_flags_str_no_filenames

def data_and_flags_master_params(master_scope_name=None):
  if(master_scope_name is not None):
    p = """\
%s
{
%s
}
"""
    return iotbx.phil.parse(p%(master_scope_name, data_and_flags_str), process_includes=True)
  else:
    return iotbx.phil.parse(data_and_flags_str, process_includes=True)

class run(object):
  """
  Encapsulates logic for extracting experimental amplitudes and R-free flags
  from the given input file(s) via using reflection_file_server. This expects
  that the standard parameter block is being used.  Determination of appropriate
  data labels will be as automatic as possible, or will give clear feedback when
  ambiguity exists. If not found in the inputs, the R-free flags can be created
  if desired.

  It is supposed to run silently.
  Error messages are accumulated in self.err, and can be rised as desired after
  the execution. Log info messages are stored in self.log and can be flushed
  after the execution as well.
  """
  def __init__(self,
               reflection_file_server,
               parameters = None,
               experimental_phases_params = None,
               working_point_group = None,
               remark_r_free_flags_md5_hexdigest = None,
               extract_r_free_flags = True,
               extract_experimental_phases = True,
               keep_going = False,
               prefer_anomalous = None,
               force_non_anomalous = False,
               allow_mismatch_flags = False,
               free_r_flags_scope = 'miller_array.labels.name',
               ):
    adopt_init_args(self, locals())
    # Buffers for error and log messages.
    self.err = []
    self.log = StringIO() # container for all log messages
    #
    if(self.parameters is None):
      self.parameters = data_and_flags_master_params().extract()
    self.r_free_flags               = None
    self.test_flag_value            = None
    self.r_free_flags_md5_hexdigest = None
    self.experimental_phases        = None
    self.raw_experimental_phases    = None
    # Get data first
    self.raw_data            = self._extract_data()
    # Apply resolution and sigma cutoffs, if requested
    self.raw_data_truncated  = self._apply_cutoffs()
    self.raw_flags_truncated = None
    # Convert to usable Fobs
    self.f_obs               = self._data_as_f_obs()
    # Then extract or generate flags
    self.raw_flags = None
    if(extract_r_free_flags):
      self.raw_flags = self._extract_flags()
      if(self.raw_flags is not None):
        flags_info = self.raw_flags.info()
        if(self.raw_flags.anomalous_flag() and
           not self.raw_data_truncated.anomalous_flag()):
          self.raw_flags = self.raw_flags.as_non_anomalous_array()
        self.raw_flags_truncated = self.raw_flags.common_set(
          self.raw_data_truncated)
    if(extract_r_free_flags and self.raw_flags is not None):
      self._get_r_free_flags()
      if(self.r_free_flags is not None):
        self.r_free_flags.set_info(flags_info)
    # Make sure they match
    if(self.r_free_flags is not None):
      f_obs_info = self.f_obs.info()
      flags_info = self.r_free_flags.info()
      self.f_obs, self.r_free_flags = self.f_obs.common_sets(self.r_free_flags)
      self.f_obs        = self.f_obs.set_info(f_obs_info)
      self.r_free_flags = self.r_free_flags.set_info(flags_info)
    # extract phases
    if extract_experimental_phases:
      self.experimental_phases = self._determine_experimental_phases(
        parameters      = experimental_phases_params,
        parameter_scope = 'miller_array.labels.name')
    # Fill in log
    self._show_summary()

  def _compose_mtz_object(self):
    f_obs_label = "F-obs"
    i_obs_label = "I-obs"
    flag_label  = "R-free-flags"
    phase_label = "HL"
    if(self.raw_data.is_xray_intensity_array()):
      column_root_label = i_obs_label
    else:
      column_root_label = f_obs_label
    mtz_dataset = self.raw_data.as_mtz_dataset(
      column_root_label = column_root_label)
    if (self.raw_flags is not None):
      mtz_dataset.add_miller_array(
        miller_array      = self.raw_flags,
        column_root_label = flag_label)
    if(self.experimental_phases is not None):
      mtz_dataset.add_miller_array(
        miller_array      = self.raw_experimental_phases,
        column_root_label = phase_label)
    labels = group_args(
      data   = column_root_label,
      flags  = flag_label,
      phases = phase_label)
    return group_args(
      mtz_dataset = mtz_dataset,
      labels      = labels)

  def result(self):
    """
    Container for:
      - raw arrays as extracted from inputs.
      - arrays that are ready to use (after aply cutoffs, map-to-asu, etc).
      - error messages (to rise if desired).
      - log messaged (to flush if desired).
      - misc.

    """
    i_obs = None
    if(self.raw_data.is_xray_intensity_array()):
      i_obs = self.raw_data
    o = group_args(
      raw_data                   = self.raw_data,
      raw_flags                  = self.raw_flags,
      raw_experimental_phases    = self.raw_experimental_phases,
      f_obs                      = self.f_obs,
      i_obs                      = i_obs,
      r_free_flags               = self.r_free_flags,
      test_flag_value            = self.test_flag_value,
      experimental_phases        = self.experimental_phases,
      r_free_flags_md5_hexdigest = self.r_free_flags_md5_hexdigest,
      err                        = self.err,
      log                        = self.log,
      mtz_object                 = self._compose_mtz_object())
    o.stop_dynamic_attributes()
    return o

  def show_summary(self, log, prefix=""):
    print(self.log.getvalue(), file=log)

  def _print(self, m): print(m, file=self.log)

  def _show_summary(self, prefix=""):
    log = self.log
    self._print("%sInput data summary:"%prefix)
    self.raw_data.show_comprehensive_summary(f=log, prefix="  "+prefix)
    if(not self.raw_data.indices().all_eq(self.f_obs.indices()) or
       type(self.raw_data.observation_type()) !=
       type(self.f_obs.observation_type())):
      print(prefix, file=log)
      print("%sData used summary:"%prefix, file=log)
      self.f_obs.show_comprehensive_summary(f=log, prefix="  "+prefix)
    if(self.r_free_flags is not None):
      print(prefix, file=log)
      self._print("%sFree-R flags summary:"%prefix)
      self.r_free_flags.show_comprehensive_summary(f=log, prefix="  "+prefix)
      print("%s  Test (R-free) flag value: %d"%(prefix, self.test_flag_value),
        file=self.log)
      self.r_free_flags.show_r_free_flags_info(out=log, prefix="  "+prefix)
    if(self.experimental_phases is not None):
      self._print("%sExperimental phases:"%prefix)
      self.experimental_phases.show_comprehensive_summary(f=log, prefix="  "+prefix)
      print("%s  Average figures of merit by resolution:"%prefix, file=log)
      figures_of_merit = abs(self.experimental_phases.phase_integrals())
      figures_of_merit.setup_binner(n_bins=10)
      legend_len = figures_of_merit.mean(use_binning=True).show(
        data_fmt="%6.3f", show_unused=False, f = log, prefix="    "+prefix)
      print("   ", ("%%%ds"%legend_len)%"overall", \
        "%6.3f"%figures_of_merit.mean(), file=log)
      print("%s  Note: Figures of merit are determined by integration of"%prefix, file=log)
      print("%s        Hendrickson-Lattman coefficients."%prefix, file=log)
      print(file=log)
      print("%s  Number and fraction of available experimental phases by resolution:"%prefix, file=log)
      self.experimental_phases.setup_binner(n_bins=10)
      self.experimental_phases.count_and_fraction_in_bins(
        data_value_to_count = (0,0,0,0),
        count_not_equal= True).show(show_unused = False, f = log, prefix="    "+prefix)
      print(file=log)

  def _determine_experimental_phases(
        self,
        parameters=None,
        parameter_scope=None,
        ignore_all_zeros = True):
    try:
      file_name = None
      labels = None
      if(parameters is not None):
        file_name = parameters.file_name
        labels    = parameters.labels
      experimental_phases = \
        self.reflection_file_server.get_experimental_phases(
          file_name        = file_name,
          labels           = labels,
          raise_no_array   = False,
          parameter_name   = "",
          ignore_all_zeros = ignore_all_zeros,
          parameter_scope  = parameter_scope)

      if(experimental_phases is None): return None
    except reflection_file_utils.Sorry_No_array_of_the_required_type:
      experimental_phases = None
    else:
      info = experimental_phases.info()
      if(parameters is not None):
        parameters.file_name = experimental_phases.info().source
        parameters.labels = [experimental_phases.info().label_string()]
      msg = miller_array_symmetry_safety_check(
        miller_array        = experimental_phases,
        data_description    = "Experimental phases",
        working_point_group = self.working_point_group)
      if(msg is not None and not self.keep_going):
        self.err.append(msg)
      experimental_phases = experimental_phases.regularize()
      self.raw_experimental_phases = experimental_phases
      if(not self.f_obs.anomalous_flag()):
        if(experimental_phases.anomalous_flag()):
          experimental_phases = experimental_phases.average_bijvoet_mates()
      elif(not experimental_phases.anomalous_flag()):
         experimental_phases = experimental_phases.generate_bijvoet_mates()
      self.experimental_phases = experimental_phases.matching_set(
        other = self.f_obs, data_substitute=(0,0,0,0))
      experimental_phases.set_info(info)
    return experimental_phases

  def _get_r_free_flags(self):
    self.r_free_flags,self.test_flag_value,self.r_free_flags_md5_hexdigest =\
      self._flags_as_r_free_flags(f_obs = self.f_obs, r_free_flags =
      self.raw_flags)
    if(self.r_free_flags is not None):
      self.r_free_flags.set_info(self.raw_flags.info())

  def _extract_data(self):
    data = self.reflection_file_server.get_xray_data(
      file_name        = self.parameters.file_name,
      labels           = self.parameters.labels,
      ignore_all_zeros = self.parameters.ignore_all_zeros,
      parameter_name   = "",
      parameter_scope  = "miller_array.labels.name",
      prefer_anomalous = self.prefer_anomalous)
    self.parameters.file_name = data.info().source
    self.parameters.labels    = [data.info().label_string()]
    if(self.working_point_group is not None):
      msg = miller_array_symmetry_safety_check(
        miller_array        = data,
        data_description    = "Reflection data",
        working_point_group = self.working_point_group)
      if(msg is not None and not self.keep_going):
        self.err.append(msg)
    return data.regularize()

  def _extract_flags(self, data_description = "R-free flags"):
    r_free_flags, test_flag_value = None, None
    params = self.parameters.r_free_flags
    # Extract
    if(not self.parameters.r_free_flags.generate):
      try:
        r_free_flags, test_flag_value = \
          self.reflection_file_server.get_r_free_flags(
            file_name                = params.file_name,
            label                    = params.label,
            test_flag_value          = params.test_flag_value,
            disable_suitability_test = params.disable_suitability_test,
            parameter_scope          = self.free_r_flags_scope)
      except reflection_file_utils.Sorry_No_array_of_the_required_type as e:
        if(self.parameters.r_free_flags.generate is not None):
          if(not self.keep_going):
            self.err.append(explain_how_to_generate_array_of_r_free_flags(
              scope = "xray_data.r_free_flags.generate"))
            self.err.append("Please try again.")
          return None
        r_free_flags, test_flag_value = None, None
      else:
        params.file_name       = r_free_flags.info().source
        params.label           = r_free_flags.info().label_string()
        params.test_flag_value = test_flag_value
        msg = miller_array_symmetry_safety_check(
          miller_array        = r_free_flags,
          data_description    = data_description,
          working_point_group = self.working_point_group)
        if(msg is not None and not self.keep_going):
          self.err.append(msg)
        info = r_free_flags.info()
        try:
          processed = r_free_flags.regularize()
        except RuntimeError as e:
          self.err.append("Bad free-r flags:\n %s"%str(e))
          return None
        if (self.force_non_anomalous):
          processed = processed.average_bijvoet_mates()
        r_free_flags = processed.set_info(info)
    # Generate or stop
    if(r_free_flags is None):
      if ((params.fraction is None) or
          (params.lattice_symmetry_max_delta is None) or
          (params.use_lattice_symmetry is None)):
        msg = """
No R-free flags are available, but one or more parameters required to generate
new flags is undefined.
"""
        self.err.append(msg)
        return None
      print("Generating a new array of R-free flags.", file=self.log)
      print(file=self.log)
      r_free_flags = self.f_obs.generate_r_free_flags(
        fraction                   = params.fraction,
        max_free                   = params.max_free,
        lattice_symmetry_max_delta = params.lattice_symmetry_max_delta,
        use_lattice_symmetry       = params.use_lattice_symmetry,
        use_dataman_shells         = params.use_dataman_shells,
        n_shells                   = params.n_shells
        ).set_info(miller.array_info(labels = ["R-free-flags"]))
      params.label           = r_free_flags.info().label_string()
      params.test_flag_value = 1
    # check if anomalous pairs are sound
    if(r_free_flags is not None):
      r_free_flags.deep_copy().as_non_anomalous_array()
    # make sure flags match anomalous flag of data
    if(self.raw_data.anomalous_flag() and not r_free_flags.anomalous_flag()):
      info = r_free_flags.info()
      observation_type = r_free_flags.observation_type()
      r_free_flags = r_free_flags.generate_bijvoet_mates()
      r_free_flags.set_observation_type(observation_type)
      r_free_flags.set_info(info)
    return r_free_flags

  def _apply_cutoffs(self):
    if(self.raw_data is None): return None
    selection = self.raw_data.all_selection()
    dd = self.raw_data.d_spacings().data()
    if(self.parameters.low_resolution is not None):
      selection &= dd <= self.parameters.low_resolution
    if(self.parameters.high_resolution is not None):
      selection &= dd >= self.parameters.high_resolution
    if(self.raw_data.sigmas() is not None):
      sigma_cutoff = self.parameters.sigma_fobs_rejection_criterion
      if(sigma_cutoff is not None and sigma_cutoff > 0):
        s = self.raw_data.data() > self.raw_data.sigmas()*sigma_cutoff
        selection &= s
    if(selection.count(True) == 0):
      self.err.append("No data left after applying resolution and sigma cutoffs.")
      return None
    return self.raw_data.select(selection)

  def _data_as_f_obs(self):
    """
    Convert input data array to amplitudes, adjusting the data type and
    applying additional filters if necessary.

    :param f_obs: selected input data
    :returns: :py:class:`cctbx.miller.array` of real numbers with observation
      type set to amplitudes
    """
    if(self.raw_data_truncated is None): return None
    f_obs = self.raw_data_truncated.deep_copy()
    # Convert to non-anomalous if requested
    if(self.force_non_anomalous):
      f_obs = f_obs.average_bijvoet_mates()
    #
    d_min = f_obs.d_min()
    if(d_min < 0.25):
      self.err.append("Resolution of data is too high: %-6.4f A"%d_min)
      return None
    if(f_obs.is_complex_array()): f_obs = abs(f_obs)
    if(f_obs.is_xray_intensity_array()):
      if(self.parameters.french_wilson_scale):
        try:
          sigI_cutoff = self.parameters.sigma_iobs_rejection_criterion
          f_obs = french_wilson.french_wilson_scale(
            miller_array                   = f_obs,
            params                         = self.parameters.french_wilson,
            sigma_iobs_rejection_criterion = sigI_cutoff,
            log                            = self.log)
        except Exception as e:
          print(str(e), file=self.log)
          print("Using alternative Iobs->Fobs conversion.", file=self.log)
          f_obs = f_obs.f_sq_as_f()
    #
    f_obs.set_observation_type_xray_amplitude()
    if(self.parameters.force_anomalous_flag_to_be_equal_to is not None):
      if(not self.parameters.force_anomalous_flag_to_be_equal_to):
        if(f_obs.anomalous_flag()):
          merged = f_obs.as_non_anomalous_array().merge_equivalents()
          f_obs = merged.array().set_observation_type( f_obs )
      elif(not f_obs.anomalous_flag()):
        observation_type = f_obs.observation_type()
        f_obs = f_obs.generate_bijvoet_mates()
        f_obs.set_observation_type(observation_type)
    else:
      f_obs = f_obs.convert_to_non_anomalous_if_ratio_pairs_lone_less_than(
        threshold=self.parameters.
          convert_to_non_anomalous_if_ratio_pairs_lone_less_than_threshold)
    f_obs.set_info(self.raw_data.info())
    return f_obs

  def _flags_as_r_free_flags(self,
        f_obs,
        r_free_flags,
        missing_show_max_lines=10):
    test_flag_value = self.parameters.r_free_flags.test_flag_value
    if(test_flag_value is None):
      msg = """
Could not determine an appropriate test flag for the data with label(s) '%s'.
This may happen if they are all a single value; please check the file to make
sure the flags are suitable for use.
"""
      self.err.append(msg % self.parameters.r_free_flags.label)
    if (isinstance(r_free_flags.data(), flex.bool)):
      r_free_flags = r_free_flags.array(
        data = r_free_flags.data() == bool(test_flag_value))
    else:
      r_free_flags = r_free_flags.array(
        data = r_free_flags.data() == test_flag_value)
    r_free_flags_md5_hexdigest = \
      r_free_flags.map_to_asu().sort(by_value="packed_indices").data() \
        .md5().hexdigest()
    if(self.remark_r_free_flags_md5_hexdigest is not None):
      self._verify_r_free_flags_md5_hexdigest(
        ignore_pdb_hexdigest = self.parameters.r_free_flags.ignore_pdb_hexdigest,
        current              = r_free_flags_md5_hexdigest,
        records              = self.remark_r_free_flags_md5_hexdigest)
    if(not f_obs.anomalous_flag()):
      if(r_free_flags.anomalous_flag()):
        print("Reducing R-free flags to non-anomalous array.", file=self.log)
        r_free_flags = r_free_flags.average_bijvoet_mates()
        print(file=self.log)
    elif(not r_free_flags.anomalous_flag()):
       print("Generating Bijvoet mates of R-free flags.", file=self.log)
       r_free_flags = r_free_flags.generate_bijvoet_mates()
       print(file=self.log)
    r_free_flags = r_free_flags.map_to_asu().common_set(f_obs)
    n_missing_r_free_flags = f_obs.indices().size() \
      - r_free_flags.indices().size()
    if(not self.allow_mismatch_flags and n_missing_r_free_flags != 0):
      msg = [
        "R-free flags not compatible with F-obs array:"
        " missing flag for %d F-obs selected for refinement"
          % n_missing_r_free_flags]
      if (missing_show_max_lines is not None and missing_show_max_lines <= 0):
        msg[0] += "."
      else:
        msg[0] += ":"
        lone = f_obs.lone_set(other=r_free_flags)
        if (missing_show_max_lines is None):
          n_not_shown = 0
        else:
          n_not_shown = lone.indices().size() - missing_show_max_lines
          if (n_not_shown > missing_show_max_lines * 0.5):
            lone = lone[:missing_show_max_lines]
          else:
            n_not_shown = 0
        if (lone.sigmas() is None):
          msg.append("    h   k   l   data")
          for hkl,f in zip(lone.indices(), lone.data()):
            msg.append("  %3d %3d %3d" % hkl + "   %.6g" % f)
        else:
          msg.append("    h   k   l   data  sigmas")
          for hkl,f,s in zip(lone.indices(), lone.data(), lone.sigmas()):
            msg.append("  %3d %3d %3d" % hkl + "   %.6g  %.6g" % (f,s))
        if (n_not_shown != 0):
          msg.append("    ... (remaining %d not shown)" % n_not_shown)
      self.err.append("\n".join(msg))
      return None,None,None
    return r_free_flags, test_flag_value, r_free_flags_md5_hexdigest

  def _verify_r_free_flags_md5_hexdigest(self,
        ignore_pdb_hexdigest,
        current,
        records):
    from_file = set()
    for record in records:
      flds = record.split()
      if (len(flds) == 3):
        from_file.add(flds[2])
    if(len(from_file) > 1):
      msg="""
Multiple conflicting REMARK r_free_flags.md5.hexdigest records found in the
input PDB file.
"""
      self.err.append(msg)
    if (len(from_file) == 1 and current not in from_file):
      log = self.log
      for i in range(2): print("*"*79, file=log)
      if(ignore_pdb_hexdigest):
        print(file=log)
        print(" ".join(["WARNING"]*9), file=log)
      print("""
The MD5 checksum for the R-free flags array summarized above is:
  %s

The corresponding MD5 checksum in the PDB file summarized above is:
  %s

These checksums should be identical but are in fact different. This is
because the R-free flags used at previous stages of refinement are
different from the R-free flags summarized above. As a consequence,
the values for R-free could be biased and misleading.

However, there is no problem if the R-free flags were just extended to
a higher resolution, or if some reflections with no data or that are
not part of the R-free set have been added or removed.""" % (
  current, sorted(from_file)[0]), end=' ', file=log)
      if (not ignore_pdb_hexdigest):
        print("""\
In this case,
simply remove the

  REMARK r_free_flags.md5.hexdigest %s

record from the input PDB file to proceed with the refinement.""" % (
  sorted(from_file)[0]), end=' ', file=log)
      print("""

Otherwise it is best to recover the previously used R-free flags
and use them consistently throughout the refinement of the model.
Run this command again with the name of the file containing the
original flags as an additional input.
""", file=log)
      if (not ignore_pdb_hexdigest):
        print("""\
If the original R-free flags are unrecoverable, remove the REMARK
record as indicated above. In this case the values for R-free will
become meaningful only after many cycles of refinement.
""", file=log)
      else:
        print("""\
If the original R-free flags are unrecoverable, the values for R-free
will become meaningful only after many cycles of refinement.
""", file=log)
      for i in range(2): print("*"*79, file=log)
      print(file=log)
      if (not ignore_pdb_hexdigest):
        if ("PHENIX_GUI_ENVIRONMENT" in os.environ):
          log.flush()
          raise Sorry("This model appears to have previously been refined "+
            "against a different set of R-free flags.  Please resolve the "+
            "mismatch; additional information and instructions are available "+
            "at the end of the log output.")
        else :
          raise Sorry("Please resolve the R-free flags mismatch.")

map_coefficients_params_str = """\
  file_name=None
    .type=path
    .short_caption=Map coefficients file
  labels=None
    .type=strings
"""

experimental_phases_params_str = """\
  file_name=None
    .type=path
    .short_caption=Experimental phase file
    .style = bold input_file file_type:hkl process_hkl child:hl_coeffs:labels
  labels=None
    .type=strings
    .input_size = 160
    .short_caption = Phase labels
    .style = renderer:draw_hl_label_widget bold
"""

experimental_phases_params = iotbx.phil.parse(
  input_string=experimental_phases_params_str)

def experimental_phases_master_params():
  return experimental_phases_params
