###
# XXX Severe duplicaiton. Remove equivalent code from mmtbx.utils.
###
from __future__ import division
from __future__ import print_function
import iotbx.phil
from libtbx import adopt_init_args
from cctbx.array_family import flex

def miller_array_symmetry_safety_check(miller_array,
                                       data_description,
                                       working_point_group,
                                       symmetry_safety_check,
                                       log):
  msg = miller_array.crystal_symmetry_is_compatible_with_symmetry_from_file(
    working_point_group = working_point_group).format_error_message(
      data_description = data_description)
  if(msg is not None):
     if(symmetry_safety_check == "warning"):
        print("*" * 79, file=log)
        print("WARNING:", msg, file=log)
        print("*" * 79, file=log)
     else:
        raise Sorry(msg + """
  The program inspects all inputs to determine the working crystal
  symmetry (unit cell & space group).
  Please check the working crystal symmetry shown above. If it is
  not correct, use the --unit_cell, --space_group, or --symmetry
  option to specify the correct unit cell parameters and space group
  symbol.
  If the working crystal symmetry is in fact correct, disable this
  error by adding
    refinement.input.symmetry_safety_check=warning
  to the command line arguments.
""")

def explain_how_to_generate_array_of_r_free_flags(log, flags_parameter_scope):
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
  print(part1 + flags_parameter_scope+""".generate=True""" + part3, file=log)

master_params_str="""
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
               If no test set is present in the reflections file, one can be \
               generated automatically, or you can use the reflection file \
               editor to combine an existing set with your X-ray or neutron data.
  {
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
        fraction = 0.1
      .type=float
      .short_caption = Fraction of reflections in test set
      .expert_level=0
    max_free = 2000
      .type=int
      .short_caption = Maximum number of reflections in test set
      .expert_level=2
    lattice_symmetry_max_delta = 5
      .type=float
      .expert_level=2
    use_lattice_symmetry = True
      .type=bool
      .short_caption = Use lattice symmetry to generate test set
      .expert_level=0
    use_dataman_shells = False
      .type = bool
      .short_caption = Assign test set in thin resolution shells
      .help = Used to avoid biasing of the test set by certain types of \
              non-crystallographic symmetry.
    n_shells = 20
      .type = int
      .short_caption = Number of resolution shells
  }
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

class run(object):
  """
  Encapsulates logic for extracting experimental amplitudes and R-free flags
  from the given input file(s).  This expects that the standard parameter block
  is being used.  Determination of appropriate data labels will be as automatic
  as possible, or will give clear feedback when ambiguity exists.  If not
  found in the inputs, the R-free flags can be created if desired.
  """
  def __init__(self, reflection_file_server,
                     parameters = None,
                     data_parameter_scope = "",
                     flags_parameter_scope = "",
                     data_description = None,
                     working_point_group = None,
                     symmetry_safety_check = None,
                     remark_r_free_flags_md5_hexdigest = None,
                     extract_r_free_flags = True,
                     keep_going = False,
                     log = None,
                     prefer_anomalous = None,
                     force_non_anomalous = False):
    adopt_init_args(self, locals())
    if(self.parameters is None):
      self.parameters = master_params().extract()
    self.intensity_flag = False
    self.f_obs = None
    self.r_free_flags = None
    self.test_flag_value = None
    self.r_free_flags_md5_hexdigest = None
    if(data_description is not None):
      print_statistics.make_header(data_description, out = log)
    self.raw_data = self.extract_data()
    data_info = self.raw_data.info()
    self.f_obs = self.data_as_f_obs(f_obs = self.raw_data)
    self.f_obs.set_info(data_info)
    if(extract_r_free_flags):
      self.raw_flags = self.extract_flags(data = self.raw_data)
      if(self.raw_flags is not None):
        flags_info = self.raw_flags.info()
    if(extract_r_free_flags and self.raw_flags is not None):
      self.get_r_free_flags()
      self.r_free_flags.set_info(flags_info)

  def get_r_free_flags(self):
    self.r_free_flags,self.test_flag_value,self.r_free_flags_md5_hexdigest =\
      self.flags_as_r_free_flags(f_obs = self.f_obs, r_free_flags =
      self.raw_flags)
    self.r_free_flags.set_info(self.raw_flags.info())

  def extract_data(self):
    data = self.reflection_file_server.get_xray_data(
      file_name        = None,
      labels           = None,
      ignore_all_zeros = self.parameters.ignore_all_zeros,
      parameter_scope  = self.data_parameter_scope,
      prefer_anomalous = self.prefer_anomalous)
    #self.parameters.file_name = data.info().source
    #self.parameters.labels = [data.info().label_string()]
    if(data.is_xray_intensity_array()):
      print("I-obs:", file=self.log)
      self.intensity_flag = True
    else:
      print("F-obs:", file=self.log)
    print(" ", data.info(), file=self.log)
    if([self.data_description, self.working_point_group,
       self.symmetry_safety_check].count(None) == 0):
      miller_array_symmetry_safety_check(
        miller_array          = data,
        data_description      = self.data_description,
        working_point_group   = self.working_point_group,
        symmetry_safety_check = self.symmetry_safety_check,
        log                   = self.log)
      print(file=self.log)
    info = data.info()
    processed = data.eliminate_sys_absent(log = self.log)
    if(processed is not data):
      info = info.customized_copy(systematic_absences_eliminated = True)
    if(not processed.is_unique_set_under_symmetry()):
      if(data.is_xray_intensity_array()):
        print("Merging symmetry-equivalent intensities:", file=self.log)
      else:
        print("Merging symmetry-equivalent amplitudes:", file=self.log)
      merged = processed.merge_equivalents()
      merged.show_summary(out = self.log, prefix="  ")
      print(file=self.log)
      processed = merged.array()
      info = info.customized_copy(merged=True)
    if (self.force_non_anomalous):
      processed = processed.average_bijvoet_mates()
    return processed.set_info(info)

  def extract_flags(self, data, data_description = "R-free flags"):
    r_free_flags, test_flag_value = None, None
    params = self.parameters.r_free_flags
    if(not self.parameters.r_free_flags.generate):
      try:
        r_free_flags, test_flag_value = \
          self.reflection_file_server.get_r_free_flags(
            file_name                = None,
            label                    = None,
            test_flag_value          = params.test_flag_value,
            disable_suitability_test = params.disable_suitability_test,
            parameter_scope          = self.flags_parameter_scope)
      except reflection_file_utils.Sorry_No_array_of_the_required_type, e:
        if(self.parameters.r_free_flags.generate is not None):
          explain_how_to_generate_array_of_r_free_flags(log = self.log,
            flags_parameter_scope = self.flags_parameter_scope)
          if(self.keep_going): return None
          raise Sorry("Please try again.")
        r_free_flags, test_flag_value = None, None
      else:
        #params.file_name = r_free_flags.info().source
        #params.label = r_free_flags.info().label_string()
        params.test_flag_value = test_flag_value
        print(data_description+":", file=self.log)
        print(" ", r_free_flags.info(), file=self.log)
        if([self.working_point_group,
           self.symmetry_safety_check].count(None) == 0):
          miller_array_symmetry_safety_check(
            miller_array          = r_free_flags,
            data_description      = data_description,
            working_point_group   = self.working_point_group,
            symmetry_safety_check = self.symmetry_safety_check,
            log                   = self.log)
          print(file=self.log)
        info = r_free_flags.info()
        processed = r_free_flags.eliminate_sys_absent(log = self.log)
        if(processed is not r_free_flags):
          info = info.customized_copy(systematic_absences_eliminated = True)
        if(not processed.is_unique_set_under_symmetry()):
           print("Checking symmetry-equivalent R-free flags for consistency:", end=' ', file=self.log)
           try:
             merged = processed.merge_equivalents()
           except RuntimeError, e:
             print(file=self.log)
             error_message = str(e)
             expected_error_message = "cctbx Error: merge_equivalents_exact: "
             assert error_message.startswith(expected_error_message)
             raise Sorry("Incompatible symmetry-equivalent R-free flags: %s" %
               error_message[len(expected_error_message):])
           else:
             print("OK", file=self.log)
             print(file=self.log)
           processed = merged.array()
           info = info.customized_copy(merged=True)
           del merged
        if (self.force_non_anomalous):
          processed = processed.average_bijvoet_mates()
        r_free_flags = processed.set_info(info)
    if(r_free_flags is None):
      if ((params.fraction is None) or
          (params.lattice_symmetry_max_delta is None) or
          (params.use_lattice_symmetry is None)):
        raise Sorry("No R-free flags are available, but one or more "+
          "parameters required to generate new flags is undefined.")
      print("Generating a new array of R-free flags.", file=self.log)
      print(file=self.log)
      libtbx.call_back(message="warn",
        data="PHENIX will generate a new array of R-free flags.  Please "+
          "check to make sure that the input data do not already contain "+
          "an R-free set; if one is present, you should cancel this job and "+
          "disable generation of new flags.  If the program you are running "+
          "outputs an MTZ file, you should be sure to use that file in all "+
          "future refinements.")
      r_free_flags = data.generate_r_free_flags(
        fraction                   = params.fraction,
        max_free                   = params.max_free,
        lattice_symmetry_max_delta = params.lattice_symmetry_max_delta,
        use_lattice_symmetry       = params.use_lattice_symmetry,
        use_dataman_shells         = params.use_dataman_shells,
        n_shells                   = params.n_shells
        ).set_info(miller.array_info(labels = ["R-free-flags"]))
      params.label = r_free_flags.info().label_string()
      params.test_flag_value = 1
    # check if anomalous pairs are sound
    if(r_free_flags is not None):
      r_free_flags.deep_copy().as_non_anomalous_array()
    return r_free_flags

  def data_as_f_obs(self, f_obs):
    """
    Convert input data array to amplitudes, adjusting the data type and
    applying additional filters if necessary.

    :param f_obs: selected input data
    :returns: :py:class:`cctbx.miller.array` of real numbers with observation
      type set to amplitudes
    """
    if(not f_obs.sigmas_are_sensible()):
      f_obs = f_obs.customized_copy(
        indices=f_obs.indices(),
        data=f_obs.data(),
        sigmas=None).set_observation_type(f_obs)
    # Delete F(0,0,0) if present
    sel = f_obs.indices()==(0,0,0)
    if(sel.count(True)>0):
      print("F(0,0,0) will be removed.", file=self.log)
      f_obs = f_obs.select(~sel)
    #
    d_min = f_obs.d_min()
    if(d_min < 0.25): # XXX what is the equivalent for neutrons ???
      raise Sorry("Resolution of data is too high: %-6.4f A"%d_min)
    f_obs.show_comprehensive_summary(f = self.log)
    f_obs_data_size = f_obs.data().size()
    print(file=self.log)
    if(f_obs.is_complex_array()): f_obs = abs(f_obs)
    f_obs_fw = None
    if(f_obs.is_xray_intensity_array()):
      if(self.parameters.french_wilson_scale):
        f_obs_fw = french_wilson.french_wilson_scale(
          miller_array=f_obs,
          params=self.parameters.french_wilson,
          sigma_iobs_rejection_criterion=\
            self.parameters.sigma_iobs_rejection_criterion,
          log=self.log)
        if f_obs_fw is not None:
          f_obs = f_obs_fw
      if (not self.parameters.french_wilson_scale or f_obs_fw is None):
        selection_by_isigma = self._apply_sigma_cutoff(
          f_obs   = f_obs,
          n       = self.parameters.sigma_iobs_rejection_criterion,
          message = "Number of reflections with |Iobs|/sigma(Iobs) < %5.2f: %d")
        if(selection_by_isigma is not None):
          f_obs = f_obs.select(selection_by_isigma)
        f_obs = f_obs.f_sq_as_f()
      print("Intensities converted to amplitudes for use in refinement.", file=self.log)
      print(file=self.log)
    #
    sigmas = f_obs.sigmas()
    if(sigmas is not None):
      selection  = sigmas > 0
      selection &= f_obs.data()>=0
      n_both_zero = selection.count(False)
      if(n_both_zero>0):
        print("Number of pairs (Fobs,sigma)=(0,0) is %s. They will be removed"%\
          n_both_zero, file=self.log)
        f_obs = f_obs.select(selection)
    #
    f_obs.set_observation_type_xray_amplitude()
    f_obs = f_obs.map_to_asu()
    selection = f_obs.all_selection()
    if(self.parameters.low_resolution is not None):
      selection &= f_obs.d_spacings().data() <= self.parameters.low_resolution
    if(self.parameters.high_resolution is not None):
      selection &= f_obs.d_spacings().data() >= self.parameters.high_resolution
    selection_positive = f_obs.data() >= 0
    print("Number of F-obs in resolution range:                  ", \
      selection.count(True), file=self.log)
    print("Number of F-obs<0 (these reflections will be rejected):", \
      selection_positive.count(False), file=self.log)
    selection_zero = f_obs.data() == 0
    print("Number of F-obs=0 (these reflections will be used in refinement):", \
      selection_zero.count(True), file=self.log)
    selection &= selection_positive
    selection_by_fsigma = self._apply_sigma_cutoff(
      f_obs   = f_obs,
      n       = self.parameters.sigma_fobs_rejection_criterion,
      message = "Number of reflections with |Fobs|/sigma(Fobs) < %5.2f: %d")
    if(selection_by_fsigma is not None): selection &= selection_by_fsigma
    selection &= f_obs.d_star_sq().data() > 0
    f_obs = f_obs.select(selection)
    rr = f_obs.resolution_range()
    print("Refinement resolution range: d_max = %8.4f" % rr[0], file=self.log)
    print("                             d_min = %8.4f" % rr[1], file=self.log)
    print(file=self.log)
    if(f_obs.indices().size() == 0):
      raise Sorry(
        "No data left after applying resolution limits and sigma cutoff.")
    if(self.parameters.force_anomalous_flag_to_be_equal_to is not None):
      if(not self.parameters.force_anomalous_flag_to_be_equal_to):
        print("force_anomalous_flag_to_be_equal_to=False", file=self.log)
        if(f_obs.anomalous_flag()):
          print("Reducing data to non-anomalous array.", file=self.log)
          merged = f_obs.as_non_anomalous_array().merge_equivalents()
          merged.show_summary(out = self.log, prefix="  ")
          f_obs = merged.array().set_observation_type( f_obs )
          del merged
          print(file=self.log)
      elif(not f_obs.anomalous_flag()):
        print("force_anomalous_flag_to_be_equal_to=True", file=self.log)
        print("Generating Bijvoet mates of X-ray data.", file=self.log)
        observation_type = f_obs.observation_type()
        f_obs = f_obs.generate_bijvoet_mates()
        f_obs.set_observation_type(observation_type)
        print(file=self.log)
    else:
      f_obs = f_obs.convert_to_non_anomalous_if_ratio_pairs_lone_less_than(
        threshold=self.parameters.
          convert_to_non_anomalous_if_ratio_pairs_lone_less_than_threshold)
    if(f_obs_data_size != f_obs.data().size()):
      print("\nFobs statistics after all cutoffs applied:\n", file=self.log)
      f_obs.show_comprehensive_summary(f = self.log)
    return f_obs

  def _apply_sigma_cutoff(self, f_obs, n, message):
    selection = None
    if(f_obs.sigmas() is not None):
      sigma_cutoff = n
      if(sigma_cutoff is not None and sigma_cutoff > 0):
        selection_by_sigma = f_obs.data() > f_obs.sigmas()*sigma_cutoff
        print(message % (sigma_cutoff,
          selection_by_sigma.count(False)), file=self.log)
        selection = selection_by_sigma
    return selection

  def flags_as_r_free_flags(self,
        f_obs,
        r_free_flags,
        missing_show_max_lines=10):
    test_flag_value = self.parameters.r_free_flags.test_flag_value
    if (test_flag_value is None):
      raise Sorry(("PHENIX could not determine an appropriate test flag "+
        "for the data with label(s) '%s'.  This may happen if they are all "+
        "a single value; please check the file to make sure the flags are "+
        "suitable for use.") % self.parameters.r_free_flags.label)
    r_free_flags.show_comprehensive_summary(f = self.log)
    print(file=self.log)
    print("Test (R-free flags) flag value:", test_flag_value, file=self.log)
    print(file=self.log)
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
      self.verify_r_free_flags_md5_hexdigest(
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
    if(n_missing_r_free_flags != 0):
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
      raise Sorry("\n".join(msg))
    r_free_flags.show_r_free_flags_info(out = self.log, prefix="")
    return r_free_flags, test_flag_value, r_free_flags_md5_hexdigest

  def verify_r_free_flags_md5_hexdigest(self,
        ignore_pdb_hexdigest,
        current,
        records):
    from_file = set()
    for record in records:
      flds = record.split()
      if (len(flds) == 3):
        from_file.add(flds[2])
    if (len(from_file) > 1):
      raise Sorry(
        "Multiple conflicting REMARK r_free_flags.md5.hexdigest records"
        " found in the input PDB file.")
    if (len(from_file) == 1 and current not in from_file):
      log = self.log
      for i in xrange(2): print("*"*79, file=log)
      if (ignore_pdb_hexdigest):
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
      for i in xrange(2): print("*"*79, file=log)
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

map_coefficents_params_str = """\
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

def determine_experimental_phases(reflection_file_server,
                                  parameters,
                                  log,
                                  parameter_scope,
                                  working_point_group,
                                  symmetry_safety_check,
                                  ignore_all_zeros = True):
  """
  Extract experimental phases from the given inputs if possible.  Returns None
  if not found.
  """
  try:
    experimental_phases = \
      reflection_file_server.get_experimental_phases(
        file_name        = parameters.file_name,
        labels           = parameters.labels,
        ignore_all_zeros = ignore_all_zeros,
        parameter_scope  = parameter_scope)
  except reflection_file_utils.Sorry_No_array_of_the_required_type:
    experimental_phases = None
  else:
    parameters.file_name = experimental_phases.info().source
    parameters.labels = [experimental_phases.info().label_string()]
    print("Experimental phases:", file=log)
    print(" ", experimental_phases.info(), file=log)
    miller_array_symmetry_safety_check(
      miller_array          = experimental_phases,
      data_description      = "Experimental phases",
      working_point_group   = working_point_group,
      symmetry_safety_check = symmetry_safety_check,
      log                   = log)
    print(file=log)
    info = experimental_phases.info()
    processed = experimental_phases.eliminate_sys_absent(log = log)
    if(processed is not experimental_phases):
       info = info.customized_copy(systematic_absences_eliminated = True)
    if(not processed.is_unique_set_under_symmetry()):
       print("Merging symmetry-equivalent Hendrickson-Lattman coefficients:", file=log)
       merged = processed.merge_equivalents()
       merged.show_summary(out = log, prefix="  ")
       print(file=log)
       processed = merged.array()
       info = info.customized_copy(merged = True)
    return processed.set_info(info)
