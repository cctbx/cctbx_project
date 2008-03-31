from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
import cctbx.xray.structure_factors
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, host_and_user, multi_out
import iotbx.phil
import libtbx.phil.command_line
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
from iotbx import pdb
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from scitbx.math import matrix
from cctbx import adptbx
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from iotbx.option_parser import iotbx_option_parser
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx.pdb import combine_unique_pdb_files
from iotbx import mtz
from scitbx.python_utils.misc import user_plus_sys_time, show_total_time
from libtbx.str_utils import show_string
from libtbx import adopt_init_args
import random, sys, os
from libtbx.test_utils import approx_equal
from mmtbx.refinement import print_statistics
import libtbx.load_env
from mmtbx.solvent import ordered_solvent

import boost.python
utils_ext = boost.python.import_ext("mmtbx_utils_ext")
from mmtbx_utils_ext import *

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
        print >> log, "*" * 79
        print >> log, "WARNING:", msg
        print >> log, "*" * 79
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
  print >> log, part1 + flags_parameter_scope+""".generate=True""" + part3

data_and_flags = iotbx.phil.parse("""\
  file_name = None
    .type=path
  labels = None
    .type=strings
  high_resolution = None
    .type=float
  low_resolution = None
    .type=float
  outliers_rejection = True
    .type=bool
  sigma_fobs_rejection_criterion = 0.0
    .type=float
  sigma_iobs_rejection_criterion = 0.0
    .type=float
  ignore_all_zeros = True
    .type=bool
  force_anomalous_flag_to_be_equal_to = None
    .type=bool
  r_free_flags {
    file_name = None
      .type=path
    label = None
      .type=str
    test_flag_value = None
      .type=int
    disable_suitability_test = False
      .type=bool
    ignore_pdb_hexdigest = False
      .type=bool
      .help=If True, disables safety check based on MD5 hexdigests stored in \
            PDB files produced by previous runs.
      .expert_level=2
    ignore_r_free_flags = False
      .type=bool
      .help = Use all reflections in refinement (work and test)
    generate = False
      .type=bool
      .help = Generate R-free flags (if not available in input files)
    fraction = 0.1
      .type=float
    max_free = 2000
      .type=int
    lattice_symmetry_max_delta = 5
      .type=float
    use_lattice_symmetry = True
      .type=bool
  }
""")

class determine_data_and_flags(object):
  def __init__(self, reflection_file_server,
                     parameters,
                     data_parameter_scope = "",
                     flags_parameter_scope = "",
                     data_description = None,
                     working_point_group = None,
                     symmetry_safety_check = None,
                     remark_r_free_flags_md5_hexdigest = None,
                     extract_r_free_flags = True,
                     log = None):
    adopt_init_args(self, locals())
    self.intensity_flag = False
    self.f_obs = None
    self.r_free_flags = None
    self.test_flag_value = None
    self.r_free_flags_md5_hexdigest = None
    if(data_description is not None):
      print_statistics.make_header(data_description, out = log)
    self.raw_data = self.extract_data()
    data_info = self.raw_data.info()
    if(extract_r_free_flags):
      self.raw_flags = self.extract_flags(data = self.raw_data)
      flags_info = self.raw_flags.info()
    self.f_obs = self.data_as_f_obs(f_obs = self.raw_data)
    if(extract_r_free_flags):
      self.r_free_flags,self.test_flag_value,self.r_free_flags_md5_hexdigest =\
        self.flags_as_r_free_flags(f_obs = self.f_obs, r_free_flags = self.raw_flags)
      self.r_free_flags.set_info(flags_info)
    self.f_obs.set_info(data_info)

  def extract_data(self):
    data = self.reflection_file_server.get_xray_data(
      file_name        = self.parameters.file_name,
      labels           = self.parameters.labels,
      ignore_all_zeros = self.parameters.ignore_all_zeros,
      parameter_scope  = self.data_parameter_scope)
    self.parameters.file_name = data.info().source
    self.parameters.labels = [data.info().label_string()]
    if(data.is_xray_intensity_array()):
      print >> self.log, "I-obs:"
      self.intensity_flag = True
    else:
      print >> self.log, "F-obs:"
    print >> self.log, " ", data.info()
    if([self.data_description, self.working_point_group,
       self.symmetry_safety_check].count(None) == 0):
      miller_array_symmetry_safety_check(
        miller_array          = data,
        data_description      = self.data_description,
        working_point_group   = self.working_point_group,
        symmetry_safety_check = self.symmetry_safety_check,
        log                   = self.log)
      print >> self.log
    info = data.info()
    processed = data.eliminate_sys_absent(log = self.log)
    if(processed is not data):
      info = info.customized_copy(systematic_absences_eliminated = True)
    if(not processed.is_unique_set_under_symmetry()):
      if(data.is_xray_intensity_array()):
        print >> self.log, "Merging symmetry-equivalent intensities:"
      else:
        print >> self.log, "Merging symmetry-equivalent amplitudes:"
      merged = processed.merge_equivalents()
      merged.show_summary(out = self.log, prefix="  ")
      print >> self.log
      processed = merged.array()
      info = info.customized_copy(merged=True)
    return processed.set_info(info)

  def extract_flags(self, data, data_description = "R-free flags"):
    r_free_flags, test_flag_value = None, None
    params = self.parameters.r_free_flags
    if(not self.parameters.r_free_flags.generate):
      try:
        r_free_flags, test_flag_value = \
          self.reflection_file_server.get_r_free_flags(
            file_name                = params.file_name,
            label                    = params.label,
            test_flag_value          = params.test_flag_value,
            disable_suitability_test = params.disable_suitability_test,
            parameter_scope          = self.flags_parameter_scope)
      except reflection_file_utils.Sorry_No_array_of_the_required_type, e:
        e.reset_tracebacklimit()
        if(self.parameters.r_free_flags.generate is not None):
          explain_how_to_generate_array_of_r_free_flags(log = self.log,
            flags_parameter_scope = self.flags_parameter_scope)
          raise Sorry("Please try again.")
        r_free_flags, test_flag_value = None, None
      else:
        params.file_name = r_free_flags.info().source
        params.label = r_free_flags.info().label_string()
        params.test_flag_value = test_flag_value
        print >> self.log, data_description+":"
        print >> self.log, " ", r_free_flags.info()
        if([self.working_point_group,
           self.symmetry_safety_check].count(None) == 0):
          miller_array_symmetry_safety_check(
            miller_array          = r_free_flags,
            data_description      = data_description,
            working_point_group   = self.working_point_group,
            symmetry_safety_check = self.symmetry_safety_check,
            log                   = self.log)
          print >> self.log
        info = r_free_flags.info()
        processed = r_free_flags.eliminate_sys_absent(log = self.log)
        if(processed is not r_free_flags):
          info = info.customized_copy(systematic_absences_eliminated = True)
        if(not processed.is_unique_set_under_symmetry()):
           print >> self.log, \
             "Checking symmetry-equivalent R-free flags for consistency:",
           try:
             merged = processed.merge_equivalents()
           except RuntimeError, e:
             print >> self.log
             error_message = str(e)
             expected_error_message = "cctbx Error: merge_equivalents_exact: "
             assert error_message.startswith(expected_error_message)
             raise Sorry("Incompatible symmetry-equivalent R-free flags: %s" %
               error_message[len(expected_error_message):])
           else:
             print >> self.log, "OK"
             print >> self.log
           processed = merged.array()
           info = info.customized_copy(merged=True)
           del merged
        r_free_flags = processed.set_info(info)
    if(r_free_flags is None):
      assert [params.fraction,
              params.max_free,
              params.lattice_symmetry_max_delta,
              params.use_lattice_symmetry].count(None) == 0
      print >> self.log, "Generating a new array of R-free flags."
      print >> self.log
      r_free_flags = data.generate_r_free_flags(
        fraction                   = params.fraction,
        max_free                   = params.max_free,
        lattice_symmetry_max_delta = params.lattice_symmetry_max_delta,
        use_lattice_symmetry       = params.use_lattice_symmetry
        ).set_info(miller.array_info(labels = ["R-free-flags"]))
      params.label = r_free_flags.info().label_string()
      params.test_flag_value = 1
    return r_free_flags

  def data_as_f_obs(self, f_obs):
    f_obs.show_comprehensive_summary(f = self.log)
    f_obs_data_size = f_obs.data().size()
    print >> self.log
    d_min = f_obs.d_min()
    if(d_min < 0.25): # XXX what is the equivalent for neutrons ???
      raise Sorry("Resolution of data is too high: %-6.4f A"%d_min)
    if(f_obs.is_complex_array()): f_obs = abs(f_obs)
    if(f_obs.is_xray_intensity_array()):
      selection_by_isigma = self._apply_sigma_cutoff(
        f_obs   = f_obs,
        n       = self.parameters.sigma_iobs_rejection_criterion,
        message = "Number of reflections with |Iobs|/sigma(Iobs) < %5.2f: %d")
      if(selection_by_isigma is not None):
        f_obs = f_obs.select(selection_by_isigma)
      print >> self.log, \
        "Intensities converted to amplitudes for use in refinement."
      f_obs = f_obs.f_sq_as_f()
      print >> self.log
    f_obs.set_observation_type_xray_amplitude()
    f_obs = f_obs.map_to_asu()
    selection = f_obs.all_selection()
    if(self.parameters.low_resolution is not None):
      selection &= f_obs.d_spacings().data() <= self.parameters.low_resolution
    if(self.parameters.high_resolution is not None):
      selection &= f_obs.d_spacings().data() >= self.parameters.high_resolution
    selection_strictly_positive = f_obs.data() > 0
    print >> self.log, \
      "Number of F-obs in resolution range:                  ", \
      selection.count(True)
    print >> self.log, \
      "Number of F-obs <= 0:                                 ", \
      selection_strictly_positive.count(False)
    selection &= selection_strictly_positive
    selection_by_fsigma = self._apply_sigma_cutoff(
      f_obs   = f_obs,
      n       = self.parameters.sigma_fobs_rejection_criterion,
      message = "Number of reflections with |Fobs|/sigma(Fobs) < %5.2f: %d")
    if(selection_by_fsigma is not None): selection &= selection_by_fsigma
    selection &= f_obs.d_star_sq().data() > 0
    f_obs = f_obs.select(selection)
    rr = f_obs.resolution_range()
    print >> self.log, "Refinement resolution range: d_max = %8.4f" % rr[0]
    print >> self.log, "                             d_min = %8.4f" % rr[1]
    print >> self.log
    if(f_obs.indices().size() == 0):
      raise Sorry(
        "No data left after applying resolution limits and sigma cutoff.")
    if(self.parameters.force_anomalous_flag_to_be_equal_to is not None):
      if(not self.parameters.force_anomalous_flag_to_be_equal_to):
        print >> self.log, "force_anomalous_flag_to_be_equal_to=False"
        if(f_obs.anomalous_flag()):
          print >> self.log, "Reducing data to non-anomalous array."
          merged = f_obs.as_non_anomalous_array().merge_equivalents()
          merged.show_summary(out = self.log, prefix="  ")
          f_obs = merged.array().set_observation_type( f_obs )
          del merged
          print >> self.log
      elif(not f_obs.anomalous_flag()):
        print >> self.log, "force_anomalous_flag_to_be_equal_to=True"
        print >> self.log, "Generating Bijvoet mates of X-ray data."
        observation_type = f_obs.observation_type()
        f_obs = f_obs.generate_bijvoet_mates()
        f_obs.set_observation_type(observation_type)
        print >> self.log
    if(f_obs_data_size != f_obs.data().size()):
      print >> self.log, "\nFobs statistics after all cutoffs applied:\n"
      f_obs.show_comprehensive_summary(f = self.log)
    return f_obs

  def _apply_sigma_cutoff(self, f_obs, n, message):
    selection = None
    if(f_obs.sigmas() is not None):
      sigma_cutoff = n
      if(sigma_cutoff is not None and sigma_cutoff > 0):
        selection_by_sigma = f_obs.data() > f_obs.sigmas()*sigma_cutoff
        print >> self.log, message % (sigma_cutoff,
          selection_by_sigma.count(False))
        selection = selection_by_sigma
    return selection

  def flags_as_r_free_flags(self, f_obs, r_free_flags):
    test_flag_value = self.parameters.r_free_flags.test_flag_value
    r_free_flags.show_comprehensive_summary(f = self.log)
    print >> self.log
    print >> self.log, "Test (R-free flags) flag value:", test_flag_value
    print >> self.log
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
        print >> self.log, "Reducing R-free flags to non-anomalous array."
        r_free_flags = r_free_flags.average_bijvoet_mates()
        print >> self.log
    elif(not r_free_flags.anomalous_flag()):
       print >> self.log, "Generating Bijvoet mates of R-free flags."
       r_free_flags = r_free_flags.generate_bijvoet_mates()
       print >> self.log
    r_free_flags = r_free_flags.map_to_asu().common_set(f_obs)
    n_missing_r_free_flags = f_obs.indices().size() \
      - r_free_flags.indices().size()
    if(n_missing_r_free_flags != 0):
      raise Sorry("R-free flags not compatible with F-obs array:"
        " missing flag for %d F-obs selected for refinement." %
        n_missing_r_free_flags)
    r_free_flags.show_r_free_flags_info(out = self.log, prefix="")
    return r_free_flags, test_flag_value, r_free_flags_md5_hexdigest

  def verify_r_free_flags_md5_hexdigest(self,
        ignore_pdb_hexdigest,
        current,
        records):
    from_file = {}
    for record in records:
      flds = record.split()
      if (len(flds) == 3):
        from_file[flds[2]] = None
    if (len(from_file) > 1):
      raise Sorry(
        "Multiple conflicting REMARK r_free_flags.md5.hexdigest records"
        " found in the input PDB file.")
    if (len(from_file) == 1 and current not in from_file):
      log = self.log
      for i in xrange(2): print >> log, "*"*79
      if (ignore_pdb_hexdigest):
        print >> log
        print >> log, " ".join(["WARNING"]*9)
      print >> log, """
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
  current, from_file.keys()[0]),
      if (not ignore_pdb_hexdigest):
        print >> log, """\
In this case,
simply remove the

  REMARK r_free_flags.md5.hexdigest %s

record from the input PDB file to proceed with the refinement.""" % (
  from_file.keys()[0]),
      print >> log, """

Otherwise it is best to recover the previously used R-free flags
and use them consistently throughout the refinement of the model.
Run this command again with the name of the file containing the
original flags as an additional input.
"""
      if (not ignore_pdb_hexdigest):
        print >> log, """\
If the original R-free flags are unrecoverable, remove the REMARK
record as indicated above. In this case the values for R-free will
become meaningful only after many cycles of refinement.
"""
      else:
        print >> log, """\
If the original R-free flags are unrecoverable, the values for R-free
will become meaningful only after many cycles of refinement.
"""
      for i in xrange(2): print >> log, "*"*79
      print >> log
      if (not ignore_pdb_hexdigest):
        raise Sorry("Please resolve the R-free flags mismatch.")


experimental_phases_params = iotbx.phil.parse("""\
  file_name=None
    .type=path
  labels=None
    .type=strings
""")

def determine_experimental_phases(reflection_file_server,
                                  parameters,
                                  log,
                                  parameter_scope,
                                  working_point_group,
                                  symmetry_safety_check,
                                  ignore_all_zeros = True):
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
    print >> log, "Experimental phases:"
    print >> log, " ", experimental_phases.info()
    miller_array_symmetry_safety_check(
      miller_array          = experimental_phases,
      data_description      = "Experimental phases",
      working_point_group   = working_point_group,
      symmetry_safety_check = symmetry_safety_check,
      log                   = log)
    print >> log
    info = experimental_phases.info()
    processed = experimental_phases.eliminate_sys_absent(log = log)
    if(processed is not experimental_phases):
       info = info.customized_copy(systematic_absences_eliminated = True)
    if(not processed.is_unique_set_under_symmetry()):
       print >> log, \
         "Merging symmetry-equivalent Hendrickson-Lattman coefficients:"
       merged = processed.merge_equivalents()
       merged.show_summary(out = log, prefix="  ")
       print >> log
       processed = merged.array()
       info = info.customized_copy(merged = True)
    return processed.set_info(info)

pdb_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .help=Model file(s) name (PDB)
    .multiple=True
""")

def get_atom_selections(all_chain_proxies,
                        xray_structure,
                        selection_strings     = None,
                        iselection            = True,
                        one_group_per_residue = False,
                        hydrogens_only        = False,
                        one_selection_array   = False,
                        log                   = None):
  if(log is None): log = sys.stdout
  # XXX Fix atom_attributes_list in-place to be able to select by element type
  # XXX in cases when it's not available.
  aal_modified_counter = 0
  aal = all_chain_proxies.stage_1.atom_attributes_list
  scatterers = xray_structure.scatterers()
  for aal_i, sc in zip(aal,scatterers):
    aal_element = aal_i.element
    if(len(aal_element.strip())==0):
      aal_i.element = sc.element_symbol()
      aal_modified_counter += 1
    if(aal_i.element is None):
      aal_i.element = "  "
  all_chain_proxies.stage_1.selection_cache(force_selection_cache_update=True)
  if(aal_modified_counter > 0):
    print >> log
    print >> log, "atom_attributes_list is modified in-place to define",\
       "missing element types when possible."
    print >> log
  #
  if(hydrogens_only):
     assert xray_structure is not None
  if(selection_strings is None or isinstance(selection_strings, str)):
     selection_strings = [selection_strings]
  elif (len(selection_strings) == 0):
     selection_strings = [None]
  n_none = selection_strings.count(None)
  ss_size = len(selection_strings)
  if((one_group_per_residue and n_none==0) or (ss_size > 1 and n_none > 0)):
     raise Sorry('Ambiguous selection.')
  selections = []
  if(ss_size == 1 and n_none == 1 and not one_group_per_residue):
     selections.append(atom_selection(all_chain_proxies = all_chain_proxies,
                                      string            = "all"))
  elif(one_group_per_residue and ss_size == 1 and n_none == 1):
     assert iselection
     assert len(selections) == 0
     residues = []
     hd_selection = None
     if(hydrogens_only):
        scat_types = xray_structure.scatterers().extract_scattering_types()
        d_selection = (scat_types == "D")
        h_selection = (scat_types == "H")
        hd_selection = d_selection | h_selection
     if(hd_selection is not None and hd_selection.count(True) == 0 and hydrogens_only):
        raise Sorry('No hydrogens to select.')
     for model in all_chain_proxies.stage_1.get_models_and_conformers():
         n_conformers = len(model.conformers)
         for conformer in model.conformers:
             residues_in_conformer = []
             for chain in conformer.get_chains():
                 for residue in chain.residues:
                     if(n_conformers == 1):
                        result = flex.size_t()
                        if(hydrogens_only):
                           for i_sel in residue.iselection:
                               if(scat_types[i_sel] in ["H", "D"]):
                                  result.append(i_sel)
                        else:
                           result = residue.iselection
                        selections.append(result)
                     else:
                        residues_in_conformer.append(residue)
             residues.append(residues_in_conformer)
         if(n_conformers > 1):
            assert len(selections) == 0
            res = []
            unique_res_ids = []
            for residuesi in residues:
              for residue in residuesi:
                res.append(residue)
                if(residue.id() not in unique_res_ids):
                  unique_res_ids.append(residue.id())
            for uri in unique_res_ids:
              s_uri = flex.size_t()
              for r in res:
                if(r.id() == uri):
                  for s in r.iselection:
                    if(s not in s_uri):
                      s_uri.append(s)
              selections.append(s_uri)
            # XXX what is this ?
            if(hydrogens_only):
               for i_seq, sel in enumerate(selections):
                 result = flex.size_t()
                 for si in sel:
                   if(scat_types[si] in ["H", "D"]):
                      result.append(si)
                 selections[i_seq] = result
  elif(ss_size != 1 or n_none == 0 and not one_group_per_residue):
     for selection_string in selection_strings:
        selections.append(atom_selection(all_chain_proxies = all_chain_proxies,
                                         string            = selection_string))
  else:
     raise Sorry('Ambiguous selection.')
  if(iselection):
     for i_seq, selection in enumerate(selections):
       if(hasattr(selection, "iselection")):
          selections[i_seq] = selections[i_seq].iselection()
  if(one_selection_array):
    s0 = selections[0]
    for s in selections[1:]:
      if(not iselection):
        s0 = s0 | s
      else:
        s0.extend(s)
    selections = s0
  return selections

def atom_selection(all_chain_proxies, string):
  try: return all_chain_proxies.selection(string = string)
  except KeyboardInterrupt: raise
  except Exception: raise Sorry("Invalid atom selection: %s" % string)

def write_pdb_file(xray_structure,
                   atom_attributes_list,
                   write_cryst1_record = True,
                   selection = None,
                   out = None):
  raw_records = []
  crystal_symmetry = xray_structure.crystal_symmetry()
  if(write_cryst1_record):
    line = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
    raw_records.append(line)
    if(out is not None):
      print >> out, line
    line = pdb.format_scale_records(unit_cell = crystal_symmetry.unit_cell())
    if(out is not None):
      print >> out, line
    for line_ in line.splitlines():
      raw_records.append(line_)
  xrs = xray_structure
  sites_cart  = xrs.sites_cart()
  scatterers  = xrs.scatterers()
  occupancies = scatterers.extract_occupancies()
  u_carts = scatterers.extract_u_cart_or_u_cart_plus_u_iso(xrs.unit_cell())
  u_isos      = xrs.extract_u_iso_or_u_equiv()
  scat_types  = scatterers.extract_scattering_types()
  if(selection is None):
     for i_seq,atom in enumerate(atom_attributes_list):
         if(atom.name is None): name = "    "
         else: name = atom.name
         if(atom.altLoc is None): altLoc = " "
         else: altLoc = atom.altLoc
         if(atom.chainID is None): chainID = " "
         else: chainID = atom.chainID
         if(atom.resSeq is None): resSeq = 1
         else: resSeq = atom.resSeq
         if(atom.iCode is None): iCode = " "
         else: iCode = atom.iCode
         if(atom.segID is None): segID = "    "
         else: segID = atom.segID
         if(atom.element is None): element = "  "
         else: element = atom.element
         if(atom.charge is None): charge = "  "
         else: charge = atom.charge
         line = pdb.format_atom_record(
                                  record_name = atom.record_name(),
                                  serial      = i_seq+1,
                                  name        = name,
                                  altLoc      = altLoc,
                                  resName     = atom.resName,
                                  chainID     = chainID,
                                  resSeq      = resSeq,
                                  iCode       = iCode,
                                  site        = sites_cart[i_seq],
                                  occupancy   = occupancies[i_seq],
                                  tempFactor  = adptbx.u_as_b(u_isos[i_seq]),
                                  segID       = segID,
                                  element     = scat_types[i_seq],#element,
                                  charge      = charge)
         if(out is not None):
           print >> out, line
         raw_records.append(line)
         if(scatterers[i_seq].flags.use_u_aniso()):
            line = pdb.format_anisou_record(
                                  serial      = i_seq+1,
                                  name        = name,
                                  altLoc      = altLoc,
                                  resName     = atom.resName,
                                  chainID     = chainID,
                                  resSeq      = resSeq,
                                  iCode       = iCode,
                                  u_cart      = u_carts[i_seq],
                                  segID       = segID,
                                  element     = scat_types[i_seq],#element,
                                  charge      = charge)
            if(out is not None):
              print >> out, line
            raw_records.append(line)
     if(out is not None):
       print >> out, "END"
  else:
     for i_seq,atom in enumerate(atom_attributes_list):
         if(selection[i_seq]):
            if(atom.name is None): name = "    "
            else: name = atom.name
            if(atom.altLoc is None): altLoc = " "
            else: altLoc = atom.altLoc
            if(atom.chainID is None): chainID = " "
            else: chainID = atom.chainID
            if(atom.resSeq is None): resSeq = 1
            else: resSeq = atom.resSeq
            if(atom.iCode is None): iCode = " "
            else: iCode = atom.iCode
            if(atom.segID is None): segID = "    "
            else: segID = atom.segID
            if(atom.element is None): element = "  "
            else: element = atom.element
            if(atom.charge is None): charge = "  "
            else: charge = atom.charge
            line = pdb.format_atom_record(
                                  record_name = atom.record_name(),
                                  serial      = i_seq+1,
                                  name        = name,
                                  altLoc      = altLoc,
                                  resName     = atom.resName,
                                  chainID     = chainID,
                                  resSeq      = resSeq,
                                  iCode       = iCode,
                                  site        = sites_cart[i_seq],
                                  occupancy   = occupancies[i_seq],
                                  tempFactor  = adptbx.u_as_b(u_isos[i_seq]),
                                  segID       = segID,
                                  element     = scat_types[i_seq],#element,
                                  charge      = charge)
            if(out is not None):
              print >> out, line
            raw_records.append(line)
            if(scatterers[i_seq].flags.use_u_aniso()):
               line = pdb.format_anisou_record(
                                  serial      = i_seq+1,
                                  name        = name,
                                  altLoc      = altLoc,
                                  resName     = atom.resName,
                                  chainID     = chainID,
                                  resSeq      = resSeq,
                                  iCode       = iCode,
                                  u_cart      = u_carts[i_seq],
                                  segID       = segID,
                                  element     = scat_types[i_seq],#element,
                                  charge      = charge)
               if(out is not None):
                 print >> out, line
               raw_records.append(line)
     if(out is not None):
       print >> out, "END"
  return raw_records

def print_programs_start_header(log, text):
  print >> log
  host_and_user().show(out= log)
  print >> log, date_and_time()
  print >> log
  print >> log, "-"*79
  print >> log, text
  print >> log, "-"*79
  print >> log

def set_log(args):
  log = multi_out()
  if(not "--quiet" in args):
     log.register(label="stdout", file_object=sys.stdout)
  string_buffer = StringIO()
  string_buffer_plots = StringIO()
  log.register(label="log_buffer", file_object=string_buffer)
  sys.stderr = log
  return log

def print_header(line, out=None):
  if (out is None): out = sys.stdout
  header_len = 80
  line_len = len(line)
  assert line_len <= header_len
  fill_len = header_len - line_len
  fill_rl = fill_len/2
  fill_r = fill_rl
  fill_l = fill_rl
  if (fill_rl*2 != fill_len):
    fill_r +=1
  out_string = "\n"+"="*(fill_l-1)+" "+line+" "+"="*(fill_r-1)+"\n"
  if(len(out_string) > 80):
    out_string = "\n"+"="*(fill_l-1)+" "+line+" "+"="*(fill_r-2)+"\n"
  print >> out, out_string
  out.flush()

def get_atom_selection(pdb_file_name, selection_string, iselection = False):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv = monomer_library.server.server(),
    ener_lib    = monomer_library.server.ener_lib(),
    file_name   = pdb_file_name,
    log         = None)
  xray_structure = processed_pdb_file.xray_structure(show_summary = False)
  result = get_atom_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    selection_strings = [selection_string],
    iselection        = iselection)
  assert len(result) == 1
  return result[0]

cif_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .help=Monomer file(s) name (CIF)
    .multiple=True
""")

class process_pdb_file_srv(object):
  def __init__(self, crystal_symmetry          = None,
                     pdb_parameters            = None,
                     pdb_interpretation_params = None,
                     stop_for_unknowns         = None,
                     log                       = None,
                     cif_objects               = None,
                     cif_parameters            = None,
                     mon_lib_srv               = None,
                     ener_lib                  = None):
    self.crystal_symmetry          = crystal_symmetry
    self.pdb_parameters            = pdb_parameters
    self.pdb_interpretation_params = pdb_interpretation_params
    self.stop_for_unknowns         = stop_for_unknowns
    self.cif_objects               = cif_objects
    self.cif_parameters            = cif_parameters
    self.log                       = log
    if(mon_lib_srv is None): self.mon_lib_srv = monomer_library.server.server()
    else: self.mon_lib_srv = mon_lib_srv
    if(ener_lib is None): self.ener_lib = monomer_library.server.ener_lib()
    else: self.ener_lib = ener_lib
    if(self.log is None): self.log = sys.stdout

  def process_pdb_files(self, pdb_file_names = None, raw_records = None):
    assert [pdb_file_names, raw_records].count(None) == 1
    if(self.cif_objects is not None):
      self._process_monomer_cif_files()
    return self._process_pdb_file(pdb_file_names = pdb_file_names,
                                  raw_records    = raw_records)

  def _process_pdb_file(self, pdb_file_names, raw_records):
    if(raw_records is None):
      pdb_combined = combine_unique_pdb_files(file_names=pdb_file_names)
      pdb_combined.report_non_unique(out=self.log)
      if (len(pdb_combined.unique_file_names) == 0):
        raise Sorry("No coordinate file given.")
      if(self.pdb_parameters is not None):
        self.pdb_parameters.file_name = [os.path.abspath(file_name)
          for file_name in pdb_combined.unique_file_names]
      raw_records = pdb_combined.raw_records
    pdb_inp = iotbx.pdb.input(source_info = None,
                              lines       = flex.std_string(raw_records))
    if(pdb_inp.atoms().size() == 0):
      msg = ["No atomic coordinates found in PDB files:"]
      if(pdb_file_names is not None):
        for file_name in pdb_file_names:
          msg.append("  %s" % show_string(file_name))
      raise Sorry("\n".join(msg))
    pdb_inp.construct_hierarchy().overall_counts() \
      .raise_duplicate_atom_labels_if_necessary()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv              = self.mon_lib_srv,
      ener_lib                 = self.ener_lib,
      params                   = self.pdb_interpretation_params,
      raw_records              = raw_records,
      strict_conflict_handling = False,
      crystal_symmetry         = self.crystal_symmetry,
      force_symmetry           = True,
      log                      = self.log)
    print >> self.log
    msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
    if(msg is not None and self.stop_for_unknowns):
       msg = "\n  ".join([msg,
         "Please edit the PDB file to resolve the problems and/or supply a",
         "CIF file with matching restraint definitions, along with",
         "apply_cif_modification and apply_cif_link parameter definitions",
         "if necessary (see phenix.refine documentation).",
         "Also note that phenix.elbow is available to create restraint",
         "definitions for unknown ligands."])
       raise Sorry(msg)
    return processed_pdb_file, pdb_inp

  def _process_monomer_cif_files(self):
    all = []
    index_dict = {}
    if(self.cif_parameters is not None):
      for file_name in self.cif_parameters.file_name:
        file_name = libtbx.path.canonical_path(file_name=file_name)
        index_dict[file_name] = len(all)
        all.append((file_name,None))
    for file_name,cif_object in self.cif_objects:
      file_name = libtbx.path.canonical_path(file_name=file_name)
      index_dict[file_name] = len(all)
      all.append((file_name,cif_object))
    unique_indices = index_dict.values()
    unique_indices.sort()
    unique = flex.select(sequence=all, permutation=unique_indices)
    if(self.cif_parameters is not None): del self.cif_parameters.file_name[:]
    for file_name,cif_object in unique:
      if(cif_object is None):
        self.mon_lib_srv.process_cif(file_name=file_name)
        self.ener_lib.process_cif(file_name=file_name)
      else:
        self.mon_lib_srv.process_cif_object(
          cif_object=cif_object, file_name=file_name)
        self.ener_lib.process_cif_object(cif_object=cif_object,
                                         file_name=file_name)
      if(self.cif_parameters is not None):
        self.cif_parameters.file_name.append(file_name)

def list_3d_as_1d(x):
  result = []
  for i in x:
    for j in i:
      for k in j:
        result.append(k)
  return result

def occupancy_selections(
      all_chain_proxies,
      xray_structure,
      add_water = False,
      other_individual_selection_strings = None,
      other_group_selection_strings = None,
      as_flex_arrays = True):
  if(other_individual_selection_strings is not None and
     len(other_individual_selection_strings) == 0):
    other_individual_selection_strings = None
  if(other_group_selection_strings is not None and
     len(other_group_selection_strings) == 0):
    other_group_selection_strings = None
  # Returns something like this:
  #
  #   One constrained group can conatain from 1 to 4 groups of atoms, and one
  #   group of atoms can conatin from 1 to N atoms, where N - arbitrary.
  #
  #   -------selections for occupancy refinement--------
  #     xxxxxxx  xxxxxxxxxxxxxxxxx  xxxxxxxxxx  xxxxxx
  #        g1        g2       g3      g4  g5      g6
  #   [ [[1,2]], [[3,4,5],[6,7,8]], [[9],[10]], [[11]] ]
  #
  #   For example:
  #     - occupancy of atom 11 is refined as individual (unconstrained group
  #       consisting of one atom)
  #     - occupancies of atoms 9 and 10 are refined as individual with the sum
  #       of them constrained to be equal to 1.0
  #     - one occupancy (q1) is refined for atoms 3,4,5 and one occupancy (q2)
  #       is refined for atoms 6,7,8 and q1+q2=1
  #     - one occupancy is refined per atoms 1,2 and no constrains.
  #
  # Behavior:
  #   always selected (constrains are determined automatically):
  #     - atoms in alternative conformations (constrained: individual if two
  #       atoms are involved or grouped otherwise)
  #     - atoms with partial non-zero occupancies
  #   plus (optional, based on parameters specified):
  #     - user defined selection for individual occupancies
  #     - user defined selection for group occupancies
  #     - add H (D) and/or water if requested
  #
  result = []
  residues = []
  selections = []
  occupancies = xray_structure.scatterers().extract_occupancies()
  for model in all_chain_proxies.stage_1.get_models_and_conformers():
    n_conformers = len(model.conformers)
    for conformer in model.conformers:
      residues_in_conformer = []
      for chain in conformer.get_chains():
        for residue in chain.residues:
          if(n_conformers == 1):
            selections.append(residue.iselection)
          else:
            residues_in_conformer.append(residue)
      residues.append(residues_in_conformer)
    if(n_conformers > 1):
      assert len(selections) == 0
      res = []
      unique_res_ids = []
      for residuesi in residues:
        for residue in residuesi:
          res.append(residue)
          if(residue.id() not in unique_res_ids):
            unique_res_ids.append(residue.id())
      for uri in unique_res_ids:
        s_uri = flex.size_t()
        for r in res:
          if(r.id() == uri):
            for s in r.iselection:
              if(s not in s_uri):
                s_uri.append(s)
        selections.append(s_uri)
  hd_selection = xray_structure.hd_selection().iselection()
  for s in selections:
    s = list(s)
    s_set = set(s)
    res = []
    for rr in residues:
      for r in rr:
        r_is = list(r.iselection)
        r_is_set = set(r_is)
        if(r_is_set.issubset(s_set)):
          res.append(r_is_set)
    zz = []
    for i, si in enumerate(res):
      for j, sj in enumerate(res):
        if(i!=j):
          r = list(si - sj)
          if(r not in zz and r != []): zz.append(r)
    if(zz != []):
      result.append(zz)
  # exclude H
  result_ = []
  for r in result:
    result__ = []
    for rr in r:
      result___ = []
      for rrr in rr:
        if(rrr not in hd_selection):
          result___.append(rrr)
      result__.append(result___)
    result_.append(result__)
  result = result_
  #
  if(other_individual_selection_strings is not None):
    sel = get_atom_selections(
      all_chain_proxies   = all_chain_proxies,
      selection_strings   = other_individual_selection_strings,
      iselection          = True,
      xray_structure      = xray_structure,
      one_selection_array = True)
    result_as_1d_array = list_3d_as_1d(x = result)
    sel_checked = []
    for i in list(sel):
      if(i not in result_as_1d_array):
        sel_checked.append([[i]])
    if(len(sel_checked) > 0):
      result.extend(sel_checked)
  if(other_group_selection_strings is not None):
    sel = get_atom_selections(
      all_chain_proxies   = all_chain_proxies,
      selection_strings   = other_group_selection_strings,
      iselection          = True,
      xray_structure      = xray_structure,
      one_selection_array = False)
    result.extend( [[list(i)] for i in sel] )
  if(add_water):
    result_as_1d_array = list_3d_as_1d(x = result)
    water_selection = get_atom_selections(
      all_chain_proxies   = all_chain_proxies,
      selection_strings   = ['water'],
      iselection          = True,
      xray_structure      = xray_structure,
      one_selection_array = True)
    for w_i_seq in water_selection:
      if(w_i_seq not in result_as_1d_array):
        result.append([[w_i_seq]])
  result_as_1d_array = list_3d_as_1d(x = result)
  for i_seq, occ in enumerate(occupancies):
    if(abs(occ-1.) > 1.e-3 and abs(occ) > 1.e-3 and
       not xray_structure.hd_selection()[i_seq]):
      if(i_seq not in result_as_1d_array):
        result.append([[i_seq]])
  result_as_1d_array = list_3d_as_1d(x=result)
  if(len(result_as_1d_array) != len(set(result_as_1d_array))):
    raise Sorry("Duplicate selection for occupancies.")
  if(as_flex_arrays):
    result_ = []
    for gsel in result:
      result__ = []
      for sel in gsel:
        result__.append(flex.size_t(sel))
      result_.append(result__)
    result = result_
  if(result == []): result = None
  return result

def assert_xray_structures_equal(x1, x2, selection = None, sites = True,
                                 adp = True, occupancies = True):
  assert x1.scatterers().size() == x2.scatterers().size()
  if(selection is not None):
    x1 = x1.select(selection)
    x2 = x2.select(selection)
  if(sites):
    assert approx_equal(x1.sites_frac(), x2.sites_frac())
  if(adp):
    assert approx_equal(x1.extract_u_iso_or_u_equiv(),
                        x2.extract_u_iso_or_u_equiv())
  if(occupancies):
    assert approx_equal(x1.scatterers().extract_occupancies(),
                        x2.scatterers().extract_occupancies())
