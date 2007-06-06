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
from iotbx import mtz
from scitbx.python_utils.misc import user_plus_sys_time, show_total_time
from libtbx.str_utils import show_string
from libtbx import adopt_init_args
import random, sys, os
from libtbx.test_utils import approx_equal
from mmtbx.refinement import print_statistics
import libtbx.load_env


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

def explain_how_to_generate_array_of_r_free_flags(neutron_flag):
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
  if(not neutron_flag):
     print >> self.log, part1 + \
                     """refinement.main.generate_r_free_flags=True""" + part3
  else:
      print >> self.log, part1 + \
             """refinement.main.generate_neutron_r_free_flags=True""" + part3

data_params = iotbx.phil.parse("""\
  file_name=None
    .type=path
  labels=None
    .type=strings
""")

def determine_data(reflection_file_server,
                   parameters,
                   parameter_scope,
                   log,
                   data_description,
                   working_point_group,
                   symmetry_safety_check,
                   ignore_all_zeros = True):
  data = reflection_file_server.get_xray_data(
    file_name        = parameters.file_name,
    labels           = parameters.labels,
    ignore_all_zeros = ignore_all_zeros,
    parameter_scope  = parameter_scope)
  parameters.file_name = data.info().source
  parameters.labels = [data.info().label_string()]
  if(data.is_xray_intensity_array()):
     print >> log, "I-obs:"
  else:
     print >> log, "F-obs:"
  print >> log, " ", data.info()
  miller_array_symmetry_safety_check(
    miller_array          = data,
    data_description      = data_description,
    working_point_group   = working_point_group,
    symmetry_safety_check = symmetry_safety_check,
    log                   = log)
  print >> log
  info = data.info()
  processed = data.eliminate_sys_absent(log = log)
  if(processed is not data):
     info = info.customized_copy(systematic_absences_eliminated = True)
  if(not processed.is_unique_set_under_symmetry()):
     if(data.is_xray_intensity_array()):
        print >> log, "Merging symmetry-equivalent intensities:"
     else:
        print >> log, "Merging symmetry-equivalent amplitudes:"
     merged = processed.merge_equivalents()
     merged.show_summary(out = log, prefix="  ")
     print >> log
     processed = merged.array()
     info = info.customized_copy(merged=True)
  return processed.set_info(info)

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

r_free_flags_params = iotbx.phil.parse("""\
  file_name=None
    .type=path
    expert_level=0
  label=None
    .type=str
    expert_level=0
  test_flag_value=None
    .type=int
    expert_level=0
  disable_suitability_test=False
    .type=bool
    .expert_level=2
""")

def determine_r_free_flags(reflection_file_server,
                           data,
                           generate_r_free_flags,
                           parameters,
                           parameter_scope,
                           working_point_group,
                           symmetry_safety_check,
                           log,
                           neutron_flag,
                           r_free_flags_fraction=None,
                           r_free_flags_max_free=None,
                           r_free_flags_lattice_symmetry_max_delta=None,
                           r_free_flags_use_lattice_symmetry=None):
  r_free_flags = None
  if(generate_r_free_flags is None or not generate_r_free_flags):
    try:
      r_free_flags, test_flag_value = \
        reflection_file_server.get_r_free_flags(
          file_name                = parameters.file_name,
          label                    = parameters.label,
          test_flag_value          = parameters.test_flag_value,
          disable_suitability_test = parameters.disable_suitability_test,
          parameter_scope          = parameter_scope)
    except reflection_file_utils.Sorry_No_array_of_the_required_type, e:
      e.reset_tracebacklimit()
      if(generate_r_free_flags is not None):
         explain_how_to_generate_array_of_r_free_flags(
                                                   neutron_flag = neutron_flag)
         raise Sorry("Please try again.")
      r_free_flags, test_flag_value = None, None
    else:
      parameters.file_name = r_free_flags.info().source
      parameters.label = r_free_flags.info().label_string()
      parameters.test_flag_value = test_flag_value
      print >> log, "R-free flags:"
      print >> log, " ", r_free_flags.info()
      miller_array_symmetry_safety_check(
        miller_array          = r_free_flags,
        data_description      = "R-free flags",
        working_point_group   = working_point_group,
        symmetry_safety_check = symmetry_safety_check,
        log                   = log)
      print >> log
      info = r_free_flags.info()
      processed = r_free_flags.eliminate_sys_absent(log = log)
      if(processed is not r_free_flags):
         info = info.customized_copy(systematic_absences_eliminated = True)
      if(not processed.is_unique_set_under_symmetry()):
         print >> log, \
           "Checking symmetry-equivalent R-free flags for consistency:",
         try:
           merged = processed.merge_equivalents()
         except RuntimeError, e:
           print >> log
           error_message = str(e)
           expected_error_message = "cctbx Error: merge_equivalents_exact: "
           assert error_message.startswith(expected_error_message)
           raise Sorry("Incompatible symmetry-equivalent R-free flags: %s" %
             error_message[len(expected_error_message):])
         else:
           print >> log, "OK"
           print >> log
         processed = merged.array()
         info = info.customized_copy(merged=True)
         del merged
      r_free_flags = processed.set_info(info)
  if(r_free_flags is None):
     assert [r_free_flags_fraction,
             r_free_flags_max_free,
             r_free_flags_lattice_symmetry_max_delta,
             r_free_flags_use_lattice_symmetry].count(None) == 0
     print >> log, "*"*79
     print >> log, "Generating a new array of R-free flags."
     print >> log, "*"*79
     print >> log
     r_free_flags = data.generate_r_free_flags(
          fraction= r_free_flags_fraction,
          max_free= r_free_flags_max_free,
          lattice_symmetry_max_delta = r_free_flags_lattice_symmetry_max_delta,
          use_lattice_symmetry = r_free_flags_use_lattice_symmetry
        ).set_info(miller.array_info(labels=["R-free-flags"]))
     parameters.label = r_free_flags.info().label_string()
     parameters.test_flag_value = True
  return r_free_flags

pdb_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .multiple=True
""")

def process_pdb_file(pdb_file_names,
                     parameters,
                     pdb_interpretation_params,
                     mon_lib_srv,
                     ener_lib,
                     crystal_symmetry,
                     log):
  if(len(pdb_file_names) == 0):
     raise Sorry("No coordinate file given.")
  parameters.file_name = [
    os.path.abspath(file_name) for file_name in pdb_file_names]
  if(len(pdb_file_names) == 1):
     pdb_file_name = pdb_file_names[0]
     raw_records = None
     pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  else:
    pdb_file_name = None
    raw_records = []
    raw_records_flex = flex.std_string()
    for file_name in pdb_file_names:
      raw_records.extend(open(file_name).readlines())
      raw_records_flex.extend(flex.split_lines(open(file_name).read()))
    pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records_flex)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv = mon_lib_srv,
    ener_lib                 = ener_lib,
    params                   = pdb_interpretation_params,
    file_name                = pdb_file_name,
    raw_records              = raw_records,
    strict_conflict_handling = False,
    crystal_symmetry         = crystal_symmetry,
    force_symmetry           = True,
    log                      = log)
  print >> log
  msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
  if(msg is not None):
     msg = "\n  ".join([msg,
       "Please edit the PDB file to resolve the problems and/or supply a",
       "CIF file with matching restraint definitions, along with",
       "apply_cif_modification and apply_cif_link parameter definitions",
       "if necessary (see phenix.refine documentation).",
       "Also note that elbow.builder is available to create restraint",
       "definitions for unknown ligands."])
     raise Sorry(msg)
  return processed_pdb_file, pdb_inp

cif_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .multiple=True
""")

def process_monomer_cif_files(cif_objects, parameters, mon_lib_srv, ener_lib):
  all = []
  index_dict = {}
  for file_name in parameters.file_name:
    file_name = libtbx.path.canonical_path(file_name=file_name)
    index_dict[file_name] = len(all)
    all.append((file_name,None))
  for file_name,cif_object in cif_objects:
    file_name = libtbx.path.canonical_path(file_name=file_name)
    index_dict[file_name] = len(all)
    all.append((file_name,cif_object))
  unique_indices = index_dict.values()
  unique_indices.sort()
  unique = flex.select(sequence=all, permutation=unique_indices)
  del parameters.file_name[:]
  for file_name,cif_object in unique:
    if (cif_object is None):
      mon_lib_srv.process_cif(file_name=file_name)
      ener_lib.process_cif(file_name=file_name)
    else:
      mon_lib_srv.process_cif_object(
        cif_object=cif_object, file_name=file_name)
      ener_lib.process_cif_object(cif_object=cif_object, file_name=file_name)
    parameters.file_name.append(file_name)

class extract_f_obs(object):
  def __init__(self, f_obs, params, neutron_data, log):
    if(params.main.scattering_table == "neutron" or neutron_data):
       header = "Neutron data"
    else:
       header = "X-ray data"
    mmtbx.refinement.print_statistics.make_header(header, out=log)
    f_obs.show_comprehensive_summary(f = log)
    print >> log
    if(f_obs.is_complex_array()): f_obs = abs(f_obs)
    if(f_obs.is_xray_intensity_array()):
       selection_by_isigma = self._apply_sigma_cutoff(
         f_obs   = f_obs,
         n       = params.main.sigma_iobs_rejection_criterion,
         message = "Number of reflections with |Iobs|/sigma(Iobs) < %5.2f: %d",
         log     = log)
       if(selection_by_isigma is not None):
          f_obs = f_obs.select(selection_by_isigma)
       print >> log, \
                   "Intensities converted to amplitudes for use in refinement."
       f_obs = f_obs.f_sq_as_f()
       print >> log
    f_obs.set_observation_type_xray_amplitude()
    f_obs = f_obs.map_to_asu()
    selection = f_obs.all_selection()
    if(neutron_data):
       d_max = params.neutron.low_resolution
       d_min = params.neutron.high_resolution
    else:
       d_max = params.main.low_resolution
       d_min = params.main.high_resolution
    if(d_max is not None):
       selection &= f_obs.d_spacings().data() <= d_max
    if(d_min is not None):
       selection &= f_obs.d_spacings().data() >= d_min
    selection_strictly_positive = f_obs.data() > 0
    print >> log, \
      "Number of F-obs in resolution range:                  ", \
      selection.count(True)
    print >> log, \
      "Number of F-obs <= 0:                                 ", \
      selection_strictly_positive.count(False)
    selection &= selection_strictly_positive
    selection_by_fsigma = self._apply_sigma_cutoff(
         f_obs   = f_obs,
         n       = params.main.sigma_fobs_rejection_criterion,
         message = "Number of reflections with |Fobs|/sigma(Fobs) < %5.2f: %d",
         log     = log)
    if(selection_by_fsigma is not None): selection &= selection_by_fsigma
    if(f_obs.sigmas() is not None):
       sigma_cutoff = params.main.sigma_fobs_rejection_criterion
       if(sigma_cutoff is not None and sigma_cutoff > 0):
          selection_by_sigma = f_obs.data() > f_obs.sigmas()*sigma_cutoff
          print >> log, \
            "Number of reflections with |Fobs|/sigma(Fobs) < %5.2f: %d" % (
              sigma_cutoff, selection_by_sigma.count(False))
          selection &= selection_by_sigma
    selection &= f_obs.d_star_sq().data() > 0
    f_obs = f_obs.select(selection)
    rr = f_obs.resolution_range()
    print >> log, "Refinement resolution range: d_max = %8.4f" % rr[0]
    print >> log, "                             d_min = %8.4f" % rr[1]
    print >> log
    if(f_obs.indices().size() == 0):
       raise Sorry(
             "No data left after applying resolution limits and sigma cutoff.")
    if(not neutron_data):
       if(params.main.force_anomalous_flag_to_be_equal_to is not None):
          if(not params.main.force_anomalous_flag_to_be_equal_to):
             print >> log, \
               "refinement.main.force_anomalous_flag_to_be_equal_to = False"
             if(f_obs.anomalous_flag()):
                print >> log, "Reducing X-ray data to non-anomalous array."
                merged = f_obs.as_non_anomalous_array().merge_equivalents()
                merged.show_summary(out=log, prefix="  ")
                f_obs = merged.array().set_observation_type( f_obs )
                del merged
                print >> log
          elif(not f_obs.anomalous_flag()):
             print >> log, \
               "refinement.main.force_anomalous_flag_to_be_equal_to = True"
             print >> log, "Generating Bijvoet mates of X-ray data."
             f_obs = f_obs.generate_bijvoet_mates()
             print >> log
       twin_analyses = mmtbx.scaling.twin_analyses.twin_analyses_brief(
                               miller_array = f_obs, cut_off = 4.0, out = None)
       if(twin_analyses   == True):  result = "twinned"
       elif(twin_analyses == False): result = "not twinned"
       elif(twin_analyses is None):  result = "unknown"
       else: raise Sorry("Twin analysis failed.")
       print >> log, "Twin analysis result: ", result
       self.twin_analyses_ = twin_analyses
    print >> log
    self.f_obs_ = f_obs
