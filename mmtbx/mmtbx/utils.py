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

def explain_how_to_generate_array_of_r_free_flags(neutron_flag, log):
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
     print >> log, part1 + \
                     """refinement.main.generate_r_free_flags=True""" + part3
  else:
      print >> log, part1 + \
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
    .help=Name of file containing free-R flags
    .expert_level=0
  label=None
    .type=str
    .expert_level=0
  test_flag_value=None
    .type=int
    .expert_level=0
  disable_suitability_test=False
    .type=bool
    .help=XXX
    .expert_level=2
  ignore_pdb_hexdigest=False
    .type=bool
    .help="If True, disables safety check based on MD5 hexdigests stored in"
          "PDB files produced by previous runs."
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
                                        neutron_flag = neutron_flag, log = log)
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
     parameters.test_flag_value = 1
  return r_free_flags

pdb_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .help=Model file(s) name (PDB)
    .multiple=True
""")

def process_pdb_file(pdb_file_names,
                     mon_lib_srv,
                     ener_lib,
                     crystal_symmetry,
                     parameters = None,
                     pdb_interpretation_params = None,
                     stop_for_unknowns = True,
                     log = None):
  if(log is None): log = sys.stdout
  if(len(pdb_file_names) == 0):
     raise Sorry("No coordinate file given.")
  if(parameters is not None):
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
    mon_lib_srv              = mon_lib_srv,
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
  if(msg is not None and stop_for_unknowns):
     msg = "\n  ".join([msg,
       "Please edit the PDB file to resolve the problems and/or supply a",
       "CIF file with matching restraint definitions, along with",
       "apply_cif_modification and apply_cif_link parameter definitions",
       "if necessary (see phenix.refine documentation).",
       "Also note that phenix.elbow is available to create restraint",
       "definitions for unknown ligands."])
     raise Sorry(msg)
  return processed_pdb_file, pdb_inp

cif_params = iotbx.phil.parse("""\
  file_name=None
    .optional=True
    .type=path
    .help=Monomer file(s) name (CIF)
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

def get_atom_selections(all_chain_proxies,
                        xray_structure,
                        selection_strings     = None,
                        iselection            = True,
                        one_group_per_residue = False,
                        hydrogens_only        = False):
  ###> Fix aal: need for proper selection
  need_to_fix_aal = False
  aal = all_chain_proxies.stage_1.atom_attributes_list
  for aal_i in aal:
    if(aal_i.element.strip() == ''):
      need_to_fix_aal = True
      break
  if(need_to_fix_aal):
    scatterers = xray_structure.scatterers()
    for aal_i, sc in zip(aal,scatterers):
      aal_i.element = sc.element_symbol()
    all_chain_proxies.stage_1.selection_cache(
                                           force_selection_cache_update = True)
  #
  if(hydrogens_only):
     assert xray_structure is not None
  if(selection_strings is None or isinstance(selection_strings, str)):
     selection_strings = [selection_strings]
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
  return selections

def atom_selection(all_chain_proxies, string):
  try: return all_chain_proxies.selection(string = string)
  except KeyboardInterrupt: raise
  except Exception: raise Sorry("Invalid atom selection: %s" % string)

def write_pdb_file(xray_structure,
                   atom_attributes_list,
                   selection = None,
                   out = None):
  if(out is None): out = sys.stdout
  crystal_symmetry = xray_structure.crystal_symmetry()
  print >> out, pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
  print >> out, pdb.format_scale_records(
                                      unit_cell = crystal_symmetry.unit_cell())
  xrs = xray_structure
  sites_cart  = xrs.sites_cart()
  scatterers  = xrs.scatterers()
  occupancies = scatterers.extract_occupancies()
  u_carts = scatterers.extract_u_cart_or_u_cart_plus_u_iso(xrs.unit_cell())
  u_isos      = xrs.extract_u_iso_or_u_equiv()
  scat_types  = scatterers.extract_scattering_types()
  #XXX high duplication
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
         print >> out, pdb.format_atom_record(
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
         if(scatterers[i_seq].flags.use_u_aniso()):
            print >> out, pdb.format_anisou_record(
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
            print >> out, pdb.format_atom_record(
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
            if(scatterers[i_seq].flags.use_u_aniso()):
               print >> out, pdb.format_anisou_record(
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
     print >> out, "END"

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
