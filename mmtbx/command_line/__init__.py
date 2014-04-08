
# XXX this module is tested implicitly in many other regression tests, but
# needs more thorough testing on its own

"""
Generic wrapper for bootstrapping high-level applications which rely on some
combination of model and data, with special attention to geometry restraints
interpretation and f_model setup.  All of the bookkeeping required to
disambiguate crystal symmetry and Miller array conventions is performed
automatically.  This is superficially similar to the setup for phenix.refine
(and re-uses many of the methods in mmtbx.utils), but is somewhat simpler
(single-dataset only) and more general-purpose.
"""

from __future__ import division
from libtbx.str_utils import make_header, make_sub_header
from libtbx.utils import Sorry, Usage, multi_out, null_out
from libtbx import Auto
from cStringIO import StringIO
import sys

cmdline_input_phil_base_str = """
input {
  include scope mmtbx.utils.xray_data_str
  %(phases)s
  %(unmerged)s
  pdb {
    include scope mmtbx.utils.pdb_params
  }
  monomers {
    include scope mmtbx.utils.cif_params
  }
  sequence = None
    .type = path
  scattering_table = wk1995  it1992  *n_gaussian  neutron
    .type = choice
  wavelength = None
    .type = float
  %(phases_flag)s
  %(twin_law)s
}
%(pdb_interpretation)s
"""

def generate_master_phil_with_inputs (
    phil_string,
    enable_twin_law=True,
    enable_automatic_twin_detection=False,
    enable_experimental_phases=False,
    enable_pdb_interpretation_params=False,
    enable_stop_for_unknowns=None,
    enable_full_geometry_params=False,
    enable_unmerged_data=False,
    as_phil_string=False) :
  import iotbx.phil
  phil_extra_dict = {
    "phases" : "",
    "unmerged" : "",
    "phases_flag" : "",
    "twin_law" : "",
    "pdb_interpretation" : "",
  }
  if (enable_automatic_twin_detection) :
    phil_extra_dict["twin_law"] = """
      skip_twin_detection = False
        .type = bool"""
  elif (enable_twin_law) :
    phil_extra_dict["twin_law"] = """
      twin_law = None
        .type = str
        .input_size = 100"""
  if (enable_experimental_phases) :
    phil_extra_dict["phases"] = """
      experimental_phases {
        include scope mmtbx.utils.experimental_phases_params_str
      }"""
    phil_extra_dict["phases_flag"] = """
      use_experimental_phases = Auto
        .type = bool
        .short_caption = Use experimental phases (Hendrickson-Lattman coefficients)
      """
  if (enable_unmerged_data) :
    phil_extra_dict["unmerged"] = """
      unmerged_data {
        file_name = None
          .type = path
        labels = None
          .type = str
      }"""
  if (enable_pdb_interpretation_params) or (enable_full_geometry_params) :
    stop_for_unknowns_params = ""
    if (enable_stop_for_unknowns is not None) :
      stop_for_unknowns_params = """
        stop_for_unknowns = %s
          .type = bool""" % enable_stop_for_unknowns
    phil_extra_dict["pdb_interpretation"] = """
      pdb_interpretation {
        include scope mmtbx.monomer_library.pdb_interpretation.master_params
        %s
      }""" % (stop_for_unknowns_params)
    if (enable_full_geometry_params) :
      phil_extra_dict["pdb_interpretation"] += """
        geometry_restraints
          .alias = refinement.geometry_restraints
        {
          edits
            .short_caption = Custom geometry restraints
          {
            include scope mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
          }
          remove {
            include scope mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str
          }
        }"""
  cmdline_input_phil_str = cmdline_input_phil_base_str % phil_extra_dict
  master_phil_str = """
    %s
    %s
  """ % (cmdline_input_phil_str, phil_string)
  if (as_phil_string) :
    return master_phil_str
  return iotbx.phil.parse(master_phil_str, process_includes=True)

def generic_simple_input_phil () :
  return generate_master_phil_with_inputs(
    phil_string="",
    enable_automatic_twin_detection=True)

class load_model_and_data (object) :
  """
  Class for processing command-line input and creating necessary objects.
  The master_phil object should include cmdline_input_phil_str above, plus
  any application-specific parameters.  Programs which use this can be invoked
  using simple file arguments or explicit parameters, e.g.

    mmtbx.some_program model.pdb data.mtz

  This class performs the following functions (mostly using other wrappers
  elsewhere in mmtbx.utils):
    1. Process all arguments and extract as Python parameters
    2. Read in data and R-free flags
    3. Filter data and flags to be consistent if necessary
    4. Read in PDB file, either using the iotbx.pdb API, or if process_pdb_file
       is True, mmtbx.monomer_library.pdb_interpretation (using any CIF files
       included in the inputs)
    5. Extract the pdb_hierarchy and xray_structure objects.
    6. Create an mmtbx.f_model.manager object using the data, flags, and
       xray_structure.
  If at any point the inputs are ambiguous, hopefully the program will stop
  and raise an interpretable error.
  """
  def __init__ (self,
      args,
      master_phil,
      update_f_part1_for="refinement",
      out=sys.stdout,
      process_pdb_file=True,
      use_conformation_dependent_library=False, # FIXME
      require_data=True,
      create_fmodel=True,
      prefer_anomalous=None,
      force_non_anomalous=False,
      set_wavelength_from_model_header=False,
      set_inelastic_form_factors=None,
      usage_string=None,
      create_log_buffer=False) :
    import mmtbx.monomer_library.pdb_interpretation
    import mmtbx.monomer_library.server
    import mmtbx.utils
    from iotbx import crystal_symmetry_from_any
    from iotbx import file_reader
    import iotbx.phil
    if isinstance(master_phil, str) :
      master_phil = iotbx.phil.parse(master_phil)
    if (usage_string is not None) :
      if (len(args) == 0) or ("--help" in args) :
        raise Usage("""%s\n\nFull parameters:\n%s""" % (usage_string,
          master_phil.as_str(prefix="  ")))
    if (force_non_anomalous) :
      assert (not prefer_anomalous)
    assert (set_inelastic_form_factors in [None, "sasaki", "henke"])
    self.args = args
    self.master_phil = master_phil
    self.processed_pdb_file = self.pdb_inp = None
    self.geometry = None
    self.sequence = None
    self.fmodel = None
    self.f_obs = None
    self.r_free_flags = None
    self.intensity_flag = None
    self.raw_data = None
    self.raw_flags = None
    self.test_flag_value = None
    self.miller_arrays = None
    self.hl_coeffs = None
    self.cif_objects = []
    self.log = out
    if ("--quiet" in args) or ("quiet=True" in args) :
      self.log = null_out()
    elif create_log_buffer :
      self.log = multi_out()
      self.log.register(label="stdout", file_object=out)
      self.log.register(label="log_buffer", file_object=StringIO())
    make_header("Collecting inputs", out=self.log)
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="input.pdb.file_name",
      reflection_file_def="input.xray_data.file_name",
      cif_file_def="input.monomers.file_name",
      seq_file_def="input.sequence")
    self.working_phil = cmdline.work
    params = self.working_phil.extract()
    if len(params.input.pdb.file_name) == 0 :
      raise Sorry("At least one PDB file is required as input.")
    self.cif_file_names = params.input.monomers.file_name
    self.pdb_file_names = params.input.pdb.file_name
    # SYMMETRY HANDLING - PDB FILES
    self.crystal_symmetry = pdb_symm = None
    for pdb_file_name in params.input.pdb.file_name :
      pdb_symm = crystal_symmetry_from_any.extract_from(pdb_file_name)
      if (pdb_symm is not None) :
        break
    # DATA INPUT
    data_and_flags = hkl_symm = hkl_in = None
    if (params.input.xray_data.file_name is None) :
      if (require_data) :
        raise Sorry("At least one reflections file is required as input.")
    else :
      # FIXME this may still require that the data file has full crystal
      # symmetry defined (although for MTZ input this will not be a problem)
      make_sub_header("Processing X-ray data", out=self.log)
      hkl_in = file_reader.any_file(params.input.xray_data.file_name)
      hkl_in.check_file_type("hkl")
      hkl_server = hkl_in.file_server
      symm = hkl_server.miller_arrays[0].crystal_symmetry()
      if ((symm is None) or
          (symm.space_group() is None) or
          (symm.unit_cell() is None)) :
        if (pdb_symm is not None) :
          from iotbx.reflection_file_utils import reflection_file_server
          print >> self.log, \
            "No symmetry in X-ray data file - using PDB symmetry:"
          pdb_symm.show_summary(f=out, prefix="  ")
          hkl_server = reflection_file_server(
            crystal_symmetry=pdb_symm,
            reflection_files=[hkl_in.file_object])
        else :
          raise Sorry("No crystal symmetry information found in input files.")
      if (hkl_server is None) :
        hkl_server = hkl_in.file_server
      data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server=hkl_server,
        parameters=params.input.xray_data,
        data_parameter_scope="input.xray_data",
        flags_parameter_scope="input.xray_data.r_free_flags",
        prefer_anomalous=prefer_anomalous,
        force_non_anomalous=force_non_anomalous,
        log=self.log)
      self.intensity_flag = data_and_flags.intensity_flag
      self.raw_data = data_and_flags.raw_data
      self.raw_flags = data_and_flags.raw_flags
      self.test_flag_value = data_and_flags.test_flag_value
      self.f_obs = data_and_flags.f_obs
      self.r_free_flags = data_and_flags.r_free_flags
      self.miller_arrays = hkl_in.file_server.miller_arrays
      hkl_symm = self.raw_data.crystal_symmetry()
    if len(self.cif_file_names) > 0 :
      for file_name in self.cif_file_names :
        cif_obj = mmtbx.monomer_library.server.read_cif(file_name=file_name)
        self.cif_objects.append((file_name, cif_obj))
    # SYMMETRY HANDLING - COMBINED
    if (hkl_symm is not None) :
      use_symmetry = hkl_symm
    from iotbx.symmetry import combine_model_and_data_symmetry
    self.crystal_symmetry = combine_model_and_data_symmetry(
      model_symmetry=pdb_symm,
      data_symmetry=hkl_symm)
    if (self.crystal_symmetry is not None) and (self.f_obs is not None) :
      self.f_obs = self.f_obs.customized_copy(
        crystal_symmetry=self.crystal_symmetry).eliminate_sys_absent().set_info(
          self.f_obs.info())
      self.r_free_flags = self.r_free_flags.customized_copy(
        crystal_symmetry=self.crystal_symmetry).eliminate_sys_absent().set_info(
          self.r_free_flags.info())
    # EXPERIMENTAL PHASES
    target_name = "ml"
    if hasattr(params.input, "experimental_phases") :
      flag = params.input.use_experimental_phases
      if (flag in [True, Auto]) :
        phases_file = params.input.experimental_phases.file_name
        if (phases_file is None) :
          phases_file = params.input.xray_data.file_name
          phases_in = hkl_in
        else :
          phases_in = file_reader.any_file(phases_file)
          phases_in.check_file_type("hkl")
        phases_in.file_server.err = self.log # redirect error output
        space_group = self.crystal_symmetry.space_group()
        point_group = space_group.build_derived_point_group()
        hl_coeffs = mmtbx.utils.determine_experimental_phases(
          reflection_file_server = phases_in.file_server,
          parameters             = params.input.experimental_phases,
          log                    = self.log,
          parameter_scope        = "input.experimental_phases",
          working_point_group    = point_group,
          symmetry_safety_check  = True)
        if (hl_coeffs is not None) :
          hl_coeffs = hl_coeffs.map_to_asu()
          if hl_coeffs.anomalous_flag() :
            if (not self.f_obs.anomalous_flag()) :
              hl_coeffs = hl_coeffs.average_bijvoet_mates()
          elif self.f_obs.anomalous_flag() :
            hl_coeffs = hl_coeffs.generate_bijvoet_mates()
          self.hl_coeffs = hl_coeffs.matching_set(other=self.f_obs,
            data_substitute=(0,0,0,0))
          target_name = "mlhl"
    # PDB INPUT
    if process_pdb_file :
      pdb_interp_params = getattr(params, "pdb_interpretation", None)
      if (pdb_interp_params is None) :
        pdb_interp_params = \
          mmtbx.monomer_library.pdb_interpretation.master_params.extract()
      make_sub_header("Processing PDB file(s)", out=self.log)
      pdb_combined = mmtbx.utils.combine_unique_pdb_files(
        file_names=params.input.pdb.file_name)
      pdb_combined.report_non_unique(out=self.log)
      pdb_raw_records = pdb_combined.raw_records
      processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
        cif_objects=self.cif_objects,
        pdb_interpretation_params=pdb_interp_params,
        crystal_symmetry=self.crystal_symmetry,
        use_neutron_distances=params.input.scattering_table=="neutron",
        stop_for_unknowns=getattr(pdb_interp_params, "stop_for_unknowns",False),
        log=self.log)
      self.processed_pdb_file, self.pdb_inp = \
        processed_pdb_files_srv.process_pdb_files(
          raw_records = pdb_raw_records,
          stop_if_duplicate_labels = False,
          allow_missing_symmetry=\
            (self.crystal_symmetry is None) and (not require_data))
      self.geometry = self.processed_pdb_file.geometry_restraints_manager(
        show_energies=False)
      assert (self.geometry is not None)
      self.xray_structure = self.processed_pdb_file.xray_structure()
      chain_proxies = self.processed_pdb_file.all_chain_proxies
      self.pdb_hierarchy = chain_proxies.pdb_hierarchy
      # FIXME note that the cdl parameter is always included in the parameters
      # for pdb_interpretation regardless of whether it is actually used by
      # the app in question
      if use_conformation_dependent_library :
        if pdb_interp_params.cdl :
          from mmtbx.conformation_dependent_library import setup_restraints
          from mmtbx.conformation_dependent_library import update_restraints
          from mmtbx import restraints
          restraints_manager = restraints.manager(geometry=self.geometry)
          cdl_proxies = setup_restraints(restraints_manager)
          update_restraints(
            hierarchy=self.pdb_hierarchy,
            restraints_manager=restraints_manager,
            current_geometry=self.xray_structure,
            cdl_proxies=cdl_proxies,
            log=self.log,
            verbose=True)
    else :
      pdb_file_object = mmtbx.utils.pdb_file(
        pdb_file_names=params.input.pdb.file_name,
        cif_objects=self.cif_objects,
        crystal_symmetry=self.crystal_symmetry,
        log=self.log)
      self.pdb_inp = pdb_file_object.pdb_inp
      self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()
      self.pdb_hierarchy.atoms().reset_i_seq()
      self.xray_structure = self.pdb_inp.xray_structure_simple(
        crystal_symmetry=self.crystal_symmetry)
    # wavelength
    if (set_wavelength_from_model_header and params.input.wavelength is None) :
      wavelength = self.pdb_inp.extract_wavelength()
      if (wavelength is not None) :
        print >> self.log, ""
        print >> self.log, "Using wavelength = %g from PDB header" % wavelength
        params.input.wavelength = wavelength
    # set scattering table
    if (data_and_flags is not None) :
      self.xray_structure.scattering_type_registry(
        d_min=self.f_obs.d_min(),
        table=params.input.scattering_table)
      if ((params.input.wavelength is not None) and
          (set_inelastic_form_factors is not None)) :
        self.xray_structure.set_inelastic_form_factors(
          photon=params.input.wavelength,
          table=set_inelastic_form_factors)
      make_sub_header("xray_structure summary", out=self.log)
      self.xray_structure.scattering_type_registry().show(out = self.log)
      self.xray_structure.show_summary(f=self.log)
    # FMODEL SETUP
    if (create_fmodel) and (data_and_flags is not None) :
      make_sub_header("F(model) initialization", out=self.log)
      skip_twin_detection = getattr(params.input, "skip_twin_detection", None)
      twin_law = getattr(params.input, "twin_law", None)
      if (twin_law is Auto) :
        if (self.hl_coeffs is not None) :
          raise Sorry("Automatic twin law determination not supported when "+
            "experimental phases are used.")
      elif (skip_twin_detection is not None) :
        twin_law = Auto
      if (twin_law is Auto) :
        print >> self.log, "Twinning will be detected automatically."
        self.fmodel = mmtbx.utils.fmodel_simple(
          update_f_part1_for=update_f_part1_for,
          xray_structures=[self.xray_structure],
          scattering_table=params.input.scattering_table,
          f_obs=self.f_obs,
          r_free_flags=self.r_free_flags,
          skip_twin_detection=skip_twin_detection,
          target_name=target_name,
          log=self.log)
      else :
        if ((twin_law is not None) and (self.hl_coeffs is not None)) :
          raise Sorry("Automatic twin law determination not supported when "+
            "experimental phases are used.")
        self.fmodel = mmtbx.utils.fmodel_manager(
          f_obs=self.f_obs,
          xray_structure=self.xray_structure,
          r_free_flags=self.r_free_flags,
          twin_law=params.input.twin_law,
          hl_coeff=self.hl_coeffs,
          target_name=target_name)
        self.fmodel.update_all_scales(
          params=None,
          log=self.log,
          optimize_mask=True,
          show=True,
          update_f_part1_for=update_f_part1_for)
      self.fmodel.info().show_rfactors_targets_scales_overall(out=self.log)
    # SEQUENCE
    if (params.input.sequence is not None) :
      seq_file = file_reader.any_file(params.input.sequence,
        force_type="seq",
        raise_sorry_if_errors=True)
      self.sequence = seq_file.file_object
    # UNMERGED DATA
    self.unmerged_i_obs = None
    if hasattr(params.input, "unmerged_data") :
      if (params.input.unmerged_data.file_name is not None) :
        self.unmerged_i_obs = load_and_validate_unmerged_data(
          f_obs=self.f_obs,
          file_name=params.input.unmerged_data.file_name,
          data_labels=params.input.unmerged_data.labels,
          log=self.log)
    self.params = params
    print >> self.log, ""
    print >> self.log, "End of input processing"

  def start_log_file (self, file_name) :
    assert type(self.log).__name__ == 'multi_out'
    log_file = open(file_name, "w")
    self.log.replace_stringio(
      old_label="log_buffer",
      new_label="log",
      new_file_object=log_file)
    return self.log

def load_and_validate_unmerged_data (f_obs, file_name, data_labels,
    log=sys.stdout) :
  from iotbx import merging_statistics
  unmerged_i_obs = merging_statistics.select_data(
    file_name=file_name,
    data_labels=data_labels,
    log=log)
  if ((unmerged_i_obs.space_group() is not None) and
      (unmerged_i_obs.unit_cell() is not None)) :
    if (not unmerged_i_obs.is_similar_symmetry(f_obs)) :
      show_symmetry_error("Data file", "Unmerged data", unmerged_i_obs, f_obs)
  return unmerged_i_obs

def show_symmetry_error (file1, file2, symm1, symm2) :
  import cStringIO
  symm_out1 = cStringIO.StringIO()
  symm_out2 = cStringIO.StringIO()
  symm1.show_summary(f=symm_out1, prefix="  ")
  symm2.show_summary(f=symm_out2, prefix="  ")
  raise Sorry("Incompatible symmetry definitions:\n%s:\n%s\n%s\n%s" %
    (file1, symm_out1.getvalue(), file2, symm_out2.getvalue()))

def validate_input_params (params) :
  if params.input.pdb.file_name is None :
    raise Sorry("No PDB file defined.")
  elif params.input.xray_data.file_name is None :
    raise Sorry("No reflection file defined.")
  elif params.input.xray_data.labels is None :
    raise Sorry("No labels chosen for reflection data.")
  elif (params.input.xray_data.r_free_flags.label is None) :
    raise Sorry("R-free flags not defined.  If you are trying to run this "+
      "program with a reflections file that is missing R-free flags, use "+
      "the reflection file editor to generate a new tests set.")
  return True
