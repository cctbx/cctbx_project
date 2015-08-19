
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
  energy = None
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
    enable_cdl=None,
    as_phil_string=False) :
  """
  Generate a complete PHIL parameter block with generic input parameters plus
  user-specified options.  The result is suitable for input for the class
  load_model_and_data.  Depending on the target application, the exact input
  options can be adjusted.

  :param phil_string: application-specific parameters
  :param enable_twin_law: allow twinned f_model calculation
  :param enable_automatic_twin_detection: allow automatic detection of twinning
      and setup of the f_model object
  :param enable_experimental_phases: use Hendrickson-Lattman coefficients
  :param enable_pdb_interpretation_params: show options for modifying the
      behavior of mmtbx.monomer_library.pdb_interpretation
  :param enable_stop_for_unknowns: modify behavior of restraint interpretation
      when unknown atoms are encountered.  Default is None; if True, the
      program will not raise an error; if False, the program will raise an
      error which may be suppressed by the user
  :param enable_full_geometry_params: include parameters for specifying custom
      geometry resraints
  :param enable_unmerged_data: accept separate unmerged intensities
  :param enable_cdl: change default setting for conformation-dependent library
  :param as_phil_string: return parameter string instead of PHIL object
  :returns: PHIL object (unless as_phil_string=True)
  """
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
    assert (enable_cdl is None)
    return master_phil_str
  master_phil = iotbx.phil.parse(master_phil_str, process_includes=True)
  if (enable_cdl is not None) :
    wp = iotbx.phil.parse("pdb_interpretation.restraints_library.cdl=%s" % enable_cdl)
    master_phil = master_phil.fetch(source=wp)
  return master_phil

def generic_simple_input_phil () :
  """
  Generate minimal PHIL input string with no additional parameters.
  """
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

  Parameters
  ----------
  args: list of command-line arguments
  master_phil: PHIL master (can optionally be an unparsed string)
  out: filehandle-like object
  process_pdb_file: run full restraints generation
  require_data: raise error if no experimental data supplied
  create_fmodel: setup mmtbx.f_model.manager object
  prefer_anomalous: preferentially use anomalous data if present
  force_non_anomalous: merge anomalous data if present
  set_wavelength_from_model_header: interpret PDB or mmCIF header to set \
    experimental wavelength
  set_inelastic_form_factors: table to use (if any) for setting anomalous \
    scattering form factors
  usage_string: console output for no arguments or --help
  create_log_buffer: store log output for later output to file
  remove_unknown_scatterers: delete atoms with scattering type 'X' (only \
    used when process_pdb_file=False)
  generate_input_phil: specifies that the master_phil object is a string \
    containing onky the app-specific options, and automatically add the \
    standard input parameters

  Attributes
  ----------
  args : list of str
  cif_file_names : libtbx.phil.scope_extract_list
  cif_objects : list of ...
  crystal_symmetry : cctbx.crystal.symmetry
  f_obs : cctbx.miller.array
  fmodel : mmtbx.f_model.manager
  geometry : cctbx.geometry_restraints.manager.manager
  hl_coeffs : ...
  intensity_flag : bool
  log : file
  master_phil : libtbx.phil.scope
  miller_arrays : list of cctbx.miller.array
  params : libtbx.phil.scope_extract
  pdb_file_names : libtbx.phil.scope_extract_list
  pdb_hierarchy : iotbx.pdb.hierarchy.root
  pdb_inp : iotbx.pdb.input
  processed_pdb_file : mmtbx.monomer_library.pdb_interpretation.process
  r_free_flags : cctbx.miller.array
  raw_data : cctbx.miller.array
  raw_flags : cctbx.miller.array
  sequence : ...
  test_flag_value : int
  unknown_residues_flag : bool
  unknown_residues_error_message : str
  unmerged_i_obs : ...
  working_phil : libtbx.phil.scope
  xray_structure : cctbx.xray.structure.structure

  Examples
  --------
  >>> from mmtbx.command_line import load_model_and_data
  >>> cmdline = load_model_and_data(
  ...  args=["model.pdb", "data.mtz"],
  ...  master_phil=master_phil,
  ...  prefer_anomalous=True,
  ...  set_wavelength_from_model_header=True,
  ...  set_inelastic_form_factors="sasaki",
  ...  )
  >>> assert cmdline.pdb_hierarchy is not None
  >>> assert cmdline.fmodel is not None
  >>> assert cmdline.params.input.wavelength is not None
  """
  def __init__ (self,
      args,
      master_phil,
      out=sys.stdout,
      process_pdb_file=True,
      require_data=True,
      create_fmodel=True,
      prefer_anomalous=None,
      force_non_anomalous=False,
      set_wavelength_from_model_header=False,
      set_inelastic_form_factors=None,
      usage_string=None,
      create_log_buffer=False,
      remove_unknown_scatterers=False,
      generate_input_phil=False) :
    import mmtbx.monomer_library.pdb_interpretation
    import mmtbx.monomer_library.server
    import mmtbx.utils
    from iotbx import crystal_symmetry_from_any
    from iotbx import file_reader
    import iotbx.phil
    if generate_input_phil :
      assert isinstance(master_phil, basestring)
      master_phil = generate_master_phil_with_inputs(phil_string=master_phil)
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
    self.pdb_hierarchy = self.xray_structure = None
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
    self.unknown_residues_flag = False
    self.unknown_residues_error_message = False
    if process_pdb_file :
      pdb_interp_params = getattr(params, "pdb_interpretation", None)
      if (pdb_interp_params is None) :
        pdb_interp_params = \
          mmtbx.monomer_library.pdb_interpretation.master_params.extract()
      make_sub_header("Processing PDB file(s)", out=self.log)
      pdb_combined = mmtbx.utils.combine_unique_pdb_files(
        file_names=params.input.pdb.file_name,)
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
      error_msg = self.processed_pdb_file.all_chain_proxies.\
        fatal_problems_message(
          ignore_unknown_scattering_types=False,
          ignore_unknown_nonbonded_energy_types=False)
      if (error_msg is not None) :
        self.unknown_residues_flag = True
        self.unknown_residues_error_message = error_msg
      self.geometry = self.processed_pdb_file.geometry_restraints_manager(
        show_energies=False)
      assert (self.geometry is not None)
      self.xray_structure = self.processed_pdb_file.xray_structure()
      chain_proxies = self.processed_pdb_file.all_chain_proxies
      self.pdb_hierarchy = chain_proxies.pdb_hierarchy
    else :
      pdb_file_object = mmtbx.utils.pdb_file(
        pdb_file_names=params.input.pdb.file_name,
        cif_objects=self.cif_objects,
        crystal_symmetry=self.crystal_symmetry,
        log=self.log)
      self.pdb_inp = pdb_file_object.pdb_inp
      self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()
      if (remove_unknown_scatterers) :
        known_sel = self.pdb_hierarchy.atom_selection_cache().selection(
          "not element X")
        if (known_sel.count(True) != len(known_sel)) :
          self.pdb_hierarchy = self.pdb_hierarchy.select(known_sel)
          self.xray_structure = self.pdb_hierarchy.extract_xray_structure(
            crystal_symmetry=self.crystal_symmetry)
      self.pdb_hierarchy.atoms().reset_i_seq()
      if (self.xray_structure is None) :
        self.xray_structure = self.pdb_inp.xray_structure_simple(
          crystal_symmetry=self.crystal_symmetry)
    # wavelength
    if (params.input.energy is not None) :
      if (params.input.wavelength is not None) :
        raise Sorry("Both wavelength and energy have been specified!")
      params.input.wavelength = 12398.424468024265 / params.input.energy
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
          show=True)
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
    """
    Open a log file and write out the existing output buffer, returning the
    multi_out pseudo-filehandle.

    :param file_name: log file to create
    :returns: libtbx.utils.multi_out object
    """
    assert type(self.log).__name__ == 'multi_out'
    log_file = open(file_name, "w")
    self.log.replace_stringio(
      old_label="log_buffer",
      new_label="log",
      new_file_object=log_file)
    return self.log

  def save_data_mtz (self, file_name) :
    """
    Write the processed amplitudes, optional Hendrickson-Lattman coefficients,
    and R-free flags to the designated MTZ file.
    """
    assert (self.f_obs is not None)
    mtz_data = self.f_obs.as_mtz_dataset(column_root_label="F")
    if (self.hl_coeffs is not None) :
      mtz_data.add_miller_array(self.hl_coeffs,
        column_root_label="HL")
    if (self.r_free_flags is not None) :
      mtz_data.add_miller_array(self.r_free_flags,
        column_root_label="FreeR_flag")
    mtz_data.mtz_object().write(file_name)

  def create_model_manager (self, log=None) :
    """
    Instantiate an mmtbx.model.manager object with the current pdb hierarchy,
    xray structure, and geometry restraints.
    """
    if (log is None) : log = self.log
    import mmtbx.restraints
    import mmtbx.model
    restraints_manager = mmtbx.restraints.manager(
      geometry=self.geometry,
      normalization=True)
    return mmtbx.model.manager(
      xray_structure=self.xray_structure,
      pdb_hierarchy=self.pdb_hierarchy,
      restraints_manager=restraints_manager,
      log=log)

def load_and_validate_unmerged_data (f_obs, file_name, data_labels,
    log=sys.stdout) :
  """
  Read in (and verify) unmerged intensities, e.g. from scalepack or XDS.
  """
  from iotbx import merging_statistics
  unmerged_i_obs = merging_statistics.select_data(
    file_name=file_name,
    data_labels=data_labels,
    log=log)
  if ((unmerged_i_obs.space_group() is not None) and
      (unmerged_i_obs.unit_cell() is not None)) :
    if (not unmerged_i_obs.is_similar_symmetry(f_obs)) :
      pg_f_obs = f_obs.space_group().build_derived_point_group()
      pg_i_obs = unmerged_i_obs.space_group().build_derived_point_group()
      if (pg_i_obs == pg_f_obs) :
        # special case: same unit cell, same point group, different space group
        if unmerged_i_obs.unit_cell().is_similar_to(f_obs.unit_cell()) :
          return unmerged_i_obs
      show_symmetry_error("Data file", "Unmerged data", unmerged_i_obs, f_obs)
  elif (unmerged_i_obs.space_group() is not None) :
    pg_f_obs = f_obs.space_group().build_derived_point_group()
    pg_i_obs = unmerged_i_obs.space_group().build_derived_point_group()
    if (pg_i_obs != pg_f_obs) :
      raise Sorry("Incompatible space groups in merged and unmerged data:"+
        "%s versus %s" % (f_obs.space_group_info(),
        unmerged_i_obs.space_group_info()))
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
  """
  Check for completeness of mandatory input parameters
  """
  if params.input.pdb.file_name is None :
    raise Sorry("No PDB file defined.")
  elif isinstance(params.input.pdb.file_name,list):
    if (len(params.input.pdb.file_name) == 0):
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
