
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
from libtbx.utils import Sorry, Usage
from libtbx import Auto
import sys

cmdline_input_phil_base_str = """
input {
  include scope mmtbx.utils.xray_data_str
  %(phases)s
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
    enable_automatic_twin_detection=True,
    enable_experimental_phases=False,
    enable_pdb_interpretation_params=False) :
  import iotbx.phil
  phil_extra_dict = {
    "phases" : "",
    "phases_flag" : "",
    "twin_law" : "",
    "pdb_interpretation" : "",
  }
  if (enable_automatic_twin_detection) :
    phil_extra_dict["twin_law"] = """
      skip_twin_detection = False
        .type = bool"""
  else :
    phil_extra_dict["twin_law"] = """
      twin_law = None
        .type = str"""
  if (enable_experimental_phases) :
    phil_extra_dict["phases"] = """
      experimental_phases {
        file_name = None
          .type = path
        labels = None
          .type = str
      }"""
    phil_extra_dict["phases_flag"] = """
      use_experimental_phases = Auto
        .type = bool
      """
  if (enable_pdb_interpretation_params) :
    phil_extra_dict["pdb_interpretation"] = """
      mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
      """
  cmdline_input_phil_str = cmdline_input_phil_base_str % phil_extra_dict
  master_phil_str = """
    %s
    %s
  """ % (cmdline_input_phil_str, phil_string)
  return iotbx.phil.parse(master_phil_str, process_includes=True)

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
      require_data=True,
      create_fmodel=True,
      prefer_anomalous=None,
      force_non_anomalous=False,
      usage_string=None) :
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
    self.args = args
    self.master_phil = master_phil
    self.processed_pdb_file = None
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
    if ("--quiet" in args) or ("quiet=True" in args) :
      out = null_out()
    make_header("Collecting inputs", out=out)
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
    # DATA INPUT
    data_and_flags = hkl_symm = hkl_in = None
    if (params.input.xray_data.file_name is None) :
      if (require_data) :
        raise Sorry("At least one reflections file is required as input.")
    else :
      # FIXME this may still require that the data file has full crystal
      # symmetry defined (although for MTZ input this will not be a problem)
      make_sub_header("Processing X-ray data", out=out)
      hkl_in = file_reader.any_file(params.input.xray_data.file_name)
      hkl_in.check_file_type("hkl")
      data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server=hkl_in.file_server,
        parameters=params.input.xray_data,
        data_parameter_scope="input.xray_data",
        flags_parameter_scope="input.xray_data.r_free_flags",
        prefer_anomalous=prefer_anomalous,
        force_non_anomalous=force_non_anomalous,
        log=out)
      self.intensity_flag = data_and_flags.intensity_flag
      self.raw_data = data_and_flags.raw_data
      self.raw_flags = data_and_flags.raw_flags
      self.test_flag_value = data_and_flags.test_flag_value
      self.f_obs = data_and_flags.f_obs
      self.r_free_flags = data_and_flags.r_free_flags
      self.miller_arrays = hkl_in.file_server.miller_arrays
      hkl_symm = self.raw_data.crystal_symmetry()
    self.cif_file_names = params.input.monomers.file_name
    self.pdb_file_names = params.input.pdb.file_name
    if len(self.cif_file_names) > 0 :
      for file_name in self.cif_file_names :
        cif_obj = mmtbx.monomer_library.server.read_cif(file_name=file_name)
        self.cif_objects.append((file_name, cif_obj))
    # SYMMETRY HANDLING
    use_symmetry = pdb_symm = None
    if (hkl_symm is not None) :
      use_symmetry = hkl_symm
    for pdb_file_name in params.input.pdb.file_name :
      pdb_symm = crystal_symmetry_from_any.extract_from(pdb_file_name)
      if (pdb_symm is not None) :
        break
    from iotbx.symmetry import combine_model_and_data_symmetry
    use_symmetry = combine_model_and_data_symmetry(
      model_symmetry=pdb_symm,
      data_symmetry=hkl_symm)
    if (use_symmetry is not None) and (self.f_obs is not None) :
      self.f_obs = self.f_obs.customized_copy(
        crystal_symmetry=use_symmetry).eliminate_sys_absent().set_info(
          self.f_obs.info())
      self.r_free_flags = self.r_free_flags.customized_copy(
        crystal_symmetry=use_symmetry).eliminate_sys_absent().set_info(
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
        point_group = use_symmetry.space_group().build_derived_point_group()
        hl_coeffs = mmtbx.utils.determine_experimental_phases(
          reflection_file_server = phases_in.file_server,
          parameters             = params.input.experimental_phases,
          log                    = out,
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
      make_sub_header("Processing PDB file(s)", out=out)
      pdb_combined = mmtbx.utils.combine_unique_pdb_files(
        file_names=params.input.pdb.file_name)
      pdb_combined.report_non_unique(out=out)
      pdb_raw_records = pdb_combined.raw_records
      processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
        cif_objects=self.cif_objects,
        pdb_interpretation_params=pdb_interp_params,
        crystal_symmetry=use_symmetry,
        use_neutron_distances=params.input.scattering_table=="neutron",
        stop_for_unknowns=getattr(pdb_interp_params, "stop_for_unknowns",False),
        log=out)
      self.processed_pdb_file, self.pdb_inp = \
        processed_pdb_files_srv.process_pdb_files(
          raw_records = pdb_raw_records,
          stop_if_duplicate_labels = False,
          allow_missing_symmetry=(use_symmetry is None) and (not require_data))
      self.geometry = self.processed_pdb_file.geometry_restraints_manager(
        show_energies=False)
      assert (self.geometry is not None)
      self.xray_structure = self.processed_pdb_file.xray_structure()
      chain_proxies = self.processed_pdb_file.all_chain_proxies
      self.pdb_hierarchy = chain_proxies.pdb_hierarchy
    else :
      pdb_file_object = pdb_file(
        pdb_file_names=params.input.pdb.file_name,
        cif_objects=self.cif_objects,
        crystal_symmetry=use_symmetry,
        log=out)
      pdb_in = pdb_file_object.pdb_inp
      self.pdb_hierarchy = pdb_in.construct_hierarchy()
      self.pdb_hierarchy.atoms().reset_i_seq()
      self.xray_structure = pdb_in.xray_structure_simple(
        crystal_symmetry=use_symmetry)
    # set scattering table
    if (data_and_flags is not None) :
      self.xray_structure.scattering_type_registry(
        d_min=self.f_obs.d_min(),
        table=params.input.scattering_table)
      make_sub_header("xray_structure summary", out=out)
      self.xray_structure.scattering_type_registry().show(out = out)
      self.xray_structure.show_summary(f=out)
    # FMODEL SETUP
    if (create_fmodel) and (data_and_flags is not None) :
      skip_twin_detection = getattr(params.input, "skip_twin_detection", None)
      twin_law = getattr(params.input, "twin_law", None)
      if (twin_law is Auto) :
        if (self.hl_coeffs is not None) :
          raise Sorry("Automatic twin law determination not supported when "+
            "experimental phases are used.")
      elif (skip_twin_detection is not None) :
        twin_law = Auto
      if (twin_law is Auto) :
        self.fmodel = mmtbx.utils.fmodel_simple(
          update_f_part1_for=update_f_part1_for,
          xray_structures=[self.xray_structure],
          scattering_table=params.input.scattering_table,
          f_obs=self.f_obs,
          r_free_flags=self.r_free_flags,
          skip_twin_detection=skip_twin_detection,
          target_name=target_name,
          log=out)
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
          log=out,
          optimize_mask=True,
          update_f_part1_for=update_f_part1_for)
    # SEQUENCE
    if (params.input.sequence is not None) :
      seq_file = file_reader.any_file(params.input.sequence,
        force_type="seq",
        raise_sorry_if_errors=True)
      self.sequence = seq_file.file_object
    self.params = params
    print >> out
    print >> out, "End of input processing"

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
