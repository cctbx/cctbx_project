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
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
from iotbx import pdb
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from scitbx.math import matrix
from cctbx import adptbx
import sys, os, math
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import iotbx.phil
import libtbx.phil.command_line
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.option_parser import iotbx_option_parser
from iotbx import crystal_symmetry_from_any
from iotbx import pdb
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import mtz
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import user_plus_sys_time, show_total_time
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, date_and_time, host_and_user, multi_out
from libtbx import adopt_init_args
import random, sys, os
from libtbx.test_utils import approx_equal
from mmtbx.refinement import print_statistics
import libtbx.load_env
from mmtbx import utils
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from iotbx.pdb import combine_unique_pdb_files


map_params_str ="""\
  map_format = *xplor
    .optional = True
    .type = choice(multi=True)
  map_coefficients_format = *mtz phs
    .optional = True
    .type = choice(multi=True)
  suppress = None
    .type = strings
    .help = List of mtz_label_amplitudes of maps to be suppressed. \
            Intended to selectively suppress computation and \
            writing of the standard maps.
    .expert_level = 1
  map
    .multiple = True
    .short_caption=Electron density map
    .style=noauto auto_align scrolled
  {
    mtz_label_amplitudes = None
      .type = str
      .short_caption=Amplitude label
    mtz_label_phases = None
      .type = str
      .short_caption=Phase label
    likelihood_weighted = None
      .type = bool
      .expert_level=1
    obs_factor = None
      .type = float
      .short_caption=Multiply Fobs by
    calc_factor = None
      .type = float
      .short_caption=Multiply Fcalc by
    kicked = False
      .type = bool
    fill_missing_f_obs_with_weighted_f_model = True
      .type = bool
  }
  map {
    mtz_label_amplitudes = 2FOFCWT
    mtz_label_phases = PH2FOFCWT
    likelihood_weighted = True
    obs_factor = 2
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = True
  }
  map {
    mtz_label_amplitudes = FOFCWT
    mtz_label_phases = PHFOFCWT
    likelihood_weighted = True
    obs_factor = 1
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = True
  }
  map {
    mtz_label_amplitudes = 2FOFCWT_no_fill
    mtz_label_phases = PH2FOFCWT_no_fill
    likelihood_weighted = True
    obs_factor = 2
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = False
  }
  map {
    mtz_label_amplitudes = FOFCWT_no_fill
    mtz_label_phases = PHFOFCWT_no_fill
    likelihood_weighted = True
    obs_factor = 1
    calc_factor = 1
    kicked = False
    fill_missing_f_obs_with_weighted_f_model = False
  }
  anomalous_difference_map
    .short_caption=Anomalous difference map
    .style = box auto_align
  {
    mtz_label_amplitudes = ANOM
      .type = str
      .short_caption=Amplitude label
    mtz_label_phases = PHANOM
      .type = str
      .short_caption=Phase label
  }

  grid_resolution_factor = 1/4
    .type = float
    .expert_level=1
  region = *selection cell
    .type = choice
    .expert_level=1
    .short_caption=Map region
  atom_selection = None
    .type = str
    .expert_level=2
    .style = selection
  atom_selection_buffer = 3
    .type = float
    .expert_level=2
  apply_sigma_scaling = True
    .type = bool
    .expert_level = 1
  apply_volume_scaling = False
    .type = bool
    .expert_level = 2
"""

map_params = iotbx.phil.parse(map_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
maps {
  %s
  crystal_symmetry
    .help = Unit cell and space group parameters
  {
    unit_cell=None
      .type=unit_cell
    space_group=None
      .type=space_group
  }
  pdb_interpretation
    .expert_level = 2
  {
    include scope mmtbx.monomer_library.pdb_interpretation.master_params
  }
  use_experimental_phases = None
    .type = bool
    .expert_level=0
  high_resolution = None
    .type = float
    .expert_level=0
  low_resolution = None
    .type = float
    .expert_level=0
  sigma_fobs_rejection_criterion = 0.0
    .type = float
    .expert_level=0
  neutron_data = False
    .type = bool
    .expert_level=0
  symmetry_safety_check = *error warning
    .type=choice
    .expert_level=2
  bulk_solvent_and_scale
    .expert_level=0
  {
   include scope mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params
  }
  input {
    pdb
    {
      include scope mmtbx.utils.pdb_params
    }
    experimental_phases {
      include scope mmtbx.utils.experimental_phases_params
    }
    monomers {
      include scope mmtbx.utils.cif_params
    }
  }
}
"""%map_params_str, process_includes=True)

def set_f_model(command_line_interpreter):
  f_obs = r_free_flags = command_line_interpreter.data()
  r_free_flags = command_line_interpreter.r_free_flags
  if(r_free_flags is None):
     r_free_flags = f_obs.array(data = flex.bool(flags.size(), False))
  else:
     r_free_flags = f_obs.array(data = flags)

def run(args, command_name="phenix.maps"):
  log = utils.set_log(args)
  utils.print_programs_start_header(
    log  = log,
    text = "  %s: tools for electron density maps calculation" % command_name)
  command_line_interpreter = interpreter(command_name  = command_name,
                                         args          = args,
                                         log           = log)
  raise Sorry("Not implemented yet.")

class interpreter:

  def __init__(self,
        command_name,
        args,
        log,
        flags=None):
    self.command_name = command_name
    self.args = args
    self.log = log
    self.mon_lib_srv = monomer_library.server.server()
    self.ener_lib = monomer_library.server.ener_lib()
    self.params = None
    self.pdb_file_names = []
    self.processed_pdb_file = None
    self.processed_pdb_file_reference = None
    self.reflection_files = []
    self.data = None
    self.neutron_data = None
    self.neutron_r_free_flags = None
    self.r_free_flags = None
    self.experimental_phases = None
    self.output_mtz_file_name = None
    self.process_args()
    self.pdb_file_names.extend(self.params.maps.input.pdb.file_name)
    self.processed_pdb_file, self.pdb_inp = utils.process_pdb_file(
               pdb_file_names            = self.pdb_file_names,
               parameters                = self.params.maps.input.pdb,
               pdb_interpretation_params = self.params.maps.pdb_interpretation,
               mon_lib_srv               = self.mon_lib_srv,
               ener_lib                  = self.ener_lib,
               crystal_symmetry          = self.crystal_symmetry,
               log                       = self.log)
    reflection_file_server = self.reflection_file_server()
    if(not self.params.maps.neutron_data):
       data_source = "X-ray data"
    else:
       data_source = "Neutron data"
    self.data = utils.determine_data(
               reflection_file_server = reflection_file_server,
               parameters             = self.params.maps.input.data,
               parameter_scope        = "maps.input.data",
               log                    = self.log,
               data_description       = data_source,
               ignore_all_zeros       = True,
               working_point_group    = self.point_group,
               symmetry_safety_check  = self.params.maps.symmetry_safety_check)
    self.r_free_flags = utils.determine_r_free_flags(
               reflection_file_server = reflection_file_server,
               data                   = self.data,
               generate_r_free_flags  = False,
               parameters             = self.params.maps.input.r_free_flags,
               parameter_scope        = "maps.input.r_free_flags",
               working_point_group    = self.point_group,
               symmetry_safety_check  = self.params.maps.symmetry_safety_check,
               log                    = self.log,
               neutron_flag           = self.params.maps.neutron_data)
    self.experimental_phases = utils.determine_experimental_phases(
               reflection_file_server = reflection_file_server,
               parameters             =
                 self.params.maps.input.experimental_phases,
               log                    = self.log,
               parameter_scope        = "refinement.input.experimental_phases",
               working_point_group    = self.point_group,
               symmetry_safety_check  = self.params.maps.symmetry_safety_check)

  def process_args(self):
    args = self.args
    if (len(args) == 0): args = ["--help"]
    description_see_also \
        = 'See also: http://www.phenix-online.org/\n' +\
          'Questions / problems: phenixbb@phenix-online.org'
    self.command_line = (iotbx_option_parser(
      usage="%s [options] [reflection_file] [pdb_file] [parameter_file]"
        % self.command_name,
      description='Example: %s data.mtz model.pdb refine.params\n\n'
        % self.command_name + description_see_also)
      .enable_show_defaults()
      .enable_symmetry_comprehensive()
      .option(None, "--unused_ok",
        action="store_true",
        default=False,
        help="Disables detection of unused parameter definitions")
      .option(None, "--quiet",
        action="store_true",
        help="Suppress output to screen")
    ).process(args=args)
    if(self.command_line.expert_level is not None):
      master_params.show(
        out = self.log,
        expert_level = self.command_line.expert_level,
        attributes_level = self.command_line.attributes_level)
      sys.exit(0)
    print >> self.log
    if (len(args) > 0):
      print >> self.log, "Command line arguments:", " ".join([show_string(arg)
        for arg in args])
      print >> self.log
    print >> self.log
    crystal_symmetries_from_coordinate_file = []
    crystal_symmetries_from_reflection_file = []
    cif_objects = []
    parameter_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil = master_params,
      home_scope  = "maps")
    parsed_params = []
    command_line_params = []
    print >> self.log, "Processing inputs. This may take a minute or two."
    for arg in self.command_line.args:
      arg_is_processed = False
      if (os.path.isfile(arg)):
        params = None
        try: params = iotbx.phil.parse(file_name=arg)
        except KeyboardInterrupt: raise
        except RuntimeError: pass
        else:
          if (len(params.objects) == 0):
            params = None
        if (params is not None):
          parsed_params.append(params)
          arg_is_processed = True
        elif (pdb.is_pdb_file(file_name=arg)):
          self.pdb_file_names.append(arg)
          arg_is_processed = True
          try:
            crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
              file_name=arg)
          except KeyboardInterrupt: raise
          except: pass
          else:
            if (crystal_symmetry is not None):
              crystal_symmetries_from_coordinate_file.append(
                crystal_symmetry)
        else:
          try:
            cif_object = mmtbx.monomer_library.server.read_cif(file_name=arg)
          except KeyboardInterrupt: raise
          except: pass
          else:
            if (len(cif_object) > 0):
              cif_objects.append((arg,cif_object))
              arg_is_processed = True
      if (not arg_is_processed):
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if (reflection_file.file_type() is not None):
          self.reflection_files.append(reflection_file)
          arg_is_processed = True
          try:
            crystal_symmetry = crystal_symmetry_from_any.extract_from(
              file_name=arg)
          except KeyboardInterrupt: raise
          except: pass
          else:
            if (crystal_symmetry is not None):
              crystal_symmetries_from_reflection_file.append(crystal_symmetry)
      if (not arg_is_processed):
        try:
          params = parameter_interpreter.process(arg=arg)
        except Sorry, e:
          if (not os.path.isfile(arg)):
            if ("=" in arg): raise
            e.reset_tracebacklimit()
            raise Sorry("File not found: %s" % show_string(arg))
          e.reset_tracebacklimit()
          raise Sorry("Unknown file format: %s" % arg)
        else:
          command_line_params.append(params)
    print >> self.log
    if (len(command_line_params) > 0):
      print >> self.log, "Command line parameter definitions:"
      for params in command_line_params:
        params.show(out=self.log, prefix="  ")
        print >> self.log
    self.params, unused_definitions = master_params.fetch(
      sources = parsed_params+command_line_params,
      track_unused_definitions = True)
    if (len(unused_definitions)):
      print >> self.log, "*"*79
      if (self.command_line.options.unused_ok):
        print >> self.log, "WARNING:",
      else:
        print >> self.log, "ERROR:",
      print >> self.log, "Unused parameter definitions:"
      for obj_loc in unused_definitions:
        print >> self.log, " ", str(obj_loc)
      print >> self.log, "*"*79
      print >> self.log
      if (not self.command_line.options.unused_ok):
        raise Sorry("""Unused parameter definitions:
  Please check the input file(s) for spelling errors and obsolete
  parameter definitions.
  To disable this error message, add
    --unused_ok
  to the command line arguments.""")
    self.params = self.params.extract()
    utils.process_monomer_cif_files(
                                cif_objects  = cif_objects,
                                parameters   = self.params.maps.input.monomers,
                                mon_lib_srv  = self.mon_lib_srv,
                                ener_lib     = self.ener_lib)
    del cif_objects
    self.crystal_symmetry = crystal.select_crystal_symmetry(
      from_command_line = self.command_line.symmetry,
      from_parameter_file=crystal.symmetry(
        unit_cell = self.params.maps.crystal_symmetry.unit_cell,
        space_group_info = self.params.maps.crystal_symmetry.space_group),
      from_coordinate_files=crystal_symmetries_from_coordinate_file,
      from_reflection_files=crystal_symmetries_from_reflection_file)
    if (   self.crystal_symmetry.unit_cell() is None
        or self.crystal_symmetry.space_group_info() is None):
      raise Sorry(
        "Crystal symmetry is not defined. Please use the --symmetry option"
        " or define the refinement.crystal_symmetry parameters.")
    print >> self.log, "Working crystal symmetry after inspecting all inputs:"
    self.crystal_symmetry.show_summary(f = self.log, prefix="  ")
    self.params.maps.crystal_symmetry.unit_cell = \
      self.crystal_symmetry.unit_cell()
    self.params.maps.crystal_symmetry.space_group = \
      self.crystal_symmetry.space_group_info()
    self.point_group = self.crystal_symmetry.space_group() \
      .build_derived_point_group()
    print >> self.log

  def reflection_file_server(self):
    return reflection_file_utils.reflection_file_server(
      crystal_symmetry=self.crystal_symmetry,
      force_symmetry=True,
      reflection_files=self.reflection_files,
      err=self.log)

### Fo-Fo map calculation code
fo_minus_fo_master_params_str = """\
f_obs_1_file_name = None
  .type = str
  .help = File with Fobs data
  .expert_level = 1
f_obs_1_label = None
  .type = str
  .expert_level = 1

f_obs_2_file_name = None
  .type = str
  .help = File with Fobs data
  .expert_level = 1
f_obs_2_label = None
  .type = str
  .expert_level = 1
high_resolution = None
  .type = float
  .help = High resolution data cutoff
  .expert_level = 1
low_resolution = None
  .type = float
  .help = Low resolution data cutoff
  .expert_level = 1
sigma_cutoff = None
  .type = float
  .help = Fobs sigma cutoff
  .expert_level = 1
phase_source = None
  .type = str
  .help = PDB file with a model or reflection file with the phases
  .expert_level = 1
scattering_table = *xray neutron
  .type = choice(multi=False)
  .help = Choices of scattering table for structure factors calculations
  .expert_level=1
"""
def fo_minus_fo_master_params():
  return iotbx.phil.parse(fo_minus_fo_master_params_str, process_includes=False)

def compute_fo_minus_fo_map(data_arrays, xray_structure, log, silent):
  fmodels = []
  for i_seq, d in enumerate(data_arrays):
    if(not silent):
      print >> log, "Data set: %d"%i_seq
    if(d.anomalous_flag()):
      d = d.average_bijvoet_mates()
    r_free_flags = d.array(data = flex.bool(d.data().size(), False))
    fmodel = mmtbx.f_model.manager(
      xray_structure = xray_structure,
      r_free_flags   = r_free_flags,
      target_name    = "ls_wunit_k1",
      f_obs          = d)
    fmodel.remove_outliers()
    fmodel.update_solvent_and_scale()
    if(not silent):
      fmodel.info().show_rfactors_targets_scales_overall()
      print >> log
    fmodels.append(fmodel)
  # prepare Fobs for map calculation (apply scaling):
  f_obss = []
  for fmodel in fmodels:
    obs = fmodel.f_obs
    fb_cart  = fmodel.fb_cart()
    scale_k2 = fmodel.scale_k2()
    f_obs_scale   = 1.0 / fb_cart * scale_k2
    obs = miller.array(miller_set = fmodel.f_model(),
                       data       = obs.data()*f_obs_scale)
    f_obss.append(obs)
  # given two Fobs sets, make them one-to-one matching, get phases and map coefficients
  # Note: f_calc below is just f_calc from atoms (no bulk solvent etc applied)
  fobs_1, f_model = f_obss[0].common_sets(other = fmodels[1].f_model())
  fobs_1, fobs_2 = fobs_1.common_sets(other = f_obss[1])
  fobs_1, f_model = fobs_1.common_sets(other = f_model)
  assert fobs_2.indices().all_eq(fobs_1.indices())
  assert f_model.indices().all_eq(fobs_1.indices())
  # scale again
  scale_k1 = 1
  den = flex.sum(flex.abs(fobs_2.data())*flex.abs(fobs_2.data()))
  if(den != 0):
    scale_k1 = flex.sum(flex.abs(fobs_1.data())*flex.abs(fobs_2.data())) / den
  # map coefficients
  diff = miller.array(
    miller_set = f_model,
    data       = fobs_1.data()-fobs_2.data()*scale_k1)
  def phase_transfer(miller_array, phase_source):
    tmp = miller.array(miller_set = miller_array,
      data = flex.double(miller_array.indices().size(), 1)
      ).phase_transfer(phase_source = phase_source)
    return miller.array(miller_set = miller_array,
      data = miller_array.data() * tmp.data() )
  map_coeff = phase_transfer(
    miller_array = diff,
    phase_source = f_model)
  # output MTZ file with map coefficients
  class map_coeffs_mtz_label_manager:
    def __init__(self, amplitudes, phases):
      self._amplitudes = amplitudes
      self._phases = phases
    def amplitudes(self):
      return self._amplitudes
    def phases(self, root_label, anomalous_sign=None):
      assert anomalous_sign is None or not anomalous_sign
      return self._phases
  mtz_history_buffer = flex.std_string()
  lbl_mgr = map_coeffs_mtz_label_manager(amplitudes = "FoFo", phases = "PHFc")
  if(map_coeff.anomalous_flag()):
    map_coeff = map_coeff.average_bijvoet_mates()
  mtz_dataset = map_coeff.as_mtz_dataset(
    column_root_label=lbl_mgr.amplitudes(),
    label_decorator=lbl_mgr)
  mtz_history_buffer.append("> column label %s = phenix %s" % (
      lbl_mgr.amplitudes(), "FoFoPHFc"))
  file_name = "FoFoPHFc.mtz"
  mtz_history_buffer.append("file name %s"%file_name)
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.add_history(mtz_history_buffer)
  mtz_object.write(file_name=file_name)

def run_fobs_minus_fobs_map(args, command_name = "phenix.fobs_minus_fobs_map"):
  if(len(args) == 0): args = ["--help"]
  examples = """Examples:

phenix.fobs_minus_fobs_map f_obs_1_file=data1.mtz f_obs_2_file=data2.sca \
f_obs_1_label=FOBS1 f_obs_2_label=FOBS2 model.pdb

phenix.fobs_minus_fobs_map f_obs_1_file=data.mtz f_obs_2_file=data.mtz \
f_obs_1_label=FOBS1 f_obs_2_label=FOBS2 phase_source=model.pdb \
high_res=2.0 sigma_cutoff=2 scattering_table=neutron"""

  command_line = (iotbx_option_parser(
    usage="%s [options]" % command_name,
    description=examples)
    .option("--silent",
      action="store_true",
      help="Suppress output to the screen.")
    ).process(args=args)
  #
  log = sys.stdout
  if(not command_line.options.silent):
    utils.print_header("phenix.fobs_minus_fobs_map", out = log)
    print >> log, "Command line arguments: "
    print >> log, args
    print >> log
  #
  processed_args = utils.process_command_line_args(args = command_line.args,
    master_params = fo_minus_fo_master_params(), log = log)
  crystal_symmetry = processed_args.crystal_symmetry
  params = processed_args.params
  #
  if(not command_line.options.silent):
    print >> log, "*** Parameters:"
    params.show(out = log)
    print >> log
  params = params.extract()
  #
  pdb_file_names = processed_args.pdb_file_names
  if(len(processed_args.pdb_file_names) == 0):
    if(params.phase_source is not None):
      pdb_file_names = [params.phase_source]
    else:
      raise Sorry("No PDB file found.")
  # Extaract Fobs1, Fobs2
  f_obss = []
  if(len(processed_args.reflection_files)==2):
    for reflection_file in processed_args.reflection_files:
      reflection_file_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        force_symmetry   = True,
        reflection_files = [reflection_file],
        err              = StringIO())
      determine_data_and_flags_result = utils.determine_data_and_flags(
        reflection_file_server  = reflection_file_server,
        keep_going              = True,
        log                     = StringIO())
      f_obss.append(determine_data_and_flags_result.f_obs)
  else:
    if([params.f_obs_1_file_name,params.f_obs_2_file_name].count(None)==2):
      raise Sorry("No reflection data file found.")
    for file_name, label in zip([params.f_obs_1_file_name,params.f_obs_2_file_name],
                                [params.f_obs_1_label,params.f_obs_2_label]):
      reflection_file = reflection_file_reader.any_reflection_file(
        file_name = file_name, ensure_read_access = False)
      reflection_file_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        force_symmetry   = True,
        reflection_files = [reflection_file],
        err              = StringIO())
      parameters = utils.data_and_flags.extract()
      if(label is not None):
        parameters.labels = [label]
      determine_data_and_flags_result = utils.determine_data_and_flags(
          reflection_file_server  = reflection_file_server,
          parameters              = parameters,
          keep_going              = True,
          log                     = StringIO())
      f_obss.append(determine_data_and_flags_result.f_obs)
  if(len(f_obss)!=2):
    raise Sorry(" ".join(errors))
  if(not command_line.options.silent):
    for ifobs, fobs in enumerate(f_obss):
      print >> log, "*** Summary for data set %d:"%ifobs
      fobs.show_comprehensive_summary(f = log)
      print >> log
  pdb_combined = combine_unique_pdb_files(file_names = pdb_file_names)
  pdb_combined.report_non_unique(out = log)
  if(len(pdb_combined.unique_file_names) == 0):
    raise Sorry("No coordinate file given.")
  xray_structure = iotbx.pdb.input(source_info = None, lines =
    pdb_combined.raw_records).xray_structure_simple()
  if(not command_line.options.silent):
    print >> log, "*** Model summary:"
    xray_structure.show_summary(f = log)
    print >> log
  f_obss[0] = f_obss[0].resolution_filter(d_min = params.high_resolution,
    d_max = params.low_resolution)
  f_obss[1] = f_obss[1].resolution_filter(d_min = params.high_resolution,
    d_max = params.low_resolution)
  if(params.sigma_cutoff is not None):
    for i in [0,1]:
      if(f_obss[i].sigmas() is not None):
        sel = f_obss[i].data() > f_obss[i].sigmas()*params.sigma_cutoff
        f_obss[i] = f_obss[i].select(sel)
  compute_fo_minus_fo_map(data_arrays = f_obss, xray_structure = xray_structure,
    log = log, silent = command_line.options.silent)
