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
import sys, os
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
from scitbx.python_utils.misc import user_plus_sys_time, show_total_time
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, date_and_time, host_and_user, multi_out
from libtbx import adopt_init_args
import random, sys, os
from libtbx.test_utils import approx_equal
from mmtbx.refinement import print_statistics
import libtbx.load_env
from mmtbx import utils

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
  {
    mtz_label_amplitudes = None
      .type = str
    mtz_label_phases = None
      .type = str
    likelihood_weighted = None
      .type = bool
    obs_factor = None
      .type = float
    calc_factor = None
      .type = float
  }
  map {
    mtz_label_amplitudes = 2FOFCWT
    mtz_label_phases = PH2FOFCWT
    likelihood_weighted = True
    obs_factor = 2
    calc_factor = 1
  }
  map {
    mtz_label_amplitudes = FOFCWT
    mtz_label_phases = PHFOFCWT
    likelihood_weighted = True
    obs_factor = 1
    calc_factor = 1
  }
  anomalous_difference_map {
    mtz_label_amplitudes = ANOM
      .type = str
    mtz_label_phases = PHANOM
      .type = str
  }
  grid_resolution_factor = 1/3
    .type = float
  region = *selection cell
    .type = choice
  atom_selection = None
    .type = str
  atom_selection_buffer = 3
    .type = float
  apply_sigma_scaling = True
    .type = bool
    .expert_level = 1
  apply_volume_scaling = False
    .type = bool
    .expert_level = 1
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
