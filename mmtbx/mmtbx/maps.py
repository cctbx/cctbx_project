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

master_params = iotbx.phil.parse("""\
maps {
  pdb_interpretation
    .expert_level = 2
  {
    include scope mmtbx.monomer_library.pdb_interpretation.master_params
  }
  use_experimental_phases = None
    .type = bool
  neutron_data = False
    .type = bool
  input {
    pdb {
      file_name=None
        .optional=True
        .type=path
        .multiple=True
    }
    data {
      file_name=None
        .type=path
      labels=None
        .type=strings
    }
    r_free_flags {
      file_name=None
        .type=path
      label=None
        .type=str
      test_flag_value=None
        .type=int
      disable_suitability_test=False
        .type=bool
    }
    experimental_phases {
      file_name=None
        .type=path
      labels=None
        .type=strings
    }
    monomers {
      file_name=None
        .optional=True
        .type=path
        .multiple=True
    }
  }
}
""", process_includes=True)

def print_header(log):
  print >> log
  host_and_user().show(out= log)
  print >> log, date_and_time()
  print >> log
  print >> log, "-"*79
  print >> log, \
   "  phenix.maps (or mmtbx.maps): tools for electron density maps calculation"
  print >> log, "-"*79
  print >> log, "Usage:"
  print >> log, "   phenix.maps [options] [reflection_file] [pdb_file] [parameter_file]"
  print >> log, "Exapmle:"
  print >> log, "   phenix.maps model.pdb data.mtz parameter_file.txt "
  print >> log, "Options:"
  print >> log, "   -h, --h, -help, --help   Show this help message and exit"
  print >> log, "   --show-defaults          Show all default parameters"
  print >> log, "   --quiet                  Suppress output to screen"
  print >> log, "See also:"
  print >> log, "   http://www.phenix-online.org/download/cci_apps/  (click on phenix.maps link)"
  print >> log, "Questions / problems: phenixbb@phenix-online.org"
  print >> log, "-"*79
  print >> log

def set_log(args):
  log = multi_out()
  if(not "--quiet" in args):
     log.register(label="stdout", file_object = sys.stdout)
  string_buffer = StringIO()
  log.register(label="log_buffer", file_object = string_buffer)
  return log

def maps(command_name, args):
  log = set_log(args = args)
  if(len(args) == 0 or "--help" in args or "--h" in args or "-h" in args or
     "-help" in args):
     print_header(log = log)
     sys.exit(0)
  elif(len(args) == 1 and "show"in args[0] and "defaults" in args[0] and
     ("-" in args[0] or "_" in args[0])):
     master_params.show(out = log)
     sys.exit(0)
  else:
    command_line_interpreter = interpreter(
      command_name  = "phenix.maps",
      master_params = master_params,
      args          = args,
      log           = log)
    raise Sorry("Not implemented yet.")

class interpreter:

  def __init__(self,
        command_name,
        master_params,
        args,
        log=None):
    self.command_name = command_name
    self.master_params = master_params
    self.args = args
    if (log is None): log = sys.stdout
    self.log = log
    self.mon_lib_srv = monomer_library.server.server()
    self.ener_lib = monomer_library.server.ener_lib()
    self.params = None
    self.pdb_file_names = []
    self.processed_pdb_file = None
    self.reflection_files = []
    self.data = None
    self.r_free_flags = None
    self.experimental_phases = None
    self.process_args()
    mmtbx.refinement.print_statistics.make_header(
                                        "Processing input PDB file", out = log)
    self.process_pdb_file()
    mmtbx.refinement.print_statistics.make_header(
                                     "Processing experimental data", out = log)
    reflection_file_server = self.reflection_file_server()
    self.determine_data(reflection_file_server)
    self.determine_r_free_flags(reflection_file_server)
    self.determine_experimental_phases(reflection_file_server)
    mmtbx.refinement.print_statistics.make_header(
                             "Complete set of effective parameters", out = log)
    print >> log, "#phil __ON__"
    self.master_params.format(self.params).show(out = log)
    print >> log, "#phil __OFF__"

  def process_args(self):
    args = self.args
    print >> self.log
    if (len(args) > 0):
      print >> self.log, "Command line arguments:", " ".join([show_string(arg)
        for arg in args])
      print >> self.log
    crystal_symmetries_from_coordinate_file = []
    crystal_symmetries_from_reflection_file = []
    cif_objects = []
    parameter_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params = self.master_params,
      home_scope = "maps")
    parsed_params = []
    command_line_params = []
    print >> self.log, "Processing inputs. This may take a minute or two."
    for arg in self.args:
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
        elif (pdb.interpretation.is_pdb_file(file_name=arg)):
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
    self.params, unused_definitions = self.master_params.fetch(
      sources=parsed_params+command_line_params,
      track_unused_definitions=True)
    if(len(unused_definitions)):
      print >> self.log, "*"*79
      print >> self.log, "Unused parameter definitions:"
      for obj_loc in unused_definitions:
        print >> self.log, " ", str(obj_loc)
      print >> self.log
      raise Sorry("""
  Please check the input file(s) for spelling errors and obsolete
  parameter definitions.""")
    self.params = self.params.extract()
    self.process_monomer_cif_files(cif_objects=cif_objects)
    del cif_objects
    self.crystal_symmetry = crystal.select_crystal_symmetry(
      from_coordinate_files = crystal_symmetries_from_coordinate_file,
      from_reflection_files = crystal_symmetries_from_reflection_file)
    if (   self.crystal_symmetry.unit_cell() is None
        or self.crystal_symmetry.space_group_info() is None):
      raise Sorry(
        "Crystal symmetry is not defined. Please make sure the model or data\n"\
        "       file has crystall symmetry information.")
    print >> self.log, "Working crystal symmetry after inspecting all inputs:"
    self.crystal_symmetry.show_summary(f=self.log, prefix="  ")
    self.point_group = self.crystal_symmetry.space_group() \
      .build_derived_point_group()
    print >> self.log

  def process_monomer_cif_files(self, cif_objects):
    all = []
    index_dict = {}
    for file_name in self.params.maps.input.monomers.file_name:
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
    del self.params.maps.input.monomers.file_name[:]
    for file_name,cif_object in unique:
      if (cif_object is None):
        self.mon_lib_srv.process_cif(file_name=file_name)
        self.ener_lib.process_cif(file_name=file_name)
      else:
        self.mon_lib_srv.process_cif_object(
          cif_object=cif_object, file_name=file_name)
        self.ener_lib.process_cif_object(
          cif_object=cif_object, file_name=file_name)
      self.params.refinement.input.monomers.file_name.append(file_name)

  def process_pdb_file(self):
    self.pdb_file_names.extend(self.params.maps.input.pdb.file_name)
    if (len(self.pdb_file_names) == 0):
      raise Sorry("No coordinate file given.")
    self.params.maps.input.pdb.file_name = [
      os.path.abspath(file_name) for file_name in self.pdb_file_names]
    if (len(self.pdb_file_names) == 1):
      pdb_file_name = self.pdb_file_names[0]
      raw_records = None
      self.pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    else:
      pdb_file_name = None
      raw_records = []
      raw_records_flex = flex.std_string()
      for file_name in self.pdb_file_names:
        raw_records.extend(open(file_name).readlines())
        raw_records_flex.extend(flex.split_lines(open(file_name).read()))
      self.pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records_flex)
    self.processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=self.mon_lib_srv,
      ener_lib=self.ener_lib,
      params = self.params.maps.pdb_interpretation,
      file_name = pdb_file_name,
      raw_records = raw_records,
      strict_conflict_handling=False,
      crystal_symmetry=self.crystal_symmetry,
      force_symmetry=True,
      log=self.log)
    print >> self.log
    msg = self.processed_pdb_file.all_chain_proxies.fatal_problems_message()
    if (msg is not None):
      msg = "\n  ".join([msg,
        "Please edit the PDB file to resolve the problems and/or supply a",
        "CIF file with matching restraint definitions, along with",
        "apply_cif_modification and apply_cif_link parameter definitions",
        "if necessary (see phenix.refine documentation).",
        "Also note that elbow.builder is available to create restraint",
        "definitions for unknown ligands."])
      raise Sorry(msg)

  def reflection_file_server(self):
    return reflection_file_utils.reflection_file_server(
      crystal_symmetry=self.crystal_symmetry,
      force_symmetry=True,
      reflection_files=self.reflection_files,
      err=self.log)

  def miller_array_symmetry_safety_check(self, miller_array, data_description):
    msg = miller_array.crystal_symmetry_is_compatible_with_symmetry_from_file(
      working_point_group=self.point_group).format_error_message(
        data_description=data_description)
    if (msg is not None):
      raise Sorry(msg + """
  %s inspects all inputs to determine the working crystal
  symmetry (unit cell & space group).
  Please check the working crystal symmetry shown above.
""" % self.command_name)

  def determine_data(self, reflection_file_server):
    p = self.params.maps.input.data
    self.data = reflection_file_server.get_xray_data(
      file_name=p.file_name,
      labels=p.labels,
      ignore_all_zeros=True,
      parameter_scope="refinement.input.data")
    p.file_name = self.data.info().source
    p.labels = [self.data.info().label_string()]
    if (self.data.is_xray_intensity_array()):
      print >> self.log, "I-obs:"
    else:
      print >> self.log, "F-obs:"
    print >> self.log, " ", self.data.info()
    self.miller_array_symmetry_safety_check(
      miller_array=self.data,
      data_description="X-ray data")
    print >> self.log
    info = self.data.info()
    processed = self.data.eliminate_sys_absent(log=self.log)
    if (processed is not self.data):
      info = info.customized_copy(systematic_absences_eliminated=True)
    if (not processed.is_unique_set_under_symmetry()):
      if (self.data.is_xray_intensity_array()):
        print >> self.log, "Merging symmetry-equivalent intensities:"
      else:
        print >> self.log, "Merging symmetry-equivalent amplitudes:"
      merged = processed.merge_equivalents()
      merged.show_summary(out=self.log, prefix="  ")
      print >> self.log
      processed = merged.array()
      info = info.customized_copy(merged=True)
    self.data = processed.set_info(info)

  def determine_r_free_flags(self, reflection_file_server):
    log = self.log
    r_free_flags = None
    params = self.params.maps
    p = params.input.r_free_flags
    try:
      parameter_scope="maps.input.r_free_flags"
      r_free_flags, test_flag_value = \
        reflection_file_server.get_r_free_flags(
          file_name=p.file_name,
          label=p.label,
          test_flag_value=p.test_flag_value,
          disable_suitability_test=p.disable_suitability_test,
          parameter_scope=parameter_scope)
    except reflection_file_utils.Sorry_No_array_of_the_required_type, e:
      e.reset_tracebacklimit()
      r_free_flags, test_flag_value = None, None
    else:
      p.file_name = r_free_flags.info().source
      p.label = r_free_flags.info().label_string()
      p.test_flag_value = test_flag_value
      print >> log, "R-free flags:"
      print >> log, " ", r_free_flags.info()
      self.miller_array_symmetry_safety_check(
        miller_array=r_free_flags,
        data_description="R-free flags")
      print >> log
      info = r_free_flags.info()
      processed = r_free_flags.eliminate_sys_absent(log=log)
      if (processed is not r_free_flags):
        info = info.customized_copy(systematic_absences_eliminated=True)
      if (not processed.is_unique_set_under_symmetry()):
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
    self.r_free_flags = r_free_flags

  def determine_experimental_phases(self, reflection_file_server):
    log = self.log
    if (self.params.maps.use_experimental_phases):
      p = self.params.maps.input.experimental_phases
      try:
        self.experimental_phases = \
          reflection_file_server.get_experimental_phases(
            file_name=p.file_name,
            labels=p.labels,
            ignore_all_zeros=True,
            parameter_scope="refinement.input.experimental_phases")
      except reflection_file_utils.Sorry_No_array_of_the_required_type:
        if (fl is not None): raise
        self.experimental_phases = None
      else:
        p.file_name = self.experimental_phases.info().source
        p.labels = [self.experimental_phases.info().label_string()]
        print >> log, "Experimental phases:"
        print >> log, " ", self.experimental_phases.info()
        self.miller_array_symmetry_safety_check(
          miller_array=self.experimental_phases,
          data_description="Experimental phases")
        print >> log
        info = self.experimental_phases.info()
        processed = self.experimental_phases.eliminate_sys_absent(log=log)
        if (processed is not self.experimental_phases):
          info = info.customized_copy(systematic_absences_eliminated=True)
        if (not processed.is_unique_set_under_symmetry()):
          print >> log, \
            "Merging symmetry-equivalent Hendrickson-Lattman coefficients:"
          merged = processed.merge_equivalents()
          merged.show_summary(out=log, prefix="  ")
          print >> log
          processed = merged.array()
          info = info.customized_copy(merged=True)
        self.experimental_phases = processed.set_info(info)

  def inputs(self):
    return inputs(
      params=self.params,
      processed_pdb_file=self.processed_pdb_file,
      pdb_inp=self.pdb_inp,
      data=self.data,
      r_free_flags=self.r_free_flags,
      experimental_phases=self.experimental_phases)

class inputs:
  def __init__(self,
        params=None,
        processed_pdb_file=None,
        pdb_inp=None,
        data=None,
        r_free_flags=None,
        experimental_phases=None):
    adopt_init_args(self, locals())


if (__name__ == "__main__" ):
  maps(command_name = sys.argv[0], args = sys.argv[1:])
