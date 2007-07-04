import iotbx.phil
from mmtbx.refinement import rigid_body
from mmtbx import utils
from cctbx.array_family import flex
from iotbx.option_parser import iotbx_option_parser
from libtbx.str_utils import show_string
import libtbx.phil.command_line
import os, sys
from iotbx import pdb
from cctbx import crystal
from libtbx.utils import Sorry
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import crystal_symmetry_from_any
from mmtbx import monomer_library

modify_params_str = """\
selection = None
  .type = str
  .multiple = True
  .help = Selection for atoms to be modified
adp
  .help = Scope of options to modify ADP of selected atoms
{
  randomize = None
    .type = bool
    .help = Randomize ADP within a certain range
  set_b_iso = None
    .type = float
    .help = Set ADP of atoms to set_b_iso
  convert_to_isotropic = None
    .type = bool
    .help = Convert atoms to isotropic
  convert_to_anisotropic = None
    .type = bool
    .help = Convert atoms to anisotropic
  shift_b_iso = None
    .type = float
    .help = Add shift_b_iso value to ADP
  scale_adp = None
    .type = float
    .help = Multiply ADP by scale_adp
}
sites
  .help = Scope of options to modify coordinates of selected atoms
{
  shake = None
    .type = float
    .help = Randomize coordinates with mean error value equal to shake
  translate = 0 0 0
    .type = strings
    # XXX FUTURE float(3)
    .help = Translational shift
  rotate = 0 0 0
    .type = strings
    # XXX FUTURE float(3)
    .help = Rotational shift
  euler_angle_convention = *xyz zyz
    .type = choice
    .expert_level=2
    .help = Euler angles convention to be used for rotation
}
occupancies
  .help = Scope of options to modify occupancies of selected atoms
{
  randomize = None
    .type = bool
    .help = Randomize occupancies within a certain range
}
output
  .help = Write out PDB file with modified model (file name is defined in \
          write_modified)
{
  pdb {
      file_name=None
        .type=str
  }
}
"""
modify_params = iotbx.phil.parse(modify_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
%s
remove {
  selection = None
    .type = str
    .help = Select atoms to be removed
}
input {
  pdb
  {
    include scope mmtbx.utils.pdb_params
  }
  crystal_symmetry
    .help = Unit cell and space group parameters
  {
    unit_cell=None
      .type=unit_cell
    space_group=None
      .type=space_group
  }
}
"""%modify_params_str, process_includes=True)

class modify(object):
  def __init__(self, xray_structure, params, all_chain_proxies, log = None):
    self.log = log
    if(self.log is None): self.log = sys.stdout
    self.xray_structure = xray_structure
    self.all_chain_proxies = all_chain_proxies
    self.params = params
    self._adp_modified = False
    self._sites_modified = False
    self._occupancies_modified = False
    self.remove_selection = None
    try:
      params_remove_selection = self.params.remove.selection
    except: params_remove_selection = None
    if(params_remove_selection is not None):
      self.remove_selection = ~utils.get_atom_selections(
                         iselection        = False,
                         all_chain_proxies = all_chain_proxies,
                         selection_strings = [params_remove_selection],
                         xray_structure    = xray_structure)[0]
    self.selection = utils.get_atom_selections(
                                     iselection        = False,
                                     all_chain_proxies = all_chain_proxies,
                                     selection_strings = self.params.selection,
                                     xray_structure    = xray_structure)[0]
    self._show_selected()
    self._convert_to_isotropic()
    self._convert_to_anisotropic()
    self._set_b_iso()
    self._scale_adp()
    self._shift_b_iso()
    self._randomize_adp()
    self._shake_sites()
    self._rb_shift()
    self._randomize_occupancies()

  def _show_selected(self):
    assert self.xray_structure.scatterers().size() == self.selection.size()
    if(self.selection.size() == self.selection.count(True)):
      print >> self.log, \
                     "All input atoms are selected for requested manipulation."
    elif(self.selection.count(True) == 0):
      raise Sorry("No atoms selected for requested manipulation. \n"\
                  "       Selection string provided: %s"%self.params.selection)
    else:
      fmt = "Number of atoms selected for requested manipulation: %d out of total: %d"
      print >> self.log, fmt%(self.selection.count(True),self.selection.size())

  def _print_action(self, text):
    print >> self.log, text

  def _assert_not_modified(self, sites=False, adp=False, occupancies=False):
    assert [sites, adp, occupancies].count(True) in [0,1]
    #if(sites):
    #  if(self._sites_modified):
    #    raise Sorry("Can't modify cooridinates (already modified).")
    #  else: self._sites_modified = True
    if(adp):
      if(self._adp_modified):
        raise Sorry("Can't modify ADP (already modified).")
      else: self._adp_modified = True
    if(occupancies):
      if(self._occupancies_modified):
        raise Sorry("Can't modify occupancies (already modified).")
      else: self._occupancies_modified = True

  def _convert_to_isotropic(self):
    if(self.params.adp.convert_to_isotropic):
       self._print_action(text = "Converting to isotropic ADP.")
       self._assert_not_modified(adp = True)
       self.xray_structure.convert_to_isotropic(
                                       selection = self.selection.iselection())

  def _convert_to_anisotropic(self):
    if(self.params.adp.convert_to_anisotropic):
       self._print_action(text = "Converting to anisotropic ADP.")
       self._assert_not_modified(adp = True)
       self.xray_structure.convert_to_anisotropic(
                                       selection = self.selection.iselection())

  def _set_b_iso(self):
    b_iso = self.params.adp.set_b_iso
    if(b_iso is not None):
       self._print_action(text = "Setting all isotropic ADP to: %8.3f"%b_iso)
       self._assert_not_modified(adp = True)
       self.xray_structure.set_b_iso(value = b_iso, selection = self.selection)

  def _scale_adp(self):
    scale = self.params.adp.scale_adp
    if(scale is not None):
       self._print_action(text = "Scaling all ADP with: %.3f"%scale)
       self._assert_not_modified(adp = True)
       self.xray_structure.scale_adp(factor = scale, selection= self.selection)

  def _shift_b_iso(self):
    shift = self.params.adp.shift_b_iso
    if(shift is not None):
       self._print_action(text = "Shift all ADP by: %8.3f"%shift)
       self._assert_not_modified(adp = True)
       self.xray_structure.shift_us(b_shift   = shift,
                                    selection = self.selection.iselection())

  def _randomize_adp(self):
    if(self.params.adp.randomize):
       self._print_action(text = "Randomizing ADP.")
       self._assert_not_modified(adp = True)
       self.xray_structure.shake_adp(selection = self.selection)

  def _shake_sites(self):
    shake_value = self.params.sites.shake
    if(shake_value is not None):
       self._print_action(
                  text = "Shaking sites (RMS difference = %8.3f)."%shake_value)
       self._assert_not_modified(sites = True)
       self.xray_structure.shake_sites_in_place(rms_difference= shake_value,
                                                selection     = self.selection)

  def _rb_shift(self):
    rb = self.params.sites
    trans = [float(i) for i in rb.translate]
    rot   = [float(i) for i in rb.rotate]
    rb_shift = (abs(trans[0]+trans[1]+trans[2]+rot[0]+rot[1]+rot[2]) > 1.e-6)
    if(rb_shift):
       self._print_action(text = "Rigid body shift.")
       self._assert_not_modified(sites = True)
       if(self.params.sites.euler_angle_convention == "zyz"):
          rot_obj = rigid_body.rb_mat_euler(phi = rot[0],
                                            psi = rot[1],
                                            the = rot[2])
       else:
          rot_obj = rigid_body.rb_mat(phi = rot[0],
                                      psi = rot[1],
                                      the = rot[2])
       self.xray_structure.apply_rigid_body_shift(
                                    rot       = rot_obj.rot_mat().as_mat3(),
                                    trans     = trans,
                                    selection = self.selection.iselection())

  def _randomize_occupancies(self):
    if(self.params.occupancies.randomize):
       self._print_action(text = "Randomizing all occupancies.")
       self._assert_not_modified(occupancies = True)
       self.xray_structure.shake_occupancies(selection = self.selection)

def pdbtools(command_name, args):
  log = utils.set_log(args)
  utils.print_programs_start_header(
  log  = log,
  text= \
   "  phenix.pdbtools (or mmtbx.pdbtools): tools for PDB model manipulations.")
  command_line_interpreter = interpreter(command_name  = "phenix.pdbtools",
                                         args          = args,
                                         log           = log)
  utils.print_header("Complete set of parameters", out = log)
  master_params.format(command_line_interpreter.params).show(out = log)
  utils.print_header("xray structure summary", out = log)
  xray_structure = command_line_interpreter.processed_pdb_file.xray_structure(
                                    show_summary = True)
  if(xray_structure is None):
    raise Sorry("Cannot extract xray_structure.")
  all_chain_proxies = \
                  command_line_interpreter.processed_pdb_file.all_chain_proxies
  utils.print_header("Performing requested model manipulations", out = log)
  result = modify(xray_structure    = xray_structure,
                  params            = command_line_interpreter.params,
                  all_chain_proxies = all_chain_proxies,
                  log               = log)
  utils.print_header("Writing output model", out = log)
  atom_attributes_list = all_chain_proxies.stage_1.atom_attributes_list
  ofn = command_line_interpreter.params.output.pdb.file_name
  ifn = command_line_interpreter.params.input.pdb.file_name
  if(ofn is None):
    if(len(ifn)==1): ofn = os.path.basename(ifn[0]) + "_modified.pdb"
    else: ofn = os.path.basename(ifn[0]) + "_et_al_modified.pdb"
  print >> log, "Output file name: ", ofn
  ofo = open(ofn, "w")
  utils.write_pdb_file(xray_structure       = xray_structure,
                       atom_attributes_list = atom_attributes_list,
                       selection            = result.remove_selection,
                       out                  = ofo)
  utils.print_header("Done", out = log)

class interpreter:
  def __init__(self,
        command_name,
        args,
        log):
    self.command_name = command_name
    self.args = args
    self.log = log
    self.mon_lib_srv = monomer_library.server.server()
    self.ener_lib = monomer_library.server.ener_lib()
    self.params = None
    self.pdb_file_names = []
    self.processed_pdb_file = None
    self.processed_pdb_file_reference = None
    self.process_args()
    self.pdb_file_names.extend(self.params.input.pdb.file_name)
    self.processed_pdb_file, self.pdb_inp = utils.process_pdb_file(
               pdb_file_names            = self.pdb_file_names,
               parameters                = self.params.input.pdb,
               pdb_interpretation_params = None,
               mon_lib_srv               = self.mon_lib_srv,
               ener_lib                  = self.ener_lib,
               crystal_symmetry          = self.crystal_symmetry,
               stop_for_unknowns         = False,
               log                       = self.log)

  def process_args(self):
    args = self.args
    if (len(args) == 0): args = ["--help"]
    description_see_also \
        = 'See also: http://www.phenix-online.org/\n' +\
          'Questions / problems: phenixbb@phenix-online.org'
    self.command_line = (iotbx_option_parser(
      usage="phenix.pdbtools [options] [pdb_file] [parameter_file]",
      description='Example: phenix.pdbtools model.pdb parameters.txt\n\n'+\
         description_see_also)
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
        expert_level = self.command_line.expert_level)
      sys.exit(0)
    if (len(args) > 0):
      utils.print_header("Getting inputs", out = self.log)
      print >> self.log, "Command line arguments:", " ".join([show_string(arg)
        for arg in args])
    crystal_symmetries_from_coordinate_file = []
    parameter_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params = master_params,
      home_scope    = "pdbtools")
    parsed_params = []
    command_line_params = []
    utils.print_header("Processing inputs", out = self.log)
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
    self.crystal_symmetry = crystal.select_crystal_symmetry(
      from_command_line = self.command_line.symmetry,
      from_parameter_file=crystal.symmetry(
        unit_cell = self.params.input.crystal_symmetry.unit_cell,
        space_group_info = self.params.input.crystal_symmetry.space_group),
      from_coordinate_files=crystal_symmetries_from_coordinate_file)
    if (   self.crystal_symmetry.unit_cell() is None
        or self.crystal_symmetry.space_group_info() is None):
      raise Sorry(
        "Crystal symmetry is not defined. Please use the --symmetry option"
        " or define the refinement.crystal_symmetry parameters.")
    print >> self.log, "Working crystal symmetry after inspecting all inputs:"
    self.crystal_symmetry.show_summary(f = self.log, prefix="  ")
    self.params.input.crystal_symmetry.unit_cell = \
      self.crystal_symmetry.unit_cell()
    self.params.input.crystal_symmetry.space_group = \
      self.crystal_symmetry.space_group_info()
    self.point_group = self.crystal_symmetry.space_group() \
      .build_derived_point_group()
    print >> self.log


if (__name__ == "__main__" ):
  pdbtools(command_name = sys.argv[0], args = sys.argv[1:])
