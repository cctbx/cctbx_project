from __future__ import generators
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
import mmtbx.restraints
import mmtbx.model
from mmtbx import model_statistics
import random
from libtbx import easy_run


fmodel_from_xray_structure_params_str = """\
k_sol = 0.0
  .type = float
  .help = Bulk solvent k_sol values
b_sol = 0.0
  .type = float
  .help = Bulk solvent b_sol values
b_cart = 0 0 0 0 0 0
  .type = strings
  # XXX FUTURE float(6)
  .help = Anisotropic scale matrix
scale = 1.0
  .type = float
  .help = Overall scale factor
structure_factors_accuracy {
  include scope mmtbx.f_model.sf_and_grads_accuracy_params
}
mask {
  include scope mmtbx.masks.mask_master_params
}
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
  .help = Choices of scattering table for structure factors calculations
"""
fmodel_from_xray_structure_params = iotbx.phil.parse(
  fmodel_from_xray_structure_params_str, process_includes=True)

fmodel_from_xray_structure_master_params_str = """\
high_resolution = None
  .type = float
low_resolution = None
  .type = float
%s
output {
  format = *mtz cns
    .type = choice
  label = FMODEL
    .type = str
  type = real *complex
    .type = choice
  file_name = None
    .type = str
}
"""%fmodel_from_xray_structure_params_str

fmodel_from_xray_structure_master_params = iotbx.phil.parse(
  fmodel_from_xray_structure_master_params_str, process_includes=True)

modify_params_str = """\
selection = None
  .type = str
  .help = Selection for atoms to be modified
adp
  .help = Scope of options to modify ADP of selected atoms
  .multiple = True
{
  selection = None
    .type = str
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
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
  .multiple = True
{
  selection = None
    .type = str
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
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
  set = None
    .type = float
    .help = Set all or selected occupancies to given value
}
output
  .help = Write out PDB file with modified model (file name is defined in \
          write_modified)
{
  file_name=None
    .type=str
}
random_seed = None
  .type = int
  .help = Random seed
"""
modify_params = iotbx.phil.parse(modify_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
%s
remove = None
  .type = str
  .help = Selection for the atoms to be removed
keep = None
  .type = str
  .help = Select atoms to keep
put_into_box_with_buffer = None
  .type = float
  .help = Move molecule into center of box.
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
f_model {
%s
}
pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
}
geometry_minimization {
  include scope mmtbx.command_line.geometry_minimization.master_params
}
"""%(modify_params_str, fmodel_from_xray_structure_master_params_str),
     process_includes=True)

class modify(object):
  def __init__(self, xray_structure, params, all_chain_proxies, log = None):
    self.log = log
    if(self.log is None): self.log = sys.stdout
    self.xray_structure = xray_structure
    self.all_chain_proxies = all_chain_proxies
    self.params = params
    self._occupancies_modified = False
    self.remove_selection = None
    self.keep_selection = None
    if(self.params.random_seed is not None):
      random.seed(self.params.random_seed)
      flex.set_random_seed(self.params.random_seed)
    try: params_remove_selection = self.params.remove
    except KeyboardInterrupt: raise
    except: params_remove_selection = None
    try: params_keep_selection = self.params.keep
    except KeyboardInterrupt: raise
    except: params_keep_selection = None
    if(params_remove_selection is not None):
      self.remove_selection = flex.smart_selection(
        flags=~utils.get_atom_selections(
          iselection        = False,
          all_chain_proxies = all_chain_proxies,
          selection_strings = [params_remove_selection],
          xray_structure    = xray_structure,
          log               = self.log)[0])
    if(params_keep_selection is not None):
      self.keep_selection = flex.smart_selection(
        flags=utils.get_atom_selections(
          iselection        = False,
          all_chain_proxies = all_chain_proxies,
          selection_strings = [params_keep_selection],
          xray_structure    = xray_structure,
          log               = self.log)[0])
    if(self.remove_selection is not None and self.keep_selection is not None):
      raise Sorry("Ambiguous selection: 'keep' and 'remove' keywords cannot"
        " be used simultaneously:\n  keep=%s\n  remove=%s" % (
          show_string(params_keep_selection),
          show_string(params_remove_selection)))
    if(self.keep_selection is not None):
      assert self.remove_selection is None
      self.remove_selection = flex.smart_selection(flags =
        self.keep_selection.flags)
    self.top_selection = flex.smart_selection(
      flags=utils.get_atom_selections(
        iselection        = False,
        all_chain_proxies = all_chain_proxies,
        selection_strings = [self.params.selection],
        xray_structure    = xray_structure,
        log               = self.log)[0])
    self._process_adp()
    self._process_sites()
    self._process_occupancies(selection = self.top_selection)
    try: self._put_in_box()
    except KeyboardInterrupt: raise
    except: pass

  def _put_in_box(self):
    if(self.params.put_into_box_with_buffer is not None):
      result = \
        self.xray_structure.orthorhombic_unit_cell_around_centered_scatterers(
          buffer_size = self.params.put_into_box_with_buffer)
      self.xray_structure.replace_scatterers(result.scatterers())

  def _print_action(self, text, selection):
    print >> self.log, "%s: selected atoms: %s" % (
      text, selection.format_summary())

  def _process_adp(self):
    for adp in self.params.adp:
      if (adp.selection is None):
        selection = self.top_selection
      else:
        selection = flex.smart_selection(
          flags=utils.get_atom_selections(
            iselection=False,
            all_chain_proxies=self.all_chain_proxies,
            selection_strings=[adp.selection],
            xray_structure=self.xray_structure,
            log = self.log)[0])
      if (adp.convert_to_isotropic):
        self._convert_to_isotropic(selection=selection)
      if (adp.convert_to_anisotropic):
        self._convert_to_anisotropic(selection=selection)
      self._set_b_iso(selection=selection, b_iso=adp.set_b_iso)
      self._scale_adp(selection=selection, factor=adp.scale_adp)
      self._shift_b_iso(selection=selection, shift=adp.shift_b_iso)
      if (adp.randomize):
        self._randomize_adp(selection=selection)

  def _convert_to_isotropic(self, selection):
    self._print_action(
      text = "Converting to isotropic ADP",
      selection = selection)
    self.xray_structure.convert_to_isotropic(selection=selection.indices)

  def _convert_to_anisotropic(self, selection):
    self._print_action(
      text = "Converting to anisotropic ADP",
      selection = selection)
    self.xray_structure.convert_to_anisotropic(selection=selection.flags)

  def _set_b_iso(self, selection, b_iso):
    if (b_iso is not None):
      self._print_action(
        text = "Setting all isotropic ADP = %.3f" % b_iso,
        selection = selection)
      self.xray_structure.set_b_iso(value=b_iso, selection=selection.flags)

  def _scale_adp(self, selection, factor):
    if (factor is not None):
      self._print_action(
        text = "Multiplying all ADP with factor = %.6g" % factor,
        selection = selection)
      self.xray_structure.scale_adp(factor=factor, selection=selection.flags)

  def _shift_b_iso(self, selection, shift):
    if (shift is not None):
      self._print_action(
        text = "Adding shift = %.2f to all ADP" % shift,
        selection = selection)
      self.xray_structure.shift_us(b_shift=shift, selection=selection.indices)

  def _randomize_adp(self, selection):
    self._print_action(
      text = "Randomizing ADP",
      selection = selection)
    self.xray_structure.shake_adp(selection=selection.flags)

  def _process_sites(self):
    for sites in self.params.sites:
      if (sites.selection is None):
        selection = self.top_selection
      else:
        selection = flex.smart_selection(
          flags=utils.get_atom_selections(
            iselection=False,
            all_chain_proxies=self.all_chain_proxies,
            selection_strings=[sites.selection],
            xray_structure=self.xray_structure,
            log = self.log)[0])
      self._shake_sites(selection=selection, rms_difference=sites.shake)
      self._rb_shift(
        selection=selection,
        translate=sites.translate,
        rotate=sites.rotate,
        euler_angle_convention=sites.euler_angle_convention)

  def _shake_sites(self, selection, rms_difference):
    if (rms_difference is not None):
      self._print_action(
        text = "Shaking sites (RMS = %.3f)" % rms_difference,
        selection = selection)
      self.xray_structure.shake_sites_in_place(
        rms_difference=rms_difference,
        selection=selection.flags)

  def _rb_shift(self, selection, translate, rotate, euler_angle_convention):
    trans = [float(i) for i in translate]
    rot   = [float(i) for i in rotate]
    if(len(trans) != 3): raise Sorry("Wrong value: translate= " + translate)
    if(len(rot) != 3): raise Sorry("Wrong value: translate= " + rotate)
    if (   trans[0] != 0 or trans[1] != 0 or trans[2] != 0
        or rot[0] != 0 or rot[1] != 0 or rot[2] != 0):
      self._print_action(
        text = "Rigid body shift",
        selection = selection)
      if (euler_angle_convention == "zyz"):
        rot_obj = rigid_body.rb_mat_zyz(phi = rot[0],
                                          psi = rot[1],
                                          the = rot[2])
      else:
        rot_obj = rigid_body.rb_mat_xyz(phi = rot[0],
                                    psi = rot[1],
                                    the = rot[2])
      self.xray_structure.apply_rigid_body_shift(
        rot       = rot_obj.rot_mat().as_mat3(),
        trans     = trans,
        selection = selection.indices)

  def _process_occupancies(self, selection):
    if(self.params.occupancies.randomize):
      self._print_action(
        text = "Randomizing occupancies",
        selection = selection)
      if (self._occupancies_modified):
        raise Sorry("Can't modify occupancies (already modified).")
      else:
        self._occupancies_modified = True
      self.xray_structure.shake_occupancies(selection=selection.flags)
    if(self.params.occupancies.set is not None):
      self._print_action(
        text = "Setting occupancies to: %8.3f"%self.params.occupancies.set,
        selection = selection)
      if (self._occupancies_modified):
        raise Sorry("Can't modify occupancies (already modified).")
      else:
        self._occupancies_modified = True
      self.xray_structure.set_occupancies(value = self.params.occupancies.set,
                                          selection = selection.flags)

  def report_number_of_atoms_to_be_removed(self):
    if (    self.remove_selection is not None
        and self.remove_selection.selected_size > 0):
      self.remove_selection.show_summary(
        out = self.log,
        label = "Atoms to be kept: ")

class fmodel_from_xray_structure(object):
  def __init__(self, xray_structure, f_obs = None, params = None,
                     target_name = "ml", r_free_flags_fraction = 0.1):
    if(params is None):
      params = fmodel_from_xray_structure_master_params.extract()
    if(f_obs is None):
      hr = None
      try: hr = params.high_resolution
      except: self.Sorry_high_resolution_is_not_defined()
      if(params.scattering_table == "neutron"):
        xray_structure.switch_to_neutron_scattering_dictionary()
      else:
        xray_structure.scattering_type_registry(
          table = params.scattering_table, d_min = hr)
      if(hr is None): self.Sorry_high_resolution_is_not_defined()
      f_obs = xray_structure.structure_factors(d_min = hr).f_calc()
      sfga = params.structure_factors_accuracy
      f_obs = f_obs.structure_factors_from_scatterers(
         xray_structure = xray_structure,
         algorithm                    = sfga.algorithm,
         cos_sin_table                = sfga.cos_sin_table,
         grid_resolution_factor       = sfga.grid_resolution_factor,
         quality_factor               = sfga.quality_factor,
         u_base                       = sfga.u_base,
         b_base                       = sfga.b_base,
         wing_cutoff                  = sfga.wing_cutoff,
         exp_table_one_over_step_size = sfga.exp_table_one_over_step_size
         ).f_calc()
      lr = None
      try: lr = params.low_resolution
      except: RuntimeError("Parameter scope does not have 'low_resolution'.")
      if(params.low_resolution is not None):
        f_obs = f_obs.resolution_filter(d_max = lr)
    else:
      hr = None
      try: hr = params.high_resolution
      except AttributeError: pass
      except: raise RuntimeError
      lr = None
      try: lr = params.low_resolution
      except AttributeError: pass
      except: raise RuntimeError
      f_obs = f_obs.resolution_filter(d_max = lr, d_min = hr)
      if(params.scattering_table == "neutron"):
        xray_structure.switch_to_neutron_scattering_dictionary()
      else:
        xray_structure.scattering_type_registry(
          table = params.scattering_table, d_min = hr)
    r_free_flags = f_obs.generate_r_free_flags(fraction = r_free_flags_fraction)
    fmodel = mmtbx.f_model.manager(
      xray_structure               = xray_structure,
      sf_and_grads_accuracy_params = params.structure_factors_accuracy,
      r_free_flags                 = r_free_flags,
      mask_params                  = params.mask,
      target_name                  = target_name,
      f_obs                        = abs(f_obs),
      k_sol                        = params.k_sol,
      b_sol                        = params.b_sol,
      b_cart                       = [float(i) for i in params.b_cart])
    f_model = fmodel.f_model()
    f_model = f_model.array(data = f_model.data()*params.scale)
    try:
      if(params.output.type == "real"):
        f_model = abs(f_model)
        f_model.set_observation_type_xray_amplitude()
    except AttributeError: pass
    except: raise RuntimeError
    self.f_model = f_model
    self.params = params
    self.fmodel = fmodel

  def Sorry_high_resolution_is_not_defined(self):
    raise Sorry("High resolution limit is not defined. "\
      "Use 'high_resolution' keyword to define it.")

  def write_to_file(self, file_name):
    assert self.params.output.format in ["mtz", "cns"]
    assert file_name is not None
    op = self.params.output
    if(self.params.output.format == "cns"):
      ofo = open(file_name, "w")
      print >> ofo, "NREFlection=%d" % self.f_model.indices().size()
      print >> ofo, "ANOMalous=%s" % {0: "FALSE"}.get(
        int(self.f_model.anomalous_flag()), "TRUE")
      for n_t in [("%s"%op.label, "%s"%op.type.upper())]:
        print >> ofo, "DECLare NAME=%s DOMAin=RECIprocal TYPE=%s END"%n_t
      if(op.type == "complex"):
        arrays = [
          self.f_model.indices(), flex.abs(self.f_model.data()),
          self.f_model.phases(deg=True).data()]
        for values in zip(*arrays):
          print >> ofo, "INDE %d %d %d" % values[0],
          print >> ofo, " %s= %.6g %.6g" % (op.label, values[1],values[2])
      else:
        arrays = [
          self.f_model.indices(), self.f_model.data()]
        for values in zip(*arrays):
          print >> ofo, "INDE %d %d %d" % values[0],
          print >> ofo, " %s= %.6g" % (op.label, values[1])
    else:
      mtz_dataset= self.f_model.as_mtz_dataset(column_root_label="%s"%op.label)
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.write(file_name = file_name)

def run(args, command_name="phenix.pdbtools"):
  log = utils.set_log(args)
  utils.print_programs_start_header(
    log  = log,
    text = "  %s tools for PDB model manipulations." % command_name)
  command_line_interpreter = interpreter(command_name  = command_name,
                                         args          = args,
                                         log           = log)
  all_chain_proxies = \
    command_line_interpreter.processed_pdb_file.all_chain_proxies
  xray_structure = command_line_interpreter.processed_pdb_file.xray_structure(
    show_summary = False)
  if(xray_structure is None):
    raise Sorry("Cannot extract xray_structure.")
### get i/o file names
  ofn = command_line_interpreter.params.output.file_name
  ifn = command_line_interpreter.params.input.pdb.file_name
  if(ofn is None):
    if(len(ifn)==1): ofn = os.path.basename(ifn[0]) + "_modified.pdb"
    else: ofn = os.path.basename(ifn[0]) + "_et_al_modified.pdb"
### show_geometry_statistics and exit
  if(command_line_interpreter.command_line.options.show_geometry_statistics):
    utils.print_header("Geometry statistics", out = log)
    command_line_interpreter.processed_pdb_file.log = None # to disable output
    geometry = command_line_interpreter.processed_pdb_file.\
      geometry_restraints_manager(show_energies = True)
    restraints_manager = mmtbx.restraints.manager(
      geometry = geometry,
      normalization = True)
    model_statistics.geometry(
      sites_cart         = xray_structure.sites_cart(),
      restraints_manager = restraints_manager).show(out = log)
    return
### show_adp_statistics and exit
  if(command_line_interpreter.command_line.options.show_adp_statistics):
    utils.print_header("ADP statistics", out = log)
    model = mmtbx.model.manager(
      xray_structure = xray_structure,
      pdb_hierarchy = all_chain_proxies.pdb_hierarchy,
      log = log)
    model.show_adp_statistics(out = log, padded = True)
    return
### add hydrogens and exit
  if(command_line_interpreter.command_line.options.add_h):
    utils.print_header("Adding hydrogen atoms", out = log)
    if(len(ifn) > 1): raise Sorry("Multiple input PDB files found.")
    ifn = ifn[0]
    if(not os.path.isfile(ifn)): raise Sorry("File %s does not exist."%ifn)
    easy_run.go("phenix.reduce %s > %s"% (ifn, ofn))
    print >> log, "Output model file name (with H added): %s\n"%ofn
    return
### show parameters
  utils.print_header("Complete set of parameters", out = log)
  master_params.format(command_line_interpreter.params).show(out = log)
### run geometry regularization
  if(command_line_interpreter.command_line.options.geometry_regularization):
    utils.print_header("Geometry regularization", out = log)
    from mmtbx.command_line import geometry_minimization
    sites_cart = geometry_minimization.run(
      params = command_line_interpreter.params.geometry_minimization,
      processed_pdb_file = command_line_interpreter.processed_pdb_file,
      log = log)
    xray_structure = xray_structure.replace_sites_cart(new_sites = sites_cart)
  utils.print_header("Performing requested model manipulations", out = log)
  result = modify(xray_structure    = xray_structure,
                  params            = command_line_interpreter.params,
                  all_chain_proxies = all_chain_proxies,
                  log               = log)
  result.report_number_of_atoms_to_be_removed()
  utils.print_header("Writing output model", out = log)
### write output file (if got to this point)
  print >> log, "Output model file name: ", ofn
  ofo = open(ofn, "w")
  utils.write_pdb_file(
    xray_structure       = xray_structure,
    pdb_hierarchy        = all_chain_proxies.pdb_hierarchy,
    pdb_atoms            = all_chain_proxies.pdb_atoms,
    selection            = getattr(result.remove_selection, "flags", None),
    write_cryst1_record  = not command_line_interpreter.fake_crystal_symmetry,
    out                  = ofo)
  ofo.close()
  if(command_line_interpreter.command_line.options.f_model):
    par = command_line_interpreter.params.f_model
    if(par.output.format == "cns"): extension = ".hkl"
    elif(par.output.format == "mtz"): extension = ".mtz"
    else: extension = ".reflections"
    ofn = par.output.file_name
    if(ofn is None):
      if(len(ifn)==1): ofn = os.path.basename(ifn[0]) + extension
      else: ofn = os.path.basename(ifn[0]) + "_et_al" + extension
    utils.print_header("Compute model structure factors", out = log)
    print >> log, "Output reflections file name: ", ofn
    fmodel_from_xray_structure(
      xray_structure = xray_structure,
      params         = command_line_interpreter.params.f_model).write_to_file(
        file_name = ofn)
  utils.print_header("Done", out = log)

class interpreter:
  def __init__(self,
        command_name,
        args,
        log):
    self.command_name = command_name
    self.args = args
    self.log = log
    self.fake_crystal_symmetry = False
    self.mon_lib_srv = monomer_library.server.server()
    self.ener_lib = monomer_library.server.ener_lib()
    self.params = None
    self.pdb_file_names = []
    self.cif_objects = []
    self.processed_pdb_file = None
    self.processed_pdb_file_reference = None
    self.process_args()
    self.pdb_file_names.extend(self.params.input.pdb.file_name)
    processed_pdb_files_srv = utils.process_pdb_file_srv(
      crystal_symmetry          = self.crystal_symmetry,
      pdb_parameters            = self.params.input.pdb,
      pdb_interpretation_params = self.params.pdb_interpretation,
      stop_for_unknowns         = True,
      log                       = self.log,
      mon_lib_srv               = self.mon_lib_srv,
      ener_lib                  = self.ener_lib,
      cif_objects               = self.cif_objects)
    self.processed_pdb_file, self.pdb_inp = \
        processed_pdb_files_srv.process_pdb_files(pdb_file_names =
          self.pdb_file_names)

  def process_args(self):
    args = self.args
    if (len(args) == 0): args = ["--help"]
    description_see_also \
        = 'See also: http://www.phenix-online.org/\n' +\
          'Questions / problems: phenixbb@phenix-online.org'
    self.command_line = (iotbx_option_parser(
      usage="%s [options] [pdb_file] [parameter_file]" % self.command_name,
      description='Example: %s model.pdb parameters.txt\n\n'
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
      .option("--show_adp_statistics",
          action="store_true",
          help="Show complete ADP (B-factors) statistics.")
      .option("--show_geometry_statistics",
          action="store_true",
          help="Show complete ADP (B-factors) statistics.")
      .option("--add_h",
          action="store_true",
          help="Add H atoms to a model using Reduce program.")
      .option("--geometry_regularization",
          action="store_true",
          help="Perform geometry regularization.")
      .option("--f_model",
          action="store_true",
          help="Compute total model structure factors (F_model) and output into a file.")
    ).process(args=args)
    if(self.command_line.expert_level is not None):
      master_params.show(
        out = self.log,
        expert_level = self.command_line.expert_level,
        attributes_level = self.command_line.attributes_level)
      sys.exit(0)
    if (len(args) > 0):
      utils.print_header("Getting inputs", out = self.log)
      print >> self.log, "Command line arguments:", " ".join([show_string(arg)
        for arg in args])
    crystal_symmetries_from_coordinate_file = []
    parameter_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil = master_params,
      home_scope  = "pdbtools")
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
          try: cif_object = mmtbx.monomer_library.server.read_cif(file_name=arg)
          except KeyboardInterrupt: raise
          except: pass
          else:
            if(len(cif_object) > 0):
              self.cif_objects.append((arg,cif_object))
              arg_is_processed = True
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
    # XXX also need to inspect input files defined via parameters
    self.crystal_symmetry = crystal.select_crystal_symmetry(
      from_command_line = self.command_line.symmetry,
      from_parameter_file=crystal.symmetry(
        unit_cell = self.params.input.crystal_symmetry.unit_cell,
        space_group_info = self.params.input.crystal_symmetry.space_group),
      from_coordinate_files=crystal_symmetries_from_coordinate_file)
    if (   self.crystal_symmetry.unit_cell() is None
        or self.crystal_symmetry.space_group_info() is None):
      self.fake_crystal_symmetry = True
      print >> self.log, "*** No input crystal symmetry found. ***"
      print >> self.log,"Functionality requiring crystal symmetry unavailable."
      print >> self.log
      self.crystal_symmetry = crystal.symmetry((1, 1, 1, 90, 90, 90), "P 1")
    print >> self.log, "Working crystal symmetry after inspecting all inputs:"
    self.crystal_symmetry.show_summary(f = self.log, prefix="  ")
    self.params.input.crystal_symmetry.unit_cell = \
      self.crystal_symmetry.unit_cell()
    self.params.input.crystal_symmetry.space_group = \
      self.crystal_symmetry.space_group_info()
    self.point_group = self.crystal_symmetry.space_group() \
      .build_derived_point_group()
    print >> self.log
