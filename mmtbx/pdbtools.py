import iotbx.phil
from mmtbx.refinement import rigid_body
from mmtbx import utils
from cctbx.array_family import flex
import iotbx.cns.miller_array
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
from libtbx import easy_run, easy_pickle
from iotbx.pdb import combine_unique_pdb_files
from mmtbx.command_line import lockit
from libtbx import runtime_utils

modify_params_str = """\
selection = None
  .type = str
  .help = Selection for atoms to be modified
  .short_caption = Modify atom selection
  .input_size=400
  .style = bold selection
adp
  .help = Scope of options to modify ADP of selected atoms
  .multiple = True
  .expert_level=2
  .short_caption=Modify ADPs
  .style = auto_align menu_item parent_submenu:model_modifications noauto
{
  atom_selection = None
    .type = str
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
    .short_caption = Modify ADPs for selection
    .input_size=400
    .style  = bold selection
  randomize = None
    .type = bool
    .help = Randomize ADP within a certain range
    .short_caption=Randomize ADPs
  set_b_iso = None
    .type = float
    .help = Set ADP of atoms to set_b_iso
    .short_caption=Set isotropic B to
  convert_to_isotropic = None
    .type = bool
    .help = Convert atoms to isotropic
  convert_to_anisotropic = None
    .type = bool
    .help = Convert atoms to anisotropic
  shift_b_iso = None
    .type = float
    .help = Add shift_b_iso value to ADP
    .short_caption=Increase B_iso by
  scale_adp = None
    .type = float
    .help = Multiply ADP by scale_adp
    .short_caption=ADP scale factor
}
sites
  .help = Scope of options to modify coordinates of selected atoms
  .multiple = True
  .short_caption=Modify coordinates
  .expert_level=2
  .style = auto_align noauto menu_item parent_submenu:model_modifications
{
  atom_selection = None
    .type = str
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
    .input_size=400
    .short_caption = Modify sites for selection
    .style = bold selection
  shake = None
    .type = float
    .help = Randomize coordinates with mean error value equal to shake
  max_rotomer_distortion = None
    .type = bool
    .help = Switch to a rotomer maximally distant from the current one
  min_rotomer_distortion = None
    .type = bool
    .help = Switch to a rotomer minimally distant from the current one
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
  .short_caption=Modify occupancies
  .expert_level=2
  .style  = noauto menu_item parent_submenu:model_modifications
{
  randomize = None
    .type = bool
    .help = Randomize occupancies within a certain range
  set = None
    .type = float
    .help = Set all or selected occupancies to given value
    .short_caption=Set occupancies to
}
truncate_to_polyala = None
  .type = bool
  .help = Truncate a model to poly-Ala.
  .style = noauto
set_chemical_element_simple_if_necessary = None
  .type = bool
  .help = Make a simple guess about what the chemical element is (based on \
          atom name and the way how it is formatted) and write it into output file.
rename_chain_id
  .help = Rename chains
{
  old_id = None
    .type = str
  new_id = None
    .type = str
}
output
  .help = Write out PDB file with modified model (file name is defined in \
          write_modified)
  .style = noauto
{
  file_name=None
    .type=path
    .input_size=400
    .short_caption = Output PDB file
    .help = Default is the original file name with the file extension \
            replaced by "_modified.pdb".
    .style = bold
}
random_seed = None
  .type = int
  .help = Random seed
  .expert_level=2
"""
modify_params = iotbx.phil.parse(modify_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
modify
  .caption = These options will only be processed if you selected "Modify atomic coordinates or properties" as the action.
  .short_caption = Modify starting model
  .style = menu_item scrolled auto_align
{
remove = None
  .type = str
  .help = Selection for the atoms to be removed
  .short_caption=Remove atom selection
  .input_size=400
  .style = bold selection
keep = None
  .type = str
  .help = Select atoms to keep
  .short_caption=Keep only atom selection
  .input_size=400
  .style = bold selection
put_into_box_with_buffer = None
  .type = float
  .help = Move molecule into center of box.
  .expert_level = 2
%s
}
input {
  pdb
  {
    include scope mmtbx.utils.pdb_params
  }
  monomer_library {
    include scope mmtbx.utils.cif_params
  }
  crystal_symmetry
    .help = Unit cell and space group parameters
    .short_caption = Crystal symmetry
    .style = auto_align menu_item
  {
    unit_cell=None
      .type=unit_cell
    space_group=None
      .type=space_group
  }
}
pdb_interpretation
  .short_caption = PDB Interpretation
  .style = menu_item
{
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
}
geometry_minimization
  .short_caption = Geometry minimization
  .caption = Note: these options will only be processed if you select \
    "Regularize model geometry" as the action.
  .style = auto_align menu_item
{
  include scope mmtbx.command_line.geometry_minimization.master_params
}
regularize_geometry = False
  .type = bool
  .short_caption = Perform geometry minimization
  .style = bold
"""%modify_params_str, process_includes=True)

class modify(object):
  def __init__(self, xray_structure, params, all_chain_proxies, log = None):
    self.log = log
    self.pdb_hierarchy = all_chain_proxies.pdb_hierarchy
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
          xray_structure    = xray_structure)[0])
    if(params_keep_selection is not None):
      self.keep_selection = flex.smart_selection(
        flags=utils.get_atom_selections(
          iselection        = False,
          all_chain_proxies = all_chain_proxies,
          selection_strings = [params_keep_selection],
          xray_structure    = xray_structure)[0])
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
        xray_structure    = xray_structure)[0])
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
      if (adp.atom_selection is None):
        selection = self.top_selection
      else:
        selection = flex.smart_selection(
          flags=utils.get_atom_selections(
            iselection=False,
            all_chain_proxies=self.all_chain_proxies,
            selection_strings=[adp.atom_selection],
            xray_structure=self.xray_structure)[0])
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
      if (sites.atom_selection is None):
        selection = self.top_selection
      else:
        selection = flex.smart_selection(
          flags=utils.get_atom_selections(
            iselection=False,
            all_chain_proxies=self.all_chain_proxies,
            selection_strings=[sites.atom_selection],
            xray_structure=self.xray_structure)[0])
      self._shake_sites(selection=selection, rms_difference=sites.shake)
      if(sites.max_rotomer_distortion is not None):
        self._max_distant_rotomer(selection=selection)
      if(sites.min_rotomer_distortion is not None):
        self._max_distant_rotomer(selection=selection, min_dist_flag=True)
      self._rb_shift(
        selection=selection,
        translate=sites.translate,
        rotate=sites.rotate,
        euler_angle_convention=sites.euler_angle_convention)

  def _max_distant_rotomer(self, selection, min_dist_flag=False):
    if(min_dist_flag):
      self._print_action(
        text = "Switching to the min distant rotomers",
        selection = selection)
    else:
      self._print_action(
        text = "Switching to the max distant rotomers",
        selection = selection)
    mon_lib_srv = mmtbx.monomer_library.server.server()
    selection = selection.flags
    sites_cart_start = self.xray_structure.sites_cart()
    sites_cart_result = sites_cart_start.deep_copy()
    xrs = self.xray_structure.deep_copy_scatterers()
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          if(len(conformers)>1): continue # XXX ignore alt conformations
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            residue_iselection = flex.size_t()
            exclude = False
            for atom in residue.atoms():
              residue_iselection.append(atom.i_seq)
            for r_i_seq in residue_iselection:
              if(not selection[r_i_seq]):
                exclude = True
                break
            if(not exclude):
              rotamer_iterator = lockit.get_rotamer_iterator(
                mon_lib_srv         = mon_lib_srv,
                residue             = residue,
                atom_selection_bool = None)
              dist_start = -1.
              if(min_dist_flag): dist_start = 1.e+6
              if(rotamer_iterator is not None):
                for rotamer, rotamer_sites_cart in rotamer_iterator:
                  sites_cart_start_ = sites_cart_start.deep_copy()
                  xray_structure_ = self.xray_structure.deep_copy_scatterers()
                  sites_cart_start_ = sites_cart_start_.set_selected(
                    residue_iselection, rotamer_sites_cart)
                  xray_structure_.set_sites_cart(sites_cart_start_)
                  dist = flex.sum(xray_structure_.distances(xrs))
                  flag = None
                  if(min_dist_flag):
                    flag = dist < dist_start and abs(dist-dist_start) > 0.3
                  else:
                    flag = dist > dist_start
                  if(flag):
                    dist_start = dist
                    sites_cart_result = sites_cart_result.set_selected(
                      residue_iselection, rotamer_sites_cart)
    self.xray_structure.set_sites_cart(sites_cart_result)

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
        text = "Setting occupancies to: %8.3f"% \
                self.params.occupancies.set,
        selection = selection)
      if (self._occupancies_modified):
        raise Sorry("Can't modify occupancies (already modified).")
      else:
        self._occupancies_modified = True
      self.xray_structure.set_occupancies(
          value = self.params.occupancies.set,
          selection = selection.flags)

  def report_number_of_atoms_to_be_removed(self):
    if (    self.remove_selection is not None
        and self.remove_selection.selected_size > 0):
      self.remove_selection.show_summary(
        out = self.log,
        label = "Atoms to be kept: ")

def set_chemical_element_simple_if_necessary(hierarchy):
  for model in hierarchy.models():
    for chain in model.chains():
      for atom in chain.atoms():
        atom.set_element(atom.determine_chemical_element_simple())

def rename_chain_id(hierarchy, params, log):
  print >> log, "old_id= '%s'"%params.old_id
  print >> log, "new_id= '%s'"%params.new_id
  counter = 0
  for model in hierarchy.models():
    for chain in model.chains():
      if(chain.id == params.old_id):
        chain.id = params.new_id
        counter += 1
  if(counter==0):
    print >> log, \
    "WARNING: no chain id renamed: check input PDB file or renaming parameters."

def truncate_to_poly_ala(hierarchy):
  import iotbx.pdb.amino_acid_codes
  aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
  ala_atom_names = set([" N  ", " CA ", " C  ", " O  ", " CB "])
  for model in hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        def have_amino_acid():
          for ag in rg.atom_groups():
            if (ag.resname in aa_resnames):
              return True
          return False
        if (have_amino_acid()):
          for ag in rg.atom_groups():
            for atom in ag.atoms():
              if (atom.name not in ala_atom_names):
                ag.remove_atom(atom=atom)

def run(args, command_name="phenix.pdbtools"):
  log = utils.set_log(args)
  utils.print_programs_start_header(
    log  = log,
    text = "  %s tools for PDB model manipulations." % command_name)
  command_line_interpreter = interpreter(command_name  = command_name,
                                         args          = args,
                                         log           = log)
  output_files = []
### Truncate to poly-Ala
  if(command_line_interpreter.params.modify.truncate_to_polyala):
    xray_structure = command_line_interpreter.pdb_inp.xray_structure_simple()
    utils.print_header("Truncating to poly-Ala", out = log)
    pdb_hierarchy = command_line_interpreter.pdb_inp.construct_hierarchy()
    truncate_to_poly_ala(hierarchy = pdb_hierarchy)
    pdbout = os.path.basename(command_line_interpreter.pdb_file_names[0])
    pdb_hierarchy.write_pdb_file(file_name = pdbout+"_modified.pdb",
      crystal_symmetry = xray_structure.crystal_symmetry())
    output_files.append(pdbout+"_modified.pdb")
    return output_files
###
  command_line_interpreter.set_ppf()
  all_chain_proxies = \
    command_line_interpreter.processed_pdb_file.all_chain_proxies
  xray_structure = command_line_interpreter.processed_pdb_file.xray_structure(
    show_summary = False)
  if(xray_structure is None):
    raise Sorry("Cannot extract xray_structure.")
### get i/o file names
  ofn = command_line_interpreter.params.modify.output.file_name
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
      hd_selection       = xray_structure.hd_selection(),
      ignore_hd          = command_line_interpreter.command_line.options.ignore_hydrogens,
      restraints_manager = restraints_manager).show(out = log)
    return None
### show_adp_statistics and exit
  if(command_line_interpreter.command_line.options.show_adp_statistics):
    utils.print_header("ADP statistics", out = log)
    model = mmtbx.model.manager(
      xray_structure = xray_structure,
      pdb_hierarchy = all_chain_proxies.pdb_hierarchy,
      log = log)
    model.show_adp_statistics(out = log, padded = True)
    return None
### add hydrogens and exit
  if(command_line_interpreter.command_line.options.add_h):
    utils.print_header("Adding hydrogen atoms", out = log)
    if(len(ifn) > 1): raise Sorry("Multiple input PDB files found.")
    ifn = ifn[0]
    if(not os.path.isfile(ifn)): raise Sorry("File %s does not exist."%ifn)
    easy_run.go("phenix.reduce %s > %s"% (ifn, ofn))
    print >> log, "Output model file name (with H added): %s\n"%ofn
    return [ofn]
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
### set_chemical_element_simple_if_necessary
  if(command_line_interpreter.params.modify.set_chemical_element_simple_if_necessary):
    xray_structure = command_line_interpreter.pdb_inp.xray_structure_simple()
    utils.print_header("Restore 77-78 column", out = log)
    print >> log, """  Make a simple guess about which chemical element should be in 77-78 column,
  and outputs and new file with this column filled."""
    pdb_hierarchy = command_line_interpreter.pdb_inp.construct_hierarchy()
    set_chemical_element_simple_if_necessary(hierarchy = pdb_hierarchy)
    pdbout = os.path.basename(command_line_interpreter.pdb_file_names[0])
    pdb_hierarchy.write_pdb_file(file_name = pdbout+"_modified.pdb",
      crystal_symmetry = xray_structure.crystal_symmetry())
    output_files.append(pdbout+"_modified.pdb")
    print >> log, "All done."
    return output_files
### rename_chain_id
  if([command_line_interpreter.params.modify.rename_chain_id.old_id,
      command_line_interpreter.params.modify.rename_chain_id.new_id].count(None)==0):
    utils.print_header("Rename chain id", out = log)
    pdb_hierarchy = command_line_interpreter.pdb_inp.construct_hierarchy()
    rename_chain_id(hierarchy = pdb_hierarchy,
      params = command_line_interpreter.params.modify.rename_chain_id, log = log)
    pdbout = os.path.basename(command_line_interpreter.pdb_file_names[0])
    pdb_hierarchy.write_pdb_file(file_name = pdbout+"_modified.pdb",
      crystal_symmetry = xray_structure.crystal_symmetry())
    output_files.append(pdbout+"_modified.pdb")
    print >> log, "All done."
    return output_files
### do other model manipulations
  utils.print_header("Performing requested model manipulations", out = log)
  result = modify(xray_structure    = xray_structure,
                  params            = command_line_interpreter.params.modify,
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
  output_files.append(ofn)
  utils.print_header("Done", out = log)
  return output_files

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    results = run(list(self.args))
    eff_file = self.args[0]
    if os.path.isfile(eff_file) :
      base, ext = os.path.splitext(eff_file)
      easy_pickle.dump("%s.pkl" % base, results)
    return results

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
    self.cif_file_names = []
    self.cif_objects = []
    self.processed_pdb_file = None
    self.processed_pdb_file_reference = None
    self.process_args()
    self.pdb_file_names.extend(self.params.input.pdb.file_name)
    if(self.fake_crystal_symmetry): # create cs for internal use
      pdb_combined = combine_unique_pdb_files(file_names = self.pdb_file_names)
      pdb_combined.report_non_unique(out = self.log)
      if(len(pdb_combined.unique_file_names) == 0):
        raise Sorry("No coordinate file given.")
      pdb_inp = iotbx.pdb.input(
        source_info = None,
        lines       = flex.std_string(pdb_combined.raw_records))
      self.crystal_symmetry = pdb_inp.xray_structure_simple().\
        cubic_unit_cell_around_centered_scatterers(
        buffer_size = 10).crystal_symmetry()
    print >> self.log, "Working crystal symmetry after inspecting all inputs:"
    self.crystal_symmetry.show_summary(f = self.log, prefix="  ")
    self.params.input.crystal_symmetry.unit_cell = \
      self.crystal_symmetry.unit_cell()
    self.params.input.crystal_symmetry.space_group = \
      self.crystal_symmetry.space_group_info()
    self.point_group = self.crystal_symmetry.space_group() \
      .build_derived_point_group()
    print >> self.log
    #
    pdb_combined = combine_unique_pdb_files(file_names = self.pdb_file_names)
    pdb_combined.report_non_unique(out = self.log)
    if(len(pdb_combined.unique_file_names) == 0):
      raise Sorry("No coordinate file given.")
    self.pdb_inp = iotbx.pdb.input(source_info = None,
      lines = flex.std_string(pdb_combined.raw_records))
    #
    for file_name in self.params.input.monomer_library.file_name :
      if os.path.isfile(file_name) and not file_name in self.cif_file_names :
        cif_object = mmtbx.monomer_library.server.read_cif(file_name=file_name)
        if(len(cif_object) > 0):
          self.cif_objects.append((file_name,cif_object))
          self.cif_file_names.append(file_name)

  def set_ppf(self):
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
      .option("--ignore_hydrogens",
          action="store_true",
          help="Do not account for H and D atoms in geometry statistics.")
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
              self.cif_file_names.append(os.path.abspath(arg))
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

def validate_params (params, callback=None) :
  if len(params.input.pdb.file_name) == 0 :
    raise Sorry("No PDB file(s) specified.")
  elif params.action is None :
    raise Sorry("You must choose an action to perform.")
  elif params.action in ["modify", "regularize"] :
    if params.modify.output.file_name is None :
      raise Sorry("Please specify an output PDB file.")
  elif params.action == "add_h" :
    raise Sorry("Please use ReadySet to add hydrogens.")
