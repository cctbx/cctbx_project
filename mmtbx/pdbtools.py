
"""
Common manipulations of (macromolecular) model files.  This module is used
as a standalone program, phenix.pdbtools, but some components are also used
in phenix.refine, and it is also a dumping ground for individual related
functions.
"""

from __future__ import division
import iotbx.phil
from mmtbx.refinement import rigid_body
from mmtbx import utils
from cctbx.array_family import flex
from iotbx.option_parser import iotbx_option_parser
from libtbx.str_utils import show_string
import libtbx.phil
from iotbx import pdb
from cctbx import crystal
from libtbx.utils import Sorry, null_out, check_if_output_directory_exists
from iotbx import crystal_symmetry_from_any
from mmtbx import monomer_library
import mmtbx.restraints
import mmtbx.model
from mmtbx import model_statistics
import random
from libtbx import easy_run
from iotbx.pdb import combine_unique_pdb_files, write_whole_pdb_file
from libtbx import runtime_utils
import scitbx.matrix
import cStringIO
import os, sys
import scitbx.rigid_body


modify_params_str = """\
selection = None
  .type = atom_selection
  .help = Selection for atoms to be modified
  .short_caption = Modify atom selection
  .input_size=400
  .style = bold noauto
adp
  .help = Scope of options to modify ADP of selected atoms
  .multiple = True
  .short_caption=Modify ADPs
  .style = auto_align menu_item parent_submenu:model_modifications noauto
{
  atom_selection = None
    .type = atom_selection
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
    .short_caption = Modify ADPs for selection
    .input_size=400
    .style  = bold
  randomize = False
    .type = bool
    .help = Randomize ADP within a certain range
    .short_caption=Randomize ADPs
  set_b_iso = None
    .type = float
    .help = Set ADP of atoms to set_b_iso
    .short_caption=Set isotropic B to
    .input_size = 64
  convert_to_isotropic = False
    .type = bool
    .help = Convert atoms to isotropic
  convert_to_anisotropic = False
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
  .style = auto_align noauto menu_item parent_submenu:model_modifications
{
  atom_selection = None
    .type = atom_selection
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
    .input_size=400
    .short_caption = Modify sites for selection
    .style = bold
  shake = None
    .type = float
    .help = Randomize coordinates with mean error value equal to shake
    .short_caption = Randomize coordinates (mean value)
  switch_rotamers = max_distant min_distant exact_match fix_outliers
    .type=choice(multi=False)
  translate = 0 0 0
    .type = floats(size=3)
    .optional = False
    .help = Translational shift (x,y,z)
  rotate = 0 0 0
    .type = floats(size=3)
    .optional = False
    .help = Rotational shift (x,y,z)
  euler_angle_convention = *xyz zyz
    .type = choice
    .help = Euler angles convention to be used for rotation
}
occupancies
  .help = Scope of options to modify occupancies of selected atoms
  .multiple = True
  .short_caption=Modify occupancies
  .style  = noauto menu_item parent_submenu:model_modifications
{
  atom_selection = None
    .type = atom_selection
    .help = Selection for atoms to be modified. \\
            Overrides parent-level selection.
    .input_size=400
    .short_caption = Modify sites for selection
    .style = bold
  randomize = False
    .type = bool
    .help = Randomize occupancies within a certain range
    .short_caption = Randomize occupancies
  set = None
    .type = float
    .help = Set all or selected occupancies to given value
    .short_caption=Set occupancies to
    .input_size = 64
  normalize = False
    .type = bool
    .help = Reset the sum of occupancies for each residue to 1.0.  Note that \
      this will be performed before any atoms are removed.
    .short_caption = Reset total occupancy to 1
}
rotate_about_axis
  .style = box auto_align
{
  axis = None
    .type = str
  angle = None
    .type = float
  atom_selection = None
    .type = str
}
change_of_basis = None
  .type = str
  .short_caption = Change of basis operator
  .help = Apply change-of-basis operator (e.g. reindexing operator) to \
    the coordinates and symmetry.  Example: 'a,c,b'.
renumber_residues = False
  .type = bool
  .help = Re-number residues
increment_resseq = None
  .type = int
  .help = Increment residue number
  .short_caption = Increment residue numbers by
truncate_to_polyala = False
  .type = bool
  .help = Truncate a model to poly-Ala.
  .short_caption = Truncate to poly-Ala
  .style = noauto
remove_alt_confs = False
  .type = bool
  .help = Deletes atoms whose altloc identifier is not blank or 'A', and \
    resets the occupancies of the remaining atoms to 1.0.
  .short_caption = Remove alternate conformers
  .style = noauto
always_keep_one_conformer = False
  .type = bool
  .help = Modifies behavior of remove_alt_confs so that residues with no \
    conformer labeled blank or A are not deleted.  Silent if remove_alt_confs \
    is False.
set_chemical_element_simple_if_necessary = None
  .type = bool
  .short_caption = Guess element field if necessary
  .help = Make a simple guess about what the chemical element is (based on \
          atom name and the way how it is formatted) and write it into output file.
set_seg_id_to_chain_id = False
  .type = bool
  .short_caption = Set segID to chain ID
  .help = Sets the segID field to the chain ID (padded with spaces).
  .style = noauto
clear_seg_id = False
  .type = bool
  .short_caption = Clear segID field
  .help = Erases the segID field.
  .style = noauto
convert_semet_to_met = False
  .type = bool
  .short_caption = Convert SeMet residues to Met
  .style = noauto
rename_chain_id
  .help = Rename chains
  .short_caption = Rename chain ID
  .style = box
{
  old_id = None
    .type = str
    .input_size = 50
    .short_caption = Old ID
  new_id = None
    .type = str
    .input_size = 50
    .short_caption = New ID
}
set_charge
  .short_caption = Set atomic charge
  .style = box auto_align
{
  charge_selection = None
    .type = atom_selection
    .short_caption = Atom selection
  charge = None
    .type = int(value_max=7,value_min=-3)
}
output
  .help = Write out file with modified model (file name is defined in \
          write_modified)
  .style = noauto
{
  file_name=None
    .type=path
    .input_size=400
    .short_caption = Output model file
    .help = Default is the original file name with the file extension \
            replaced by "_modified.pdb".
    .style = bold new_file file_type:pdb
  format = *pdb mmcif
    .type = choice
    .help = Choose the output format of coordinate file (PDB or mmCIF)
}
remove_first_n_atoms_fraction = None
  .short_caption = Remove first N atoms (fraction)
  .type = float
random_seed = None
  .type = int
  .help = Random seed
"""
modify_params = iotbx.phil.parse(modify_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
use_neutron_distances = False
  .type = bool
  .help = Use neutron X-H distances (which are longer than X-ray ones)
modify
  .short_caption = Modify starting model
  .style = menu_item scrolled auto_align
{
remove = None
  .type = atom_selection
  .help = Selection for the atoms to be removed
  .short_caption=Remove atom selection
  .input_size=400
  .style = bold noauto
keep = None
  .type = atom_selection
  .help = Select atoms to keep
  .short_caption=Keep only atom selection
  .input_size=400
  .style = bold noauto
put_into_box_with_buffer = None
  .type = float
  .help = Move molecule into center of box.
%s
move_waters_last = False
  .type = bool
  .short_caption = Move waters to end of model
  .help = Transfer waters to the end of the model.  Addresses some \
    limitations of water picking in phenix.refine.
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
    .style = hidden
  {
    unit_cell=None
      .type=unit_cell
      .style = noauto
    space_group=None
      .type=space_group
      .style = noauto
  }
}
model_statistics = None
  .type = bool
  .style = hidden
stop_for_unknowns = True
  .type = bool
  .short_caption = Stop for residues missing geometry restraints
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
include scope libtbx.phil.interface.tracking_params
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
    self.xray_structure_was_replaced = False
    if(self.params.random_seed is not None):
      random.seed(self.params.random_seed)
      flex.set_random_seed(self.params.random_seed)
    try: params_remove_selection = self.params.remove
    except KeyboardInterrupt: raise
    except Exception: params_remove_selection = None
    try: params_keep_selection = self.params.keep
    except KeyboardInterrupt: raise
    except Exception: params_keep_selection = None
    if (params.change_of_basis is not None) :
      print >> log, "Applying change-of-basis operator '%s'" % \
        params.change_of_basis
      from cctbx import sgtbx
      cb_op = sgtbx.change_of_basis_op(params.change_of_basis)
      self.xray_structure = self.xray_structure.change_basis(cb_op)
      self.pdb_hierarchy.atoms().set_xyz(self.xray_structure.sites_cart())
      print >> log, "New symmetry:"
      self.xray_structure.crystal_symmetry().show_summary(f=log, prefix="  ")
      self.xray_structure_was_replaced = True
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
    self._process_remove_first_n_atoms_fraction()
    self._rotate_about_axis()
    self._process_adp()
    self._process_sites()
    self._process_occupancies()
    try: self._put_in_box()
    except KeyboardInterrupt: raise
    except Exception: pass

  def _process_remove_first_n_atoms_fraction(self):
    if(self.params.remove_first_n_atoms_fraction is not None):
      scatterers = self.xray_structure.scatterers()
      n_remove = int(scatterers.size()*self.params.remove_first_n_atoms_fraction)
      sel = flex.bool(n_remove, False).concatenate(
        flex.bool(scatterers.size()-n_remove, True))
      if(self.remove_selection is not None):
        sel = sel & self.remove_selection.flags
      self.remove_selection = flex.smart_selection(flags=sel)

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
      self._switch_rotamers(selection=selection, mode=sites.switch_rotamers)
      self._rb_shift(
        selection=selection,
        translate=sites.translate,
        rotate=sites.rotate,
        euler_angle_convention=sites.euler_angle_convention)

  def _switch_rotamers(self, selection, mode):
    if(mode is None): return
    self._print_action(
      text = "Switching rotamers; mode = %s"%mode,
      selection = selection)
    mon_lib_srv = mmtbx.monomer_library.server.server()
    self.pdb_hierarchy.atoms().set_xyz(self.xray_structure.sites_cart())
    self.pdb_hierarchy = utils.switch_rotamers(
      pdb_hierarchy=self.pdb_hierarchy,
      mode=mode,
      selection=selection.flags)
    self.xray_structure.set_sites_cart(
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz())

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
        rot_obj = scitbx.rigid_body.rb_mat_zyz(
          phi = rot[0],
          psi = rot[1],
          the = rot[2])
      else:
        rot_obj = scitbx.rigid_body.rb_mat_xyz(
          phi = rot[0],
          psi = rot[1],
          the = rot[2])
      self.xray_structure.apply_rigid_body_shift(
        rot       = rot_obj.rot_mat().as_mat3(),
        trans     = trans,
        selection = selection.indices)

  def _process_occupancies(self):
    def check_if_already_modified() :
      if(self.top_selection): return
      if (self._occupancies_modified):
        raise Sorry("Can't modify occupancies (already modified).")
      else:
        self._occupancies_modified = True
    for occ in self.params.occupancies:
      if(occ.atom_selection is None):
        selection = self.top_selection
      else:
        selection = flex.smart_selection(
          flags=utils.get_atom_selections(
            iselection=False,
            all_chain_proxies=self.all_chain_proxies,
            selection_strings=[occ.atom_selection],
            xray_structure=self.xray_structure)[0])
      if(occ.randomize):
        self._print_action(
          text = "Randomizing occupancies",
          selection = selection)
        check_if_already_modified()
        self.xray_structure.shake_occupancies(selection=selection.flags)
      if(occ.set is not None):
        self._print_action(
          text = "Setting occupancies to: %8.3f"%occ.set, selection = selection)
        check_if_already_modified()
        self.xray_structure.set_occupancies(
            value = occ.set,
            selection = selection.flags)
      # FIXME since this is done *before* any atoms are removed from the
      # structure, it may not always have the desired effect
      if (occ.normalize) :
        self._print_action(
          text = "Resetting total occupancy to 1.0",
          selection = selection)
        check_if_already_modified()
        normalize_occupancies(
          hierarchy=self.pdb_hierarchy,
          selection=selection.flags,
          xray_structure=self.xray_structure,
          log=self.log)

  def report_number_of_atoms_to_be_removed(self):
    if (    self.remove_selection is not None
        and self.remove_selection.selected_size > 0):
      self.remove_selection.show_summary(
        out = self.log,
        label = "Atoms to be kept: ")

  def _rotate_about_axis(self):
    raap = self.params.rotate_about_axis
    sites_cart = self.xray_structure.sites_cart()
    if([raap.axis, raap.atom_selection, raap.angle].count(None)==0):
      axis = []
      try:
        for a in raap.axis.split():
          axis.append(float(a))
      except Exception:
        sel = utils.get_atom_selections(
          iselection=False,
          all_chain_proxies=self.all_chain_proxies,
          selection_strings=raap.axis,
          xray_structure=self.xray_structure)[0]
        axis = [i for i in sites_cart.select(sel).as_double()]
      if(len(axis)!=6):
        raise Sorry("Bad selection rotate_about_axis.axis: %s"%str(raap.axis))
      p1 = scitbx.matrix.col(axis[:3])
      p2 = scitbx.matrix.col(axis[3:])
      raa = p1.rt_for_rotation_around_axis_through(
        point=p2, angle=raap.angle, deg=True)
      #
      sel = utils.get_atom_selections(
          iselection=False,
          all_chain_proxies=self.all_chain_proxies,
          selection_strings=raap.atom_selection,
          xray_structure=self.xray_structure)[0]
      if(sel.count(True)==0):
        raise Sorry(
          "Empty selection rotate_about_axis.selection: %s"%str(raap.atom_selection))
      sites_cart_rotated = raa * sites_cart.select(sel)
      self.xray_structure.set_sites_cart(
        sites_cart.set_selected(sel, sites_cart_rotated))

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

def set_atomic_charge (
    hierarchy,
    xray_structure,
    pdb_atoms,
    selection,
    charge,
    log=None) :
  """
  Set the charge for selected atoms.  Note that both the hierarchy and the
  xray_structure must be passed here - in the context of the pdbtools app,
  the scattering_type attributes will override the atom's element and charge
  when the structure is written out.
  """
  assert isinstance(charge, int)
  if (selection is None) :
    raise Sorry("You must specify an atom selection to apply a charge to.")
  if (abs(charge) >= 10) :
    raise Sorry("The charge must be in the range from -9 to 9.")
  if (log is None) : log = null_out()
  sel_cache = hierarchy.atom_selection_cache()
  isel = sel_cache.selection(selection).iselection()
  if (len(isel) == 0) :
    raise Sorry("No atoms selected for charge modification (selection = %s)." %
      selection)
  if (charge == 0) :
    charge = "  "
  elif (charge < 0) :
    charge = "%1d-" % abs(charge)
  else :
    charge = "%1d+" % charge
  scatterers = xray_structure.scatterers()
  for i_seq in isel :
    atom = pdb_atoms[i_seq]
    atom.set_charge(charge)
    elem_symbol = scatterers[i_seq].element_and_charge_symbols()[0]
    scatterers[i_seq].scattering_type = elem_symbol + charge
    print >> log, "  %s : set charge to %s" % (atom.id_str(), charge)

def truncate_to_poly_ala(hierarchy):
  """
  Remove all protein sidechain atoms beyond C-alpha (as well as all hydrogens).
  Does not change the chemical identity of the residues.
  """
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

def remove_alt_confs (hierarchy, always_keep_one_conformer=False) :
  """
  Remove all alternate conformers from the hierarchy.  Depending on the
  value of always_keep_one_conformer, this will either remove any atom_group
  with altloc other than blank or 'A', or it will remove any atom_group
  beyond the first conformer found.
  """
  for model in hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        atom_groups = residue_group.atom_groups()
        assert (len(atom_groups) > 0)
        if always_keep_one_conformer :
          if (len(atom_groups) == 1) and (atom_groups[0].altloc == '') :
            continue
          atom_groups_and_occupancies = []
          for atom_group in atom_groups :
            if (atom_group.altloc == '') :
              continue
            mean_occ = flex.mean(atom_group.atoms().extract_occ())
            atom_groups_and_occupancies.append((atom_group, mean_occ))
          atom_groups_and_occupancies.sort(lambda a,b: cmp(b[1], a[1]))
          for atom_group, occ in atom_groups_and_occupancies[1:] :
            residue_group.remove_atom_group(atom_group=atom_group)
          single_conf, occ = atom_groups_and_occupancies[0]
          single_conf.altloc = ''
        else :
          for atom_group in atom_groups :
            if (not atom_group.altloc in ["", "A"]) :
              residue_group.remove_atom_group(atom_group=atom_group)
            else :
              atom_group.altloc = ""
          if (len(residue_group.atom_groups()) == 0) :
            chain.remove_residue_group(residue_group=residue_group)
      if (len(chain.residue_groups()) == 0) :
        model.remove_chain(chain=chain)
  atoms = hierarchy.atoms()
  new_occ = flex.double(atoms.size(), 1.0)
  atoms.set_occ(new_occ)

def normalize_occupancies (hierarchy, selection=None,
    xray_structure=None, log=None) :
  """
  Reset the sum of occupancies for each residue_group to 1.0.  This will fail
  if any (selected) atom_group has non-uniform occupancies.
  """
  if (log is None) : log = null_out()
  for model in hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        i_seqs = residue_group.atoms().extract_i_seq()
        if (selection is not None) :
          selection_rg = selection.select(i_seqs)
          if (not selection_rg.all_eq(True)) :
            continue
        atom_groups = residue_group.atom_groups()
        occupancy_sum = 0
        for atom_group in atom_groups :
          atom_occ = atom_group.atoms().extract_occ()
          if (len(atom_occ) == 0) :
            residue_group.remove_atom_group(atom_group)
            continue
          if (atom_group.altloc.strip() == '') :
            continue
          if (not atom_occ.all_eq(atom_occ[0])) :
            raise Sorry("Non-uniform occupancies for atom_group '%s'" %
              atom_group.id_str())
          occupancy_sum += atom_occ[0]
        if (occupancy_sum != 1.0) :
          print >> log, "  residue_group '%s': occupancy=%.2f" % \
            (residue_group.id_str(), occupancy_sum)
        gap = 1.0 - occupancy_sum
        reset = False
        n_groups = len(atom_groups)
        for atom_group in atom_groups :
          atoms = atom_group.atoms()
          if (atom_group.altloc.strip() == '') or (n_groups == 1.0) :
            print >> log, "    altloc ' ': set occupancy to 1.0"
            atoms.set_occ(flex.double(len(atoms), 1.0))
          elif (not reset) :
            occ = atoms.extract_occ()
            atoms.set_occ(occ + flex.double(len(atoms), gap))
            print >> log, "    altloc '%s': set occupancy to %.2f" % \
              (atom_group.altloc, atoms[0].occ)
            reset = True
  occ = hierarchy.atoms().extract_occ()
  if (xray_structure is not None) :
    xray_structure.scatterers().set_occupancies(occ)

def convert_semet_to_met (pdb_hierarchy, xray_structure) :
  """
  Change the chemical identity of MSE residues to MET, without changing
  coordinates.
  """
  n_mse = 0
  scatterers = xray_structure.scatterers()
  for i_seq, atom in enumerate(pdb_hierarchy.atoms()) :
    if (atom.name.strip() == "SE") and (atom.element.upper() == "SE") :
      atom_group = atom.parent()
      if (atom_group.resname == "MSE") :
        n_mse += 1
        atom_group.resname = "MET"
        atom.name = " SD "
        atom.element = " S"
        scatterer = scatterers[i_seq]
        scatterer.scattering_type = 'S'
        for ag_atom in atom_group.atoms() :
          ag_atom.hetero = False
  return n_mse

def renumber_residues(pdb_hierarchy, renumber_from=None,
    atom_selection=None, log=None):
  if (log is None) : log = null_out()
  selected_i_seqs = None
  if (atom_selection is not None) :
    sel_cache = pdb_hierarchy.atom_selection_cache()
    selected_i_seqs = sel_cache.selection(atom_selection).iselection()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      if (selected_i_seqs is not None) :
        chain_i_seqs = chain.atoms().extract_i_seq()
        intersection = selected_i_seqs.intersection(chain_i_seqs)
        if (len(intersection) == 0) :
          continue
        elif (len(intersection) != len(chain_i_seqs)) :
          print >> log, "Warning: chain '%s' is only partially selected (%d out of %d) - will not renumber." % (chain.id, len(intersection), len(chain_i_seqs))
          continue
      if (renumber_from is None) :
        counter = 1
        for rg in chain.residue_groups():
          rg.resseq=counter
          counter += 1
      else :
        for rg in chain.residue_groups() :
          resseq = rg.resseq_as_int()
          resseq += renumber_from
          rg.resseq = "%4d" % resseq

def write_cif_file(pdb_hierarchy, crystal_symmetry, file_name):
  pdb_hierarchy.write_mmcif_file(
    file_name=file_name, crystal_symmetry=crystal_symmetry)

# XXX only works for one MODEL at present
def move_waters (pdb_hierarchy, xray_structure, out) :
  """
  Re-arranges the PDB hierarchy so that waters are located at the end of the
  atom list.  This is done automatically in phenix.refine if water picking is
  enabled, but the routine here is more flexible.  (Note that there is some
  overlap in functionality with mmtbx.command_line.sort_hetatms, but the
  latter module is significantly more complicated.)
  """
  if (len(pdb_hierarchy.models()) > 1) :
    raise Sorry("Rearranging water molecules is not supported for "+
      "multi-MODEL structures.")
  sel_cache = pdb_hierarchy.atom_selection_cache()
  water_sel = sel_cache.selection("resname HOH or resname WAT")
  n_waters = water_sel.count(True)
  if (n_waters == 0) :
    print >> out, "No waters found, skipping"
    return pdb_hierarchy, xray_structure
  else :
    print >> out, "%d atoms will be moved." % n_waters
    hierarchy_water = pdb_hierarchy.select(water_sel)
    hierarchy_non_water = pdb_hierarchy.select(~water_sel)
    xrs_water = xray_structure.select(water_sel)
    xrs_non_water = xray_structure.select(~water_sel)
    for chain in hierarchy_water.only_model().chains() :
      hierarchy_non_water.only_model().append_chain(chain.detached_copy())
    xrs_non_water.add_scatterers(xrs_water.scatterers())
    assert len(xrs_non_water.scatterers()) == \
           len(hierarchy_non_water.atoms()) == len(pdb_hierarchy.atoms())
    for atom, scat in zip(hierarchy_non_water.atoms(),
                          xrs_non_water.scatterers()) :
      assert (atom.id_str() == scat.label)
    return hierarchy_non_water, xrs_non_water

def run(args, command_name="phenix.pdbtools", out=sys.stdout,
    replace_stderr=True):
  log = utils.set_log(args, out=out, replace_stderr=replace_stderr)
  utils.print_programs_start_header(
    log  = log,
    text = "  %s tools for PDB model manipulations." % command_name)
  command_line_interpreter = interpreter(command_name  = command_name,
                                         args          = args,
                                         log           = log)
  output_files = []
  ### get i/o file names
  params = command_line_interpreter.params
  ofn = params.modify.output.file_name
  output_format = params.modify.output.format
  ifn = params.input.pdb.file_name
  if(ofn is None):
    if output_format == "pdb": ext = "pdb"
    elif output_format == "mmcif": ext = "cif"
    if(len(ifn)==1): ofn = os.path.basename(ifn[0]) + "_modified."+ext
    elif(len(ifn)>1): ofn = os.path.basename(ifn[0]) + "_et_al_modified"+ext
    else:
      pdbout = os.path.basename(command_line_interpreter.pdb_file_names[0])
      ofn = pdbout+"_modified."+ext
  command_line_interpreter.set_ppf()
  all_chain_proxies = \
    command_line_interpreter.processed_pdb_file.all_chain_proxies
  xray_structure = command_line_interpreter.processed_pdb_file.xray_structure(
    show_summary = False)
  pdb_hierarchy = all_chain_proxies.pdb_hierarchy
  pdb_atoms = all_chain_proxies.pdb_atoms
  if(xray_structure is None):
    raise Sorry("Cannot extract xray_structure.")
### show_geometry_statistics and exit
  if(params.model_statistics):
    use_molprobity = True
    rotarama_dir = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data",
      test=os.path.isdir)
    if (rotarama_dir is None) or (not libtbx.env.has_module("probe")) :
      use_molprobity = False
    utils.print_header("Geometry statistics", out = log)
    geometry = command_line_interpreter.processed_pdb_file.\
      geometry_restraints_manager(show_energies = True,
        params_edits=params.geometry_restraints.edits)
    restraints_manager = mmtbx.restraints.manager(
      geometry = geometry,
      normalization = True)
    pdb_hierarchy = command_line_interpreter.pdb_inp.construct_hierarchy()
    ph = pdb_hierarchy
    rm = restraints_manager
    not_hd_sel = ~xray_structure.hd_selection()
    if(command_line_interpreter.command_line.options.ignore_hydrogens):
      ph = ph.select(not_hd_sel)
      rm = rm.select(not_hd_sel)
    mso = model_statistics.geometry(
      pdb_hierarchy      = ph,
      molprobity_scores  = use_molprobity,
      restraints_manager = rm)
    print >> log
    mso.show(out = log, prefix="")
    utils.print_header("ADP statistics", out = log)
    model = mmtbx.model.manager(
      xray_structure     = xray_structure,
      pdb_hierarchy      = all_chain_proxies.pdb_hierarchy,
      restraints_manager = restraints_manager,
      log                = log)
    model.show_adp_statistics(out = log, padded = True)
    print >> log
    print >> log, "Mean(|Bi-Bj|): %-6.2f"%model.rms_b_iso_or_b_equiv_bonded()
    print >> log, "  Bi and Bj are B-factors of bonded atoms i and j"
    print >> log
    return None
### add hydrogens and exit
  if(command_line_interpreter.command_line.options.add_h):
    if output_format == "mmcif":
      raise Sorry("mmcif format not currently supported with the option add_h")
    utils.print_header("Adding hydrogen atoms", out = log)
    if(len(ifn) > 1): raise Sorry("Multiple input PDB files found.")
    ifn = command_line_interpreter.pdb_file_names[0]
    if(not os.path.isfile(ifn)): raise Sorry("File %s does not exist."%ifn)
    easy_run.go("phenix.reduce %s > %s"% (ifn, ofn))
    print >> log, "Output model file name (with H added): %s\n"%ofn
    return [ofn]
### show parameters
  utils.print_header("Complete set of parameters", out = log)
  master_params.format(params).show(out = log)
### set_chemical_element_simple_if_necessary
  if(params.modify.set_chemical_element_simple_if_necessary):
    utils.print_header("Restore 77-78 column", out = log)
    print >> log, """\
  Make a simple guess about which chemical element should be in 77-78 column,
  and outputs and new file with this column filled."""
    pdb_hierarchy.atoms().set_chemical_element_simple_if_necessary()
### rename_chain_id
  if([params.modify.rename_chain_id.old_id,
      params.modify.rename_chain_id.new_id].count(None)==0):
    utils.print_header("Rename chain id", out = log)
    rename_chain_id(hierarchy = pdb_hierarchy,
      params=params.modify.rename_chain_id,
      log = log)
### Renumber residues
  if (params.modify.increment_resseq) or (params.modify.renumber_residues) :
    utils.print_header("Re-numbering residues", out = log)
    renumber_residues(pdb_hierarchy = pdb_hierarchy,
      renumber_from=params.modify.increment_resseq,
      atom_selection=params.modify.selection,
      log=log)
### segID manipulations
  if (params.modify.set_seg_id_to_chain_id) :
    if (params.modify.clear_seg_id) :
      raise Sorry("Parameter conflict - set_seg_id_to_chain_id=True and "+
        "clear_seg_id=True.  Please choose only one of these options.")
    for atom in pdb_hierarchy.atoms() :
      labels = atom.fetch_labels()
      atom.segid = "%-4s" % labels.chain_id
  elif (params.modify.clear_seg_id) :
    for atom in pdb_hierarchy.atoms() :
      atom.segid = "    "
### convert MSE to MET
  if (params.modify.convert_semet_to_met) :
    convert_semet_to_met(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure)
### set atomic charge
  if (params.modify.set_charge.charge_selection is not None) :
    set_atomic_charge(
      hierarchy=pdb_hierarchy,
      pdb_atoms=pdb_atoms,
      xray_structure=xray_structure,
      selection=params.modify.set_charge.charge_selection,
      charge=params.modify.set_charge.charge,
      log=log)
### do other model manipulations
  utils.print_header("Performing requested model manipulations", out = log)
  result = modify(xray_structure    = xray_structure,
                  params            = params.modify,
                  all_chain_proxies = all_chain_proxies,
                  log               = log)
  result.report_number_of_atoms_to_be_removed()
  if (result.xray_structure_was_replaced) :
    xray_structure = result.xray_structure
  selection = getattr(result.remove_selection, "flags", None)
  if selection is not None:
    xray_structure = xray_structure.select(selection)
    pdb_hierarchy = pdb_hierarchy.select(selection)
### Truncate to poly-ALA (keeping original residue names)
  if(params.modify.truncate_to_polyala):
    utils.print_header("Truncating to poly-Ala", out = log)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    truncate_to_poly_ala(hierarchy=pdb_hierarchy)
    pdb_atoms = pdb_hierarchy.atoms()
    xray_structure = xray_structure.select(pdb_atoms.extract_i_seq())
### Remove alt. confs.
  if (params.modify.remove_alt_confs) :
    utils.print_header("Removing alternate conformations", out = log)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    remove_alt_confs(hierarchy = pdb_hierarchy,
      always_keep_one_conformer=params.modify.always_keep_one_conformer)
    pdb_atoms = pdb_hierarchy.atoms()
    xray_structure = xray_structure.select(pdb_atoms.extract_i_seq())
    xray_structure.set_occupancies(1.0)
    print >> log, "All occupancies reset to 1.0."
  if (params.modify.move_waters_last) :
    utils.print_header("Moving waters to end of model", out=log)
    pdb_hierarchy, xray_structure = move_waters(pdb_hierarchy, xray_structure,
      out=log)
  utils.print_header("Writing output model", out = log)
### write output file (if got to this point)
  print >> log, "Output model file name: ", ofn
  pdb_hierarchy.adopt_xray_structure(xray_structure,
    assert_identical_id_str=False)
  crystal_symmetry = xray_structure.crystal_symmetry()
  if command_line_interpreter.fake_crystal_symmetry:
    crystal_symmetry = None
  ss_ann = None
  if hasattr(command_line_interpreter.pdb_inp, "extract_secondary_structure"):
    ss_ann = command_line_interpreter.pdb_inp.extract_secondary_structure()
  if output_format == "pdb":
    write_whole_pdb_file(
        file_name=ofn,
        pdb_hierarchy=pdb_hierarchy,
        crystal_symmetry=crystal_symmetry,
        append_end=True,
        atoms_reset_serial_first_value=1,
        ss_annotation=ss_ann)
  elif output_format =="mmcif":
    write_cif_file(pdb_hierarchy, crystal_symmetry, file_name=ofn)
  output_files.append(ofn)
  utils.print_header("Done", out = log)
  return output_files

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
    self.process_args()
    self.pdb_file_names.extend(self.params.input.pdb.file_name)
    if(self.fake_crystal_symmetry): # create cs for internal use
      pdb_combined = combine_unique_pdb_files(file_names = self.pdb_file_names)
      pdb_combined.report_non_unique(out = self.log)
      if(len(pdb_combined.unique_file_names) == 0):
        raise Sorry("No coordinate file given.")
      pdb_inp = iotbx.pdb.input(
        source_info = None,
        lines       = flex.std_string(pdb_combined.raw_records),
        raise_sorry_if_format_error = True)
      #XXXself.crystal_symmetry = pdb_inp.xray_structure_simple().\
      #XXX  cubic_unit_cell_around_centered_scatterers(
      #XXX  buffer_size = 10).crystal_symmetry()
      from cctbx import uctbx
      self.crystal_symmetry = \
        uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
          sites_cart=pdb_inp.atoms().extract_xyz(),
          buffer_layer=5).crystal_symmetry()
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
      lines = flex.std_string(pdb_combined.raw_records),
      raise_sorry_if_format_error=True)
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
      stop_for_unknowns         = self.params.stop_for_unknowns,
      log                       = self.log,
      mon_lib_srv               = self.mon_lib_srv,
      use_neutron_distances     = self.params.use_neutron_distances,
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
    phil_out = cStringIO.StringIO()
    master_params.show(out=phil_out, prefix="  ")
    phil_help_str = phil_out.getvalue()
    self.command_line = (iotbx_option_parser(
      usage="%s [options] [pdb_file] [parameter_file]" % self.command_name,
      description='Example: %s model.pdb parameters.txt\n\nFull parameters:\n%s'
        % (self.command_name + description_see_also, phil_help_str))
      .enable_show_defaults()
      .enable_symmetry_comprehensive()
      .option(None, "--unused_ok",
        action="store_true",
        default=False,
        help="Disables detection of unused parameter definitions")
      .option(None, "--quiet",
        action="store_true",
        help="Suppress output to screen")
      .option("--add_h",
          action="store_true",
          help="Add H atoms to a model using Reduce program.")
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
    parameter_interpreter = master_params.command_line_argument_interpreter(
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
        elif (pdb.is_pdb_file(file_name=arg) or
              pdb.is_pdb_mmcif_file(file_name=arg)):
          # XXX this is an inefficient way to do this - a CIF file could get
          # read 3 times, once by pdb.is_pdb_mmcif_file, once by
          # crystal_symmetry_from_cif.extract_from and finally once again when
          # we actually try to read the file properly
          def extract_from(file_name):
            from iotbx.cif import crystal_symmetry_from_cif
            from iotbx.pdb import crystal_symmetry_from_pdb
            for extract_function in (crystal_symmetry_from_pdb.extract_from,
                                     crystal_symmetry_from_cif.extract_from):
              try: return extract_function(file_name)
              except Exception: continue
          self.pdb_file_names.append(arg)
          arg_is_processed = True
          crystal_symmetry_from_any.extract_and_append(
            file_names=[arg],
            target_list=crystal_symmetries_from_coordinate_file,
            extract_function=extract_from)
        else:
          try: cif_object = mmtbx.monomer_library.server.read_cif(file_name=arg)
          except KeyboardInterrupt: raise
          except Exception: pass
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
            raise Sorry("File not found: %s" % show_string(arg))
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
  if (len(params.input.pdb.file_name) == 0) :
    raise Sorry("No PDB file(s) specified.")
  elif (params.modify.output.file_name is not None) :
    if (os.path.isdir(params.modify.output.file_name)) :
      raise Sorry("The specified output file is a currently existing "+
        "directory.")
    check_if_output_directory_exists(params.modify.output.file_name)
  return True

def finish_job (result) :
  output_files = []
  if isinstance(result, list) :
    for file_name in result :
      file_desc = None
      base, ext = os.path.splitext(file_name)
      if ext == ".pdb" :
        file_desc = "Modified PDB file"
      else :
        file_desc = "Unknown file"
      output_files.append((file_desc, file_name))
  return (output_files, [])

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(args=list(self.args))
