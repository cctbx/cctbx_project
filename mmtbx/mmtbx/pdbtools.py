import iotbx.phil
from mmtbx.refinement import rigid_body
from mmtbx import utils
from cctbx.array_family import flex

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
  convert_to_isotropic = False
    .type = bool
    .help = Convert atoms to isotropic
  convert_to_anisotropic = False
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
write_modified = None
  .type = str
  .help = Write out PDB file with modified model (file name is defined in \
          write_modified)
"""
modify_params = iotbx.phil.parse(modify_params_str, process_includes=True)

master_params = iotbx.phil.parse("""\
%s
remove {
  selection = None
    .type = str
    .help = Select atoms to be removed
}
"""%modify_params_str, process_includes=True)

class modify(object):
  def __init__(self, xray_structure, params, all_chain_proxies, log = None):
    self.log = log
    if(self.log is None): self.log = sys.stdout
    self.xray_structure = xray_structure
    self.all_chain_proxies = all_chain_proxies
    self.params = params
    self.remove_selection = None
    try:
      params_remove_selection = self.params.remove.selection
    except: params_remove_selection = "None"
    if(params_remove_selection != "None"):
      self.remove_selection = utils.get_atom_selections(
                         iselection        = False,
                         all_chain_proxies = all_chain_proxies,
                         selection_strings = [self.params.remove.selection])[0]
    self.selection = utils.get_atom_selections(
                                       iselection        = False,
                                       all_chain_proxies = all_chain_proxies,
                                       selection_strings = params.selection)[0]
    self._convert_to_isotropic()
    self._convert_to_anisotropic()
    self._set_b_iso()
    self._scale_adp()
    self._shift_b_iso()
    self._randomize_adp()
    self._shake_sites()
    self._rb_shift()
    self._randomize_occupancies()

  def _print_action(self, text):
    print >> self.log, text

  def _convert_to_isotropic(self):
    if(self.params.adp.convert_to_isotropic):
       self._print_action(text = "Converting to isotropic ADP.")
       self.xray_structure.convert_to_isotropic(
                                       selection = self.selection.iselection())

  def _convert_to_anisotropic(self):
    if(self.params.adp.convert_to_anisotropic):
       self._print_action(text = "Converting to anisotropic ADP.")
       self.xray_structure.convert_to_anisotropic(
                                       selection = self.selection.iselection())

  def _set_b_iso(self):
    b_iso = self.params.adp.set_b_iso
    if(b_iso is not None):
       self._print_action(text = "Setting all isotropic ADP to: %8.3f"%b_iso)
       self.xray_structure.set_b_iso(value = b_iso, selection = self.selection)

  def _scale_adp(self):
    scale = self.params.adp.scale_adp
    if(scale is not None):
       self._print_action(text = "Scaling all ADP with: %.3f"%scale)
       self.xray_structure.scale_adp(factor = scale, selection= self.selection)

  def _shift_b_iso(self):
    shift = self.params.adp.shift_b_iso
    if(shift is not None):
       self._print_action(text = "Shift all ADP by: %8.3f"%shift)
       self.xray_structure.shift_us(b_shift   = shift,
                                    selection = self.selection.iselection())

  def _randomize_adp(self):
    if(self.params.adp.randomize):
       self._print_action(text = "Randomizing ADP.")
       self.xray_structure.shake_adp(selection = self.selection)

  def _shake_sites(self):
    shake_value = self.params.sites.shake
    if(shake_value is not None):
       self._print_action(
                  text = "Shaking sites (RMS difference = %8.3f)."%shake_value)
       self.xray_structure.shake_sites_in_place(rms_difference= shake_value,
                                                selection     = self.selection)

  def _rb_shift(self):
    rb = self.params.sites
    trans = [float(i) for i in rb.translate]
    rot   = [float(i) for i in rb.rotate]
    rb_shift = (abs(trans[0]+trans[1]+trans[2]+rot[0]+rot[1]+rot[2]) > 1.e-6)
    if(rb_shift):
       self._print_action(text = "Rigid body shift.")
       if(self.params.sites.euler_angle_convention == "zyz"):
          rot_obj = rigid_body.rb_mat_euler(phi = rot[0],
                                            psi = rot[1],
                                            the = rot[2])
       else:
          rot_obj = rigid_body.rb_mat(phi = rot[0],
                                      psi = rot[1],
                                      the = rot[2])
       print self.selection
       self.xray_structure.apply_rigid_body_shift(
                                    rot       = rot_obj.rot_mat().as_mat3(),
                                    trans     = trans,
                                    selection = self.selection.iselection())

  def _randomize_occupancies(self):
    if(self.params.occupancies.randomize):
       self._print_action(text = "Randomizing all occupancies.")
       self.xray_structure.shake_occupancies(selection = self.selection)
