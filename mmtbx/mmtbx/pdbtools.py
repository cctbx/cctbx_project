
master_params = iotbx.phil.parse("""\
selection = None
  .type = str
  .multiple = True
  .help = Selection of atoms to be modified
adp
  .help = Scope of options to modify ADP of selected atoms
{
  randomize = None
    .type = bool
    .help = Randomize ADP within a certain range
  set_b_iso_to = None
    .type = float
    .help = Set ADP of atoms to set_b_iso_to
  convert_to_isotropic = False
    .type = bool
    .help = Convert atoms to isotropic
  convert_to_anisotropic = False
    .type = bool
    .help = Convert atoms to isotropic
  shift_b_iso_by = None
    .type = float
    .help = Add shift_b_iso_by values to ADP
  scale_adp = None
    .type = float
    .help = Multiply all ADP by scale_adp
}
sites
  .help = Scope of options to modify coordinates of selected atoms
{
  shake = None
    .type = float
    .help = Randomize coordinates with mean error value equal to shake
  rigid_body
    .help = Shift model as a rigid body
  {
    t = 0 0 0
      .type = strings
      # XXX FUTURE float(3)
      .help = Translational shift
    r = 0 0 0
      .type = strings
      # XXX FUTURE float(3)
      .help = Rotational shift
  }
}
occupancies
  .help = Scope of options to modify occupancies of selected atoms
{
  randomize = None
    .type = bool
    .help = Randomize occupancies within a certain range
}
remove_hydrogens = True
  .type = bool
  .help = Remove hydrogen atoms from imput model
remove_ordered_solvent = False
  .type = bool
  .help = Remove ordered solvent (water) from imput model
write_modified = None
  .type = str
  .help = Write out PDB file with modified model (file name is defined in \
          write_modified)
""")


class modify(object):
  def __init__(self, xray_structure, params, all_chain_proxies, log):
    self.xray_structure = xray_structure
    self.all_chain_proxies = all_chain_proxies
    self.params = params
    selection = get_atom_selections(
                    iselection        = False,
                    all_chain_proxies = all_chain_proxies,
                    selection_strings = params.modify_start_model.selection)[0]
    if(self.params.adp.convert_to_isotropic):
       print >> log, "Converting to isotropic ADP."
       self.xray_structure.convert_to_isotropic(
                                            selection = selection.iselection())
    if(self.params.adp.set_b_iso_to is not None):
       biso = adp.set_b_iso_to
       print >> log, "Setting all ADP to: %8.3f"%biso
       self.xray_structure.set_b_iso(value = biso, selection = selection)
    if(self.params.adp.scale_adp is not None):
       print >> log, "Scaling all ADP with: %.3f"%adp.scale_adp
       self.xray_structure.scale_adp(factor = adp.scale_adp)
    if(self.params.adp.shift_b_iso_by is not None):
       biso = adp.shift_b_iso_by
       print >> log, "Shift all ADP by: %8.3f"%biso
       self.xray_structure.shift_us(b_shift = biso)
    if(self.params.adp.set_biso_to_wilson_b and wilson_b is not None):
       biso = wilson_b
       print >> log, "Setting ADP to Wilson B: %8.3f"%biso
       self.xray_structure.set_b_iso(value = biso, selection = selection)
    if(self.params.adp.randomize):
       print >> log, "Randomizing ADP."
       self.xray_structure.shake_adp(selection = selection)
    if(self.params.sites.shake is not None):
       print >> log, \
         "Shaking sites (RMS difference = %8.3f)." % sites.shake
       self.xray_structure.shake_sites_in_place(
         rms_difference=sites.shake,
         selection=selection)
    rb = self.params.sites.rigid_body
    trans = [float(i) for i in rb.t]
    rot   = [float(i) for i in rb.r]
    rb_shift = (abs(trans[0]+trans[1]+trans[2]+rot[0]+rot[1]+rot[2]) > 1.e-6)
    if(rb_shift):
       print >> log, "Rigid body shift."
       if(params.rigid_body.euler_angle_convention == "zyz"):
          rot_obj = mmtbx.refinement.rigid_body.rb_mat_euler(phi = rot[0],
                                                             psi = rot[1],
                                                             the = rot[2])
       else:
          rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = rot[0],
                                                       psi = rot[1],
                                                       the = rot[2])
       dim = self.xray_structure.scatterers().size()
       self.xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                          xray_structure      = self.xray_structure,
                          rotation_matrices   = [rot_obj.rot_mat()],
                          translation_vectors = [(trans[0],trans[1],trans[2])],
                          selections          = [flex.bool(dim, True)])
    if(self.params.occupancies.randomize):
       print >> log, "Randomizing all occupancies."
       self.xray_structure.shake_occupancies(selection = selection)
