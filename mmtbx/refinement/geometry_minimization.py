from __future__ import absolute_import, division, print_function
from cctbx import geometry_restraints
import cctbx.geometry_restraints.lbfgs
import mmtbx.refinement.minimization_ncs_constraints
import scitbx.lbfgs
import sys
import mmtbx.utils
from scitbx.array_family import flex
from mmtbx import monomer_library
from six.moves import range

class lbfgs(geometry_restraints.lbfgs.lbfgs):

  def __init__(self,
        sites_cart,
        geometry_restraints_manager,
        geometry_restraints_flags,
        lbfgs_termination_params,
        correct_special_position_tolerance,
        riding_h_manager=None,
        sites_cart_selection=None,
        lbfgs_exception_handling_params=None,
        rmsd_bonds_termination_cutoff=0,
        rmsd_angles_termination_cutoff=0,
        states_collector=None,
        site_labels=None):
    self.rmsd_bonds_termination_cutoff = rmsd_bonds_termination_cutoff
    self.rmsd_angles_termination_cutoff = rmsd_angles_termination_cutoff
    self.states_collector = states_collector
    geometry_restraints.lbfgs.lbfgs.__init__(self,
      sites_cart=sites_cart,
      riding_h_manager = riding_h_manager,
      correct_special_position_tolerance=correct_special_position_tolerance,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=lbfgs_termination_params,
      sites_cart_selection=sites_cart_selection,
      lbfgs_exception_handling_params=lbfgs_exception_handling_params,
      site_labels=site_labels,
      states_collector=states_collector)

  def callback_after_step(self, minimizer):
    self.apply_shifts()
    if([self.rmsd_angles, self.rmsd_bonds].count(None) == 0):
      if(self.rmsd_angles < self.rmsd_angles_termination_cutoff and
         self.rmsd_bonds < self.rmsd_bonds_termination_cutoff):
        return True

def add_rotamer_restraints(
      pdb_hierarchy,
      restraints_manager,
      selection,
      sigma,
      mode=None,
      accept_allowed=True,
      mon_lib_srv=None,
      rotamer_manager=None):
  if(mode is not None):
    pdb_hierarchy_for_proxies = mmtbx.utils.switch_rotamers(
      pdb_hierarchy  = pdb_hierarchy.deep_copy(),
      mode           = mode,
      accept_allowed = accept_allowed,
      selection      = selection,
      mon_lib_srv    = mon_lib_srv,
      rotamer_manager= rotamer_manager)
  else:
    pdb_hierarchy_for_proxies = pdb_hierarchy.deep_copy()
  sc = pdb_hierarchy_for_proxies.atoms().extract_xyz()
  if selection is not None:
    sc = sc.select(selection)
  restraints_manager.geometry.add_chi_torsion_restraints_in_place(
    pdb_hierarchy   = pdb_hierarchy_for_proxies,
    sites_cart      = sc,
    selection       = selection,
    chi_angles_only = True,
    sigma           = sigma)
  return pdb_hierarchy_for_proxies, restraints_manager

class run2(object):
  def __init__(self,
               restraints_manager,
               pdb_hierarchy,
               correct_special_position_tolerance,
               riding_h_manager               = None,
               ncs_restraints_group_list      = [], # These are actually for NCS CONSTRAINTS!
               max_number_of_iterations       = 500,
               number_of_macro_cycles         = 5,
               selection                      = None,
               bond                           = False,
               nonbonded                      = False,
               angle                          = False,
               dihedral                       = False,
               chirality                      = False,
               planarity                      = False,
               parallelity                    = False,
               rmsd_bonds_termination_cutoff  = 0,
               rmsd_angles_termination_cutoff = 0,
               alternate_nonbonded_off_on     = False,
               cdl                            = False,
               rdl                            = False,
               correct_hydrogens              = False,
               fix_rotamer_outliers           = True,
               allow_allowed_rotamers         = True,
               states_collector               = None,
               log                            = None,
               mon_lib_srv                    = None,
               ias_selection                  = None,
               ):
    self.log = log
    if self.log is None:
      self.log = sys.stdout
    self.pdb_hierarchy = pdb_hierarchy
    self.ias_selection = ias_selection
    self.minimized = None
    self.mon_lib_srv = mon_lib_srv
    if self.mon_lib_srv is None:
      self.mon_lib_srv = monomer_library.server.server()
    self.restraints_manager = restraints_manager
    assert max_number_of_iterations+number_of_macro_cycles > 0
    assert [bond,nonbonded,angle,dihedral,chirality,planarity,
            parallelity].count(False) < 7
    self.cdl_proxies = None
    self.rdl_proxies = None
    self.rotamer_manager = None
    if fix_rotamer_outliers:
      from mmtbx.rotamer.rotamer_eval import RotamerEval
      self.rotamer_manager = RotamerEval(mon_lib_srv=self.mon_lib_srv)
    if(cdl):
      from mmtbx.conformation_dependent_library.cdl_setup import setup_restraints
      self.cdl_proxies = setup_restraints(self.restraints_manager.geometry)
    self.correct_hydrogens = correct_hydrogens
    if(alternate_nonbonded_off_on and number_of_macro_cycles % 2 != 0):
      number_of_macro_cycles += 1
    import scitbx.lbfgs
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_number_of_iterations)
    exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True)
    geometry_restraints_flags = geometry_restraints.flags.flags(
      bond               = bond,
      nonbonded          = nonbonded,
      angle              = angle,
      dihedral           = dihedral,
      chirality          = chirality,
      planarity          = planarity,
      parallelity        = parallelity,
      reference_coordinate = True,
      reference_dihedral = True,
      bond_similarity    = True,
      ramachandran_restraints = True)
    self.update_cdl_restraints()
    self.show()
    for i_macro_cycle in range(number_of_macro_cycles):
      print("  macro-cycle:", i_macro_cycle, file=self.log)
      self.restraints_manager.geometry.update_ramachandran_restraints_phi_psi_targets(
        hierarchy=self.pdb_hierarchy)
      if(alternate_nonbonded_off_on and i_macro_cycle<=number_of_macro_cycles/2):
        geometry_restraints_flags.nonbonded = bool(i_macro_cycle % 2)
      self.update_cdl_restraints(macro_cycle=i_macro_cycle)
      if(fix_rotamer_outliers):
        junk, self.restraints_manager = add_rotamer_restraints(
          pdb_hierarchy      = self.pdb_hierarchy,
          restraints_manager = self.restraints_manager,
          selection          = selection,
          sigma              = 10,
          mode               = "fix_outliers",
          accept_allowed     = allow_allowed_rotamers,
          mon_lib_srv        = self.mon_lib_srv,
          rotamer_manager    = self.rotamer_manager)
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      if rdl:
        self.updaterdl(prefix="Update RDL restraints")
      if (ncs_restraints_group_list is not None
          and len(ncs_restraints_group_list)) > 0:
        # do ncs minimization
        print("Using NCS constraints.", file=self.log)
        xrs = self.pdb_hierarchy.extract_xray_structure().deep_copy_scatterers()
        refine_selection = flex.size_t(range(xrs.scatterers().size()))
        tfg_obj = mmtbx.refinement.minimization_ncs_constraints.\
            target_function_and_grads_geometry_minimization(
                xray_structure=xrs,
                ncs_restraints_group_list=ncs_restraints_group_list, # CONSTRAINTS
                refine_selection=refine_selection,
                restraints_manager=self.restraints_manager.geometry,
                refine_sites=True,
                refine_transformations=False,
                )
        minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
          target_and_grads_object      = tfg_obj,
          xray_structure               = xrs,
          ncs_restraints_group_list    = ncs_restraints_group_list, # CONSTRAINTS
          refine_selection             = refine_selection,
          finite_grad_differences_test = False,
          max_iterations               = max_number_of_iterations,
          refine_sites                 = True,
          refine_transformations       = False)
        self.pdb_hierarchy.adopt_xray_structure(xrs)
      else:
        sites_cart_orig = sites_cart.deep_copy()
        if ias_selection is not None and ias_selection.count(True) > 0:
          sites_cart = sites_cart.select(~ias_selection)
        self.minimized = lbfgs(
          sites_cart                      = sites_cart,
          riding_h_manager                = riding_h_manager,
          correct_special_position_tolerance=correct_special_position_tolerance,
          geometry_restraints_manager     = restraints_manager.geometry,
          geometry_restraints_flags       = geometry_restraints_flags,
          lbfgs_termination_params        = lbfgs_termination_params,
          lbfgs_exception_handling_params = exception_handling_params,
          sites_cart_selection            = selection,
          rmsd_bonds_termination_cutoff   = rmsd_bonds_termination_cutoff,
          rmsd_angles_termination_cutoff  = rmsd_angles_termination_cutoff,
          states_collector                = states_collector,
          site_labels                     = None)
        if(ias_selection is not None):
          for i_seq, ias_s in enumerate(ias_selection): # assumes that IAS appended to the back
            if(not ias_s):
              sites_cart_orig[i_seq] = sites_cart[i_seq]
        else:
          sites_cart_orig = sites_cart
        self.pdb_hierarchy.atoms().set_xyz(sites_cart_orig)
      self.show()
      self.log.flush()
      geometry_restraints_flags.nonbonded = nonbonded
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = max_number_of_iterations)

  def update_cdl_restraints(self, macro_cycle=None):
    if(self.cdl_proxies is not None):
      from mmtbx.conformation_dependent_library import update_restraints
      if(macro_cycle is None):
        rc = update_restraints(
          self.pdb_hierarchy,
          self.restraints_manager.geometry,
          cdl_proxies=self.cdl_proxies,
          log=self.log,
          verbose=False)
      elif(macro_cycle>0):
        rc = update_restraints(
          self.pdb_hierarchy,
          self.restraints_manager.geometry,
          sites_cart=self.pdb_hierarchy.atoms().extract_xyz(),
          cdl_proxies=self.cdl_proxies,
          log=self.log,
          verbose=False)

  def updaterdl(self, prefix):
    if self.restraints_manager is None: return
    from mmtbx.conformation_dependent_library import rotamers
    from mmtbx.refinement import print_statistics
    print_statistics.make_header(prefix, out=self.log)
    self.rdl_proxies = None #rotamers.setup_restraints(result)
    rc = rotamers.update_restraints(
      self.pdb_hierarchy,
      self.restraints_manager.geometry,
      current_geometry=self.pdb_hierarchy.extract_xray_structure(),
      rdl_proxies=self.rdl_proxies,
      log=self.log,
      verbose=False,
      )
    print("="*79, file=self.log)
    return rc

  def show(self):
    es = self.restraints_manager.geometry.energies_sites(
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
      compute_gradients = False)
    es.show(prefix="    ", f=self.log)


def minimize_wrapper_for_ramachandran(
    model,
    original_pdb_h,
    excl_string_selection,
    processed_pdb_file = None,
    log=None,
    reference_rotamers = True,
    number_of_cycles=1,
    run_first_minimization_without_reference=False,
    oldfield_weight_scale=3,
    oldfield_plot_cutoff=0.03,
    nonbonded_weight=500,
    reference_sigma=0.7):
  """ Wrapper around geometry minimization specifically tuned for eliminating
  Ramachandran outliers.
  probably not working anymore... no processed_pdb_file available.
  WARNING: no setting sites_cart at the end...
  """
  grm = model.get_restraints_manager()
  assert grm is not None
  from mmtbx.geometry_restraints import reference
  from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
      reference_model, reference_model_params
  from libtbx.utils import null_out
  from scitbx.array_family import flex
  if log is None:
    log = null_out()
  # assert hierarchy.atoms_size()==xrs.scatterers().size(), "%d %d" % (
  #     hierarchy.atoms_size(), xrs.scatterers().size())

  ncs_restraints_group_list = model.get_ncs_groups()
  if ncs_restraints_group_list is None:
    ncs_restraints_group_list = []

  grm.geometry.pair_proxies(
      sites_cart=model.get_sites_cart())
  grm.geometry.update_ramachandran_restraints_phi_psi_targets(
      hierarchy=model.get_hierarchy())

  if reference_rotamers and original_pdb_h is not None:
    # make selection excluding rotamer outliers
    from mmtbx.rotamer.rotamer_eval import RotamerEval
    # print "Excluding rotamer outliers"
    rotamer_manager = model.get_rotamer_manager()
    non_rot_outliers_selection = flex.bool(model.get_number_of_atoms(), False)
    for m in original_pdb_h.models():
      for chain in m.chains():
        for conf in chain.conformers():
          for res in conf.residues():
            ev = rotamer_manager.evaluate_residue_2(res)
            if ev != "OUTLIER" or ev is None:
              for a in res.atoms():
                non_rot_outliers_selection[a.i_seq] = True
            # else:
            #   print "  ", res.id_str()

    if processed_pdb_file is not None:
      rm_params = reference_model_params.extract()
      rm_params.reference_model.enabled=True
      rm_params.reference_model.strict_rotamer_matching=False
      rm_params.reference_model.main_chain=False
      rm = reference_model(
        processed_pdb_file=processed_pdb_file,
        reference_file_list=None,
        reference_hierarchy_list=[original_pdb_h],
        mon_lib_srv=model.get_mon_lib_srv(),
        ener_lib=model.get_ener_lib(),
        has_hd=None,
        params=rm_params.reference_model,
        selection=non_rot_outliers_selection,
        log=log)
      rm.show_reference_summary(log=log)
      grm.geometry.adopt_reference_dihedral_manager(rm)

  # dealing with SS
  if model.get_ss_annotation() is not None:
    from mmtbx.secondary_structure import manager
    ss_manager = manager(
        pdb_hierarchy=model.get_hierarchy(),
        geometry_restraints_manager=grm.geometry,
        sec_str_from_pdb_file=model.get_ss_annotation(),
        params=None,
        mon_lib_srv=model.get_mon_lib_srv(),
        verbose=-1,
        log=log)
    grm.geometry.set_secondary_structure_restraints(
        ss_manager=ss_manager,
        hierarchy=model.get_hierarchy(),
        log=log)

  if run_first_minimization_without_reference:
    obj = run2(
      restraints_manager=grm,
      pdb_hierarchy=model.get_hierarchy(),
      correct_special_position_tolerance=1.0,
      ncs_restraints_group_list=ncs_restraints_group_list,
      max_number_of_iterations=300,
      number_of_macro_cycles=number_of_cycles,
      bond=True,
      nonbonded=True,
      angle=True,
      dihedral=True,
      chirality=True,
      planarity=True,
      fix_rotamer_outliers=True,
      log=log)

  if original_pdb_h is not None:
    if not excl_string_selection or len(excl_string_selection) == 0:
      excl_string_selection = "all"
    asc = original_pdb_h.atom_selection_cache()
    sel = asc.selection("(%s) and (name CA or name C or name N or name O)" % excl_string_selection)

    grm.geometry.append_reference_coordinate_restraints_in_place(
        reference.add_coordinate_restraints(
            sites_cart = original_pdb_h.atoms().extract_xyz().select(sel),
            selection  = sel,
            sigma      = reference_sigma,
            top_out_potential=True))

  obj = run2(
      restraints_manager       = grm,
      pdb_hierarchy            = model.get_hierarchy(),
      correct_special_position_tolerance = 1.0,
      ncs_restraints_group_list=ncs_restraints_group_list,
      max_number_of_iterations = 300,
      number_of_macro_cycles   = number_of_cycles,
      bond                     = True,
      nonbonded                = True,
      angle                    = True,
      dihedral                 = True,
      chirality                = True,
      planarity                = True,
      fix_rotamer_outliers     = True,
      log                      = log)
  grm.geometry.reference_dihedral_manager=None
