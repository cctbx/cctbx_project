from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import iotbx.phil
from cctbx import geometry_restraints
import cctbx.geometry_restraints.lbfgs
import mmtbx.refinement.minimization_ncs_constraints
import scitbx.lbfgs
import sys
import mmtbx.utils
from scitbx.array_family import flex
from mmtbx import monomer_library

master_params_str = """\
  alternate_nonbonded_off_on=False
    .type = bool
  restrain_c_alpha_positions=False
    .type = bool
  max_iterations=500
    .type = int
  macro_cycles=1
    .type = int
  show_geometry_restraints=False
    .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

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

def run(processed_pdb_file, params=master_params().extract(), log=sys.stdout):
  co = params
  geometry_restraints_flags = geometry_restraints.flags.flags(default=True)
  all_chain_proxies = processed_pdb_file.all_chain_proxies
  reference_manager = None
  if (co.restrain_c_alpha_positions):
    assert 0, "Broken! - rewrite this part first."
    from mmtbx.geometry_restraints import reference
    ca_selection=all_chain_proxies.pdb_hierarchy.get_peptide_c_alpha_selection()
    ca_sites_cart = \
      all_chain_proxies.sites_cart.deep_copy().select(ca_selection)
    reference_manager = reference.manager()
    reference_manager.add_coordinate_restraints(
      sites_cart=ca_sites_cart,
      selection=ca_selection)
  geometry_restraints_manager = processed_pdb_file.\
    geometry_restraints_manager(show_energies = False,
                                reference_manager=\
                                  reference_manager)
  special_position_settings = all_chain_proxies.special_position_settings
  sites_cart = all_chain_proxies.sites_cart_exact().deep_copy()
  atom_labels = [atom.id_str() for atom in all_chain_proxies.pdb_atoms]
  geometry_restraints_manager.site_symmetry_table \
    .show_special_position_shifts(
      special_position_settings=special_position_settings,
      site_labels=atom_labels,
      sites_cart_original=all_chain_proxies.sites_cart,
      sites_cart_exact=sites_cart,
      out=log,
      prefix="  ")
  if (co.show_geometry_restraints):
    geometry_restraints_manager.show_sorted(
      flags=geometry_restraints_flags,
      sites_cart=sites_cart,
      site_labels=atom_labels)
  pair_proxies =  geometry_restraints_manager.pair_proxies(
    sites_cart=all_chain_proxies.sites_cart,
    flags=geometry_restraints_flags)
  pair_proxies.bond_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    site_labels=atom_labels,
    f=log,
    max_items=10)
  if (pair_proxies.nonbonded_proxies is not None):
    pair_proxies.nonbonded_proxies.show_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      site_labels=atom_labels,
      f=log,
      max_items=10)
  del pair_proxies
  print(file=log)
  log.flush()
  if (co.alternate_nonbonded_off_on and co.macro_cycles % 2 != 0):
    co.macro_cycles += 1
    print("INFO: Number of macro cycles increased by one to ensure use of", file=log)
    print("      nonbonded interactions in last macro cycle.", file=log)
    print(file=log)
  for i_macro_cycle in range(co.macro_cycles):
    if (co.alternate_nonbonded_off_on):
      geometry_restraints_flags.nonbonded = bool(i_macro_cycle % 2)
      print("Use nonbonded interactions this macro cycle:", \
        geometry_restraints_flags.nonbonded, file=log)
    minimized = lbfgs(
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=co.max_iterations),
      lbfgs_exception_handling_params=
        scitbx.lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound=True))
    print("Energies at start of minimization:", file=log)
    minimized.first_target_result.show(f=log)
    print(file=log)
    print("Number of minimization iterations:", minimized.minimizer.iter(), file=log)
    print("Root-mean-square coordinate difference: %.3f" % (
      all_chain_proxies.sites_cart.rms_difference(sites_cart)), file=log)
    print(file=log)
    print("Energies at end of minimization:", file=log)
    minimized.final_target_result.show(f=log)
    print(file=log)
    geometry_restraints_manager.pair_proxies(
      sites_cart=sites_cart,
      flags=geometry_restraints_flags) \
        .bond_proxies.show_sorted(
          by_value="residual",
          sites_cart=sites_cart,
          site_labels=atom_labels,
          f=log,
          max_items=10)
    print(file=log)
  assert geometry_restraints_flags.nonbonded
  return sites_cart

def add_rotamer_restraints(
      pdb_hierarchy,
      restraints_manager,
      selection,
      sigma,
      mode,
      accept_allowed=True,
      mon_lib_srv=None,
      rotamer_manager=None):
  pdb_hierarchy_for_proxies = mmtbx.utils.switch_rotamers(
    pdb_hierarchy  = pdb_hierarchy.deep_copy(),
    mode           = "exact_match",
    accept_allowed = accept_allowed,
    selection      = selection,
    mon_lib_srv    = mon_lib_srv,
    rotamer_manager= rotamer_manager)
  mmtbx.utils.switch_rotamers(
    pdb_hierarchy  = pdb_hierarchy,
    mode           = mode,
    accept_allowed = accept_allowed,
    selection      = selection,
    mon_lib_srv    = mon_lib_srv,
    rotamer_manager= rotamer_manager)
  restraints_manager.geometry.add_chi_torsion_restraints_in_place(
      pdb_hierarchy   = pdb_hierarchy_for_proxies,
      sites_cart      = pdb_hierarchy_for_proxies.atoms().extract_xyz(),
      chi_angles_only = True,
      sigma           = sigma)
  return pdb_hierarchy, restraints_manager

class run2(object):
  def __init__(self,
               restraints_manager,
               pdb_hierarchy,
               correct_special_position_tolerance,
               riding_h_manager               = None,
               ncs_restraints_group_list      = [],
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
               mon_lib_srv                    = None
               ):
    self.log = log
    if self.log is None:
      self.log = sys.stdout
    self.pdb_hierarchy = pdb_hierarchy
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
        sites_cart=self.pdb_hierarchy.atoms().extract_xyz())
      if(alternate_nonbonded_off_on and i_macro_cycle<=number_of_macro_cycles/2):
        geometry_restraints_flags.nonbonded = bool(i_macro_cycle % 2)
      self.update_cdl_restraints(macro_cycle=i_macro_cycle)
      if(fix_rotamer_outliers):
        self.pdb_hierarchy, self.restraints_manager = add_rotamer_restraints(
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
                ncs_restraints_group_list=ncs_restraints_group_list,
                refine_selection=refine_selection,
                restraints_manager=self.restraints_manager.geometry,
                refine_sites=True,
                refine_transformations=False,
                )
        minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
          target_and_grads_object      = tfg_obj,
          xray_structure               = xrs,
          ncs_restraints_group_list    = ncs_restraints_group_list,
          refine_selection             = refine_selection,
          finite_grad_differences_test = False,
          max_iterations               = max_number_of_iterations,
          refine_sites                 = True,
          refine_transformations       = False)
        self.pdb_hierarchy.adopt_xray_structure(xrs)
      else:
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
        self.pdb_hierarchy.atoms().set_xyz(sites_cart)
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
    hierarchy,
    xrs,
    original_pdb_h,
    excl_string_selection,
    grm,
    log=None,
    ncs_restraints_group_list=[],
    ss_annotation = None,
    mon_lib_srv=None,
    ener_lib=None,
    rotamer_manager=None,
    reference_rotamers = True,
    number_of_cycles=1,
    run_first_minimization_without_reference=False,
    oldfield_weight_scale=3,
    oldfield_plot_cutoff=0.03,
    nonbonded_weight=500,
    reference_sigma=0.7):
  """ Wrapper around geometry minimization specifically tuned for eliminating
  Ramachandran outliers.
  """
  try:
    import pickle as pickle
  except ImportError:
    import pickle
  assert grm is not None
  from time import time
  from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
  from mmtbx.geometry_restraints import reference
  from mmtbx.command_line.geometry_minimization import \
      get_geometry_restraints_manager
  from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
      reference_model, reference_model_params
  from libtbx.utils import null_out
  from scitbx.array_family import flex
  if log is None:
    log = null_out()
  # assert hierarchy.atoms_size()==xrs.scatterers().size(), "%d %d" % (
  #     hierarchy.atoms_size(), xrs.scatterers().size())

  grm.geometry.pair_proxies(
      sites_cart=hierarchy.atoms().extract_xyz())
  if grm.geometry.ramachandran_manager is not None:
    grm.geometry.ramachandran_manager.update_phi_psi_targets(
        sites_cart=hierarchy.atoms().extract_xyz())

  if reference_rotamers and original_pdb_h is not None:
    # make selection excluding rotamer outliers
    from mmtbx.rotamer.rotamer_eval import RotamerEval
    # print "Excluding rotamer outliers"
    if rotamer_manager is None:
      rotamer_manager = RotamerEval(mon_lib_srv=mon_lib_srv)
    non_rot_outliers_selection = flex.bool(hierarchy.atoms_size(), False)
    for model in original_pdb_h.models():
      for chain in model.chains():
        for conf in chain.conformers():
          for res in conf.residues():
            ev = rotamer_manager.evaluate_residue_2(res)
            if ev != "OUTLIER" or ev is None:
              for a in res.atoms():
                non_rot_outliers_selection[a.i_seq] = True
            # else:
            #   print "  ", res.id_str()


    rm_params = reference_model_params.extract()
    rm_params.reference_model.enabled=True
    rm_params.reference_model.strict_rotamer_matching=False
    rm_params.reference_model.main_chain=False
    rm = reference_model(
      processed_pdb_file=processed_pdb_file,
      reference_file_list=None,
      reference_hierarchy_list=[original_pdb_h],
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      has_hd=None,
      params=rm_params.reference_model,
      selection=non_rot_outliers_selection,
      log=log)
    rm.show_reference_summary(log=log)
    grm.geometry.adopt_reference_dihedral_manager(rm)

  # dealing with SS
  if ss_annotation is not None:
    from mmtbx.secondary_structure import manager
    ss_manager = manager(
        pdb_hierarchy=hierarchy,
        geometry_restraints_manager=grm.geometry,
        sec_str_from_pdb_file=ss_annotation,
        params=None,
        mon_lib_srv=mon_lib_srv,
        verbose=-1,
        log=log)
    grm.geometry.set_secondary_structure_restraints(
        ss_manager=ss_manager,
        hierarchy=hierarchy,
        log=log)

  # grm pickle-unpickle
  # t0 = time()
  # prefix="grm"
  # pklfile = open("%s.pkl" % prefix, 'wb')
  # pickle.dump(grm.geometry, pklfile)
  # pklfile.close()
  # t1 = time()
  # pklfile = open("%s.pkl" % prefix, 'rb')
  # grm_from_file = pickle.load(pklfile)
  # pklfile.close()
  # t2 = time()
  # print "Time pickling/unpickling: %.4f, %.4f" % (t1-t0, t2-t1)
  # grm.geometry=grm_from_file


  if run_first_minimization_without_reference:
    obj = run2(
      restraints_manager=grm,
      pdb_hierarchy=hierarchy,
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
  # grm.geometry.write_geo_file(
  #     sites_cart=hierarchy.atoms().extract_xyz(),
  #     site_labels=[atom.id_str() for atom in hierarchy.atoms()],
  #     file_name="last_gm.geo")
  obj = run2(
      restraints_manager       = grm,
      pdb_hierarchy            = hierarchy,
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
