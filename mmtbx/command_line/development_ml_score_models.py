from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME phenix.development.ml_score_models

import iotbx.pdb
from cctbx.array_family import flex
import sys
import mmtbx.f_model
from libtbx.utils import null_out
import mmtbx.utils
from iotbx import reflection_file_utils
import sys, random
from libtbx.test_utils import approx_equal
from libtbx import group_args
import mmtbx.model
from mmtbx.refinement import refinement_flags as reffl
from phenix.refinement import weight_xray_chem
from mmtbx.refinement import adp_refinement
from mmtbx.refinement import minimization
import scitbx.lbfgs
from cctbx import adptbx
from libtbx import adopt_init_args
from cctbx import xray
import os

random.seed(0)
flex.set_random_seed(0)

class lbfgs(object):

  def __init__(self, fmodels,
                     restraints_manager       = None,
                     model                    = None,
                     target_weights           = None,
                     refine_xyz               = False,
                     refine_adp               = False,
                     lbfgs_termination_params = None,
                     verbose                  = 0,
                     correct_special_position_tolerance = 1.0,
                     iso_restraints           = None,
                     u_min                    = adptbx.b_as_u(-5.0),
                     u_max                    = adptbx.b_as_u(999.99),
                     log                      = None):
    adopt_init_args(self, locals())
    self.f=None
    self.xray_structure = self.fmodels.fmodel_xray().xray_structure
    if(self.refine_adp):
      self.xray_structure.tidy_us()
      self.fmodels.update_xray_structure(
        xray_structure = self.xray_structure,
        update_f_calc  = True)
    self.weights = target_weights
    self.correct_special_position_tolerance = correct_special_position_tolerance
    if(refine_xyz and target_weights is not None):
      self.weights = target_weights.xyz_weights_result
    elif(refine_adp and target_weights is not None):
      self.weights = target_weights.adp_weights_result
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    self.compute_target(compute_gradients = False,u_iso_refinable_params = None)
    if(self.refine_adp):
      self.xray_structure.tidy_us()
      self.fmodels.update_xray_structure(
        xray_structure = self.xray_structure,
        update_f_calc  = True)

  def apply_shifts(self):
    if(self.refine_adp):
      xray.ext.truncate_shifts(
        shifts    = self.x,
        min_value = self.u_min,
        max_value = self.u_max)
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    if(self.refine_xyz):
      site_symmetry_table = self.xray_structure.site_symmetry_table()
      for i_seq in site_symmetry_table.special_position_indices():
        try:
          scatterers_shifted[i_seq].site = crystal.correct_special_position(
            crystal_symmetry = self.xray_structure,
            special_op       = site_symmetry_table.get(i_seq).special_op(),
            site_frac        = scatterers_shifted[i_seq].site,
            site_label       = scatterers_shifted[i_seq].label,
            tolerance        = self.correct_special_position_tolerance)
        except Exception, e:
          print >> self.log, str(e)
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)
    if(self.refine_adp):
      return apply_shifts_result.u_iso_refinable_params
    else:
      return None

  def compute_target(self, compute_gradients, u_iso_refinable_params):
    self.stereochemistry_residuals = None
    self.fmodels.update_xray_structure(
      xray_structure = self.xray_structure,
      update_f_calc  = True)
    fmodels_target_and_gradients = self.fmodels.target_and_gradients(
      weights                = self.weights,
      compute_gradients      = compute_gradients,
      u_iso_refinable_params = u_iso_refinable_params)
    self.f = fmodels_target_and_gradients.target()
    self.g = fmodels_target_and_gradients.gradients()
    if(self.refine_xyz and self.restraints_manager is not None and
       self.weights.w > 0.0):
      self.stereochemistry_residuals = \
        self.model.restraints_manager_energies_sites(
          compute_gradients = compute_gradients)
      er = self.stereochemistry_residuals.target
      self.f += er * self.weights.w
      if(compute_gradients):
        sgc = self.stereochemistry_residuals.gradients
        # ias do not participate in geometry restraints
        if(self.model is not None and self.model.ias_selection is not None and
           self.model.ias_selection.count(True) > 0):
          sgc.extend(flex.vec3_double(
            self.model.ias_selection.count(True),[0,0,0]))
        xray.minimization.add_gradients(
          scatterers     = self.xray_structure.scatterers(),
          xray_gradients = self.g,
          site_gradients = sgc*self.weights.w)
    if(self.refine_adp and self.restraints_manager is not None and
       self.restraints_manager.geometry is not None
       and self.weights.w > 0.0 and self.iso_restraints is not None):
      energies_adp = self.model.energies_adp(
        iso_restraints    = self.iso_restraints,
        use_hd            = False,
        compute_gradients = compute_gradients)
      self.f += energies_adp.target * self.weights.w
      if(compute_gradients):
        if(energies_adp.u_aniso_gradients is None):
          xray.minimization.add_gradients(
            scatterers      = self.xray_structure.scatterers(),
            xray_gradients  = self.g,
            u_iso_gradients = energies_adp.u_iso_gradients * self.weights.w)
        else:
          energies_adp.u_aniso_gradients *= self.weights.w
          if(energies_adp.u_iso_gradients is not None):
            energies_adp.u_iso_gradients *= self.weights.w
          xray.minimization.add_gradients(
            scatterers        = self.xray_structure.scatterers(),
            xray_gradients    = self.g,
            u_aniso_gradients = energies_adp.u_aniso_gradients,
            u_iso_gradients   = energies_adp.u_iso_gradients)
          energies_adp.u_aniso_gradients = None # just for safety
          energies_adp.u_iso_gradients = None

  def compute_functional_and_gradients(self):
    u_iso_refinable_params = self.apply_shifts()
    self.compute_target(compute_gradients     = True,
                        u_iso_refinable_params = u_iso_refinable_params)
    return self.f, self.g

def macro_cycle(m, fmodels, weights, alpha_beta, log):
  iso_restraints = adp_refinement.adp_restraints_master_params.extract().iso
  for it in xrange(3):
    # XYZ refinement
    m.model.set_refine_individual_sites()
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = 25)
    lbfgs(
      fmodels                  = fmodels,
      restraints_manager       = m.restraints_manager,
      model                    = m.model,
      target_weights           = weights,
      refine_xyz               = True,
      refine_adp               = False,
      lbfgs_termination_params = lbfgs_termination_params,
      correct_special_position_tolerance = 1.0,
      iso_restraints           = iso_restraints,
      log                      = log)
    rwm, rfm, twm, tfm = get_nums(fmodels.fmodel_xray(), alpha_beta)
    print "  XYZ: r_work, r_free: %6.4f %6.4f"%(rwm, rfm), twm, tfm
    # ADP refinement
    m.model.set_refine_individual_adp()
    lbfgs(
      fmodels                  = fmodels,
      restraints_manager       = m.restraints_manager,
      model                    = m.model,
      target_weights           = weights,
      refine_xyz               = False,
      refine_adp               = True,
      lbfgs_termination_params = lbfgs_termination_params,
      correct_special_position_tolerance = 1.0,
      iso_restraints           = iso_restraints,
      log                      = log)
    rwm, rfm, twm, tfm = get_nums(fmodels.fmodel_xray(), alpha_beta)
    print "  ADP: r_work, r_free: %6.4f %6.4f"%(rwm, rfm), twm, tfm

def get_model_objects(pdb_file_name, processed_args, log):
  # Remove H
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  ph = pdb_inp.construct_hierarchy()
  sel = ph.atom_selection_cache().selection("not (element H or element D)")
  ph = ph.select(sel)
  ph.write_pdb_file(file_name="tmp.pdb",
    crystal_symmetry=pdb_inp.crystal_symmetry())
  #
  mmtbx_pdb_file_1 = mmtbx.utils.pdb_file(
    pdb_file_names        = ["tmp.pdb"],
    cif_objects           = processed_args.cif_objects,
    crystal_symmetry      = processed_args.crystal_symmetry,
    use_neutron_distances = False,
    log                   = log)
  os.remove("tmp.pdb")
  mmtbx_pdb_file_1.set_ppf(stop_if_duplicate_labels = True)
  processed_pdb_file = mmtbx_pdb_file_1.processed_pdb_file
  ph = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5.0,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  rflags = flex.bool(xrs.scatterers().size(), True)
  refinement_flags = reffl.manager(
    individual_sites     = True,
    individual_adp       = True,
    sites_individual     = rflags,
    adp_individual_iso   = rflags,
    adp_individual_aniso = None)
  refinement_flags.check_all()
  model = mmtbx.model.manager(
    xray_structure     = xrs,
    pdb_hierarchy      = ph,
    refinement_flags   = refinement_flags,
    restraints_manager = restraints_manager)
  return group_args(
    pdb_hierachy       = ph,
    xray_structure     = xrs,
    restraints_manager = restraints_manager,
    model              = model)

def get_data(hkl_file, crystal_symmetry, log):
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = hkl_file,
    err              = log)
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    keep_going              = True,
    force_non_anomalous     = True,
    log                     = log)
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  return f_obs, r_free_flags

def get_nums(fmodel, alpha_beta):
  return fmodel.r_work(), fmodel.r_free(), \
    fmodel.target_unscaled_w(alpha_beta), \
    fmodel.target_unscaled_t(alpha_beta)

def run(args, log=null_out()):
  # All CL inputs
  processed_args = mmtbx.utils.process_command_line_args(args=args, log=log)
  # PDB files
  if(len(processed_args.pdb_file_names)!=2):
    raise Sorry("Two PDB files required.")
  pdb_file_1, pdb_file_2 = processed_args.pdb_file_names
  print "PDB files 1 (M1) and 2 (M2):", pdb_file_1, pdb_file_2
  m1 = get_model_objects(
    pdb_file_name  = pdb_file_1,
    processed_args = processed_args,
    log            = log)
  m2 = get_model_objects(
    pdb_file_name  = pdb_file_2,
    processed_args = processed_args,
    log            = log)
  mmtbx.utils.assert_xray_structures_equal(
    x1               = m1.xray_structure,
    x2               = m2.xray_structure,
    sites            = False,
    adp              = True,
    occupancies      = True,
    elements         = True,
    scattering_types = True,
    eps              = 1.e-6)
  # Data
  f_obs, r_free_flags = get_data(
    hkl_file         = processed_args.reflection_files,
    crystal_symmetry = processed_args.crystal_symmetry,
    log              = log,
    )
  # Fmodel
  print "Set up fmodel..."
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    xray_structure = m1.xray_structure,
    r_free_flags   = r_free_flags)
  #
  print "START: Update all scales..."
  fmodel.update_all_scales(update_f_part1=False)
  fmodels = mmtbx.fmodels(fmodel_xray = fmodel)
  alpha_beta = fmodel.alpha_beta()
  #
  #
  rwm1, rfm1, twm1, tfm1 = get_nums(fmodel, alpha_beta)
  print "  M1: r_work, r_free: %6.4f %6.4f"%(rwm1, rfm1), twm1, tfm1

  fmodels.update_xray_structure(xray_structure=m2.xray_structure,
    update_f_calc=True, update_f_mask=True)
  rwm2, rfm2, twm2, tfm2 = get_nums(fmodel, alpha_beta)
  print "  M2: r_work, r_free: %6.4f %6.4f"%(rwm2, rfm2), twm2, tfm2

  fmodels.update_xray_structure(xray_structure=m1.xray_structure,
    update_f_calc=True, update_f_mask=True)
  rwm1_, rfm1_, twm1_, tfm1_ = get_nums(fmodel, alpha_beta)

  fmodels.update_xray_structure(xray_structure=m2.xray_structure,
    update_f_calc=True, update_f_mask=True)
  rwm2_, rfm2_, twm2_, tfm2_ = get_nums(fmodel, alpha_beta)

  assert approx_equal(rwm1, rwm1_)
  assert approx_equal(rwm2, rwm2_)
  assert approx_equal(rfm1, rfm1_)
  assert approx_equal(rfm2, rfm2_)

  assert approx_equal(twm1, twm1_)
  assert approx_equal(twm2, twm2_)
  assert approx_equal(tfm1, tfm1_)
  assert approx_equal(tfm2, tfm2_)
  # Compute weights
  log.flush()
  print "Compute weights..."
  iso_restraints = adp_refinement.adp_restraints_master_params.extract().iso
  weights = weight_xray_chem.weight(
    fmodels                            = fmodels,
    model                              = m1.model,
    correct_special_position_tolerance = 2.,
    target_weights_params              = weight_xray_chem.master_params.extract(),
    iso_restraints                     = iso_restraints,
    macro_cycle                        = 0,
    show_summary                        = False)
  # Compute alpha and beta and target functors
  fmodels.create_target_functors(alpha_beta = alpha_beta)
  # Refine M1 and M2
  print "Refine M1, then M2 (independently)..."
  for i, m in enumerate([m1,m2]):
    i+=1
    fmodels.update_xray_structure(xray_structure=m.xray_structure,
      update_f_calc=True, update_f_mask=True)
    rwm, rfm, twm, tfm = get_nums(fmodel, alpha_beta)
    print "Refine M%s, start r_work, r_free: %6.4f %6.4f"%(
      i, rwm, rfm), twm, tfm
    macro_cycle(m=m, fmodels=fmodels, weights=weights, alpha_beta=alpha_beta,
      log=log)


if(__name__ == "__main__"):
  run(sys.argv[1:])





                                                                                                        
