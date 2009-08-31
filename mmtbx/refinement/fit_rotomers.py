from cctbx.array_family import flex
from libtbx import adopt_init_args
from mmtbx.command_line import lockit
from mmtbx import real_space_correlation
from libtbx import adopt_init_args
from cctbx import maptbx
import scitbx.lbfgs
import iotbx.pdb
import mmtbx.monomer_library
import mmtbx.model
import mmtbx.refinement.real_space
from cctbx import miller

class select_map(object):
  def __init__(self, unit_cell, map_data_1, map_data_2):
    adopt_init_args(self, locals())
    assert map_data_1.focus() == map_data_2.focus()
    assert map_data_1.all() == map_data_2.all()
    self.fft_n_real = map_data_1.focus()
    self.fft_m_real = map_data_1.all()

  def select(self, sites_cart, atom_radius=2.0):
    return maptbx.grid_indices_around_sites(
      unit_cell  = self.unit_cell,
      fft_n_real = self.fft_n_real,
      fft_m_real = self.fft_m_real,
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), atom_radius))

  def get_cc(self, sites_cart, residue_iselection):
    sel_map = self.select(sites_cart = sites_cart)
    m1 = self.map_data_1.select(sel_map)
    m2 = self.map_data_2.select(sel_map)
    return flex.linear_correlation(x = m1, y = m2).coefficient()

class refiner(object):
  def __init__(self,
               pdb_hierarchy,
               target_map,
               geometry_restraints_manager,
               real_space_target_weight,
               real_space_gradients_delta,
               max_iterations,
               min_iterations):
    adopt_init_args(self, locals())
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations, min_iterations = min_iterations)

  def refine_constrained(self, residue):
    return lockit.residue_refine_constrained(
      pdb_hierarchy               = self.pdb_hierarchy,
      residue                     = residue,
      density_map                 = self.target_map,
      geometry_restraints_manager = self.geometry_restraints_manager,
      real_space_target_weight    = self.real_space_target_weight,
      real_space_gradients_delta  = self.real_space_gradients_delta,
      lbfgs_termination_params    = self.lbfgs_termination_params)

  def refine_restrained(self, residue):
    return lockit.residue_refine_restrained(
      pdb_hierarchy               = self.pdb_hierarchy,
      residue                     = residue,
      density_map                 = self.target_map,
      geometry_restraints_manager = self.geometry_restraints_manager,
      real_space_target_weight    = self.real_space_target_weight,
      real_space_gradients_delta  = self.real_space_gradients_delta,
      lbfgs_termination_params    = self.lbfgs_termination_params)

def iterate_rotamers(pdb_hierarchy,
                     xray_structure,
                     map_data_1,
                     map_data_2,
                     map_data_3,
                     mon_lib_srv,
                     rsr_manager,
                     f_obs, fft_map,
                     poor_cc_threshold,
                     log):
  assert map_data_1.focus() == map_data_2.focus()
  assert map_data_1.all() == map_data_2.all()
  fmt="  %s %s %s: cc_start=%6.4f score: start=%6.2f final=%6.2f diff=%6.2f best_rotomer_id: %5s"
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    map_data_1 = map_data_1,
    map_data_2 = map_data_2)
  get_class = iotbx.pdb.common_residue_names_get_class
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  sites_cart_start = xray_structure.sites_cart()
  #
  def target(sites_cart_residue):
    sites_frac_residue = xray_structure.unit_cell().fractionalize(
      sites_cart_residue)
    result = 0
    for rsf in sites_frac_residue:
      result += map_data_1.eight_point_interpolation(rsf)
    return result
  def target_diff(sites_cart_residue):
    sites_frac_residue = xray_structure.unit_cell().fractionalize(
      sites_cart_residue)
    result = 0
    for rsf in sites_frac_residue:
      result += map_data_3.eight_point_interpolation(rsf)
    return result
  #
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        residue_iselection = flex.size_t()
        for atom in residue.atoms():
          residue_iselection.append(atom.i_seq)
        sites_cart = xray_structure.select(residue_iselection).sites_cart()
        cc_start = map_selector.get_cc(
          sites_cart         = sites_cart,
          residue_iselection = residue_iselection)
        cc = cc_start
        t_start = target(sites_cart)
        t3_start = target_diff(sites_cart)+t_start
        if(cc_start < poor_cc_threshold):
          #t_best = 0
          t_best = 0#t_start
          t3_best = t3_start
          if(get_class(residue.resname) != "common_amino_acid"):
            n_other_residues += 1
          else:
            rotamer_iterator = lockit.get_rotamer_iterator(
              mon_lib_srv         = mon_lib_srv,
              residue             = residue,
              atom_selection_bool = None)
            if(rotamer_iterator is None):
              n_amino_acids_ignored += 1
            else:
              n_amino_acids_scored += 1
              t_start = rsr_manager.refine_restrained(residue = residue).rf_only
              t_best = t_start
              rotamer_id_best = None
              residue_sites_best = None
              for rotamer, rotamer_sites_cart in rotamer_iterator:
                #XXX
                if 1:
                  T = target(rotamer_sites_cart)
                  T3_ = target_diff(rotamer_sites_cart)
                  T3 = T3_+T
                  if(T3 > t3_best and T3_ > 0):
                  #if(T > t_best):
                    #t_best = T
                    t3_best = T3
                    rotamer_id_best = rotamer.id
                    residue_sites_best = rotamer_sites_cart.deep_copy()
                    residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
                    sites_cart_start = sites_cart_start.set_selected(
                      residue_iselection, rotamer_sites_cart)
                #XXX
                else:
                  residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
                  refined = rsr_manager.refine_restrained(residue = residue)
                  rotamer_sites_cart_refined = \
                    refined.sites_cart_residue.deep_copy()
                  if(refined.rf_only > t_best):
                    t_best = refined.rs_f_final
                    rotamer_id_best = rotamer.id
                    residue_sites_best = rotamer_sites_cart_refined.deep_copy()
                    residue.atoms().set_xyz(new_xyz=rotamer_sites_cart_refined)
                    sites_cart_start = sites_cart_start.set_selected(
                      residue_iselection, rotamer_sites_cart_refined)
              if(residue_sites_best is not None):
                residue.atoms().set_xyz(new_xyz=residue_sites_best)
                ##
                refined = rsr_manager.refine_restrained(residue = residue)
                T = target(refined.sites_cart_residue)
                T3_ = target_diff(refined.sites_cart_residue)
                T3 = T3_+T
                if(T3 > t3_best):
                  #print t3_best, T3
                  # XXX use CC to accept the result!!! ...hm, this doesn't work...
                  #XXXxrs_dc = xray_structure.deep_copy_scatterers()
                  #XXXsc_dc = sites_cart_start.deep_copy()
                  #XXXsc_dc.set_selected(residue_iselection, refined.sites_cart_residue)
                  #XXXxrs_dc.set_sites_cart(sc_dc)
                  #XXXf_calc = f_obs.structure_factors_from_scatterers(
                  #XXX  xray_structure = xrs_dc).f_calc()
                  #XXXfft_map_r = miller.fft_map(
                  #XXX  crystal_gridding     = fft_map,
                  #XXX  fourier_coefficients = f_calc)
                  #XXXfft_map_r.apply_sigma_scaling()
                  #XXXmdr = fft_map_r.real_map_unpadded()
                  #XXXms = select_map(
                  #XXX  unit_cell  = xray_structure.unit_cell(),
                  #XXX  map_data_1 = map_data_1,
                  #XXX  map_data_2 = mdr)
                  #XXXcc = map_selector.get_cc(
                  #XXX  sites_cart         = refined.sites_cart_residue,
                  #XXX  residue_iselection = residue_iselection)
                  # XXX
                  if 1:#(cc > cc_start):
                    rotamer_sites_cart_refined = \
                        refined.sites_cart_residue.deep_copy()
                    sites_cart_start = sites_cart_start.set_selected(
                          residue_iselection, rotamer_sites_cart_refined)
                    residue.atoms().set_xyz(new_xyz=rotamer_sites_cart_refined)
                ##
          if(t3_best is not None and t3_best-t3_start > 0):
            print >> log, fmt%(chain.id,residue.resname,residue.resseq,cc_start,
              t3_start, t3_best, t3_best-t3_start, rotamer_id_best), \
              "%6.4f %6.4f"%(cc_start, cc)
          #if(t_best is not None and t_best-t_start > 0):
          #  print >> log, fmt%(chain.id,residue.resname,residue.resseq,cc_start,
          #    t_start, t_best, t_best-t_start, rotamer_id_best)
  xray_structure.set_sites_cart(sites_cart_start)

def run(fmodel,
        geometry_restraints_manager,
        pdb_hierarchy,
        number_of_macro_cycles,
        do_not_use_dihedrals,
        solvent_selection,
        poor_cc_threshold,
        log):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  def get_map_data(map_type, resolution_factor=1./4):
    map_obj = fmodel.electron_density_map()
    fft_map = map_obj.fft_map(resolution_factor = resolution_factor,
      map_type = map_type)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded(),fft_map
  if(do_not_use_dihedrals):
    sel = flex.bool(fmodel.xray_structure.scatterers().size(), True)
    geometry_restraints_manager = geometry_restraints_manager.select(sel)
    geometry_restraints_manager.remove_dihedrals_in_place(sel)
    #XXX geometry_restraints_manager = geometry_restraints_manager.select(
    #XXX   ~solvent_selection)
  # XXX
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints_manager,
    normalization = True)
  model = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = fmodel.xray_structure,
    pdb_hierarchy = pdb_hierarchy)
  # XXX
  fmt = "Macro-cycle %2d: r_work=%6.4f r_free=%6.4f"
  print >> log, fmt%(0, fmodel.r_work(), fmodel.r_free())
  for macro_cycle in range(1,number_of_macro_cycles+1):
    map_data_1,fft_map_1 = get_map_data(map_type = "2mFo-DFc")
    map_data_2,fft_map_2 = get_map_data(map_type = "Fc")
    map_data_3,fft_map_3 = get_map_data(map_type = "mFo-DFc")
    rsr_manager = refiner(
      pdb_hierarchy               = pdb_hierarchy,
      target_map                  = map_data_1,
      geometry_restraints_manager = geometry_restraints_manager,
      real_space_target_weight    = 100,
      real_space_gradients_delta  = fmodel.f_obs.d_min()/4,
      max_iterations              = 55,
      min_iterations              = 50)
    iterate_rotamers(
      pdb_hierarchy     = pdb_hierarchy,
      xray_structure    = fmodel.xray_structure,
      map_data_1        = map_data_1,
      map_data_2        = map_data_2,
      map_data_3        = map_data_3,
      mon_lib_srv       = mon_lib_srv,
      rsr_manager       = rsr_manager,
      f_obs=fmodel.f_obs,fft_map = fft_map_1,
      poor_cc_threshold = poor_cc_threshold,
      log               = log)
    fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
    print >> log, "1:", fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    if 1:
      assert model.xray_structure is fmodel.xray_structure
      params = mmtbx.refinement.real_space.master_params().extract()
      params.real_space_refinement.mode="diff_map"
      mmtbx.refinement.real_space.run(
        fmodel = fmodel, # XXX neutron ?
        model  = model,
        params = params.real_space_refinement,
        log    = log)
      assert model.xray_structure is fmodel.xray_structure
      fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
      print >> log, "2:",fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
