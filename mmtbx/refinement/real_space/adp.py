from __future__ import absolute_import, division, print_function
import mmtbx.refinement.real_space.utils
import mmtbx.refinement.utils
from scitbx.array_family import flex
from cctbx import adptbx
from libtbx import easy_mp
from mmtbx import bulk_solvent
from libtbx.test_utils import approx_equal
from six.moves import range
from cctbx import crystal
from mmtbx.refinement import adp_refinement
from cctbx import adp_restraints
from libtbx import group_args
from libtbx.str_utils import format_value
import scitbx.minimizers

import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

import random
random.seed(0)
flex.set_random_seed(0)

def map_and_model_to_fmodel(map_data, xray_structure, atom_radius, d_min,
                            reset_adp=True):
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure = xray_structure,
    map_data       = map_data,
    box_cushion    = atom_radius)
  box.apply_mask_inplace(atom_radius = atom_radius)
  f_obs_complex = box.box_map_coefficients(d_min = d_min)
  f_obs = abs(f_obs_complex)
  if(flex.mean(f_obs.data())<1.e-6): return None
  xrs = box.xray_structure_box.deep_copy_scatterers()
  if(reset_adp):
    vals_init = xrs.extract_u_iso_or_u_equiv()
    xrs = xrs.set_b_iso(value=0)
    assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    sc = flex.sum(abs(f_obs).data()*abs(f_calc).data())/ \
         flex.sum(abs(f_calc).data()*abs(f_calc).data())
    f_calc = f_calc.array(data = f_calc.data()*sc)
    o = bulk_solvent.complex_f_kb_scaled(
      f1      = f_obs_complex.data(),
      f2      = f_calc.data(),
      b_range = flex.double(range(5,505,5)),
      ss      = 1./flex.pow2(f_calc.d_spacings().data()) / 4.)
    xrs = xrs.set_b_iso(value=o.b())
    k_isotropic = flex.double(f_calc.data().size(), o.k())
    if(o.k()<1.e-6):
      k_isotropic = flex.double(f_calc.data().size(), 1)
      xrs.set_u_iso(values = vals_init)
  fmodel = mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs)
  if(reset_adp):
    fmodel.update_core(k_isotropic = k_isotropic)
  fmodel.update(target_name="ls_wunit_k1")
  fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
    remove_outliers=False)
  return fmodel

def get_plain_pair_sym_table(crystal_symmetry, sites_frac, plain_pairs_radius=5):
  asu_mappings = crystal.symmetry.asu_mappings(crystal_symmetry,
    buffer_thickness = plain_pairs_radius)
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry = crystal_symmetry)
  sites_cart = crystal_symmetry.unit_cell().orthogonalize(sites_frac)
  site_symmetry_table = special_position_settings.site_symmetry_table(
    sites_cart = sites_cart)
  asu_mappings.process_sites_frac(
    original_sites      = sites_frac,
    site_symmetry_table = site_symmetry_table)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff = plain_pairs_radius)
  return pair_asu_table.extract_pair_sym_table()

class tg(object):
  def __init__(self, fmodel, x, restraints_weight, bound_flags,
               lower_bound, upper_bound):
    self.bound_flags = bound_flags
    self.lower_bound = lower_bound
    self.upper_bound = upper_bound
    self.restraints_weight = restraints_weight
    self.fmodel = fmodel
    self.plain_pair_sym_table = None
    self.adp_iso_params = None
    if(self.restraints_weight is not None):
      self.plain_pair_sym_table = get_plain_pair_sym_table(
        crystal_symmetry = self.fmodel.xray_structure.crystal_symmetry(),
        sites_frac       = self.fmodel.xray_structure.sites_frac())
      self.adp_iso_params = \
        adp_refinement.adp_restraints_master_params.extract().iso
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
      iselection = self.fmodel.xray_structure.all_selection().iselection())
    # required fields
    self.x = x
    self.t = None
    self.g = None
    self.d = None
    self.use_curvatures=False
    #
    self.weight = None
    if(self.restraints_weight is not None):
      self.weight = self._weight()
    self.tgo = self._compute(x = self.x)
    self.update_target_and_grads(x=x)

  def _weight(self):
    num = self._restraints().gradients.norm()
    den = self._data().gradient_xray.norm()
    if(den==0): return 1
    return num/den

  def _restraints(self):
    return adp_restraints.energies_iso(
      plain_pair_sym_table = self.plain_pair_sym_table,
      xray_structure       = self.fmodel.xray_structure,
      parameters           = self.adp_iso_params,
      compute_gradients    = True,
      use_u_local_only     = self.adp_iso_params.use_u_local_only,
      use_hd               = False)

  def _data(self):
    fmodels = mmtbx.fmodels(fmodel_xray = self.fmodel)
    return fmodels.target_and_gradients(compute_gradients=True)

  def _compute(self, x):
    self.fmodel.xray_structure.set_b_iso(values = x)
    self.fmodel.update_xray_structure(update_f_calc = True)
    D = self._data()
    if(self.restraints_weight is not None):
      R = self._restraints()
      self.tgo = group_args(
        target   = D.target()*self.weight + R.target*self.restraints_weight,
        gradient = D.gradient_xray*self.weight + R.gradients*self.restraints_weight)
    else:
      self.tgo = group_args(
        target   = D.target(),
        gradient = D.gradient_xray)
    return self.tgo

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    self.tgo = self._compute(x=self.x)
    self.t = self.tgo.target
    self.g = self.tgo.gradient

  def target(self): return self.t

  def gradients(self): return self.g

  def gradient(self): return self.gradients()

class ncs_aware_refinement(object):
  def __init__(self, map_model_manager, d_min, atom_radius, nproc=1,
               log = None, individual = True, restraints_weight = 1,
               group_mode = None):
    self.mmm               = map_model_manager
    self.nproc             = nproc
    self.d_min             = d_min
    self.atom_radius       = atom_radius
    self.log               = log
    self.individual        = individual
    self.restraints_weight = restraints_weight
    self.group_mode        = group_mode
    #
    ncs_groups = self.mmm.model().get_ncs_groups()
    if(ncs_groups is None or len(ncs_groups)==0):
      values = self.run_one()
      self.mmm.model().set_b_iso(values = values)
    else:
      values = self.mmm.model().get_b_iso()
      for i, g in enumerate(ncs_groups):
        values_g = self.run_one(selection = g.master_iselection)
        values = values.set_selected(g.master_iselection, values_g)
        for j, c in enumerate(g.copies):
          values = values.set_selected(c.iselection, values_g)
      self.mmm.model().set_b_iso(values = values)

  def run_one(self, selection=None):
    model = self.mmm.model()
    if(selection is not None): model = model.select(selection)
    values = model.get_b_iso()
    model.get_hierarchy().atoms().reset_i_seq()
    if(self.nproc==1):
      args = [model,]
      return self.run_one_one(args = args)
    else:
      argss = []
      selections = []
      for sel in model.macromolecule_plus_hetatms_by_chain_selections():
        model_i = model.select(sel)
        if(model_i.size()==1):
          chain_ids = " ".join([c.id for c in model_i.get_hierarchy().chains()])
          print("Skip one atom model, chains: (%s)"%chain_ids, file=self.log)
          continue
        argss.append([model_i,])
        selections.append(sel) # XXX CAN BE BIG
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.run_one_one,
        args         = argss,
        func_wrapper = "buffer_stdout_stderr")
      for i, result in enumerate(stdout_and_results):
        values = values.set_selected(selections[i], result[1])
      model.set_b_iso(values = values)
      return values

  def run_one_one(self, args):
    model = args[0]
    log = self.log
    if(self.nproc>1): log = None
    #
    fmodel = map_and_model_to_fmodel(
      map_data       = self.mmm.map_data().deep_copy(),
      xray_structure = model.get_xray_structure(),
      atom_radius    = self.atom_radius,
      d_min          = self.d_min)
    if(fmodel is None):
      return model.get_xray_structure().extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    # selections for group ADP
    ph_box = model.get_hierarchy()
    ph_box.atoms().reset_i_seq()
    group_adp_sel = []
    if(self.group_mode=="one_per_residue"):
      for rg in ph_box.residue_groups():
        group_adp_sel.append(rg.atoms().extract_i_seq())
    elif(self.group_mode=="one_per_chain"):
      for chain in ph_box.chains():
        group_adp_sel.append(chain.atoms().extract_i_seq())
    elif(self.group_mode=="two_per_residue"):
      group_adp_sel = \
        mmtbx.refinement.utils.two_group_selections_per_residue(model)
    else:
      assert self.group_mode is None
    if self.group_mode is not None:
      #
      b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
      if(flex.max(b_isos)<1.e-2):
        b_isos = flex.random_double(model.size())*10
        model.set_b_iso(values = b_isos)
        fmodel.xray_structure.set_b_iso(values = b_isos)
        fmodel.update_xray_structure(xray_structure = fmodel.xray_structure,
          update_f_calc = True)
      #
      number_of_macro_cycles = 5
      if(self.individual): number_of_macro_cycles = 1
      group_b_manager = mmtbx.refinement.group.manager(
        fmodel                   = fmodel,
        selections               = group_adp_sel,
        convergence_test         = False,
        max_number_of_iterations = 50,
        number_of_macro_cycles   = number_of_macro_cycles,
        run_finite_differences_test = False,
        use_restraints           = True,
        refine_adp               = True,
        log                      = log)
      fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
        remove_outliers=False)
      b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
      model.set_b_iso(values = b_isos)
    #
    if(self.individual):
      if(log is not None):
        rmsb = format_value("%4.2f", model.rms_b_iso_or_b_equiv()).strip()
        print("r_work (start): %6.4f rms_B_bonded: %s"%(fmodel.r_work(),
          rmsb), file=log)
      rw = self.restraints_weight
      if(rw is not None):
        flipped = False
        for it in range(1,20):
          x = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          lower = flex.double(x.size(), 0)
          upper = flex.double(x.size(), flex.max(x)*2)
          calculator = tg(
            fmodel            = fmodel,
            x                 = x,
            restraints_weight = rw,
            bound_flags       = flex.int(x.size(), 2),
            lower_bound       = lower,
            upper_bound       = upper)
          rw_prev = rw
          b_isos_prev = b_isos
          rms_b_prev = model.rms_b_iso_or_b_equiv()
          m = scitbx.minimizers.lbfgs(
            mode='lbfgsb', max_iterations=100, calculator=calculator)
          b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          model.set_b_iso(values = b_isos)
          if(rms_b_prev is not None):
            rms_b = model.rms_b_iso_or_b_equiv()
            if(rms_b<5):
              rw = rw/2
              if(flipped):
                b_isos = b_isos_prev
                model.set_b_iso(values = b_isos)
                break
            else:
              if(rms_b > rms_b_prev):
                b_isos = b_isos_prev
                model.set_b_iso(values = b_isos)
                break
              rw = rw*2
              flipped = True
            if(log is not None):
              print(
                "r_work: %6.4f rms_B_bonded: %4.2f restraints_weight: %6.4f"%(
                fmodel.r_work(), rms_b, rw), file=log)
          if(rms_b_prev is None): break
      else:
        for it in range(1,2):
          x = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          ma = min(999,flex.max(x)*2)
          if ma<1: ma = 100
          lower = flex.double(x.size(), 0)
          upper = flex.double(x.size(), ma )
          calculator = tg(
            fmodel            = fmodel,
            x                 = x,
            restraints_weight = None,
            bound_flags       = flex.int(x.size(), 2),
            lower_bound       = lower,
            upper_bound       = upper)
          m = scitbx.minimizers.lbfgs(
            mode='lbfgsb', max_iterations=100, calculator=calculator)
          b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
          model.set_b_iso(values = b_isos)
          if(log is not None):
            print("r_work: %6.4f"%fmodel.r_work(), file=log)
    #
    fmodel.xray_structure.set_b_iso(values = b_isos)
    fmodel.update_xray_structure(xray_structure = fmodel.xray_structure,
      update_f_calc = True)
    if(log is not None):
      print("r_work (final): %6.4f"%fmodel.r_work(), file=log)
    #
    return b_isos
