from __future__ import absolute_import, division, print_function
import mmtbx.refinement.real_space.utils
import mmtbx.refinement.utils
from scitbx.array_family import flex
from cctbx import adptbx
from libtbx import easy_mp
from mmtbx import bulk_solvent
from libtbx.test_utils import approx_equal
from six.moves import range

import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")

def map_and_model_to_fmodel(map_data, xray_structure, atom_radius, d_min,
                            reset_adp=True):
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure = xray_structure,
    map_data       = map_data,
    box_cushion    = atom_radius)
  box.apply_mask_inplace(atom_radius = atom_radius)
  f_obs_complex = box.box_map_coefficients(d_min = d_min)
  f_obs = abs(f_obs_complex)
  xrs = box.xray_structure_box.deep_copy_scatterers()
  if(reset_adp):
    xrs = xrs.set_b_iso(value=0)
    assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    o = bulk_solvent.complex_f_kb_scaled(
      f1      = f_obs_complex.data(),
      f2      = f_calc.data(),
      b_range = flex.double(range(5,505,5)),
      ss      = 1./flex.pow2(f_calc.d_spacings().data()) / 4.)
    xrs = xrs.set_b_iso(value=o.b())
    k_isotropic = flex.double(f_calc.data().size(), o.k())
  fmodel = mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs)
  if(reset_adp):
    fmodel.update_core(k_isotropic = k_isotropic)
  fmodel.update(target_name="ls_wunit_k1")
  fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
    remove_outliers=False)
  return fmodel

class ncs_aware_refinement(object):
  def __init__(self, map_model_manager, d_min, atom_radius, nproc=1, log = None):
    self.mmm   = map_model_manager
    self.nproc = nproc
    self.d_min = d_min
    self.atom_radius = atom_radius
    self.log         = log
    #
    if(self.nproc>1): self.log = None
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
    if(self.nproc==1):
      args = [model,]
      return self.run_one_one(args = args)
    else:
      argss = []
      selections = []
      for c in model.get_hierarchy().chains():
        sel = c.atoms().extract_i_seq()
        argss.append([model.select(sel),])
        selections.append(sel) # XXX CAN BE BIG
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.run_one_one,
        args         = argss,
        func_wrapper = "buffer_stdout_stderr")
      values = model.get_b_iso()
      for i, result in enumerate(stdout_and_results):
        values = values.set_selected(selections[i], result[1])
      model.set_b_iso(values = values)

  def run_one_one(self, args):
    model = args[0]
    fmodel = map_and_model_to_fmodel(
      map_data       = self.mmm.map_data().deep_copy(),
      xray_structure = model.get_xray_structure(),
      atom_radius    = self.atom_radius,
      d_min          = self.d_min)
    #
    from mmtbx.refinement import adp_refinement
    adp_iso_params = adp_refinement.adp_restraints_master_params.extract().iso
    energies_adp_iso = self.restraints_manager.energies_adp_iso(
        xray_structure    = fmodel.xray_structure,
        parameters        = adp_iso_params,
        use_u_local_only  = adp_iso_params.use_u_local_only,
        use_hd            = False,
        compute_gradients = True)
    #
    # selections for group ADP
    ph_box = model.get_hierarchy()
    ph_box.atoms().reset_i_seq()
    group_adp_sel = []
    for rg in ph_box.residue_groups():
      group_adp_sel.append(rg.atoms().extract_i_seq())
    #
    group_b_manager = mmtbx.refinement.group.manager(
      fmodel                   = fmodel,
      selections               = group_adp_sel,
      convergence_test         = False,
      max_number_of_iterations = 50,
      number_of_macro_cycles   = 3,
      run_finite_differences_test = False,
      use_restraints           = True,
      refine_adp               = True,
      log                      = self.log)
    return fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
