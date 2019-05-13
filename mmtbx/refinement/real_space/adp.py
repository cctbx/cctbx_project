from __future__ import division
from __future__ import print_function
import mmtbx.refinement.real_space.utils
import mmtbx.refinement.utils
from scitbx.array_family import flex
from libtbx import adopt_init_args
from cctbx import adptbx
from libtbx import easy_mp
import mmtbx.secondary_structure
from mmtbx import bulk_solvent
from libtbx.test_utils import approx_equal
from mmtbx import masks

class real_space_group_adp_refinery_via_reciprocal_space(object):
  def __init__(self,
               map,
               d_min,
               pdb_hierarchy,
               crystal_symmetry,
               atom_radius,
               use_adp_restraints,
               nproc,
               log=None):
    adopt_init_args(self, locals())
    self.xray_structure = self.pdb_hierarchy.extract_xray_structure(
      crystal_symmetry = self.crystal_symmetry)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    self.chain_selections = []
    for chain in self.pdb_hierarchy.chains():
      self.chain_selections.append(chain.atoms().extract_i_seq())
    # XXX Probably better is to use this...
    # self.chain_selections = mmtbx.secondary_structure.contiguous_ss_selections(
    #     pdb_hierarchy = self.pdb_hierarchy)
    # XXX ...but splitting needs to be improved (can result to tiny fragments).
    #
    self.chain_selections = [pdb_hierarchy.atoms().extract_i_seq()]
    # Soft mask out the map
    rad_smooth = min(5., d_min)
    mask_object = masks.smooth_mask(
      xray_structure = self.xray_structure,
      n_real         = self.map.all(),
      rad_smooth     = rad_smooth)
    self.map = self.map * mask_object.mask_smooth

  def refine_box_with_selected(self, selection=None):
    if(selection is None): selections = self.xray_structure.all_selection()
    ph_box = self.pdb_hierarchy.select(selection)
    ph_box.atoms().reset_i_seq()
    box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure         = self.xray_structure.select(selection),
      map_data               = self.map,
      box_cushion            = self.atom_radius)
    ph_box.adopt_xray_structure(box.xray_structure_box)
    group_adp_sel = []
    for rg in ph_box.residue_groups():
      group_adp_sel.append(rg.atoms().extract_i_seq())
    f_obs_box_complex =  box.box_map_coefficients(d_min = self.d_min)
    f_obs_box = abs(f_obs_box_complex)
    #
    xrs = box.xray_structure_box.deep_copy_scatterers().set_b_iso(value=0)
    assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
    f_calc = f_obs_box.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    # Get overall B estimate
    o = bulk_solvent.complex_f_kb_scaled(
      f1      = f_obs_box_complex.data(),
      f2      = f_calc.data(),
      b_range = flex.double(range(5,505,5)),
      ss      = 1./flex.pow2(f_calc.d_spacings().data()) / 4.)
    #
    xrs = xrs.set_b_iso(value=o.b())
    f_calc = f_obs_box.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    k_isotropic = flex.double(f_calc.data().size(), o.k())
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs_box,
      xray_structure = xrs)
    fmodel.update_core(k_isotropic = k_isotropic)
    fmodel.update(target_name="ls_wunit_k1") # Risky?
    #
    fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
      remove_outliers=False)
    #
    if(self.nproc>1): log = None
    else:             log = self.log
    group_b_manager = mmtbx.refinement.group.manager(
      fmodel                   = fmodel,
      selections               = group_adp_sel,
      convergence_test         = False,
      max_number_of_iterations = 50,
      number_of_macro_cycles   = 3,
      run_finite_differences_test = False,
      use_restraints           = self.use_adp_restraints,
      refine_adp               = True,
      log                      = log)
    return fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)

  def refine(self):
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    if(self.nproc==1):
      for sel in self.chain_selections:
        b_isos_refined = self.refine_box_with_selected(selection=sel)
        b_isos = b_isos.set_selected(sel, b_isos_refined)
    else:
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.refine_box_with_selected,
        args         = self.chain_selections,
        func_wrapper = "buffer_stdout_stderr")
      for i, it in enumerate(stdout_and_results):
        so, b_isos_refined = it
        b_isos = b_isos.set_selected(self.chain_selections[i], b_isos_refined)
        print(so, file=self.log)
    self.xray_structure = self.xray_structure.set_b_iso(values = b_isos)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
