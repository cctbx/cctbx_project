from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import adptbx
from mmtbx import bulk_solvent
from cctbx.array_family import flex
from cctbx import adptbx
import mmtbx
from libtbx import group_args
import mmtbx.arrays
import mmtbx.bulk_solvent.scaler
from libtbx.test_utils import approx_equal
from libtbx.math_utils import ifloor, iceil
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from six.moves import zip, range

class run(mmtbx.f_model.manager):
  """
  This is a very specialized routine to perform complex protocols of updating
  all scales of fmodel, including case of twininng, presence of H and lileky
  more. Inside it pretends to be fmodel proper (done by dictionary updates
  before and after - any better ideas of how to do it nicer are welcome!).
  """

  def __init__(self,
               fmodel,
               apply_back_trace,
               remove_outliers,
               fast,
               params,
               refine_hd_scattering,
               log):
    ### Must be first thing here
    self.__dict__.update(fmodel.__dict__)
    # From this point on: self = fmodel
    ###
    russ = self.compute(apply_back_trace = apply_back_trace, remove_outliers =
      remove_outliers, fast = fast, params = params,
      refine_hd_scattering = refine_hd_scattering, log = log)
    ### Must be next to last...
    fmodel.__dict__.update(self.__dict__)
    ### ...and this one is last
    self.russ = russ

  def compute(self, apply_back_trace, remove_outliers, fast,
              params, refine_hd_scattering, log):
    assert [self.arrays.core_twin, self.twin_law].count(None) in [0,2]
    self.show(prefix = "start", log = log)
    self.reset_all_scales()
    self.show(prefix = "re-set all scales", log = log)
    if(remove_outliers and not self.twinned()):
      for iii in range(5):
        self.remove_outliers(use_model = False, log = None) # XXX
      self.show(prefix = "remove outliers", log = log)
    result = None
    if(self.twinned()):
      for cycle in range(2):
        if(log is not None): print("cycle %d:"%cycle, file=log)
        self.update_twin_fraction()
        self.show(prefix = "update twin fraction", log = log)
        result = self.update_solvent_and_scale_twin(log = log,
          refine_hd_scattering = refine_hd_scattering)
    else:
      result = self.update_solvent_and_scale_2(
        fast                 = fast,
        params               = params,
        apply_back_trace     = apply_back_trace,
        refine_hd_scattering = refine_hd_scattering,
        log                  = log)
    #XXX if(remove_outliers and not self.twinned()):
    #XXX   self.remove_outliers(use_model = True, log = None) # XXX
    if(remove_outliers and not self.twinned()):
      for iii in range(5):
        self.remove_outliers(use_model = True, log = None) # XXX
      self.show(prefix = "remove outliers", log = log)
    return result

  def reset_all_scales(self):
    size = self.f_obs().data().size()
    zero_c = flex.complex_double(size,0)
    zero_d = flex.double(size,0)
    one_d  = flex.double(size,1)
    f_part1_twin = self.f_calc_twin()
    f_part2_twin = self.f_calc_twin()
    if(f_part1_twin is not None):
      f_part1_twin = self.f_calc_twin().array(data=zero_c)
      f_part2_twin = self.f_calc_twin().array(data=zero_c)
    self.update_core(
      f_part1       = self.f_calc().array(data=zero_c),
      f_part2       = self.f_calc().array(data=zero_c),
      f_part1_twin  = f_part1_twin,
      f_part2_twin  = f_part2_twin,
      k_isotropic   = one_d,
      k_anisotropic = one_d,
      k_mask        = [zero_d]*len(self.k_masks()))

  def show(self, prefix, log, r=None):
    if(log is None): return
    if(r is None): r = self.r_all()
    m = "%24s: r(all,work,free)=%6.4f %6.4f %6.4f n_refl.: %d"%(prefix, r,
      self.r_work(), self.r_free(), self.f_obs().data().size())
    if(not self.twinned()):
      print(m, file=log)
    else:
      print(m+" twin_fraction=%4.2f"%self.twin_fraction, file=log)

  def need_to_refine_hd_scattering_contribution(self):
    if(self.xray_structure is None): return False
    refine_hd_scattering = True
    hd_selection = self.xray_structure.hd_selection()
    occ_h_all_zero = self.xray_structure.select(
      hd_selection).scatterers().extract_occupancies().all_eq(0.0) # riding H
    if(self.xray_structure.guess_scattering_type_neutron() or
       hd_selection.count(True)==0 or
       not occ_h_all_zero):
      refine_hd_scattering = False
    return refine_hd_scattering

  def update_solvent_and_scale_2(self, fast, params, apply_back_trace,
                                 refine_hd_scattering, log):
    if(params is None): params = bss.master_params.extract()
    if(self.xray_structure is not None):
      # Figure out Fcalc and Fmask based on presence of H
      hd_selection = self.xray_structure.hd_selection()
      xrs_no_h = self.xray_structure.select(~hd_selection)
      xrs_h    = self.xray_structure.select(hd_selection)
    # Create data container for scalers. If H scattering is refined then it is
    # assumed that self.f_calc() does not contain H contribution at all.
    fmodel_kbu = mmtbx.f_model.manager_kbu(
      f_obs   = self.f_obs(),
      f_calc  = self.f_calc(),
      f_masks = self.f_masks(),
      ss      = self.ss)
    # Compute k_total and k_mask using one of the two methods (anal or min).
    # Note: this intentionally ignores previously existing f_part1 and f_part2.
    #
    k_sol, b_sol, b_cart, b_adj = [None,]*4
    if(fast): # analytical
      assert len(fmodel_kbu.f_masks)==1
      result = mmtbx.bulk_solvent.scaler.run_simple(
        fmodel_kbu     = fmodel_kbu,
        r_free_flags   = self.r_free_flags(),
        bulk_solvent   = params.bulk_solvent,
        aniso_scale    = params.anisotropic_scaling,
        bin_selections = self.bin_selections)
      r_all_from_scaler = result.r_all() # must be here, before apply_back_trace
    else: # using minimization: exp solvent and scale model (k_sol,b_sol,b_cart)
      result = bss.bulk_solvent_and_scales(
        fmodel_kbu = fmodel_kbu,
        params     = params)
      k_sol, b_sol, b_cart = result.k_sols(), result.b_sols(), result.b_cart()
      r_all_from_scaler = result.r_all() # must be here, before apply_back_trace
    if(apply_back_trace and len(fmodel_kbu.f_masks)==1 and
       self.xray_structure is not None):
      o = result.apply_back_trace_of_overall_exp_scale_matrix(
        xray_structure = self.xray_structure)
      b_adj = o.b_adj
      if(not fast): b_sol, b_cart = [o.b_sol], o.b_cart
      self.update_xray_structure(
        xray_structure = o.xray_structure,
        update_f_calc  = True)
      fmodel_kbu = fmodel_kbu.update(f_calc = self.f_calc())
      self.show(prefix = "overall B=%s to atoms"%str("%7.2f"%o.b_adj).strip(),
        log = log)
    # Update self with new arrays so that H correction knows current R factor.
    # If no H to account for, then this is the final result.
    k_masks       = result.k_masks()
    k_anisotropic = result.k_anisotropic()
    k_isotropic   = result.k_isotropic()
    self.update_core(
      k_mask        = k_masks,
      k_anisotropic = k_anisotropic,
      k_isotropic   = k_isotropic)
    self.show(prefix = "bulk-solvent and scaling", log = log)
    # Consistency check
    if(not apply_back_trace):
      assert approx_equal(self.r_all(), r_all_from_scaler)
    # Add contribution from H (if present and riding). This goes to f_part2.
    kh, bh = 0, 0
    if(refine_hd_scattering and
       self.need_to_refine_hd_scattering_contribution()):
      # Obsolete previous contribution f_part2
      f_part2 = fmodel_kbu.f_calc.array(data=fmodel_kbu.f_calc.data()*0)
      self.update_core(f_part2 = f_part2)
      xrs_h = xrs_h.set_occupancies(value=1).set_b_iso(value = 0)
      f_h = self.compute_f_calc(xray_structure = xrs_h)
      # Accumulate all mask contributions: Fcalc_atoms+Fbulk_1+...+Fbulk_N
      data = fmodel_kbu.f_calc.data()
      for k_mask_, f_mask_ in zip(k_masks, fmodel_kbu.f_masks):
        data = data + k_mask_*f_mask_.data()
      f_calc_plus_f_bulk_no_scales = fmodel_kbu.f_calc.array(data = data)
      # Consistency check
      assert approx_equal(self.f_model().data(),
        f_calc_plus_f_bulk_no_scales.data()*k_isotropic*k_anisotropic)
      assert approx_equal(self.f_model_no_scales().data(),
        f_calc_plus_f_bulk_no_scales.data())
      #
      # Compute contribution from H (F_H)
      #
      # Coarse sampling
      b_mean = flex.mean(xrs_no_h.extract_u_iso_or_u_equiv())*adptbx.u_as_b(1.)
      b_min = int(max(0,b_mean)*0.5)
      b_max = int(b_mean*1.5)
      sc = 1000.
      kr=[i/sc for i in range(ifloor(0*sc), iceil(1.5*sc)+1, int(0.1*sc))]
      br=[i/sc for i in range(ifloor(b_min*sc), iceil(b_max*sc)+1, int(5.*sc))]
      o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs       = fmodel_kbu.f_obs.data(),
        f_calc      = f_calc_plus_f_bulk_no_scales.data(),
        f_mask      = f_h.data(),
        k_total     = k_isotropic*k_anisotropic,
        ss          = fmodel_kbu.ss,
        k_sol_range = flex.double(kr),
        b_sol_range = flex.double(br),
        r_ref       = self.r_work())
      if(o.updated()):
        f_part2 = f_h.array(data = o.k_mask()*f_h.data())
        kh, bh = o.k_sol(), o.b_sol()
        self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log, r=o.r())
      # Fine sampling
      k_min = max(0,o.k_sol()-0.1)
      k_max = o.k_sol()+0.1
      b_min = max(0,o.b_sol()-5.)
      b_max = o.b_sol()+5.
      kr=[i/sc for i in range(ifloor(k_min*sc),iceil(k_max*sc)+1,int(0.01*sc))]
      br=[i/sc for i in range(ifloor(b_min*sc),iceil(b_max*sc)+1,int(1.*sc))]
      o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs       = fmodel_kbu.f_obs.data(),
        f_calc      = f_calc_plus_f_bulk_no_scales.data(),
        f_mask      = f_h.data(),
        k_total     = k_isotropic*k_anisotropic,
        ss          = fmodel_kbu.ss,
        k_sol_range = flex.double(kr),
        b_sol_range = flex.double(br),
        r_ref       = o.r())
      if(o.updated()):
        f_part2 = f_h.array(data = o.k_mask()*f_h.data())
        kh, bh = o.k_sol(), o.b_sol()
        self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log, r=o.r())
      # THIS HELPS if fast=true is used, see how it works in reality
      #
      if(fast):
        fmodel_kbu_ = mmtbx.f_model.manager_kbu(
          f_obs   = self.f_obs(),
          f_calc  = f_calc_plus_f_bulk_no_scales,
          f_masks = [f_part2],
          ss      = self.ss)
        result = mmtbx.bulk_solvent.scaler.run_simple(
          fmodel_kbu     = fmodel_kbu_,
          r_free_flags   = self.r_free_flags(),
          bulk_solvent   = params.bulk_solvent,
          aniso_scale    = params.anisotropic_scaling,
          bin_selections = self.bin_selections)
        f_part2 = f_part2.array(data = result.core.k_mask()*f_part2.data())
        k_isotropic   = result.core.k_isotropic*result.core.k_isotropic_exp
        k_anisotropic = result.core.k_anisotropic
      # Update self with final scales
      self.update_core(
        k_mask        = k_masks,
        k_anisotropic = k_anisotropic,
        k_isotropic   = k_isotropic,
        f_part2       = f_part2)
      # Make sure what came out of scaling matches what self thinks it really is
      # It must match at least up to 1.e-6.
      self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log)
      if(fast):
        assert approx_equal(result.r_work(), self.r_work(), 1.e-4)
      else:
        assert approx_equal(self.r_all(), o.r()), [self.r_all(), o.r()]
    return group_args(
      k_sol  = k_sol,
      b_sol  = b_sol,
      b_cart = b_cart,
      k_h    = kh,
      b_h    = bh,
      b_adj  = b_adj)

  def update_solvent_and_scale_twin(self, refine_hd_scattering, log):
    if(not self.twinned()): return
    assert len(self.f_masks()) == 1
    # Re-set all scales to unit or zero
    self.show(prefix = "update scales twin start", log = log)
    self.reset_all_scales()
    self.show(prefix = "reset f_part, k_(total,mask)", log = log)
    f_calc_data      = self.f_calc().data()
    f_calc_data_twin = self.f_calc_twin().data()
    # Initial trial set
    sc = 1000.
    ksr = [i/sc for i in range(ifloor(0*sc), iceil(0.6*sc)+1,  int(0.05*sc))]
    bsr = [i/sc for i in range(ifloor(0*sc), iceil(150.*sc)+1, int(10.*sc))]
    o_kbu_sol = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
      f_obs          = self.f_obs().data(),
      f_calc_1       = f_calc_data,
      f_calc_2       = f_calc_data_twin,
      f_mask_1       = self.arrays.core.f_masks[0].data(),
      f_mask_2       = self.arrays.core_twin.f_masks[0].data(),
      ss             = self.ss,
      twin_fraction  = self.twin_fraction,
      k_sol_range    = flex.double(ksr),
      b_sol_range    = flex.double(bsr),
      miller_indices = self.f_obs().indices(), #XXX ??? What about twin-related?
      unit_cell      = self.f_obs().unit_cell(),
      r_ref          = self.r_all())
    if(o_kbu_sol.updated()):
      self.update(
        k_mask        = o_kbu_sol.k_mask(),
        k_anisotropic = o_kbu_sol.k_anisotropic())
    # Second (finer) trial set
    k_min = max(o_kbu_sol.k_sol()-0.05, 0)
    k_max = min(o_kbu_sol.k_sol()+0.05, 0.6)
    ksr = [i/sc for i in range(ifloor(k_min*sc), iceil(k_max*sc)+1, int(0.01*sc))]
    b_min = max(o_kbu_sol.b_sol()-10, 0)
    b_max = min(o_kbu_sol.b_sol()+10, 150)
    bsr = [i/sc for i in range(ifloor(b_min*sc), iceil(b_max*sc)+1, int(1.*sc))]
    o_kbu_sol = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
      f_obs          = self.f_obs().data(),
      f_calc_1       = f_calc_data,
      f_calc_2       = f_calc_data_twin,
      f_mask_1       = self.arrays.core.f_masks[0].data(),
      f_mask_2       = self.arrays.core_twin.f_masks[0].data(),
      ss             = self.ss,
      twin_fraction  = self.twin_fraction,
      k_sol_range    = flex.double(ksr),
      b_sol_range    = flex.double(bsr),
      miller_indices = self.f_obs().indices(), #XXX ??? What about twin-related?
      unit_cell      = self.f_obs().unit_cell(),
      r_ref          = o_kbu_sol.r())
    if(o_kbu_sol.updated()):
      self.update(
        k_mask        = o_kbu_sol.k_mask(),
        k_anisotropic = o_kbu_sol.k_anisotropic())
      # Disable due to rare failures. Technically they should always match. But
      # since different routines are used tiny disagreements are possible.
      # See examples in : /net/anaconda/raid1/afonine/work/bugs/twin_refinement
      #assert approx_equal(self.r_all(), o_kbu_sol.r(), 1.e-5)
      ##############
      # use apply_back_trace in if below
      if(self.xray_structure is not None):
        o = mmtbx.bulk_solvent.scaler.tmp(
          xray_structure = self.xray_structure,
          k_anisotropic  = o_kbu_sol.k_anisotropic(),
          k_masks        = [o_kbu_sol.k_mask()],
          ss             = self.ss)
        self.update_xray_structure(
          xray_structure = o.xray_structure,
          update_f_calc  = True)
      #############
        self.update(
          k_mask        = o.k_masks,
          k_anisotropic = o.k_anisotropic)

    self.show(prefix = "bulk-solvent and scaling", log = log)
    #
    # Add contribution from H (if present and riding). This goes to f_part2.
    #
    kh, bh = 0, 0
    if(refine_hd_scattering and
       self.need_to_refine_hd_scattering_contribution()):
      hd_selection = self.xray_structure.hd_selection()
      xrs_no_h = self.xray_structure.select(~hd_selection)
      xrs_h    = self.xray_structure.select(hd_selection)
      # Accumulate all mask contributions: Fcalc_atoms+Fbulk_1+...+Fbulk_N
      data = self.f_calc().data()+self.f_masks()[0].data()*self.k_masks()[0]
      f_calc_plus_f_bulk_no_scales = self.f_calc().array(data = data)
      data = self.f_calc_twin().data()+\
        self.f_masks_twin()[0].data()*self.k_masks_twin()[0]
      f_calc_plus_f_bulk_no_scales_twin = self.f_calc_twin().array(data = data)
      # Initial FH contribution
      xrs_h = xrs_h.set_occupancies(value=1).set_b_iso(value = 0)
      f_h = self.compute_f_calc(xray_structure = xrs_h)
      f_h_twin = self.compute_f_calc(xray_structure = xrs_h,
        miller_array = self.f_calc_twin())
      # Coarse sampling
      b_mean = flex.mean(xrs_no_h.extract_u_iso_or_u_equiv())*adptbx.u_as_b(1.)
      b_min = int(max(0,b_mean)*0.5)
      b_max = int(b_mean*1.5)
      sc = 1000.
      kr=[i/sc for i in range(ifloor(0*sc), iceil(1.5*sc)+1, int(0.1*sc))]
      br=[i/sc for i in range(ifloor(b_min*sc), iceil(b_max*sc)+1, int(5.*sc))]
      obj = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs          = self.f_obs().data(),
        f_calc_1       = f_calc_plus_f_bulk_no_scales.data(),
        f_calc_2       = f_calc_plus_f_bulk_no_scales_twin.data(),
        f_mask_1       = f_h.data(),
        f_mask_2       = f_h_twin.data(),
        ss             = self.ss,
        twin_fraction  = self.twin_fraction,
        k_sol_range    = flex.double(kr),
        b_sol_range    = flex.double(br),
        miller_indices = self.f_obs().indices(), # XXX What about twin-related?
        unit_cell      = self.f_obs().unit_cell(),
        r_ref          = self.r_work())
      if(obj.updated()):
        f_part2      = f_h.array(     data = obj.k_mask()*f_h.data())
        f_part2_twin = f_h_twin.array(data = obj.k_mask()*f_h_twin.data())
        kh, bh = obj.k_sol(), obj.b_sol()
      # Fine sampling
      k_min = max(0,obj.k_sol()-0.1)
      k_max = obj.k_sol()+0.1
      b_min = max(0,obj.b_sol()-5.)
      b_max = obj.b_sol()+5.
      kr=[i/sc for i in range(ifloor(k_min*sc),iceil(k_max*sc)+1,int(0.01*sc))]
      br=[i/sc for i in range(ifloor(b_min*sc),iceil(b_max*sc)+1,int(5.*sc))]
      obj = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs          = self.f_obs().data(),
        f_calc_1       = f_calc_plus_f_bulk_no_scales.data(),
        f_calc_2       = f_calc_plus_f_bulk_no_scales_twin.data(),
        f_mask_1       = f_h.data(),
        f_mask_2       = f_h_twin.data(),
        ss             = self.ss,
        twin_fraction  = self.twin_fraction,
        k_sol_range    = flex.double(kr),
        b_sol_range    = flex.double(br),
        miller_indices = self.f_obs().indices(), # XXX What about twin-related?
        unit_cell      = self.f_obs().unit_cell(),
        r_ref          = obj.r())
      if(obj.updated()):
        f_part2      = f_h.array(     data = obj.k_mask()*f_h.data())
        f_part2_twin = f_h_twin.array(data = obj.k_mask()*f_h_twin.data())
        kh, bh = obj.k_sol(), obj.b_sol()
      self.update_core(
        f_part2       = f_part2,
        f_part2_twin  = f_part2_twin,
        k_anisotropic = obj.k_anisotropic())
      self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log)
    b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(
                             self.f_obs().unit_cell(), o_kbu_sol.u_star()))
    return group_args(
      k_sol  = o_kbu_sol.k_sol(),
      b_sol  = o_kbu_sol.b_sol(),
      b_cart = b_cart,
      k_h    = kh,
      b_h    = bh)
