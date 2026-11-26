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
               log,
               mosaic = False,
               squeeze = False,
               water_selection = None):
    ### Must be first thing here
    self.__dict__.update(fmodel.__dict__)
    # From this point on: self = fmodel
    ###
    russ = self.compute(apply_back_trace = apply_back_trace, remove_outliers =
      remove_outliers, fast = fast, params = params,
      refine_hd_scattering = refine_hd_scattering, log = log,
      mosaic = mosaic, squeeze = squeeze, water_selection = water_selection)
    ### Must be next to last...
    fmodel.__dict__.update(self.__dict__)
    ### ...and this one is last
    self.russ = russ

  def compute(self, apply_back_trace, remove_outliers, fast,
              params, refine_hd_scattering, log, mosaic = False,
              squeeze = False, water_selection = None):
    assert [self.arrays.core_twin, self.twin_law].count(None) in [0,2]
    if(params is None): params = bss.master_params.extract()
    bss.validate_params(
      params  = params,
      fast    = fast,
      twinned = self.twinned(),
      n_masks = len(self.f_masks()))
    self._show_parameters(params = params, fast = fast,
      apply_back_trace = apply_back_trace, log = log)
    self.show(prefix = "start", log = log)
    self.reset_all_scales()
    self.show(prefix = "re-set all scales", log = log)
    if(remove_outliers and not self.twinned()):
      for iii in range(5):
        self.remove_outliers(use_model = False, log = None) # XXX
      self.show(prefix = "remove outliers", log = log)
    result = None
    if(mosaic):
      from mmtbx.bulk_solvent import mosaic as mosaic_model
      reason = mosaic_model.not_applicable(self)
      if(reason is not None):
        mosaic = False
        self.remove_mosaic()
        if(log is not None): print("mosaic: skipped (%s)"%reason, file=log)
    if(squeeze):
      # Determined below, from the model without it
      from mmtbx.bulk_solvent import squeeze as squeeze_model
      self.remove_squeeze()
      reason = squeeze_model.not_applicable(fmodel = self,
        params = self._squeeze_params(params))
      if(reason is not None):
        squeeze = False
        if(log is not None): print("squeeze: skipped (%s)"%reason, file=log)
    elif(self._squeeze_result is not None and not self.squeeze_is_current()):
      self.remove_squeeze()
      if(log is not None):
        print("squeeze: stored contribution dropped (the model changed)",
          file=log)
    if(self.twinned()):
      for cycle in range(2):
        if(log is not None): print("cycle %d:"%cycle, file=log)
        self.update_twin_fraction()
        self.show(prefix = "update twin fraction", log = log)
        result = self.update_solvent_and_scale_twin(log = log,
          refine_hd_scattering = refine_hd_scattering,
          back_trace = (cycle == 1))
    else:
      kwds = dict(
        fast                 = fast,
        params               = params,
        apply_back_trace     = apply_back_trace,
        refine_hd_scattering = refine_hd_scattering,
        log                  = log)
      if(mosaic):
        # Flat model first: the mosaic model starts from its scales
        result = self.update_solvent_and_scale_2(**kwds)
        self._update_mosaic(params = params, log = log)
      elif(self.mosaic_correction is not None and
           not self.mosaic_correction_is_current()):
        self.remove_mosaic()
        if(log is not None):
          print("mosaic: stored correction dropped (the model changed), "
            "flat bulk-solvent model is used", file=log)
      if(self.mosaic_correction is not None):
        result = self.update_solvent_and_scale_2(
          mask_correction = self.mosaic_correction, **kwds)
        self.show(prefix = "mosaic bulk solvent", log = log)
      elif(not mosaic):
        result = self.update_solvent_and_scale_2(**kwds)
      if(squeeze):
        result_squeeze = self._update_squeeze(params = params,
          water_selection = water_selection, log = log, kwds = kwds)
        if(result_squeeze is not None): result = result_squeeze
    #XXX if(remove_outliers and not self.twinned()):
    #XXX   self.remove_outliers(use_model = True, log = None) # XXX
    if(remove_outliers and not self.twinned()):
      for iii in range(5):
        self.remove_outliers(use_model = True, log = None) # XXX
      self.show(prefix = "remove outliers", log = log)
    self._show_scaling_b_values(result = result, log = log)
    return result

  def _update_mosaic(self, params, log):
    """
    Regions of the solvent mask and their occupancies, with the scales of the
    flat bulk-solvent model that self has now. Stores the correction to the
    mask structure factors, or removes a stored one if there is nothing to
    correct.
    """
    from mmtbx.bulk_solvent import mosaic as mosaic_model
    result = mosaic_model.run(
      fmodel = self,
      params = getattr(params, "mosaic", None),
      log    = log)
    self.remove_mosaic()
    if(result is None): return
    f_obs = self.f_obs()
    self._set_mosaic_arrays(
      correction = f_obs.array(data = result.correction),
      signature  = self.mosaic_signature(),
      result     = result)
    if(log is not None): print(result.summary(), file=log)

  def _squeeze_params(self, params):
    result = getattr(params, "squeeze", None)
    if(result is None):
      from mmtbx.bulk_solvent import squeeze as squeeze_model
      result = squeeze_model.default_params()
    return result

  def _update_squeeze(self, params, water_selection, log, kwds):
    """
    Squeeze contribution from the difference map of self, which has the
    scales of the bulk-solvent model (flat or mosaic) and no squeeze
    contribution. Then rounds of: store it, all scales with it, its scale
    and B with the new scales. Returns the result of the last scaling, None
    if there is no contribution.
    """
    from mmtbx.bulk_solvent import squeeze as squeeze_model
    sp = self._squeeze_params(params)
    sq = squeeze_model.run(
      fmodel          = self,
      params          = sp,
      water_selection = water_selection,
      log             = log)
    if(sq is None): return None
    result = None
    for cycle in range(sp.number_of_cycles):
      if(cycle > 0):
        sq = squeeze_model.refit(fmodel = self,
          result = self._squeeze_result, params = sp)
      self._set_squeeze(result = sq, signature = self.mosaic_signature())
      result = self.update_solvent_and_scale_2(
        mask_correction = self.mosaic_correction, **kwds)
      self.show(prefix = "squeeze (k=%.2f B=%.1f)"%(
        self._squeeze_result.k, self._squeeze_result.b), log = log)
    if(log is not None):
      print(self._squeeze_result.summary(), file=log)
    return result

  def _squeeze_follows_overall_b(self, b_adj):
    """
    Overall B moved to atoms by b_adj: the squeeze contribution gets the same
    factor as the atoms and the bulk solvent do, so that F_model is the same.
    """
    sq = self._squeeze_result
    if(sq is None or b_adj is None or b_adj == 0): return
    self._set_squeeze(
      result    = sq.customized_copy(b = sq.b + b_adj),
      signature = self.squeeze_signature)

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

  def _show_parameters(self, params, fast, apply_back_trace, log):
    """
    Tell which protocol is used and with which parameters. apply_back_trace
    is what is going to be done: the parameter and the caller both have to
    ask for it.
    """
    if(log is None): return
    if(self.twinned()):
      print("Bulk-solvent and scaling: twin protocol (mode and slow.* not "
        "used)", file=log)
      self._show_mask_parameters(log = log)
      return
    mode = "slow"
    if(fast): mode = "fast"
    print("Bulk-solvent and scaling: mode=%s"%mode, file=log)
    print("  bulk_solvent=%s anisotropic_scaling=%s apply_back_trace=%s"%(
      str(params.bulk_solvent), str(params.anisotropic_scaling),
      str(bool(apply_back_trace))), file=log)
    self._show_mask_parameters(log = log)
    if(fast): return
    s = params.slow
    if(s.k_sol_b_sol_grid_search):
      print("  grid search: k_sol %.2f-%.2f step %.2f, "
        "b_sol %.1f-%.1f step %.1f"%(s.k_sol_grid_search_min,
        s.k_sol_grid_search_max, s.k_sol_step, s.b_sol_grid_search_min,
        s.b_sol_grid_search_max, s.b_sol_step), file=log)
    else:
      print("  grid search: off", file=log)
    print("  minimization: k_sol,b_sol=%s b_cart=%s "
      "(symmetry constraints=%s)"%(str(s.minimization_k_sol_b_sol),
      str(s.minimization_b_cart), str(s.symmetry_constraints_on_b_cart)),
      file=log)
    print("  bounds: k_sol %.2f-%.2f, b_sol %.1f-%.1f"%(s.k_sol_min,
      s.k_sol_max, s.b_sol_min, s.b_sol_max), file=log)
    fixed = []
    if(s.fix_k_sol is not None): fixed.append("k_sol=%.4f"%s.fix_k_sol)
    if(s.fix_b_sol is not None): fixed.append("b_sol=%.4f"%s.fix_b_sol)
    fbc = s.fix_b_cart
    fbc = [fbc.b11, fbc.b22, fbc.b33, fbc.b12, fbc.b13, fbc.b23]
    if(fbc.count(None) == 0):
      fixed.append("b_cart=%s"%",".join(["%.4f"%b for b in fbc]))
    if(len(fixed) == 0): fixed = ["none"]
    print("  fixed: %s"%" ".join(fixed), file=log)

  def _show_mask_parameters(self, log):
    """
    Mask parameters as they are used (neutron data always with H in the
    mask).
    """
    import mmtbx.masks
    p = self.mask_params
    if(p is None): return
    print("  mask: solvent_radius=%.2f shrink_truncation_radius=%.2f "
      "use_asu_masks=%s"%(p.solvent_radius, p.shrink_truncation_radius,
      str(p.use_asu_masks)), file=log)
    print("  mask: gridding from %s"%mmtbx.masks.gridding_source(
      mask_params = p, d_min = self.f_obs().d_min()), file=log)
    note = ""
    if(self.data_type() == "neutron" and not p.ignore_hydrogens):
      note = " (always False for neutron data)"
    print("  mask: ignore_zero_occupancy_atoms=%s ignore_hydrogens=%s%s"%(
      str(p.ignore_zero_occupancy_atoms), str(p.ignore_hydrogens), note),
      file=log)
    line = "  mask: Fmask_res_high=%.2f n_radial_shells=%d"%(
      p.Fmask_res_high, p.n_radial_shells)
    if(p.n_radial_shells > 1):
      line += " radial_shell_width=%.2f"%p.radial_shell_width
    print(line, file=log)

  def _capture_scaling_b_values(self, fast, result):
    # Capture the anisotropic (b_cart from k_anisotropic) and isotropic
    # (b_isotropic from k_isotropic) B-values of the bulk-solvent scaling for
    # whichever mode actually ran (fast = analytical, slow = minimization).
    # Returned via russ and printed later, after the block of r-factor lines.
    import scitbx.math
    if(fast):
      b_cart = getattr(result, "scale_matrices", None)
      if(b_cart is not None and len(b_cart) > 6): b_cart = None # polynomial scale
    else:
      b_cart = result.b_cart()
    b_isotropic = None
    k_isotropic = result.k_isotropic()
    if(k_isotropic is not None and k_isotropic.size() == self.ss.size()):
      b_isotropic = scitbx.math.gaussian_fit_1d_analytical(
        x = flex.sqrt(self.ss), y = k_isotropic).b
    k_sol, b_sol = None, None
    if(not fast):
      k_sol, b_sol = list(result.k_sols()), list(result.b_sols())
    return group_args(
      mode        = "fast" if fast else "slow",
      k_sol       = k_sol,
      b_sol       = b_sol,
      b_cart      = b_cart,
      b_isotropic = b_isotropic)

  def _show_scaling_b_values(self, result, log):
    i = getattr(result, "scaling_b_info", None)
    if(log is None or i is None): return
    print(file=log)
    print("bulk-solvent scaling B-values (mode: %s):" % i.mode, file=log)
    if(i.k_sol is not None):
      print("  k_sol:       %s" % " ".join(["%.4f" % k for k in i.k_sol]),
        file=log)
      print("  b_sol:       %s" % " ".join(["%.4f" % b for b in i.b_sol]),
        file=log)
    if(i.b_cart is not None):
      print("  b_cart:      %s" % " ".join(["%.4f" % b for b in i.b_cart]),
        file=log)
    if(i.b_isotropic is not None):
      print("  b_isotropic: %.4f" % i.b_isotropic, file=log)

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

  def _h_scattering_grid_search(self, f_obs, f_calc_plus_f_bulk, f_h, k_total,
                                ss, b_mean, r_ref, zero, log):
    """
    Two-step (coarse, fine) grid search of the scale k_h and B-factor b_h of
    the H contribution k_h*exp(-b_h*s^2/4)*F_H added to
    k_total*(F_calc+F_bulk). Returns k_h, b_h, the contribution (zero array if
    no trial point lowers the R-factor r_ref) and the last search object.
    """
    kh, bh, f_part2 = 0, 0, zero
    b_min = int(max(0,b_mean)*0.5)
    b_max = int(b_mean*1.5)
    sc = 1000.
    # Coarse sampling
    kr=[i/sc for i in range(ifloor(0*sc), iceil(1.5*sc)+1, int(0.1*sc))]
    br=[i/sc for i in range(ifloor(b_min*sc), iceil(b_max*sc)+1, int(5.*sc))]
    o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
      f_obs       = f_obs,
      f_calc      = f_calc_plus_f_bulk,
      f_mask      = f_h.data(),
      k_total     = k_total,
      ss          = ss,
      k_sol_range = flex.double(kr),
      b_sol_range = flex.double(br),
      r_ref       = r_ref)
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
      f_obs       = f_obs,
      f_calc      = f_calc_plus_f_bulk,
      f_mask      = f_h.data(),
      k_total     = k_total,
      ss          = ss,
      k_sol_range = flex.double(kr),
      b_sol_range = flex.double(br),
      r_ref       = o.r())
    if(o.updated()):
      f_part2 = f_h.array(data = o.k_mask()*f_h.data())
      kh, bh = o.k_sol(), o.b_sol()
      self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log, r=o.r())
    return kh, bh, f_part2, o

  def update_solvent_and_scale_2(self, fast, params, apply_back_trace,
                                 refine_hd_scattering, log,
                                 mask_correction = None):
    """
    Bulk-solvent and overall scales (fast: per resolution bin; slow: k_sol,
    b_sol, b_cart by minimization), then the contribution of riding H
    (zero occupancy in the model) as k_h*exp(-b_h*s^2/4)*F_H by grid search.

    The scales and the H contribution depend on each other: scales fitted
    with H scattering missing are biased, and so is the H contribution found
    on top of them. In slow mode the two steps are therefore repeated, the
    scales being fitted with the current H contribution included, until k_h
    and b_h stop changing. The overall B is moved to the atoms
    (apply_back_trace) after that, which changes b_sol, b_cart and b_h by
    b_adj but not the model structure factors.

    mask_correction (cctbx.miller.array, optional): added to the mask
    structure factors for the purpose of scaling (mosaic bulk-solvent model:
    F_mask + sum (q_n-1) F_n). F_mask of self is not changed; the scaled
    correction k_mask*mask_correction is stored as f_part1.

    A stored squeeze contribution (self.f_squeeze()) is a fixed part of the
    model during scaling, as F_calc is; it is added to f_part1.
    """
    f_masks = self.f_masks()
    if(mask_correction is not None):
      assert len(f_masks) == 1
      f_masks = [f_masks[0].array(
        data = f_masks[0].data() + mask_correction.data())]
    def f_part1(k_masks):
      data = self.f_squeeze().data()
      if(mask_correction is not None):
        data = data + k_masks[0]*mask_correction.data()
      return zero.array(data = data)
    refine_h = (refine_hd_scattering and
      self.need_to_refine_hd_scattering_contribution())
    if(self.xray_structure is not None):
      # Figure out Fcalc and Fmask based on presence of H
      hd_selection = self.xray_structure.hd_selection()
      xrs_no_h = self.xray_structure.select(~hd_selection)
      xrs_h    = self.xray_structure.select(hd_selection)
    back_trace = (apply_back_trace and len(self.f_masks())==1 and
      self.xray_structure is not None)
    if(apply_back_trace and len(self.f_masks())>1 and log is not None):
      print("  apply_back_trace skipped: more than one mask", file=log)
    alternate = (refine_h and not fast)
    n_cycles = 1
    if(alternate): n_cycles = 5
    # If H scattering is refined then it is assumed that self.f_calc() does
    # not contain H contribution at all.
    f_calc_no_h = self.f_calc()
    zero = f_calc_no_h.array(data=f_calc_no_h.data()*0)
    f_part2 = zero
    f_h = None
    kh, bh = 0, 0
    k_sol, b_sol, b_cart, b_adj = [None,]*4
    for cycle in range(n_cycles):
      # Create data container for scalers. Note: this intentionally ignores
      # previously existing f_part1 and f_part2 (f_part2 is the H contribution
      # found in the previous cycle, zero in the first one).
      fmodel_kbu = mmtbx.f_model.manager_kbu(
        f_obs   = self.f_obs(),
        f_calc  = f_calc_no_h.array(data = f_calc_no_h.data()+f_part2.data()+
          self.f_squeeze().data()),
        f_masks = f_masks,
        ss      = self.ss)
      # Compute k_total and k_mask using one of the two methods (anal or min).
      if(fast): # analytical
        assert len(fmodel_kbu.f_masks)==1
        result = mmtbx.bulk_solvent.scaler.run_simple(
          fmodel_kbu     = fmodel_kbu,
          r_free_flags   = self.r_free_flags(),
          bulk_solvent   = params.bulk_solvent,
          aniso_scale    = params.anisotropic_scaling,
          bin_selections = self.bin_selections)
      else: # using minimization: exp solvent and scale model (k_sol,b_sol,b_cart)
        result = bss.bulk_solvent_and_scales(
          fmodel_kbu = fmodel_kbu,
          params     = params)
        k_sol, b_sol, b_cart = result.k_sols(), result.b_sols(), result.b_cart()
      r_all_from_scaler = result.r_all() # must be here, before apply_back_trace
      if(back_trace and not alternate):
        o = result.apply_back_trace_of_overall_exp_scale_matrix(
          xray_structure = self.xray_structure)
        b_adj = o.b_adj
        if(not fast): b_sol, b_cart = [o.b_sol], o.b_cart
        self.update_xray_structure(
          xray_structure = o.xray_structure,
          update_f_calc  = True)
        f_calc_no_h = self.f_calc()
        self._squeeze_follows_overall_b(b_adj = b_adj)
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
        k_isotropic   = k_isotropic,
        f_part1       = f_part1(k_masks),
        f_part2       = f_part2)
      self.show(prefix = "bulk-solvent and scaling", log = log)
      scaling_b_info = self._capture_scaling_b_values(fast=fast, result=result)
      # Consistency check
      if(not apply_back_trace or alternate):
        assert approx_equal(self.r_all(), r_all_from_scaler)
      if(not refine_h): break
      #
      # Add contribution from H (riding, zero occupancy). This goes to f_part2.
      #
      kh_prev, bh_prev = kh, bh
      # Obsolete previous contribution f_part2
      self.update_core(f_part2 = zero)
      if(f_h is None or not alternate):
        f_h = self.compute_f_calc(
          xray_structure = xrs_h.set_occupancies(value=1).set_b_iso(value = 0))
      # Accumulate all mask contributions: Fcalc_atoms+Fbulk_1+...+Fbulk_N
      data = f_calc_no_h.data() + self.f_squeeze().data()
      for k_mask_, f_mask_ in zip(k_masks, fmodel_kbu.f_masks):
        data = data + k_mask_*f_mask_.data()
      f_calc_plus_f_bulk_no_scales = f_calc_no_h.array(data = data)
      # Consistency check
      assert approx_equal(self.f_model().data(),
        f_calc_plus_f_bulk_no_scales.data()*k_isotropic*k_anisotropic)
      assert approx_equal(
        self.f_model_no_scales().data() + self.f_part1().data(),
        f_calc_plus_f_bulk_no_scales.data())
      # Compute contribution from H (F_H)
      b_mean = flex.mean(xrs_no_h.extract_u_iso_or_u_equiv())*adptbx.u_as_b(1.)
      kh, bh, f_part2, o = self._h_scattering_grid_search(
        f_obs              = fmodel_kbu.f_obs.data(),
        f_calc_plus_f_bulk = f_calc_plus_f_bulk_no_scales.data(),
        f_h                = f_h,
        k_total            = k_isotropic*k_anisotropic,
        ss                 = fmodel_kbu.ss,
        b_mean             = b_mean,
        # R of all reflections: the search computes R with all reflections
        # it is given, and returns r_ref if it finds no H contribution
        r_ref              = self.r_all(),
        zero               = zero,
        log                = log)
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
        f_part1       = f_part1(k_masks),
        f_part2       = f_part2)
      # Make sure what came out of scaling matches what self thinks it really is
      # It must match at least up to 1.e-6.
      self.show(prefix = "add H (%4.2f, %6.2f)"%(kh, bh), log = log)
      if(fast):
        assert approx_equal(result.r_work(), self.r_work(), 1.e-4)
      else:
        assert approx_equal(self.r_all(), o.r()), [self.r_all(), o.r()]
      # Scales in self were fitted with the H contribution of the previous
      # cycle: done once it is the same as the current one.
      if(alternate and kh == kh_prev and bh == bh_prev): break
    if(alternate and back_trace):
      # Move overall B to atoms. H contribution follows the atoms.
      o = result.apply_back_trace_of_overall_exp_scale_matrix(
        xray_structure = self.xray_structure)
      b_adj = o.b_adj
      b_sol, b_cart = [o.b_sol], o.b_cart
      self.update_xray_structure(
        xray_structure = o.xray_structure,
        update_f_calc  = True)
      f_part2 = f_part2.array(data = f_part2.data()*flex.exp(-self.ss*b_adj))
      if(kh != 0): bh = bh + b_adj
      self._squeeze_follows_overall_b(b_adj = b_adj)
      self.update_core(
        k_mask        = result.k_masks(),
        k_anisotropic = result.k_anisotropic(),
        k_isotropic   = result.k_isotropic(),
        f_part1       = f_part1(result.k_masks()),
        f_part2       = f_part2)
      scaling_b_info = self._capture_scaling_b_values(fast=fast, result=result)
      self.show(prefix = "overall B=%s to atoms"%str("%7.2f"%o.b_adj).strip(),
        log = log)
    return group_args(
      k_sol  = k_sol,
      b_sol  = b_sol,
      b_cart = b_cart,
      k_h    = kh,
      b_h    = bh,
      b_adj  = b_adj,
      scaling_b_info = scaling_b_info)

  def update_solvent_and_scale_twin(self, refine_hd_scattering, log,
                                    back_trace = True):
    """
    Grid search of k_sol, b_sol (with analytical anisotropic scale), then of
    the scale and B of the riding H contribution. The H contribution found by
    the previous call is part of the model during the k_sol, b_sol search.
    back_trace: move the overall B to the atoms once all is done (changes
    b_h by b_adj, not the model structure factors).
    """
    if(not self.twinned()): return
    assert len(self.f_masks()) == 1
    # Re-set all scales to unit or zero
    self.show(prefix = "update scales twin start", log = log)
    # H contribution found in the previous cycle (zero in the first one) is
    # part of the model while k_sol and b_sol are searched for: without it
    # they are biased.
    f_calc_data      = self.f_calc().data()      + self.f_part2().data()
    f_calc_data_twin = self.f_calc_twin().data() + self.f_part2_twin().data()
    self.reset_all_scales()
    self.show(prefix = "reset f_part, k_(total,mask)", log = log)
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
      br=[i/sc for i in range(ifloor(b_min*sc),iceil(b_max*sc)+1,int(1.*sc))]
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
    #
    # Move overall B to atoms; solvent and H contributions follow the atoms.
    #
    b_adj = 0
    if(back_trace and o_kbu_sol.updated() and self.xray_structure is not None):
      o = mmtbx.bulk_solvent.scaler.tmp(
        xray_structure = self.xray_structure,
        k_anisotropic  = self.k_anisotropic(),
        k_masks        = [self.k_masks()[0]],
        ss             = self.ss)
      b_adj = o.b_adj
      f_part2, f_part2_twin = self.f_part2(), self.f_part2_twin()
      self.update_xray_structure(
        xray_structure = o.xray_structure,
        update_f_calc  = True)
      self.update(
        k_mask        = o.k_masks,
        k_anisotropic = o.k_anisotropic)
      scale = flex.exp(-self.ss*b_adj)
      self.update_core(
        f_part2      = f_part2.array(data = f_part2.data()*scale),
        f_part2_twin = f_part2_twin.array(data = f_part2_twin.data()*scale))
      if(kh != 0): bh = bh + b_adj
      self.show(prefix = "overall B=%s to atoms"%str("%7.2f"%b_adj).strip(),
        log = log)
    b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(
                             self.f_obs().unit_cell(), o_kbu_sol.u_star()))
    return group_args(
      k_sol  = o_kbu_sol.k_sol(),
      b_sol  = o_kbu_sol.b_sol(),
      b_cart = b_cart,
      k_h    = kh,
      b_h    = bh,
      b_adj  = b_adj)
