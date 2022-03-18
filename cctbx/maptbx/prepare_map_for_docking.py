from __future__ import print_function
from __future__ import division

import math
from scitbx.array_family import flex
from scitbx.dtmin.minimizer import Minimizer
from scitbx.dtmin.refinebase import RefineBase
from scitbx.dtmin.reparams import Reparams
from scitbx.dtmin.bounds import Bounds
from cctbx import adptbx
import sys

class RefineCryoemErrors(RefineBase):
  # Set up refinement class for dtmin minimiser (based on Phaser minimiser)
  def __init__(self, mc1, mc2, ssqr_bins, target_spectrum, start_x):
    RefineBase.__init__(self)
    # Precompute data that will be used repeatedly for fgh evaluation
    f1 = flex.abs(mc1.data())
    f2 = flex.abs(mc2.data())
    p1 = mc1.phases().data()
    p2 = mc2.phases().data()
    self.sumfsqr_miller = mc1.customized_copy(data=flex.pow2(f1) + flex.pow2(f2))
    self.sumfsqr_miller.use_binner_of(mc1)
    self.f1f2cos_miller = mc1.customized_copy(data=f1 * f2 * flex.cos(p2 - p1))

    self.n_bins = mc1.binner().n_bins_used()  # Assume consistent binning
    self.ssqr_bins = ssqr_bins
    self.unit_cell = mc1.unit_cell()
    assert (self.n_bins == len(ssqr_bins))
    self.target_spectrum = target_spectrum
    self.start_x = start_x
    self.x = start_x[:]         # Full set of parameters
    recip_params = self.unit_cell.reciprocal_parameters()
    astar = recip_params[0]
    bstar = recip_params[1]
    cstar = recip_params[2]
    self.large_shifts_beta = [astar * astar, bstar * bstar, cstar * cstar, astar * bstar, astar * cstar, bstar * cstar]
    ssqr_max = flex.max(mc1.d_star_sq().data())
    d_min = math.sqrt(1./ssqr_max)
    sigmaSphericity = 100*(d_min/2.5)**2 # 50-200?
    self.sigmaSphericityBeta = flex.double(self.large_shifts_beta) * sigmaSphericity/4

  def sphericity_restraint(self, anisoBeta, do_gradient=True, do_hessian=True,
      sigma_factor = 1.):
    sigmaSphericityBeta = self.sigmaSphericityBeta * sigma_factor
    f = 0.
    g = flex.double(6, 0)
    h = flex.double(6 * 6, 0)
    h.reshape(flex.grid(6, 6))
    unit_cell = self.unit_cell
    aniso_u_iso = adptbx.beta_as_u_iso(unit_cell, anisoBeta)
    aniso_delta_beta = flex.double(adptbx.u_iso_as_beta(unit_cell, aniso_u_iso))
    anisoRemoveIso = flex.double(anisoBeta) - aniso_delta_beta
    dBetaIso_by_dBIso = list(adptbx.u_iso_as_beta(unit_cell, adptbx.b_as_u(1.)))
    mm = unit_cell.metrical_matrix()
    dBIso_by_dBetaAno = list(4./3.*flex.double(mm))

    for ni in range(6):
      f += math.pow(anisoRemoveIso[ni]/sigmaSphericityBeta[ni],2)/2.
      if (do_gradient):
        for nj in range(6):
          dBetaIso_by_dBetaAnoI = dBetaIso_by_dBIso[nj]*dBIso_by_dBetaAno[ni]
          if ni == nj:
            dBetaAno_by_dBetaAnoI = 1.
          else:
            dBetaAno_by_dBetaAnoI = 0.
          g[ni] += (anisoRemoveIso[nj]*(dBetaAno_by_dBetaAnoI-dBetaIso_by_dBetaAnoI) /
                        math.pow(sigmaSphericityBeta[nj],2) )
          if (do_hessian):
            for nk in range(6):
              dBetaIso_by_dBetaAnoI = dBetaIso_by_dBIso[nk]*dBIso_by_dBetaAno[ni]
              if nk == ni:
                dBetaAno_by_dBetaAnoI = 1.
              else:
                dBetaAno_by_dBetaAnoI = 0.
              dBetaIso_by_dBetaAnoJ = dBetaIso_by_dBIso[nk]*dBIso_by_dBetaAno[nj]
              if nk == nj:
                dBetaAno_by_dBetaAnoJ = 1.
              else:
                dBetaAno_by_dBetaAnoJ = 0.
              h[ni,nj] += ( (dBetaAno_by_dBetaAnoJ - dBetaIso_by_dBetaAnoJ) *
                            (dBetaAno_by_dBetaAnoI - dBetaIso_by_dBetaAnoI) /
                            math.pow(sigmaSphericityBeta[nk],2) )
    return (f, g, h)

  def target_gradient_hessian(self, do_gradient=True, do_hessian=True):
    if do_hessian:
      assert (do_gradient)
    # Extract parameters into variables with sensible names
    n_bins = self.n_bins
    i_par = 0
    asqr_scale = self.x[i_par]
    i_par += 1
    sigmaT_bins = self.x[i_par:i_par + n_bins]
    i_par += n_bins
    asqr_beta = tuple(self.x[i_par:i_par + 6])
    i_par += 6
    sigmaE_scale = self.x[i_par]
    i_par += 1
    sigmaE_bins = self.x[i_par:i_par + n_bins]
    i_par += n_bins
    sigmaE_beta = tuple(self.x[i_par:i_par + 6])
    i_par += 6
    assert (i_par == len(self.x))

    # Initialise function, gradient and Hessian with zeros
    f = 0.
    g = flex.double(self.nmp, 0)
    h = flex.double(self.nmp * self.nmp, 0)
    h.reshape(flex.grid(self.nmp, self.nmp))

    # Loop over bins to accumulate target, gradient, Hessian
    ma = self.sumfsqr_miller # Miller array holding associated information
    i_bin_used = 0 # Keep track in case full range of bins not used
    for i_bin in ma.binner().range_used():
      sel = ma.binner().selection(i_bin)
      ma_sel = ma.select(sel)

      # Make Miller array as basis for computing aniso corrections in bin
      # Let u = A^2*sigmaT to simplify computation of derivatives
      ones_array = flex.double(ma_sel.size(), 1)
      all_ones = ma_sel.customized_copy(data=ones_array)
      beta_miller_A = all_ones.apply_debye_waller_factors(
        u_star=adptbx.beta_as_u_star(asqr_beta))
      u_terms = (asqr_scale * sigmaT_bins[i_bin_used]
        * self.target_spectrum[i_bin_used]) * beta_miller_A.data()
      beta_miller_E = all_ones.apply_debye_waller_factors(
        u_star=adptbx.beta_as_u_star(sigmaE_beta))
      sigmaE_terms = ((sigmaE_scale * sigmaE_bins[i_bin_used])
        * beta_miller_E.data())

      # Variance term per reflection is function of these terms
      u2sigE = 2 * u_terms + sigmaE_terms
      var_terms = u2sigE * sigmaE_terms

      f1f2cos = self.f1f2cos_miller.data().select(sel)
      sumfsqr = self.sumfsqr_miller.data().select(sel)
      # Leave out constant twologpi*ma_sel.size()
      minusLL_terms = (sumfsqr * (u_terms + sigmaE_terms)
        - 2 * u_terms * f1f2cos) / var_terms + flex.log(var_terms)
      f += flex.sum(minusLL_terms)

      fgh_a_restraint = self.sphericity_restraint(asqr_beta,
          do_gradient, do_hessian, sigma_factor = 2.) # Looser for amplitude-squared
      f += fgh_a_restraint[0]
      fgh_e_restraint = self.sphericity_restraint(sigmaE_beta,
          do_gradient, do_hessian)
      f += fgh_e_restraint[0]

      if do_gradient:
        # Define some intermediate results needed below
        u2sigE2 = flex.pow2(u2sigE)
        sigmaE_sqr = flex.pow2(sigmaE_terms)
        sumsqrcos = sumfsqr + 2 * f1f2cos
        hyposqr = sumfsqr - 2 * f1f2cos
        if (self.refine_Asqr_scale or self.refine_sigmaT_bins
              or self.refine_Asqr_beta):
          dmLL_by_du_terms = (2 * u2sigE - sumsqrcos) / u2sigE2
        if (self.refine_sigmaE_scale or self.refine_sigmaE_bins
              or self.refine_sigmaE_beta):
          dmLL_by_dsigE_terms = (
            (2 * sigmaE_terms - hyposqr) / (2 * sigmaE_sqr)
            - sumsqrcos / (2 * u2sigE2) + 1. / (u2sigE))
        if self.refine_Asqr_beta or self.refine_sigmaE_beta:
          h_as_double, k_as_double, l_as_double = (
            ma_sel.indices().as_vec3_double().parts())
          hh = flex.pow2(h_as_double)
          kk = flex.pow2(k_as_double)
          ll = flex.pow2(l_as_double)
          hk = h_as_double * k_as_double
          hl = h_as_double * l_as_double
          kl = k_as_double * l_as_double

        i_par = 0 # Keep track of index for unrefined parameters
        i_ref = 0 # Keep track of refined parameters
        if self.refine_Asqr_scale: # Only affects U
          du_by_dAsqr_scale = u_terms / asqr_scale
          g[i_ref] += flex.sum(dmLL_by_du_terms * du_by_dAsqr_scale)
          i_ref += 1
        i_par += 1
        if self.refine_sigmaT_bins: # Only affects U, just current bin
          du_by_dsigmaT_bin = (asqr_scale
            * self.target_spectrum[i_bin_used]) * beta_miller_A.data()
          i_sigmaT_bin = i_ref+i_bin_used # Save for restraint terms below
          g[i_sigmaT_bin] += flex.sum(dmLL_by_du_terms*du_by_dsigmaT_bin)
          i_ref += self.n_bins
        i_par += self.n_bins
        if self.refine_Asqr_beta:  # Only affects U
          hh_factors = [-hh, -kk, -ll, -2*hk, -2*hl, -2*kl]
          du_by_dbetaA = []
          for i_beta in range(6):
            du_by_dbetaA.append(hh_factors[i_beta]*u_terms)
          for i_beta in range(6):
            g[i_ref+i_beta] += flex.sum(dmLL_by_du_terms * du_by_dbetaA[i_beta])
            g[i_ref+i_beta] += fgh_a_restraint[1][i_beta] # Restraint term
          i_ref += 6
        i_par += 6

        # Note that sigmaE_scale is fixed if sigmaE_bins and/or
        # sigmaE_beta are refined
        if self.refine_sigmaE_scale:
          dsigE_by_dscaleE = sigmaE_terms/sigmaE_scale
          g[i_ref] += flex.sum(dmLL_by_dsigE_terms * dsigE_by_dscaleE)
          i_ref += 1
        i_par += 1
        if self.refine_sigmaE_bins: # Just current bin
          dsigE_by_dsigmaE_bin = sigmaE_scale*beta_miller_E.data()
          g[i_ref+i_bin_used] += flex.sum(
            dmLL_by_dsigE_terms * dsigE_by_dsigmaE_bin)
          i_ref += self.n_bins
        i_par += self.n_bins
        if self.refine_sigmaE_beta: # Only affects SigmaE
          dsigE_by_dbetaE = []
          for i_beta in range(6):
            dsigE_by_dbetaE.append(hh_factors[i_beta]*sigmaE_terms)
          for i_beta in range(6):
            g[i_ref+i_beta] += flex.sum(dmLL_by_dsigE_terms * dsigE_by_dbetaE[i_beta])
            g[i_ref+i_beta] += fgh_e_restraint[1][i_beta] # Restraint term

          i_ref += 6
        i_par += 6
        assert (i_par == len(self.x))
        assert (i_ref == self.nmp)

        if do_hessian:
          u2sigE3 = u2sigE * u2sigE2
          if (self.refine_Asqr_scale or self.refine_sigmaT_bins
                or self.refine_Asqr_beta):
            d2mLL_by_du2_terms = 4 * (
              sumsqrcos - 2 * u_terms - sigmaE_terms) / u2sigE3
          if (self.refine_sigmaE_scale or self.refine_sigmaE_bins
                or self.refine_sigmaE_beta):
            d2mLL_by_dsigE2_terms = ( (hyposqr - sigmaE_terms)
              / (sigmaE_sqr * sigmaE_terms) + sumsqrcos/u2sigE3 - 1./u2sigE2)

          i_par = 0 # Keep track of index for unrefined parameters
          i_ref = 0  # Keep track of refined parameters
          # Note that various second derivatives are zero, i.e. of:
          #   u wrt Asqr_scale and sigmaT_bins
          #   sigmaE wrt sigmaE_bins and sigmaE_scale
          if self.refine_Asqr_scale: # Only affects U
            h[i_ref,i_ref] += (flex.sum(d2mLL_by_du2_terms
              * flex.pow2(du_by_dAsqr_scale)))
            i_ref += 1
          i_par += 1
          if self.refine_sigmaT_bins: # Only affects U, current bin
            h[i_sigmaT_bin,i_sigmaT_bin] += flex.sum(
              d2mLL_by_du2_terms * flex.pow2(du_by_dsigmaT_bin))
            i_ref += self.n_bins
          i_par += self.n_bins
          if self.refine_Asqr_beta:  # Only affects U
            for i_beta in range(6):
              for j_beta in range(6):
                h[i_ref+i_beta, i_ref+j_beta] += (
                  flex.sum(d2mLL_by_du2_terms * du_by_dbetaA[i_beta]*du_by_dbetaA[j_beta])
                  + flex.sum(dmLL_by_du_terms * hh_factors[i_beta]*hh_factors[j_beta]*u_terms) )
                # Also add restraint term
                h[i_ref+i_beta,i_ref+j_beta] += fgh_a_restraint[2][i_beta, j_beta]
            i_ref += 6
          i_par += 6

          # Note that sigmaE_scale and either sigmaE_bins or sigmaE_beta are
          # mutually exclusive in practice. If scale and bins were refined
          # simultaneously with no restraints, there would be one redundant parameter
          if self.refine_sigmaE_scale:
            h[i_ref, i_ref] += flex.sum(
              d2mLL_by_dsigE2_terms * flex.pow2(dsigE_by_dscaleE))
            i_ref += 1
          i_par += 1
          if self.refine_sigmaE_bins:
            h[i_ref+i_bin_used, i_ref+i_bin_used] += flex.sum(
              d2mLL_by_dsigE2_terms * flex.pow2(dsigE_by_dsigmaE_bin))
            i_ref += self.n_bins
          i_par += self.n_bins
          if self.refine_sigmaE_beta: # Only affects SigmaE
            for i_beta in range(6):
              for j_beta in range(6):
                h[i_ref+i_beta, i_ref+j_beta] += (
                  flex.sum(d2mLL_by_dsigE2_terms * dsigE_by_dbetaE[i_beta]*dsigE_by_dbetaE[j_beta])
                  + flex.sum(dmLL_by_dsigE_terms * hh_factors[i_beta]*hh_factors[j_beta]*sigmaE_terms) )
                h[i_ref+i_beta,i_ref+j_beta] += fgh_e_restraint[2][i_beta, j_beta]
            i_ref += 6
          i_par += 6
          assert (i_par == len(self.x))
          assert (i_ref == self.nmp)

      # Add other restraint terms
      # Restrain log of sigmaT_bins to 0, but downweighting low resolution
      d_bin = math.sqrt(self.ssqr_bins[i_bin_used])
      sigmascale = 0.15 + 0.001 / d_bin ** 3
      stbin = sigmaT_bins[i_bin_used]
      logbin = math.log(stbin)
      f += (logbin/sigmascale)**2 / 2
      if do_gradient and self.refine_sigmaT_bins:
        g[i_sigmaT_bin] += logbin / (stbin * sigmascale**2)
        if do_hessian:
          h[i_sigmaT_bin,i_sigmaT_bin] += (1.-logbin)/(stbin*sigmascale)**2

      i_bin_used += 1

    return (f, g, h, False)

  def target(self):
    f_g_h = self.target_gradient_hessian(do_gradient=False, do_hessian=False)
    return f_g_h[0]

  def target_gradient(self):
    f_g_h = self.target_gradient_hessian(do_hessian=False)
    f = f_g_h[0]
    g = f_g_h[1]
    return (f, g)

  def get_macrocycle_parameters(self):
    if len(self.refine_mask) == 0: # All parameters being refined
      return self.x

    mp = []  # Parameters for this macrocycle
    for i in range(len(self.x)):
      if self.refine_mask[i]:
        mp.append(self.x[i])
    assert (len(mp) == self.nmp)
    return mp

  def set_macrocycle_parameters(self, newx):
    if len(self.refine_mask) == 0:  # All parameters being refined
      self.x = newx
    else:
      npref = 0
      for i in range(len(self.x)):
        if self.refine_mask[i]:
          self.x[i] = newx[npref]
          npref += 1
      assert (npref == self.nmp)

  def macrocycle_large_shifts(self):
    i_par = 0 # Keep track of index for unrefined parameters
    large_shifts = []
    if self.refine_Asqr_scale:
      large_shifts.append(self.start_x[i_par]/10.)
    i_par += 1
    if self.refine_sigmaT_bins:
      for i_bin in range(self.n_bins):
        large_shifts.append(0.05)
    i_par += self.n_bins
    if self.refine_Asqr_beta:
      large_shifts.extend(self.large_shifts_beta)
    i_par += 6
    if self.refine_sigmaE_scale:
      large_shifts.append(0.01)
    i_par += 1
    if self.refine_sigmaE_bins:
      for i_bin in range(self.n_bins):
        large_shifts.append(self.start_x[i_par+i_bin]/30.)
    i_par += self.n_bins
    if self.refine_sigmaE_beta:
      large_shifts.extend(self.large_shifts_beta)
    i_par += 6
    assert (i_par == len(self.x))
    assert (len(large_shifts) == self.nmp)
    return large_shifts

  def set_macrocycle_protocol(self, macrocycle_protocol):
    # Possible parameters include overall scale of signal,
    # bin parameters for signal (BEST-like curve), anisotropy tensor for signal,
    # bin parameters for error, anisotropy tensor for error
    # Start with everything being refined, turn some things off for different protocols
    self.refine_mask = []  # Indicates "all" if left empty
    self.refine_Asqr_scale = True
    self.refine_sigmaT_bins = True
    self.refine_Asqr_beta = True
    self.refine_sigmaE_scale = True
    self.refine_sigmaE_bins = True
    self.refine_sigmaE_beta = True

    # For each protocol, define variables that aren't refined
    if macrocycle_protocol == ["default"]:
      self.refine_sigmaE_scale = False

    elif macrocycle_protocol == ["Eprior"]:
      self.refine_sigmaE_bins = False
      self.refine_sigmaE_beta = False

    else:
      print("Macrocycle protocol", macrocycle_protocol, " not recognised")
      sys.stdout.flush()
      exit

    # Now accumulate mask
    self.nmp = 0

    if self.refine_Asqr_scale:
      self.refine_mask.append(True)
      self.nmp += 1
    else:
      self.refine_mask.append(False)

    if self.refine_sigmaT_bins:
      self.refine_mask.extend([True for i in range(self.n_bins)])
      self.nmp += self.n_bins
    else:
      self.refine_mask.extend([False for i in range(self.n_bins)])

    if self.refine_Asqr_beta:
      self.refine_mask.extend([True for i in range(6)])
      self.nmp += 6
    else:
      self.refine_mask.extend([False for i in range(6)])

    if self.refine_sigmaE_scale:
      self.refine_mask.append(True)
      self.nmp += 1
    else:
      self.refine_mask.append(False)

    if self.refine_sigmaE_bins:
      self.refine_mask.extend([True for i in range(self.n_bins)])
      self.nmp += self.n_bins
    else:
      self.refine_mask.extend([False for i in range(self.n_bins)])

    if self.refine_sigmaE_beta:
      self.refine_mask.extend([True for i in range(6)])
      self.nmp += 6
    else:
      self.refine_mask.extend([False for i in range(6)])

    assert (len(self.refine_mask) == len(self.x))

  def macrocycle_parameter_names(self, full_list=False):
    parameter_names = []
    if full_list or self.refine_Asqr_scale:
      parameter_names.append("Asqr_scale")
    if full_list or self.refine_sigmaT_bins:
      for i in range(self.n_bins):
        parameter_names.append("SigmaT_bin#" + str(i + 1))
    if full_list or self.refine_Asqr_beta:
      parameter_names.append("Asqr_beta11")
      parameter_names.append("Asqr_beta22")
      parameter_names.append("Asqr_beta33")
      parameter_names.append("Asqr_beta12")
      parameter_names.append("Asqr_beta13")
      parameter_names.append("Asqr_beta23")
    if full_list or self.refine_sigmaE_scale:
      parameter_names.append("sigmaE_scale")
    if full_list or self.refine_sigmaE_bins:
      for i in range(self.n_bins):
        parameter_names.append("sigmaE_bin#" + str(i + 1))
    if full_list or self.refine_sigmaE_beta:
      parameter_names.append("sigmaE_beta11")
      parameter_names.append("sigmaE_beta22")
      parameter_names.append("sigmaE_beta33")
      parameter_names.append("sigmaE_beta12")
      parameter_names.append("sigmaE_beta13")
      parameter_names.append("sigmaE_beta23")

    if not full_list:
      assert (len(parameter_names) == self.nmp)
    else:
      assert (len(parameter_names) == len(self.x))
    return parameter_names

  def reparameterize(self):
    i_par = 0 # Keep track of index for unrefined parameters
    repar = []

    if self.refine_Asqr_scale:
      repar.append(Reparams(True,0.))
    i_par += 1

    if self.refine_sigmaT_bins:
      repar.extend([Reparams(True, 0.) for i in range(self.n_bins)])
    i_par += self.n_bins

    if self.refine_Asqr_beta:
      repar.extend([Reparams(False) for i in range(6)])
    i_par += 6

    if self.refine_sigmaE_scale:
      repar.append(Reparams(True,0.))
    i_par += 1

    if self.refine_sigmaE_bins:
      repar.extend([Reparams(True, 0.) for i in range(self.n_bins)])
    i_par += self.n_bins

    if self.refine_sigmaE_beta:
      repar.extend([Reparams(False) for i in range(6)])
    i_par += 6

    assert (i_par == len(self.x))
    assert (len(repar) == self.nmp)

    return repar

  def bounds(self):
    i_par = 0
    bounds_list = []

    if self.refine_Asqr_scale:
      this_bound = Bounds()
      this_bound.lower_on(0.001*self.start_x[i_par])
      bounds_list.append(this_bound)
    i_par += 1

    if self.refine_sigmaT_bins:
      this_bound = Bounds()
      this_bound.lower_on(0.001)
      for i in range(self.n_bins):
        bounds_list.append(this_bound)
    i_par += self.n_bins

    if self.refine_Asqr_beta:
      this_bound = Bounds()
      this_bound.off()
      for i in range(6):
        bounds_list.append(this_bound)
    i_par += 6

    if self.refine_sigmaE_scale:
      this_bound = Bounds()
      this_bound.lower_on(0.01)
      bounds_list.append(this_bound)
    i_par += 1

    if self.refine_sigmaE_bins:
      for i_bin in range(self.n_bins):
        this_bound = Bounds()
        this_bound.lower_on(0.001*self.start_x[i_par+i_bin])
        bounds_list.append(this_bound)
    i_par += self.n_bins

    if self.refine_sigmaE_beta:
      this_bound = Bounds()
      this_bound.off()
      for i in range(6):
        bounds_list.append(this_bound)
    i_par += 6

    assert (i_par == len(self.x))
    assert (len(bounds_list) == self.nmp)
    return bounds_list

  def current_statistics(self, level=3, full_list=False):
    self.log_tab_printf(1, level, "Log-likelihood: %10.6g\n", -self.target())
    self.log_blank(level)

    parameter_names = self.macrocycle_parameter_names(full_list=full_list)
    if full_list:
      self.log_tab(1, level, "All parameters")
    else:
      self.log_tab(1, level, "Refined parameters")
    list_all = (full_list or len(self.refine_mask) == 0)
    iref = 0
    for i in range(len(self.x)):
      if list_all or self.refine_mask[i]:
        self.log_tab_printf(2, level, "%-15s %10.5g\n", (parameter_names[iref], self.x[i]))
        iref += 1

  def initial_statistics(self):
    level=2
    self.log_blank(level)
    self.log_tab(1, level, "Initial statistics")
    self.current_statistics(level=level, full_list=True)

  def final_statistics(self):
    level=2
    self.log_blank(level)
    self.log_tab(1, level, "Final statistics")
    self.current_statistics(level=level, full_list=True)

  def cleanup(self):
    # Take out overall scale and isotropic B from sigmaT_bins, put into asqr_scale and asqr_beta
    # Take out overall isotropic B from anisotropy in errors, put into bins

    n_bins = self.n_bins
    # Asqr_scale = self.x[0] # Unneeded parameters listed for completeness
    sigmaT_bins = self.x[1:n_bins + 1]
    # Asqr_beta = self.x[n_bins + 1 : n_bins + 7]
    sigmaE_scale = self.x[n_bins + 7]
    # sigmaE_bins = self.x[n_bins + 8 : 2*n_bins + 8]
    sigmaE_beta = tuple(self.x[2*n_bins + 8 : 2*n_bins + 14])
    sumw = sumwx = sumwa = sumwx2 = sumwxa = 0.
    for i_bin in range(n_bins):
      x = self.ssqr_bins[i_bin]
      d_bin = math.sqrt(x)
      a = math.log(sigmaT_bins[i_bin])
      sigmascale = 0.15 + 0.001 / d_bin**3  # Downweight low resolution as in refinement
      w = 1./sigmascale**2
      sumw += w
      sumwx += w * x
      sumwa += w * a
      sumwx2 += w * x ** 2
      sumwxa += w * x * a

    if self.refine_sigmaT_bins:
      # Make sigmaT_bins values as close as possible to 1 by taking out overall scale
      # and B, and putting them into asqr_scale and asqr_beta terms
      slope_a = (sumw * sumwxa - (sumwx * sumwa)) / (sumw * sumwx2 - sumwx ** 2)
      intercept_a = (sumwa - slope_a * sumwx) / sumw
      scale_a = math.exp(intercept_a)
      deltaB_a = -4 * slope_a
      self.x[0] = self.x[0] * scale_a  # Update overall scale
      for i_bin in range(n_bins): # Take slope out of sigmaT_bins
        self.x[1 + i_bin] = self.x[1 + i_bin] / scale_a * math.exp(deltaB_a * self.ssqr_bins[i_bin] / 4)
      delta_beta_a = list(adptbx.u_iso_as_beta(self.unit_cell,adptbx.b_as_u(deltaB_a)))
      for i_beta in range(6):  # Then put slope into asqr_beta
        self.x[1 + n_bins + i_beta] = self.x[1 + n_bins + i_beta] + delta_beta_a[i_beta]

    if not (sigmaE_scale == 1.): # Put change into bins before next step
      assert (sigmaE_scale > 0.)
      for i_bin in range(n_bins):
        self.x[8 + n_bins + i_bin] = self.x[8 + n_bins + i_bin] / sigmaE_scale
      self.x[7 + n_bins] = 1.
      sigmaE_scale = 1.

    if self.refine_sigmaE_beta:
      # Extract isotropic B from sigmaE_beta, put it into sigmaE_bins
      sigmaE_u_iso = adptbx.beta_as_u_iso(self.unit_cell, sigmaE_beta)
      sigmaE_delta_beta = list(adptbx.u_iso_as_beta(self.unit_cell, sigmaE_u_iso))
      sigmaE_b_iso = adptbx.u_as_b(sigmaE_u_iso)
      for i_bin in range(n_bins):  # Put isotropic B into bins
        self.x[8 + n_bins + i_bin] = (self.x[8 + n_bins + i_bin]
          * math.exp(-sigmaE_b_iso * self.ssqr_bins[i_bin] / 4))
      for i_beta in range(6):  # Remove isotropic B from sigmaE_beta
        self.x[8 + 2*n_bins + i_beta] = (
          self.x[8 + 2*n_bins + i_beta] - sigmaE_delta_beta[i_beta])

def default_target_spectrum(ssqr):
  # Placeholder for something better based on analysis of cryoEM reconstructions
  # Scaled data from BEST curve. Original data obtained from Sasha Popov, then
  # rescaled to correspond at higher resolution to the average X-ray scattering
  # factor from proteins atoms (with average atomic composition)
  best_data = ((0.009, 3.40735),
              (0.013092, 2.9006),
              (0.0171839, 2.33083),
              (0.0212759, 1.80796),
              (0.0253679, 1.65133),
              (0.0294599, 1.75784),
              (0.0335518, 2.06865),
              (0.0376438, 2.57016),
              (0.0417358, 3.13121),
              (0.0458278, 3.62596),
              (0.0499197, 3.92071),
              (0.0540117, 3.98257),
              (0.0581037, 3.91846),
              (0.0621956, 3.80829),
              (0.0662876, 3.69517),
              (0.0703796, 3.59068),
              (0.0744716, 3.44971),
              (0.0785635, 3.30765),
              (0.0826555, 3.16069),
              (0.0867475, 2.98656),
              (0.0908395, 2.77615),
              (0.0949314, 2.56306),
              (0.0990234, 2.37314),
              (0.103115, 2.22874),
              (0.107207, 2.09477),
              (0.111299, 1.98107),
              (0.115391, 1.8652),
              (0.119483, 1.75908),
              (0.123575, 1.67093),
              (0.127667, 1.59257),
              (0.131759, 1.52962),
              (0.135851, 1.48468),
              (0.139943, 1.45848),
              (0.144035, 1.43042),
              (0.148127, 1.40953),
              (0.152219, 1.37291),
              (0.156311, 1.34217),
              (0.160403, 1.3308),
              (0.164495, 1.32782),
              (0.168587, 1.30862),
              (0.172679, 1.31319),
              (0.176771, 1.30907),
              (0.180863, 1.31456),
              (0.184955, 1.31055),
              (0.189047, 1.31484),
              (0.193139, 1.31828),
              (0.197231, 1.32321),
              (0.201323, 1.30853),
              (0.205415, 1.30257),
              (0.209507, 1.2851),
              (0.213599, 1.26912),
              (0.217691, 1.24259),
              (0.221783, 1.24119),
              (0.225875, 1.2382),
              (0.229967, 1.21605),
              (0.234059, 1.17269),
              (0.23815, 1.13909),
              (0.242242, 1.1165),
              (0.246334, 1.08484),
              (0.250426, 1.0495),
              (0.254518, 1.01289),
              (0.25861, 0.974819),
              (0.262702, 0.940975),
              (0.266794, 0.900938),
              (0.270886, 0.861657),
              (0.274978, 0.830192),
              (0.27907, 0.802167),
              (0.283162, 0.780746),
              (0.287254, 0.749194),
              (0.291346, 0.720884),
              (0.295438, 0.694409),
              (0.29953, 0.676239),
              (0.303622, 0.650672),
              (0.307714, 0.632438),
              (0.311806, 0.618569),
              (0.315898, 0.605762),
              (0.31999, 0.591398),
              (0.324082, 0.579308),
              (0.328174, 0.572076),
              (0.332266, 0.568138),
              (0.336358, 0.559537),
              (0.34045, 0.547927),
              (0.344542, 0.539319),
              (0.348634, 0.529009),
              (0.352726, 0.516954),
              (0.356818, 0.512218),
              (0.36091, 0.511836),
              (0.365002, 0.511873),
              (0.369094, 0.506957),
              (0.373186, 0.502738),
              (0.377278, 0.50191),
              (0.38137, 0.492422),
              (0.385462, 0.488461),
              (0.389553, 0.483436),
              (0.393645, 0.481468),
              (0.397737, 0.473786),
              (0.401829, 0.468684),
              (0.405921, 0.468291),
              (0.410013, 0.46645),
              (0.414105, 0.4643),
              (0.418197, 0.45641),
              (0.422289, 0.450462),
              (0.426381, 0.444678),
              (0.430473, 0.443807),
              (0.434565, 0.441158),
              (0.438657, 0.441303),
              (0.442749, 0.437144),
              (0.446841, 0.428504),
              (0.450933, 0.420459),
              (0.455025, 0.413754),
              (0.459117, 0.412064),
              (0.463209, 0.406677),
              (0.467301, 0.40253),
              (0.471393, 0.396454),
              (0.475485, 0.393192),
              (0.479577, 0.390452),
              (0.483669, 0.38408),
              (0.487761, 0.379456),
              (0.491853, 0.373123),
              (0.495945, 0.374026),
              (0.500037, 0.373344),
              (0.504129, 0.377639),
              (0.508221, 0.374029),
              (0.512313, 0.374691),
              (0.516405, 0.371632),
              (0.520497, 0.370724),
              (0.524589, 0.366095),
              (0.528681, 0.369447),
              (0.532773, 0.369043),
              (0.536865, 0.368967),
              (0.540956, 0.36583),
              (0.545048, 0.370593),
              (0.54914, 0.371047),
              (0.553232, 0.372723),
              (0.557324, 0.371915),
              (0.561416, 0.372882),
              (0.565508, 0.371052),
              (0.5696, 0.36775),
              (0.573692, 0.369884),
              (0.577784, 0.374098),
              (0.581876, 0.374169),
              (0.585968, 0.37261),
              (0.59006, 0.372356),
              (0.594152, 0.377055),
              (0.598244, 0.3817),
              (0.602336, 0.381867),
              (0.606428, 0.377746),
              (0.61052, 0.377157),
              (0.614612, 0.376604),
              (0.618704, 0.37532),
              (0.622796, 0.372488),
              (0.626888, 0.373312),
              (0.63098, 0.377505),
              (0.635072, 0.381011),
              (0.639164, 0.379326),
              (0.643256, 0.380193),
              (0.647348, 0.381122),
              (0.65144, 0.387213),
              (0.655532, 0.391928),
              (0.659624, 0.398986),
              (0.663716, 0.402951),
              (0.667808, 0.405893),
              (0.6719, 0.40217),
              (0.675992, 0.401806),
              (0.680084, 0.404238),
              (0.684176, 0.409404),
              (0.688268, 0.413486),
              (0.692359, 0.413167),
              (0.696451, 0.414008),
              (0.700543, 0.417128),
              (0.704635, 0.420275),
              (0.708727, 0.423617),
              (0.712819, 0.42441),
              (0.716911, 0.426445),
              (0.721003, 0.429012),
              (0.725095, 0.430132),
              (0.729187, 0.42992),
              (0.733279, 0.425202),
              (0.737371, 0.423159),
              (0.741463, 0.423913),
              (0.745555, 0.425542),
              (0.749647, 0.426682),
              (0.753739, 0.431186),
              (0.757831, 0.433959),
              (0.761923, 0.433839),
              (0.766015, 0.428679),
              (0.770107, 0.425968),
              (0.774199, 0.426528),
              (0.778291, 0.427093),
              (0.782383, 0.426848),
              (0.786475, 0.424549),
              (0.790567, 0.423785),
              (0.794659, 0.419892),
              (0.798751, 0.417391),
              (0.802843, 0.413128),
              (0.806935, 0.408498),
              (0.811027, 0.402764),
              (0.815119, 0.404852),
              (0.819211, 0.405915),
              (0.823303, 0.392919),
              (0.827395, 0.384632),
              (0.831487, 0.382626),
              (0.835579, 0.379891),
              (0.839671, 0.376414),
              (0.843762, 0.372915),
              (0.847854, 0.375089),
              (0.851946, 0.371918),
              (0.856038, 0.36652),
              (0.86013, 0.358529),
              (0.864222, 0.356496),
              (0.868314, 0.354707),
              (0.872406, 0.348802),
              (0.876498, 0.343693),
              (0.88059, 0.34059),
              (0.884682, 0.342432),
              (0.888774, 0.345099),
              (0.892866, 0.344524),
              (0.896958, 0.342489),
              (0.90105, 0.328009),
              (0.905142, 0.323685),
              (0.909234, 0.321378),
              (0.913326, 0.318832),
              (0.917418, 0.314999),
              (0.92151, 0.311775),
              (0.925602, 0.30844),
              (0.929694, 0.30678),
              (0.933786, 0.303484),
              (0.937878, 0.301197),
              (0.94197, 0.296788),
              (0.946062, 0.295353),
              (0.950154, 0.298028),
              (0.954246, 0.298098),
              (0.958338, 0.295081),
              (0.96243, 0.289337),
              (0.966522, 0.286116),
              (0.970614, 0.284319),
              (0.974706, 0.280972),
              (0.978798, 0.28015),
              (0.98289, 0.279016),
              (0.986982, 0.277532),
              (0.991074, 0.276013),
              (0.995165, 0.270923),
              (0.999257, 0.269446),
              (1.00335, 0.266567),
              (1.00744, 0.263561),
              (1.01153, 0.261002),
              (1.01563, 0.255349),
              (1.01972, 0.258644),
              (1.02381, 0.254974),
              (1.0279, 0.2523),
              (1.03199, 0.244489),
              (1.03609, 0.249418),
              (1.04018, 0.249519),
              (1.04427, 0.249316),
              (1.04836, 0.249197),
              (1.05245, 0.24415),
              (1.05655, 0.244556),
              (1.06064, 0.241169),
              (1.06473, 0.238484),
              (1.06882, 0.2392),
              (1.07291, 0.240651),
              (1.077, 0.243724),
              (1.0811, 0.243174),
              (1.08519, 0.239545),
              (1.08928, 0.239106),
              (1.09337, 0.238763),
              (1.09746, 0.238971),
              (1.10156, 0.229925),
              (1.10565, 0.225123),
              (1.10974, 0.226932),
              (1.11383, 0.23118),
              (1.11792, 0.228654),
              (1.12202, 0.225084),
              (1.12611, 0.225866),
              (1.1302, 0.227717),
              (1.13429, 0.229508),
              (1.13838, 0.227977),
              (1.14248, 0.226799),
              (1.14657, 0.228456),
              (1.15066, 0.22383),
              (1.15475, 0.22188),
              (1.15884, 0.219986),
              (1.16294, 0.217418),
              (1.16703, 0.214356),
              (1.17112, 0.211027),
              (1.17521, 0.210011),
              (1.1793, 0.210609),
              (1.1834, 0.210893),
              (1.18749, 0.212583),
              (1.19158, 0.208415),
              (1.19567, 0.204557),
              (1.19976, 0.198068),
              (1.20386, 0.197603),
              (1.20795, 0.196691),
              (1.21204, 0.200617),
              (1.21613, 0.199803),
              (1.22022, 0.199199),
              (1.22432, 0.196859),
              (1.22841, 0.197471),
              (1.2325, 0.19799))
  # 300 data points from 0.009 to 1.2325, so separated by 0.004091973
  s1 = (ssqr - 0.009) / 0.004091973
  is1 = int(math.floor(s1))
  if is1 < 0:
    return best_data[0][1] # Below low-res limit for BEST data
  elif is1 >= 299:
    return best_data[0][299]  # Above high-res limit, about 0.9A
  else:
    ds = s1 - is1
    is2 = is1 + 1
    best_val = (1.-ds)*best_data[is1][1] + ds*best_data[is2][1]
    return best_val

def sphere_enclosing_model(model):
  sites_cart = model.get_sites_cart()
  cart_min = flex.double(sites_cart.min())
  cart_max = flex.double(sites_cart.max())
  model_midpoint = (cart_max + cart_min) / 2
  dsqrmax = flex.max((sites_cart - tuple(model_midpoint)).norms()) ** 2
  model_radius = math.sqrt(dsqrmax)
  return model_midpoint, model_radius

def flatten_model_region(mmm, d_min):
  # Flatten the region covered by the model
  # For map_manager, replace it by the mask-weighted mean within this region
  # For half-maps, replace by the mean and put back the original half-map
  # map difference to preserve the error signal.
  mm = mmm.map_manager()
  mm1 = mmm.map_manager_1()
  mm2 = mmm.map_manager_2()
  delta_mm = mm1.customized_copy(map_data = mm1.map_data() - mm2.map_data())
  model = mmm.model()
  mmm.create_mask_around_atoms(model=model, soft_mask=True, soft_mask_radius=d_min/2)
  working_mmm = mmm.deep_copy() # Save a copy before masking to work on later

  working_mmm.add_map_manager_by_id(delta_mm, map_id = 'delta_map')

  # Invert original mask and apply to original maps to flatten density under model
  mask_mm_inverse = mmm.get_map_manager_by_id('mask')
  mask_mm_inverse.set_map_data(map_data = 1. - mask_mm_inverse.map_data())
  mmm.apply_mask_to_maps(set_outside_to_mean_inside = False)

  # Apply original mask to working_mmm to get density and difference density under model
  mask_mm = working_mmm.get_map_manager_by_id('mask')
  working_mmm.apply_mask_to_maps(set_outside_to_mean_inside = False)

  # Get mean density for part of each map under model, to add back to
  # flattened region of each map
  mask_info = working_mmm.mask_info()
  weighted_points = mask_info.size*mask_info.mean
  wmm = working_mmm.map_manager()
  mean_map = flex.sum(wmm.map_data()) / weighted_points
  mask_data = mask_mm.map_data()
  mm.set_map_data(map_data = mm.map_data() + mean_map * mask_data)
  wmm1 = working_mmm.map_manager_1()
  mean_map1 = flex.sum(wmm1.map_data()) / weighted_points
  mm1.set_map_data(map_data = mm1.map_data() +
      mean_map1 * mask_data + delta_mm.map_data()/2)
  wmm2 = working_mmm.map_manager_2()
  mean_map2 = flex.sum(wmm2.map_data()) / weighted_points
  mm2.set_map_data(map_data = mm2.map_data() +
      mean_map2 * mask_data - delta_mm.map_data()/2)
  # mm.write_map('masked_map.map')
  # mm1.write_map('masked_half_map_1.map')
  # mm2.write_map('masked_half_map_2.map')
  mmm.remove_map_manager_by_id('mask')

def add_local_squared_deviation_map(
    mmm, radius, d_min, map_id_in='map_manager', map_id_out='variance_map'):
  """
  Add spherically-averaged squared map to map_model_manager

  Compulsory arguments:
  mmm:    map_model_manager containing input map
  radius: radius of sphere over which squared deviation is averaged
  d_min:  estimate of best resolution for map

  Optional arguments:
  map_id_in:  identifier of input map, defaults to map_manager
  map_id_out: identifier of output map, defaults to 'variance map'
  """

  mm_in = mmm.get_map_manager_by_id(map_id = map_id_in)
  coeffs_in = mm_in.map_as_fourier_coefficients(d_min=d_min)
  map_out = coeffs_in.local_standard_deviation_map(radius=radius, d_min=d_min)
  # map_out is an fft_map object, which can't easily be added to
  # map_model_manager as a similar map_manager object, so cycle through FT
  mm_out = map_out.as_map_manager()
  mean_square_map_coeffs = mm_out.map_as_fourier_coefficients(d_min=d_min)
  mmm.add_map_from_fourier_coefficients(mean_square_map_coeffs,
      map_id=map_id_out)
  # All map values should be positive, but round trip through FT might change
  # this. Check and add an offset if required to make minimum non-negative.
  mm_out = mmm.get_map_manager_by_id(map_id=map_id_out)
  min_map_value = flex.min(mm_out.map_data())
  if min_map_value < 0:
    offset = -min_map_value
    mm_out.set_map_data(map_data = mm_out.map_data() + offset)

def auto_sharpen_isotropic(mmm, d_min):
  working_mmm = mmm.deep_copy()

  mc1 = working_mmm.map_as_fourier_coefficients(d_min=d_min,
      map_id='map_manager_1')
  mc2 = working_mmm.map_as_fourier_coefficients(d_min=d_min,
      map_id='map_manager_2')

  nref = mc1.size()
  num_per_bin = 1000
  max_bins = 50
  min_bins = 6
  n_bins = int(round(max(min(nref / num_per_bin, max_bins), min_bins)))
  mc1.setup_binner(n_bins=n_bins)
  mc2.use_binner_of(mc1)

  sumw = 0
  sumwx = 0.
  sumwy = 0.
  sumwx2 = 0.
  sumwxy = 0.
  for i_bin in mc1.binner().range_used():
    sel = mc1.binner().selection(i_bin)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    mapCC = mc1sel.map_correlation(other=mc2sel)
    assert (mapCC < 1.) # Ensure these are really independent half-maps
    mapCC = max(mapCC,0.001) # Avoid zero or negative values
    FSCref = math.sqrt(2./(1.+1./mapCC))
    ssqr = mc1sel.d_star_sq().data()
    x = flex.mean_default(ssqr, 0) # Mean 1/d^2 for bin
    fsq = flex.pow2(flex.abs(mc1sel.data()))
    meanfsq = flex.mean_default(fsq, 0)
    y = math.log(meanfsq/(FSCref*FSCref))
    w = fsq.size()
    sumw += w
    sumwx += w * x
    sumwy += w * y
    sumwx2 += w * x**2
    sumwxy += w * x * y

  assert (nref == sumw) # Check no Fourier terms lost outside bins
  slope = (sumw * sumwxy - (sumwx * sumwy)) / (sumw * sumwx2 - sumwx**2)
  b_sharpen = 2 * slope
  all_ones = mc1.customized_copy(data = flex.double(mc1.size(), 1))
  b_terms_miller = all_ones.apply_debye_waller_factors(b_iso = b_sharpen)
  mc1 = mc1.customized_copy(data = mc1.data()*b_terms_miller.data())
  working_mmm.add_map_from_fourier_coefficients(mc1, map_id='map_manager_1')
  mc2 = mc2.customized_copy(data = mc2.data()*b_terms_miller.data())
  working_mmm.add_map_from_fourier_coefficients(mc2, map_id='map_manager_2')
  return working_mmm

def add_ordered_volume_mask(
    mmm, d_min, rad_factor=2, protein_mw=None, nucleic_mw=None,
    map_id_out='ordered_volume_mask'):
  """
  Add map defining mask covering the volume of most ordered density required
  to contain the specified content of protein and nucleic acid, judged by ratio
  of local map variance and local noise variance.

  Compulsory arguments:
  mmm: map_model_manager containing input half-maps in default map_managers
  d_min: estimate of best resolution for map

  Optional arguments:
  rad_factor: factor by which d_min is multiplied to get radius for averaging
    sphere, defaults to 2
  protein_mw*: molecular weight of protein expected in map, if any
  nucleic_mw*: molecular weight of nucleic acid expected in map
  map_id_out: identifier of output map, defaults to ordered_volume_mask

  * Note that at least one of protein_mw and nucleic_mw must be specified
  """

  assert((protein_mw is not None) or (nucleic_mw is not None))

  # Compute local average of squared density, using a sphere that will cover a
  # sufficient number of independent points and extending at least 5 A to
  # sample non-bonded contact distances. A rad_factor of 2 should yield
  # 4*Pi/3 * (2*2)^3 or about 270 independent points for the average; fewer
  # if the higher resolution data barely contribute. Larger values give less
  # noise but lower resolution for producing a mask. A minimum radius of 5
  # is enforced to explore next-nearest-neighbour density.
  radius = max(d_min*rad_factor, 5.)

  working_mmm = auto_sharpen_isotropic(mmm,d_min)
  wmm1 = working_mmm.map_manager_1()
  wmm2 = working_mmm.map_manager_2()
  # Input map_manager may contain arbitrarily filtered map that is not the
  # simple average of the two half-maps
  mm_mean = wmm1.customized_copy(map_data = (wmm1.map_data() + wmm2.map_data()) / 2)
  working_mmm.set_map_manager(mm_mean)
  delta_mm = wmm1.customized_copy(map_data = wmm1.map_data() - wmm2.map_data())
  working_mmm.add_map_manager_by_id(delta_mm, map_id = 'delta_map')

  add_local_squared_deviation_map(working_mmm, radius, d_min,
      map_id_in='map_manager', map_id_out='map_variance')
  mvmm = working_mmm.get_map_manager_by_id('map_variance')
  add_local_squared_deviation_map(working_mmm, radius, d_min,
      map_id_in='delta_map', map_id_out='noise_variance')
  # When maps are masked near the corners and edges, both signal and noise can
  # approach zero. Avoid getting close to dividing zero by zero, by adding a
  # small offset to the noise variance.
  # At the same time, correct noise variance by factor of two for averaging
  nvmm = working_mmm.get_map_manager_by_id('noise_variance')
  nvmm_min = min(flex.min(nvmm.map_data()) , 0.)
  nvmm_max = flex.max(nvmm.map_data())
  nvmm.set_map_data(map_data = (nvmm.map_data() - nvmm_min + nvmm_max/100.) / 2.)

  # Compute Z-score where values much greater than 1 indicate signal
  mm_Zscore = mvmm.customized_copy(
      map_data = flex.sqrt(mvmm.map_data()/nvmm.map_data()))

  # Choose enough points in averaged squared density map to covered expected
  # ordered structure. An alternative that could be implemented is to assign
  # any parts of map unlikely to arise from noise as ordered density, without
  # reference to expected content. This might require figuring out
  # the effective number of independent points in the averaging sphere to
  # calibrate the chi-square distribution.
  map_volume = working_mmm.map_manager().unit_cell().volume()
  Zscore_map_data = mm_Zscore.map_data()
  numpoints = Zscore_map_data.size()
  # Convert content into volume using partial specific volumes
  target_volume = 0.
  if protein_mw is not None:
    target_volume += protein_mw*1.229
  if nucleic_mw is not None:
    target_volume += nucleic_mw*0.945
  # Expand volume by amount needed for sphere expanded by d_min in radius
  equivalent_radius = math.pow(target_volume / (4*math.pi/3.),1./3)
  volume_factor = ((equivalent_radius+d_min)/equivalent_radius)**3
  expanded_target_volume = target_volume*volume_factor
  target_points = int(expanded_target_volume/map_volume * numpoints)

  # Find threshold for target number of masked points
  from cctbx.maptbx.segment_and_split_map import find_threshold_in_map
  threshold = find_threshold_in_map(target_points = target_points,
      map_data = Zscore_map_data)
  temp_bool_3D = (Zscore_map_data >=  threshold)
  # as_double method doesn't work for multidimensional flex.bool
  mask_shape = temp_bool_3D.all()
  overall_mask = temp_bool_3D.as_1d().as_double()
  overall_mask.reshape(flex.grid(mask_shape))
  new_mm = mmm.map_manager().customized_copy(map_data=overall_mask)
  mmm.add_map_manager_by_id(new_mm,map_id=map_id_out)

def get_grid_spacings(unit_cell, unit_cell_grid):
  assert unit_cell.parameters()[3:] == (90,90,90) # Required for this method
  sp = []
  for a,n in zip(unit_cell.parameters()[:3], unit_cell_grid):
    sp.append(a/n)
  return sp

def get_distance_from_center(c, unit_cell, unit_cell_grid = None,
    center = None):
  """
  Return a 3D flex array containing the distance of each grid point from the
  center of the map.
  Code provided by Tom Terwilliger
  """
  acc = c.accessor()
  if not unit_cell_grid: # Assume c contains complete unit cell
    unit_cell_grid = acc.all()
  nu,nv,nw = unit_cell_grid
  dx,dy,dz = get_grid_spacings(unit_cell, unit_cell_grid)
  if not center:
    center = (nu//2,nv//2,nw//2)

  # d is initially going to be squared distance from center
  d = flex.double(nu*nv*nw, 0)
  d.reshape(acc)

  # sum over x,y,z in slices
  dx2 = dx**2
  for i in range(nu):
    dist_sqr = dx2 * (i - center[0])**2
    d[i:i+1,0:nv,0:nw] += dist_sqr
  dy2 = dy**2
  for j in range(nv):
    dist_sqr = dy2 * (j - center[1])**2
    d[0:nu,j:j+1,0:nw] += dist_sqr
  dz2 = dz**2
  for k in range(nw):
    dist_sqr = dz2 * (k - center[2])**2
    d[0:nu,0:nv,k:k+1] += dist_sqr

  # Take square root to get distances
  d = flex.sqrt(d)

  return d

def get_maximal_mask_radius(mm_ordered_mask):
  unit_cell = mm_ordered_mask.unit_cell()
  om_data = mm_ordered_mask.map_data()
  d_from_c = get_distance_from_center(om_data, unit_cell = unit_cell)
  sel = om_data > 0
  selected_grid_indices = sel.iselection()
  mask_distances = d_from_c.select(selected_grid_indices)
  maximal_radius = mask_distances.min_max_mean().max
  return maximal_radius

def get_mask_radius(mm_ordered_mask,frac_coverage):
  """
  Get radius of sphere around map center enclosing desired fraction of
  ordered density
  """
  unit_cell = mm_ordered_mask.unit_cell()
  om_data = mm_ordered_mask.map_data()
  d_from_c = get_distance_from_center(om_data, unit_cell = unit_cell)
  sel = om_data > 0
  selected_grid_indices = sel.iselection()
  mask_distances = d_from_c.select(selected_grid_indices)
  mask_distances = mask_distances.select(flex.sort_permutation(data=mask_distances))
  masked_points = mask_distances.size()
  mask_radius = mask_distances[math.floor(frac_coverage*masked_points)-1]
  return mask_radius

def run_refine_cryoem_errors(
    mmm, d_min,
    map_1_id="map_manager_1", map_2_id="map_manager_2",
    ordered_mask_id='ordered_volume_mask',
    sphere_cent=None, radius=None, verbosity=1, prior_params=None,
    shift_map_origin=True, keep_full_map=False):
  """
  Refine error parameters from half-maps, make weighted map coeffs for region.

  Compulsory arguments:
  mmm: map_model_manager object containing two half-maps from reconstruction
  d_min: target resolution, either best resolution for map or resolution for
    target region

  Optional arguments:
  map_1_id: identifier of first half-map, if different from default of
    map_manager_1
  map_2_id: same for second half-map
  sphere_cent: center of sphere defining target region for analysis
    default is center of map
  radius: radius of sphere
    default (when sphere center not defined either) is 1/4 narrowest map width
  prior_params: refined parameters from previous call, usually from the
    whole reconstruction before focusing on a target region
  shift_map_origin: should map coefficients be shifted to correspond to
    original origin, rather than the origin being the corner of the box,
    default True
  keep_full_map: don't mask or box the input map
    default False, use with caution
  verbosity: 0/1/2/3/4 for mute/log/verbose/debug/testing
  """

  from scipy import interpolate
  from libtbx import group_args
  from iotbx.map_model_manager import map_model_manager

  if verbosity > 0:
    print("\nPrepare map for docking by analysing signal and errors")

  # Start from two half-maps and ordered volume mask in map_model_manager
  # Keep track of ordered volume in whole reconstruction
  ordered_mm = mmm.get_map_manager_by_id(map_id=ordered_mask_id)
  total_ordered_points = flex.mean(ordered_mm.map_data()) * ordered_mm.map_data().size()
  # Get map coefficients for maps after spherical masking
  ucpars = mmm.map_manager().unit_cell().parameters()
  if sphere_cent is None:
    # Default to sphere in center of cell extending halfway to nearest edge
    sphere_cent = flex.double((ucpars[0], ucpars[1], ucpars[2]))/2.
    if radius is None:
      radius = min(ucpars[0], ucpars[1], ucpars[2])/4.
  else:
    assert radius is not None
    sphere_cent = flex.double(sphere_cent)

  # Define box big enough to hold sphere plus soft masking
  boundary_to_smoothing_ratio = 2
  soft_mask_radius = d_min
  padding = soft_mask_radius * boundary_to_smoothing_ratio
  cushion = flex.double(3,radius+padding)
  cart_min = flex.double(sphere_cent) - cushion
  cart_max = flex.double(sphere_cent) + cushion
  for i in range(3): # Keep within unit cell
    cart_min[i] = max(cart_min[i],0)
    cart_max[i] = min(cart_max[i],ucpars[i])

  cs = mmm.crystal_symmetry()
  uc = cs.unit_cell()
  # working_mmm = map_model_manager()
  # Set some parameters that will be overwritten if masked and cut out
  over_sampling_factor = 1.
  weighted_points = mmm.map_data().size()
  box_volume = uc.volume()
  masked_volume = box_volume
  d_max = max(ucpars[0], ucpars[1], ucpars[2]) + d_min

  if keep_full_map:
    working_mmm = mmm.deep_copy()
  else:
    # Box the map within xyz bounds, converted to map grid units
    lower_frac = uc.fractionalize(tuple(cart_min))
    upper_frac = uc.fractionalize(tuple(cart_max))
    map_data = mmm.map_data()
    all_orig = map_data.all()
    lower_bounds = [int(math.floor(f * n)) for f, n in zip(lower_frac, all_orig)]
    upper_bounds = [int(math.ceil( f * n)) for f, n in zip(upper_frac, all_orig)]
    working_mmm = mmm.extract_all_maps_with_bounds(
        lower_bounds=lower_bounds, upper_bounds=upper_bounds)
    # Make and apply spherical mask
    working_mmm.create_spherical_mask(
      soft_mask_radius=soft_mask_radius,
      boundary_to_smoothing_ratio=boundary_to_smoothing_ratio)
    working_mmm.apply_mask_to_maps()
    mask_info = working_mmm.mask_info()

    # Keep track of weighted volume of map (in voxels) for determining relative
    # scale of sigmaE in different subvolumes: includes disordered volume
    weighted_points = mask_info.size*mask_info.mean # Weighted volume
    box_volume = uc.volume() * (
        working_mmm.map_data().size()/mmm.map_data().size())
    masked_volume = box_volume * mask_info.mean

    # Use weighted volume of ordered region of map to compute fraction of total
    # scattering contained in the working map
    cutout_ordered_mm = working_mmm.get_map_manager_by_id(map_id=ordered_mask_id)
    cutout_ordered_points = (flex.mean(cutout_ordered_mm.map_data()) *
          cutout_ordered_mm.map_data().size())
    fraction_scattering = cutout_ordered_points / total_ordered_points
    # Calculate an oversampling ratio defined as the ratio between the size of
    # the cut-out cube and the size of a cube that could contain a sphere
    # big enough to hold the volume of ordered density. Because the ratio of
    # the size of a cube to the sphere inscribed in it is about 1.9 (which is
    # close to the factor of two between typical protein volume and unit cell
    # volume in a crystal), this should yield likelihood scores on a similar
    # scale to crystallographic ones.
    over_sampling_factor = (cutout_ordered_mm.map_data().size() /
        (cutout_ordered_points * 6./math.pi) )

    d_max = 2*(radius+padding) + d_min # Size of sphere plus a bit

  mc1 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id=map_1_id)
  mc2 = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id=map_2_id)

  # Default binner may be preferable, but probably needs wider tests.
  # Could use bins of equal width in d_star_sq instead.
  # mc1.setup_binner_d_star_sq_bin_size()
  nref = mc1.size()
  num_per_bin = 1000
  max_bins = 50
  min_bins = 6
  n_bins = int(round(max(min(nref / num_per_bin, max_bins), min_bins)))
  mc1.setup_binner(n_bins=n_bins)
  mc2.use_binner_of(mc1)
  ssqmin = flex.min(mc1.d_star_sq().data())
  ssqmax = flex.max(mc1.d_star_sq().data())

  # Initialise parameters.  This requires slope and intercept of Wilson plot,
  # plus mapCC per bin.
  ssqr_bins = flex.double()
  target_spectrum = flex.double()
  meanfsq_bins = flex.double()
  mapCC_bins = flex.double()
  sumw = 0
  sumwx = 0.
  sumwy = 0.
  sumwx2 = 0.
  sumwxy = 0.
  for i_bin in mc1.binner().range_used():
    sel = mc1.binner().selection(i_bin)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    mapCC = mc1sel.map_correlation(other=mc2sel)
    assert (mapCC < 1.) # Ensure these are really independent half-maps
    mapCC = max(mapCC,0.001) # Avoid zero or negative values
    mapCC_bins.append(mapCC)
    ssqr = mc1sel.d_star_sq().data()
    x = flex.mean_default(ssqr, 0) # Mean 1/d^2 for bin
    ssqr_bins.append(x)  # Save for later
    fsq = flex.pow2(flex.abs(mc1sel.data()))
    meanfsq = flex.mean_default(fsq, 0)
    meanfsq_bins.append(meanfsq)
    y = math.log(meanfsq)
    w = fsq.size()
    sumw += w
    sumwx += w * x
    sumwy += w * y
    sumwx2 += w * x**2
    sumwxy += w * x * y
    target_power = default_target_spectrum(x) # Could have a different target
    target_spectrum.append(target_power)

  assert (nref == sumw) # Check no Fourier terms lost outside bins
  slope = (sumw * sumwxy - (sumwx * sumwy)) / (sumw * sumwx2 - sumwx**2)
  intercept = (sumwy - slope * sumwx) / sumw
  wilson_scale_intensity = math.exp(intercept)
  wilson_b_intensity = -4 * slope
  n_bins = ssqr_bins.size()

  if prior_params is not None:
    if d_min < 0.99*math.sqrt(1./prior_params['ssqmax']):
      print("Requested resolution is higher than prior parameters support")
      sys.stdout.flush()
      exit
    ssqr_prior = tuple(prior_params['ssqr_bins'])
    sigmaE_prior = tuple(prior_params['sigmaE_bins'])
    sEinterp = interpolate.interp1d(ssqr_prior,sigmaE_prior,fill_value="extrapolate")
    # Set sigmaE_scale to 1 after rescaling sigmaE_bins by volume comparison.
    # In principle just this could be refined instead of the bins.
    sigmaE_scale = 1.
    sigmaE_bins = flex.double(sEinterp(ssqr_bins))*(weighted_points/prior_params['weighted_points'])
    sigmaE_baniso = prior_params['sigmaE_baniso']
    sigmaE_beta = adptbx.u_star_as_beta(adptbx.u_cart_as_u_star(mc1.unit_cell(),adptbx.b_as_u(sigmaE_baniso)))
  else:
    sigmaE_scale = 1. # Fix at 1
    sigmaE_bins = []
    for i_bin in range(n_bins):
      sigmaE = meanfsq_bins[i_bin] * (1.-mapCC_bins[i_bin])
      sigmaE_bins.append(sigmaE)  # Error bin parameter
    sigmaE_beta = list(adptbx.u_iso_as_beta(mc1.unit_cell(), 0.))

  sigmaT_bins = [1.]*n_bins
  start_params = []
  start_params.append(wilson_scale_intensity/3.5) # Asqr_scale, factor out low-res BEST value
  start_params.extend(sigmaT_bins)
  wilson_u=adptbx.b_as_u(wilson_b_intensity)
  asqr_beta=list(adptbx.u_iso_as_beta(mc1.unit_cell(), wilson_u))
  start_params.extend(asqr_beta)
  start_params.append(sigmaE_scale)
  start_params.extend(sigmaE_bins)
  start_params.extend(sigmaE_beta)

  # create inputs for the minimizer's run method
  # Constrained refinement using prior parameters could be revisited later.
  # However, this would require all cryo-EM maps to obey the assumption
  # that half-maps are completely unmasked.
  if prior_params is not None:
    macro1 = ["Eprior"]
  else:
    macro1 = ["default"]        # protocol: fix error terms using prior
  macro2 = ["default"]       # protocol: refine sigmaE terms too
  protocol = [macro1, macro2]   # overall minimization protocol
  ncyc = 100                  # maximum number of microcycles per macrocycle
  minimizer_type = "bfgs"     # minimizer, bfgs or newton
  study_params = False        # flag for calling studyparams procedure
  output_level=verbosity      # 0/1/2/3/4 for mute/log/verbose/debug/testing

  # create instances of refine and minimizer
  refine_cryoem_errors = RefineCryoemErrors(
    mc1=mc1, mc2=mc2,
    ssqr_bins = ssqr_bins, target_spectrum = target_spectrum,
    start_x = start_params)
  minimizer = Minimizer(output_level=output_level)

  # Run minimizer
  minimizer.run(refine_cryoem_errors, protocol, ncyc, minimizer_type, study_params)
  refined_params=refine_cryoem_errors.x

  # Extract and report refined parameters
  i_par = 0
  asqr_scale = refined_params[i_par] # Not used for correction: leave map on original scale
  i_par += 1
  sigmaT_bins = refined_params[i_par:i_par + n_bins]
  i_par += n_bins
  asqr_beta = tuple(refined_params[i_par:i_par + 6])
  i_par += 6
  sigmaE_scale = refined_params[i_par] # Used here but not saved later
  i_par += 1
  sigmaE_bins = list(sigmaE_scale * flex.double(refined_params[i_par:i_par + n_bins]))
  i_par += n_bins
  sigmaE_beta = tuple(refined_params[i_par:i_par + 6])
  i_par += 6
  assert (i_par == len(refined_params))

  # Convert asqr_beta to a_beta for application in weights
  a_beta = tuple(flex.double(asqr_beta)/2)

  # Convert beta parameters to Baniso for (optional) use and output
  a_baniso = adptbx.u_as_b(adptbx.beta_as_u_cart(mc1.unit_cell(), a_beta))
  sigmaE_baniso = adptbx.u_as_b(adptbx.beta_as_u_cart(mc1.unit_cell(), sigmaE_beta))

  if verbosity > 0:
    print("\nRefinement of scales and error terms completed\n")
    print("\nParameters for A and BEST curve correction")
    print("  A overall scale: ",math.sqrt(asqr_scale))
    for i_bin in range(n_bins):
      print("  Bin #", i_bin + 1, "BEST curve correction: ", sigmaT_bins[i_bin])
    print("  A tensor as beta:", a_beta)
    print("  A tensor as Baniso: ", a_baniso)
    es = adptbx.eigensystem(a_baniso)
    a_beta_ev = es.vectors
    print("  Eigenvalues and eigenvectors:")
    for iv in range(3):
      print("  ",es.values()[iv],es.vectors(iv))

    print("\nParameters for SigmaE")
    if prior_params is not None:
      print("  SigmaE scale applied to prior bins:", sigmaE_scale)
    for i_bin in range(n_bins):
      print("  Bin #", i_bin + 1, "SigmaE base: ", sigmaE_bins[i_bin])
    print("  SigmaE tensor as beta:", sigmaE_beta)
    print("  SigmaE tensor as Baniso (intensity scale): ", sigmaE_baniso)
    es = adptbx.eigensystem(sigmaE_baniso)
    sigmaE_beta_ev = es.vectors
    print("  Eigenvalues and eigenvectors:")
    for iv in range(3):
      print("  ",es.values()[iv],es.vectors(iv))
    sys.stdout.flush()

  # Loop over bins to compute expectedE and Dobs for each Fourier term
  # Start with mean of half-map Fourier terms and make Miller array for Dobs
  expectE = mc1.customized_copy(data = (mc1.data() + mc2.data())/2)
  expectE.use_binner_of(mc1)
  dobs = expectE.customized_copy(data=flex.double(expectE.size(),0))
  i_bin_used = 0 # Keep track in case full range of bins not used
  weighted_map_noise = 0.
  if verbosity > 0:
    print("MapCC before and after rescaling as a function of resolution")
    print("Bin   <ssqr>   mapCC_before   mapCC_after")
  for i_bin in mc1.binner().range_used():
    sel = expectE.binner().selection(i_bin)
    eEsel = expectE.select(sel)

    # Make Miller array as basis for computing aniso corrections in this bin
    ones_array = flex.double(eEsel.size(), 1)
    all_ones = eEsel.customized_copy(data=ones_array)
    beta_a_miller = all_ones.apply_debye_waller_factors(
      u_star=adptbx.beta_as_u_star(a_beta))
    beta_sE_miller = all_ones.apply_debye_waller_factors(
      u_star=adptbx.beta_as_u_star(sigmaE_beta))

    # SigmaT is target_spectrum times sigmaT_bins correction factor
    sigmaT = sigmaT_bins[i_bin_used] * target_spectrum[i_bin_used]
    abeta_terms = beta_a_miller.data() # Anisotropy correction per reflection
    a2beta_terms = flex.pow2(abeta_terms)
    asqrSigmaT = asqr_scale * sigmaT * a2beta_terms
    sigmaE_terms = sigmaE_bins[i_bin_used] * beta_sE_miller.data()

    scale_terms = 1./flex.sqrt(asqrSigmaT + sigmaE_terms/2.)
    dobs_terms = 1./flex.sqrt(1. + sigmaE_terms/(2*asqrSigmaT))
    expectE.data().set_selected(sel, expectE.data().select(sel) * scale_terms)
    dobs.data().set_selected(sel, dobs_terms)
    weighted_map_noise += flex.sum(sigmaE_terms/(sigmaE_terms + 2*asqrSigmaT))

    # Apply corrections to mc1 and mc2 to compute mapCC after rescaling
    # SigmaE variance is twice as large for half-maps before averaging
    scale_terms_12 = 1./(abeta_terms + sigmaE_terms / (asqr_scale * sigmaT * abeta_terms))
    mc1.data().set_selected(sel, mc1.data().select(sel) * scale_terms_12)
    mc2.data().set_selected(sel, mc2.data().select(sel) * scale_terms_12)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    mapCC = mc1sel.map_correlation(other=mc2sel)
    if verbosity > 0:
      print(i_bin_used+1, ssqr_bins[i_bin_used], mapCC_bins[i_bin_used], mapCC)
      sys.stdout.flush()
    mapCC_bins[i_bin_used] = mapCC # Update for returned output
    i_bin_used += 1

  # At this point, weighted_map_noise is the sum of the noise variance for a
  # weighted half-map. In the following, this sum could be multiplied by two
  # for Friedel symmetry, but then divided by two for effects of averaging.
  weighted_map_noise = math.sqrt(weighted_map_noise) / masked_volume

  if verbosity > 0:
    print("Fraction of full map scattering: ",fraction_scattering)
    print("Over-sampling factor: ",over_sampling_factor)
    print("Weighted map noise: ",weighted_map_noise)

  # The following code could be used if we wanted to return a map_model_manager
  # wEmean = dobs*expectE
  # working_mmm.add_map_from_fourier_coefficients(map_coeffs=wEmean, map_id='map_manager_wtd')
  # working_mmm.write_map(map_id='map_manager_wtd',file_name='prepmap.map')
  # working_mmm.remove_map_manager_by_id(map_1_id)
  # working_mmm.remove_map_manager_by_id(map_2_id)

  shift_cart = working_mmm.shift_cart()
  if shift_map_origin:
    ucwork = expectE.crystal_symmetry().unit_cell()
    # shift_cart is position of original origin in boxed-map coordinate system
    # shift_frac should correspond to what has to be done to a model to put it
    # into the map, i.e. move it in the opposite direction
    shift_frac = ucwork.fractionalize(shift_cart)
    shift_frac = tuple(-flex.double(shift_frac))
    expectE = expectE.translational_shift(shift_frac)

  resultsdict = dict(
    n_bins = n_bins,
    ssqr_bins = ssqr_bins,
    ssqmin = ssqmin,
    ssqmax = ssqmax,
    weighted_points = weighted_points,
    asqr_scale = asqr_scale,
    sigmaT_bins = sigmaT_bins,
    asqr_beta = asqr_beta,
    a_baniso = a_baniso,
    sigmaE_bins = sigmaE_bins,
    sigmaE_baniso = sigmaE_baniso,
    mapCC_bins = mapCC_bins)
  return group_args(
    shift_cart = shift_cart,
    expectE = expectE, dobs = dobs,
    over_sampling_factor = over_sampling_factor,
    fraction_scattering = fraction_scattering,
    weighted_map_noise = weighted_map_noise,
    resultsdict = resultsdict)

# Command-line interface using argparse
def run():
  """
  Prepare cryo-EM map for docking by preparing weighted MTZ file.

  Obligatory command-line arguments (no keywords):
  half_map_1: name of file containing the first half-map from a reconstruction
  half_map_2: name of file containing the second half-map
  d_min: desired resolution, either best for whole map or for local region

  Optional command-line arguments (keyworded):
  --protein_mw*: molecular weight expected for protein component of ordered density
  --nucleic_mw*: same for nucleic acid component
  --model: PDB file for model that can either be used to flatten the map around
          the model, to search for the next component, or to cut out a sphere
          of density to refine and score the fit of the model
  --flatten_model: Use model to define region where map should be flattened
  --cutout_model: Use model to define sphere to process for refining and
          scoring the docking of this model
  --sphere_cent: Centre of sphere defining target map region (3 floats)
          defaults to centre of map
  --radius: radius of sphere (1 float)
          must be given if sphere_cent defined, otherwise
          defaults to narrowest extent of input map divided by 4
  --shift_map_origin: shift output mtz file to match input map on its origin?
          default True
  --file_root: root name for output files
  --write_params: write out refined parameters as a pickle file
  --read_params: start with refined parameters from earlier run
  --mute (or -m): mute output
  --verbose (or -v): verbose output
  --testing: extra verbose output for debugging
  * NB: At least one of protein_mw or nucleic_mw must be given
        Either a model or a sphere can be used to specify a cutout region, but
        not both
        The flatten_model option cannot be combined with cutting out a sphere.
  """
  import argparse
  import pickle
  from iotbx.map_model_manager import map_model_manager
  from iotbx.data_manager import DataManager
  dm = DataManager()
  dm.set_overwrite(True)

  parser = argparse.ArgumentParser(
          description='Prepare cryo-EM map for docking')
  parser.add_argument('map1',help='Map file for half-map 1')
  parser.add_argument('map2', help='Map file for half-map 2')
  parser.add_argument('d_min', help='d_min for maps', type=float)
  parser.add_argument('--protein_mw',
                      help='Molecular weight of protein component of map',
                      type=float)
  parser.add_argument('--nucleic_mw',
                      help='Molecular weight of nucleic acid component of map',
                      type=float)
  parser.add_argument('--model',help='Placed model')
  parser.add_argument('--flatten_model',help='Flatten map around model',
                      action='store_true')
  parser.add_argument('--cutout_model',help='Cut out sphere around model',
                      action='store_true')
  parser.add_argument('--sphere_cent',help='Centre of sphere for docking',
                      nargs=3, type=float)
  parser.add_argument('--radius',help='Radius of sphere for docking', type=float)
  parser.add_argument('--file_root',
                      help='Root of filenames for output')
  parser.add_argument('--read_params', help='Filename for prior parameters')
  parser.add_argument('--write_params', help='Write out refined parameters',
                      action='store_true')
  parser.add_argument('--shift_map_origin', dest='shift_map_origin', type=bool)
  parser.set_defaults(shift_map_origin=True)
  parser.add_argument('-m', '--mute', dest = 'mute',
                      help = 'Mute output', action = 'store_true')
  parser.add_argument('-v', '--verbose', dest = 'verbose',
                      help = 'Set output as verbose', action = 'store_true')
  parser.add_argument('--testing', dest = 'testing',
                      help='Set output as testing', action='store_true')
  args = parser.parse_args()
  d_min = args.d_min
  verbosity = 1
  if args.mute: verbosity = 0
  if args.verbose: verbosity = 2
  if args.testing: verbosity = 4
  shift_map_origin = args.shift_map_origin

  cutout_specified = False
  sphere_cent = None
  radius = None
  model = None

  protein_mw = None
  nucleic_mw = None
  if (args.protein_mw is None) and (args.nucleic_mw is None):
    print("At least one of protein_mw or nucleic_mw must be given")
    sys.stdout.flush()
    exit
  if args.protein_mw is not None:
    protein_mw = args.protein_mw
  if args.nucleic_mw is not None:
    nucleic_mw = args.nucleic_mw

  if args.model is not None:
    if not (args.cutout_model or args.flatten_model):
      print('Use for model must be specified (flatten or cut out map')
      sys.stdout.flush()
      exit
    model_file = args.model
    model = dm.get_model(model_file)

  if (args.sphere_cent is not None) and args.cutout_model:
    print("Only one method to define region to cut out (sphere or model) can be given")
    sys.stdout.flush()
    exit
  if args.sphere_cent is not None:
    assert args.radius is not None
    sphere_cent = tuple(args.sphere_cent)
    radius = args.radius
    cutout_specified = True
  if args.cutout_model:
    assert args.model is not None
    sphere_cent, radius = sphere_enclosing_model(model)
    radius = radius + d_min # Expand to allow width for density
    cutout_specified = True

  # Get prior parameters if provided
  if args.read_params is not None:
    infile = open(args.read_params,"rb")
    prior_params = pickle.load(infile)
    infile.close()
  else:
    prior_params = None

  # Create map_model_manager containing half-maps
  map1_filename = args.map1
  mm1 = dm.get_real_map(map1_filename)
  map2_filename = args.map2
  mm2 = dm.get_real_map(map2_filename)
  # delta_mm = mm1.customized_copy(map_data = mm1.map_data() - mm2.map_data())
  mmm = map_model_manager(model=model, map_manager_1=mm1, map_manager_2=mm2) #,
      # extra_map_manager_list=[delta_mm], extra_map_manager_id_list=['delta_map'])

  # Prepare maps by flattening model volume if desired
  if args.flatten_model:
    assert args.model is not None
    flatten_model_region(mmm, d_min)

  # Add mask map for ordered component of map.
  mask_id = 'ordered_volume_mask'
  add_ordered_volume_mask(mmm, d_min,
      protein_mw=protein_mw, nucleic_mw=nucleic_mw,
      map_id_out=mask_id)
  if verbosity>1:
    ordered_mm = mmm.get_map_manager_by_id(map_id=mask_id)
    if args.file_root is not None:
      map_file_name = args.file_root + "_ordered_volume_mask.map"
    else:
      map_file_name = "ordered_volume_mask.map"
    ordered_mm.write_map(map_file_name)

  if prior_params is None:
    # Initial refinement to get overall error parameters
    results = run_refine_cryoem_errors(mmm, d_min, verbosity=verbosity,
      shift_map_origin=shift_map_origin)
    prior_params = results.resultsdict
    if args.write_params:
      if args.file_root is not None:
        paramsfile = args.file_root + ".pickle"
      else:
        paramsfile = "prior_params.pickle"
      outf = open(paramsfile,"wb")
      pickle.dump(prior_params,outf,2)
      outf.close()

  # The following could loop over different regions
  if cutout_specified:
    # Refine to get scale and error parameters for docking region
    results = run_refine_cryoem_errors(mmm, d_min, verbosity=verbosity,
      sphere_cent=sphere_cent, radius=radius, prior_params=prior_params,
      shift_map_origin=shift_map_origin)

  expectE = results.expectE
  mtz_dataset = expectE.as_mtz_dataset(column_root_label='Emean')
  dobs = results.dobs
  mtz_dataset.add_miller_array(
      dobs,column_root_label='Dobs',column_types='W')
  mtz_object=mtz_dataset.mtz_object()

  if args.file_root is not None:
    mtzout_file_name = args.file_root + ".mtz"
  else:
    mtzout_file_name = "weighted_map_data.mtz"
  print ("Writing mtz for docking as",mtzout_file_name)
  if not shift_map_origin:
    shift_cart = results.shift_cart
    print ("Origin of full map relative to mtz:", shift_cart)
  dm.write_miller_array_file(mtz_object, filename=mtzout_file_name)
  over_sampling_factor = results.over_sampling_factor
  fraction_scattering = results.fraction_scattering
  print ("Over-sampling factor for Fourier terms:",over_sampling_factor)
  print ("Fraction of total scattering:",fraction_scattering)
  sys.stdout.flush()

if __name__ == "__main__":
  run()
