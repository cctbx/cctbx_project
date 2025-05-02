from __future__ import print_function
from __future__ import division

import math
from scitbx.array_family import flex
from scitbx.dtmin.minimizer import Minimizer
from scitbx.dtmin.refinebase import RefineBase
from scitbx.dtmin.reparams import Reparams
from scitbx.dtmin.bounds import Bounds
from cctbx import adptbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import sys

class RefineCryoemSignal(RefineBase):
  # Set up refinement class for dtmin minimiser (based on Phaser minimiser)
  def __init__(self, sumfsqr_lm, f1f2cos_lm, deltafsqr_lm,
               r_star, ssqr_bins, target_spectrum, start_x):
    RefineBase.__init__(self)
    # Prepare data that will be used repeatedly for fgh evaluation
    # Sample local mean with spacing derived from averaging radius
    # Assume box is at least close to being a cube
    self.unit_cell = sumfsqr_lm.unit_cell()
    recip_params = self.unit_cell.reciprocal_parameters()
    (astar, bstar, cstar) = recip_params[:3]
    spacing = int(round(r_star/(2*max(astar,bstar,cstar))))
    self.subsample = spacing**3
    h,k,l = sumfsqr_lm.indices().as_vec3_double().parts()
    ih = h.iround()
    ik = k.iround()
    il = l.iround()
    # Modulus of negative element in flex array is negative or zero, but we need
    # it to be positive. Subtract 1 so that lowest resolution reflection is 1,1,1
    ihmod = (ih % spacing + spacing - 1) % spacing
    ikmod = (ik % spacing + spacing - 1) % spacing
    ilmod = (il % spacing + spacing - 1) % spacing
    modsum = ihmod + ikmod + ilmod
    sel = (modsum == 0)
    self.sumfsqr_miller = sumfsqr_lm.select(sel)
    self.sumfsqr_miller.use_binning_of(sumfsqr_lm) # Matches subset
    self.f1f2cos_miller = f1f2cos_lm.select(sel)
    self.f1f2cos_miller.use_binner_of(self.sumfsqr_miller)
    self.sigmaE_miller = (deltafsqr_lm.select(sel)) / 2. # Half of deltafsqr power
    self.sigmaE_miller.use_binner_of(self.sumfsqr_miller)

    self.n_bins = self.sumfsqr_miller.binner().n_bins_used()  # Assume consistent binning
    self.ssqr_bins = ssqr_bins
    assert (self.n_bins == len(ssqr_bins))
    self.target_spectrum = target_spectrum
    self.start_x = start_x
    self.x = start_x[:]         # Full set of parameters
    d_min = math.sqrt(1./flex.max(sumfsqr_lm.d_star_sq().data()))
    cell_tensor = flex.double((astar*astar, bstar*bstar, cstar*cstar, astar*bstar, astar*cstar, bstar*cstar))

    # Maximum beta should usually give an exponential argument < 50 at resolution limit
    self.max_beta = cell_tensor * 45 * d_min**2
    # But make sure that maximum is higher than initial estimate
    # This can be an issue for subvolumes with unusually low signal
    asqr_beta = flex.double(self.x[self.n_bins + 1:self.n_bins + 7])
    if (self.max_beta[0] < 3 * asqr_beta[0]):
      rescale = 3 * asqr_beta[0]/self.max_beta[0]
      self.max_beta = rescale * self.max_beta
    self.max_beta = list(self.max_beta)

    # Choose maximum shift in beta to change by less than factor of 2 at resolution limit
    self.large_shifts_beta = list(cell_tensor*math.log(2.)*d_min**2)

    # Apply a relatively weak restraint to sphericity of beta
    self.sigmaSphericityBeta = cell_tensor * 2*d_min**2
    # sigmaSphericityBiso = adptbx.u_as_b(adptbx.beta_as_u_cart(self.unit_cell, tuple(self.sigmaSphericityBeta)))
    # print("sigmaSphericityBiso: ",tuple(sigmaSphericityBiso))


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
    subsample = self.subsample
    n_bins = self.n_bins
    i_par = 0
    asqr_scale = self.x[i_par]
    i_par += 1
    sigmaT_bins = self.x[i_par:i_par + n_bins]
    i_par += n_bins
    asqr_beta = tuple(self.x[i_par:i_par + 6])
    i_par += 6
    assert (i_par == len(self.x))

    # Initialise function, gradient and Hessian with zeros
    f = 0.
    g = flex.double(self.nmp, 0)
    h = flex.double(self.nmp * self.nmp, 0)
    h.reshape(flex.grid(self.nmp, self.nmp))

    # Loop over bins to accumulate target, gradient, Hessian
    i_bin_used = 0 # Keep track in case full range of bins not used
    for i_bin in self.sumfsqr_miller.binner().range_used():
      sel = self.sumfsqr_miller.binner().selection(i_bin)
      sumfsqr_miller_sel = self.sumfsqr_miller.select(sel)
      sumfsqr = sumfsqr_miller_sel.data()
      f1f2cos = self.f1f2cos_miller.data().select(sel)
      sigmaE_terms = self.sigmaE_miller.data().select(sel)

      # Make Miller array as basis for computing aniso corrections in bin
      # Let u = A^2*sigmaT to simplify computation of derivatives
      ones_array = flex.double(sumfsqr.size(), 1)
      all_ones = sumfsqr_miller_sel.customized_copy(data=ones_array)
      beta_miller_Asqr = all_ones.apply_debye_waller_factors(
        u_star=adptbx.beta_as_u_star(asqr_beta))
      u_terms = (asqr_scale * sigmaT_bins[i_bin_used]
        * self.target_spectrum[i_bin_used]) * beta_miller_Asqr.data()

      # Use local sphere fsc to compute relative weight of local vs global fitting
      # Steepness of sigmoid controlled by factor applied to fsc.
      # Mean value of f1f2cos should dominate sigmaS_terms calculation when signal
      # is good but it should be completely determined by the anisotropic error model
      # when there is little signal.
      # Keep some influence of anisotropic error model throughout for training.
      # wt_terms weight the anisotropy model in sigmaS calculation, with this
      # weight falling off with increasing local fsc.
      # Behaviour of wt as a function of fsc determined by steepness parameter for
      # the sigmoid function and the maximum fractional weight for the local
      # statistics. Values of local_weight near 1 seem to be better than lower
      # values.
      steep = 9.
      local_weight = 0.95
      fsc = f1f2cos/(sumfsqr/2.)
      # Make sure fsc is in range 0 to 1
      fsc = (fsc + flex.abs(fsc)) / 2
      fsc = (fsc + 1 - flex.abs(fsc - 1)) / 2
      wt_terms = flex.exp(steep*fsc)
      wt_terms = 1. - local_weight * (wt_terms/(math.exp(steep/2) + wt_terms))
      meanwt = flex.mean_default(wt_terms,0.)
      sigmaS_terms = wt_terms*u_terms + (1.-wt_terms)*f1f2cos
      # Make sure sigmaS is non-negative
      sigmaS_terms = (sigmaS_terms + flex.abs(sigmaS_terms))/2

      s2sigE = 2*sigmaS_terms + sigmaE_terms
      var_terms = s2sigE*sigmaE_terms

      # Leave out constant twologpi * number of reflections involved
      assert (flex.min(var_terms) > 0.)
      minusLL_terms = (sumfsqr * (sigmaS_terms + sigmaE_terms)
        - 2 * sigmaS_terms * f1f2cos) / var_terms + flex.log(var_terms)
      f += subsample*flex.sum(minusLL_terms)

      fgh_a_restraint = self.sphericity_restraint(asqr_beta,
          do_gradient, do_hessian, sigma_factor = 2.) # Looser for amplitude-squared
      f += meanwt * fgh_a_restraint[0]

      if do_gradient:
        s2sigE2 = flex.pow2(s2sigE)
        sumsqrcos = sumfsqr + 2*f1f2cos
        if (self.refine_Asqr_scale or self.refine_sigmaT_bins
              or self.refine_Asqr_beta):
          dmLL_by_dsigmaS_terms = (2 * s2sigE - sumsqrcos) / s2sigE2
          dmLL_by_du_terms = wt_terms * dmLL_by_dsigmaS_terms
        if self.refine_Asqr_beta:
          h_as_double, k_as_double, l_as_double = (
            sumfsqr_miller_sel.indices().as_vec3_double().parts())
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
          i_Asqr_scale = i_ref # Save for mixed second derivatives
          g[i_Asqr_scale] += subsample*flex.sum(dmLL_by_du_terms*du_by_dAsqr_scale)
          i_ref += 1
        i_par += 1
        if self.refine_sigmaT_bins: # Only affects U, just current bin
          du_by_dsigmaT_bin = u_terms / sigmaT_bins[i_bin_used]
          i_sigmaT_bin = i_ref+i_bin_used # Save for restraint terms below
          g[i_sigmaT_bin] += subsample*flex.sum(dmLL_by_du_terms*du_by_dsigmaT_bin)
          i_ref += self.n_bins
        i_par += self.n_bins
        if self.refine_Asqr_beta:  # Only affects U
          hh_factors = [-hh, -kk, -ll, -2*hk, -2*hl, -2*kl]
          du_by_dbetaA = []
          for i_beta in range(6):
            du_by_dbetaA.append(hh_factors[i_beta]*u_terms)
          for i_beta in range(6):
            g[i_ref+i_beta] += subsample*flex.sum(dmLL_by_du_terms * du_by_dbetaA[i_beta])
            g[i_ref+i_beta] += meanwt * fgh_a_restraint[1][i_beta] # Restraint term
          i_ref += 6
        i_par += 6

        assert (i_par == len(self.x))
        assert (i_ref == self.nmp)

        if do_hessian:
          s2sigE3 = s2sigE * s2sigE2
          if (self.refine_Asqr_scale or self.refine_sigmaT_bins
                or self.refine_Asqr_beta):

            d2mLL_by_dsigmaS2_terms = 4 * (
              sumsqrcos - 2 * sigmaS_terms - sigmaE_terms) / s2sigE3
            d2mLL_by_du2_terms = flex.pow2(wt_terms) * d2mLL_by_dsigmaS2_terms

          i_par = 0 # Keep track of index for unrefined parameters
          i_ref = 0  # Keep track of refined parameters
          # Note that second derivatives u wrt Asqr_scale and sigmaT_bins are 0
          if self.refine_Asqr_scale: # Only affects U
            h[i_ref,i_ref] += subsample*(flex.sum(d2mLL_by_du2_terms
              * flex.pow2(du_by_dAsqr_scale)))
            i_ref += 1
          i_par += 1
          if self.refine_sigmaT_bins: # Only affects U, current bin
            h[i_sigmaT_bin,i_sigmaT_bin] += subsample*flex.sum(
              d2mLL_by_du2_terms * flex.pow2(du_by_dsigmaT_bin))
            if self.refine_Asqr_scale:
              d2u_by_dAsqr_scale_by_dsigmaT_bin = du_by_dsigmaT_bin / asqr_scale
              cross_term = subsample*flex.sum(
                d2mLL_by_du2_terms * du_by_dAsqr_scale * du_by_dsigmaT_bin +
                dmLL_by_du_terms * d2u_by_dAsqr_scale_by_dsigmaT_bin)
              h[i_Asqr_scale,i_sigmaT_bin] += cross_term
              h[i_sigmaT_bin,i_Asqr_scale] += cross_term
            i_ref += self.n_bins
          i_par += self.n_bins
          if self.refine_Asqr_beta:  # Only affects U
            for i_beta in range(6):
              if self.refine_Asqr_scale: # Add Asqr_scale mixed derivatives
                d2u_by_dAsqr_scale_by_dbetaA = du_by_dbetaA[i_beta] / asqr_scale
                cross_term = subsample*flex.sum(
                  d2mLL_by_du2_terms * du_by_dAsqr_scale * du_by_dbetaA[i_beta] +
                  dmLL_by_du_terms * d2u_by_dAsqr_scale_by_dbetaA)
                h[i_Asqr_scale,i_ref+i_beta] += cross_term
                h[i_ref+i_beta,i_Asqr_scale] += cross_term
              for j_beta in range(6):
                h[i_ref+i_beta, i_ref+j_beta] += subsample*(
                  flex.sum(d2mLL_by_du2_terms * du_by_dbetaA[i_beta]*du_by_dbetaA[j_beta])
                  + flex.sum(dmLL_by_du_terms * hh_factors[i_beta]*hh_factors[j_beta]*u_terms) )
                # Also add restraint term
                h[i_ref+i_beta,i_ref+j_beta] += meanwt * fgh_a_restraint[2][i_beta, j_beta]
            i_ref += 6
          i_par += 6

          assert (i_par == len(self.x))
          assert (i_ref == self.nmp)

      # Restrain log of sigmaT_bins to 0, but downweighting low resolution
      d_bin = 1./math.sqrt(self.ssqr_bins[i_bin_used])
      sigmascale = 0.15 + 0.0005 * d_bin**3
      stbin = sigmaT_bins[i_bin_used]
      logbin = math.log(stbin)
      f += meanwt * (logbin/sigmascale)**2 / 2
      if do_gradient and self.refine_sigmaT_bins:
        g[i_sigmaT_bin] += meanwt * logbin / (stbin * sigmascale**2)
        if do_hessian:
          h[i_sigmaT_bin,i_sigmaT_bin] += meanwt * (1.-logbin)/(stbin*sigmascale)**2

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
      large_shifts.append(self.x[i_par]/4.)
    i_par += 1
    if self.refine_sigmaT_bins:
      for i_bin in range(self.n_bins):
        large_shifts.append(self.x[i_par+i_bin]/4.)
    i_par += self.n_bins
    if self.refine_Asqr_beta:
      large_shifts.extend(self.large_shifts_beta)
    i_par += 6
    assert (i_par == len(self.x))
    assert (len(large_shifts) == self.nmp)
    return large_shifts

  def set_macrocycle_protocol(self, macrocycle_protocol):
    # Possible parameters include overall scale of signal,
    # bin parameters for signal (BEST-like curve), anisotropy tensor for signal
    self.refine_mask = []  # Indicates "all" if left empty
    self.refine_Asqr_scale = True
    self.refine_sigmaT_bins = True
    self.refine_Asqr_beta = True

    # For each protocol, define variables that aren't refined
    # Currently only have default protocol to refine everything.
    if macrocycle_protocol == ["default"]:
      pass

    else:
      raise Sorry("Macrocycle protocol", macrocycle_protocol, " not recognised")

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

    if not full_list:
      assert (len(parameter_names) == self.nmp)
    else:
      assert (len(parameter_names) == len(self.x))
    return parameter_names

  def reparameterize(self):
    i_par = 0 # Keep track of index for unrefined parameters
    repar = []

    if self.refine_Asqr_scale:
      repar.append(Reparams(False))
    i_par += 1

    if self.refine_sigmaT_bins:
      repar.extend([Reparams(False) for i in range(self.n_bins)])
    i_par += self.n_bins

    if self.refine_Asqr_beta:
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
      this_bound.lower_on(0.0001*self.start_x[i_par])
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
      # this_bound.off()
      for i in range(6):
        this_bound.on(-self.max_beta[i],self.max_beta[i])
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
    pass

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
  model_centre = sites_cart.mean()
  dsqrmax = flex.max((sites_cart - tuple(model_centre)).norms()) ** 2
  model_radius = math.sqrt(dsqrmax)
  return model_centre, model_radius

def mask_fixed_model_region(mmm, d_min,
                                 fixed_model=None,
                                 ordered_mask_id=None,
                                 fixed_mask_id='mask_around_atoms'):
  # Flatten the region of the ordered volume mask covered by the fixed model
  # and keep the mask for the fixed model.

  # Do nothing if no fixed model
  if fixed_model is None:
    return

  radius = 3.
  if d_min is not None:
    radius = max(radius, d_min)
  mmm.create_mask_around_atoms(model=fixed_model, mask_atoms_atom_radius=radius, mask_id=fixed_mask_id)
  fixed_model_mask = mmm.get_map_manager_by_id(map_id=fixed_mask_id)

  if ordered_mask_id is not None:
    # Invert fixed_model_mask before flattening
    fixed_model_mask.set_map_data(map_data = 1. - fixed_model_mask.map_data())
    mmm.apply_mask_to_map(map_id=ordered_mask_id, mask_id=fixed_mask_id)
    # Invert back so we can use to assess fraction of fixed model in masked sphere
    fixed_model_mask.set_map_data(map_data = 1. - fixed_model_mask.map_data())

def add_local_squared_deviation_map(
    mmm, coeffs_in, radius, d_min, map_id_out):
  """
  Add spherically-averaged squared map to map_model_manager

  Compulsory arguments:
  mmm:        map_model_manager to add new map
  coeffs_in:  Miller array with coefficients for input map
  radius:     radius of sphere over which squared deviation is averaged
  d_min:      resolution to use for calculation
  map_id_out: identifier of output map
  """

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

def auto_sharpen_by_FSC(mc1, mc2):
  # Simply sharpen by bin-wise <F^2> and multiply by FSC
  # This didn't work as well as auto_sharpen_isotropic, but may be worth trying
  # again with a smooth curve over resolution instead of bins

  mapCC = mc1.map_correlation(other=mc2)
  if mapCC >= 1.:
    raise Sorry("Half-maps must be non-identical")
  mc1.setup_binner_d_star_sq_bin_size()
  mc2.use_binner_of(mc1)
  mc1s = mc1.customized_copy(data = mc1.data())
  mc2s = mc2.customized_copy(data = mc2.data())

  for i_bin in mc1.binner().range_used():
    sel = mc1.binner().selection(i_bin)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    mapCC = mc1sel.map_correlation(other=mc2sel)
    mapCC = max(mapCC,0.001) # Avoid zero or negative values
    FSCref = math.sqrt(2./(1.+1./mapCC))
    fsq = flex.pow2(flex.abs(mc1sel.data()))
    meanfsq = flex.mean_default(fsq, 1.e-10)
    mc1s.data().set_selected(sel, mc1.data().select(sel) * FSCref/math.sqrt(meanfsq))
    mc2s.data().set_selected(sel, mc2.data().select(sel) * FSCref/math.sqrt(meanfsq))

  return mc1s, mc2s

def auto_sharpen_isotropic(mc1, mc2):
  # Use Wilson plot in which <F^2> is divided by FSC^2, so that sharpening
  # downweights data with poor FSC. Results of fit to Wilson plot could be
  # odd if stated d_min goes much beyond real signal.

  mapCC = mc1.map_correlation(other=mc2)
  if mapCC >= 1.:
    raise Sorry("Half-maps must be non-identical")
  mc1.setup_binner_d_star_sq_bin_size()
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
    if mc1sel.size() < 2:
      continue
    mapCC = mc1sel.map_correlation(other=mc2sel)
    mapCC = max(mapCC,0.001) # Avoid zero or negative values
    FSCref = math.sqrt(2./(1.+1./mapCC))
    ssqr = mc1sel.d_star_sq().data()
    x = flex.mean_default(ssqr, 0) # Mean 1/d^2 for bin
    fsq = flex.pow2(flex.abs(mc1sel.data()))
    meanfsq = flex.mean_default(fsq, 0)
    y = math.log(meanfsq/(FSCref*FSCref))
    w = fsq.size() * FSCref
    sumw += w
    sumwx += w * x
    sumwy += w * y
    sumwx2 += w * x**2
    sumwxy += w * x * y

  slope = (sumw * sumwxy - (sumwx * sumwy)) / (sumw * sumwx2 - sumwx**2)
  b_sharpen = 2 * slope # Divided by 2 to apply to amplitudes
  all_ones = mc1.customized_copy(data = flex.double(mc1.size(), 1))
  b_terms_miller = all_ones.apply_debye_waller_factors(b_iso = b_sharpen)
  mc1s = mc1.customized_copy(data = mc1.data()*b_terms_miller.data())
  mc2s = mc2.customized_copy(data = mc2.data()*b_terms_miller.data())

  return mc1s, mc2s

def add_ordered_volume_mask(
    mmm, d_min, rad_factor=2, protein_mw=None, nucleic_mw=None,
    ordered_mask_id='ordered_volume_mask'):
  """
  Add map defining mask covering the volume of most ordered density required
  to contain the specified content of protein and nucleic acid, judged by
  local map variance.

  Compulsory arguments:
  mmm: map_model_manager containing input half-maps in default map_managers
  d_min: estimate of best resolution for map

  Optional arguments:
  rad_factor: factor by which d_min is multiplied to get radius for averaging
    sphere, defaults to 2
  protein_mw*: molecular weight of protein expected in map, if any
  nucleic_mw*: molecular weight of nucleic acid expected in map
  ordered_mask_id: identifier of output map, defaults to ordered_volume_mask

  * Note that at least one of protein_mw and nucleic_mw must be specified,
    and that the map must be complete
  """

  if (protein_mw is None) and (nucleic_mw is None):
    raise Sorry("At least one of protein_mw and nucleic_mw must be defined")
  if (not mmm.map_manager_1().is_full_size()) or  (
      not mmm.map_manager_2().is_full_size()):
    raise Sorry("Please supply half-maps that are full-size, not cut out from a larger box")

  if d_min is None or d_min <= 0:
    spacings = get_grid_spacings(mmm.map_manager().unit_cell(),
                                 mmm.map_manager().unit_cell_grid)
    d_min = 2.5 * max(spacings)

  # Compute local average of squared density, using a sphere that will cover a
  # sufficient number of independent points. A rad_factor of 2 should yield
  # 4*Pi/3 * (2*2)^3 or about 270 independent points for the average; fewer
  # if the higher resolution data barely contribute. Larger values give less
  # noise but lower resolution for producing a mask. A minimum radius of 5
  # is enforced to explore next-nearest-neighbour density. A maximum of 20
  # is enforced to keep some data at very low resolution
  radius = min(max(d_min*rad_factor, 5.),20.)

  d_work = (d_min + radius) # Save some time by lowering resolution
  mm1 = mmm.map_manager_1()
  mc1_in = mm1.map_as_fourier_coefficients(d_min=d_work)
  mm2 = mmm.map_manager_2()
  mc2_in = mm2.map_as_fourier_coefficients(d_min=d_work)
  if d_min < 10.: # Auto-sharpen unless very low resolution
    mc1s, mc2s = auto_sharpen_isotropic(mc1_in, mc2_in)
  else:
    mc1s = mc1_in
    mc2s = mc2_in
  mcs_mean  = mc1s.customized_copy(data = (mc1s.data() + mc2s.data())/2)

  add_local_squared_deviation_map(mmm, mcs_mean, radius, d_work,
      map_id_out='map_variance')
  mvmm = mmm.get_map_manager_by_id('map_variance')
  # mvmm.write_map("mvmm.map") # Uncomment to check intermediate result

  # Choose enough points in averaged squared density map to covered expected
  # ordered structure. An alternative that could be implemented is to assign
  # any parts of map unlikely to arise from noise as ordered density, without
  # reference to expected content, possibly like the false discovery rate approach.
  map_volume = mmm.map_manager().unit_cell().volume()
  mvmm_map_data = mvmm.map_data()
  numpoints = mvmm_map_data.size()
  # Convert content into volume using partial specific volumes
  target_volume = 0.
  if protein_mw is not None:
    target_volume += protein_mw*1.229
  if nucleic_mw is not None:
    target_volume += nucleic_mw*0.945
  # Expand volume by amount needed for sphere expanded by 1.5*d_min in radius
  equivalent_radius = math.pow(target_volume / (4*math.pi/3.),1./3)
  volume_factor = ((equivalent_radius+1.5*d_min)/equivalent_radius)**3
  expanded_target_volume = target_volume*volume_factor
  target_points = int(expanded_target_volume/map_volume * numpoints)

  # Find threshold for target number of masked points
  from cctbx.maptbx.segment_and_split_map import find_threshold_in_map
  threshold = find_threshold_in_map(target_points = target_points,
      map_data = mvmm_map_data)
  temp_bool_3D = (mvmm_map_data >=  threshold)
  # as_double method doesn't work for multidimensional flex.bool
  mask_shape = temp_bool_3D.all()
  overall_mask = temp_bool_3D.as_1d().as_double()
  overall_mask.reshape(flex.grid(mask_shape))
  new_mm = mmm.map_manager().customized_copy(map_data=overall_mask)

  # Clean up temporary map, then add ordered volume mask
  mmm.remove_map_manager_by_id('map_variance')
  mmm.add_map_manager_by_id(new_mm,map_id=ordered_mask_id)

def get_grid_spacings(unit_cell, unit_cell_grid):
  if unit_cell.parameters()[3:] != (90,90,90):
    raise Sorry("Unit cell must be orthogonal") # Required for this method
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
  if unit_cell.parameters()[3:] != (90,90,90):
    raise Sorry("Unit cell must be orthogonal") # Required for this method

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

  # Check assumption that this is a full map
  if not mm_ordered_mask.is_full_size():
    raise Sorry("Provided map or mask must be full-size")

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

  # Check assumption that this is a full map
  if not mm_ordered_mask.is_full_size():
    raise Sorry("Provided map or mask must be full-size")

  unit_cell = mm_ordered_mask.unit_cell()
  om_data = mm_ordered_mask.map_data()
  d_from_c = get_distance_from_center(om_data, unit_cell = unit_cell)
  sel = om_data > 0
  selected_grid_indices = sel.iselection()
  mask_distances = d_from_c.select(selected_grid_indices)
  mask_distances = mask_distances.select(flex.sort_permutation(data=mask_distances))
  masked_points = mask_distances.size()
  mask_radius = mask_distances[int(math.floor(frac_coverage*masked_points))-1]
  return mask_radius

def get_ordered_volume_exact(mm_ordered_mask,sphere_center,sphere_radius):
  """
  Get volume of density flagged as ordered inside sphere.
  Useful as a followup to fast spherical average approach, which can be deceived
  by counting ordered volume in a neighbouring cell if the cryoEM map is cropped
  tightly compared to the sphere radius.
  mm_ordered_mask: mask defining ordered volume
  sphere_center: centre of sphere in orthogonal coordinates
  sphere_radius: radius of sphere
  """

  # Check assumption that this is a full map
  if not mm_ordered_mask.is_full_size():
    raise Sorry("Provided map or mask must be full-size")

  unit_cell = mm_ordered_mask.unit_cell()
  om_data = mm_ordered_mask.map_data()
  sphere_center_frac = unit_cell.fractionalize(tuple(sphere_center))
  sphere_center_grid = [round(n * f) for n,f in zip(mm_ordered_mask.map_data().all(),
      sphere_center_frac)]
  d_from_c = get_distance_from_center(om_data, unit_cell = unit_cell,
      center = sphere_center_grid)
  sel = d_from_c <= sphere_radius
  selected_grid_indices = sel.iselection()
  om_data_sel = om_data.select(selected_grid_indices)
  sel = om_data_sel > 0
  ordered_in_sphere = om_data_sel.select(sel)
  ordered_volume = (ordered_in_sphere.size()/om_data.size()) * unit_cell.volume()
  return ordered_volume

def get_flex_max(a,b):
  c = a.deep_copy()
  sel = (b>a)
  c.set_selected(sel,b.select(sel))
  return c

def write_mtz(miller, file_name, root):
  mtz_dataset = miller.as_mtz_dataset(column_root_label=root)
  mtz_object=mtz_dataset.mtz_object()
  dm = DataManager()
  dm.set_overwrite(True)
  dm.write_miller_array_file(mtz_object, filename=file_name)

def largest_prime_factor(i):
  from libtbx.math_utils import prime_factors_of
  pf = prime_factors_of(i)
  if len(pf) == 0: # 0 or 1
    return 1
  else:
    return pf[-1]

def next_allowed_grid_size(i, largest_prime=5):
  if i<2:
    j = 2
  elif i%2 == 1:
    j = i+1
  else:
    j = i
  while (largest_prime_factor(j) > largest_prime):
    j += 2
  return j

def get_sharpening_b(intensities):

  # Initialise data for Wilson plot
  sumw = 0.
  sumwx = 0.
  sumwy = 0.
  sumwxy = 0.
  sumwx2 = 0.

  # Assume that binning has been defined and get data in bins
  int_binned = intensities.customized_copy(data = intensities.data())
  int_binned.setup_binner_d_star_sq_bin_size()
  for i_bin in int_binned.binner().range_used():
    sel = int_binned.binner().selection(i_bin)
    int_sel = int_binned.select(sel)
    ssqr = int_sel.d_star_sq().data()
    x = flex.mean_default(ssqr, 0) # Mean 1/d^2 for bin
    mean_int_sel_data = flex.mean_default(int_sel.data(),0)
    if mean_int_sel_data > 0.:
      y = math.log(mean_int_sel_data)
      w = int_sel.size()
      sumw += w
      sumwx += w * x
      sumwy += w * y
      sumwxy += w * x * y
      sumwx2 += w * x * x

  slope = (sumw * sumwxy - (sumwx * sumwy)) / (sumw * sumwx2 - sumwx**2)
  b_sharpen = 4 * slope

  return b_sharpen

def intensities_as_expanded_map(mm,marray):
  '''
  Take miller array (assumed here to be real-valued and in P1, but could be
  generalised), place values in a map, offset to the center.
  '''
  h_indices, k_indices, l_indices = marray.indices().as_vec3_double().parts()
  h_indices = h_indices.iround()
  k_indices = k_indices.iround()
  l_indices = l_indices.iround()
  hmin = flex.min(h_indices)
  hmax = flex.max(h_indices)
  kmin = flex.min(k_indices)
  kmax = flex.max(k_indices)
  lmin = flex.min(l_indices)
  lmax = flex.max(l_indices)
  assert ((hmin < 0) and (hmax >= -hmin))
  assert ((kmin < 0) and (kmax >= -kmin))
  assert (lmin == 0)
  data = marray.data()

  # Shift hkl to center, add buffer on edges to get to sensible grid for FFT
  # Apply Friedel symmetry to get full sphere of Fourier coefficients
  h_range = next_allowed_grid_size(2*hmax+2)
  h_offset = h_range//2
  h_map_indices = h_offset + h_indices
  h_inverse_indices = h_offset - h_indices

  k_range = next_allowed_grid_size(2*kmax+2)
  k_offset = k_range//2
  k_map_indices = k_offset + k_indices
  k_inverse_indices = k_offset - k_indices

  l_range = next_allowed_grid_size(2*lmax+2)
  l_offset = l_range//2
  l_map_indices = l_offset + l_indices
  l_inverse_indices = l_offset - l_indices

  miller_as_map_data = flex.double(h_range*k_range*l_range,0.)
  miller_as_map_data.reshape(flex.grid(h_range,k_range,l_range))
  for ih,ik,il,coeff in zip(h_map_indices,k_map_indices,l_map_indices,data):
    miller_as_map_data[ih,ik,il] = coeff
  for ih,ik,il,coeff in zip(h_inverse_indices,k_inverse_indices,l_inverse_indices,data):
    miller_as_map_data[ih,ik,il] = coeff

  # Unit cell for reciprocal space sphere of data is obtained from the reciprocal cell
  # parameters for the coefficients computed from the boxed map, times the
  # extents in h,k,l
  from cctbx import crystal
  abcstar = mm.crystal_symmetry().unit_cell().reciprocal_parameters()[:3]
  ucr = tuple(list(flex.double(abcstar)*flex.double((h_range,k_range,l_range)))+[90,90,90])
  crystal_symmetry=crystal.symmetry(ucr,1)
  from iotbx.map_manager import map_manager
  mm_data = map_manager(map_data=miller_as_map_data,unit_cell_grid=miller_as_map_data.all(),
      unit_cell_crystal_symmetry=crystal_symmetry,wrapping=False)
  return group_args(mm_data = mm_data,
                    h_map_indices = h_map_indices,
                    k_map_indices = k_map_indices,
                    l_map_indices = l_map_indices)

def expanded_map_as_intensities(rmm, marray, h_map_indices, k_map_indices, l_map_indices):
  '''
  Fetch miller array data back from map. Only take unique values, ignoring
  Friedel mates.
  '''
  rmap_data = rmm.map_data()
  miller_data = flex.double()
  for ih,ik,il in zip(h_map_indices,k_map_indices,l_map_indices):
    miller_data.append(rmap_data[ih,ik,il])
  miller_array = marray.customized_copy(data = miller_data)
  return miller_array

def local_mean_intensities(mm, d_min, intensities, r_star, b_sharpen):
  """
  Compute local means of input intensities (or amplitudes) using a convolution
  followed optionally by a resolution-dependent renormalisation

  Compulsory argument:
  mm: map_manager corresponding to system from which intensities are obtained
  d_min: target resolution of locally-averaged intensities
  intensities: Terms to be locally averaged, extended by r_star from 1/d_min
  r_star: radius in reciprocal space for averaging
  """

  # Check that resolution range has been extended from d_min
  d_star_sq_max = flex.max(intensities.d_star_sq().data())
  assert d_star_sq_max > 1/d_min**2

  # Sharpen intensities to reduce dynamic range for numerical stability
  # Put the overall B back at end
  int_sharp = intensities.customized_copy(data = intensities.data())
  all_ones = int_sharp.customized_copy(data = flex.double(int_sharp.size(), 1))
  b_terms_miller = all_ones.apply_debye_waller_factors(b_iso = b_sharpen)
  int_sharp = int_sharp.customized_copy(data = intensities.data()*b_terms_miller.data())

  # Turn intensity values into a spherical map inside a cube
  results = intensities_as_expanded_map(mm,int_sharp)
  coeffs_as_map = results.mm_data # map_manager with intensities
  h_map_indices = results.h_map_indices # grid positions for original hkl in order
  k_map_indices = results.k_map_indices
  l_map_indices = results.l_map_indices

  # Put map into map_model_manager to ensure matching grid for spherically-averaged map
  rmmm = map_model_manager(map_manager = coeffs_as_map)

  # Cell for "map" of intensities is reciprocal of real-space cell
  # Compute reciprocal d_min as twice the grid spacing between intensities
  rcp = coeffs_as_map.unit_cell().parameters()
  nu, nv, nw = coeffs_as_map.map_data().all()
  rdmin = 2*flex.min(flex.double(rcp[:3])/flex.double((nu,nv,nw)))
  inverse_intensities = coeffs_as_map.map_as_fourier_coefficients(d_min=rdmin)

  # Compute G-function, avoiding divide by zero and numerical precision issues
  # for the origin term
  stol = flex.sqrt(inverse_intensities.sin_theta_over_lambda_sq().data())
  w = 4 * stol * math.pi * r_star
  sel = (w < 0.001)
  w.set_selected(sel,0.001)
  sphere_reciprocal = 3 * (flex.sin(w) - w * flex.cos(w))/flex.pow(w, 3)
  sphere_reciprocal.set_selected(sel,1.)

  # Spherically-averaged intensities from FT of product of FTs
  prod_coeffs = inverse_intensities.customized_copy(
      data = inverse_intensities.data()*sphere_reciprocal)
  rmmm.add_map_from_fourier_coefficients(prod_coeffs, map_id='local_mean_map')
  local_mean_as_map = rmmm.get_map_manager_by_id(map_id='local_mean_map')

  # Associate map values with hkl using saved indices, then restrict to dmin
  extended_mean = expanded_map_as_intensities(local_mean_as_map,int_sharp,
      h_map_indices,k_map_indices,l_map_indices)

  min_in = flex.min(int_sharp.data())
  if min_in >= 0: # Make sure strictly non-negative remains non-negative
    extended_mean = extended_mean.customized_copy(data =
      (extended_mean.data()+min_in + flex.abs(extended_mean.data()-min_in))/2 )

  # Select data to desired resolution, remove sharpening from above
  local_mean = extended_mean.select(extended_mean.d_spacings().data() >= d_min)
  all_ones = local_mean.customized_copy(data = flex.double(local_mean.size(), 1))
  b_terms_miller = all_ones.apply_debye_waller_factors(b_iso = -b_sharpen)
  local_mean = local_mean.customized_copy(data = local_mean.data()*b_terms_miller.data())

  return local_mean

def local_mean_density(mm, radius):
  """
  Compute spherically-averaged map

  Compulsory argument:
  mm: map_manager containing map to be averaged
  radius: radius of sphere for averaging
  """

  # Put map into map_model_manager to ensure matching grid for spherically-averaged map
  mmm = map_model_manager(map_manager = mm)
  mc = mm.map_as_fourier_coefficients(d_min=radius/5.)

  # Compute G-function, avoiding divide by zero and numerical precision issues
  # for the origin term
  stol = flex.sqrt(mc.sin_theta_over_lambda_sq().data())
  w = 4 * stol * math.pi * radius
  sel = (w == 0.)
  w.set_selected(sel,0.0001)
  sphere_reciprocal = 3 * (flex.sin(w) - w * flex.cos(w))/flex.pow(w, 3)
  sphere_reciprocal.set_selected(sel,1.)

  # Spherically-averaged map values from FT of product of FTs
  prod_coeffs = mc.customized_copy(data = mc.data()*sphere_reciprocal)
  mmm.add_map_from_fourier_coefficients(prod_coeffs, map_id='local_mean_map')
  local_mean_mm = mmm.get_map_manager_by_id(map_id='local_mean_map')

  return local_mean_mm

def sigmaA_from_model_in_map(expectE, dobs, Fcalc, over_sampling_factor,
                             verbosity=0, log=sys.stdout):
  """
  Compute sigmaA values needed to explain agreement of model with map given map errors

  Compulsory arguments:
  expectE: Emean giving expected E for map
  dobs:    weighting factor for map coefficient accuracy
  Fcalc:   calculated structure factors for placed model
  over_sampling_factor: account for oversampling in LLG calculation
                        relevant here for precision of sigmaA estimate
  Optional arguments:
  verbosity: levels greater than default of 0 progressively more verbose
  log:       destination for printing if any
  """
  model_sigmaA = expectE.customized_copy(data=flex.double(expectE.size(),0))
  Ecalc = Fcalc.customized_copy(data=Fcalc.data())
  xdat = []
  ydat = []
  wdat = []
  if expectE.binner() is None:
    expectE.setup_binner_d_star_sq_bin_size()
  if verbosity > 1:
    print("SigmaA curve for placed model contribution",file=log)
    print(" # of terms   mean(ssqr)  sigmaA   sigma(sigmaA)")
  for i_bin in expectE.binner().range_used():
    sel = expectE.binner().selection(i_bin)
    eEsel = expectE.select(sel)
    abseE = eEsel.amplitudes().data()
    p1 = eEsel.phases().data()
    Fcalc_sel = Fcalc.select(sel)
    absFc = Fcalc_sel.amplitudes().data()
    sigmaP = flex.mean_default(flex.pow2(absFc),0.)
    absEc = absFc/math.sqrt(sigmaP)
    Ecalc.data().set_selected(sel,
                  Fcalc.data().select(sel) / math.sqrt(sigmaP))
    p2 = Fcalc_sel.phases().data()
    cosdphi = flex.cos(p2 - p1)
    dobssel = dobs.data().select(sel)

    # Solve for sigmaA yielding 1st derivative of LLG approximately equal to zero (assuming only numerator matters)
    sum0 = flex.sum(2.*dobssel*abseE*absEc*cosdphi)
    sum1 = flex.sum(2*flex.pow2(dobssel) * (1. - flex.pow2(abseE) - flex.pow2(absEc)))
    sum2 = flex.sum(2*flex.pow(dobssel,3)*abseE*absEc*cosdphi)
    sum3 = flex.sum(- 2*flex.pow(dobssel,4))
    u1 = -2*sum2**3 + 9*sum1*sum2*sum3 - 27*sum0*sum3**2
    sqrt_arg = u1**2 - 4*(sum2**2 - 3*sum1*sum3)**3
    if sqrt_arg < 0: # Paranoia: haven't run into this in a wide variety of calculations
      raise Sorry("Argument of square root in sigmaA calculation is negative")
    third = 1./3
    x1 = (u1 + math.sqrt(sqrt_arg))**third
    sigmaA = ( (2 * 2.**third * sum2**2 - 6 * 2.**third*sum1*sum3
                - 2*sum2*x1 + 2.**(2*third) * x1**2) /
              (6 * sum3 * x1) )
    sigmaA = max(min(sigmaA,0.999),1.E-6)
    # Now work out approximate e.s.d. of sigmaA estimate from 2nd derivative of LLG at optimal sigmaA
    denom = flex.pow((1-flex.pow2(dobssel)*sigmaA**2),3)
    sum0 = flex.sum((2*flex.pow2(dobssel) * (1. - flex.pow2(abseE) - flex.pow2(absEc)))/denom)
    sum1 = flex.sum((12*flex.pow(dobssel,3)*abseE*absEc*cosdphi)/denom)
    sum2 = flex.sum((-6*flex.pow(dobssel,4)*(flex.pow2(abseE) + flex.pow2(absEc)))/denom)
    sum3 = flex.sum((4*flex.pow(dobssel,5)*abseE*absEc*cosdphi)/denom)
    sum4 = flex.sum((-2*flex.pow(dobssel,6))/denom)
    d2LLGbydsigmaA = sum0 + sum1*sigmaA + sum2*sigmaA**2 + sum3*sigmaA**3 + sum4*sigmaA**4
    d2LLGbydsigmaA /= over_sampling_factor
    sigma_sigmaA = 1./math.sqrt(abs(d2LLGbydsigmaA))
    ssqr = flex.mean_default(eEsel.d_star_sq().data(),0)
    ndat = dobssel.size()
    if verbosity > 1:
      print(ndat,ssqr,sigmaA,sigma_sigmaA,file=log)
    xdat.append(ssqr)
    ydat.append(math.log(sigmaA))
    wdat.append(abs(d2LLGbydsigmaA)) # Inverse variance estimate

  xdat = flex.double(xdat)
  ydat = flex.double(ydat)
  wdat = flex.double(wdat)
  W = flex.sum(wdat)
  Wx = flex.sum(wdat * xdat)
  Wy = flex.sum(wdat * ydat)
  Wxx = flex.sum(wdat * xdat * xdat)
  Wxy = flex.sum(wdat * xdat * ydat)
  slope = (W * Wxy - Wx * Wy) / (W * Wxx - Wx * Wx)
  intercept = (Wy * Wxx - Wx * Wxy) / (W * Wxx - Wx * Wx)
  frac_complete = math.exp(2.*intercept)
  if verbosity > 0:
    print("Slope and intercept of sigmaA plot for placed model contribution: ",slope,intercept,file=log)
    print("Approximate fraction of scattering: ",frac_complete,file=log)
  linlogsiga = flex.double()
  logsiga_combined = flex.double()
  sigma_linlog = math.log(1.5) # Allow 50% error in linear fit
  i_bin_used = 0
  for i_bin in expectE.binner().range_used():
    sigmaA = math.exp(ydat[i_bin_used])
    # For complete model, slope should be negative, but can be positive if the
    # modelled part is much better ordered than the larger unmodelled part
    linlog = min(math.log(0.999),(intercept + slope*xdat[i_bin_used]))
    linlogsiga.append(linlog)
    sigma_sigmaA = 1./math.sqrt(wdat[i_bin_used])
    sigma_lnsigmaA = math.exp(-ydat[i_bin_used])*sigma_sigmaA
    combined_logsiga = ((ydat[i_bin_used]/sigma_lnsigmaA**2 + linlog/sigma_linlog**2) /
                        (1./sigma_lnsigmaA**2 + 1./sigma_linlog**2) )
    logsiga_combined.append(combined_logsiga)
    combined_siga = math.exp(combined_logsiga)

    sel = expectE.binner().selection(i_bin)
    model_sigmaA.data().set_selected(sel, combined_siga)
    i_bin_used += 1

  if verbosity > 1:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,1)
    ax[0].plot(xdat,ydat,label="ln(sigmaA)")
    ax[0].plot(xdat,linlogsiga,label="line fit")
    ax[0].plot(xdat,logsiga_combined,label="combined")
    ax[0].legend(loc="upper right")
    ax[1].plot(xdat,flex.exp(ydat),label="sigmaA")
    ax[1].plot(xdat,flex.exp(linlogsiga),label="line fit")
    ax[1].plot(xdat,flex.exp(logsiga_combined),label="combined")
    ax[1].legend(loc="lower left")
    plt.savefig("fixed_model_sigmaA_plot.png")

  return(group_args(model_sigmaA=model_sigmaA, Ecalc=Ecalc))

def assess_halfmap_errors(
          working_mmm, d_min, sphere_points=500, verbosity=1, log=sys.stdout):
  """
  Refine error parameters from half-maps, make weighted map coeffs for region.

  Compulsory arguments:
  working_mmm: map_model_manager object containing two half-maps from reconstruction
  d_min: target resolution, either best resolution for map or resolution for
    target region. Determined automatically if None.

  Optional arguments:
  sphere_points: number of reflections to use for local spherical average
  verbosity: 0/1/2/3/4 for mute/log/verbose/debug/testing
  log: destination for standard output
  """

  # Compute local averages needed for initial covariance terms
  # Do these calculations at extended resolution to avoid edge effects up to
  # desired resolution
  wuc = working_mmm.crystal_symmetry().unit_cell()
  unit_cell_grid = working_mmm.map_manager().unit_cell_grid
  ucpars = wuc.parameters()
  d_max = max(ucpars[0], ucpars[1], ucpars[2])
  work_mm = working_mmm.map_manager()
  v_star = 1./wuc.volume()
  r_star = math.pow(3*sphere_points*v_star/(4*math.pi),1./3.)
  spacings = get_grid_spacings(wuc,unit_cell_grid)

  # Check and potentially re-evaluate d_min
  guess_d_min = False
  minimum_possible_d_min = 2.*max(spacings)
  if d_min is None or d_min <= 0.:
    # Set first guess of d_min if no value provided
    # The majority of cryo-EM deposits have a ratio of nominal resolution to pixel
    # size between 2.5 and 3.5.
    guess_d_min = True
    d_min = 2.5 * max(spacings)
  elif d_min < minimum_possible_d_min:
    d_min = minimum_possible_d_min
    guess_d_min = True

  # Compute map coefficients from half-maps, then either proceed or evaluate
  # sensible d_min from different choices
  d_min_extended = 1./(1./d_min + r_star)
  mc1 = working_mmm.map_manager_1().map_as_fourier_coefficients(d_min=d_min_extended, d_max=d_max)
  mc2 = working_mmm.map_manager_2().map_as_fourier_coefficients(d_min=d_min_extended, d_max=d_max)

  while guess_d_min:
    # Check whether signal goes to higher resolution than initial guess
    d_min_from_fsc = mc1.d_min_from_fsc(other=mc2, fsc_cutoff=0.05).d_min
    # d_min_from_fsc returns None if the fsc_cutoff has not been passed by resolution limit
    if d_min_from_fsc is not None:
      d_min = max(d_min_from_fsc,minimum_possible_d_min)
      guess_d_min = False
    else:
      if d_min > minimum_possible_d_min: # Set d_min to highest possible value and repeat
        d_min = minimum_possible_d_min
      else:
        guess_d_min = False # d_min is already highest possible value
    d_min_extended = 1./(1./d_min + r_star)
    mc1 = working_mmm.map_manager_1().map_as_fourier_coefficients(d_min=d_min_extended, d_max=d_max)
    mc2 = working_mmm.map_manager_2().map_as_fourier_coefficients(d_min=d_min_extended, d_max=d_max)

  f1 = flex.abs(mc1.data())
  f2 = flex.abs(mc2.data())
  p1 = mc1.phases().data()
  p2 = mc2.phases().data()
  sumfsqr = mc1.customized_copy(data = flex.pow2(f1) + flex.pow2(f2))
  f1f2cos = mc1.customized_copy(data = f1 * f2 * flex.cos(p2 - p1))
  deltafsqr = mc1.customized_copy(data = flex.pow2(flex.abs(mc1.data()-mc2.data())))

  # Local mean calculations of covariance terms use data to extended resolution
  # but return results to desired resolution
  # Calculating sumfsqr_local_mean from f1f2cos and deltafsqr local means turns out to have
  # improved numerical stability, i.e. don't end up with negative or zero deltafsqr_local_mean
  b_sharpen = get_sharpening_b(sumfsqr)
  f1f2cos_local_mean = local_mean_intensities(work_mm, d_min, f1f2cos,
                                              r_star, b_sharpen)
  deltafsqr_local_mean = local_mean_intensities(work_mm, d_min, deltafsqr,
                                                r_star, b_sharpen)
  sumfsqr_local_mean = f1f2cos_local_mean.customized_copy(data = 2*f1f2cos_local_mean.data() + deltafsqr_local_mean.data())

  # Trim starting sumfsqr back to desired resolution limit, setup binning
  # NB: f1f2cos and deltafsqr aren't used any more
  mc1 = mc1.select(mc1.d_spacings().data() >= d_min)
  assert (mc1.size() == sumfsqr_local_mean.size())
  mc2 = mc2.select(mc2.d_spacings().data() >= d_min)
  sumfsqr = sumfsqr.select(sumfsqr.d_spacings().data() >= d_min)
  mc1.setup_binner_d_star_sq_bin_size()
  mc2.use_binner_of(mc1)
  sumfsqr.use_binner_of(mc1)
  sumfsqr_local_mean.use_binner_of(mc1)
  f1f2cos_local_mean.use_binner_of(mc1)

  # Prepare starting parameters for signal refinement
  mapCC = mc1.map_correlation(other=mc2)
  if mapCC >= 1.:
    raise Sorry("Half-maps must be non-identical")
  mapCC_bins = flex.double()
  ssqr_bins = flex.double()
  target_spectrum = flex.double()
  sumw = 0
  sumwx = 0.
  sumwy = 0.
  sumwx2 = 0.
  sumwxy = 0.
  nkept = 0
  for i_bin in sumfsqr.binner().range_used():
    sel = sumfsqr.binner().selection(i_bin)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    if mc1sel.size() == 0:
      continue
    mapCC = mc1sel.map_correlation(other=mc2sel)
    mapCC_bins.append(mapCC) # Store for before/after comparison
    sumfsqsel = sumfsqr.select(sel)
    meanfsq = flex.mean_default(sumfsqsel.data(),0.) / 2
    ssqr = mc1sel.d_star_sq().data()
    x = flex.mean_default(ssqr, 0) # Mean 1/d^2 for bin
    ssqr_bins.append(x)  # Save for later
    target_power = default_target_spectrum(x) # Could have a different target
    target_spectrum.append(target_power)
    if mapCC < 0.143 and nkept > 3: # Only fit initial signal estimate where clear overall
      continue
    signal = meanfsq * mapCC / target_power
    if signal <= 0.:
      continue
    y = math.log(signal)
    w = mc1sel.size()
    nkept += 1
    sumw += w
    sumwx += w * x
    sumwy += w * y
    sumwx2 += w * x**2
    sumwxy += w * x * y

  denominator = sumw * sumwx2 - sumwx**2
  if denominator == 0.:
    raise Sorry("Problem with fitting signal to mapCC curve")
  slope = (sumw * sumwxy - (sumwx * sumwy)) / denominator
  intercept = (sumwy - slope * sumwx) / sumw
  wilson_scale_signal = math.exp(intercept)
  wilson_b_signal = -4 * slope
  if verbosity > 0:
    print("\n  Wilson B for signal power: ", wilson_b_signal, file=log)
  n_bins = ssqr_bins.size()

  sigmaT_bins = [1.]*n_bins
  start_params = []
  start_params.append(wilson_scale_signal/3.5) # Asqr_scale, factor out low-res BEST value
  start_params.extend(sigmaT_bins)
  wilson_u=adptbx.b_as_u(wilson_b_signal)
  asqr_beta=list(adptbx.u_iso_as_beta(mc1.unit_cell(), wilson_u))
  start_params.extend(asqr_beta)

  # create inputs for the minimizer's run method
  # Constrained refinement using prior parameters could be revisited later.
  # However, this would require all cryo-EM maps to obey the assumption
  # that half-maps are completely unmasked.
  macro1 = ["default"]
  macro2 = ["default"]
  protocol = [macro1, macro2] # overall minimization protocol
  ncyc = 100                  # maximum number of microcycles per macrocycle
  minimizer_type = "bfgs"     # minimizer, bfgs or newton
  if verbosity < 4:
    study_params = False      # flag for calling studyparams procedure
  else:
    study_params = True
  level=verbosity      # 0/1/2/3/4 for mute/log/verbose/debug/testing

  # create instances of refine and minimizer
  refine_cryoem_signal = RefineCryoemSignal(
    sumfsqr_lm = sumfsqr_local_mean, f1f2cos_lm = f1f2cos_local_mean,
    deltafsqr_lm = deltafsqr_local_mean, r_star = r_star, ssqr_bins = ssqr_bins,
    target_spectrum = target_spectrum, start_x = start_params)
  minimizer = Minimizer(output_level=level)

  # Run minimizer
  minimizer.run(refine_cryoem_signal, protocol, ncyc, minimizer_type, study_params)
  refined_params=refine_cryoem_signal.x

  # Extract and report refined parameters
  i_par = 0
  asqr_scale = refined_params[i_par] # Not used for correction: leave map on original scale
  i_par += 1
  sigmaT_bins = refined_params[i_par:i_par + n_bins]
  i_par += n_bins
  asqr_beta = tuple(refined_params[i_par:i_par + 6])
  i_par += 6
  assert (i_par == len(refined_params))

  # Convert asqr_beta to a_beta for reporting
  a_beta = tuple(flex.double(asqr_beta)/2)

  # Convert beta parameters to Baniso for (optional) use and output
  a_baniso = adptbx.u_as_b(adptbx.beta_as_u_cart(mc1.unit_cell(), a_beta))

  # For debugging and data inspection, set make_intermediate_files true, but false normally
  make_intermediate_files = False

  if verbosity > 0:
    print("\nRefinement of scales and error terms completed\n", file=log)
    print("\nParameters for A and BEST curve correction", file=log)
    print("  A overall scale: ",math.sqrt(asqr_scale), file=log)
    for i_bin in range(n_bins):
      print("  Bin #", i_bin + 1, "BEST curve correction: ", sigmaT_bins[i_bin], file=log)
    print("  A tensor as beta:", a_beta, file=log)
    print("  A tensor as Baniso: ", a_baniso, file=log)
    es = adptbx.eigensystem(a_baniso)
    print("  Eigenvalues and eigenvectors:", file=log)
    for iv in range(3):
      print("  ",es.values()[iv],es.vectors(iv), file=log)

    sys.stdout.flush()

  # Loop over bins to compute expectedE and Dobs for each Fourier term
  # Start with mean of half-map Fourier terms and make Miller array for Dobs
  expectE = mc1.customized_copy(data = (mc1.data() + mc2.data())/2)
  expectE.use_binner_of(mc1)
  dobs = expectE.customized_copy(data=flex.double(expectE.size(),0))
  if make_intermediate_files:
    sigmaS = expectE.customized_copy(data=flex.double(expectE.size(),0))
    sigmaE = expectE.customized_copy(data=flex.double(expectE.size(),0))
  i_bin_used = 0 # Keep track in case full range of bins not used
  if verbosity > 0:
    print("MapCC before and after rescaling as a function of resolution", file=log)
    print("Bin   <ssqr>   Nterms   mapCC_before   mapCC_after", file=log)
  for i_bin in mc1.binner().range_used():
    sel = expectE.binner().selection(i_bin)
    eEsel = expectE.select(sel)

    # Make Miller array as basis for computing aniso corrections in this bin
    ones_array = flex.double(eEsel.size(), 1)
    all_ones = eEsel.customized_copy(data=ones_array)
    beta_miller_Asqr = all_ones.apply_debye_waller_factors(
      u_star=adptbx.beta_as_u_star(asqr_beta))
    u_terms = (asqr_scale * sigmaT_bins[i_bin_used]
      * target_spectrum[i_bin_used]) * beta_miller_Asqr.data()

    # Noise is half of deltafsqr power
    sigmaE_terms = (deltafsqr_local_mean.data().select(sel)) / 2.
    f1f2cos_terms = f1f2cos_local_mean.data().select(sel)
    sumfsqr_terms = sumfsqr_local_mean.data().select(sel)

    # Compute the relative weight between the local statistics and the overall
    # anisotropic signal model, as in refinement code.
    steep = 9.
    local_weight = 0.95
    fsc = f1f2cos_terms/(sumfsqr_terms/2.)
    # Make sure fsc is in range 0 to 1
    fsc = (fsc + flex.abs(fsc)) / 2
    fsc = (fsc + 1 - flex.abs(fsc - 1)) / 2
    wt_terms = flex.exp(steep*fsc)
    wt_terms = 1. - local_weight * (wt_terms/(math.exp(steep/2) + wt_terms))
    sigmaS_terms = wt_terms*u_terms + (1.-wt_terms)*f1f2cos_terms
    # Make sure sigmaS is non-negative
    sigmaS_terms = (sigmaS_terms + flex.abs(sigmaS_terms))/2
    scale_terms = 1./flex.sqrt(sigmaS_terms + sigmaE_terms/2.)
    # Arrange Dobs calculation so sigmaS can be zero
    dobs_terms = flex.sqrt(sigmaS_terms / (sigmaS_terms + sigmaE_terms/2.))
    expectE.data().set_selected(sel, expectE.data().select(sel) * scale_terms)
    dobs.data().set_selected(sel, dobs_terms)

    if make_intermediate_files:
      sigmaS.data().set_selected(sel, sigmaS_terms)
      sigmaE.data().set_selected(sel, sigmaE_terms)

    # Apply corrections to mc1 and mc2 to compute mapCC after rescaling
    # SigmaE variance is twice as large for half-maps before averaging
    scale_terms_12 = 1./flex.sqrt(sigmaS_terms + sigmaE_terms)
    # Overwrite mc1 and mc2, but not needed later.
    mc1.data().set_selected(sel, mc1.data().select(sel) * scale_terms_12)
    mc2.data().set_selected(sel, mc2.data().select(sel) * scale_terms_12)
    mc1sel = mc1.select(sel)
    mc2sel = mc2.select(sel)
    mapCC = mc1sel.map_correlation(other=mc2sel)
    if verbosity > 0:
      print(i_bin_used+1, ssqr_bins[i_bin_used], mc1sel.size(), mapCC_bins[i_bin_used], mapCC, file=log)
    mapCC_bins[i_bin_used] = mapCC # Update for returned output
    i_bin_used += 1

  # Compensate for numerical instability when sigmaS is extremely small, in which
  # case dobs is very small and expectE can be very large
  sel = dobs.data() < 0.00001
  expectE.data().set_selected(sel,1.)

  if make_intermediate_files: # For debugging and investigating signal/noise
    # Write out sigmaS, sigmaE and Dobs both as mtz files and intensities-as-maps

    write_mtz(sigmaS,"sigmaS.mtz","sigmaS")
    results = intensities_as_expanded_map(work_mm,sigmaS)
    coeffs_as_map = results.mm_data # map_manager with intensities
    coeffs_as_map.write_map("sigmaS.map")

    write_mtz(sigmaE,"sigmaE.mtz","sigmaE")
    results = intensities_as_expanded_map(work_mm,sigmaE)
    coeffs_as_map = results.mm_data # map_manager with intensities
    coeffs_as_map.write_map("sigmaE.map")

    write_mtz(dobs,"Dobs.mtz","Dobs")
    results = intensities_as_expanded_map(work_mm,dobs)
    coeffs_as_map = results.mm_data # map_manager with intensities
    coeffs_as_map.write_map("Dobs.map")

  return group_args(expectE=expectE, dobs=dobs)

def map_errors_from_d_min(working_mmm, d_min, verbosity=1, log=sys.stdout):
  """
  Use generic relationship between d_min and shape of FSC curve to compute
  radially-symmetry Dobs and Emean values from a single full map.
  Make weighted map coeffs for region.

  Compulsory arguments:
  working_mmm: map_model_manager object containing single map from reconstruction
  d_min: target resolution, either best resolution for map or resolution for
    target region. Legal value must be provided

  Optional arguments:
  verbosity: 0/1/2/3/4 for mute/log/verbose/debug/testing
  log: destination for standard output
  """

  wuc = working_mmm.crystal_symmetry().unit_cell()
  ucpars = wuc.parameters()
  d_max = max(ucpars[0], ucpars[1], ucpars[2])
  map_coeffs = working_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max)
  e_array = map_coeffs.amplitudes()
  e_array.setup_binner(auto_binning = True)
  e_array = e_array.quasi_normalize_structure_factors()
  expectE = e_array.phase_transfer(phase_source = map_coeffs.phases(deg=False), deg = False)
  # The signal_ratio and deltaB are computed from a database investigation of nearly
  # 500 EMDB entries, to match the FSC curves associated with different limits
  signal_ratio_0 = 34.468*math.exp(3.00533/d_min + 4.27895/d_min**2)
  deltaB = 21.3225*d_min**2 + 12.0213*d_min+17.1158
  signal_ratio = expectE.customized_copy(data=flex.double(expectE.size(), signal_ratio_0))
  signal_ratio = signal_ratio.apply_debye_waller_factors(b_iso=deltaB)
  dobs = expectE.customized_copy(data=signal_ratio.data()/(1.+signal_ratio.data())) # FSC
  dobs = dobs.customized_copy(data=flex.sqrt(2*dobs.data()/(1.+dobs.data()))) # FSC to Dobs
  return group_args(expectE=expectE, dobs=dobs)

def assess_cryoem_errors(
        mmm=None,
        d_min=None,
        half_maps_provided=False,
        determine_ordered_volume=True,
        sphere_cent=None,
        radius=None,
        double_map_box=False,
        sphere_points=500,
        shift_map_origin=True,
        ordered_mask_id=None,
        fixed_mask_id=None,
        verbosity=1,
        keep_full_map=False,
        log=sys.stdout):
  """
  Determine error parameters from half-maps (preferred) or full map, and make
  weighted map coeffs for region.

  Compulsory arguments:
  mmm: map_model_manager object containing map or half-maps from reconstruction
    and optionally fixed model plus fixed model mask
  d_min: target resolution, either best resolution for map or resolution for
    target region

  Optional arguments:
  determine_ordered_volume: flag for whether ordered volume should be assessed
    (only applicable for half-maps)
  ordered_mask_id: identifier for mask defining ordered volume
  sphere_cent: center of sphere defining target region for analysis
    default is center of map
  radius: radius of sphere
    default (when sphere center not defined either) is 1/3 narrowest map width
  sphere_points: number of reflections to use for local spherical average
  shift_map_origin: should map coefficients be shifted to correspond to
    original origin, rather than the origin being the corner of the box,
    default True
  define_ordered_volume: should ordered volume over full map be evaluated, to
    compare with cutout volume,
    default True
  keep_full_map: don't mask or box the input map
    default False, use with caution
  verbosity: 0/1/2/3/4 for mute/log/verbose/debug/testing
  """

  if verbosity > 0:
    print("\nPrepare map for docking by assessing signal and errors", file=log)

  # Start from two half-maps and ordered volume mask in map_model_manager
  # Determine ordered volume in whole reconstruction and fraction that will be in sphere
  unit_cell = mmm.map_manager().unit_cell()
  unit_cell_grid = mmm.map_manager().unit_cell_grid
  map_grid = mmm.map_manager().map_data().all()
  if (not mmm.map_manager().is_full_size()):
    print("\nERROR: Unit cell grid does not match map dimensions:")
    print("Unit cell grid: ", unit_cell_grid)
    print("Map dimensions: ", map_grid)
    raise Sorry("Provided maps must be full-size")
  spacings = get_grid_spacings(unit_cell,unit_cell_grid)
  ucpars = unit_cell.parameters()
  # Determine voxel count for fixed model mask
  mm_model_mask = mmm.get_map_manager_by_id(fixed_mask_id)
  voxels_in_fixed_model = 0
  if mm_model_mask is not None:
    model_mask_data = mm_model_mask.map_data()
    mask_sel = model_mask_data > 0
    voxels_in_fixed_model = mask_sel.count(True)

  # Keep track of shifted origin
  origin_shift = mmm.map_manager().origin_shift_grid_units
  if flex.sum(flex.abs(flex.double(origin_shift))) > 0:
    shifted_origin = True
  else:
    shifted_origin = False
  if sphere_cent is None:
    # Default to sphere in center of cell extending 2/3 distance to nearest edge
    # sphere_cent_grid and sphere_cent_map are relative to map origin
    sphere_cent_grid = [round(n/2) for n in unit_cell_grid]
    sphere_cent = [(g * s) for g,s in zip(sphere_cent_grid, spacings)]
    if radius is None:
      radius = min(ucpars[0], ucpars[1], ucpars[2])/3.
  else:
    if radius is None:
      raise Sorry("Radius must be provided if sphere_cent is defined")
    # Force sphere center to be on a map grid point for easier comparison of different spheres
    # Account for origin shift if present in input maps
    sphere_cent_frac = unit_cell.fractionalize(sphere_cent)
    sphere_cent_grid = [round(n * f) for n,f in zip(unit_cell_grid, sphere_cent_frac)]
    if shifted_origin:
      sphere_cent_grid = [(c - o) for c,o in zip(sphere_cent_grid, origin_shift)]
  for i in range(3):
    if sphere_cent_grid[i] > map_grid[i] or sphere_cent_grid[i] < 0:
      raise Sorry("Sphere centre is outside map grid")
  sphere_cent_map = [(g * s) for g,s in zip(sphere_cent_grid, spacings)]
  sphere_cent_map = flex.double(sphere_cent_map)

  if determine_ordered_volume:
    ordered_mm = mmm.get_map_manager_by_id(map_id=ordered_mask_id)
    total_ordered_volume = flex.mean(ordered_mm.map_data()) * unit_cell.volume()
    ordered_volume_in_sphere = get_ordered_volume_exact(ordered_mm, sphere_cent_map, radius)
    fraction_scattering = ordered_volume_in_sphere / total_ordered_volume
  else:
    fraction_scattering = None

  # Get map coefficients for maps after spherical masking
  # Define box big enough to hold sphere plus soft masking
  # Warn if sphere too near edge of full map for soft mask to reach zero
  # or if less than about 85% of sphere would be within full map
  # (outside by >= half-radius)
  boundary_to_smoothing_ratio = 2
  # d_min may be waiting to be evaluated, if there are half-maps
  sensible_d_min = 2.5*max(spacings) # Lower end of typical range from EMDB
  minimum_possible_d_min = 2.*max(spacings)
  if d_min is None or d_min < minimum_possible_d_min:
    if not half_maps_provided:
      print("Given map spacing, minimum possible value for d_min is ",
            minimum_possible_d_min)
      raise Sorry("Sensible value for d_min must be given for single map")
    soft_mask_radius = sensible_d_min
  else:
    soft_mask_radius = d_min
  padding = soft_mask_radius * boundary_to_smoothing_ratio
  cushion = flex.double(3,radius+padding)
  cart_min = flex.double(sphere_cent_map) - cushion
  cart_max = flex.double(sphere_cent_map) + cushion
  max_outside = 0.
  for i in range(3):
    if cart_min[i] < 0:
      max_outside = max(max_outside, -cart_min[i])
    if cart_max[i] > ucpars[i]-spacings[i]:
      max_outside = max(max_outside, cart_max[i]-(ucpars[i]-spacings[i]))
  if max_outside > padding + radius/2:
    print("\nWARNING: substantial fraction of sphere is outside map volume", file=log)
  elif max_outside > padding:
    print("\nWARNING: sphere is partially outside map volume", file=log)
  elif max_outside > 0:
    print("\nWARNING: sphere too near map edge to allow full extent of smooth masking", file=log)

  cs = mmm.crystal_symmetry()
  uc = cs.unit_cell()

  # Set some parameters that will be overwritten if masked and cut out
  over_sampling_factor = 1.
  box_volume = uc.volume()

  if keep_full_map: # Rearrange later so this choice comes first?
    working_mmm = mmm.deep_copy()
  else:
    # Box the map within xyz bounds, converted to map grid units
    lower_frac = uc.fractionalize(tuple(cart_min))
    upper_frac = uc.fractionalize(tuple(cart_max))
    all_orig = mmm.map_data().all()
    lower_bounds = [int(math.floor(f * n)) for f, n in zip(lower_frac, all_orig)]
    upper_bounds = [int(math.ceil( f * n)) for f, n in zip(upper_frac, all_orig)]
    working_mmm = mmm.extract_all_maps_with_bounds(
        lower_bounds=lower_bounds, upper_bounds=upper_bounds)
    # Make and apply spherical mask
    working_mmm.create_spherical_mask(
      soft_mask_radius=soft_mask_radius,
      boundary_to_smoothing_ratio=boundary_to_smoothing_ratio)
    working_mmm.apply_mask_to_maps()
    if double_map_box:
      delta_cushion = 2*radius - cushion[0]
      delta_cushion_grid = int(math.ceil(delta_cushion/spacings[0]))
      # Note that extract_all_maps_with_bounds works in terms of the lower
      # and upper ranges of the map being operated on. Here we are padding the
      # map out with zeroes on all sides.
      new_lower_bounds = tuple(-flex.int(3,delta_cushion_grid))
      new_upper_bounds = tuple(flex.int(upper_bounds) - flex.int(lower_bounds) + flex.int(3,delta_cushion_grid))
      working_mmm = working_mmm.extract_all_maps_with_bounds(
          lower_bounds=new_lower_bounds, upper_bounds=new_upper_bounds,
          stay_inside_current_map = False)
    mask_info = working_mmm.mask_info()

    working_map_size = working_mmm.map_data().size()
    box_volume = uc.volume() * (
        working_map_size/mmm.map_data().size())

    voxels_in_masked_model = 0
    if voxels_in_fixed_model > 0:
      mask_info = working_mmm.mask_info(mask_id=fixed_mask_id)
      mm_model_mask = working_mmm.get_map_manager_by_id(fixed_mask_id)
      model_mask_data = mm_model_mask.map_data()
      voxels_in_masked_model = working_map_size * mask_info.mean

    # Calculate an oversampling ratio defined as the ratio between the size of
    # the cut-out cube and the size of a cube that could contain a sphere
    # big enough to hold the volume of ordered density. Because the ratio of
    # the size of a cube to the sphere inscribed in it is about 1.9 (which is
    # close to the factor of two between typical protein volume and unit cell
    # volume in a crystal), this should yield likelihood scores on a similar
    # scale to crystallographic ones.
    if determine_ordered_volume:
      assert(ordered_volume_in_sphere > 0.)
      over_sampling_factor = box_volume / (ordered_volume_in_sphere * 6./math.pi)
    else:
      over_sampling_factor = 4. # Guess!

  if verbosity > 0:
    if determine_ordered_volume:
      print("Fraction of full map scattering: ",fraction_scattering, file=log)
    print("Over-sampling factor: ",over_sampling_factor, file=log)
    sys.stdout.flush()

  if half_maps_provided:
    results = assess_halfmap_errors(working_mmm, d_min,
                                    sphere_points=sphere_points,
                                    verbosity=verbosity,
                                    log=log)
  else:
    results = map_errors_from_d_min(working_mmm, d_min,
                                    verbosity=verbosity,
                                    log=log)
  expectE = results.expectE
  dobs = results.dobs
  # d_min may have been evaluated from half-maps
  eps = 1.E-6 # Add a little slop for rare rounding errors leading to mismatch of limits
  d_min = flex.min(expectE.d_spacings().data()) - eps
  d_max = flex.max(expectE.d_spacings().data()) + eps

  # If there is a fixed model and it occupies a significant part of the sphere,
  # compute the structure factors corresponding to the cutout density for the fixed
  # model.

  masked_model_contributes = False
  masked_model_E = None
  masked_model_sigmaA = None
  fraction_masked = voxels_in_masked_model/working_map_size
  if fraction_masked > 0.01:
    if verbosity > 1:
      print("Masked voxels from fixed model, total in box, fraction: ",
            voxels_in_masked_model, working_map_size, fraction_masked, file=log)
    masked_model_contributes = True
    mm_masked_model = working_mmm.get_map_manager_by_id(map_id='fixed_atom_map')
    # mm_masked_model.write_map('masked_model.map')
    masked_model_F = mm_masked_model.map_as_fourier_coefficients(d_min=d_min,d_max=d_max)

  if masked_model_contributes:
    results = sigmaA_from_model_in_map(expectE, dobs, masked_model_F, over_sampling_factor, verbosity)
    masked_model_sigmaA = results.model_sigmaA
    masked_model_E = results.Ecalc
    # Possibly temporary choice: modify expectE instead of adding fixed contribution in phasertng
    expectE_mod = expectE.customized_copy(data=expectE.data()
                  - masked_model_sigmaA.data()*masked_model_E.data())
  else:
    expectE_mod = expectE.deep_copy()

  # Return map_model_managers with weighted maps
  wEmean = dobs*expectE_mod
  working_mmm.add_map_from_fourier_coefficients(map_coeffs=wEmean, map_id='map_manager_wtd')
  sigmaA = 0.9 # Default for likelihood-weighted map
  dosa = sigmaA*dobs.data()
  lweight = dobs.customized_copy(data=(2*dosa)/(1.-flex.pow2(dosa)))
  wEmean2 = lweight*expectE_mod
  working_mmm.add_map_from_fourier_coefficients(map_coeffs=wEmean2, map_id='map_manager_lwtd')
  if half_maps_provided:
    working_mmm.remove_map_manager_by_id('map_manager_1')
    working_mmm.remove_map_manager_by_id('map_manager_2')

  shift_cart = working_mmm.shift_cart()
  if shift_map_origin:
    ucwork = expectE.crystal_symmetry().unit_cell()
    # shift_cart is position of original origin in boxed-map coordinate system
    # shift_frac should correspond to what has to be done to a model to put it
    # into the map, i.e. move it in the opposite direction
    shift_frac = ucwork.fractionalize(shift_cart)
    shift_frac = tuple(-flex.double(shift_frac))
    expectE = expectE.translational_shift(shift_frac)
    expectE_mod = expectE_mod.translational_shift(shift_frac)
    if masked_model_contributes:
      masked_model_E = masked_model_E.translational_shift(shift_frac)
      # write_mtz(expectE,"expectE.mtz","eE")
      # write_mtz(expectE_mod,"expectE_mod.mtz","expectEmod")
      # write_mtz(masked_model_E,"masked_model.mtz","masked_model")

  return group_args(
    new_mmm = working_mmm,
    shift_cart = shift_cart,
    sphere_center = list(sphere_cent_map),
    expectE = expectE_mod, dobs = dobs, # Temporary kludge
    # expectE = expectE, dobs = dobs,
    masked_model_E = masked_model_E, masked_model_sigmaA = masked_model_sigmaA,
    over_sampling_factor = over_sampling_factor,
    fraction_scattering = fraction_scattering)

# Command-line interface using argparse
def run():
  """
  Prepare cryo-EM map for docking by preparing weighted MTZ file.

  Command-line arguments (keyworded):
  Either two half-maps or a single full map must be provided
    --map:   name of file containing the full final reconstructed map
    --map1: name of file containing the first half-map from a reconstruction
    --map2: name of file containing the second half-map

    --d_min: desired resolution, either best for whole map or for local region

  Content of full reconstruction should be provided and is compulsory if the sphere
  parameters are not defined, either explicitly or in the form of a model
    --protein_mw: molecular weight expected for protein component of ordered density
          in full map
    --nucleic_mw: same for nucleic acid component

  Optional parameters specifying centre and radius of cut-out sphere
    --model: PDB file for placed model used to define the centre and radius of
          the cut-out sphere, used as an alternative to sphere_cent and radius
    --sphere_cent: Centre of sphere defining target map region (3 floats)
          defaults to centre of map
    --radius: radius of sphere (1 float)
          must be given if sphere_cent defined, otherwise
          defaults to either radius enclosing 98% of ordered volume (if determined)
          or narrowest extent of input map divided by 4

  Option to create output map box twice as wide as sphere, oversampling the
  transform of the sphere so that only alternate hkl terms are needed in likelihood
    --double_map_box: make box with double diameter of sphere
          default False

  Optional fixed model, contribution of which is removed from output data
    --fixed_model: PDB file for fixed background model

  Optional parameters with sensible defaults
    --sphere_points: number of reflections to be used for averaging in Fourier space
          default 500
    --no_shift_map_origin: don't shift output mtz file to match input map on its origin
          default False
    --no_determine_ordered_volume: don't determine ordered volume
          default False
    --file_root: root name for output files
    --write_maps: write sharpened and likelihood-weighted maps
  Control level of output. Default is logfile level of output
    --mute (or -m): mute output
    --verbose (or -v): verbose output
    --testing: extra verbose output for debugging
  """
  import argparse
  dm = DataManager()
  dm.set_overwrite(True)

  parser = argparse.ArgumentParser(
          description='Prepare cryo-EM map for docking')
  parser.add_argument('--map',
                      help='Map file for final reconstruction, instead of half-maps')
  parser.add_argument('--map1', help='Map file for half-map 1')
  parser.add_argument('--map2', help='Map file for half-map 2')
  parser.add_argument('--d_min', help='d_min for maps', type=float)
  parser.add_argument('--protein_mw',
                      help='Molecular weight of protein component of full map',
                      type=float)
  parser.add_argument('--nucleic_mw',
                      help='Molecular weight of nucleic acid component of full map',
                      type=float)
  parser.add_argument('--model',
                      help='''Optional roughly-placed model defining size and
                              position of cut-out sphere for docking.
                              Orientation of model is ignored.''')
  parser.add_argument('--sphere_cent',
                      help='Centre of sphere for docking, if no model provided',
                      nargs=3, type=float)
  parser.add_argument('--radius',
                      help='Radius of sphere for docking, if no model provided',
                      type=float)
  parser.add_argument('--double_map_box',
                      help='Make map box twice as wide as sphere',
                      action = 'store_true')
  parser.add_argument('--fixed_model',
                      help='Optional fixed model accounting for explained map features')
  parser.add_argument('--sphere_points',help='Target nrefs in averaging sphere',
                      type=float, default=500.)
  parser.add_argument('--no_shift_map_origin', action='store_true')
  parser.add_argument('--no_determine_ordered_volume', action='store_true')
  parser.add_argument('--file_root',
                      help='Root of filenames for output')
  parser.add_argument('--write_maps', help = 'Write sharpened and likelihood-weighted maps',
                      action = 'store_true')
  parser.add_argument('-m', '--mute', help = 'Mute output', action = 'store_true')
  parser.add_argument('-v', '--verbose', help = 'Set output as verbose',
                      action = 'store_true')
  parser.add_argument('--testing', help='Set output as testing level', action='store_true')
  args = parser.parse_args()
  d_min = args.d_min
  verbosity = 1
  if args.mute: verbosity = 0
  if args.verbose: verbosity = 2
  if args.testing: verbosity = 4
  shift_map_origin = not(args.no_shift_map_origin)
  determine_ordered_volume = not(args.no_determine_ordered_volume)

  sphere_points = args.sphere_points

  cutout_specified = False
  sphere_cent = None
  radius = None
  double_map_box = args.double_map_box

  model = None
  if args.model is not None:
    model_file = args.model
    model = dm.get_model(model_file)
    if args.sphere_cent is not None:
      raise Sorry('Either model or sphere specification may be given but not both')
    sphere_cent, radius = sphere_enclosing_model(model)
    radius = radius + d_min # Expand to allow width for density
    cutout_specified = True

  fixed_model = None
  if args.fixed_model is not None:
    model_file = args.fixed_model
    fixed_model = dm.get_model(model_file)

  if args.sphere_cent is not None:
    if args.radius is None:
      raise Sorry("Radius must be provided if sphere_cent is defined")
    sphere_cent = tuple(args.sphere_cent)
    radius = args.radius
    cutout_specified = True
  if (not determine_ordered_volume and not cutout_specified):
    raise Sorry("Must determine ordered volume if cutout not specified")

  protein_mw = None
  nucleic_mw = None
  if (args.protein_mw is None) and (args.nucleic_mw is None):
    if cutout_specified:
      determine_ordered_volume = False
    if determine_ordered_volume and not cutout_specified:
      raise Sorry("At least one of protein_mw or nucleic_mw must be given")
  if args.protein_mw is not None:
    protein_mw = args.protein_mw
  if args.nucleic_mw is not None:
    nucleic_mw = args.nucleic_mw

  # Create map_model_manager containing half-maps, or full map if no half-maps
  half_maps_provided = True
  if args.map1 is not None:
    if args.map2 is None:
      raise Sorry("Half-maps must be provided in a pair")
    map1_filename = args.map1
    mm1 = dm.get_real_map(map1_filename)
    map2_filename = args.map2
    mm2 = dm.get_real_map(map2_filename)
    mmm = map_model_manager(map_manager_1=mm1, map_manager_2=mm2)
  else:
    if args.map is None:
      raise Sorry("Either half-maps or full map must be provided")
    half_maps_provided = False
    map_filename = args.map
    mm = dm.get_real_map(map_filename)
    mmm = map_model_manager(map_manager=mm)

  if (fixed_model):
    mmm.add_model_by_id(model=fixed_model, model_id='fixed_model')
    mmm.generate_map(model=fixed_model, d_min=d_min, map_id='fixed_atom_map')

  ordered_mask_id = None
  if determine_ordered_volume:
    ordered_mask_id = 'ordered_volume_mask'
    add_ordered_volume_mask(mmm, d_min,
        protein_mw=protein_mw, nucleic_mw=nucleic_mw,
        ordered_mask_id=ordered_mask_id)
    if verbosity>1:
      ordered_mm = mmm.get_map_manager_by_id(map_id=ordered_mask_id)
      if args.file_root is not None:
        map_file_name = args.file_root + "_ordered_volume_mask.map"
      else:
        map_file_name = "ordered_volume_mask.map"
      ordered_mm.write_map(map_file_name)

  fixed_mask_id = None
  if fixed_model:
    fixed_mask_id = 'mask_around_atoms'
    mask_fixed_model_region(mmm, d_min,
                            fixed_model=fixed_model,
                            ordered_mask_id=ordered_mask_id,
                            fixed_mask_id=fixed_mask_id)

  # Default to sphere containing most of ordered volume
  if not cutout_specified:
    ucpars = mm1.unit_cell().parameters()
    sphere_cent = (ucpars[0]/2,ucpars[1]/2,ucpars[2]/2)
    if determine_ordered_volume:
      target_completeness = 0.98
      radius = get_mask_radius(mmm.get_map_manager_by_id(ordered_mask_id),target_completeness)
    else:
      radius = min(ucpars[0],ucpars[1],ucpars[2])/4.

  # Refine to get scale and error parameters for docking region
  results = assess_cryoem_errors(
                        mmm=mmm,
                        d_min=d_min,
                        half_maps_provided=half_maps_provided,
                        determine_ordered_volume=determine_ordered_volume,
                        sphere_cent=sphere_cent,
                        radius=radius,
                        double_map_box=double_map_box,
                        sphere_points=sphere_points,
                        shift_map_origin=shift_map_origin,
                        ordered_mask_id=ordered_mask_id,
                        fixed_mask_id=fixed_mask_id,
                        verbosity=verbosity)

  expectE = results.expectE
  mtz_dataset = expectE.as_mtz_dataset(column_root_label='Emean')
  dobs = results.dobs
  mtz_dataset.add_miller_array(
      dobs,column_root_label='Dobs',column_types='W')
  mtz_object=mtz_dataset.mtz_object()

  if args.file_root is not None:
    mtzout_file_name = args.file_root + ".mtz"
    mapout1_file_name = args.file_root + "_weighted.map"
    mapout2_file_name = args.file_root + "_likelihood_weighted.map"
  else:
    mtzout_file_name = "weighted_map_data.mtz"
    mapout1_file_name = "weighted.map"
    mapout2_file_name = "likelihood_weighted.map"
  print ("Writing mtz for docking as", mtzout_file_name)
  if not shift_map_origin:
    shift_cart = results.shift_cart
    print ("Origin of full map relative to mtz:", shift_cart)
  dm.write_miller_array_file(mtz_object, filename=mtzout_file_name)
  if args.write_maps:
    print ("Writing weighted map as", mapout1_file_name)
    print ("Writing likelihood-weighted map as", mapout2_file_name)
    new_mmm = results.new_mmm
    new_mmm.write_map(map_id='map_manager_wtd', file_name = mapout1_file_name)
    new_mmm.write_map(map_id='map_manager_lwtd', file_name = mapout2_file_name)
  over_sampling_factor = results.over_sampling_factor
  fraction_scattering = results.fraction_scattering
  print ("Over-sampling factor for Fourier terms:",over_sampling_factor)
  print ("Fraction of total scattering:",fraction_scattering)
  sys.stdout.flush()

if __name__ == "__main__":
  run()
