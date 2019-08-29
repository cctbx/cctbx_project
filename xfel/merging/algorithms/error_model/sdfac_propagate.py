from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
import math
from rstbx.symmetry.constraints.parameter_reduction \
    import symmetrize_reduce_enlarge
from scitbx.matrix import sqr, col

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base
from xfel.merging.algorithms.error_model.sdfac_refine_lbfgs import finite_difference

from libtbx import group_args
from six.moves import zip

"""
Classes to support propagating erros after postrefinement in cxi.merge
"""

def r2d(radians):
  return 180*radians/math.pi

# Bucket to hold refinable error terms
class error_terms(group_args):
  @staticmethod
  def from_x(x):
    return error_terms(sigma_thetax  = x[0],
                       sigma_thetay  = x[1],
                       sigma_lambda  = x[2],
                       sigma_deff    = x[3],
                       sigma_gstar   = x[4:])

  def to_x(self):
    x = flex.double([self.sigma_thetax,
                     self.sigma_thetay,
                     self.sigma_lambda,
                     self.sigma_deff])
    x.extend(self.sigma_gstar)
    return x

class sdfac_propagate(error_modeler_base):
  def __init__(self, scaler, error_terms = None, verbose = True):
    super(sdfac_propagate, self).__init__(scaler)

    # Note, since the rs algorithm doesn't explicitly refine eta and deff separately, but insteads refines RS,
    # assume rs only incorporates information from deff and set eta to zero.
    ct = scaler.crystal_table
    ct['deff'] = 1/ct['RS']
    ct['eta'] = flex.double(len(ct), 0)

    self.error_terms = error_terms
    self.verbose = verbose

  def finite_difference(self, parameter_name, table, DELTA = 1.E-7):
    """ Compute finite difference given a parameter name """
    refls = self.scaler.ISIGI

    def target():
      r = self.compute_intensity_parameters()
      return refls['iobs'] / r['D']

    functional = target()

    if parameter_name.startswith('c'):
      # Handle the crystal parameters
      from scitbx.matrix import sqr
      cryst_param = int(parameter_name.lstrip('c'))
      parameter_name = 'b_matrix'
      current = table[parameter_name]*1 # make a copy

      sre = symmetrize_reduce_enlarge(self.scaler.params.target_space_group.group())
      for i in range(len(table)):
        sre.set_orientation(orientation=table['b_matrix'][i])
        vals = list(sre.forward_independent_parameters())
        vals[cryst_param] += DELTA
        newB = sqr(sre.backward_orientation(vals).reciprocal_matrix())
        table['b_matrix'][i] = newB
    else:
      current = table[parameter_name]
      table[parameter_name] = current + DELTA

    dfunctional = target()
    #calculate by finite_difference
    finite_g = (dfunctional-functional )/DELTA

    table[parameter_name] = current

    return finite_g

  def compute_intensity_parameters(self):
    """ Create a new reflection table with all the derived parameters needed
    to apply corrections from RS postrefinement """
    refls = self.scaler.ISIGI
    ct = self.scaler.crystal_table

    rx = flex.mat3_double() # crystal rotation around x
    ry = flex.mat3_double() # crystal rotation around y
    u = flex.mat3_double()  # U matrix (orientation)
    b = flex.mat3_double()  # B matrix (cell parameters)
    wavelength = flex.double()
    G = flex.double()       # scaling gfactor
    B = flex.double()       # wilson B factor
    s0 = flex.vec3_double() # beam vector
    deff = flex.double()    # effective domain size
    eta = flex.double()     # effective mosaic domain misorientation angle

    ex = col((1,0,0))       # crystal rotation x axis
    ey = col((0,1,0))       # crystal rotation y axis

    for i in range(len(ct)):
      # Need to copy crystal specific terms for each reflection. Equivalent to a JOIN in SQL.
      n_refl = ct['n_refl'][i]
      rx.extend(flex.mat3_double(n_refl, ex.axis_and_angle_as_r3_rotation_matrix(ct['thetax'][i])))
      ry.extend(flex.mat3_double(n_refl, ey.axis_and_angle_as_r3_rotation_matrix(ct['thetay'][i])))
      u.extend(flex.mat3_double(n_refl, ct['u_matrix'][i]))
      b.extend(flex.mat3_double(n_refl, ct['b_matrix'][i]))
      wavelength.extend(flex.double(n_refl, ct['wavelength'][i]))
      G.extend(flex.double(n_refl, ct['G'][i]))
      B.extend(flex.double(n_refl, ct['B'][i]))
      s0.extend(flex.vec3_double(n_refl, (0,0,-1)) * (1/ct['wavelength'][i]))
      deff.extend(flex.double(n_refl, ct['deff'][i]))
      eta.extend(flex.double(n_refl, ct['eta'][i]))

    h          = refls['miller_index_original'].as_vec3_double()
    q          = ry * rx * u * b * h                  # vector pointing from origin of reciprocal space to RLP
    qlen       = q.norms()                            # length of q
    d          = 1/q.norms()                          # resolution
    #rs         = (1/deff)+(eta/(2*d))                # proper formulation of RS
    rs         = 1/deff                               # assumes eta is zero
    rs_sq      = rs*rs                                # square of rs
    s          = (s0+q)                               # vector from center of Ewald sphere to RLP
    slen       = s.norms()                            # length of s
    rh         = slen-(1/wavelength)                  # distance from RLP to Ewald sphere
    p_n        = rs_sq                                # numerator of partiality lorenzian expression
    p_d        = (2. * (rh * rh)) + rs_sq             # denominator of partiality lorenzian expression
    partiality = p_n/p_d
    theta      = flex.asin(wavelength/(2*d))
    epsilon    = -8*B*(flex.sin(theta)/wavelength)**2 # exponential term in partiality
    eepsilon   = flex.exp(epsilon)                    # e^epsilon
    D          = partiality * G * eepsilon            # denominator of partiality lorenzian expression
    thetah     = flex.asin(wavelength/(2*d))          # reflecting angle
    sinthetah  = flex.sin(thetah)
    er         = sinthetah/wavelength                 # ratio term in epsilon

    # save all the columns
    r = flex.reflection_table()
    r['rx']         = rx
    r['ry']         = ry
    r['u']          = u
    r['b']          = b
    r['h']          = h
    r['q']          = q
    r['qlen']       = qlen
    r['D']          = D
    r['rs']         = rs
    r['eta']        = eta
    r['deff']       = deff
    r['d']          = d
    r['s']          = s
    r['slen']       = slen
    r['wavelength'] = wavelength
    r['p_n']        = p_n
    r['p_d']        = p_d
    r['partiality'] = partiality
    r['G']          = G
    r['B']          = B
    r['eepsilon']   = eepsilon
    r['thetah']     = thetah
    r['sinthetah']  = sinthetah
    r['er']         = er

    return r

  def initial_estimates(self):
    # Compute errors by examining distributions of parameters
    ct = self.scaler.crystal_table
    stats_thetax = flex.mean_and_variance(ct['thetax'])
    stats_thetay = flex.mean_and_variance(ct['thetay'])
    stats_lambda = flex.mean_and_variance(ct['wavelength'])
    stats_deff   = flex.mean_and_variance(ct['deff'])
    stats_rs     = flex.mean_and_variance(ct['RS'])
    sigma_thetax = stats_thetax.unweighted_sample_standard_deviation()
    sigma_thetay = stats_thetay.unweighted_sample_standard_deviation()
    sigma_lambda = stats_lambda.unweighted_sample_standard_deviation()
    sigma_deff   = stats_deff.unweighted_standard_error_of_mean()
    sigma_rs     = stats_rs.unweighted_sample_standard_deviation()
    print("ThetaX %.4f +/- %.4f"    %(r2d(stats_thetax.mean()), r2d(sigma_thetax)), file=self.log)
    print("Thetay %.4f +/- %.4f"    %(r2d(stats_thetay.mean()), r2d(sigma_thetay)), file=self.log)
    print("Wavelength %.4f +/- %.4f"%(    stats_lambda.mean(),      sigma_lambda), file=self.log)
    print("DEFF %.4f +/- %.4f"      %(    stats_deff.mean(),        sigma_deff), file=self.log)
    print("RS %.6f +/- %.6f"        %(    stats_rs.mean(),          sigma_rs), file=self.log)

    sre = symmetrize_reduce_enlarge(self.scaler.params.target_space_group.group())
    c_gstar_params = None

    for i in range(len(ct)):
      # Get the G* unit cell parameters from cctbx
      sre.set_orientation(orientation=ct['b_matrix'][i])
      p = sre.forward_independent_parameters()
      if c_gstar_params is None:
        c_gstar_params = [flex.double() for j in range(len(p))]

      for j in range(len(p)):
        c_gstar_params[j].append(p[j])

    # Compute the error in the unit cell terms from the distribution of unit cell parameters provided
    print("Free G* parameters", file=self.log)
    sigma_gstar = flex.double()
    for j in range(len(c_gstar_params)):
      stats  = flex.mean_and_variance(c_gstar_params[j])
      print("G* %d %.4f *1e-5 +/- %.4f *1e-5"%(j, stats.mean()*1e5, stats.unweighted_sample_standard_deviation()*1e5), file=self.log)
      sigma_gstar.append(stats.unweighted_sample_standard_deviation())

    self.error_terms = error_terms(sigma_thetax  = sigma_thetax,
                                   sigma_thetay  = sigma_thetay,
                                   sigma_lambda  = sigma_lambda,
                                   sigma_deff    = sigma_deff,
                                   sigma_gstar   = sigma_gstar)

  def dI_derrorterms(self):
    refls = self.scaler.ISIGI
    ct = self.scaler.crystal_table

    # notation: dP1_dP2 is derivative of parameter 1 with respect to parameter 2. Here,
    # for example, is the derivative of rx wrt thetax
    drx_dthetax = flex.mat3_double()
    dry_dthetay = flex.mat3_double()
    s0hat = flex.vec3_double(len(refls), (0,0,-1))

    ex = col((1,0,0))
    ey = col((0,1,0))

    # Compute derivatives
    sre = symmetrize_reduce_enlarge(self.scaler.params.target_space_group.group())
    gstar_params = None
    gstar_derivatives = None

    for i in range(len(ct)):
      n_refl = ct['n_refl'][i]

      # Derivatives of rx/y wrt thetax/y come from cctbx
      drx_dthetax.extend(flex.mat3_double(n_refl, ex.axis_and_angle_as_r3_derivative_wrt_angle(ct['thetax'][i])))
      dry_dthetay.extend(flex.mat3_double(n_refl, ey.axis_and_angle_as_r3_derivative_wrt_angle(ct['thetay'][i])))

      # Derivatives of the B matrix wrt to the unit cell parameters also come from cctbx
      sre.set_orientation(orientation=ct['b_matrix'][i])
      p = sre.forward_independent_parameters()
      dB_dp = sre.forward_gradients()
      if gstar_params is None:
        assert gstar_derivatives is None
        gstar_params = [flex.double() for j in range(len(p))]
        gstar_derivatives = [flex.mat3_double() for j in range(len(p))]
      assert len(p) == len(dB_dp) == len(gstar_params) == len(gstar_derivatives)
      for j in range(len(p)):
        gstar_params[j].extend(flex.double(n_refl, p[j]))
        gstar_derivatives[j].extend(flex.mat3_double(n_refl, tuple(dB_dp[j])))

    # Compute the scalar terms used while computing derivatives
    self.r = r = self.compute_intensity_parameters()

    # Begin computing derivatives
    sigma_Iobs = refls['scaled_intensity']/refls['isigi']
    dI_dIobs = 1/r['D']

    def compute_dI_dp(dq_dp):
      """ Deriviatives of the scaled intensity I wrt to thetax, thetay and the unit cell parameters
      are computed the same, starting with the deriviatives of those parameters wrt to q """
      dqlen_dp = r['q'].dot(dq_dp)/r['qlen']
      dd_dp    = -(1/(r['qlen']**2)) * dqlen_dp
      drs_dp   = -(r['eta']/(2 * r['d']**2)) * dd_dp
      dslen_dp = r['s'].dot(dq_dp)/r['slen']
      drhsq_dp = 2 * (r['slen'] - (1/r['wavelength'])) * dslen_dp
      dPn_dp   = 2 * r['rs'] * drs_dp
      dPd_dp   = 2 * ((r['rs'] * drs_dp) + drhsq_dp)
      dP_dp    = ((r['p_d'] * dPn_dp)-(r['p_n'] * dPd_dp))/(r['p_d']**2)
      dI_dp    = -(refls['iobs']/(r['partiality']**2 * r['G'] * r['eepsilon'])) * dP_dp
      return dI_dp

    # Derivatives wrt the unit cell parameters
    dI_dgstar = []
    for j in range(len(gstar_params)):
      dI_dgstar.append(compute_dI_dp(r['ry'] * r['rx'] * r['u'] * gstar_derivatives[j] * r['h']))

    # Derivatives wrt the crystal orientation
    dI_dthetax = compute_dI_dp(r['ry'] * drx_dthetax * r['u'] * r['b'] * r['h'])
    dI_dthetay = compute_dI_dp(dry_dthetay * r['rx'] * r['u'] * r['b'] * r['h'])

    # Derivatives wrt to the wavelength
    dthetah_dlambda  = 1/(flex.sqrt(1 - ((r['wavelength']/(2 * r['d']))**2)) * 2 * r['d'])
    den_dlambda      = flex.cos(r['thetah']) * dthetah_dlambda
    der_dlambda      = ((r['wavelength'] * den_dlambda) - r['sinthetah'])/r['wavelength']**2
    depsilon_dlambda = -16 * r['B'] * r['er'] * der_dlambda
    ds0_dlambda      = s0hat*(-1/r['wavelength']**2)
    dslen_dlambda    = r['s'].dot(ds0_dlambda)/r['slen']
    drhsq_dlambda    = 2*(r['slen']-(1/r['wavelength']))*(dslen_dlambda+(1/r['wavelength']**2))
    dP_dlambda       = -2*(r['p_n']/r['p_d']**2) * drhsq_dlambda
    dD_dlambda       = (r['G'] * r['eepsilon'] * dP_dlambda) + (r['partiality'] * r['G'] * r['eepsilon'] * depsilon_dlambda)
    dI_dlambda       = -(refls['iobs']/r['D']**2) * dD_dlambda

    # Derivatives wrt to the deff
    drs_deff = -1/(r['deff']**2)
    dPn_deff = 2 * r['rs'] * drs_deff
    dPd_deff = 2 * r['rs'] * drs_deff
    dP_deff  = ((r['p_d'] * dPn_deff)-(r['p_n'] * dPd_deff))/(r['p_d']**2)
    dI_deff  = -(refls['iobs']/(r['partiality']**2 * r['G'] * r['eepsilon'])) * dP_deff

    # Derivatives wrt to eta (unused for RS refinement)
    # drs_deta = 1/(2*r['d'])
    # dPn_deta = 2 * r['rs'] * drs_deta
    # dPd_deta = 2 * r['rs'] * drs_deta
    # dP_deta  = ((r['p_d']*dPn_deta)-(r['p_n']*dPd_deta))/(r['p_d']**2)
    # dI_deta  = -(refls['iobs']/(r['partiality']**2 * r['G'] * r['eepsilon'])) * dP_deta

    if self.verbose:
      # Show comparisons to finite differences
      n_cryst_params = sre.constraints.n_independent_params()
      print("Showing finite differences and derivatives for each parameter (first few reflections only)")
      for parameter_name, table, derivatives, delta, in zip(['iobs', 'thetax', 'thetay', 'wavelength', 'deff'] + ['c%d'%cp for cp in range(n_cryst_params)],
                                                    [refls, ct, ct, ct, ct] + [ct]*n_cryst_params,
                                                    [dI_dIobs, dI_dthetax, dI_dthetay, dI_dlambda, dI_deff] + dI_dgstar,
                                                    [1e-7]*5 + [1e-11]*n_cryst_params):
        finite_g = self.finite_difference(parameter_name, table, delta)
        print(parameter_name)
        for refl_id in range(min(10, len(refls))):
          print("%d % 21.1f % 21.1f"%(refl_id, finite_g[refl_id], derivatives[refl_id]))
        stats = flex.mean_and_variance(finite_g-derivatives)
        stats_finite = flex.mean_and_variance(finite_g)
        percent = 0 if stats_finite.mean() == 0 else 100*stats.mean()/stats_finite.mean()
        print("Mean difference between finite and analytical: % 24.4f +/- % 24.4f (%8.3f%% of finite d.)"%( \
            stats.mean(), stats.unweighted_sample_standard_deviation(), percent))
        print()

    return [dI_dIobs, dI_dthetax, dI_dthetay, dI_dlambda, dI_deff] + dI_dgstar

  def adjust_errors(self, dI_derrorterms = None, compute_sums = True):
    """ Propagate errors to the scaled and merged intensity errors based on statistical error propagation.
    This uses 1) and estimate of the errors in the post-refined parametes from the observed population
    and 2) partial derivatives of the scaled intensity with respect to each of the post-refined parameters.
    """
    assert self.scaler.params.postrefinement.algorithm == 'rs'

    refls = self.scaler.ISIGI
    ct = self.scaler.crystal_table

    if self.error_terms is None:
      self.initial_estimates()

    if dI_derrorterms is None:
      dI_derrorterms = self.dI_derrorterms()
    dI_dIobs, dI_dthetax, dI_dthetay, dI_dlambda, dI_deff = dI_derrorterms[0:5]
    dI_dgstar = dI_derrorterms[5:]
    sigma_Iobs = refls['scaled_intensity']/refls['isigi']
    r = self.r

    # Propagate errors
    refls['isigi'] = refls['scaled_intensity'] / \
                     flex.sqrt((sigma_Iobs**2 * dI_dIobs**2) +
                     sum([self.error_terms.sigma_gstar[j]**2 * dI_dgstar[j]**2 for j in range(len(self.error_terms.sigma_gstar))]) +
                     (self.error_terms.sigma_thetax**2 * dI_dthetax**2) +
                     (self.error_terms.sigma_thetay**2 * dI_dthetay**2) +
                     (self.error_terms.sigma_lambda**2 * dI_dlambda**2) +
                     (self.error_terms.sigma_deff**2 * dI_deff**2))
    if self.verbose:
      # Show results of propagation
      from scitbx.math import five_number_summary
      all_data = [(refls['iobs'], "Iobs"),
                  (sigma_Iobs, "Original errors"),
                  (1/r['D'], "Total scale factor"),
                  (refls['iobs']/r['D'], "Inflated intensities"),
                  (refls['scaled_intensity']/refls['isigi'], "Propagated errors"),
                  (flex.sqrt(sigma_Iobs**2 * dI_dIobs**2), "Iobs term"),
                  (flex.sqrt(self.error_terms.sigma_thetax**2 * dI_dthetax**2), "Thetax term"),
                  (flex.sqrt(self.error_terms.sigma_thetay**2 * dI_dthetay**2), "Thetay term"),
                  (flex.sqrt(self.error_terms.sigma_lambda**2 * dI_dlambda**2), "Wavelength term"),
                  (flex.sqrt(self.error_terms.sigma_deff**2 * dI_deff**2), "Deff term")] + \
                 [(flex.sqrt(self.error_terms.sigma_gstar[j]**2 * dI_dgstar[j]**2), "Gstar term %d"%j) for j in range(len(self.error_terms.sigma_gstar))]
      print("%20s % 20s % 20s % 20s"%("Data name","Quartile 1", "Median", "Quartile 3"), file=self.log)
      for data, title in all_data:
        fns = five_number_summary(data)
        print("%20s % 20d % 20d % 20d"%(title, fns[1], fns[2], fns[3]), file=self.log)

    if compute_sums:
      # Final terms for cxi.merge
      self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
      self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

      Intensity = refls['scaled_intensity']
      sigma = Intensity / refls['isigi']
      variance = sigma * sigma

      for i in range(len(refls)):
        j = refls['miller_id'][i]
        self.scaler.summed_wt_I[j] += Intensity[i] / variance[i]
        self.scaler.summed_weight[j] += 1 / variance[i]
