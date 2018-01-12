from __future__ import division
from dials.array_family import flex
import math
from rstbx.symmetry.constraints.parameter_reduction \
    import symmetrize_reduce_enlarge
from scitbx.matrix import sqr, col

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base
from xfel.merging.algorithms.error_model.sdfac_refine_lbfgs import finite_difference

def r2d(radians):
  return 180*radians/math.pi

class sdfac_propagate(error_modeler_base):
  def finite_difference(self, parameter_name, table, DELTA = 1.E-7):
    """ Compute finite difference given a parameter name """
    refls = self.scaler.ISIGI

    def target():
      r = self.compute_intensity_parameters()
      return refls['iobs'] / r['D']

    functional = target()

    if parameter_name.startswith('c'):
      from scitbx.matrix import sqr
      cryst_param = int(parameter_name.lstrip('c'))
      parameter_name = 'b_matrix'
      current = table[parameter_name]*1 # make a copy

      sre = symmetrize_reduce_enlarge(self.scaler.params.target_space_group.group())
      for i in xrange(len(table)):
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

    rx = flex.mat3_double()
    ry = flex.mat3_double()
    u = flex.mat3_double()
    b = flex.mat3_double()
    wavelength = flex.double()
    G = flex.double()
    B = flex.double()
    s0 = flex.vec3_double()
    deff = flex.double()
    eta = flex.double()

    ex = col((1,0,0))
    ey = col((0,1,0))

    for i in xrange(len(ct)):
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

    iobs       = refls['iobs']
    h          = refls['miller_index_original'].as_vec3_double()
    q          = ry * rx * u * b * h
    qlen       = q.norms()
    d          = 1/q.norms()
    #rs         = (1/deff)+(eta/(2*d))
    rs         = 1/deff # assumes eta is zero
    #eta        = flex.double(len(refls), 0)
    #deff       = 1/rs # assumes eta is zero
    s          = (s0+q)
    slen       = s.norms()
    rh         = slen-(1/wavelength)
    rs_sq      = rs*rs
    p_n        = rs_sq
    p_d        = (2. * (rh * rh)) + rs_sq
    partiality = p_n/p_d
    theta      = flex.asin(wavelength/(2*d))
    epsilon    = -8*B*(flex.sin(theta)/wavelength)**2
    eepsilon   = flex.exp(epsilon)
    D          = partiality * G * eepsilon
    thetah     = flex.asin(wavelength/(2*d))
    sinthetah  = flex.sin(thetah)
    er         = sinthetah/wavelength


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

  def adjust_errors(self):
    """ Propagate errors based on statistical error propagation and the observed population of
    images to the scaled and merged intensity errors """
    assert self.scaler.params.postrefinement.algorithm == 'rs'

    refls = self.scaler.ISIGI
    ct = self.scaler.crystal_table

    # this version doesn't post-refine deff and eta directly so compute these values from the refined rs value
    ct['deff'] = 1/ct['RS'] # assumes eta is zero
    ct['eta'] = flex.double(len(ct), 0)

    # Compute errors by examining distributions of parameters
    stats_thetax = flex.mean_and_variance(ct['thetax'])
    stats_thetay = flex.mean_and_variance(ct['thetay'])
    stats_lambda = flex.mean_and_variance(ct['wavelength'])
    #stats_eta    = flex.mean_and_variance(ct['ETA'])
    stats_deff   = flex.mean_and_variance(ct['deff'])
    stats_rs     = flex.mean_and_variance(ct['RS'])
    sigma_thetax = stats_thetax.unweighted_sample_standard_deviation()
    sigma_thetay = stats_thetay.unweighted_sample_standard_deviation()
    sigma_lambda = stats_lambda.unweighted_sample_standard_deviation()
    sigma_eta    = 0 #stats_eta.unweighted_sample_standard_deviation()
    sigma_deff   = stats_deff.unweighted_sample_standard_deviation()
    sigma_rs     = stats_rs.unweighted_sample_standard_deviation()
    print >> self.log, "ThetaX %.4f +/- %.4f"    %(r2d(stats_thetax.mean()), r2d(sigma_thetax))
    print >> self.log, "Thetay %.4f +/- %.4f"    %(r2d(stats_thetay.mean()), r2d(sigma_thetay))
    print >> self.log, "Wavelength %.4f +/- %.4f"%(    stats_lambda.mean(),      sigma_lambda)
    #print "ETA %.4f +/- %.4f"       %(    stats_eta.mean(),         sigma_eta)
    print >> self.log, "DEFF %.4f +/- %.4f"      %(    stats_deff.mean(),        sigma_deff)
    print >> self.log, "RS %.6f +/- %.6f"        %(    stats_rs.mean(),          sigma_rs)

    drx_dthetax = flex.mat3_double()
    dry_dthetay = flex.mat3_double()
    s0hat = flex.vec3_double(len(refls), (0,0,-1))

    ex = col((1,0,0))
    ey = col((0,1,0))

    # Compute derivatives
    sre = symmetrize_reduce_enlarge(self.scaler.params.target_space_group.group())
    c_gstar_params = None
    c_gstar_derivatives = None
    gstar_params = None
    gstar_derivatives = None

    for i in xrange(len(ct)):
      n_refl = ct['n_refl'][i]
      drx_dthetax.extend(flex.mat3_double(n_refl, ex.axis_and_angle_as_r3_derivative_wrt_angle(ct['thetax'][i])))
      dry_dthetay.extend(flex.mat3_double(n_refl, ey.axis_and_angle_as_r3_derivative_wrt_angle(ct['thetay'][i])))

      sre.set_orientation(orientation=ct['b_matrix'][i])
      p = sre.forward_independent_parameters()
      dB_dp = sre.forward_gradients()
      if gstar_params is None:
        assert gstar_derivatives is None and c_gstar_params is None and c_gstar_derivatives is None
        c_gstar_params = [flex.double() for j in xrange(len(p))]
        c_gstar_derivatives = [flex.mat3_double() for j in xrange(len(p))]
        gstar_params = [flex.double() for j in xrange(len(p))]
        gstar_derivatives = [flex.mat3_double() for j in xrange(len(p))]
      assert len(p) == len(dB_dp) == len(gstar_params) == len(gstar_derivatives) == len(c_gstar_params) == len(c_gstar_derivatives)
      for j in xrange(len(p)):
        c_gstar_params[j].append(p[j])
        c_gstar_derivatives[j].append(tuple(dB_dp[j]))
        gstar_params[j].extend(flex.double(n_refl, p[j]))
        gstar_derivatives[j].extend(flex.mat3_double(n_refl, tuple(dB_dp[j])))

    r = self.compute_intensity_parameters()

    print >> self.log, "Free G* parameters"
    sigma_gstar = []
    sI_dgstar = []
    for j in xrange(len(gstar_params)):
      stats  = flex.mean_and_variance(c_gstar_params[j])
      print >> self.log, "G* %d %.4f *1e-5 +/- %.4f *1e-5"%(j, stats.mean()*1e5, stats.unweighted_sample_standard_deviation()*1e5)
      sigma_gstar.append(stats.unweighted_sample_standard_deviation())

    sigma_Iobs = refls['scaled_intensity']/refls['isigi']
    dI_dIobs = 1/r['D']

    def compute_dI_dp(dq_dp):
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

    dI_dgstar = []
    for j in xrange(len(gstar_params)):
      dI_dgstar.append(compute_dI_dp(r['ry'] * r['rx'] * r['u'] * gstar_derivatives[j] * r['h']))

    dI_dthetax = compute_dI_dp(r['ry'] * drx_dthetax * r['u'] * r['b'] * r['h'])
    dI_dthetay = compute_dI_dp(dry_dthetay * r['rx'] * r['u'] * r['b'] * r['h'])

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

    drs_deff = -1/(r['deff']**2)
    dPn_deff = 2 * r['rs'] * drs_deff
    dPd_deff = 2 * r['rs'] * drs_deff
    dP_deff  = ((r['p_d'] * dPn_deff)-(r['p_n'] * dPd_deff))/(r['p_d']**2)
    dI_deff  = -(refls['iobs']/(r['partiality']**2 * r['G'] * r['eepsilon'])) * dP_deff

    drs_deta = 1/(2*r['d'])
    dPn_deta = 2 * r['rs'] * drs_deta
    dPd_deta = 2 * r['rs'] * drs_deta
    dP_deta  = ((r['p_d']*dPn_deta)-(r['p_n']*dPd_deta))/(r['p_d']**2)
    dI_deta  = -(refls['iobs']/(r['partiality']**2 * r['G'] * r['eepsilon'])) * dP_deta

    if True:
      # Finite differences
      n_cryst_params = sre.constraints.n_independent_params()
      print "Showing finite differences and derivatives for each parameter (first few reflections only)"
      for parameter_name, table, derivatives, delta, in zip(['iobs', 'thetax', 'thetay', 'wavelength', 'deff', 'eta'] + ['c%d'%cp for cp in xrange(n_cryst_params)],
                                                    [refls, ct, ct, ct, ct, ct] + [ct]*n_cryst_params,
                                                    [dI_dIobs, dI_dthetax, dI_dthetay, dI_dlambda, dI_deff, dI_deta] + dI_dgstar,
                                                    [1e-7]*6 + [1e-11]*n_cryst_params):
        finite_g = self.finite_difference(parameter_name, table, delta)
        print parameter_name
        for refl_id in xrange(min(10, len(refls))):
          print "%d % 21.1f % 21.1f"%(refl_id, finite_g[refl_id], derivatives[refl_id])
        stats = flex.mean_and_variance(finite_g-derivatives)
        stats_finite = flex.mean_and_variance(finite_g)
        percent = 0 if stats_finite.mean() == 0 else 100*stats.mean()/stats_finite.mean()
        print "Mean difference between finite and analytical: % 24.4f +/- % 24.4f (%8.3f%% of finite d.)"%( \
            stats.mean(), stats.unweighted_sample_standard_deviation(), percent)
        print

    # Propagate errors
    refls['isigi'] = refls['scaled_intensity'] / flex.sqrt(((sigma_Iobs**2 * dI_dIobs**2) +
                                                            sum([sigma_gstar[j]**2 * dI_dgstar[j]**2 for j in xrange(len(sigma_gstar))]) +
                                                            (sigma_thetax**2 * dI_dthetax**2) +
                                                            (sigma_thetay**2 * dI_dthetay**2) +
                                                            (sigma_lambda**2 * dI_dlambda**2) +
                                                            (sigma_deff**2 * dI_deff**2) +
                                                            (sigma_eta**2 * dI_deta**2)))

    # Show results of propagation
    from scitbx.math import five_number_summary
    all_data = [(refls['iobs'], "Iobs"),
                (sigma_Iobs, "Original errors"),
                (1/r['D'], "Total scale factor"),
                (refls['iobs']/r['D'], "Inflated intensities"),
                (refls['scaled_intensity']/refls['isigi'], "Propagated errors"),
                (flex.sqrt(sigma_Iobs**2 * dI_dIobs**2), "Iobs term"),
                (flex.sqrt(sigma_thetax**2 * dI_dthetax**2), "Thetax term"),
                (flex.sqrt(sigma_thetay**2 * dI_dthetay**2), "Thetay term"),
                (flex.sqrt(sigma_lambda**2 * dI_dlambda**2), "Wavelength term"),
                (flex.sqrt(sigma_deff**2 * dI_deff**2), "Deff term"),
                (flex.sqrt(sigma_eta**2 * dI_deta**2), "Eta term")] + \
               [(flex.sqrt(sigma_gstar[j]**2 * dI_dgstar[j]**2), "Gstar term %d"%j) for j in xrange(len(sigma_gstar))]
    print >> self.log, "%20s % 20s % 20s % 20s"%("Data name","Quartile 1", "Median", "Quartile 3")
    for data, title in all_data:
      fns = five_number_summary(data)
      print >> self.log, "%20s % 20d % 20d % 20d"%(title, fns[1], fns[2], fns[3])

    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    Intensity = refls['scaled_intensity']
    sigma = Intensity / refls['isigi']
    variance = sigma * sigma

    for i in xrange(len(refls)):
      j = refls['miller_id'][i]
      self.scaler.summed_wt_I[j] += Intensity[i] / variance[i]
      self.scaler.summed_weight[j] += 1 / variance[i]
