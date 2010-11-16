import math
from cctbx.array_family import flex
from cctbx.eltbx.xray_scattering import wk1995
from scitbx.math import bessel_i1_over_i0
from mmtbx import max_lik
import math
from cctbx.array_family import flex
from cctbx import miller
from libtbx import adopt_init_args
import iotbx.phil

alpha_beta_params = iotbx.phil.parse("""\
  free_reflections_per_bin = 140
    .type = int
  number_of_macromolecule_atoms_absent = 225
    .type = int
  n_atoms_included = 0
    .type = int
  bf_atoms_absent = 15.0
    .type = float
  final_error = 0.0
    .type = float
  absent_atom_type = "O"
    .type = str
  method = *est calc
    .type = choice
  estimation_algorithm = *analytical iterative
    .type = choice
  verbose = -1
    .type = int
  interpolation = True
    .type = bool
  fix_scale_for_calc_option = None
    .type = float
  number_of_waters_absent = 613
    .type = float
""")

def fo_fc_alpha_over_eps_beta(f_obs, f_model, alpha, beta):
  # Parameter "t". The coefficient 2 is already included, so for example
  # fom = th(t) or I1(t)/I0(t) depending on cf.
  return max_lik.fo_fc_alpha_over_eps_beta(
    f_obs          = f_obs.data(),
    f_model        = flex.abs(f_model.data()),
    alpha          = alpha.data(),
    beta           = beta.data(),
    space_group    = f_obs.space_group(),
    miller_indices = f_obs.indices())

def figures_of_merit_(f_obs,
                      f_calc,
                      alpha,
                      beta):
  return max_lik.fom_and_phase_error(f_obs          = f_obs.data(),
                                     f_model        = f_calc.data(),
                                     alpha          = alpha,
                                     beta           = beta,
                                     space_group    = f_obs.space_group(),
                                     miller_indices = f_obs.indices()).fom()

def phase_error(f_obs,
                f_calc,
                alpha,
                beta):
  return max_lik.fom_and_phase_error(
                                f_obs          = f_obs.data(),
                                f_model        = f_calc.data(),
                                alpha          = alpha,
                                beta           = beta,
                                space_group    = f_obs.space_group(),
                                miller_indices = f_obs.indices()).phase_error()

def fom(t,lcent):
    #
    # calculates the ratio I1(2t)/I0(2t) if lcent=0 and th(t) in other cases
    # I1() - Modified Bessel Function of 1 order
    # I0() - Modified Bessel Function of 0 order
    #
    if lcent == 0:
       result = bessel_i1_over_i0(2.*t)
    else:
       result = math.tanh(t)
    return result
###############################################################################
class alpha_beta_calc(object):

  def __init__(self,f,
                    n_atoms_absent,
                    n_atoms_included,
                    bf_atoms_absent,
                    final_error,
                    absent_atom_type):
    #
    #   ss=s**2 = (2*sin(teta)/lambda)**2 for the given reflection;
    #   final_error - desired mean error in atomic positions (in A);
    #                 it must be specified as 0., if the user has
    #                 no idea about its other value;
    #   n_atoms_included - an approximate number of non-hydrogen atoms in
    #                      the ASYMMETRIC PART OF THE UNIT CELL, which
    #                      are INCLUDED into the current model for refinement;
    #   n_atoms_absent - an approximate number of non-hydrogen atoms in
    #                    the ASYMMETRIC PART OF THE UNIT CELL, which are
    #                    NOT INCLUDED into the current model for refinement;
    #.....................................................................
    # P.Afonine, V.Lunin & A.Urzhumtsev.(2003).J.Appl.Cryst.36,158-159
    #
    self.f = f
    assert n_atoms_absent >= 0
    assert n_atoms_included >= 0
    assert f.size() > 0
    self.ss = 1./flex.pow2(f.d_spacings().data())
    assert self.ss.size() == f.data().size()
    self.nsym = f.space_group().order_z()
    assert self.nsym >= 1
    self.n_atoms_absent   = n_atoms_absent
    self.n_atoms_included = n_atoms_included
    self.bf_atoms_absent = bf_atoms_absent
    if final_error is None : final_error = 0.0
    self.final_error = final_error
    assert final_error >= 0.0
    if absent_atom_type is None : absent_atom_type="C"
    self.absent_atom_type = absent_atom_type

  def alpha_beta(self):
    #
    # alpha, beta by formulas
    #.........................................................................
    # V.Lunin & T.Skovoroda. Acta Cryst. (1995). A51, 880-887
    # A.Urzhumtsev, T.Skovoroda & V.Lunin. J. Appl. Cryst. (1996). 29, 741-744
    # V.Lunin, P.Afonine & A.Urzhumtsev. Acta Cryst. (2002). A58, 270-282
    #
    alpha=[] ; beta=[]
    n_part = self.nsym * self.n_atoms_included
    n_lost = self.nsym * self.n_atoms_absent
    for ssi in self.ss:
       ak = math.exp(-0.25*ssi*(self.final_error**2)*(math.pi**3))
       alpha.append( ak )
       fact = self.form_factor(ssi)*math.exp(-self.bf_atoms_absent/4.0*ssi)
       beta.append(((1.0-ak**2)*n_part+n_lost)*fact**2)
    alpha_data = flex.double(alpha)
    beta_data = flex.double(beta)
    alpha = miller.array(miller_set = self.f, data = alpha_data)
    beta  = miller.array(miller_set = self.f, data = beta_data)
    return alpha, beta

  def form_factor(self,ss):
    #
    # W & K form-factor of atom C
    #
    table=wk1995(self.absent_atom_type).fetch()
    a_wk=table.array_of_a()
    b_wk=table.array_of_b()
    c_wk=table.c()
    result_wk=c_wk
    for i in xrange(5):
       result_wk += a_wk[i]*math.exp(-b_wk[i]*ss/4.0)
    return result_wk

###############################################################################
class alpha_beta_est_manager(object):

  def __init__(self,f_obs,
                    f_calc,
                    free_reflections_per_bin,
                    flags,
                    interpolation):
    adopt_init_args(self, locals())
    #
    # icent - array contains 0 for acentric reflections and >0 integer
    #         for centric reflections
    # epsilon  - array contains the correction factors for intensity
    #         they are equal to how many times the transposed symmetry
    #         matrixes leaves the reciprocal space point at the same place
    # icont - array contains the control set; the value 1 means that
    #         this reflection will be used in the calculation of
    #         likelihood function;
    # f_calc - array of calculated magnitude values
    # f_obs  - array of experimental magnitude values
    # free_reflections_per_bin - minimal number of reflections in given
    #                    resolution zone used for alpha,beta calculation
    # V.Lunin & T.Skovoroda. Acta Cryst. (1995). A51, 880-887
    # P.Afonine, V.Lunin & A.Urzhumtsev.(2003).J.Appl.Cryst.36,158-159
    #
    assert len(self.flags) == self.f_obs.data().size()
    assert self.f_obs.data().size() == self.f_calc.data().size()
    assert self.f_calc.indices().all_eq(self.f_obs.indices()) == 1
    self.f_calc = abs(self.f_calc)
    if(self.flags.count(True) > 0):
      if free_reflections_per_bin > flags.count(True):
         self.free_reflections_per_bin = flags.count(True)
      self.f_obs_test  = self.f_obs.select(self.flags)
      self.f_calc_test = self.f_calc.select(self.flags)
    if(self.flags.count(True) == 0):
      self.f_obs_test  = self.f_obs.select(~self.flags)
      self.f_calc_test = self.f_calc.select(~self.flags)
    self.f_obs_test.setup_binner_counting_sorted(
      reflections_per_bin= self.free_reflections_per_bin)
    self.fo_test_sets = []
    self.fm_test_sets = []
    self.indices_sets = []
    for i_bin in self.f_obs_test.binner().range_used():
       sel = self.f_obs_test.binner().selection(i_bin)
       sel_f_obs_test = self.f_obs_test.select(sel)
       sel_f_calc_test = self.f_calc_test.select(sel)
       if(sel.count(True) > 0): # XXX I do not understand why it can be 0 (in rare cases)
         self.fo_test_sets.append(sel_f_obs_test.data())
         self.fm_test_sets.append(sel_f_calc_test.data())
         self.indices_sets.append(sel_f_obs_test.indices())
    for a,b,c in zip(self.fo_test_sets, self.fm_test_sets, self.indices_sets):
      assert a.size() == b.size() == c.size() != 0
    obj = max_lik.alpha_beta_est(fo_test     = self.fo_test_sets,
                                 fm_test     = self.fm_test_sets,
                                 indices     = self.indices_sets,
                                 space_group = self.f_obs_test.space_group())
    self.alpha_in_zones, self.beta_in_zones = obj.alpha(), obj.beta()
    alpha, beta = self.alpha_beta_for_each_reflection()
    self.alpha = miller.array(miller_set = self.f_obs, data = alpha)
    self.beta  = miller.array(miller_set = self.f_obs, data = beta)

  def alpha_beta(self):
    return self.alpha, self.beta

  def smooth(self,x):
    if len(x) > 1:
      x1=x[0]
      x2=x[1]
      for i in range(1,len(x)-1,1):
         x3=x[i+1]
         tmp = (x1+x2+x3)/3.0
         x[i]=tmp
         x1=x2
         x2=x3
      for i in range(0,len(x),1):
        if(x[i] < 0.01):
          try: x[i] = x[i-1]
          except: x[i] = x[i+1]

    return x

  def alpha_beta_for_each_reflection(self):
    alpha = flex.double(self.f_obs.size())
    beta = flex.double(self.f_obs.size())
    self.f_obs.setup_binner(n_bins= len(self.alpha_in_zones))
    binner = self.f_obs.binner()
    if(self.interpolation == True):
      az = flex.double(self.smooth(self.alpha_in_zones))
      bz = flex.double(self.smooth(self.beta_in_zones) )
      alpha = binner.interpolate(az, 0)
      beta  = binner.interpolate(bz, 0)
    elif(self.interpolation == False):
      for i_bin, az, bz in zip(binner.range_used(),self.alpha_in_zones,
                               self.beta_in_zones):
        sel = binner.selection(i_bin)
        alpha.set_selected(sel, az)
        beta.set_selected(sel, bz)
    return alpha, beta


class alpha_beta(object):
  def __init__(self, f_obs            = None,
                     f_calc           = None,
                     free_reflections_per_bin = None,
                     flags            = None,
                     verbose          = None,
                     n_atoms_absent   = None,
                     n_atoms_included = None,
                     bf_atoms_absent  = None,
                     final_error      = None,
                     absent_atom_type = None,
                     method           = None,
                     interpolation    = None):
    adopt_init_args(self, locals())
    assert self.method == "calc" or self.method == "est" or \
           self.method == "calc_and_est"
    assert self.verbose is not None
    if (self.method == "est"):
      assert self.interpolation is not None
      assert self.f_obs.data().size() == self.f_calc.data().size()
      assert self.flags.size() == self.f_obs.data().size()
      assert self.f_calc.indices().all_eq(self.f_obs.indices()) == 1
      self.alpha, self.beta = alpha_beta_est_manager(
        f_obs           = self.f_obs,
        f_calc          = abs(self.f_calc),
        free_reflections_per_bin = self.free_reflections_per_bin,
        flags           = self.flags,
        interpolation   = self.interpolation).alpha_beta()
    if (self.method == "calc"):
      assert self.f_obs is not None or self.f_calc is not None
      if (self.f_obs is not None): f = self.f_obs
      else: f = self.f_calc
      assert self.n_atoms_absent is not None
      assert self.n_atoms_included is not None
      assert self.bf_atoms_absent is not None
      assert self.absent_atom_type is not None
      self.alpha, self.beta = alpha_beta_calc(
                                    f                = f,
                                    n_atoms_absent   = self.n_atoms_absent,
                                    n_atoms_included = self.n_atoms_included,
                                    bf_atoms_absent  = self.bf_atoms_absent,
                                    final_error      = self.final_error,
                                    absent_atom_type = self.absent_atom_type).alpha_beta()
    if (self.method == "calc_and_est"):
      assert self.interpolation    is not None
      assert self.f_obs            is not None
      assert self.f_calc           is not None
      assert self.free_reflections_per_bin is not None
      assert self.flags            is not None
      assert self.n_atoms_absent   is not None
      assert self.n_atoms_included is not None
      assert self.bf_atoms_absent  is not None
      assert self.final_error      is not None
      assert self.absent_atom_type is not None
      assert self.f_obs.data().size() == self.f_calc.data().size() == \
             self.flags.size()
      assert self.f_calc.indices().all_eq(self.f_obs.indices()) == 1
      self.alpha_calc, self.beta_calc = alpha_beta_calc(
                                    f                = self.f_obs,
                                    n_atoms_absent   = self.n_atoms_absent,
                                    n_atoms_included = self.n_atoms_included,
                                    bf_atoms_absent  = self.bf_atoms_absent,
                                    final_error      = self.final_error,
                                    absent_atom_type = self.absent_atom_type).alpha_beta()
      self.alpha_est, self.beta_est = alpha_beta_est(
        f_obs           = self.f_obs,
        f_calc          = abs(self.f_calc),
        free_reflections_per_bin = self.free_reflections_per_bin,
        flags           = self.flags,
        interpolation   = self.interpolation).alpha_beta()
      alpha_calc_ma = miller.array(miller_set= self.f_obs,data= self.alpha_calc)
      beta_calc_ma  = miller.array(miller_set= self.f_obs,data= self.beta_calc)
      alpha_est_ma  = miller.array(miller_set= self.f_obs,data= self.alpha_est)
      beta_est_ma   = miller.array(miller_set= self.f_obs,data= self.beta_est)

      ss = 1./flex.pow2(alpha_calc_ma.d_spacings().data())
      omega_calc = []
      omega_est  = []
      for ac,ae,ssi in zip(self.alpha_calc, self.alpha_est,ss):
        if(ac > 1.0): ac = 1.0
        if(ae > 1.0): ae = 1.0
        if(ac <= 0.0): ac = 1.e-6
        if(ae <= 0.0): ae = 1.e-6
        coeff = -4./(math.pi**3*ssi)
        omega_calc.append( math.sqrt( math.log(ac) * coeff ) )
        omega_est.append( math.sqrt( math.log(ae) * coeff ) )
      omega_calc_ma = miller.array(miller_set= self.f_obs,data= flex.double(omega_calc))
      omega_est_ma  = miller.array(miller_set= self.f_obs,data= flex.double(omega_est))

      if(self.flags.count(True) > 0):
        omega_calc_ma_test = omega_calc_ma.select(self.flags)
        omega_est_ma_test  = omega_est_ma.select(self.flags)
        alpha_calc_ma_test= alpha_calc_ma.select(self.flags)
        beta_calc_ma_test= beta_calc_ma.select(self.flags)
        alpha_est_ma_test= alpha_est_ma.select(self.flags)
        beta_est_ma_test= beta_est_ma.select(self.flags)
      if(self.flags.count(True) == 0):
        omega_calc_ma_test = omega_calc_ma.select(~self.flags)
        omega_est_ma_test  = omega_est_ma.select(~self.flags)
        alpha_calc_ma_test= alpha_calc_ma.select(~self.flags)
        beta_calc_ma_test= beta_calc_ma.select(~self.flags)
        alpha_est_ma_test= alpha_est_ma.select(~self.flags)
        beta_est_ma_test= beta_est_ma.select(~self.flags)

      alpha_calc_ma_test.setup_binner(
        reflections_per_bin = self.free_reflections_per_bin)
      beta_calc_ma_test.use_binning_of(alpha_calc_ma_test)
      alpha_est_ma_test.use_binning_of(alpha_calc_ma_test)
      beta_est_ma_test.use_binning_of(alpha_calc_ma_test)
      omega_calc_ma_test.use_binning_of(alpha_calc_ma_test)
      omega_est_ma_test.use_binning_of(alpha_calc_ma_test)

      print "    Resolution           Estimated and calculated alpha, beta and model error"
      print "   d1        d2    nref alpha_e    beta_e    err_e alpha_c   beta_c     err_c"
      for i_bin in alpha_calc_ma_test.binner().range_used():
        sel = alpha_calc_ma_test.binner().selection(i_bin)
        sel_alpha_calc_ma_test = alpha_calc_ma_test.select(sel)
        sel_beta_calc_ma_test  = beta_calc_ma_test.select(sel)
        sel_alpha_est_ma_test  = alpha_est_ma_test.select(sel)
        sel_beta_est_ma_test   = beta_est_ma_test.select(sel)
        sel_omega_calc_ma_test = omega_calc_ma_test.select(sel)
        sel_omega_est_ma_test  = omega_est_ma_test.select(sel)
        size = sel_alpha_calc_ma_test.data().size()
        if(self.interpolation == False):
          i=0
          while i < size:
            v_alpha_est = sel_alpha_est_ma_test.data()[i]
            if(sel_alpha_est_ma_test.data().count(v_alpha_est) > size/2):
              v_beta_est = sel_beta_est_ma_test.data()[i]
              if(sel_beta_est_ma_test.data().count(v_beta_est) > size/2): break
            i+=1
        elif(self.interpolation == True):
          v_alpha_est = flex.mean(sel_alpha_est_ma_test.data())
          v_beta_est  = flex.mean(sel_beta_est_ma_test.data())
        alpha_calc = flex.mean(sel_alpha_calc_ma_test.data())
        beta_calc  = flex.mean(sel_beta_calc_ma_test.data())
        omega_c = flex.mean(sel_omega_calc_ma_test.data())
        omega_e = flex.mean(sel_omega_est_ma_test.data())
        d1 = alpha_calc_ma_test.binner().bin_d_range(i_bin)[0]
        d2 = alpha_calc_ma_test.binner().bin_d_range(i_bin)[1]
        print "%8.4f %8.4f %5d %6.5f %12.3f %5.3f %6.5f %12.3f %5.3f" % \
          (d1,d2,size,v_alpha_est,\
           v_beta_est,omega_e,alpha_calc,beta_calc,omega_c)

  def alpha_beta(self):
    return self.alpha, self.beta

def sigma_miss(miller_array, n_atoms_absent, bf_atoms_absent, absent_atom_type):
  result = flex.double()
  if(n_atoms_absent == 0): return flex.double(miller_array.indices().size(), 0)
  def form_factor(ssi, absent_atom_type):
    table=wk1995(absent_atom_type).fetch()
    a_wk=table.array_of_a()
    b_wk=table.array_of_b()
    c_wk=table.c()
    result_wk=c_wk
    for i in xrange(5):
      result_wk += a_wk[i]*math.exp(-b_wk[i]*ssi/4.0)
    return result_wk
  ss = 1./flex.pow2(miller_array.d_spacings().data())
  nsym = miller_array.space_group().order_z()
  #
  for ssi in ss:
     fact = form_factor(ssi, absent_atom_type)*math.exp(-bf_atoms_absent/4.0*ssi)
     result.append(fact * nsym * n_atoms_absent)
  return result
