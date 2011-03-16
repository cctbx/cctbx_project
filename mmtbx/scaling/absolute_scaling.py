## Peter Zwart, April 18, 2005
from cctbx.array_family import flex
from mmtbx import scaling
from cctbx import uctbx
from cctbx import adptbx
from cctbx.eltbx import xray_scattering
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from scitbx.linalg import eigensystem
from libtbx.utils import show_exception_info_if_full_testing
import scitbx.lbfgs
import math
import sys
from libtbx.str_utils import format_value


class gamma_protein:
  def __init__(self, d_star_sq):
    ## Coefficient for gamma as a fuynction of resolution
    ## described by a chebyshev polynome.
    coefs_mean = \
            [-0.24994838652402987,    0.15287426147680838,
              0.068108692925184011,   0.15780196907582875,
             -0.07811375753346686,    0.043211175909300889,
             -0.043407219965134192,   0.024613271516995903,
              0.0035146404613345932, -0.064118486637211411,
              0.10521875419321854,   -0.10153928782775833,
              0.0335706778430487,    -0.0066629477818811282,
             -0.0058221659481290031,  0.0136026246654981,
             -0.013385834361135244,   0.022526368996167032,
             -0.019843844247892727,   0.018128145323325774,
             -0.0091740188657759101,  0.0068283902389141915,
             -0.0060880807366142566,  0.0004002124110802677,
             -0.00065686973991185187,-0.0039358839200389316,
              0.0056185833386634149, -0.0075257168326962913,
             -0.0015215201587884459, -0.0036383549957990221,
             -0.0064289154284325831,  0.0059080442658917334,
             -0.0089851215734611887,  0.0036488156067441039,
             -0.0047375008148055706, -0.00090999496111171302,
              0.00096986728652170276,-0.0051006830761911011,
              0.0046838536228956777, -0.0031683076118337885,
              0.0037866523617167236,  0.0015810274077361975,
              0.0011030841357086191,  0.0015715596895281762,
             -0.0041354783162507788]

    coefs_mean = flex.double(coefs_mean)

    ## Coefficients descriubing the standard deviation of the above term
    coefs_sigma =\
            [ 0.040664777671929754,   0.015978463897495004,
              0.013068746157907499,  -0.00022301829884293671,
             -3.6158446516473002e-05,-0.0038078038613494048,
             -0.0016909798393289835, -0.003513220109509628,
             -0.00097404245360874664,-0.0037187631008071421,
              0.00031757305576918596,-0.0045169213215130082,
             -0.00086517495110945088,-0.0028285478264746034,
              0.00079199093295738952, 0.00062164093803723265,
              0.0018573920701968332,  0.0012407439272070957,
              0.0002221210281800356, -0.00026366883273910315,
             -0.0011200111831549365, -0.00083301832564164021,
             -0.0011493805850995207, -0.00083852924805887018,
             -0.00019326928638570406,-0.00039548849332659371,
             -0.00022423676327422079,-0.00093964924566378681,
             -0.00063394326307545398,-0.00033546289190621823,
              0.00040784160433128666, 0.001443822356327713,
              0.0013545273776501506,  0.0016735119787882374,
              0.0014189539521390023,  0.0013332965706491645,
              0.00087782212397582277, 0.00043943299411748545,
              0.00063740734290143163, 0.0007244345027539082,
              0.00099161845846209391, 0.00083124705069858489,
              0.00013275127763292434,-0.00023456844928215886,
             -0.001094278715119489]

    coefs_sigma = flex.double( coefs_sigma )

    self.low = 0.008
    self.high = 0.69
    self.gamma_cheb=chebyshev_polynome(
      45, self.low, self.high, coefs_mean)
    self.gamma_sigma_cheb = chebyshev_polynome(
      45, self.low, self.high, coefs_sigma)

    ## Make sure that the d_star_sq array
    ## does not have any elements that fall outside
    ## the allowed range (determined by the range
    ## on which the chebyshev polynome was computed;
    ## i.e: self.low<=d_star_sq<=self.high
    self.d_star_sq = d_star_sq.deep_copy()

    lower_then_low_limit = self.d_star_sq <= self.low
    self.d_star_sq = self.d_star_sq.set_selected(
      lower_then_low_limit, self.low)

    higher_then_high_limit = self.d_star_sq >= self.high
    self.d_star_sq = self.d_star_sq.set_selected(
      higher_then_high_limit, self.high)
    ## If everything is okai, this assertion should be
    ## passed withgout any problem
    assert ( flex.min(self.d_star_sq)>=self.low  )
    assert ( flex.max(self.d_star_sq)<=self.high )

    self.gamma = self.gamma_cheb.f( self.d_star_sq )
    self.sigma_gamma = self.gamma_sigma_cheb.f( self.d_star_sq )

class gamma_nucleic:
  def __init__(self, d_star_sq):
    ## Coefficient for gamma as a fuynction of resolution
    ## described by a chebyshev polynome.
    coefs_mean = \
            [-0.30244055359995714,    0.14551540259035137,
              0.06404885418364728,    0.14126884888694674,
             -0.076010848339056358,   0.10445808072140797,
             -0.13817185173869803,    0.03162832786042917,
              0.041599262771740309,  -0.088562816354662108,
              0.063708058411421117,  -0.044796037515868393,
              0.0088130627883259444,  0.02692514601906355,
             -0.050931900034622016,   0.019443590649642444,
             -0.0011195556039252301, -0.01644343506476071,
              0.0065957064914017524, -0.018596655500261718,
              0.0096270346410321905,  0.016307576063048754,
             -0.02680646640174009,    0.020734177937708331,
              0.0028123353629064345,  0.0045005299107411609,
              0.0076229925628943053, -0.008362403313550976,
              0.0034163962268388306,  0.001904748909797396,
             -0.013325099196913409,   0.0048138529463141863,
              0.0037576434237086738, -0.011440938719878148,
              0.0070463203562045043, -0.014417892444775739,
              0.00051623208479814814,-0.007030834594537072,
             -0.010592510032603445,   0.0099794223029419579,
             -0.0042803299088959388,  0.0018056147902035455,
              9.1732385471614747e-05, 0.0048087303990040917,
              0.0033924291685209331]

    coefs_mean = flex.double(coefs_mean)

    ## Coefficients descriubing the standard deviation of the above term
    coefs_sigma =\
            [ 0.13066942262051989,    0.02540993472514427,
              0.022640055258519923,   0.010682155584811278,
              0.0055933901688389942,  0.010202224633747257,
             -0.0068126652213876008,  0.00074050873381034524,
              0.0043056775404382332, -0.0068162235999210587,
              0.0043883564931143154, -0.0046223069963272981,
             -0.0021799388224634842,  0.0018700994720378869,
             -0.0051883494911385414, -0.00036639670195559728,
             -0.0018731351222522098, -0.005641953724585742,
             -0.0021296034270177015, -0.0037654091933662288,
             -0.0031915331246228089,  0.0017569392295630887,
             -0.0023581953932665491,  0.0043374380859762538,
              0.003490459547329672,   0.0030620317182512053,
              0.0037626939824912907, -0.0014184248052271247,
             -0.0032475452005936508, -0.0053177954201788511,
             -0.0085157840734136816, -0.0057322608003856712,
             -0.0051182987317167803, -0.0052003177422633084,
             -0.001085721076048506,  -0.00072459199543249329,
              0.0010209328663554243,  0.00076695099249463397,
              0.00034115347063572426, 0.0021264541997130233,
             -0.00031955842674212867,-0.00148958769833968,
              0.0003181991857060145, -0.00069586514533741132,
             -0.00046211335387235546]

    coefs_sigma = flex.double( coefs_sigma )

    self.low = 0.008
    self.high = 0.69
    self.gamma_cheb=chebyshev_polynome(
      45, self.low, self.high, coefs_mean)
    self.gamma_sigma_cheb = chebyshev_polynome(
      45, self.low, self.high, coefs_sigma)

    ## Make sure that the d_star_sq array
    ## does not have any elements that fall outside
    ## the allowed range (determined by the range
    ## on which the chebyshev polynome was computed;
    ## i.e: self.low<=d_star_sq<=self.high
    self.d_star_sq = d_star_sq.deep_copy()
    lower_then_low_limit = self.d_star_sq <= self.low
    self.d_star_sq = self.d_star_sq.set_selected(
      lower_then_low_limit, self.low)

    higher_then_high_limit = self.d_star_sq >= self.high
    self.d_star_sq = self.d_star_sq.set_selected(
      higher_then_high_limit, self.high)
    ## If everything is okai, this assertion should be
    ## passed withgout any problem
    assert ( flex.min(self.d_star_sq)>=self.low  )
    assert ( flex.max(self.d_star_sq)<=self.high )

    self.gamma = self.gamma_cheb.f( self.d_star_sq )
    self.sigma_gamma = self.gamma_sigma_cheb.f( self.d_star_sq )


## This class computes the expected intensity
## for a given d_star_sq_array given some basic info
class expected_intensity:
  def __init__(self,
               scattering_info,
               d_star_sq_array,
               p_scale=0.0,
               b_wilson=0.0,
               magic_fudge_factor=2.0):
    ## First recompute some parameters
    scattering_info.scat_data(d_star_sq_array)
    ## mean intensity
    self.mean_intensity = scattering_info.sigma_tot_sq
    self.mean_intensity = self.mean_intensity*(1.0+scattering_info.gamma_tot)
    ## I am missing a factor 2 somewhere
    self.mean_intensity/=magic_fudge_factor
    self.mean_intensity=self.mean_intensity*flex.exp(
      -d_star_sq_array*b_wilson/2.0)
    self.mean_intensity*=math.exp(-p_scale)

    ## the associated standard deviation
    self.sigma_intensity = scattering_info.gamma_tot_sigma
    self.sigma_intensity = scattering_info.sigma_tot_sq*self.sigma_intensity
    self.sigma_intensity = self.sigma_intensity*flex.exp(
      -d_star_sq_array*b_wilson/2.0)
    self.sigma_intensity*= math.exp(-p_scale)




class scattering_information(object):
  def __init__(self,
               n_residues=None,
               n_bases=None,
               asu_contents=None,
               fraction_protein=None,
               fraction_nucleic=None):
    """ Returns scattering info for specified structure"""
    ## Preference is given to a supplied asu_content dictionairy
    ## If this is not available, one will made up given the
    ## number of residues or bases.

    self.sigma_tot_sq = None
    self.gamma_tot = None
    self.gamma_tot_sigma = None

    ## if aus_contents is not none no resiudes should be specified
    if asu_contents is not None:
      assert (n_residues==None)
      assert (n_bases==None)
    ## and vise versa
    if ( (n_residues is not None) or (n_bases is not None) ):
      assert (asu_contents==None)
    ## if aus_contents is not none, fractions need to bve specified
    if asu_contents is not None:
      assert( fraction_protein is not None)
      assert( fraction_nucleic is not None)

    ## determin the fractions first
    if asu_contents is None:
      if n_residues is None:
        n_residues=0.0
      fraction_protein = (8.0*1.0*1.0+
                          5.0*6.0*6.0+
                          1.5*7.0*7.0+
                          1.2*8.0*8.0)*float(n_residues)
      if n_bases is None:
        n_bases=0.0
      fraction_nucleic = (16.0*1.0*1.0+
                          9.7*6.0*6.0+
                          3.8*7.0*7.0+
                          5.9*8.0*8.0+
                          1.0*15.0*15.0)*float(n_bases)
      tot = fraction_protein+fraction_nucleic
      fraction_protein/=tot
      fraction_nucleic/=tot


    if fraction_protein is not None:
      self.fraction_protein = fraction_protein
    else:
      self.fraction_protein = 0.0

    if fraction_nucleic is not None:
      self.fraction_nucleic = fraction_nucleic
    else:
      self.fraction_nucleic = 0.0
    assert ( self.fraction_protein+self.fraction_nucleic<=1.0 )


    if asu_contents is None:
      asu_contents = None
      if n_residues > 0:
        asu_contents = {"H":8.0*float(n_residues),
                        "C":5.0*float(n_residues),
                        "N":1.5*float(n_residues),
                        "O":1.2*float(n_residues)
                       }
      if n_bases > 0:
        ## These values are rather approximate
        asu_contents = {"H":float(n_bases)*16.0,
                        "C":float(n_bases)*9.7,
                        "N":float(n_bases)*3.8,
                        "O":float(n_bases)*5.9,
                        "P":float(n_bases)*1.0
                        }
      if n_bases > 0:
        if n_residues > 0:
          asu_contents = {"H":float(n_bases)*16.0 + 8.0*float(n_residues),
                          "C":float(n_bases)*9.7  + 5.0*float(n_residues),
                          "N":float(n_bases)*3.8  + 1.5*float(n_residues),
                          "O":float(n_bases)*5.9  + 1.2*float(n_residues),
                          "P":float(n_bases)*1.0
                          }
    self.asu_contents =  asu_contents


  def scat_data(self, d_star_sq=None):

    if d_star_sq is None:
      self.sigma_tot_sq=None
      self.gamma_tot_sigma=None
      self.gamma_tot=None

    if d_star_sq is not None:
      self.sigma_tot_sq = flex.double( d_star_sq.size() )
      gaussians = {}
      for chemical_type, n_atoms in self.asu_contents.items():
        gaussians[chemical_type] = xray_scattering.wk1995(
          chemical_type).fetch()
        f0 = gaussians[chemical_type].at_d_star_sq(d_star_sq)
        self.sigma_tot_sq += f0*f0*n_atoms

      if(d_star_sq.size()>0):
        ## Protein part
        gamma_prot = gamma_protein(d_star_sq)
        self.gamma_prot = gamma_prot.gamma*self.fraction_protein
        ## Nucleotide part; needs to be completed
        gamma_nuc = gamma_nucleic(d_star_sq)
        self.gamma_nuc = gamma_nuc.gamma*self.fraction_nucleic ##
        ## Totals
        self.gamma_tot = self.gamma_prot*self.fraction_protein +\
                         self.gamma_nuc*self.fraction_nucleic
        self.gamma_tot_sigma = (gamma_prot.sigma_gamma*self.fraction_protein)*\
                               (gamma_prot.sigma_gamma*self.fraction_protein)+\
                               (gamma_nuc.sigma_gamma*self.fraction_nucleic)*\
                               (gamma_nuc.sigma_gamma*self.fraction_nucleic)
        self.gamma_tot_sigma = flex.sqrt(  self.gamma_tot_sigma )




def anisotropic_correction(cache_0,
                           p_scale,
                           u_star,
                           b_add=None):
  ## Make sure that u_star is not rwgk scaled, i.e. like you get it from
  ## the ml_absolute_scale_aniso routine (!which is !!NOT!! scaled!)
  work_array = None
  try:
    work_array = cache_0.input.select( cache_0.input.data() > 0 )
  except KeyboardInterrupt: raise
  except: pass
  if work_array is None:
    work_array = cache_0.select( cache_0.data() > 0 )

  change_back_to_intensity=False
  if work_array.is_xray_intensity_array():
    work_array = work_array.f_sq_as_f()
    change_back_to_intensity=True

  assert not work_array.is_xray_intensity_array()

  if b_add is not None:
    u_star_add =  adptbx.b_iso_as_u_star( work_array.unit_cell(),
                                          b_add )
    u_star = u_star+u_star_add



  corrected_amplitudes = scaling.ml_normalise_aniso( work_array.indices(),
                                                     work_array.data(),
                                                     p_scale,
                                                     work_array.unit_cell(),
                                                     u_star )
  if work_array.sigmas() is not None:
    corrected_sigmas = scaling.ml_normalise_aniso( work_array.indices(),
                                                   work_array.sigmas(),
                                                   p_scale,
                                                   work_array.unit_cell(),
                                                   u_star )
  else:
    corrected_sigmas = None


  work_array = work_array.customized_copy(
    data = corrected_amplitudes,
    sigmas = corrected_sigmas ).set_observation_type(work_array)
  if change_back_to_intensity:
    work_array = work_array.f_as_f_sq()
  return work_array



class ml_iso_absolute_scaling(object):
  def __init__(self,
              miller_array,
              n_residues=None,
              n_bases=None,
              asu_contents=None,
              prot_frac = 1.0,
              nuc_frac= 0.0):
    """ Maximum likelihood isotropic wilson scaling.
    input:  a miller array (miller_array)
            number of protein residues (n_residues)
            number of nucleotides (n_bases)
               OR
            a scattering dictionairy (asu_contents)
              ( i.e. {'Au':1, 'C':2.5, 'O':1', 'H':3 } )
            The keys prot_frac and nuc_frac will be
            will be used at a later stage.
     """
    self.p_scale, self.b_wilson = None, None
    ## Checking input combinations
    if (n_residues is None):
      if (n_bases is None):
        assert asu_contents is not None
        assert (type(asu_contents) == type({}) )
    if asu_contents is None:
      assert ( (n_residues is not None) or (n_bases is not None) )

    assert (prot_frac+nuc_frac<=1.0)
    if ( miller_array.is_xray_intensity_array() ):
      miller_array = miller_array.f_sq_as_f()
    if ( miller_array.is_real_array() ):
      ## Save the information of the file name etc
      self.info = miller_array.info()
      work_array = miller_array.resolution_filter(
        1.0/math.sqrt(  scaling.get_d_star_sq_low_limit() ),
        1.0/math.sqrt( scaling.get_d_star_sq_high_limit() )
        )
      if work_array.data().size() > 0:
        work_array = work_array.select(work_array.data()>0)

        self.d_star_sq = flex.double(work_array.d_spacings().data()*
                                     work_array.d_spacings().data())
        self.d_star_sq = 1.0/self.d_star_sq
        self.scat_info =  None
        if asu_contents is None:
          self.scat_info= scattering_information(
                                            n_residues=n_residues,
                                            n_bases = n_bases,
                                            fraction_protein = prot_frac,
                                            fraction_nucleic = nuc_frac)
        else:
          self.scat_info = scattering_information(
                                             asu_contents = asu_contents,
                                             fraction_protein = prot_frac,
                                             fraction_nucleic = nuc_frac)
        if (work_array.size() > 0 ):
          ## Compute the terms
          self.scat_info.scat_data(self.d_star_sq)
          self.f_obs = work_array.data()
          ## Make sure sigma's are used when available
          if (work_array.sigmas() is not None):
            self.sigma_f_obs = work_array.sigmas()
          else:
            self.sigma_f_obs = flex.double(self.f_obs.size(),0.0)
          if (flex.min( self.sigma_f_obs ) < 0):
            self.sigma_f_obs = self.sigma_f_obs*0.0
          ## multiplicities and d_star_sq
          self.epsilon = work_array.epsilons().data().as_double()
          ## centric flags
          self.centric = flex.bool(work_array.centric_flags().data())
          ## Wilson parameters come from scattering_information class
          self.gamma_prot = self.scat_info.gamma_tot
          self.sigma_prot_sq = self.scat_info.sigma_tot_sq
          ## Optimisation stuff
          self.x = flex.double(2,0.0)
          self.x[0]=0.0
          self.x[1]=50.0
          self.f=0
          term_parameters = scitbx.lbfgs.termination_parameters( max_iterations = 1e6 ) # just for safety
          self.minimizer = scitbx.lbfgs.run(target_evaluator=self, termination_params=term_parameters)
          self.p_scale = self.x[0]
          self.b_wilson = self.x[1]
          ## this we do not need anymore
          del self.x
          del self.f_obs
          del self.sigma_f_obs
          del self.epsilon
          del self.gamma_prot
          del self.sigma_prot_sq
          del self.d_star_sq
          del self.centric

  def compute_functional_and_gradients(self):

    f = scaling.wilson_total_nll(self.d_star_sq,
                                 self.f_obs,
                                 self.sigma_f_obs,
                                 self.epsilon,
                                 self.sigma_prot_sq,
                                 self.gamma_prot,
                                 self.centric,
                                 self.x[0],
                                 self.x[1])
    g = flex.double( scaling.wilson_total_nll_gradient(self.d_star_sq,
                                          self.f_obs,
                                          self.sigma_f_obs,
                                          self.epsilon,
                                          self.sigma_prot_sq,
                                          self.gamma_prot,
                                          self.centric,
                                          self.x[0],
                                          self.x[1]) )
    self.f = f
    return f, g

  def show(self,
           out=None,
           verbose=0):
    if out is None:
      out = sys.stdout
    if verbose>0:
      print >> out, "ML estimate of overall B value of %s:" \
            % str(self.info)
      print >> out, "%s"%format_value("%5.2f", self.b_wilson), "A**(-2)"
      print >> out, "Estimated -log of scale factor of %s:" \
            % str(self.info)
      print >> out, "%s"%format_value("%5.2f", self.p_scale)




class ml_aniso_absolute_scaling(object):
  def __init__(self,
               miller_array,
               n_residues=None,
               n_bases=None,
               asu_contents=None,
               prot_frac = 1.0,
               nuc_frac= 0.0):
    """ Maximum likelihood anisotropic wilson scaling"""
    #Checking input
    if (n_residues is None):
      if (n_bases is None):
        assert asu_contents is not None
        assert (type(asu_contents) == type({}) )
    if asu_contents is None:
      assert ( (n_residues is not None) or (n_bases is not None) )
    assert (prot_frac+nuc_frac<=1.0)

    if ( miller_array.is_xray_intensity_array() ):
      miller_array = miller_array.f_sq_as_f()
    if ( miller_array.is_real_array() ):
      self.info = miller_array.info()

      work_array = miller_array.resolution_filter(
        1.0/math.sqrt(  scaling.get_d_star_sq_low_limit() ),
        1.0/math.sqrt( scaling.get_d_star_sq_high_limit() )
        )
      work_array = work_array.select(work_array.data()>0)

      self.d_star_sq = flex.double(work_array.d_spacings().data()*
                                   work_array.d_spacings().data())
      self.d_star_sq = 1.0/self.d_star_sq
      self.scat_info =  None
      if asu_contents is None:
        self.scat_info= scattering_information(
                                          n_residues=n_residues,
                                          n_bases = n_bases)
      else:
        self.scat_info = scattering_information(
                                           asu_contents = asu_contents,
                                           fraction_protein = prot_frac,
                                           fraction_nucleic = nuc_frac)
      self.scat_info.scat_data(self.d_star_sq)

      if (work_array.size() > 0 ):
        self.hkl = work_array.indices()
        self.f_obs = work_array.data()
        self.unit_cell =  uctbx.unit_cell(
          miller_array.unit_cell().parameters() )
        ## Make sure sigma's are used when available
        if (work_array.sigmas() is not None):
          self.sigma_f_obs = work_array.sigmas()
        else:
          self.sigma_f_obs = flex.double(self.f_obs.size(),0.0)
        if (flex.min( self.sigma_f_obs ) < 0):
          self.sigma_f_obs = self.sigma_f_obs*0.0

        ## multiplicities
        self.epsilon = work_array.epsilons().data().as_double()
        ## Determine Wilson parameters
        self.gamma_prot = self.scat_info.gamma_tot
        self.sigma_prot_sq = self.scat_info.sigma_tot_sq
        ## centric flags
        self.centric = flex.bool(work_array.centric_flags().data())
        ## Symmetry stuff
        self.sg = work_array.space_group()
        self.adp_constraints = self.sg.adp_constraints()
        self.dim_u = self.adp_constraints.n_independent_params
        ## Setup number of parameters
        assert self.dim_u()<=6
        ## Optimisation stuff
        self.x = flex.double(self.dim_u()+1, 0.0) ## B-values and scale factor!
        exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = False,
          ignore_line_search_failed_step_at_upper_bound = False,
          ignore_line_search_failed_maxfev              = False)
        term_parameters = scitbx.lbfgs.termination_parameters(
          max_iterations = 50)

        self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
          termination_params=term_parameters,
          exception_handling_params=exception_handling_params)

        ## Done refining
        Vrwgk = math.pow(self.unit_cell.volume(),2.0/3.0)
        self.p_scale = self.x[0]
        self.u_star = self.unpack()
        self.u_star = list( flex.double(self.u_star) / Vrwgk )
        self.b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(self.unit_cell,
                                         self.u_star))
        self.u_cif = adptbx.u_star_as_u_cif(self.unit_cell,
                                            self.u_star)
        #get eigenvalues of B-cart
        eigen = eigensystem.real_symmetric( self.b_cart )
        self.eigen_values = eigen.values()
        self.eigen_vectors = eigen.vectors()

        self.work_array  = work_array # i need this for further analyses

        del self.x
        del self.f_obs
        del self.sigma_f_obs
        del self.epsilon
        del self.gamma_prot
        del self.sigma_prot_sq
        del self.centric
        del self.hkl
        del self.d_star_sq

  def pack(self,g):
    ## generate a set of reduced parameters for the minimizer
    g_independent = [g[0]]
    g_independent = g_independent + \
      list(
       self.adp_constraints.independent_gradients(list(g[1:]))
      )
    return flex.double(g_independent)

  def unpack(self):
    ## generate all parameters from the reduced set
    ## for target function computing and so forth
    u_star_full = self.adp_constraints.all_params(list(self.x[1:]))
    return u_star_full


  def compute_functional_and_gradients(self):
    u = self.unpack()
    f = scaling.wilson_total_nll_aniso(self.hkl,
                                       self.f_obs,
                                       self.sigma_f_obs,
                                       self.epsilon,
                                       self.sigma_prot_sq,
                                       self.gamma_prot,
                                       self.centric,
                                       self.x[0],
                                       self.unit_cell,
                                       u)
    self.f=f
    g_full_exact = flex.double( scaling.wilson_total_nll_aniso_gradient(
      self.hkl,
      self.f_obs,
      self.sigma_f_obs,
      self.epsilon,
      self.sigma_prot_sq,
      self.gamma_prot,
      self.centric,
      self.x[0],
      self.unit_cell,
      u ))

    g = self.pack(g_full_exact)
    return f, g

  def format_it(self,x,format="%3.2f"):
    xx = format%(x)
    if x > 0:
      xx = " "+xx
    return(xx)


  def aniso_ratio_p_value(self,rat):
    return -3
    coefs = flex.double( [-1.7647171873040273, -3.4427008004789115, -1.097150249786379, 0.17303317520973829, 0.35955513268118661, 0.066276397961476205, -0.064575726062529232, -0.0063025873711609016, 0.0749945566688624, 0.14803702885155121, 0.154284467861286])
    fit_e = scitbx.math.chebyshev_polynome(11,0,1.0,coefs)
    x = flex.double( range(1000) )/999.0
    start = int(rat*1000)
    norma = flex.sum(flex.exp(fit_e.f(x)))/x[1]
    x = x*(1-rat)+rat
    norma2 = flex.sum(flex.exp(fit_e.f(x)))/(x[1]-x[0])
    return -math.log(norma2/norma )


  def analyze_aniso_correction(self, n_check=2000, p_check=0.25, level=3, z_level=9):

    correction_factors = self.work_array.customized_copy(
                 data=self.work_array.data()*0.0+1.0, sigmas=None )

    correction_factors = anisotropic_correction(
      correction_factors,0.0,self.u_star ).data()

    self.work_array = self.work_array.f_as_f_sq()

    isigi = self.work_array.data() / (
              self.work_array.sigmas()+max(1e-8,flex.min(self.work_array.sigmas()))
            )
    d_spacings = self.work_array.d_spacings().data().as_double()

    d_sort   = flex.sort_permutation( d_spacings )
    d_select = d_sort[0:n_check]

    min_d = d_spacings[ d_select[0] ]
    max_d = d_spacings[ d_select[ n_check-1] ]

    isigi = isigi.select( d_select )
    mean_isigi = flex.mean( isigi )
    observed_count = flex.bool( isigi > level ).as_double()
    mean_count = flex.mean( observed_count )
    correction_factors = correction_factors.select( d_select )

    isigi_rank      = flex.sort_permutation(isigi)
    correction_rank = flex.sort_permutation(correction_factors, reverse=True)

    n_again = int(correction_rank.size()*p_check )
    sel_hc = correction_rank[0:n_again]
    sel_lc = correction_rank[n_again:]
    mean_isigi_low_correction_factor  = flex.mean( isigi.select(sel_lc) )
    mean_isigi_high_correction_factor = flex.mean( isigi.select(sel_hc) )
    frac_below_low_correction         = flex.mean( observed_count.select(sel_lc) )
    frac_below_high_correction        = flex.mean( observed_count.select(sel_hc) )
    mu = flex.mean( observed_count )
    var = math.sqrt(mu*(1.0-mu)/n_again)
    z_low  = abs(frac_below_low_correction-mean_count)/max(1e-8,var)
    z_high = abs(frac_below_high_correction-mean_count)/max(1e-8,var)
    z_tot  = math.sqrt( (z_low*z_low + z_high*z_high) )

    message = """indicates that there probably is no significant systematic
noise amplification."""
    if z_tot > z_level:
      if mean_isigi_high_correction_factor < level:
        message =  """indicates that there probably is significant
systematic noise amplification that could possibly lead to artefacts in the
maps or difficulties in refinement"""
      else:
        message =  """indicates that there probably is some systematic dependence
between the anisotropy and not-so-well-defined  intensities. Because the signal
to noise for the most affected intensities is relatively good, the affect on maps
or refinement behavoir is most likely not very serious."""

    txt = """
For the resolution shell spanning between %4.2f - %4.2f Angstrom,
the mean I/sigI is equal to %5.2f. %4.1f %% of these intensities have
an I/sigI > 3. When sorting these intensities by their anisotropic
correction factor and analysing the I/sigI behavoir for this ordered
list, we can gauge the presence of 'anisotropy induced noise amplification'
in reciprocal space.

  The quarter of Intensities *least* affected by the anisotropy correction show
    <I/sigI>                 :   %5.2e
    Fraction of I/sigI > 3   :   %5.2e     ( Z = %8.2f )

  The quarter of Intensities *most* affected by the anisotropy correction show
    <I/sigI>                 :   %5.2e
    Fraction of I/sigI > 3   :   %5.2e     ( Z = %8.2f )

The combined Z-score of %8.2f %s

Z-scores are computed on the basis of a Bernoulli model assuming independence of weak reflections wrst anisotropy.

    """%(max_d, min_d,
         mean_isigi, 100.0*mean_count,
         mean_isigi_low_correction_factor, frac_below_low_correction,  z_low,
         mean_isigi_high_correction_factor, frac_below_high_correction, z_high,
         z_tot, message
        )


    #assert 3 ==0
    return txt










  def show(self,
           out=None,
           verbose=1):
    if out is None:
      out = sys.stdout
    try: b_cart = self.b_cart
    except AttributeError, e:
      print >> out, "*** ERROR ***"
      print >> out, str(e)
      show_exception_info_if_full_testing()
      return

    if verbose>0:

      print >> out
      print >> out
      print >> out,"ML estimate of overall B_cart value of %s:" \
            % str(self.info)
      # XXX: for GUI
      self.overall_b_cart = """\
%5.2f, %5.2f, %5.2f
%12.2f, %5.2f
%19.2f
""" % (self.b_cart[0], self.b_cart[3], self.b_cart[4],
                       self.b_cart[1], self.b_cart[5],
                                       self.b_cart[2])
      #print >> out,"%5.2f," %(self.b_cart[0]), "%5.2f,"\
      #      %(self.b_cart[3]),"%5.2f" %(self.b_cart[4])
      #print >> out,"%12.2f," %(self.b_cart[1]),"%5.2f" %(self.b_cart[5])
      #print >> out,"%19.2f"  %(self.b_cart[2])
      print >> out, self.overall_b_cart

      print >> out,"Equivalent representation as U_cif: "
      print >> out,"%5.2f," %(self.u_cif[0]), "%5.2f," \
            %(self.u_cif[3]), "%5.2f" %(self.u_cif[4])
      print >> out,"%12.2f," %(self.u_cif[1]),"%5.2f" %(self.u_cif[5])
      print >> out,"%19.2f"  %(self.u_cif[2])
      print >> out
      print >> out, "Eigen analyses of B-cart:"
      print >> out, "                 Value   Vector"
      print >> out, "Eigenvector 1 : %s  (%s, %s, %s)"%(self.format_it(self.eigen_values[0],"%5.3f"),
                                                        self.format_it(self.eigen_vectors[0]),
                                                        self.format_it(self.eigen_vectors[1]),
                                                        self.format_it(self.eigen_vectors[2]))
      print >> out, "Eigenvector 2 : %s  (%s, %s, %s)"%(self.format_it(self.eigen_values[1],"%5.3f"),
                                                        self.format_it(self.eigen_vectors[3]),
                                                        self.format_it(self.eigen_vectors[4]),
                                                        self.format_it(self.eigen_vectors[5]))
      print >> out, "Eigenvector 3 : %s  (%s, %s, %s)"%(self.format_it(self.eigen_values[2],"%5.3f"),
                                                        self.format_it(self.eigen_vectors[6]),
                                                        self.format_it(self.eigen_vectors[7]),
                                                        self.format_it(self.eigen_vectors[8]))
      print >> out
      print >> out,"ML estimate of  -log of scale factor of %s:" \
             % str(self.info)
      print >> out,"%5.2f" %(self.p_scale)

      print >> out
      print >> out
      print >> out, "----------------    Anisotropicity analyses     ----------------"
      print >> out
      anirat = abs(self.eigen_values[0]-self.eigen_values[2])/self.eigen_values[0]
      self.anirat = anirat
      ani_rat_p = self.aniso_ratio_p_value(anirat)
      if ani_rat_p < 0:
        ani_rat_p = 0.0
      print >> out, "Anisotropicity   ( [MaxAnisoB-MinAnisoB]/[MaxAnisoB] ) :  %7.3e"%(anirat)
      print >> out, "                             Anisotropic ratio p-value :  %7.3e"%(ani_rat_p)
      print >> out
      print >> out, "     The p-value is a measure of of the severity of anisotropy as observed in the PDB."
      print >> out, "     The p-value of %5.3e indicates that roughly %4.1f %% of dataset available in the PDB have"%(ani_rat_p,100.0*math.exp(-ani_rat_p))
      print >> out, "     an anisotropy equal or worse as compared to this dataset."
      print >> out

      print >> out, self.analyze_aniso_correction()
      print >> out




class kernel_normalisation(object):
  def __init__(self,
               miller_array,
               kernel_width=None,
               n_bins=23,
               n_term=13,
               d_star_sq_low=None,
               d_star_sq_high=None,
               auto_kernel=False,
               number_of_sorted_reflections_for_auto_kernel=50):
    ## Autokernel is either False, true or a specific integer
    if kernel_width is None:
      assert (auto_kernel is not False)
    if auto_kernel is not False:
      assert (kernel_width==None)

    assert miller_array.size()>0


    ## intensity arrays please
    work_array = None
    if not miller_array.is_real_array():
      raise RuntimeError("Please provide real arrays only")
      ## I might have to change this upper condition

    if miller_array.is_xray_amplitude_array():
      work_array = miller_array.f_as_f_sq()
    if miller_array.is_xray_intensity_array():
      work_array = miller_array.deep_copy()
      work_array = work_array.set_observation_type(miller_array)
    ## If type is not intensity or amplitude
    ## raise an execption please
    if not miller_array.is_xray_intensity_array():
      if not miller_array.is_xray_amplitude_array():
        raise RuntimeError("Observation type unknown")

    ## declare some shorthands
    I_obs = work_array.data()
    epsilons = work_array.epsilons().data().as_double()
    d_star_sq_hkl = work_array.d_spacings().data()
    d_star_sq_hkl = 1.0/(d_star_sq_hkl*d_star_sq_hkl)
    ## Set up some limits
    if d_star_sq_low is None:
      d_star_sq_low = flex.min(d_star_sq_hkl)

    if d_star_sq_high is None:
      d_star_sq_high = flex.max(d_star_sq_hkl)

    ## A feeble attempt to determine an appropriate kernel width
    ## that seems to work reasonable in practice
    self.kernel_width=kernel_width
    if auto_kernel is not False:
      ## get the d_star_sq_array and sort it
      sort_permut = flex.sort_permutation(d_star_sq_hkl)
      ##
      if auto_kernel==True:
        number=number_of_sorted_reflections_for_auto_kernel
      else:
        number=int(auto_kernel)

      if number > d_star_sq_hkl.size():
        number = d_star_sq_hkl.size()-1

      self.kernel_width = d_star_sq_hkl[sort_permut[number]]-d_star_sq_low
      assert self.kernel_width > 0

    ## Making the d_star_sq_array
    assert (n_bins>1) ## assure that there are more then 1 bins for interpolation

    self.d_star_sq_array = chebyshev_lsq_fit.chebyshev_nodes(
          n=n_bins,
          low=d_star_sq_low,
          high=d_star_sq_high,
          include_limits=True)

    ## Now get the average intensity please
    ##
    ## This step can be reasonably time consuming

    self.mean_I_array = scaling.kernel_normalisation(
      d_star_sq_hkl = d_star_sq_hkl,
      I_hkl = I_obs,
      epsilon = epsilons,
      d_star_sq_array = self.d_star_sq_array,
      kernel_width = self.kernel_width
      )

    self.var_I_array = scaling.kernel_normalisation(
      d_star_sq_hkl = d_star_sq_hkl,
      I_hkl = I_obs*I_obs,
      epsilon = epsilons*epsilons,
      d_star_sq_array = self.d_star_sq_array,
      kernel_width = self.kernel_width
      )

    self.var_I_array = self.var_I_array - self.mean_I_array*self.mean_I_array

    self.weight_sum = self.var_I_array = scaling.kernel_normalisation(
      d_star_sq_hkl = d_star_sq_hkl,
      I_hkl = I_obs*0.0+1.0,
      epsilon = epsilons*0.0+1.0,
      d_star_sq_array = self.d_star_sq_array,
      kernel_width = self.kernel_width
      )


    #assert flex.min( self.mean_I_array ) > 0
    sel_pos = self.mean_I_array > 0
    self.mean_I_array = self.mean_I_array.select(sel_pos)
    self.d_star_sq_array = self.d_star_sq_array.select(sel_pos)
    self.var_I_array = flex.log( self.var_I_array.select( sel_pos ) )


    self.mean_I_array = flex.log( self.mean_I_array )
    ## Fit a chebyshev polynome please
    normalizer_fit_lsq = chebyshev_lsq_fit.chebyshev_lsq_fit(
      n_term,
      self.d_star_sq_array,
      self.mean_I_array )
    self.normalizer = chebyshev_polynome(
      n_term,
      d_star_sq_low,
      d_star_sq_high,
      normalizer_fit_lsq.coefs)

    var_lsq_fit = chebyshev_lsq_fit.chebyshev_lsq_fit(
      n_term,
      self.d_star_sq_array,
      self.var_I_array )
    self.var_norm = chebyshev_polynome(
      n_term,
      d_star_sq_low,
      d_star_sq_high,
      var_lsq_fit.coefs)

    ws_fit = chebyshev_lsq_fit.chebyshev_lsq_fit(
      n_term,
      self.d_star_sq_array,
      self.weight_sum )
    self.weight_sum = chebyshev_polynome(
      n_term,
      d_star_sq_low,
      d_star_sq_high,
      ws_fit.coefs)




    ## The data wil now be normalised using the
    ## chebyshev polynome we have just obtained
    self.mean_I_array = flex.exp( self.mean_I_array)
    self.normalizer_for_miller_array =  flex.exp( self.normalizer.f(d_star_sq_hkl) )

    self.var_I_array = flex.exp( self.var_I_array )
    self.var_norm = flex.exp( self.var_norm.f(d_star_sq_hkl) )
    self.weight_sum = flex.exp( self.weight_sum.f(d_star_sq_hkl))


    self.normalised_miller = None
    self.normalised_miller_dev_eps = None

    if work_array.sigmas() is not None:
      self.normalised_miller = work_array.customized_copy(
        data = work_array.data()/self.normalizer_for_miller_array,
        sigmas = work_array.sigmas()/self.normalizer_for_miller_array
        ).set_observation_type(work_array)
      self.normalised_miller_dev_eps = self.normalised_miller.customized_copy(
        data = self.normalised_miller.data()/epsilons,
        sigmas = self.normalised_miller.sigmas()/epsilons)\
        .set_observation_type(work_array)

    else:
      self.normalised_miller = work_array.customized_copy(
        data = work_array.data()/self.normalizer_for_miller_array
        ).set_observation_type(work_array)
      self.normalised_miller_dev_eps = self.normalised_miller.customized_copy(
        data = self.normalised_miller.data()/epsilons)\
        .set_observation_type(work_array)

    ## all done
