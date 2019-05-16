from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx.array_family import flex
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import pair_analyses
import sys
from six.moves import range


class combined_scaling(object):
  def __init__(self,
               miller_array_x1,
               miller_array_x2,
               options=None,
               out=None):
    ## miller_array_x1 : 'reference'
    ## miller_array_x2 : 'derivative' or something similar
    ## options : a phil file with scaling options
    ## convergence : rescale after rejection?
    ##
    ##
    ## These are the tasks ahead
    ## 0) map to asu, common sets
    ## 1) least squares scaling of x1 and x2
    ## 2) local scaling of x1 and x2
    ## 3) outlier rejections
    ## 4) least squares scaling of x1 and x2
    ## 5) local scaling of x1 and x2
    ## 6) generation of delta F's
    ##
    ## we will perform operations on this arrays
    self.x1 = miller_array_x1.deep_copy().map_to_asu()
    self.x2 = miller_array_x2.deep_copy().map_to_asu()
    ## Get the scaling and rejection options
    self.lsq_options = None
    self.loc_options = None
    self.out_options = None
    self.overall_protocol = None
    self.cut_level_rms = None
    if options is not None:
      self.lsq_options = options.least_squares_options
      self.loc_options = options.local_scaling_options
      self.out_options = options.outlier_rejection_options
      self.overall_protocol = options.target
      self.cut_level_rms = None
    ## cycling options
    self.auto_cycle=False
    if options is not None:
      if options.iterations=='auto':
        self.auto_cycle=True

    self.max_iterations = 10
    if options is not None:
      self.max_iterations = options.max_iterations
    if not self.auto_cycle:
      assert self.max_iterations > 0

    ## output options
    self.out=out
    if self.out==None:
      self.out=sys.stdout
    ## get common sets, and make a reference copy please
    self.x1, self.x2 = self.x1.common_sets( self.x2 )
    self.s1 = self.x1.deep_copy().map_to_asu()
    self.s2 = self.x2.deep_copy().map_to_asu()

    scaling_tasks={'lsq':False, 'local':False }

    if self.overall_protocol=='ls':
      scaling_tasks['lsq']=True
      scaling_tasks['local']=False

    if self.overall_protocol=='loc':
      scaling_tasks['lsq']=False
      scaling_tasks['local']=True

    if self.overall_protocol=='ls_and_loc':
      scaling_tasks['lsq']=True
      scaling_tasks['local']=True

    print(scaling_tasks)
    print(self.overall_protocol)
    #assert ( scaling_tasks['lsq'] or scaling_tasks['local'] )

    self.convergence=False
    counter = 0

    print(file=self.out)
    print("==========================================", file=self.out)
    print("=             Relative scaling           =", file=self.out)
    print("==========================================", file=self.out)
    print(file=self.out)


    while not self.convergence:
      print(file=self.out)
      print("--------------------------", file=self.out)
      print("    Scaling cycle %i   "%(counter), file=self.out)
      print("--------------------------", file=self.out)

      if counter == 0:
        self.cut_level_rms = 3
        if options is not None:
          self.cut_level_rms = self.out_options.cut_level_rms_primary
      if counter > 0:
        self.cut_level_rms = 3
        if options is not None:
          self.cut_level_rms = self.out_options.cut_level_rms_secondary

      ## take the common sets
      self.s1 = self.s1.common_set(  self.x1 )
      self.s2 = self.s2.common_set(  self.x2 )
      self.x1, self.x2 = self.s1.common_sets( self.s2 )

      if scaling_tasks['lsq']:
        self.perform_least_squares_scaling()
      if scaling_tasks['local']:
        self.perform_local_scaling()
      num_reject = self.perform_outlier_rejection()
      if num_reject==0:
        self.convergence=True
      if not self.auto_cycle:
        if counter==self.max_iterations:
          self.convergence=True
      counter+=1

    ## Now the datasets have been scaled, we can return bnotyh datasets
    ## They can be used by the usewr to computer either
    ## - isomorphous differences
    ## - be used in an F1-F2 FFT
    del self.s1
    del self.s2




  def perform_least_squares_scaling(self):
    print(file=self.out)
    print("Least squares scaling", file=self.out)
    print("---------------------", file=self.out)
    print(file=self.out)

    ##-----Get options-----
    use_exp_sigmas=self.lsq_options.use_experimental_sigmas

    use_int = True
    if self.lsq_options.scale_data=='amplitudes':
      use_int=False

    use_wt=True
    if self.lsq_options.scale_target == 'basic':
      use_wt=False
    ##-------------------------------------------
    ls_scaling = relative_scaling.ls_rel_scale_driver(
        self.x1,
        self.x2,
        use_intensities=use_int,
        scale_weight=use_wt,
        use_weights=use_exp_sigmas)
    ls_scaling.show(out=self.out)
    ##----- Update the miller arrays please-------
    self.x1 = ls_scaling.native.deep_copy()
    self.x2 = ls_scaling.derivative.deep_copy()

  def perform_local_scaling(self):
    print(file=self.out)
    print("Local scaling", file=self.out)
    print("-------------", file=self.out)
    print(file=self.out)

    ##-----Get options-----
    use_exp_sigmas=self.loc_options.use_experimental_sigmas
    use_int = True

    if self.loc_options.scale_data=='amplitudes':
      use_int=False

    local_scaling_target_dictionary ={
      'local_moment':False,
      'local_lsq':False,
      'local_nikonov':False}
    local_scaling_target_dictionary[self.loc_options.scale_target]=True
    print(local_scaling_target_dictionary)

    ##--------------------
    local_scaling = relative_scaling.local_scaling_driver(
      self.x1,
      self.x2,
      local_scaling_target_dictionary,
      use_intensities=use_int,
      use_weights=use_exp_sigmas,
      max_depth=self.loc_options.max_depth,
      target_neighbours=self.loc_options.target_neighbours,
      sphere=self.loc_options.neighbourhood_sphere,
      out=self.out)

    self.x1 = local_scaling.native.deep_copy()
    self.x2 = local_scaling.derivative.deep_copy()


  def perform_outlier_rejection(self):
    print(file=self.out)
    print("Outlier rejections", file=self.out)
    print("------------------", file=self.out)
    print(file=self.out)
    print(" sigma criterion : %4.1f "%(self.out_options.cut_level_sigma), file=self.out)
    print(" rms criterion   : %4.1f "%(self.cut_level_rms), file=self.out)
    print(" protocol        :", self.out_options.protocol, file=self.out)
    print(file=self.out)

    outlier_protocol ={'solve':False,'rms':False, 'rms_and_sigma':False }
    ## Please set the protocol to what is specified
    outlier_protocol[ self.out_options.protocol ]=True

    outlier_rejection = pair_analyses.outlier_rejection(
      self.x1,
      self.x2,
      cut_level_rms=self.cut_level_rms,
      cut_level_sigma=self.out_options.cut_level_sigma
    )

    number_of_outliers = self.x1.size() - outlier_rejection.nat.size()
    number_of_outliers += self.x2.size() - outlier_rejection.der.size()

    self.x1 = outlier_rejection.nat.deep_copy()
    self.x2 = outlier_rejection.der.deep_copy()

    self.x1, self.x2 = self.x1.common_sets( self.x2 )
    return( number_of_outliers )




class ano_scaling(object):
  def __init__(self,
               miller_array_x1,
               options=None,
               out=None):
    ## These are the tasks ahead
    ##
    ## 1) splitting up x1 in hemsispheres x1p x1n
    ## 2) sumbiut the two halve data sets to the combined scaler
    ##

    assert miller_array_x1.anomalous_flag()
    assert miller_array_x1.indices().size() > 0

    self.options = options

    self.s1p, self.s1n = miller_array_x1.hemispheres_acentrics()

    self.s1p = self.s1p.set_observation_type( miller_array_x1 )

    self.s1n = self.s1n.customized_copy( indices=-self.s1n.indices() )
    self.s1n = self.s1n.set_observation_type( miller_array_x1 )


    assert self.s1p.indices().size() == self.s1n.indices().size()

    self.x1p = self.s1p.deep_copy()
    self.x1n = self.s1n.deep_copy()

    ## Now we have a 'native' and a 'derivative'
    ##
    ## Submit these things to the combined scaler
    if self.options is not None:
      ano_scaler=combined_scaling(
        self.x1p,
        self.x1n,
        options,
        out)
      self.s1p = ano_scaler.x1.deep_copy()
      self.s1n = ano_scaler.x2.deep_copy()
      del ano_scaler



class naive_fa_estimation(object):
  def __init__(self,
               ano,
               iso,
               options,
               out=None):
    if out == None:
      out = sys.stdout

    ## get stuff
    self.options = options
    self.iso = iso.deep_copy().map_to_asu()
    self.ano = ano.deep_copy().map_to_asu()
    ## get common sets
    self.iso, self.ano = self.iso.common_sets( self.ano )

    ## perform normalisation
    normalizer_iso = absolute_scaling.kernel_normalisation(
      self.iso, auto_kernel=True, n_term=options.number_of_terms_in_normalisation_curve)
    normalizer_ano = absolute_scaling.kernel_normalisation(
      self.ano, auto_kernel=True, n_term=options.number_of_terms_in_normalisation_curve)

    self.fa = self.iso.customized_copy(
      data = flex.sqrt( self.iso.data()*self.iso.data()\
               /normalizer_iso.normalizer_for_miller_array
               +
               self.ano.data()*self.ano.data()\
               /normalizer_ano.normalizer_for_miller_array
              ),
      sigmas = flex.sqrt( self.iso.sigmas()*self.iso.sigmas()\
               /(normalizer_iso.normalizer_for_miller_array*
                 normalizer_iso.normalizer_for_miller_array
                 )
               +
               self.ano.sigmas()*self.ano.sigmas()\
               /(normalizer_ano.normalizer_for_miller_array
                 *normalizer_ano.normalizer_for_miller_array)
              ))





class cns_fa_driver(object):
  def __init__(self,
               lambdas):
    ## first generate all anomalous differences
    self.ano_and_iso = []
    self.na_ano=0

    for set in lambdas:
      assert set.is_xray_amplitude_array()
      if set.anomalous_flag():
        self.na_ano+=1
        plus, minus = set.hemispheres_acentrics()
        d_ano = plus.customized_copy(
          data = flex.abs( plus.data() - minus.data() ),
          sigmas = flex.sqrt( plus.sigmas()*plus.sigmas() +
                              minus.sigmas()*minus.sigmas() )
          ).set_observation_type( set )
        self.ano_and_iso.append( d_ano )
    #now generate all isomorphous differences
    self.n_iso=0
    for set1 in range(len(lambdas)):
      for set2 in range(set1+1,len(lambdas)):
        self.n_iso+=1
        t1 = lambdas[set1].average_bijvoet_mates().set_observation_type(
          lambdas[set1] )
        t2 = lambdas[set2].average_bijvoet_mates().set_observation_type(
          lambdas[set2] )
        tmp1,tmp2 = t1.common_sets(t2)
        tmp1 = tmp1.customized_copy(
          data = flex.abs( tmp1.data() - tmp2.data() ),
          sigmas = flex.sqrt( tmp1.sigmas()*tmp1.sigmas() +
                              tmp2.sigmas()*tmp2.sigmas() )
        ).set_observation_type( tmp1 )
        self.ano_and_iso.append( tmp1 )


    self.normalise_all()
    self.average_all()

  def normalise_all(self):
    ## normalise all difference data please
    for set in self.ano_and_iso:
      tmp_norm = absolute_scaling.kernel_normalisation(
        set,
        auto_kernel=True)
      set = tmp_norm.normalised_miller.deep_copy().set_observation_type(
        tmp_norm.normalised_miller)


  def average_all(self):
    ## get started quickly please
    mean_index = self.ano_and_iso[0].indices()
    mean_data = self.ano_and_iso[0].data()
    mean_sigmas = self.ano_and_iso[0].sigmas()

    ## loop over the remaining arrays
    for set_no in range( 1,len(self.ano_and_iso) ):
      mean_index = mean_index.concatenate(
        self.ano_and_iso[ set_no ].indices() )

      mean_data= mean_data.concatenate(
        self.ano_and_iso[ set_no ].data() )

      mean_sigmas = mean_sigmas.concatenate(
        self.ano_and_iso[ set_no ].sigmas() )

    final_miller = self.ano_and_iso[0].customized_copy(
      indices= mean_index,
      data = mean_data,
      sigmas = mean_sigmas).set_observation_type( self.ano_and_iso[0] )

    final_miller = final_miller.f_as_f_sq()

    merged = final_miller.merge_equivalents()
    merged.show_summary()
    self.fa = merged.array().set_observation_type(final_miller).f_sq_as_f()



class mum_dad(object):
  def __init__(self,
               lambda1,
               lambda2,
               k1=1.0):
    ## assumed is of course that the data are scaled.
    ## lambda1 is the 'reference'
    self.w1=lambda1.deep_copy()
    self.w2=lambda2.deep_copy()

    if not self.w1.is_xray_amplitude_array():
      self.w1 = self.w1.f_sq_as_f()
    if not self.w2.is_xray_amplitude_array():
      self.w2 = self.w2.f_sq_as_f()

    self.w1, self.w2 = self.w1.common_sets( self.w2 )

    l1p, l1n = self.w1.hemispheres_acentrics()
    self.mean1 = l1p.data()+l1n.data()
    self.diff1 = l1p.data()-l1n.data()
    self.v1 = ( l1p.sigmas()*l1p.sigmas() +
                l1n.sigmas()*l1n.sigmas() )

    l2p, l2n = self.w2.hemispheres_acentrics()
    self.mean2 = l2p.data()+l2n.data()
    self.diff2 = l2p.data()-l2n.data()
    self.v2 = ( l2p.sigmas()*l2p.sigmas() +
                l2n.sigmas()*l2n.sigmas() )

    self.new_diff = flex.abs( (self.diff1 + k1*self.diff2)/2.0 )
    self.new_sigma_mean = flex.sqrt( (self.v1+k1*k1*self.v2)/2.0 )

    self.dad = l1p.customized_copy(
      data = self.new_diff,
      sigmas = self.new_sigma_mean ).set_observation_type( self.w1 )


class singh_ramasheshan_fa_estimate(object):
  def __init__(self,
               w1,
               w2,
               k1,
               k2):
    self.w1=w1.deep_copy()
    self.w2=w2.deep_copy()

    if self.w1.is_xray_amplitude_array():
      self.w1 = self.w1.f_as_f_sq()
    if self.w2.is_xray_amplitude_array():
      self.w2 = self.w2.f_as_f_sq()

    ## common sets please
    self.w1,self.w2 = self.w1.common_sets( self.w2 )

    ## get differences and sums please
    self.p1, self.n1 = self.w1.hemispheres_acentrics()
    self.p2, self.n2 = self.w2.hemispheres_acentrics()

    self.diff1 = self.p1.data() - self.n1.data()
    self.diff2 = self.p2.data() - self.n2.data()

    self.s1 =   self.p1.sigmas()*self.p1.sigmas()\
              + self.n1.sigmas()*self.n1.sigmas()
    self.s1 =  flex.sqrt( self.s1 )

    self.s2 =   self.p2.sigmas()*self.p2.sigmas()\
              + self.n2.sigmas()*self.n2.sigmas()
    self.s2 =  flex.sqrt( self.s2 )

    self.sum1 = self.p1.data() + self.n1.data()
    self.sum2 = self.p2.data() + self.n2.data()

    self.k1_sq = k1*k1
    self.k2_sq = k2*k2

    self.determinant=None
    self.fa=None
    self.sigfa=None


    self.selector=None
    self.iselector=None

    self.a=None
    self.b=None
    self.c=None

    self.compute_fa_values()

    self.fa = self.p1.customized_copy(
      data = self.fa,
      sigmas = self.sigfa).set_observation_type( self.p1 )

  def set_sigma_ratio(self):
    tmp = self.s2/flex.abs( self.diff2 +1e-6 ) +\
          self.s1/flex.abs( self.diff1 +1e-6 )
    tmp = tmp/2.0
    self.sigfa = tmp

  def compute_coefs(self):
    self.a = ( self.k2_sq*self.k2_sq
               + self.k2_sq*(1 + self.k1_sq)
               + (self.k1_sq-1)*(self.k1_sq-1)
             )
    self.b = -self.k2_sq*(self.sum1 + self.sum2) \
             -(self.k1_sq-1)*(self.sum1 - self.sum2)
    self.c = 0.25*(self.sum1 - self.sum2)*(self.sum1 - self.sum2)\
            +(1.0/8.0)*self.k2_sq*(  self.diff2*self.diff2
                                   + self.diff1*self.diff1/self.k1_sq)

  def compute_determinant(self):
    self.determinant = self.b*self.b - 4.0*self.a*self.c
    self.selector = (self.determinant>0)
    self.iselector = self.selector.iselection()

  def compute_fa_values(self):
    self.compute_coefs()
    self.compute_determinant()
    reset_selector = (~self.selector).iselection()
    self.determinant = self.determinant.set_selected(  reset_selector, 0 )

    choice1 = -self.b + flex.sqrt( self.determinant )
    choice1 /= 2*self.a
    choice2 = -self.b - flex.sqrt( self.determinant )
    choice2 /= 2*self.a

    select1 = choice1 > choice2
    select2 = ~select1

    choice1 = choice1.set_selected( select1.iselection(), 0 )
    choice2 = choice2.set_selected( select2.iselection(), 0 )

    choice1 = choice1+choice2
    select1 = (choice1<0).iselection()
    choice1 = choice1.set_selected( select1 , 0 )
    self.fa =  choice1 + choice2

    self.set_sigma_ratio()

    self.sigfa = self.sigfa*self.fa



class twmad_fa_driver(object):
  def __init__(self,
               lambda1,
               lambda2,
               k1,
               k2,
               options,
               out=None):
    self.out=out
    if self.out==None:
      self.out=sys.stdout

    self.options = options
    print("FA estimation", file=self.out)
    print("=============", file=self.out)

    if k1 is None:
      raise Sorry("f\"(w1)/f\"(w2) ratio is not defined. Please provide f\" values upon input")

    if k2 is None:
      if self.options.protocol=='algebraic':
        raise Sorry("""
delta f' f\" ratio is not defined.
Either provide f' and f\" values upon input,
or chose different Fa estimation protocol.
               """)

    self.options = options

    protocol = {'algebraic': False,
                'cns': False,
                'combine_ano': False}
    protocol[ self.options.protocol ] = True

    self.fa_values = None

    if protocol['algebraic']:
      print(" Using algebraic approach to estimate FA values ", file=self.out)
      print(file=self.out)
      tmp = singh_ramasheshan_fa_estimate(
        lambda1,
        lambda2,
        k1,
        k2)
      self.fa_values = tmp.fa.f_sq_as_f()

    if protocol['cns']:
      print(" Using CNS approach to estimate FA values ", file=self.out)
      print(file=self.out)

      tmp = cns_fa_driver( [lambda1, lambda2] )
      self.fa_values = tmp.fa

    if protocol['combine_ano']:
      print(" Combining anomalous data only", file=self.out)
      print(file=self.out)

      tmp = mum_dad(
        lambda1,
        lambda2,
        k1)
      self.fa_values = tmp.dad

    norma = absolute_scaling.kernel_normalisation(
      self.fa_values,
      auto_kernel=True)

    self.fa_values = norma.normalised_miller.f_sq_as_f()
