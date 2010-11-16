from cctbx import miller
from cctbx import crystal
from libtbx import table_utils
import cctbx.sgtbx.cosets
from cctbx.crystal import reindex
from cctbx.array_family import flex
from libtbx.utils import Sorry
import mmtbx.scaling
from scitbx.python_utils import robust_statistics
import sys, math
import scitbx.lbfgs



class reindexing(object):
  __doc__=""" Reindexing matrices """
  def __init__(self,
               set_a,
               set_b,
               out=None,
               relative_length_tolerance=0.05,
               absolute_angle_tolerance=3.0,
               lattice_symmetry_max_delta=3.0,
               file_name=None):

    self.relative_length_tolerance=relative_length_tolerance
    self.absolute_angle_tolerance=absolute_angle_tolerance
    self.lattice_symmetry_max_delta=lattice_symmetry_max_delta
    self.out = out
    if self.out is None:
      self.out = sys.stdout
    self.file_name=file_name

    ## make deep copy for safety, and make it non-anomalous
    self.set_b_ori = set_b.deep_copy().set_observation_type( set_b )
    self.set_a = set_a.deep_copy().set_observation_type( set_a ).average_bijvoet_mates()
    self.set_b = set_b.deep_copy().set_observation_type( set_b ).average_bijvoet_mates()

    self.set_a_xs = crystal.symmetry( unit_cell=self.set_a.unit_cell(),
                                      space_group=self.set_a.space_group() )
    self.set_b_xs = crystal.symmetry( unit_cell=self.set_b.unit_cell(),
                                      space_group=self.set_b.space_group() )

    self.cosets = reindex.reindexing_operators(
      xs1=self.set_a_xs,
      xs2=self.set_b_xs,
      relative_length_tolerance=self.relative_length_tolerance,
      absolute_angle_tolerance=self.absolute_angle_tolerance,
      max_delta=self.lattice_symmetry_max_delta,
      anomalous_flag=self.set_a.anomalous_flag(),
      out=self.out
      )

    self.set_a_to_niggli = self.cosets.cb_op_to_n_1
    self.set_b_to_niggli = self.cosets.cb_op_to_n_2

    self.set_a_nig = self.set_a.change_basis( self.set_a_to_niggli  ).map_to_asu()
    self.set_b_nig = self.set_b.change_basis( self.set_b_to_niggli  ).map_to_asu()

    self.ops_in_niggli_setting = self.cosets.cb_ops_in_niggli_setting()
    self.nice_cb_ops = self.cosets.combined_cb_ops()

    self.cc_values = []
    self.matches = []
    if len(self.nice_cb_ops)>0:
      self.analyse()
      self.select_and_transform()
    else:
      print >> self.out, "No simple reindexing relation found between unit cells"

  def analyse(self):
    table_data=[]
    for (cb_op, comb_cb_op) in zip(self.ops_in_niggli_setting,
                                   self.nice_cb_ops):
      tmp_set_b_nig = self.set_b_nig.change_basis( cb_op ).map_to_asu()
      tmp_set_a_nig, tmp_set_b_nig = self.set_a_nig.common_sets(
        tmp_set_b_nig, assert_is_similar_symmetry=False )
      tmp_cc = tmp_set_a_nig.correlation(tmp_set_b_nig,
                                         assert_is_similar_symmetry=False )
      self.cc_values.append(  tmp_cc.coefficient()  )
      self.matches.append( 100.0*float(   tmp_set_a_nig.indices().size()  )/
                           float(  self.set_a_nig.indices().size()  )   )



  def select_and_transform(self,
                           matches_cut_off=0.75
                           ):
    ## hopsa
    max_cc=-1.0
    location = 0
    table_data=[]
    for ii in range(len(self.nice_cb_ops)):
      table_data.append(
        [self.nice_cb_ops[ii].as_hkl(),
         "%4.3f"%(self.cc_values[ii]),
         "%4.1f"%(self.matches[ii]),
         '   ']
        )

      if self.matches[ii]>=matches_cut_off:
        if max_cc<self.cc_values[ii]:
          max_cc = self.cc_values[ii]
          location = ii

    legend = ('Operator', 'Correlation', 'matches (%)', 'choice')
    table_data[location][3]=' <--- '
    self.table = table_utils.format([legend]+table_data,
                                       comments=None,
                                       has_header=True,
                                       separate_rows=False,
                                       prefix='| ',
                                       postfix=' |')
    print >> self.out
    print >> self.out, "Reference analyses"
    print >> self.out, "   The following reindexing operators have been found:"
    print >> self.out
    print >> self.out, self.table
    print >> self.out
    if str(self.nice_cb_ops[location].as_hkl()) == "h,k,l":
      print >> self.out, "The data doesn't need to be reindexed! Indexing is consistent between these datasets."
    else:
      print >> self.out, "If the data is reindexed with operator (%s), the correlation of"%( self.nice_cb_ops[location].as_hkl() )
      print >> self.out, "the intensities to the reference data is maximized. "
      if self.file_name is not None:
        print self.file_name
        print >> self.out, "This can be done for instance with:"
        print >> self.out, "  phenix.reflection_file_converter %s --change_of_basis=\"%s\" <output_options> "%(self.file_name, self.nice_cb_ops[location].as_hkl() )
    print >> self.out, "-------------------------------------------------------------------------------"
    ##  change things in primitive setting

    transform_b = self.set_b_ori.change_basis( self.set_b_to_niggli ).change_basis(
      self.ops_in_niggli_setting[location]  ).change_basis(
      self.set_b_to_niggli.inverse() ).map_to_asu().set_observation_type(
      self.set_b_ori)
    return ( transform_b )



class delta_generator(object):
  def __init__(self,
               nat,
               der,
               nsr_bias=1.0):
    self.nat=nat.deep_copy()
    self.der=der.deep_copy()
    self.nsr_bias=1.0/nsr_bias

    assert self.nat.is_real_array()
    assert self.nat.is_real_array()

    if self.nat.is_xray_intensity_array():
      self.nat.f_sq_as_f()
    if self.der.is_xray_intensity_array():
      self.der.f_sq_as_f()

    self.nat,self.der = self.nat.common_sets(self.der)

    self.der = self.der.customized_copy(
      data = self.der.data()*self.nsr_bias,
      sigmas = self.der.sigmas()*self.nsr_bias).set_observation_type(
        self.der)

    self.delta_f=self.nat.customized_copy(
      data = ( self.der.data() - self.nat.data() ),
      sigmas = flex.sqrt( self.der.sigmas()*self.der.sigmas()+
                          self.nat.sigmas()*self.nat.sigmas() )
      ).set_observation_type( self.nat )

    self.abs_delta_f=self.nat.customized_copy(
      data = flex.abs( self.der.data() - self.nat.data() ),
      sigmas = flex.sqrt( self.der.sigmas()*self.der.sigmas()+
                          self.nat.sigmas()*self.nat.sigmas() )
      ).set_observation_type( self.der )

    if not self.nat.is_xray_intensity_array():
      self.nat.f_as_f_sq()
    if not self.der.is_xray_intensity_array():
      self.der.f_as_f_sq()

    self.delta_i=self.nat.customized_copy(
      data = ( self.der.data() - self.nat.data() ),
      sigmas = flex.sqrt( self.der.sigmas()*self.der.sigmas()+
                          self.nat.sigmas()*self.nat.sigmas() )
      ).set_observation_type( self.nat )

    self.abs_delta_i=self.nat.customized_copy(
      data = flex.abs( self.der.data() - self.nat.data() ),
      sigmas = flex.sqrt( self.der.sigmas()*self.der.sigmas()+
                          self.nat.sigmas()*self.nat.sigmas() )
      ).set_observation_type( self.der )



class outlier_rejection(object):
  def __init__(self,
               nat,
               der,
               cut_level_rms=3,
               cut_level_sigma=0,
               method={'solve':False,'rms':False, 'rms_and_sigma':True},
               out=None
               ):
    self.out=out
    self.method = method
    if self.out == None:
      self.out = sys.stdout

    self.cut_level_rms=cut_level_rms
    self.cut_level_sigma=cut_level_sigma
    ## Just make sure that we have the common sets
    self.nat = nat.deep_copy()
    self.der = der.deep_copy()

    self.nat, self.der = self.nat.common_sets(self.der)
    #Make sure we have amplitudes
    assert self.nat.is_real_array()
    assert self.nat.is_real_array()

    if self.nat.is_xray_intensity_array():
      self.nat=self.nat.f_sq_as_f()
    if self.der.is_xray_intensity_array():
      self.der=self.der.f_sq_as_f()

    ## Construct delta f's
    delta_gen = delta_generator( self.nat, self.der )
    self.delta_f = delta_gen.abs_delta_f
    ## Set up a binner please
    self.delta_f.setup_binner_d_star_sq_step(auto_binning=True)
    ## for each bin, I would like to compute
    ## mean dF**2
    ## mean sigma**2
    self.mean_df2 = self.delta_f.mean_of_intensity_divided_by_epsilon(
      use_binning=True,
      return_fail=1e12)
    self.mean_sdf2 = self.delta_f.mean_of_squared_sigma_divided_by_epsilon(
      use_binning=True,
      return_fail=1e12)

    self.result = flex.bool(  self.delta_f.indices().size(), True )

    self.detect_outliers()
    self.remove_outliers()

  def detect_outliers(self):
    count_true = 0
    if (self.method['solve']):
      count_true+=1
    if (self.method['rms']):
      count_true+=1
    if (self.method['rms_and_sigma']):
      count_true+=1


    if not count_true==1:
      raise Sorry("Outlier removal protocol not specified properly")

    if (self.method['solve']):
      self.detect_outliers_solve()
    if (self.method['rms']):
      self.detect_outliers_rms()
    if (self.method['rms_and_sigma']):
      self.detect_outliers_sigma()




  def detect_outliers_solve(self):
    """
    TT says:
    I toss everything > 3 sigma in the scaling,
    where sigma comes from the rms of everything being scaled:

    sigma**2 = <delta**2>- <experimental-sigmas**2>

    Then if a particular
    delta**2 > 3 sigma**2 + experimental-sigmas**2
    then I toss it.
    """
    terwilliger_sigma_array = flex.double(self.mean_df2.data) -\
                              flex.double(self.mean_sdf2.data)

    for bin_number in self.delta_f.binner().range_all():
      ## The selection tells us wether or not somthing is in the correct bin
      selection =  self.delta_f.binner().selection( bin_number ).iselection()
      ## Now just make a global check to test for outlierness:
      tmp_sigma_array =  terwilliger_sigma_array[bin_number] -\
                         self.delta_f.sigmas()*self.delta_f.sigmas()
      tmp_sigma_array = flex.sqrt(tmp_sigma_array)*self.cut_level_rms

      potential_outliers = ( self.delta_f.data()  >  tmp_sigma_array )
      potential_outliers =  potential_outliers.select( selection )

      self.result = self.result.set_selected( selection, potential_outliers )

    print >> self.out
    print >> self.out, " %8i potential outliers detected" %(
      self.result.count(True) )
    print >> self.out, " They will be removed from the data set"
    print >> self.out


  def detect_outliers_rms(self):
    for bin_number in self.delta_f.binner().range_all():
      selection =  self.delta_f.binner().selection( bin_number ).iselection()
      potential_outliers = (
        self.delta_f.data()  >  self.cut_level_rms*math.sqrt(
        self.mean_df2.data[bin_number])  )
      potential_outliers =  potential_outliers.select( selection )
      self.result = self.result.set_selected( selection, potential_outliers )

    print >> self.out
    print >> self.out, " %8i potential outliers detected" %(
      self.result.count(True) )
    print >> self.out, " They will be removed from the data set"
    print >> self.out


  def detect_outliers_sigma(self):
    ## Locate outliers in native
    potential_outlier_nat = (self.nat.data()/self.nat.sigmas()
                               < self.cut_level_sigma)
    nat_select = potential_outlier_nat.iselection()

    ## Locate outliers in derivative
    potential_outlier_der = (self.der.data()/self.der.sigmas()
                               <self.cut_level_sigma)
    der_select = potential_outlier_der.iselection()

    for bin_number in self.delta_f.binner().range_all():
      ## RMS outlier removal
      selection =  self.delta_f.binner().selection( bin_number ).iselection()
      potential_outliers = (
        self.delta_f.data()  >  self.cut_level_rms*math.sqrt(
        self.mean_df2.data[bin_number])  )
      potential_outliers =  potential_outliers.select( selection )
      self.result = self.result.set_selected( selection, potential_outliers )

    self.result = self.result.set_selected( nat_select, True )
    self.result = self.result.set_selected( der_select, True )

    print >> self.out
    print >> self.out, " %8i potential outliers detected" %(
      self.result.count(True) )
    print >> self.out, " They will be removed from the data set"
    print >> self.out


  def remove_outliers(self):
    potential_outliers = self.nat.select( self.result )

    matches = miller.match_indices( self.nat.indices(),
                                    potential_outliers.indices()  )

    self.nat = self.nat.select( matches.single_selection(0) )

    self.nat, self.der = self.nat.common_sets(self.der)


class f_double_prime_ratio(object):
  def __init__(self,
               lambda1,
               lambda2):
    ## make sure we have anomalous data
    assert ( lambda1.anomalous_flag() )
    assert ( lambda2.anomalous_flag() )
    ## make temporary copies
    tmp_l1 = lambda1.deep_copy()
    tmp_l2 = lambda2.deep_copy()
    ##make sure we have intensities
    assert tmp_l1.is_real_array()
    if tmp_l1.is_xray_amplitude_array():
      tmp_l1=tmp_l1.f_as_f_sq()

    assert tmp_l2.is_real_array()
    if tmp_l2.is_xray_amplitude_array():
      tmp_l2=tmp_l2.f_as_f_sq()

    ## we only need the anomalous diffs
    l1p, l1n = tmp_l1.hemispheres_acentrics()
    self.diff1 = l1p.data()-l1n.data()
    self.v1 = ( l1p.sigmas()*l1p.sigmas() +
                l1n.sigmas()*l1n.sigmas() )

    l2p, l2n = tmp_l2.hemispheres_acentrics()
    self.diff2 = l2p.data()-l2n.data()
    self.v2 = ( l2p.sigmas()*l2p.sigmas() +
                l2n.sigmas()*l2n.sigmas() )

    self.x=flex.double([1.0])
    scitbx.lbfgs.run(target_evaluator=self)
    self.ratio=self.x[0]

  def compute_functional_and_gradients(self):
    f = self.compute_functional()
    g = self.compute_gradient()
    ##g2 =  self.compute_gradient_fd()
    return f,g

  def compute_functional(self):
    top = self.diff1 - self.diff2*self.x[0]
    top= top*top
    bottom = self.v1+self.v2*self.x[0]*self.x[0]
    result = top/bottom
    result=flex.sum(result)
    return result

  def compute_gradient(self):
    tmp_bottom = self.v1+self.v2*self.x[0]*self.x[0]
    tmp_top = self.diff1 - self.diff2*self.x[0]
    part1 = -2.0*self.x[0]*tmp_top*tmp_top*self.v2/( tmp_bottom*tmp_bottom )
    part2 = -2.0*self.diff2*tmp_top/tmp_bottom
    result=flex.sum( part1+part2)
    return(flex.double([result]))

  def compute_gradient_fd(self):
    h=0.0000001
    current = self.compute_functional()
    self.x[0]+=h
    new = self.compute_functional()
    self.x[0]-=h
    result = current - new
    result/=-h
    return flex.double([result])

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "Estimated ratio of fdp(w1)/fwp(w2): %3.2f"%(self.x[0])
    print >> out



class delta_f_prime_f_double_prime_ratio(object):
  def __init__(self,
               lambda1,
               lambda2,
               level=1.0):
    self.tmp_l1 = lambda1.deep_copy()
    self.tmp_l2 = lambda2.deep_copy()

    ## make sure we have amplitudes please
    if self.tmp_l1.is_xray_intensity_array():
      self.tmp_l1 = self.tmp_l1.f_sq_as_f()
    if self.tmp_l2.is_xray_intensity_array():
      self.tmp_l2 = self.tmp_l2.f_sq_as_f()

    ## Now this is done, swe have to set up a binner
    self.tmp_l1, self.tmp_l2 = self.tmp_l1.common_sets( self.tmp_l2 )
    self.tmp_l1.setup_binner_d_star_sq_step( auto_binning=True )
    self.tmp_l2.use_binner_of( self.tmp_l1 )

    self.tmp_l1_no_ano = self.tmp_l1.average_bijvoet_mates()
    self.tmp_l2_no_ano = self.tmp_l2.average_bijvoet_mates()
    self.tmp_l1_no_ano.setup_binner_d_star_sq_step( auto_binning=True )
    self.tmp_l2_no_ano.setup_binner_d_star_sq_step( auto_binning=True )

    self.plus2, self.minus2 = self.tmp_l2.hemispheres_acentrics()
    self.plus2.use_binning_of( self.tmp_l1_no_ano )
    self.minus2.use_binning_of( self.tmp_l1_no_ano )

    ## we assume that the data is properly scaled of course
    ## Loop over all bins, and in each bin,
    ## compute <diso>, <df> and their ratio
    estimates = flex.double()
    count = 0
    for bin in self.tmp_l1.binner().range_all():
      selection =  self.tmp_l1_no_ano.binner().selection( bin ).iselection()
      tmp1 = self.tmp_l1_no_ano.select( selection ).data()
      tmp2 = self.tmp_l2_no_ano.select( selection ).data()
      stmp1 = self.tmp_l1_no_ano.select( selection ).sigmas()
      stmp2 = self.tmp_l2_no_ano.select( selection ).sigmas()

      selection = self.plus2.binner().selection( bin ).iselection()
      tmp3 = self.plus2.select( selection ).data()
      tmp4 = self.minus2.select( selection ).data()
      stmp3 = self.plus2.select( selection ).sigmas()
      stmp4 = self.minus2.select( selection ).sigmas()

      tmpiso=None
      tmpsiso=None
      tmpano=None
      tmpsano=None

      if tmp1.size() > 0:
        tmpiso = flex.mean( (tmp1 - tmp2)*(tmp1 - tmp2) )
        tmpsiso = flex.mean( stmp1*stmp1 +  stmp2*stmp2  )
      if tmp3.size() > 0:
        tmpano = flex.mean( (tmp3 - tmp4)*(tmp3 - tmp4) )
        tmpsano = flex.mean( stmp3*stmp3 +  stmp4*stmp4  )

      if tmp1.size() > 0: ## make sure something is there
        if tmp3.size() > 0: ## make sure something is there
          if math.sqrt(tmpiso/tmpsiso)>=level:  ## s/n is okai
            if math.sqrt(tmpano/tmpsano)>=level: ## s/n is okai
              delta_iso_sq = math.sqrt( tmpiso  )
              delta_ano_sq = math.sqrt( tmpano  )
              tmp=2.0*delta_iso_sq/delta_ano_sq
              print bin, tmp
              estimates.append( tmp )
              count+=1.0

    ## compute the trimean please
    self.ratio = None
    self.number = count
    if (self.number>0):
      self.ratio=robust_statistics.trimean( estimates )



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
