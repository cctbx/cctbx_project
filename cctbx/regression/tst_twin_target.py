from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import miller
from cctbx import xray
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.development import random_structure as rs
from scitbx import differential_evolution as de
import scitbx.lbfgs
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random
import math
from six.moves import range
from six.moves import zip

random.seed(0)
flex.set_random_seed(0)

def tst_ls_on_i():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P4' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = tmp.structure_factors( False, 3.5  ).f_calc()
  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data(),
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=0.0,
                        u_sol=0.0,
                        f_part=sfs.data(),
                        k_part=0.0,
                        u_part=0.0 )

  target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
    hkl_obs=sfs.indices(),
    i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0.0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  target = target_evaluator.target( sfs.data() )
  # The target vaslue should be zero
  assert approx_equal( target, 0, eps=1e-6)

  # the derivatives as well
  derivs_complex = target_evaluator.d_target_d_fmodel( sfs.data()  )
  derivs_ab = target_evaluator.d_target_d_ab( sfs.data()  )
  for cmplx,da,db in zip(derivs_complex,
                         derivs_ab[0],
                         derivs_ab[1]):
    assert approx_equal( cmplx.real, da, eps=1e-5)
    assert approx_equal( cmplx.imag, db, eps=1e-5)

  for alpha in flex.double(range(50))/100.0:
    #----------------------------------------------------------------
    # use fin diffs to check the derivatives to a and b
    old_target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
      hkl_obs=sfs.indices(),
      i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.1,
      w_obs=None,
      hkl_calc=sfs.indices(),
      space_group=sfs.space_group(),
      anomalous_flag=False,
      alpha=alpha,
      twin_law=[-1,0,0, 0,1,0, 0,0,-1]
      )
    old_target_value = old_target_evaluator.target( sfs.data() )
    old_derivs = old_target_evaluator.d_target_d_ab( sfs.data() )

    new_data =  sfs.data()
    h=0.0001

    for N_test in range(sfs.data().size() ):
      ori = complex( sfs.data()[N_test] )
      #print "----------------"
      #print alpha
      #print sfs.indices()[N_test]
      #print sfs.data()[N_test]
      new_data[N_test] = ori+complex(h,0)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_real = float((new_target_value-old_target_value)/h)
      new_data[N_test] = ori+complex(0,h)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_imag = float( (new_target_value-old_target_value)/h )
      # only use 'large' first derivatives.
      if old_derivs[0][N_test]>0:
        #print "real", N_test, fdif_real,old_derivs[0][N_test], (fdif_real-old_derivs[0][N_test])/old_derivs[0][N_test]
        if old_derivs[0][N_test]>2500:
          assert approx_equal( (fdif_real-old_derivs[0][N_test])/fdif_real,0, eps=1e-2)
      if old_derivs[1][N_test]>0:
        #print  "Imag", N_test, fdif_imag,old_derivs[1][N_test], (fdif_imag-old_derivs[1][N_test])/old_derivs[1][N_test]
        if old_derivs[1][N_test]>2500:
          assert approx_equal( (fdif_imag-old_derivs[1][N_test])/fdif_imag,0, eps=1e-2)
      new_data[N_test] = ori

  #-------------------------------------
  # use fin diffs to test derivatives wrst alpha, the twin fraction
  h=0.0000000001

  target_evaluator = xray.least_squares_hemihedral_twinning_on_i(
    hkl_obs=sfs.indices(),
    i_obs=( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  tst_alpha = [0.1, 0.2, 0.3, 0.4, 0.5]

  for ii in tst_alpha:
    target_evaluator.alpha(ii)
    old_target = target_evaluator.target( sfs.data()*1.0 )
    target_evaluator.alpha( ii + h )
    new_target = target_evaluator.target( sfs.data()*1.0 )
    fd = (new_target-old_target)/h
    target_evaluator.alpha(ii)
    an = target_evaluator.d_target_d_alpha(sfs.data()*1.0)
    assert approx_equal( fd/an , 1.0, eps=1e-2 )






def tst_ls_on_f():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P4' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = tmp.structure_factors( False, 3.5  ).f_calc()
  f_mod = xray.f_model_core_data( hkl = sfs.indices(),
                        f_atoms= sfs.data(),
                        f_mask = sfs.data(),
                        unit_cell = sfs.unit_cell(),
                        k_overall=1.0,
                        u_star=(0,0,0,0,0,0),
                        k_sol=0.0,
                        u_sol=0.0,
                        f_part=sfs.data(),
                        k_part=0.0,
                        u_part=0.0 )

  target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
    hkl_obs=sfs.indices(),
    f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0.0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  target = target_evaluator.target( sfs.data() )
  # The target vaslue should be zero
  assert approx_equal( target, 0, eps=1e-6)

  # the derivatives as well
  derivs_complex = target_evaluator.d_target_d_fmodel( sfs.data()  )
  derivs_ab = target_evaluator.d_target_d_ab( sfs.data()  )
  for cmplx,da,db in zip(derivs_complex,
                         derivs_ab[0],
                         derivs_ab[1]):
    assert approx_equal( cmplx.real, da, eps=1e-5)
    assert approx_equal( cmplx.imag, -db, eps=1e-5)

  for alpha in flex.double(range(50))/100.0:
    #----------------------------------------------------------------
    # use fin diffs to check the derivatives to a and b
    old_target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs=sfs.indices(),
      f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.1,
      w_obs=None,
      hkl_calc=sfs.indices(),
      space_group=sfs.space_group(),
      anomalous_flag=False,
      alpha=alpha,
      twin_law=[-1,0,0, 0,1,0, 0,0,-1]
      )
    old_target_value = old_target_evaluator.target( sfs.data() )
    old_derivs = old_target_evaluator.d_target_d_ab( sfs.data() )

    new_data =  sfs.data()
    h=0.0001
    checked=0
    for N_test in range(sfs.data().size() ):
      ori = complex( sfs.data()[N_test] )
      #print "----------------"
      #print alpha
      #print sfs.indices()[N_test]
      #print sfs.data()[N_test]
      new_data[N_test] = ori+complex(h,0)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_real = float((new_target_value-old_target_value)/h)
      new_data[N_test] = ori+complex(0,h)
      new_target_value = old_target_evaluator.target( new_data )
      fdif_imag = float( (new_target_value-old_target_value)/h )
      # only use 'large' first derivative
      if 1: #old_derivs[0][N_test]>0:
        #print "real", N_test, fdif_real,old_derivs[0][N_test], (fdif_real-old_derivs[0][N_test])/old_derivs[0][N_test]
        if old_derivs[0][N_test]>1:
          checked+=1
          assert approx_equal( (fdif_real-old_derivs[0][N_test])/fdif_real,0, eps=1e-3)
      if abs(old_derivs[1][N_test])>0:
        #print  "Imag", N_test, fdif_imag,old_derivs[1][N_test], (fdif_imag-old_derivs[1][N_test])/old_derivs[1][N_test]
        if old_derivs[1][N_test]>1:
          checked+=1
          assert approx_equal( (fdif_imag-old_derivs[1][N_test])/fdif_imag,0, eps=1e-3)
      new_data[N_test] = ori
    assert checked>0
  #-------------------------------------
  # use fin diffs to test derivatives wrst alpha, the twin fraction
  h=0.00001

  target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
    hkl_obs=sfs.indices(),
    f_obs=flex.sqrt( flex.abs(sfs.data())*flex.abs(sfs.data()))*1.0,
    w_obs=None,
    hkl_calc=sfs.indices(),
    space_group=sfs.space_group(),
    anomalous_flag=False,
    alpha=0,
    twin_law=[-1,0,0, 0,1,0, 0,0,-1]
    )

  tst_alpha = [0.1, 0.2, 0.3, 0.4, 0.5]

  for ii in tst_alpha:
    target_evaluator.alpha(ii)
    old_target = target_evaluator.target( sfs.data()*1.0 )
    target_evaluator.alpha( ii + h )
    new_target = target_evaluator.target( sfs.data()*1.0 )
    fd = (new_target-old_target)/h
    target_evaluator.alpha(ii)
    an = target_evaluator.d_target_d_alpha(sfs.data()*1.0)
    assert approx_equal( fd/an , 1.0, eps=1e-2 )


class find_max(object):
  def __init__(self,
               tot):
    self.to=tot
    self.n=2
    self.x=[]
    self.domain= [ (0,5), (0,5) ]
    self.eps=1e-16
    self.optim =  de.differential_evolution_optimizer(
      self,
      population_size=100,
      f=0.3,
      cr=0.9,
      n_cross=1,
      show_progress=True,
      show_progress_nth_cycle=1,eps=1e-12)


  def print_status(self, mean, min, best,count):
    None
    #print count, " \t:::  ", mean, min, list(best), list(
    #  self.to.d_log_p_d_f( best[0],best[1] ) )


  def target(self, vector):
    f1=vector[0]
    if f1<self.eps:
      f1=self.eps

    f2=vector[1]
    if f2<self.eps:
      f2=self.eps

    return -self.to.log_p(f1,f2)


class minim(object):
  def __init__(self,
               tot, fm1, fm2):
    self.tot=tot
    self.x = flex.double([fm1,fm2])
    scitbx.lbfgs.run(target_evaluator=self)

  def g(self,f1,f2,h=0.00001):
    fo = self.tot.log_p( f1, f2 )
    fp1= self.tot.log_p( f1+h, f2 )
    fp2= self.tot.log_p( f1, f2+h )
    tmp = [fp1-fo, fp2-fo]
    tmp = flex.double( tmp )/h
    return tmp

  def compute_functional_and_gradients(self):
    f = self.tot.log_p( self.x[0], self.x[1] )
    g = self.tot.d_log_p_d_f( self.x[0], self.x[1] )
    g1= self.g(self.x[0], self.x[1] )
    # test gradients wrst fin diffs
    div = int(math.log(abs(g[0])))
    if div<0:
      div=0
    div = math.exp(div)
    if div<1e5: # if this divisor is very large, we have numerical issues .....
      assert approx_equal( g[0]/div, g1[0]/div, eps=5e-1)

    div = int(math.log(abs(g[1])))
    if div<0:
      div=0
    div = math.exp(div)
    if div<1e5:
      assert approx_equal( g[1]/div, g1[1]/div, eps=5e-1)
    #print list(g), list(g1)
    return -f,-flex.double(g)



def tst_snd_der(tot, f1, f2 , h=0.01):
  # see AMS55, 25.3.22 and 25.3.26
  #
  snd_der = tot.dd_log_p_dd_f(f1,f2)
  tmp0 = tot.log_p( f1,   f2   )
  tmp1 = tot.log_p( f1-h, f2   )
  tmp2 = tot.log_p( f1+h, f2   )

  ddf1 = (tmp1-2.0*tmp0 + tmp2)/(h*h)

  tmp0 = tot.log_p( f1,   f2   )
  tmp1 = tot.log_p( f1, f2-h   )
  tmp2 = tot.log_p( f1, f2+h   )

  ddf2 = (tmp1-2.0*tmp0+tmp2)/(h*h)

  div = int(math.log(abs( snd_der[0] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( ddf1/div, snd_der[0]/div, eps=5e-2 )

  div = int(math.log(abs( snd_der[1] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( ddf2/div, snd_der[1]/div, eps=5e-2 )

  tmp0 = tot.log_p( f1-h,   f2-h )
  tmp1 = tot.log_p( f1-h,   f2+h )
  tmp2 = tot.log_p( f1+h,   f2-h )
  tmp3 = tot.log_p( f1+h,   f2+h )

  tmp = (tmp3-tmp2-tmp1+tmp0)
  tmp = tmp/(h*h*4.0)

  div = int(math.log(abs( snd_der[2] )))
  if div<0:
    div=0
  div = math.exp(div)
  assert approx_equal( tmp/div, snd_der[2]/div, eps=5e-2 )


def tst_integration(twin_lik, range_value=7 ):
  #first find the maximum
  tmp_1_lbfgs = minim( twin_lik, 1, 1 )
  fm1=tmp_1_lbfgs.x[0]
  fm2=tmp_1_lbfgs.x[1]
  snd_der = twin_lik.dd_log_p_dd_f(fm1, fm2)
  is11 = snd_der[0]
  is22 = snd_der[1]
  is12 = snd_der[2]
  det = is11*is22-is12*is12
  s11= is22/det
  s22= is11/det
  s12= -is12/det
  s1= math.sqrt( abs(s11)*2.0 )
  s2= math.sqrt( abs(s22)*2.0 )
  quad_value = twin_lik.num_integrate(fm1, s1, fm2, s2, range_value)
  laplace_value= twin_lik.laplace_integrate(fm1,fm2)
  assert approx_equal(quad_value, laplace_value,eps=5e-2)
  #print "--->", quad_value, laplace_value, 100*(quad_value-laplace_value)/laplace_value


def plot_function(twin_lik, range_value=7):
  tmp_1_lbfgs = minim( twin_lik, 1, 1 )
  fm1=tmp_1_lbfgs.x[0]
  fm2=tmp_1_lbfgs.x[1]
  snd_der = twin_lik.dd_log_p_dd_f(fm1, fm2)
  is11 = snd_der[0]
  is22 = snd_der[1]
  is12 = snd_der[2]
  det = is11*is22-is12*is12
  s11= is22/det
  s22= is11/det
  s12= -is12/det
  s1= math.sqrt( abs(s11)*2.0 )
  s2= math.sqrt( abs(s22)*2.0 )
  for ii in range(-25,26):
    for jj in range(-25,26):
      x=range_value*ii/25.0 + fm1
      y=range_value*jj/25.0 + fm2
      print(x,y, twin_lik.log_p(x,y)-twin_lik.log_p(fm1,fm2))

def tst_single_likelihood(centric=False,
                          twin_fraction=0.1,
                          sigmaa=0.7,
                          plot=False):
  io1=0.6
  io2=0.6
  so1=1.0/100.0
  so2=1.0/100.0
  fc1=0.3
  fc2=0.3
  a  =sigmaa
  b  =1.0-a*a
  centric=centric
  tf=twin_fraction
  tmp_1 = xray.single_twin_likelihood(io1,so1,io2,so2,
                                      fc1,fc2,1.0,1.0,
                                      centric,centric,
                                      a,b,
                                      tf, 90)
  tmp_1_de = find_max( tmp_1 )

  tmp_1_lbfgs = minim( tmp_1 ,tmp_1_de.x[0], tmp_1_de.x[1] )
  #Check if similar optimia are found
  #assert approx_equal( tmp_1_de.x[0], tmp_1_lbfgs.x[0], eps=5e-1 )
  #assert approx_equal( tmp_1_de.x[1], tmp_1_lbfgs.x[1], eps=5e-1 )

  #test second derivatives
  tst_snd_der(tmp_1, 3.0, 2.0 )
  tst_snd_der(tmp_1, 2.0, 3.0 )

  #test the integration
  tst_integration( tmp_1 )
  if plot:
    plot_function( tmp_1 )


def tst_twin_completion():
  uc = uctbx.unit_cell( "40,40,70,90,90,90" )
  xs = crystal.symmetry( unit_cell=uc, space_group="P1" )
  miller_set = miller.build_set( crystal_symmetry=xs,
                                 anomalous_flag=False,
                                 d_min=3.0 ).map_to_asu()
  select = flex.bool( miller_set.indices().size(), True )
  select[300]=False
  miller_set_mod = miller_set.select( select )
  # make sure we threw away a reflection
  assert not miller_set_mod.indices().all_eq(  miller_set.indices() )
  checker = xray.twin_completion(
    miller_set_mod.indices(),
    xs.space_group(),
    False,
    [0,1,0,1,0,0,0,0,-1]
    )
  new_hkl = checker.twin_complete()
  miller_set_mod = miller_set_mod.customized_copy(
    indices=new_hkl)
  miller_set_mod = miller_set_mod.map_to_asu()
  a,b = miller_set_mod.common_sets(  miller_set )
  assert  a.indices().size() == miller_set_mod.indices().size()
  assert  a.indices().size() == miller_set.indices().size()
  assert  miller_set_mod.is_unique_set_under_symmetry()

  checker = xray.twin_completion(
    miller_set.indices(),
    xs.space_group(),
    False,
    [0,1,0,1,0,0,0,0,-1]
    )

  basic_flags = miller_set.generate_r_free_flags_basic()
  lattice_flags = miller_set.generate_r_free_flags_on_lattice_symmetry()

  assert not checker.check_free_flags( basic_flags.data() ) # this should give False
  assert checker.check_free_flags( lattice_flags.data() )   # this should give True

  selection_array = checker.get_free_model_selection( lattice_flags.indices(), lattice_flags.data() )
  assert selection_array.all_eq( lattice_flags.data()  )


def tst_detwin():
  tmp = rs.xray_structure(sgtbx.space_group_info( 'P4' ),
                          elements=['C']*310,
                          n_scatterers=310)

  sfs = abs( tmp.structure_factors( False, 2.5  ).f_calc() )
  tmp_detwin = xray.hemihedral_detwinner(
     hkl_obs = sfs.indices() ,
     hkl_calc= sfs.indices(),
     space_group = sfs.space_group(),
     anomalous_flag = False,
     twin_law = [-1,0,0,  0,1,0,  0,0,-1] )#"-h,k,-l" )

  sfs = sfs.f_as_f_sq()

  tf = [0.1,0.2,0.3,0.4,0.49]
  for t in tf:
    i,s = tmp_detwin.twin_with_twin_fraction( i_obs = sfs.data(),
                                              sigma_obs = sfs.sigmas(),
                                              twin_fraction = t )
    diff = sfs.data() - i
    diff = flex.sum( flex.abs(diff) )
    assert diff>1e-3
    #check what happen when no sigmas ae passed in
    dti,dts = tmp_detwin.detwin_with_twin_fraction( i_obs = i,
                                                    sigma_obs = None,
                                                    twin_fraction= t )


    dti,dts = tmp_detwin.detwin_with_twin_fraction( i_obs = i,
                                                    sigma_obs = s,
                                                    twin_fraction= t )
    diff = sfs.data() - dti
    diff = flex.sum( flex.abs(diff) ) # / flex.sum( sfs.data() )
    assert approx_equal( diff, 0, eps=1e-5 )
    permut = tmp_detwin.obs_to_twin_obs()
    ind = range( permut.size() )
    for ii, jj, kk, pp, mm  in zip( sfs.data(), i, dti, ind, permut ):
      no = (1-t)*sfs.data()[pp] + t*sfs.data()[mm]
      assert approx_equal( jj-no, 0, eps=1e-5)
      assert approx_equal( ii-kk, 0, eps=1e-5)
      #print sfs.indices()[pp], sfs.indices()[ mm ] , ii, jj, kk, pp,mm, jj-no


def run():
  tst_detwin()
  tst_twin_completion()

  tst_ls_on_i()
  tst_ls_on_f()

  tst_single_likelihood(True,0.15,0.0)
  tst_single_likelihood(True,0.15,0.1)
  tst_single_likelihood(True,0.15,0.4)
  tst_single_likelihood(True,0.15,0.6)
  tst_single_likelihood(True,0.15,0.8)

  tst_single_likelihood(True,0.25,0.0)
  tst_single_likelihood(True,0.25,0.1)
  tst_single_likelihood(True,0.25,0.4)
  tst_single_likelihood(True,0.25,0.6)
  tst_single_likelihood(True,0.25,0.8)

  tst_single_likelihood(True,0.40,0.0)
  tst_single_likelihood(True,0.40,0.1)
  tst_single_likelihood(True,0.40,0.4)
  tst_single_likelihood(True,0.40,0.6)
  tst_single_likelihood(True,0.40,0.8)

  tst_single_likelihood(True,0.44,0.0)
  tst_single_likelihood(True,0.44,0.1)
  tst_single_likelihood(True,0.44,0.4)
  tst_single_likelihood(True,0.44,0.6)
  tst_single_likelihood(True,0.44,0.8)

  tst_single_likelihood(False,0.15,0.0)
  tst_single_likelihood(False,0.15,0.1)
  tst_single_likelihood(False,0.15,0.4)
  tst_single_likelihood(False,0.15,0.6)
  tst_single_likelihood(False,0.15,0.8)

  tst_single_likelihood(False,0.25,0.0)
  tst_single_likelihood(False,0.25,0.1)
  tst_single_likelihood(False,0.25,0.4)
  tst_single_likelihood(False,0.25,0.6)
  tst_single_likelihood(False,0.25,0.8)

  tst_single_likelihood(False,0.40,0.0)
  tst_single_likelihood(False,0.40,0.1)
  tst_single_likelihood(False,0.40,0.4)
  tst_single_likelihood(False,0.40,0.6)
  tst_single_likelihood(False,0.40,0.8)

  tst_single_likelihood(False,0.44,0.0)
  tst_single_likelihood(False,0.44,0.1)
  tst_single_likelihood(False,0.44,0.4)
  tst_single_likelihood(False,0.44,0.6)
  tst_single_likelihood(False,0.44,0.8)

if (__name__ == "__main__"):
  run()
  print(format_cpu_times())
