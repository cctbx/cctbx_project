from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import adptbx
from scitbx import lbfgs
from mmtbx import masks
import cctbx.xray.structure_factors
from cctbx.eltbx.xray_scattering import wk1995
from libtbx import adopt_init_args
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
from libtbx.math_utils import iround
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import mmtbx.scaling
import scitbx.math as sm
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import sigmaa_estimation
from mmtbx import masks
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
import mmtbx.f_model
from libtbx import table_utils
from libtbx.utils import Sorry, user_plus_sys_time
import scitbx.lbfgs
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from scitbx import differential_evolution
import sys, os, math, time, string
import mmtbx.f_model
from libtbx.str_utils import format_value, show_string
from mmtbx.scaling import outlier_rejection
from mmtbx.scaling import absolute_scaling
import mmtbx.scaling.twin_analyses
from libtbx import Auto

master_params =  iotbx.phil.parse("""
  twin_law = None
  .type=str
  detwin{
    mode = algebraic proportional *auto
    .type= choice
    map_types{
      twofofc = *two_m_dtfo_d_fc two_dtfo_fc
      .type = choice
      fofc = *m_dtfo_d_fc gradient m_gradient
      .type = choice
      aniso_correct = False
      .type=bool

    }
  }
  """)


class twin_fraction_object(object):
  """provides methods for derivatives and
  transformastion of twin fraction"""
  def __init__(self, twin_fraction = 0):
    self.min_frac = 0.001
    self.max_frac = 0.999
    self.twin_fraction = float(twin_fraction)
    if (self.twin_fraction<=self.min_frac):
      self.twin_fraction = self.min_frac + self.min_frac/10.0
    if (self.twin_fraction>=self.max_frac):
      self.twin_fraction = self.max_frac - self.min_frac/10.0

  def twin_fraction_to_ref( self ):
    tmp = self.twin_fraction - self.min_frac
    tmp = (self.max_frac-self.min_frac)/tmp -1.0
    if tmp < 1e-70:
      tmp = 1e-70
    tmp = -math.log( tmp )
    return tmp

  def ref_to_twin_fraction(self, x):
    if (x<-10):
      x=-10
    tmp = self.min_frac + (self.max_frac-self.min_frac)/(1+math.exp(-x) )
    self.twin_fraction = tmp

  def d_t_d_twin_fraction_ref(self, dtdp ):
    tmp = self.twin_fraction_to_ref()
    tmp2 = 1.0+math.exp(-tmp )
    tmp = (self.max_frac - self.min_frac)*math.exp( -tmp )/(tmp2*tmp2)
    # (d target)/(d twin_fraction)* (d twin_fraction)/(d refinable parameter)
    # |--------------------------|  |---------------------------------------|
    #         from outside                       calculated above
    return dtdp*tmp

  def show(self,out=None):
    if out is None:
      out = sys.stdout
    print >> out, "twin fraction: %4.3f" %( self.twin_fraction  )



class scaling_parameters_object(object):
  """Object holds a set of parameters needed for f model.
  provides tranformations for parameter optimisation"""
  def __init__(self,
               xs=None,
               k_overall=1.0,
               u_star=(0,0,0,0,0,0),
               k_sol=0,
               u_sol=0,
               k_part=0,
               u_part=0,
               object=None,
               max_u_sol  = 5.0,
               max_u_part = 5.0 ):

    if object is not None:
      k_overall = object.k_overall
      u_star    = object.u_star
      k_sol     = object.k_sol
      u_sol     = object.u_sol
      k_part    = object.k_part
      u_part    = object.u_part
      u_star    = object.u_star
    # this is complete paranoia,. Trying to ensure that one always obtains 2 unique objects.
    self.k_overall = float(k_overall)
    self.u_part = float(u_part)
    self.k_sol = float(k_sol)
    self.u_sol = float(u_sol)
    self.k_part = float(k_part)
    self.u_star = ( float(u_star[0]),
                    float(u_star[1]),
                    float(u_star[2]),
                    float(u_star[3]),
                    float(u_star[4]),
                    float(u_star[5])
                  )
    self.max_u_sol = float(max_u_sol)
    self.max_u_part= float(max_u_part)
    if xs is None:
      self.xs=object.xs
    else:
      self.xs=xs

    assert self.xs is not None

    self.adp_constraints = self.xs.space_group().adp_constraints()
    self.vrwgk =  math.pow(self.xs.unit_cell().volume(),-2.0/3.0)
    self.n_u_indep = self.xs.space_group().adp_constraints(
      ).n_independent_params()

    # make sure that the supplied adp follows the symmetry constraints
    self.u_star = self.xs.space_group().average_u_star( self.u_star )

  def ref_to_k_overall(self,x):
    self.k_overall = math.exp( x )

  def ref_to_k_sol(self,x):
    if x>10:
      self.k_sol = math.exp( 10 )
    else:
      self.k_sol = math.exp( x )

  def ref_to_u_sol(self, x):
    if x>10:
      self.u_sol = math.exp(10.0)
    else:
      self.u_sol = math.exp( x )
    if self.u_sol > self.max_u_sol:
      self.u_sol = self.max_u_sol

  def ref_to_k_part(self, x):
    if x > 10:
      self.k_part = math.exp(10)
    else:
      self.k_part = math.exp( x )

  def ref_to_u_part(self, x):
    self.u_part = math.exp( x )
    if self.u_part > self.max_u_part:
      self.u_part = self.max_u_part

  def ref_to_u_star(self, x ):
    # first we need to expand the bugger to the full size
    tmp =  self.adp_constraints.all_params( x )
    # now it needs to be scaled back to something
    # physical
    tmp =list( flex.double(tmp) * self.vrwgk )
    # done
    self.u_star = tmp

  def k_overall_to_ref(self):
    if self.k_overall > 0:
      return math.log( self.k_overall )
    else:
      return None
  def k_sol_to_ref(self):
    if self.k_sol>0:
      return math.log( self.k_sol )
    else:
      return 0
  def k_part_to_ref(self):
    if self.k_part > 0:
      return math.log( self.k_part )
    else:
      return 0
  def u_sol_to_ref(self):
    if self.u_sol > 0:
      return math.log( self.u_sol )
    else:
     return -1000.0

  def u_part_to_ref(self):
    if self.u_part>0:
      return math.log( self.u_part )
    else:
      return -1000.0

  def u_star_to_ref(self):
    # first we pick the independent parameters
    tmp = self.xs.space_group().adp_constraints(
      ).independent_params(all_params=self.u_star)
    # now do the scaling please
    tmp =  list( flex.double(tmp)/self.vrwgk )
    return tmp

  # derivatives of refinable parameter wrst to target
  def d_t_d_k_overall_ref(self,dtdp):
    return self.k_overall*dtdp
  def d_t_d_k_sol_ref(self,dtdp):
    return self.k_sol*dtdp
  def d_t_d_k_part_ref(self,dtdp):
    return self.k_part*dtdp
  def d_t_d_u_sol_ref(self, dtdp):
    return self.u_sol*dtdp
  def d_t_d_u_part_ref(self, dtdp):
    return self.u_part*dtdp
  def d_t_d_u_star_ref(self, dtdp):
    # first introduce the scaling
    tmp = list( flex.double(dtdp) * self.vrwgk )
    #now do the symmetry completion please
    tmp = list( self.adp_constraints.independent_gradients(
      list( tmp ) ) )
    return tmp

  def show(self,out=None):
    if out is None:
      out=sys.stdout
    print >> out
    print >> out, "F-model scaling parameters"
    print >> out, "k_overall : %5.2e"%(self.k_overall)
    print >> out, "u_star    : %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e"%(
      self.u_star[0], self.u_star[1], self.u_star[2],
      self.u_star[3], self.u_star[4], self.u_star[5])
    print >> out, "   (%i independent parameters)"%(self.n_u_indep)
    print >> out, "k_sol     : %5.2e"%(self.k_sol)
    print >> out, "u_sol     : %5.2e"%(self.u_sol)
    print >> out, "    B_sol : %5.2f"%(self.u_sol*79.0)
    print >> out, "k_part    : %5.2e"%(self.k_part)
    print >> out, "u_part    : %5.2e"%(self.u_part)
    print >> out


  def deep_copy(self):
    new = scaling_parameters_object(object=self)
    return new



def get_initial_scale(miller_obs,
                      f_atoms):
  tmp_calc = f_atoms.deep_copy().map_to_asu()
  tmp_obs = miller_obs.deep_copy().map_to_asu()
  tmp_calc, tmp_obs = tmp_obs.common_sets(
    abs(tmp_calc)  )
  init_scale = flex.sum( tmp_calc.data()*tmp_obs.data() )/ \
               flex.sum( tmp_calc.data()*tmp_calc.data() )
  return init_scale

class de_bulk_solvent_scaler(object):
  def __init__(self,
               scaling_parameters,
               twin_fraction_obj,
               target_evaluator,
               f_model_core_data,
               fix_sol      = False,
               fix_aniso    = False,
               fix_twin     = False,
               fix_scale    = False,
               best_sol     = None,
               best_aniso   = None,
               best_twin    = None,
               best_scale   = None,
               twin_domain  = (-1,1),
               scale_domain = (-4,4),
               sol_domain   = [(-2,0),(-2,0)],
               aniso_domain = (-1,1) ,
               out          = None):
    self.out = out
    if self.out is None:
      self.out = sys.stdout

    self.fix_sol = fix_sol
    self.fix_aniso = fix_aniso
    self.fix_twin = fix_twin

    self.scaling_parameters=scaling_parameters
    self.target_evaluator=target_evaluator
    self.f_model_core_data=f_model_core_data
    self.twin_fraction_object = twin_fraction_obj
    self.best_score = None
    #first determin the number of parameters please
    self.n = 1+2+1+self.scaling_parameters.n_u_indep
    # parameter map
    # 0  : twin fraction
    # 1  : scale factor
    # 2  : ksol
    # 3  : usol
    # 4-N: ustar

    if fix_sol:
      assert len(best_sol)==2
      # this should be a list [ksol,usol]
    if fix_scale:
      best_scale >= 0

    if self.fix_aniso:
      assert len(self.best_aniso) == self.scaling_parameters.n_u_indep
      # these should be independent parameters for aniso u star!

    if self.fix_twin:
      assert self.best_twin >= 0.0
      assert self.best_twin <= 1.0

    if best_scale is not None:
      scale_domain = (math.log(best_scale)+scale_domain[0], math.log(best_scale)+scale_domain[1])
    x = 0
    if best_twin is not None:
      x = self.twin_fraction_to_ref( best_twin )
      twin_domain = ( x+twin_domain[0] , x+twin_domain[1] )

    if best_aniso is not None:
      tmp_aniso_domain = []
      for ii in xrange( self.scaling_parameters.n_u_indep ):
        tmp_aniso_domain.append( (best_aniso[ii]+aniso_domain[0], best_aniso[ii]+aniso_domain[1]) )
      aniso_domain = tmp_aniso_domain
    else:
      tmp_aniso_domain = []
      for ii in xrange( self.scaling_parameters.n_u_indep ):
        tmp_aniso_domain.append( (aniso_domain[0],aniso_domain[1]) )
      aniso_domain = tmp_aniso_domain

    if best_sol is not None:
      sol_domain = [ (math.log(best_sol[0]+1e-12)+sol_domain[0][0], math.log(best_sol[0]+1e-12)+sol_domain[0][1] ),
                     (math.log(best_sol[1]+1e-12)+sol_domain[1][0], math.log(best_sol[1]+1e-12)+sol_domain[1][1]) ]

    self.domain = [twin_domain]  + [scale_domain] + sol_domain + aniso_domain
    self.x = flex.double([0]*self.n)

    # amke a default solution vector
    solution_vector = [ 0, 0, math.log( 0.30), math.log(46.0/80.0) ]
    for ii in xrange(self.scaling_parameters.n_u_indep):
      solution_vector.append( 0 )
    if best_twin is not None:
      solution_vector[0]=x
    if best_scale is not None:
      solution_vector[1]=math.log(best_scale+1e-12)
    if best_sol is not None:
      solution_vector[2]=math.log(best_sol[0]+1e-12)
      solution_vector[3]=math.log(best_sol[1]+1e-12)
    if best_aniso is not None:
      for ii in xrange(self.scaling_parameters.n_u_indep):
        solution_vector[ii+4]=best_aniso[ii]

    self.de = differential_evolution.differential_evolution_optimizer(
     self,
     population_size=min(self.n+2,8),
     f=0.8,
     cr=0.7,
     n_cross=2,
     eps=1e-8,
     show_progress=False,
     max_iter = 500,
     insert_solution_vector = flex.double(solution_vector )
    )

  def twin_fraction_to_ref( self, tf ):
    tmp = tf - 0.001
    tmp = (0.999-0.001)/tmp -1.0
    if tmp < 1e-70:
      tmp = 1e-70
    tmp = -math.log( tmp )
    return tmp



  def update_parameters(self, vector):
    self.scaling_parameters.ref_to_k_overall( vector[1] )
    if self.fix_sol is not False:
      self.scaling_parameters.k_sol = self.fix_sol[0]
      self.scaling_parameters.u_sol = self.fix_sol[1]
    else:
      self.scaling_parameters.ref_to_k_sol( vector[2] )
      self.scaling_parameters.ref_to_u_sol( vector[3] )

    if self.fix_aniso is not False:
      self.scaling_parameters.u_star = list( self.fix_aniso)
    else:
      self.scaling_parameters.ref_to_u_star(   list(vector[4:])   )

    if self.fix_twin is not False:
      self.twin_fraction_object.twin_fraction = self.fix_twin
    else:
      self.twin_fraction_object.ref_to_twin_fraction( vector[0] )

  def target(self, vector):
    #first make sure our place holder for scaling parameters is up to date!
    self.update_parameters(vector)
    #self.scaling_parameters.show()
    #self.twin_fraction_object.show()

    #make the core data summat aware of the changes
    self.f_model_core_data.koverall(
      self.scaling_parameters.k_overall )
    self.f_model_core_data.ustar(
      self.scaling_parameters.u_star)
    self.f_model_core_data.ksol(
      self.scaling_parameters.k_sol )
    self.f_model_core_data.usol(
      self.scaling_parameters.u_sol )
    #do the same thing for the twin fraction please
    self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )

    # we can get the target value!
    f = self.target_evaluator.target( self.f_model_core_data.f_model() )
    if self.best_score is None:
      self.best_score = f
    else:
      if self.best_score > f:
        self.best_score = f

    return f

  def print_status(self,
                   best_score,
                   mean_score,
                   vector,
                   count=None):
    self.update_parameters(vector)
    b_cart = adptbx.u_star_as_u_cart(self.scaling_parameters.xs.unit_cell(),
                                     self.scaling_parameters.u_star)
    b_cart = adptbx.u_as_b( b_cart )
    print >> self.out, "#--------------------------------------------#"
    if count is not None:
      print >> self.out, "| Generation   :     %8s                |"%(count)
    print >> self.out, "| best score   :     %8.6e            |"%(best_score)
    print >> self.out, "| mean score   :     %8.6e            |"%(mean_score)
    print >> self.out, "| k_overall    :     %8.6e            |"%(
      self.scaling_parameters.k_overall)
    print >> self.out, "| b_cart(ii)   : %8s %8s %8s  |"%("%5.2f"%(b_cart[0]),
                                                       "%5.2f"%(b_cart[1]),
                                                       "%5.2f"%(b_cart[2]) )
    print >> self.out, "| b_cart(ij)   : %8s %8s %8s  |"%("%5.2f"%(b_cart[3]),
                                                       "%5.2f"%(b_cart[4]),
                                                       "%5.2f"%(b_cart[5]) )
    print >> self.out, "| k_sol, b_sol :     %4.3f  %5.2f            |"%( self.scaling_parameters.k_sol,
                                                            adptbx.u_as_b(self.scaling_parameters.u_sol)
                                                            )
    print >> self.out, "| twin fraction:     %4.3f                   |"%(
      self.twin_fraction_object.twin_fraction)
    print >> self.out, "#--------------------------------------------#"


class bulk_solvent_scaler(object):
  def __init__(self,
               scaling_parameters,
               twin_fraction_obj,
               target_evaluator,
               f_model_core_data,
               parameter_mask):
    self.parameter_mask = parameter_mask
    self.scaling_parameters=scaling_parameters
    self.target_evaluator=target_evaluator
    self.f_model_core_data=f_model_core_data
    self.twin_fraction_object = twin_fraction_obj

    #first determin the number of parameters please
    self.n = 1+2+1+self.scaling_parameters.n_u_indep
    self.x = flex.double([0]*self.n)
    self.f = None
    self.update_x()
    self.compute_functional_and_gradients()
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations = 30)
    exception_handling_parameters = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      ignore_line_search_failed_step_at_upper_bound=True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters,
      exception_handling_params=exception_handling_parameters)

  def update_x(self):
    # this is only need when starting,
    # all other updates are done by the minimizer
    self.x[0] = self.twin_fraction_object.twin_fraction_to_ref()
    self.x[1] = self.scaling_parameters.k_overall_to_ref()
    self.x[2] = self.scaling_parameters.k_sol_to_ref()
    self.x[3] = self.scaling_parameters.u_sol_to_ref()
    tmp       = self.scaling_parameters.u_star_to_ref()
    for item,ii in zip(tmp,xrange(4,self.n)) :
      self.x[ii]=item

  def update_parameter_objects(self):
    """ from refinebale parameters to physical meaning full params"""
    self.scaling_parameters.ref_to_k_overall( self.x[1] )
    self.scaling_parameters.ref_to_k_sol( self.x[2] )
    self.scaling_parameters.ref_to_u_sol( self.x[3] )
    self.scaling_parameters.ref_to_u_star(   list(self.x[4:])   )
    self.twin_fraction_object.ref_to_twin_fraction( self.x[0] )

  def update_core_data(self):
    #make the core data summat aware of the changes
    self.f_model_core_data.koverall(
      self.scaling_parameters.k_overall )
    self.f_model_core_data.ustar(
      self.scaling_parameters.u_star)
    self.f_model_core_data.ksol(
      self.scaling_parameters.k_sol )
    self.f_model_core_data.usol(
      self.scaling_parameters.u_sol )
    #do the same thing for the twin fraction please
    self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )

  def update_parameters(self):
    self.update_parameter_objects()
    self.update_core_data()

  def compute_functional_and_gradients(self):
    #first make sure our place holder for scaling parameters is up to date!
    self.update_parameters()

    f = self.compute_functional()
    g = flex.double(self.compute_gradients())
    #gfd = self.compute_gradients_fin_diff()
    self.f = f
    return f, g

  def compute_functional(self):
    target = self.target_evaluator.target( self.f_model_core_data.f_model() )
    return target

  def compute_gradients(self):
    dtdab =  self.target_evaluator.d_target_d_ab(
      self.f_model_core_data.f_model() )
    gradient_flags = [True,True,True,False,True,False]
    gradient_object =  self.f_model_core_data.d_target_d_all(
      dtdab[0], dtdab[0], flex.bool(gradient_flags) )
    dtdalpha = self.target_evaluator.d_target_d_alpha(
      self.f_model_core_data.f_model())
    g = [0,0,0,0]
    g[0] = self.twin_fraction_object.d_t_d_twin_fraction_ref( dtdalpha )*self.parameter_mask.twin_fraction
    g[1] = self.scaling_parameters.d_t_d_k_overall_ref(
      gradient_object.koverall() )*self.parameter_mask.k_overall
    g[2] = self.scaling_parameters.d_t_d_k_sol_ref(
      gradient_object.ksol() )*self.parameter_mask.k_sol
    g[3] = self.scaling_parameters.d_t_d_u_sol_ref(
      gradient_object.usol() )*self.parameter_mask.u_sol
    tmp = self.scaling_parameters.d_t_d_u_star_ref(
      gradient_object.ustar() )
    for aa in tmp:
      g.append(aa*self.parameter_mask.u_star)
    return g

class scaling_parameter_mask(object):
  def __init__(self,
               twin_fraction=True,
               k_overall=True,
               u_star=True,
               k_sol=True,
               u_sol=True):
    self.twin_fraction = 0.0
    self.k_overall     = 0.0
    self.u_star        = 0.0
    self.k_sol         = 0.0
    self.u_sol         = 0.0
    if twin_fraction:
      self.twin_fraction = 1.0
    if k_overall:
      self.k_overall     = 1.0
    if u_star:
      self.u_star        = 1.0
    if k_sol:
      self.k_sol         = 1.0
    if u_sol:
      self.u_sol         = 1.0



class bulk_solvent_scaling_manager(object):
  def __init__(self,
               target_evaluator,
               f_model_core_data,
               crystal_symmetry,
               scaling_parameters=None,
               twin_fraction=None,
               k_sol_limits=(0.1, 0.8),
               u_sol_limits=(10.0/80.0,80.0/80.0),
               n_trials=100,
               out=None):

    self.out = out
    if out is None:
      self.out = sys.stdout
    self.crystal_symmetry = crystal_symmetry

    # value used currently by minimizers
    self.scaling_parameters = scaling_parameters
    self.twin_fraction = twin_fraction_object( twin_fraction.twin_fraction )

    # cached best values
    self.best_scaling_parameters = scaling_parameters
    self.best_twin_fraction = twin_fraction_object( twin_fraction.twin_fraction )
    self.best_score_until_now = None # will be filled in later

    self.target_evaluator = target_evaluator
    self.f_model_core_data = f_model_core_data
    self.n_trials = n_trials

    # setting up a grid for sampling
    self.k_u_grid = sm.square_halton_sampling(
      k_sol_limits[0], k_sol_limits[1],
      u_sol_limits[0], u_sol_limits[1] )

  def print_it(self,
               score,
               scaling_params,
               twin_fractions):
    b_cart = adptbx.u_star_as_u_cart(scaling_params.xs.unit_cell(),
                                     scaling_params.u_star)
    b_cart = adptbx.u_as_b( b_cart )

    print >> self.out
    print >> self.out, "#--------------------------------------------#"
    if score is not None:
      print >> self.out, "| score        :     %8.6e            |"%(score)
    print >> self.out, "| k_overall    :     %8.6e            |"%(
      scaling_params.k_overall)
    print >> self.out, "| b_cart(ii)   : %8s %8s %8s  |"%("%5.2f"%(b_cart[0]),
                                                       "%5.2f"%(b_cart[1]),
                                                       "%5.2f"%(b_cart[2]) )
    print >> self.out, "| b_cart(ij)   : %8s %8s %8s  |"%("%5.2f"%(b_cart[3]),
                                                       "%5.2f"%(b_cart[4]),
                                                       "%5.2f"%(b_cart[5]) )
    print >> self.out, "| k_sol, b_sol :     %4.3f  %5.2f            |"%( scaling_params.k_sol,
                                                            adptbx.u_as_b(scaling_params.u_sol)
                                                            )
    print >> self.out, "| twin fraction:     %4.3f                   |"%(
      twin_fractions.twin_fraction)
    print >> self.out, "#--------------------------------------------#"
    print >> self.out

  def print_best(self):
    self.print_it( self.best_score_until_now, self.best_scaling_parameters, self.best_twin_fraction)


  def initial_scale_and_twin_fraction(self):
    # ksol and usol are fixed
    scaling_parameters = scaling_parameters_object(
      xs    = self.crystal_symmetry,
      k_sol = 0.40,
      u_sol = 50.0/80.0)
    twin_fraction = twin_fraction_object(twin_fraction=0.10)
    # first refine the scale fcator
    parameter_mask = scaling_parameter_mask(twin_fraction=False,
                                            k_overall=True,
                                            u_star=False,
                                            k_sol=False,
                                            u_sol=False)
    scaler = bulk_solvent_scaler(
      scaling_parameters,
      twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)
    # now please refine the twin fraction
    parameter_mask = scaling_parameter_mask(twin_fraction=True,
                                            k_overall=False,
                                            u_star=False,
                                            k_sol=False,
                                            u_sol=False)
    scaler = bulk_solvent_scaler(
      scaling_parameters,
      twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)
    # refine twin and scale factor simulataneously
    parameter_mask = scaling_parameter_mask(twin_fraction=True,
                                            k_overall=True,
                                            u_star=False,
                                            k_sol=False,
                                            u_sol=False)
    scaler = bulk_solvent_scaler(
      scaling_parameters,
      twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)

    # set the scaling parameters please
    self.scaling_parameters.k_overall = scaling_parameters.k_overall
    self.twin_fraction.twin_fraction = twin_fraction.twin_fraction

    #check if a 'best score' is allready in place
    if self.best_score_until_now is None:
      self.best_score_until_now = scaler.f
      self.best_scaling_parameters.k_overall = scaling_parameters.k_overall
      self.best_twin_fraction.twin_fraction = twin_fraction.twin_fraction
    else:
      if self.best_score_until_now < scaler.f:
        self.best_score_until_now = scaler.f
        self.best_scaling_parameters.k_overall = scaling_parameters.k_overall
        self.best_twin_fraction.twin_fraction = twin_fraction.twin_fraction

  def update_best(self):
    self.best_scaling_parameters = scaling_parameters_object( object = self.scaling_parameters )
    self.best_twin_fraction.twin_fraction = float(self.twin_fraction.twin_fraction)

  def grid_search(self,n_cycle="Auto"):

    converged=False
    cycle_count=0
    last_score = None
    while not converged:
      self.setup_next_trial(reset=True)
      score = self.target_evaluator.target(
        self.f_model_core_data.f_model() )
      min_score = score
      k_b = (self.scaling_parameters.k_sol,
             self.scaling_parameters.u_sol)
      for ii in xrange(self.n_trials-1):
        self.setup_next_trial()
        score = self.target_evaluator.target(
          self.f_model_core_data.f_model() )
        if score < min_score:
          min_score = score
          k_b = (self.scaling_parameters.k_sol,
                 self.scaling_parameters.u_sol)
          if min_score < self.best_score_until_now:
            #self.print_it(score,
            #              self.scaling_parameters,
            #              self.twin_fraction)
            self.best_score_until_now = score
            self.update_best()
      self.scale_it()
      cycle_count += 1

      if n_cycle is not "Auto":
        if cycle_count == n_cycle:
          converged=True
      else:
        if last_score is not None:
          if self.best_score_until_now <= last_score:
            converged = True
          else:
            last_score = self.best_score_until_now
        else:
          last_score = self.best_score_until_now



    #finsih up
    self.scale_it()
    score = self.target_evaluator.target(
      self.f_model_core_data.f_model() )
    if score < self.best_score_until_now:
      # a new global minimum is found
      self.print_it(score,
                    self.scaling_parameters,
                    self.twin_fraction)
      # update the best parameters please


  def de_search(self):
    best_sol = None
    if min(self.scaling_parameters.k_sol, self.scaling_parameters.u_sol)<=1E-5:
      best_sol = None
    else:
      best_sol = (self.scaling_parameters.k_sol, self.scaling_parameters.u_sol)
    de = de_bulk_solvent_scaler(
               self.scaling_parameters,
               self.twin_fraction,
               self.target_evaluator,
               self.f_model_core_data,
               out = self.out,
               best_twin = self.twin_fraction.twin_fraction,
               best_scale = self.scaling_parameters.k_overall,
               best_aniso = self.scaling_parameters.u_star_to_ref(),
               best_sol   = best_sol,
               scale_domain=(-0.1, 0.1 ),
               sol_domain = [ (-0.2, 0.2), (-0.1,0.1) ],
               aniso_domain = (-0.01, 0.01 ),
               #fix_twin=0.05,
               #fix_sol=[0.2,50/80.0]
               )
    self.scaling_parameters = de.scaling_parameters
    self.twin_fraction = de.twin_fraction_object
    if (de.best_score < self.best_score_until_now) or (self.best_score_until_now is None):
      self.update_best()

  def scale_it(self):
    parameter_mask = scaling_parameter_mask(twin_fraction=True,
                                            k_overall=True,
                                            u_star=False,
                                            k_sol=False,
                                            u_sol=False)
    scaler = bulk_solvent_scaler(
      self.scaling_parameters,
      self.twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)
    self.scaling_parameters=scaler.scaling_parameters
    self.twin_fraction=scaler.twin_fraction_object

    parameter_mask = scaling_parameter_mask(twin_fraction=False,
                                            k_overall=False,
                                            u_star=True,
                                            k_sol=False,
                                            u_sol=False)
    scaler = bulk_solvent_scaler(
      self.scaling_parameters,
      self.twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)
    self.scaling_parameters=scaler.scaling_parameters
    self.twin_fraction=scaler.twin_fraction_object

    parameter_mask = scaling_parameter_mask(twin_fraction=True,
                                            k_overall=True,
                                            u_star=True,
                                            k_sol=True,
                                            u_sol=True)
    scaler = bulk_solvent_scaler(
      self.scaling_parameters,
      self.twin_fraction,
      self.target_evaluator,
      self.f_model_core_data,
      parameter_mask)
    self.scaling_parameters=scaler.scaling_parameters
    self.twin_fraction=scaler.twin_fraction_object



  def setup_next_trial(self,reset=False):
    # 1a. set the overall scale factors from the best guess we have until now
    self.f_model_core_data.koverall(
      self.best_scaling_parameters.k_overall )
    self.f_model_core_data.ustar(
      self.best_scaling_parameters.u_star)
    # 1b. same goes for the twin fraction
    self.target_evaluator.alpha(self.best_twin_fraction.twin_fraction)


    # 2. get the values for k_sol and u_sol from the halton grid
    k_sol = None
    u_sol = None
    if not reset:
      k_sol, u_sol = self.k_u_grid.next()
    else:
      k_sol, u_sol = self.k_u_grid.start()

    self.f_model_core_data.ksol( k_sol )
    self.f_model_core_data.usol( u_sol )

    # some copying is needed unfortunately
    self.scaling_parameters.k_overall = self.best_scaling_parameters.k_overall
    self.scaling_parameters.u_star    = self.best_scaling_parameters.u_star
    self.scaling_parameters.k_sol     = k_sol
    self.scaling_parameters.u_sol     = u_sol
    self.twin_fraction.twin_fraction  = self.best_twin_fraction.twin_fraction

class target_attributes(mmtbx.f_model.target_attributes):

  def __init__(self):
    mmtbx.f_model.target_attributes.__init__(self, family="ls")
    self.twin = "amplitudes"
    self.pseudo_ml = False

class twin_model_manager(mmtbx.f_model.manager_mixin):
  def __init__(self,
               f_obs              = None,
               f_mask             = None,
               free_array         = None,
               xray_structure     = None,
               scaling_parameters = None,
               sf_and_grads_accuracy_params =
                           mmtbx.f_model.sf_and_grads_accuracy_params.extract(),
               mask_params        = None,
               out                = None,
               twin_law           = None,
               start_fraction     = 0.1,
               n_refl_bin         = 2000,
               max_bins           = 20,
               detwin_mode = "auto",
               map_types = master_params.extract().detwin.map_types
                ):
    self.alpha_beta_params=None
    self.twin = True
    self.sfg_params = sf_and_grads_accuracy_params
    self.target_name="twin_lsq_f"
    self._target_attributes = target_attributes()
    self.out = out
    self.did_search = 0

    if self.out is None:
      self.out = sys.stdout
    self.twin_fraction_object = twin_fraction_object(twin_fraction=start_fraction)
    self.twin_law=twin_law
    self.twin_fraction=start_fraction

    self.possible_detwin_modes = ["proportional",
                                  "algebraic",
                                  "gradient",
                                  "auto"
                                  ]
    self.detwin_mode = detwin_mode
    if self.detwin_mode is Auto:
      self.detwin_mode="auto"
    assert self.detwin_mode in self.possible_detwin_modes
    self.detwin_switch_twin_fraction = 0.45

    self.map_types = map_types

    assert (self.twin_law is not None)
    self.f_obs = f_obs.map_to_asu()
    self.free_array = free_array.map_to_asu()

    assert self.f_obs.indices().all_eq( self.free_array.indices() )

    self.f_obs_w = self.f_obs.select( ~self.free_array.data() )
    self.f_obs_f = self.f_obs.select( self.free_array.data() )

    #setup the binners if this has not been done yet
    self.max_bins = max_bins
    self.n_refl_bin = n_refl_bin
    if (self.n_refl_bin>self.f_obs.data().size() ) or (self.n_refl_bin is None):
      self.n_refl_bin = self.f_obs.data().size()
    if self.f_obs.binner() is None:
      if self.f_obs.indices().size()/float(n_refl_bin) > max_bins:
        self.f_obs.setup_binner(n_bins = max_bins)
      else:
        self.f_obs.setup_binner( reflections_per_bin=self.n_refl_bin )

    self.f_obs_w.use_binning_of( self.f_obs )
    self.f_obs_f.use_binning_of( self.f_obs )

    self.xray_structure = xray_structure
    self.xs = crystal.symmetry( unit_cell=f_obs.unit_cell(),
                                space_group=f_obs.space_group() )
    self.scaling_parameters = scaling_parameters_object(
      xs = self.xs,
      object = scaling_parameters)
    if self.scaling_parameters is None:
      self.scaling_parameters = scaling_parameters_object(self.xs)

    self.mask_params=None
    if mask_params is not None:
      self.mask_params = mask_params
    else:
      self.mask_params = mmtbx.masks.mask_master_params.extract()


    self.norma_sum_f_sq = flex.sum( self.f_obs.data() * self.f_obs.data() )
    self.norma_sum_f_sq_w = flex.sum( self.f_obs_w.data() * self.f_obs_w.data() )
    self.norma_sum_f_sq_f = flex.sum( self.f_obs_f.data() * self.f_obs_f.data() )
    #-------------------
    self.miller_set = None
    self.f_atoms = None
    self.free_flags_for_f_atoms = None
    self.miller_set = None
    self.f_atoms = self.compute_f_atoms()

    #-------------------
    self.f_mask_array = None
    if f_mask is not None:
      if  f_mask.data().size() == self.f_atoms.data().size():
        self.f_mask_array = f_mask
      else:
        self.update_f_mask()
    else:
      self.update_f_mask()

    #-------------------
    self.f_partial_array = None

    #-------------------
    self.data_core = xray.f_model_core_data(
      hkl = self.f_atoms.indices(),
      f_atoms= self.f_atoms.data(),
      f_mask = self.f_mask_array.data(),
      unit_cell = self.f_atoms.unit_cell(),
      k_overall=self.scaling_parameters.k_overall,
      u_star=self.scaling_parameters.u_star,
      k_sol=self.scaling_parameters.k_sol,
      u_sol=self.scaling_parameters.u_sol,
      f_part=None,
      k_part=self.scaling_parameters.k_part,
      u_part=self.scaling_parameters.u_part )


    self.target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs       = self.f_obs_w.indices(),
      f_obs         = self.f_obs_w.data(),
      w_obs         = self.f_obs_w.sigmas(),
      hkl_calc      = self.f_atoms.indices(),
      space_group   = self.f_obs.space_group(),
      anomalous_flag= self.f_obs.anomalous_flag(),
      alpha         = self.twin_fraction,
      twin_law      = self.twin_law.as_double_array()[0:9] )

    self.free_target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs        = self.f_obs_f.indices(),
      f_obs          = self.f_obs_f.data(),
      w_obs          = self.f_obs_f.sigmas(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs.space_group(),
      anomalous_flag = self.f_obs.anomalous_flag(),
      alpha          = self.twin_fraction,
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.bss=bulk_solvent_scaling_manager(
      self.target_evaluator,
      self.data_core,
      self.xs,
      self.scaling_parameters,
      self.twin_fraction_object,
      n_trials=1000,
      out=self.out)
    self.scaling_parameters = self.bss.best_scaling_parameters
    self.twin_fraction_object = self.bss.best_twin_fraction
    ###

    self.r_all_object = xray.hemihedral_r_values(
      hkl_obs        = self.f_obs.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs.space_group(),
      anomalous_flag = self.f_obs.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.r_work_object = xray.hemihedral_r_values(
      hkl_obs        = self.f_obs_w.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_w.space_group(),
      anomalous_flag = self.f_obs.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.r_free_object = xray.hemihedral_r_values(
      hkl_obs        = self.f_obs_f.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_f.space_group(),
      anomalous_flag = self.f_obs_f.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.work_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs_w.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_w.space_group(),
      anomalous_flag = self.f_obs_w.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )
    self.free_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs_f.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_f.space_group(),
      anomalous_flag = self.f_obs_f.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )
    self.full_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs.space_group(),
      anomalous_flag = self.f_obs.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.structure_factor_gradients_w = cctbx.xray.structure_factors.gradients(
      miller_set                   = self.miller_set,
      cos_sin_table                = self.sfg_params.cos_sin_table,
      grid_resolution_factor       = self.sfg_params.grid_resolution_factor,
      quality_factor               = self.sfg_params.quality_factor,
      u_base                       = self.sfg_params.u_base,
      b_base                       = self.sfg_params.b_base,
      wing_cutoff                  = self.sfg_params.wing_cutoff,
      exp_table_one_over_step_size = self.sfg_params.exp_table_one_over_step_size)

    self.sigmaa_object_cache = None
    self.update_sigmaa_object = True

    self.xray_structure_mask_cache = None
    if self.xray_structure is not None:
      self.xray_structure_mask_cache = self.xray_structure.deep_copy_scatterers()

    self.epsilons_w = self.f_obs_w.epsilons().data().as_double()
    self.epsilons_f = self.f_obs_f.epsilons().data().as_double()

  def info(self, free_reflections_per_bin = 140, max_number_of_bins = 30):
    return mmtbx.f_model.info(
      fmodel                   = self,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins       = max_number_of_bins)

  def outlier_selection(self, show = False, log = None):
    # XXX
    return None

  def remove_outliers(self, show = False, log = None):
    # XXX
    return self

  def wilson_b(self, force_update = False):
    # XXX
    return None

  def scale_k1(self):
    # XXX is it true scale_k1 ?
    # XXX name k_overall is ambigous.
    return self.scaling_parameters.k_overall

  def show_parameter_summary(self, manager=None):
    # print bulk solvent parameters
    print >> self.out, "Usol           ", self.scaling_parameters.u_sol, self.data_core.usol()
    print >> self.out, "Ksol           ", self.scaling_parameters.k_sol, self.data_core.ksol()
    print >> self.out, "Koverall       ", self.scaling_parameters.k_overall, self.data_core.koverall()
    print >> self.out, "Ustar          ", self.scaling_parameters.u_star, self.data_core.ustar()
    print >> self.out, "Twin fraction  ", self.twin_fraction_object.twin_fraction, self.twin_fraction, self.target_evaluator.alpha()
    print >> self.out, "mask step      ", self.mask_params.grid_step_factor
    print >> self.out, "mask shift     ", self.mask_params.mean_shift_for_mask_update
    print >> self.out, "mask trunk rad ", self.mask_params.shrink_truncation_radius
    print >> self.out, "mask solv rad  ", self.mask_params.solvent_radius

    if manager is not None:
      x = self.f_model().data()
      y = manager.f_model().data()
      print >> self.out, "Fmodel delta  " , flex.sum( flex.abs(x - y) )
      x = self.f_calc().data()
      y = manager.f_calc().data()
      print >> self.out, "Fatoms delta ", flex.sum( flex.abs(x - y) )
      x = self.f_mask_array.data()
      y = manager.f_mask_array.data()
      print >> self.out, "Fmask delta  ", flex.sum( flex.abs(x - y) )
      x = flex.abs( self.bulk_solvent_mask().data - manager.bulk_solvent_mask().data )
      print >> self.out, "Bit wise diff mask ", flex.sum( x )



  def deep_copy(self):
    new_object = twin_model_manager(
      f_obs              = self.f_obs.deep_copy(),
      f_mask             = self.f_mask_array.deep_copy(),
      free_array         = self.free_array.deep_copy(),
      xray_structure     = self.xray_structure.deep_copy_scatterers(),
      scaling_parameters = self.scaling_parameters.deep_copy(),
      mask_params        = self.mask_params,
      out                = self.out,
      twin_law           = self.twin_law,
      start_fraction     = self.twin_fraction,
      n_refl_bin         = self.n_refl_bin,
      max_bins           = self.max_bins,
      detwin_mode        = self.detwin_mode,
      map_types          = self.map_types,
      sf_and_grads_accuracy_params = self.sfg_params,
      )
    new_object.twin_fraction_object.twin_fraction = float(self.twin_fraction_object.twin_fraction)
    new_object.twin_fraction = float(self.twin_fraction_object.twin_fraction)
    new_object.update()
    new_object.did_search = self.did_search
    return new_object

  def resolution_filter(self,d_max=None,d_min=None):
    dc = self.deep_copy()
    # make a dummy copy of the observed data please
    dummy_obs  = dc.f_obs.resolution_filter(d_max,d_min)
    twin_complete = dc.construct_miller_set(external_miller_array = dummy_obs )
    appropriate_f_mask_array = dc.f_mask_array.common_set( twin_complete )
    new_object = twin_model_manager(
      f_obs        = dummy_obs,
      f_mask             = appropriate_f_mask_array, #dc.f_mask_array.resolution_filter(d_max,d_min),
      free_array         = dc.free_array.resolution_filter(d_max,d_min),
      xray_structure     = dc.xray_structure,
      scaling_parameters = dc.scaling_parameters.deep_copy(),
      mask_params        = dc.mask_params,
      out                = dc.out,
      twin_law           = dc.twin_law,
      start_fraction     = dc.twin_fraction,
      n_refl_bin         = dc.n_refl_bin,
      max_bins           = dc.max_bins,
      detwin_mode        = dc.detwin_mode,
      map_types          = dc.map_types,
      sf_and_grads_accuracy_params = dc.sfg_params
      )

    new_object.update()
    new_object.did_search = self.did_search
    return new_object

  def select(self, selection):
    dc = self.deep_copy()
    new_object = twin_model_manager(
      f_obs        = dc.f_obs.select(selection) ,
      f_mask             = dc.f_mask_array.selection(selection),
      free_array         = dc.free_array.selection(selection),
      xray_structure     = dc.xray_structure,
      scaling_parameters = dc.scaling_parameters.deep_copy(),
      mask_params        = dc.mask_params,
      out                = dc.out,
      twin_law           = dc.twin_law,
      start_fraction     = dc.twin_fraction,
      n_refl_bin         = dc.n_refl_bin,
      max_bins           = dc.max_bins,
      detwin_mode        = dc.detwin_mode,
      map_types          = dc.map_types,
      sf_and_grads_accuracy_params = dc.sfg_params
      )
    new_object.did_search = self.did_search
    return new_object


  def f_model(self):
    tmp_f_model = self.f_atoms.customized_copy(
      data = self.data_core.f_model()
    )
    return tmp_f_model

  def f_model_w(self):
    tmp = self.f_model()
    return tmp.select(~self.free_flags_for_f_atoms )

  def f_model_t(self):
    tmp = self.f_model()
    return tmp.select( self.free_flags_for_f_atoms )

  def f_calc(self):
    if self.f_atoms is None:
      self.f_atoms = self.compute_f_atoms()
    return self.f_atoms

  def f_calc_w(self):
    tmp = self.f_calc()
    return tmp.select(~self.free_flags_for_f_atoms )

  def f_calc_t(self):
    tmp = self.f_calc()
    return tmp.select( self.free_flags_for_f_atoms )

  def target_attributes(self):
    return self._target_attributes

  def update_solvent_and_scale(self,
                               params=None,
                               bulk_solvent_parameters=None,
                               twin_fraction_parameters=None,
                               refine=False,
                               grid_search=False,
                               initialise=False,
                               de_search=True,
                               out=None,
                               verbose=None
                               ):
    if initialise:
      self.bss.initial_scale_and_twin_fraction()
      self.scaling_parameters = self.bss.best_scaling_parameters
      self.twin_fraction_object = self.bss.best_twin_fraction
    if bulk_solvent_parameters is not None:
      self.scaling_parameters = bulk_solvent_parameters
      self.bss.best_scaling_parameters = self.scaling_parameters
    if twin_fraction_parameters is not None:
      self.bss.best_twin_fraction = self.twin_fraction_object
      self.twin_fraction_object = twin_fraction_parameters
    if de_search:
      self.did_search += 1
      self.bss.de_search()

      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters
      refine = True

    if grid_search:
      self.bss.grid_search(n_cycle=40)
      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters
    if refine:
      self.bss.scale_it()
      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters

    self.twin_fraction = self.twin_fraction_object.twin_fraction
    self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    self.free_target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    self.data_core.koverall( self.scaling_parameters.k_overall )
    self.data_core.ustar( self.scaling_parameters.u_star )
    self.data_core.ksol( self.scaling_parameters.k_sol )
    self.data_core.usol( self.scaling_parameters.u_sol )
    self.apply_back_b_iso()
    self.update_xray_structure(update_f_calc=True)


  def update_core(self,
                  f_calc        = None,
                  f_mask        = None,
                  f_part        = None,
                  b_cart        = None,
                  k_sol         = None,
                  b_sol         = None,
                  u_sol         = None,
                  k_overall     = None,
                  twin_fraction = None,
                  r_free_flags  = None):
    if f_calc is not None:
      self.data_core.renew_fatoms( f_calc.data() )
      self.f_atoms = f_calc
    else:
      assert self.f_atoms.indices().all_eq( self.miller_set.indices() )
      self.data_core.renew_fatoms( self.f_atoms.data() )

    if f_mask is not None:
      self.data_core.renew_fmask( f_mask.data() )
      self.f_mask_array = f_mask
    else:
      self.data_core.renew_fmask( self.f_mask_array.data() )

    if f_part is not None:
      self.data_core.renew_fpart( f_part.calc() )
      self.f_partial_array = f_part
    else:
      if self.f_partial_array is not None:
        self.data_core.renew_fpart( self.f_partial_array.data() )

    assert ([u_sol,b_sol]).count(None)>1

    if b_sol is not None:
      u_sol = adptbx.b_as_u( b_sol )
    if u_sol is not None:
      self.data_core.usol( u_sol )
      self.scaling_parameters.u_sol = u_sol
    if u_sol is None:
       self.data_core.usol( self.scaling_parameters.u_sol )

    if k_sol is not None:
      self.data_core.ksol( k_sol )
      self.scaling_parameters.k_sol = k_sol
    else:
      self.data_core.ksol( self.scaling_parameters.k_sol )

    if k_overall is not None:
      self.scaling_parameters.k_overall = k_overall
      self.data_core.koverall( self.scaling_parameters.k_overall )
    else:
      self.data_core.koverall( self.scaling_parameters.k_overall )

    if b_cart is not None:
      u_star = adptbx.u_cart_as_u_star( self.xs.unit_cell(), adptbx.b_as_u( list(b_cart) ) )
      self.data_core.ustar(u_star)
      self.scaling_parameters.u_star = u_star

    if twin_fraction is None:
      self.twin_fraction = self.twin_fraction_object.twin_fraction
      self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
      self.free_target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    else:
      self.twin_fraction = twin_fraction
      self.target_evaluator.alpha( twin_fraction )
      self.free_target_evaluator.alpha( twin_fraction )


  def update(self, f_calc              = None,
                   f_obs               = None,
                   f_mask              = None,
                   f_ordered_solvent   = None,
                   r_free_flags        = None,
                   b_cart              = None,
                   k_sol               = None,
                   b_sol               = None,
                   sf_and_grads_accuracy_params = None,
                   target_name         = None,
                   abcd                = None,
                   alpha_beta_params   = None,
                   xray_structure      = None,
                   mask_params         = None,
                   overall_scale       = None,
                   twin_fraction       = None ):
    if(sf_and_grads_accuracy_params is not None):
      self.sfg_params = sf_and_grads_accuracy_params
      self.update_xray_structure(update_f_calc  = True)

    if(f_calc is not None):
       assert f_calc.indices().all_eq(self.f_model.indices())
       self.update_core(f_calc = f_calc)

    if(mask_params is not None):
       self.mask_params = mask_params
    if(f_obs is not None):
       assert f_obs.data().size() == self.f_obs.data().size()
       self.f_obs = f_obs
       self.f_obs_w = self.f_obs.select(~self.free_array.data() )
       self.f_obs_f = self.f_obs.select( self.free_array.data() )
    if(f_mask is not None):
      assert f_mask.indices().all_eq( self.f_mask_array().indices() )
      assert f_mask.data().size() == self.f_mask_array().data().size()
      self.update_core(f_mask = f_mask)

    if(r_free_flags is not None):
      self.update_r_free_flags(r_free_flags)
      self.update_core(r_free_flags = r_free_flags)
    if(b_cart is not None):
      try: assert b_cart.size() == 6
      except: assert len(b_cart) == 6
      self.update_core(b_cart = b_cart)
    if overall_scale is not None:
      self.scaling_parameters.k_overall = overall_scale
      self.update_core()

    if twin_fraction is None:
      self.twin_fraction = self.twin_fraction_object.twin_fraction
      self.target

    if twin_fraction is None:
      self.twin_fraction = self.twin_fraction_object.twin_fraction
      self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
      self.free_target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    else:
      self.twin_fraction = twin_fraction
      self.target_evaluator.alpha( twin_fraction )
      self.free_target_evaluator.alpha( twin_fraction )



  def construct_miller_set(self, return_free_f_atoms_array=False, external_miller_array=None):
    completion = None
    tmp_miller = external_miller_array
    if tmp_miller is None:
      tmp_miller = self.f_obs

    completion = xray.twin_completion( tmp_miller.indices(),
                                       self.xs.space_group(),
                                       tmp_miller.anomalous_flag(),
                                       self.twin_law.as_double_array()[0:9] )

    indices = completion.twin_complete()
    miller_set = miller.set(
      crystal_symmetry = self.xs,
      indices =indices,
      anomalous_flag = tmp_miller.anomalous_flag() ).map_to_asu()

    assert miller_set.is_unique_set_under_symmetry()
    if not return_free_f_atoms_array:
      return miller_set
    else:
      free_array_for_f_atoms = completion.get_free_model_selection(
        miller_set.indices(),
        self.free_array.data() )
      return miller_set, free_array_for_f_atoms




  def compute_f_atoms(self):
    """Get f calc from the xray structure"""
    if self.miller_set is None:
      self.miller_set, self.free_flags_for_f_atoms = self.construct_miller_set(True)
    tmp = self.miller_set.structure_factors_from_scatterers(
      xray_structure = self.xray_structure,
      algorithm                    = self.sfg_params.algorithm,
      cos_sin_table                = self.sfg_params.cos_sin_table,
      grid_resolution_factor       = self.sfg_params.grid_resolution_factor,
      quality_factor               = self.sfg_params.quality_factor,
      u_base                       = self.sfg_params.u_base,
      b_base                       = self.sfg_params.b_base,
      wing_cutoff                  = self.sfg_params.wing_cutoff,
      exp_table_one_over_step_size =
                      self.sfg_params.exp_table_one_over_step_size
    )
    f_atoms = tmp.f_calc()
    return f_atoms


  def apply_back_b_iso(self):
    eps = math.pi**2*8
    unit_cell = self.xray_structure.unit_cell()
    b_min = min(self.b_sol(), self.xray_structure.min_u_cart_eigenvalue())
    if(b_min < 0):
      self.xray_structure.tidy_us(u_min = 1.e-6)
    b_iso = self.b_iso()
    b_test = b_min+b_iso
    if(b_test < 0.0): b_adj = b_iso + abs(b_test) + 0.001
    else: b_adj = b_iso
    if(abs(b_adj) <= 300.0):
      b_cart = self.b_cart()
      b_cart_new = [b_cart[0]-b_adj,b_cart[1]-b_adj,b_cart[2]-b_adj,
                    b_cart[3],      b_cart[4],      b_cart[5]]
      self.update(b_cart = b_cart_new)
      self.update(b_sol = self.k_sol_b_sol()[1] + b_adj)
      self.xray_structure.shift_us(b_shift = b_adj)
      b_min = min(self.b_sol(), self.xray_structure.min_u_cart_eigenvalue())
      assert b_min >= 0.0
      self.xray_structure.tidy_us(u_min = 1.e-6)
      self.update_xray_structure(
        xray_structure           = self.xray_structure,
        update_f_calc            = True,
        update_f_mask            = False,
        update_f_ordered_solvent = False,
        out                      = None)




  def _get_step(self, update_f_ordered_solvent = False):
    step = self.f_obs.d_min()/self.mask_params.grid_step_factor
    if(step < 0.3): step = 0.3
    step = min(0.8, step)
    if(update_f_ordered_solvent): step = 0.3
    return step

  def _update_f_mask_flag(self, xray_structure, mean_shift):
    if(self.xray_structure_mask_cache is None):
       self.xray_structure_mask_cache = xray_structure.deep_copy_scatterers()
       return True
    else:
       sites_cart_1 = self.xray_structure_mask_cache.sites_cart()
       sites_cart_2 = xray_structure.sites_cart()
       self.xray_structure_mask_cache = xray_structure.deep_copy_scatterers()
       if(sites_cart_1.size() != sites_cart_2.size()): return True
       atom_atom_distances = flex.sqrt((sites_cart_1 - sites_cart_2).dot())
       mean_shift_ = flex.mean(atom_atom_distances)
       update_f_mask = False
       if(mean_shift_ >= mean_shift):
          update_f_mask = True
       return update_f_mask

  def print_diffs(self):
    sites_cart_1 = self.xray_structure_mask_cache.sites_cart()
    sites_cart_2 = self.xray_structure.sites_cart()
    atom_atom_distances = flex.sqrt((sites_cart_1 - sites_cart_2).dot())
    mean_shift_ = flex.mean(atom_atom_distances)
    print >> self.out, "MEAN SHIFT: ", mean_shift_


  def update_xray_structure(self,
                            xray_structure           = None,
                            update_f_calc            = False,
                            update_f_mask            = False,
                            update_f_ordered_solvent = False,
                            force_update_f_mask      = True,
                            out                      = None,
                            k_sol                    = None,
                            b_sol                    = None,
                            b_cart                   = None):
    if (xray_structure is not None):
      self.xray_structure = xray_structure
    if(update_f_mask):
       if(force_update_f_mask):
          consider_mask_update = True
       else:
          consider_mask_update = self._update_f_mask_flag(
                  xray_structure = self.xray_structure,
                  mean_shift     = self.mask_params.mean_shift_for_mask_update)
    f_calc = None
    if(update_f_calc):
       timer = user_plus_sys_time()
       assert self.xray_structure is not None
       f_calc = self.miller_set.structure_factors_from_scatterers(
         xray_structure = self.xray_structure,
         algorithm                    = self.sfg_params.algorithm,
         cos_sin_table                = self.sfg_params.cos_sin_table,
         grid_resolution_factor       = self.sfg_params.grid_resolution_factor,
         quality_factor               = self.sfg_params.quality_factor,
         u_base                       = self.sfg_params.u_base,
         b_base                       = self.sfg_params.b_base,
         wing_cutoff                  = self.sfg_params.wing_cutoff,
         exp_table_one_over_step_size =
                         self.sfg_params.exp_table_one_over_step_size).f_calc()
    f_mask = None
    set_core_flag=True
    if(update_f_mask and consider_mask_update):
       bulk_solvent_mask_obj = self.bulk_solvent_mask()
       f_mask = bulk_solvent_mask_obj.structure_factors(miller_set= self.miller_set)
    if([f_calc, f_mask].count(None) == 2): set_core_flag = False
    else: set_core_flag = True
    if(f_calc is None): f_calc = self.f_calc()
    if(f_mask is None): f_mask = self.f_mask()
    if(set_core_flag):
       self.update_core(f_calc = f_calc,
                        f_mask = f_mask,
                        b_cart = b_cart,
                        k_sol  = k_sol,
                        b_sol  = b_sol)


  def bulk_solvent_mask(self):
    step = self._get_step()
    result = masks.bulk_solvent(
      xray_structure           = self.xray_structure,
      grid_step                = step,
      ignore_zero_occupancy_atoms=self.mask_params.ignore_zero_occupancy_atoms,
      solvent_radius           = self.mask_params.solvent_radius,
      shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    return result

  def update_f_mask(self):
    mask = self.bulk_solvent_mask()
    self.f_mask_array = mask.structure_factors( self.miller_set )


  def r_values(self, table=True, rows=False, d_min=None, d_max=None, again=False):
    if rows:
      table=True

    additional_selection_w = flex.bool(self.f_obs_w.data().size(), True)
    d_w = self.f_obs_w.d_spacings().data()
    if d_max is not None:
      exclude_low_w  = flex.bool(d_w<d_max)
      additional_selection_w = additional_selection_w&exclude_low_w
    if d_min is not None:
      exclude_high_w = flex.bool(d_w>d_min)
      additional_selection_w = additional_selection_w&exclude_high_w

    additional_selection_f = flex.bool(self.f_obs_f.data().size(), True)
    d_f = self.f_obs_f.d_spacings().data()
    if d_max is not None:
      exclude_low_f  = flex.bool(d_f<d_max)
      additional_selection_f = additional_selection_f&exclude_low_f
    if d_min is not None:
      exclude_high_f = flex.bool(d_f>d_min)
      additional_selection_f = additional_selection_f&exclude_high_f

    r_abs_work_f_overall = self.r_work_object.r_amplitude_abs(
      f_obs         = self.f_obs_w.data(),
      f_model       = self.data_core.f_model(),
      selection     = additional_selection_w,
      twin_fraction = self.twin_fraction_object.twin_fraction)

    r_abs_free_f_overall = self.r_free_object.r_amplitude_abs(
      self.f_obs_f.data(),
      self.data_core.f_model(),
      additional_selection_f,
      self.twin_fraction_object.twin_fraction)

    #make a sigmaa object
    tmp_sigmaa_object = self.sigmaa_object()

    if table:
      r_abs_work_f_bin = []
      r_abs_free_f_bin = []
      bin_low = []
      bin_high= []
      n_free = []
      n_work = []
      rows = []
      bins = []
      for i_bin in self.f_obs_f.binner().range_used():
        selection = flex.bool( self.f_obs_w.binner().bin_indices() == i_bin )
        #combine selection
        n_work = selection.count(True)
        tmp_work = self.r_work_object.r_amplitude_abs(
          f_obs         = self.f_obs_w.data(),
          f_model       = self.data_core.f_model(),
          selection     = selection,
          twin_fraction = self.twin_fraction_object.twin_fraction)

        mean_f_obs_w = flex.mean_default( self.f_obs_w.data().select( selection ), None )


        selection = flex.bool( self.f_obs_f.binner().bin_indices() == i_bin )
        selection = selection&additional_selection_f
        n_free = selection.count(True)
        tmp_free = self.r_free_object.r_amplitude_abs(
          f_obs         = self.f_obs_f.data(),
          f_model       = self.data_core.f_model(),
          selection     = selection,
          twin_fraction = self.twin_fraction_object.twin_fraction)

        r_abs_work_f_bin.append(tmp_work)
        r_abs_free_f_bin.append(tmp_free)
        d_max,d_min = self.f_obs_w.binner().bin_d_range( i_bin )
        d_range = self.f_obs_w.binner().bin_legend(
          i_bin=i_bin, show_bin_number=False, show_counts=False)

        alpha_w, beta_w = self.alpha_beta_w()
        alpha_f, beta_f = self.alpha_beta_f()
        n_additional_selection_w = flex.bool(alpha_w.data().size(), True)
        d_w = alpha_w.d_spacings().data()
        if d_max is not None:
          exclude_low_w  = flex.bool(d_w<d_max)
          n_additional_selection_w = n_additional_selection_w&exclude_low_w
        if d_min is not None:
          exclude_high_w = flex.bool(d_w>d_min)
          n_additional_selection_w = n_additional_selection_w&exclude_high_w

        n_additional_selection_f = flex.bool(alpha_f.data().size(), True)
        d_f = alpha_f.d_spacings().data()
        if d_max is not None:
          exclude_low_f  = flex.bool(d_f<d_max)
          n_additional_selection_f = n_additional_selection_f&exclude_low_f
        if d_min is not None:
          exclude_high_f = flex.bool(d_f>d_min)
          n_additional_selection_f = n_additional_selection_f&exclude_high_f

        alpha_f = flex.mean_default( alpha_f.select( n_additional_selection_f ).data(), None )
        beta_f  = flex.mean_default( beta_f.select( n_additional_selection_f ).data(), None )
        alpha_w = flex.mean_default( alpha_w.select( n_additional_selection_w ).data(), None )
        beta_w  = flex.mean_default( beta_w.select( n_additional_selection_w ).data(), None )

        phase_error_w = flex.mean_default( self.phase_errors_work().select( n_additional_selection_w), None )
        phase_error_f = flex.mean_default( self.phase_errors_test().select( n_additional_selection_f), None )

        fom_w = flex.mean_default( self.figures_of_merit_work().select( n_additional_selection_w ), None )

        tmp = [ str( "%3i"%(i_bin)    ),
                str( "%5.2f"%(d_max)  ),
                str( "%5.2f"%(d_min)  ),
                str( "%5i"%(n_work)   ),
                str( "%3.2f"%(tmp_work) ),
                str( "%5i"%(n_free)   ),
                str( "%3.2f"%(tmp_free) ) ]
        """
               i_bin        = None,
               d_range      = None,
               completeness = None,
               alpha_work   = None,
               beta_work    = None,
               r_work       = None,
               r_free       = None,
               target_work  = None,
               target_free  = None,
               n_work       = None,
               n_free       = None,
               mean_f_obs   = None,
               fom_work     = None,
               pher_work    = None,
               pher_free    = None
        """
        rows.append( tmp )
        bin = resolution_bin(i_bin=i_bin,
                             d_range=d_range,
                             completeness=0.98,
                             alpha_work=alpha_w,
                             beta_work=beta_w,
                             r_work=tmp_work,
                             r_free=tmp_free,
                             target_work=None,
                             target_free=None,
                             n_work=n_work,
                             n_free=n_free,
                             scale_k1_work= None, # XXX for Peter to fix.
                             mean_f_obs = mean_f_obs_w,
                             fom_work = fom_w,
                             pher_work = phase_error_w,
                             pher_free = phase_error_f )

        bins.append( bin )

      if not rows:
        header = ("bin", "d_max", "d_min", "n_work", "r_work", "n_free", "r_free")
        comments = """
Overall r values
R Work : %4.3f
R Free : %4.3f

R  = \sum_h( |Ft-Fo| )/ \sum_h(Fo)
Ft = Sqrt(tf*F1^2 + (1-tf)F2^2)
F1,F2 are twin related model amplitudes.
tf is the twin fraction and Fo is an observed amplitude."""%(r_abs_work_f_overall, r_abs_free_f_overall)

        table_txt = table_utils.format( [header]+rows,
                                        comments=comments,
                                        has_header=True,
                                        separate_rows=False,
                                        prefix='| ',
                                        postfix=' |')
        print >> self.out, "------------------------  R values ------------------------"
        print >> self.out, "  twin law      : %s"%( sgtbx.change_of_basis_op( self.twin_law ).as_hkl() )
        print >> self.out, "  twin fraction : %4.3f"%( self.twin_fraction_object.twin_fraction)
        print >> self.out, table_txt
        print >> self.out, "-----------------------------------------------------------"
        print >> self.out
        self.r_work_in_lowest_resolution_bin(show=True)
        self.r_overall_low_high(show=True)
      else:
        return bins

    else:
      return r_abs_work_f_overall, r_abs_free_f_overall

  def r_work(self):
    w,f = self.r_values(False)
    return w

  def r_free(self):
    w,f = self.r_values(False)
    return f

  def r_all(self):
    selection = flex.bool( self.f_obs.data().size(), True )
    overall_r = self.r_all_object.r_amplitude_abs(
          f_obs         = self.f_obs.data(),
          f_model       = self.data_core.f_model(),
          selection     = selection,
          twin_fraction = self.twin_fraction_object.twin_fraction)
    return overall_r


  def r_work_in_lowest_resolution_bin(self, reflections_per_bin=200, show=False):
    d_star_sq = self.f_obs_w.d_star_sq().data()
    sort_permut = flex.sort_permutation( d_star_sq )
    if sort_permut.size() < reflections_per_bin:
      reflections_per_bin = sort_permut.size()
    i_select = sort_permut[:reflections_per_bin-1]
    b_select = flex.bool(sort_permut.size(), False )
    b_select = b_select.set_selected( i_select, True )

    tmp_work = self.r_work_object.r_amplitude_abs(
      f_obs         = self.f_obs_w.data(),
      f_model       = self.data_core.f_model(),
      selection     = b_select,
      twin_fraction = self.twin_fraction_object.twin_fraction)

    if not show:
      return tmp_work
    else:
      print >> self.out, "-----------------------------------------------------------"
      print >> self.out, "  R-value for the %i lowest resolution reflections:"%(reflections_per_bin)
      print >> self.out, "    %4.3f" %(self.r_work_in_lowest_resolution_bin(reflections_per_bin))
      print >> self.out, "-----------------------------------------------------------"


  def r_overall_low_high(self, d = 6.0, show=False):
    r_work = self.r_work()
    d_max, d_min = self.f_obs_w.d_max_min()
    if(d_max < d): d = d_max
    if(d_min > d): d = d_min

    n_low = self.f_obs_w.resolution_filter(d_min = d, d_max = 999.9).data().size()
    if(n_low > 0):
       r_work_l = self.r_values(d_min = d, d_max = 999.9, table=False )[0]
    else:
       r_work_l = None


    n_high = self.f_obs_w.resolution_filter(d_min = 0.0, d_max = d).data().size()
    if(n_high > 0):
       r_work_h = self.r_values(d_min = 0.0, d_max = d,table=False)[0]
    else:
       r_work_h = None


    if(r_work_l is not None):
       r_work_l = r_work_l
    else:
       r_work_l = 0.0

    if(r_work_h is not None):
       r_work_h = r_work_h
    else:
       r_work_h = 0.0
    if not show:
      return r_work, r_work_l, r_work_h, n_low, n_high
    else:
      print >> self.out, "----------------------------------------------------------"
      print >> self.out, "Overall, low and high resolution R-work values"
      print >> self.out
      print >> self.out, "Limits: Overall: %6.2f -- %6.2f"%(d_max,d_min)
      print >> self.out, "        Low    : %6.2f -- %6.2f"%(d_max,d)
      print >> self.out, "        High   : %6.2f -- %6.2f"%(d,d_min)
      print >> self.out
      print >> self.out, "R values    : Overall    low    high"
      print >> self.out, "              %6.3f   %6.3f  %6.3f"%(r_work,r_work_l,r_work_h)
      print >> self.out, "Contributors:%7i  %7i %7i"%(n_low+n_high, n_low,n_high)
      print >> self.out, "----------------------------------------------------------"



  def twin_fraction_scan(self, n=10):
    """for each twin fraction, compute the target value and r value"""
    print >> self.out
    print >> self.out
    print >> self.out, "------------------------ Twin fraction scan ----------------------"
    print >> self.out
    print >> self.out, " R-values and target values for various twin fractions are listed."
    print >> self.out
    current_twin_fraction = twin_fraction_object(self.twin_fraction_object.twin_fraction)
    trail_twin_fractions = list( flex.double( range(n+1) )/(2.0*n) )
    rows = []
    for tf in trail_twin_fractions:
      tmp_twin_fraction = twin_fraction_object( tf )
      self.update_solvent_and_scale( twin_fraction_parameters =  tmp_twin_fraction )
      rw,rf = self.r_values(table=False)
      ttw,ttf = self.target(print_it=False)
      tmp = [ "%4.3f"%(tf),
              "%4.3f"%(rw),
              "%4.3f"%(rf),
              "%5.4e"%(ttw),
              "%5.4e"%(ttf)
              ]
      rows.append( tmp )

    legend = ( "Twin fraction", "R-work", "R-free", "Target-work", "Target-free" )
    table_txt = table_utils.format( [legend]+rows,
                                    comments=None,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')
    print >> self.out, table_txt
    print >> self.out
    print >> self.out,  "------------------------------------------------------------------"
    print >> self.out
    print >> self.out
    self.update_solvent_and_scale( twin_fraction_parameters =  current_twin_fraction )



  def target(self, print_it=True):
    tmp_w=self.target_evaluator.target( self.data_core.f_model() )/self.norma_sum_f_sq_w
    tmp_f=self.free_target_evaluator.target( self.data_core.f_model() )/self.norma_sum_f_sq_f
    if print_it:
      print >> self.out
      print >> self.out, "----------------- Target values -----------------"
      print >> self.out, "   working set  : %8.6e "%(tmp_w)
      print >> self.out, "   free set     : %8.6e "%(tmp_f)
      print >> self.out, "-------------------------------------------------"
    else:
      return(tmp_w,tmp_f)

  def target_functor(self):
    return target_functor(manager=self)

  def target_f(self):
    return self.target_t()

  def detwin_data(self, mode=None):
    #determine how to detwin
    if mode is None:
      mode = "auto"
    if mode == "auto":
      if self.twin_fraction_object.twin_fraction > self.detwin_switch_twin_fraction:
        mode = "proportional"
        if self.twin_fraction_object.twin_fraction > 1.0 - self.detwin_switch_twin_fraction:
          mode = "algebraic"
      else:
        mode = "algebraic"

    if mode == "algebraic":
      if  abs(self.twin_fraction_object.twin_fraction-0.5)<1e-3:
        print >> self.out, "Automatic adjustment: detwinning mode set to proportional"

    assert mode in self.possible_detwin_modes
    assert mode != "auto"

    #detwinning should be done against
    tmp_i_obs = self.f_obs.deep_copy().f_as_f_sq()
    untouched = self.f_obs.deep_copy().f_as_f_sq()
    dt_f_obs = None
    tmp_free = self.free_array.deep_copy()

    # now please detwin the data
    if mode == "proportional":
      if (tmp_i_obs.sigmas() is None):
        raise Sorry("Detwinning requires experimental sigmas.")
      dt_iobs, dt_isigma = self.full_detwinner.detwin_with_model_data(
        tmp_i_obs.data(),
        tmp_i_obs.sigmas(),
        self.data_core.f_model(),
        self.twin_fraction_object.twin_fraction )
      tmp_i_obs = tmp_i_obs.customized_copy(
        data = dt_iobs, sigmas = dt_isigma).set_observation_type( tmp_i_obs )
      dt_f_obs = tmp_i_obs.f_sq_as_f()

    if mode == "algebraic":
      dt_iobs, dt_isigma = self.full_detwinner.detwin_with_twin_fraction(
        tmp_i_obs.data(),
        tmp_i_obs.sigmas(),
        self.twin_fraction_object.twin_fraction )
      # find out which intensities are zero or negative, they will be eliminated later on
      zeros = flex.bool( dt_iobs <= 0 )
      x = dt_iobs.select( ~zeros )
      y = tmp_i_obs.data().select( ~zeros )

      tmp_i_obs = tmp_i_obs.customized_copy(
        data = dt_iobs, sigmas = dt_isigma).set_observation_type( tmp_i_obs )
      dt_f_obs = tmp_i_obs.select( ~zeros ).f_sq_as_f()
      tmp_free = tmp_free.select( ~zeros )
      untouched = untouched.select( ~zeros )

    #we can now quickly scale the two and see what hapens.
    abs_tmp_f_model = abs( self.f_model() ).common_set( dt_f_obs ).set_observation_type( dt_f_obs )
    tmp_f_model = self.f_model().common_set( dt_f_obs )
    scaler = relative_scaling.ls_rel_scale_driver(
      miller_native = abs_tmp_f_model,
      miller_derivative = dt_f_obs,
      use_intensities = False,
      scale_weight = False,
      use_weights = False )

    if dt_f_obs.data().size()==0:
      if mode == "algebraic":
        raise Sorry("Algebraic detwinning of data resulted in a dataset without data! \n Please try to use detwinning using proportionality rules (detwin.mode=proportional) ")
      else:
        raise Sorry("This should never have happend. Please contact authors")
    return dt_f_obs, tmp_f_model, tmp_free

  def sigmaa_object(self, detwinned_data=None, f_model_data=None, forced_update=False):
    if self.sigmaa_object_cache is None:
      forced_update = True
    assert ( [detwinned_data,f_model_data] ).count(None) != 1
    tmp_free = self.free_array
    if (detwinned_data is None):
      if forced_update or self.update_sigmaa_object is True:
        self.update_sigmaa_object = True
        detwinned_data,f_model_data,tmp_free = self.detwin_data(mode=self.detwin_mode)
    if forced_update:
      self.update_sigmaa_object = True
      detwinned_data,f_model_data,tmp_free = self.detwin_data(mode=self.detwin_mode)
    if self.update_sigmaa_object:
      self.update_sigmaa_object = False
      self.sigmaa_object_cache = sigmaa_estimation.sigmaa_estimator(
        miller_obs   = detwinned_data,
        miller_calc  = f_model_data,
        r_free_flags = tmp_free,
        kernel_width_free_reflections=200,
        )
    return self.sigmaa_object_cache

  def model_error_ml(self):
    return None

  def alpha_beta(self, external_sigmaa_object=None):
    sigmaa_object = self.sigmaa_object()
    return sigmaa_object.alpha_beta()

  def alpha_beta_w(self, only_if_required_by_target=False):
    a,b = self.alpha_beta()
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    a = a.select( tmp_free.data() )
    b = b.select( tmp_free.data() )
    return a,b

  def alpha_beta_f(self,only_if_required_by_target=False):
    a,b = self.alpha_beta()
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    a = a.select( ~tmp_free.data() )
    b = b.select( ~tmp_free.data() )
    return a,b

  def figures_of_merit(self):
    sigmaa_object = self.sigmaa_object()
    return sigmaa_object.fom().data()

  def figures_of_merit_work(self):
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    fom = self.figures_of_merit().select(
      tmp_free.data())
    return fom

  def figures_of_merit_t(self):
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    fom = self.figures_of_merit().select(
      ~tmp_free.data())
    return fom

  def phase_errors(self):
    sigmaa_object = self.sigmaa_object()
    return sigmaa_object.phase_errors().data()

  def phase_errors_work(self):
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    pher = self.phase_errors().select(tmp_free.data())
    return pher

  def phase_errors_test(self):
    tmp_sigmaa = self.sigmaa_object()
    tmp_free = tmp_sigmaa.r_free_flags
    pher = self.phase_errors().select(~tmp_free.data())
    return pher

  def determine_n_bins(self,
        free_reflections_per_bin,
        max_n_bins=None,
        min_n_bins=1,
        min_refl_per_bin=100):
    assert free_reflections_per_bin > 0
    n_refl = self.free_array.data().size()
    n_free = self.free_array.data().count(True)
    n_refl_per_bin = free_reflections_per_bin
    if (n_free != 0):
      n_refl_per_bin *= n_refl / n_free
    n_refl_per_bin = min(n_refl, iround(n_refl_per_bin))
    result = max(1, iround(n_refl / max(1, n_refl_per_bin)))
    if (min_n_bins is not None):
      result = max(result, min(min_n_bins, iround(n_refl / min_refl_per_bin)))
    if (max_n_bins is not None):
      result = min(max_n_bins, result)
    return result



  def _map_coeff(self, f_obs, f_model, f_obs_scale, f_model_scale):
    d_obs = miller.array(miller_set = f_model,
                         data       = f_obs.data()*f_obs_scale
                        ).phase_transfer(phase_source = f_model)
    return miller.array(miller_set = f_model,
                        data       = d_obs.data()-f_model.data()*f_model_scale)



  def map_coefficients(self,
                       map_type = None,
                       k        = None,
                       n        = None,
                       w1       = None,
                       w2       = None
                       ):
    assert map_type in ("Fo-Fc", "Fobs-Fmodel",
                        "2mFo-DFc", "2mFobs-DFmodel",
                        "mFo-DFc", "mFobs-DFmodel",
                        "gradient",
                        "m_gradient"
                        )
    # this is to modify default behavoir of phenix.refine
    if (map_type == "mFo-DFc") or (map_type == "mFobs-DFmodel") :
      if self.map_types.fofc == "gradient":
        map_type = "gradient"
      if self.map_types.fofc == "m_gradient":
        map_type = "m_gradient"
      if self.map_types.fofc == "m_dtfo_d_fc":
        map_type = "m_dtfo_d_fc"

    if (map_type == "2mFo-DFc") or (map_type=="2mFobs-DFmodel") :
      if self.map_types.twofofc == "two_m_dtfo_d_fc":
        map_type = "two_m_dtfo_d_fc"
      if  self.map_types.twofofc == "two_dtfo_fc":
        map_type = "two_dtfo_fc"

    #detwin
    dt_f_obs, tmp_f_model, tmp_free = self.detwin_data(mode=self.detwin_mode)
    #for aniso correction
    aniso_scale = 1.0/self.data_core.overall_scale() # anisotropy correction
    aniso_scale = self.f_atoms.customized_copy(
      data = aniso_scale ).common_set( dt_f_obs )
    aniso_scale = aniso_scale.data()

    if map_type not in ["gradient","m_gradient"]:
      result = None
      if map_type == "Fo-Fc":
        if ([k,n]).count(None) > 0:
          raise Sorry("Map coefficient multipliers (k and n) must be provided to generate detwinned maps")
        result = self._map_coeff( f_obs         = dt_f_obs,
                                  f_model       = tmp_f_model,
                                  f_obs_scale   = k,
                                  f_model_scale = n )

        assert result is not None
      else:
        sigmaa_object = self.sigmaa_object()
        dt_f_obs, tmp_f_model = dt_f_obs.common_sets( tmp_f_model )
        m, dt_f_obs = sigmaa_object.fom().common_sets( dt_f_obs)
        d, dt_f_obs = sigmaa_object.alpha_beta()[0].common_sets( dt_f_obs )
        m = m.data()
        d = d.data()
        dt_f_obs, tmp_f_model = dt_f_obs.common_sets( tmp_f_model )

        if map_type == "m_dtfo_d_fc":
          result = self._map_coeff( f_obs   = dt_f_obs,
                                    f_model = tmp_f_model,
                                    f_obs_scale  = m ,
                                    f_model_scale =d )
        if map_type == "dtfo_fc":
          result = self._map_coeff( f_obs   = dt_f_obs,
                                    f_model = tmp_f_model,
                                    f_obs_scale  = 1.0,
                                    f_model_scale =1.0 )
        if map_type == "two_m_dtfo_d_fc":
          result = self._map_coeff( f_obs   = dt_f_obs,
                                    f_model = tmp_f_model,
                                    f_obs_scale  = 2*m,
                                    f_model_scale = d )
        if map_type == "two_dtfo_fc":
          result = self._map_coeff( f_obs   = dt_f_obs,
                                    f_model = tmp_f_model,
                                    f_obs_scale  = 2,
                                    f_model_scale = 1 )

        assert result is not None
      assert result != None

      if self.map_types.aniso_correct:
        result = result.customized_copy( data = result.data()*aniso_scale )

      return result

    else:
      # get coefficients for a gradient map please
      gradients = self.target_evaluator.d_target_d_fmodel(
        self.data_core.f_model() )
      gradients = self.f_atoms.customized_copy(
        data = -gradients).common_set( self.f_obs )
      if map_type == "m_gradient":
        # get the FOMs please
        m = self.sigmaa_object().fom().common_set(self.f_obs).data()
        gradients = gradients.customized_copy(
          data = gradients.data()*m )

      if self.map_types.aniso_correct:
        gradients = gradients.customized_copy( data = gradients.data()*aniso_scale )
      return gradients



  def electron_density_map(self,
                           map_type          = "Fo-Fc",
                           k                 = 1,
                           n                 = 1,
                           w1                = None,
                           w2                = None,
                           resolution_factor = 1/3.,
                           symmetry_flags = None):
    assert map_type in ("Fo-Fc", "Fobs-Fmodel"
                        "2mFo-DFc", "2mFobs-DFmodel",
                        "mFo-DFc", "mFobs-DFmodel",
                        "gradient",
                        "m_gradient")

    return self.map_coefficients(
      map_type          = map_type,
      k                 = k,
      n                 = n,
      w1                = w1,
      w2                = w2).fft_map(
         resolution_factor = resolution_factor,
         symmetry_flags    = symmetry_flags)

  def u_star(self):
    return self.data_core.ustar()

  def u_cart(self):
    tmp = self.u_star()
    tmp = adptbx.u_star_as_u_cart(self.xs.unit_cell(),tmp)
    return tmp

  def b_cart(self):
    b_cart = adptbx.u_as_b( self.u_cart() )
    return b_cart

  def b_iso(self):
    b_cart = self.b_cart()
    return (b_cart[0]+b_cart[1]+b_cart[2])/3.0

  def u_iso(self):
    u_cart = self.u_cart()
    return (u_cart[0]+u_cart[1]+u_cart[2])/3.0

  def u_iso_as_u_cart(self):
    ui = self.u_iso()
    return [ui,ui,ui,0.0,0.0,0.0]

  def k_sol(self):
    return self.data_core.ksol()

  def u_sol(self):
    return self.data_core.usol()

  def b_sol(self):
    return adptbx.u_as_b( self.u_sol() )

  def k_sol_b_sol(self):
    return self.k_sol(), self.b_sol()

  def k_sol_u_sol(self):
    return self.k_sol(), self.u_sol()


  def f_mask(self):
    return self.f_mask_array

  def f_mask_w(self):
    return self.f_mask().select(~self.free_flags_for_f_atoms )

  def f_mask_t(self):
    return self.f_mask().select( self.free_flags_for_f_atoms )

  def f_bulk(self):
    tmp = self.data_core.f_bulk()
    tmp = self.f_mask_array.customized_copy(
      data = tmp ).set_observation_type( self.f_mask_array )
    return tmp

  def f_bulk_t(self):
    tmp = self.f_bulk()
    return tmp.select( self.free_flags_for_f_atoms )

  def f_bulk_w(self):
    tmp = self.f_bulk()
    return tmp.select(~self.free_flags_for_f_atoms )

  def fb_bulk(self):
    tmp = self.data_core.f_bulk()
    multi = self.data_core.overall_scale()
    tmp = self.f_mask_array.customized_copy(
      data = tmp*multi ).set_observation_type( self.f_mask_array )
    return tmp

  def fb_bulk_t(self):
    tmp = self.f_bulk()
    return tmp.select( self.free_flags_for_f_atoms )

  def fb_bulk_w(self):
    tmp = self.f_bulk()
    return tmp.select(~self.free_flags_for_f_atoms )

  def scale_k1_w(self):
    return self.data_core.koverall()
  def scale_k1_t(self):
    return self.data_core.koverall()
  def scale_k3_t(self):
    return self.data_core.koverall()
  def scale_k3_w(self):
    return self.data_core.koverall()




  def fft_vs_direct(self, reflections_per_bin = 250,
                          n_bins              = 0,
                          out                 = None):
    print >> self.out, "Direct vs FFT comparison not yet implemented. "

  def r_work_scale_k1_completeness_in_bins(self, reflections_per_bin = 500,
                                                 n_bins              = 0,
                                                 prefix              = "",
                                                 out                 = None):
    self.r_values(table=True)


  def show_k_sol_b_sol_b_cart_target(self, header=None,target=None,out=None):
    if(out is None): out = self.out
    p = " "
    if(header is None): header = ""
    line_len = len("|-"+"|"+header)
    fill_len = 80-line_len-1
    print >> out, "|-"+header+"-"*(fill_len)+"|"
    k_sol = self.k_sol()
    b_sol = self.b_sol()
    u0,u1,u2,u3,u4,u5 = self.b_cart()

    target_w=self.target_w()

    alpha, beta = self.alpha_beta_w()
    alpha_d = alpha.data()
    a_mean = flex.mean(alpha_d)
    a_zero = (alpha_d <= 0.0).count(True)
    r_work = self.r_work()
    u_isos = self.xray_structure.extract_u_iso_or_u_equiv()
    b_iso_mean = flex.mean(u_isos * math.pi**2*8)
    print >> out, "| k_sol=%5.2f b_sol=%7.2f target_w =%20.6f r_work=%7.4f" % \
                  (k_sol, b_sol, target_w, r_work) + 5*p+"|"
    print >> out, "| B(11,22,33,12,13,23)=%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f |" % \
                  (u0,u1,u2,u3,u4,u5)
    print >> out, "| trace(B) = (B11 + B22 + B33)/3 = %-10.3f                                 |"%self.u_iso()
    print >> out, "| mean alpha:%8.4f  number of alpha <= 0.0:%7d" % \
                  (a_mean, a_zero)+25*p+"|"
    print >> out, "|"+"-"*77+"|"
    out.flush()


  def show_essential(self, header = None, out=None):
    if(out is None): out = self.out
    out.flush()
    p = " "
    if(header is None): header = ""
    d_max, d_min = self.f_obs.d_max_min()
    line1 = "---(resolution: "
    line2 = n_as_s("%6.2f",d_min)
    line3 = n_as_s("%6.2f",d_max)
    line4 = " - "
    line5 = " A)"
    tl = header+line1+line2+line4+line3+line5
    line_len = len("|-"+"|"+tl)
    fill_len = 80-line_len-1
    print >> out, "|-"+tl+"-"*(fill_len)+"|"
    print >> out, "| "+"  "*38+"|"
    r_work = n_as_s("%6.4f",self.r_work()    )
    r_free = n_as_s("%6.4f",self.r_free()    )
    scale  = n_as_s("%6.3f",self.scale_k1_w())
    k_sol  = n_as_s("%4.2f",self.k_sol())
    b_sol  = n_as_s("%6.2f",self.b_sol())
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.2f",self.b_cart())
    b_iso  = n_as_s("%7.2f",self.b_iso())

    #XXXX Model error analyses required
    #err    = n_as_s("%6.2f",self.model_error_ml())
    err=" None "
    try:    target_work = n_as_s("%.4g",self.target_w())
    except: target_work = str(None)

    line = "| r_work= "+r_work+"   r_free= "+r_free+"   ksol= "+k_sol+\
           "   Bsol= "+b_sol+"   scale= "+scale
    np = 79 - (len(line) + 1)
    if(np < 0): np = 0
    print >> out, line + p*np + "|"
    print >> out, "| "+"  "*38+"|"
    print >> out, "| overall anisotropic scale matrix (Cartesian basis):    "\
                  "                     |"
    c = ","
    line4 = "| (B11,B22,B33,B12,B13,B23)= ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+")"
    np = 79 - (len(line4) + 1)
    line4 = line4 + " "*np + "|"
    print >> out, line4
    line5 = "| (B11+B22+B33)/3 = "+b_iso
    np = 79 - (len(line5) + 1)
    line5 = line5 + " "*np + "|"
    print >> out, line5
    print >> out, "| "+"  "*38+"|"
    line5_and_a_half = "| Twin law : %s     Twin fraction: %4.3f"%(self.twin_law.r().as_hkl(),self.twin_fraction_object.twin_fraction)
    np = 79 - (len(line5_and_a_half) + 1)
    line5_and_a_half = line5_and_a_half + " "*np + "|"
    print >> out, line5_and_a_half
    print >> out, "| "+"  "*38+"|"
    line6="| Target ("+self.target_name+")= "+target_work+\
          " | ML estimate for coordinates error: "+err+" A"
    np = 79 - (len(line6) + 1)
    line6 = line6 + " "*np + "|"
    print >> out, line6
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_comprehensive(self,
                         header = "",
                         free_reflections_per_bin = 140,
                         max_number_of_bins  = 30,
                         out=None):
    self.r_values(table=True)
    self.twin_fraction_scan(n=15)
    self.sigmaa_object().show(out=self.out)




  def statistics_in_resolution_bins(self,
                                    free_reflections_per_bin = 200,
                                    max_number_of_bins  = 30,
                                    out=None):
    results = self.r_values(table=True, rows=True)
    return results
    # first detemine number of bins please
    # self.twin_fraction_scan(n=15)
    # self.sigmaa_object().show(out=self.out)

  def r_factors_in_resolution_bins(self,
                                   free_reflections_per_bin = 200,
                                   max_number_of_bins  = 30,
                                   out=None):
    results = self.r_values(table=True, rows=True)
    return results
    #self.twin_fraction_scan(n=15)


  def r_work_scale_k1_completeness_in_bins(self, reflections_per_bin = 500,
                                                 n_bins              = 0,
                                                 prefix              = "",
                                                 out                 = None):
    #actively ignoring input
    self.r_values(table=True)



  def show_fom_phase_error_alpha_beta_in_bins(self,
                                              free_reflections_per_bin = 200,
                                              max_number_of_bins  = 30,
                                              out=None):
    self.sigmaa_object().show(out=self.out)


  def export(self, out=None, format="mtz"):
    assert format in ["mtz", "cns"]
    file_name = None
    if (out is None):
      out = sys.stdout
    elif (hasattr(out, "name")):
      file_name = libtbx.path.canonical_path(file_name=out.name)
    warning = [
      "DO NOT USE THIS FILE AS INPUT FOR REFINEMENT!",
      "Resolution and sigma cutoffs may have been applied to FOBS."]
    width = max([len(line) for line in warning])
    warning.insert(0, "*" * width)
    warning.append(warning[0])


    if (format == "cns"):
      for line in warning:
        print >> out, "{ %s%s }" % (line, " "*(width-len(line)))
      print >> out, "{ %s }" % date_and_time()
      if (file_name is not None):
        print >> out, "{ file name: %s }" % os.path.basename(file_name)
        print >> out, "{ directory: %s }" % os.path.dirname(file_name)
      self.explain_members(out=out, prefix="{ ", suffix=" }")
      crystal_symmetry_as_cns_comments(
        crystal_symmetry=self.f_obs, out=out)
      print >> out, "NREFlection=%d" % self.f_obs.indices().size()
      print >> out, "ANOMalous=%s" % {0: "FALSE"}.get(
        int(self.f_obs.anomalous_flag()), "TRUE")
      have_sigmas = self.f_obs.sigmas() is not None
      for n_t in [("FOBS", "REAL"),
                  ("SIGFOBS", "REAL"),
                  ("R_FREE_FLAGS", "INTEGER"),
                  ("FMODEL", "COMPLEX"),
                  ("FCALC", "COMPLEX"),
                  ("FMASK", "COMPLEX"),
                  ("FBULK", "COMPLEX"),
                  ("FOM", "REAL"),
                  ("ALPHA", "REAL"),
                  ("BETA", "REAL")]:
        if (not have_sigmas and n_t[0] == "SIGFOBS"): continue
        print >> out, "DECLare NAME=%s DOMAin=RECIprocal TYPE=%s END" % n_t
      f_model            = self.f_model_scaled_with_k1()
      f_model_amplitudes = f_model.amplitudes().data()
      f_model_phases     = f_model.phases(deg=True).data()
      f_calc_amplitudes  = self.f_calc().amplitudes().data()
      f_calc_phases      = self.f_calc().phases(deg=True).data()
      f_mask_amplitudes  = self.f_mask().amplitudes().data()
      f_mask_phases      = self.f_mask().phases(deg=True).data()
      f_bulk_amplitudes  = self.f_bulk().amplitudes().data()
      f_bulk_phases      = self.f_bulk().phases(deg=True).data()
      alpha, beta        = [item.data() for item in self.alpha_beta()]
      arrays = [
        self.f_obs.indices(), self.f_obs.data(), self.f_obs.sigmas(),
        self.free_array.data(),
        f_model_amplitudes, f_model_phases,
        f_calc_amplitudes, f_calc_phases,
        f_mask_amplitudes, f_mask_phases,
        f_bulk_amplitudes, f_bulk_phases,
        self.figures_of_merit(),
        alpha, beta]
      if (not have_sigmas):
        del arrays[2]
        i_r_free_flags = 2
      else:
        i_r_free_flags = 3
      for values in zip(*arrays):
        print >> out, "INDE %d %d %d" % values[0]
        print >> out, " FOBS= %.6g" % values[1],
        if (have_sigmas):
          print >> out, " SIGFOBS= %.6g" % values[2],
        print >> out, \
          " R_FREE_FLAGS= %d FMODEL= %.6g %.6g\n" \
          " FCALC= %.6g %.6g FMASK= %.6g %.6g FBULK= %.6g %.6g\n" \
          " FB_CART= %.6g FOM= %.6g ALPHA= %.6g BETA= %.6g"  \
            % values[i_r_free_flags:]
      if (file_name is not None):
        out.close()

    else:
      assert file_name is not None
      mtz_dataset = self.f_obs.as_mtz_dataset(column_root_label="FOBS")
      mtz_dataset.add_miller_array(
        miller_array=self.free_array, column_root_label="R_FREE_FLAGS")
      mtz_dataset.add_miller_array(
        miller_array=self.f_model(), column_root_label="FMODEL")
      mtz_dataset.add_miller_array(
        miller_array=self.f_calc(), column_root_label="FCALC")
      mtz_dataset.add_miller_array(
        miller_array=self.f_mask(), column_root_label="FMASK")
      mtz_dataset.add_miller_array(
        miller_array=self.f_bulk(), column_root_label="FBULK")
      mtz_dataset.add_miller_array(
        miller_array= self.sigmaa_object().fom(),
        column_root_label="FOM", column_types="W")
      alpha, beta = self.alpha_beta()
      mtz_dataset.add_miller_array(
        miller_array=alpha, column_root_label="ALPHA", column_types="W")
      mtz_dataset.add_miller_array(
        miller_array=beta, column_root_label="BETA", column_types="W")
      mtz_history_buffer = flex.std_string(warning)
      ha = mtz_history_buffer.append
      ha(date_and_time())
      ha("file name: %s" % os.path.basename(file_name))
      ha("directory: %s" % os.path.dirname(file_name))
      s = StringIO()
      self.explain_members(out=s)
      for line in s.getvalue().splitlines():
        ha(line)
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.add_history(lines=mtz_history_buffer)
      out.close()
      mtz_object.write(file_name=file_name)



  def explain_members(self, out=None, prefix="", suffix=""):
    if (out is None): out = sys.stdout
    def zero_if_almost_zero(v, eps=1.e-6):
      if (abs(v) < eps): return 0
      return v
    for line in [
          "Fmodel   = scale_k1 * fb_cart * (Fcalc + Fbulk)",
          "Fcalc    = structure factors calculated from atomic model",
          "Fbulk    = k_sol * exp(-b_sol*s**2/4) * Fmask",
          "A        = orthogonalization matrix",
          "k_sol    = %.6g" % self.k_sol(),
          "b_sol    = %.6g" % zero_if_almost_zero(self.b_sol()),
          "B_cart = (B11, B22, B33, B12, B13, B23)",
          "       = (%s)" % ", ".join(
            ["%.6g" % zero_if_almost_zero(v) for v in self.b_cart()])]:
      print >> out, prefix + line + suffix



  def show_targets(self, out=None, text=""):
    if(out is None): out = self.out
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print >> out, part1 + "-"*n + part2
    part3 = "| target_work(%s"%self.target_name+") = %.6e  r_work = %6.4f  r_free = %6.4f"%\
                                (self.target_w(), self.r_work(), self.r_free())
    n = 78 - len(str(part3)+"|")
    print >> out, part3, " "*n +"|"
    print >> out, "|" +"-"*77+"|"
    out.flush()

  def _header_resolutions_nreflections(self, header, out):
    out.flush()
    if(header is None): header = ""
    d_max, d_min = self.f_obs.d_max_min()
    line1 = "(resolution: "
    line2 = n_as_s("%6.2f",d_min)
    line3 = n_as_s("%6.2f",d_max)
    line4 = " - "
    line5 = " A; n_refl. = "
    line6 = n_as_s("%d",self.f_obs.data().size())
    tl = header+"-"+line1+line2+line4+line3+line5+line6
    line_len = len("|-"+"|"+tl)
    fill_len = 80-line_len-1
    print >> out, "|-"+tl+"-"*(fill_len)+"|"
    out.flush()

  def _rfactors_and_bulk_solvent_and_scale_params(self, out):
    out.flush()
    r_work = n_as_s("%6.4f",self.r_work()    )
    r_free = n_as_s("%6.4f",self.r_free()    )
    scale  = n_as_s("%6.3f",self.scale_k1_w())
    k_sol  = n_as_s("%4.2f",self.k_sol())
    b_sol  = n_as_s("%6.2f",self.b_sol())
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.2f",self.b_cart())
    b_iso  = n_as_s("%7.2f",self.b_iso())
    line = "| r_work= "+r_work+"   r_free= "+r_free+"   ksol= "+k_sol+\
           "   Bsol= "+b_sol+"   scale= "+scale
    np = 79 - (len(line) + 1)
    if(np < 0): np = 0
    print >> out, line + " "*np + "|"
    print >> out, "| "+"  "*38+"|"
    print >> out, "| overall anisotropic scale matrix (Cartesian basis; B11,B22,B33,B12,B13,B23):|"
    c = ","
    line4 = "| ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+"); trace/3= "+b_iso
    np = 79 - (len(line4) + 1)
    line4 = line4 + " "*np + "|"
    print >> out, line4
    out.flush()


def ls_ff_weights(f_obs, atom, B):
  d_star_sq_data = f_obs.d_star_sq().data()
  table = wk1995(atom).fetch()
  ff = table.at_d_star_sq(d_star_sq_data) * flex.exp(-B/4.0*d_star_sq_data)
  weights = 1.0/flex.pow2(ff)
  return weights

class target_functor(object):

  def __init__(self, manager):
    self.manager = manager

  def __call__(self, compute_gradients=False):
    return target_result(manager=self.manager)

class target_result(mmtbx.f_model.target_result_mixin):

  def __init__(self, manager):
    self.manager = manager

  def target_work(self):
    return self.manager.target(False)[0]

  def target_test(self):
    return self.manager.target(False)[1]

  def d_target_d_f_model_work(self):
    manager = self.manager
    return manager.miller_set.array(
      data=manager.target_evaluator.d_target_d_fmodel(
             manager.data_core.f_model()))

  def d_target_d_f_calc_work(self):
    manager = self.manager
    d_t_d_f_m = self.d_target_d_f_model_work()
    return d_t_d_f_m.array(
      data=d_t_d_f_m.data() * manager.data_core.d_f_model_core_data_d_f_atoms()
                            / manager.norma_sum_f_sq)

def ls_sigma_weights(f_obs):
  if(f_obs.sigmas() is not None):
     sigmas_squared = flex.pow2(f_obs.sigmas())
  else:
     sigmas_squared = flex.double(f_obs.data().size(), 1.0)
  assert sigmas_squared.all_gt(0)
  weights = 1 / sigmas_squared
  return weights

def kb_range(x_max, x_min, step):
  x_range = []
  x = x_min
  while x <= x_max + 0.0001:
    x_range.append(x)
    x += step
  return x_range

def n_as_s(format, value):
  vt = type(value).__name__
  if(vt in ["int","float"]):
     return str(format%value).strip()
  else:
     new = []
     for item in value:
       new.append( str(format%item).strip() )
  return new

class resolution_bin(object):
  def __init__(self,
               i_bin        = None,
               d_range      = None,
               completeness = None,
               alpha_work   = None,
               beta_work    = None,
               r_work       = None,
               r_free       = None,
               target_work  = None,
               target_free  = None,
               n_work       = None,
               n_free       = None,
               mean_f_obs   = None,
               fom_work     = None,
               scale_k1_work= None,
               pher_work    = None,
               pher_free    = None):
    adopt_init_args(self, locals())


def n_as_s(format, value):
  if value == "none":
    return "None"
  if value == "None":
    return "None"
  if ( value is None ):
    return format_value(format=format, value=value)
  if (isinstance(value, (int, float))):
    return (format % value).strip()
  return [(format % v).strip() for v in value]
