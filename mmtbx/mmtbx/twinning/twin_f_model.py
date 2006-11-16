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
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import mmtbx.scaling
import scitbx.math as sm
from mmtbx.scaling import absolute_scaling, relative_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx.scaling import twin_detwin_data, sigmaa_estimation
from mmtbx import masks
from libtbx import table_utils
import scitbx.lbfgs
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
from scitbx import differential_evolution
import sys, os, math, time

class twin_fraction_object(object):
  """provides methods for derivatives and
  transformastion of twin fraction"""
  def __init__(self, twin_fraction = 0):
    self.min_frac = 0.001
    self.max_frac = 0.499
    self.twin_fraction = float(twin_fraction)
    if (self.twin_fraction<=self.min_frac):
      self.twin_fraction = self.min_frac + self.min_frac/10.0
    if (self.twin_fraction>=self.max_frac):
      self.twin_fraction = self.max_frac - self.min_frac/10.0

  def twin_fraction_to_ref( self ):
    tmp = self.twin_fraction - self.min_frac
    tmp = (self.max_frac-self.min_frac)/tmp -1.0
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
               object=None):

    if object is not None:
      k_overall = object.k_overall
      u_star    = object.u_star
      k_sol     = object.k_sol
      u_sol     = object.u_sol
      k_part    = object.k_part
      u_part    = object.u_part
      u_star    = object.u_star
    # this is complete paranoia, making casts just to ensure that one always obtains 2 unique objects
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

  def ref_to_k_part(self, x):
    if x > 10:
      self.k_part = math.exp(10)
    else:
      self.k_part = math.exp( x )

  def ref_to_u_part(self, x):
    self.u_part = math.exp( x )

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
      return None
  def k_part_to_ref(self):
    if self.k_part > 0:
      return math.log( self.k_part )
    else:
      return None
  def u_sol_to_ref(self):
    if self.u_sol > 0:
      return math.log( self.u_sol )
    else:
     return None
  def u_part_to_ref(self):
    if self.u_part>0:
      return math.log( self.u_part )
    else:
      return None
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

class mask_parameter_object(object):
  def __init__(self,
               gridding_n_real=None,
               grid_step=None,
               solvent_radius=1.0,
               shrink_truncation_radius=1.0):
    assert ( [gridding_n_real, grid_step] ).count(None) is 1
    self.gridding_n_real=gridding_n_real
    self.grid_step=grid_step
    self.solvent_radius=solvent_radius
    self.shrink_truncation_radius=shrink_truncation_radius


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
               out=None):
    self.out = out
    if self.out is None:
      self.out = sys.stdout

    self.scaling_parameters=scaling_parameters
    self.target_evaluator=target_evaluator
    self.f_model_core_data=f_model_core_data
    self.twin_fraction_object = twin_fraction_obj
    self.best_score = None
    #first determin the number of parameters please
    self.n = 1+2+1+self.scaling_parameters.n_u_indep
    self.domain = [ (-1,1), (-4,0), (-2,0), (-2,0) ] + [ (-1,1) ]*self.scaling_parameters.n_u_indep
    self.x = flex.double([0]*self.n)
    self.de = differential_evolution.differential_evolution_optimizer(
     self,
     population_size=10,
     f=0.8,
     cr=0.7,
     n_cross=2,
     eps=1e-12,
     show_progress=True
    )

  def update_parameters(self, vector):
    self.scaling_parameters.ref_to_k_overall( vector[1] )
    self.scaling_parameters.ref_to_k_sol( vector[2] )
    self.scaling_parameters.ref_to_u_sol( vector[3] )
    self.scaling_parameters.ref_to_u_star(   list(vector[4:])   )
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
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters)



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

  def compute_gradients_fin_diff(self, h=0.001):
    # obsolete code, will be removed at one stage
    fo = self.compute_functional()
    fd = []
    for ii in xrange( self.x.size() ):
      old = float(self.x[ii])
      self.x[ii] = self.x[ii]+h
      self.update_parameters()
      tmp = self.compute_functional()
      fd.append( tmp )
      self.x[ii] = float(old)
    fd = flex.double(fd)
    fd = (fd-fo)/h
    #print list(fd)
    fd[0] = fd[0]*self.parameter_mask.twin_fraction
    fd[1] = fd[1]*self.parameter_mask.k_overall
    fd[2] = fd[2]*self.parameter_mask.k_sol
    fd[3] = fd[3]*self.parameter_mask.u_sol
    for ii in xrange(3,self.x.size()):
      fd[ii] = fd[ii]*self.parameter_mask.u_star
    #print list(fd)
    return fd

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
        if cycle_count == n_cylce:
          converged=True
      else:
        if last_score is not None:
          if self.best_score_until_now <= last_score:
            converged = True
          else:
            last_score = self.best_score_until_now
        else:
          last_score = self.best_score_until_now
      print cycle_count, converged, self.best_score_until_now , last_score,
      if  last_score is not None:
        print self.best_score_until_now - last_score
      else:
        print
      self.print_it(score,
                    self.scaling_parameters,
                    self.twin_fraction)



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
    de = de_bulk_solvent_scaler(
               self.scaling_parameters,
               self.twin_fraction,
               self.target_evaluator,
               self.f_model_core_data,
               out = self.out)
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


class twin_model_manager(object):
  def __init__(self,
               f_obs_array=None,
               free_array=None,
               xray_structure=None,
               scaling_parameters=None,
               mask_parameters=None,
               out=None,
               twin_law=None,
               start_fraction=0.1,
               n_refl_bin=2000,
               max_bins=20):
    self.out = out
    if self.out is None:
      self.out = sys.stdout
    self.twin_fraction_object = twin_fraction_object(twin_fraction=start_fraction)
    self.twin_law=twin_law
    self.twin_fraction=start_fraction
    assert (self.twin_law is not None)
    self.f_obs_array = f_obs_array.map_to_asu()
    self.free_array = free_array.map_to_asu()

    self.f_obs_array_work = self.f_obs_array.select( ~self.free_array.data() )
    self.f_obs_array_free = self.f_obs_array.select( self.free_array.data() )

    #setup the binners if this has not been done yet
    if self.f_obs_array.binner() is None:
      if self.f_obs_array.indices().size()/float(n_refl_bin) > max_bins:
        self.f_obs_array.setup_binner(n_bins = max_bins)
      else:
        self.f_obs_array.setup_binner( reflections_per_bin=n_refl_bin )

    self.f_obs_array_work.use_binning_of( self.f_obs_array )
    self.f_obs_array_free.use_binning_of( self.f_obs_array )

    self.xray_structure = xray_structure
    self.xs = crystal.symmetry( unit_cell=f_obs_array.unit_cell(),
                                space_group=f_obs_array.space_group() )
    self.scaling_parameters = scaling_parameters
    if self.scaling_parameters is None:
      self.scaling_parameters = scaling_parameters_object(self.xs)

    self.mask_parameters = mask_parameters
    if self.mask_parameters is None:
      self.mask_parameters = mask_parameter_object(
        grid_step=f_obs_array.d_min()/4.0 )

    #-------------------
    self.miller_set = None
    self.f_atoms = None
    self.miller_set = None
    print "updating f atoms"
    self.update_f_atoms()

    #-------------------
    self.mask = None
    self.f_mask = None
    print "updating mask"
    self.update_f_mask()
    #-------------------
    self.f_partial = None

    #-------------------
    print "make data core"
    self.data_core = xray.f_model_core_data(
      hkl = self.f_atoms.indices(),
      f_atoms= self.f_atoms.data(),
      f_mask = self.f_mask.data(),
      unit_cell = self.f_atoms.unit_cell(),
      k_overall=self.scaling_parameters.k_overall,
      u_star=self.scaling_parameters.u_star,
      k_sol=self.scaling_parameters.k_sol,
      u_sol=self.scaling_parameters.u_sol,
      f_part=None,
      k_part=self.scaling_parameters.k_part,
      u_part=self.scaling_parameters.u_part )


    print "Making target evaluators (free and work)"
    self.target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs       = self.f_obs_array_work.indices(),
      f_obs         = self.f_obs_array_work.data(),
      w_obs         = self.f_obs_array_work.sigmas(),
      hkl_calc      = self.f_atoms.indices(),
      space_group   = self.f_obs_array.space_group(),
      anomalous_flag= self.f_obs_array.anomalous_flag(),
      alpha         = self.twin_fraction,
      twin_law      = self.twin_law.as_double_array()[0:9] )

    self.free_target_evaluator = xray.least_squares_hemihedral_twinning_on_f(
      hkl_obs        = self.f_obs_array_free.indices(),
      f_obs          = self.f_obs_array_free.data(),
      w_obs          = self.f_obs_array_free.sigmas(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array.space_group(),
      anomalous_flag = self.f_obs_array.anomalous_flag(),
      alpha          = self.twin_fraction,
      twin_law       = self.twin_law.as_double_array()[0:9] )

    print "Making bulk solvent scaler"
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
    print "Making r-value objects"
    self.r_work_object = xray.hemihedral_r_values(
      hkl_obs        = self.f_obs_array_work.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array_work.space_group(),
      anomalous_flag = self.f_obs_array.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    self.r_free_object = xray.hemihedral_r_values(
      hkl_obs        = self.f_obs_array_free.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array_free.space_group(),
      anomalous_flag = self.f_obs_array_free.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

    print "Making detwinning objects"
    self.work_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs_array_work.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array_work.space_group(),
      anomalous_flag = self.f_obs_array_work.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )
    self.free_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs_array_free.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array_free.space_group(),
      anomalous_flag = self.f_obs_array_free.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )
    self.full_detwinner = xray.hemihedral_detwinner(
      hkl_obs        = self.f_obs_array.indices(),
      hkl_calc       = self.f_atoms.indices(),
      space_group    = self.f_obs_array.space_group(),
      anomalous_flag = self.f_obs_array.anomalous_flag(),
      twin_law       = self.twin_law.as_double_array()[0:9] )

  def f_model(self):
    tmp_f_model = self.f_atoms.customized_copy(
      data = self.data_core.f_model()
    )
    return tmp_f_model


  def update_bulk_solvent_parameters(self,
                                     bulk_solvent_parameters=None,
                                     twin_fraction_parameters=None,
                                     refine=False,
                                     grid_search=False,
                                     initialise=False,
                                     de_search=False,
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
      self.bss.de_search()
      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters
    if grid_search:
      self.bss.grid_search()
      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters
    if refine:
      self.bss.scale_it()
      self.twin_fraction_object = self.bss.best_twin_fraction
      self.scaling_parameters = self.bss.best_scaling_parameters

    self.target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    self.free_target_evaluator.alpha( self.twin_fraction_object.twin_fraction )
    self.data_core.koverall( self.scaling_parameters.k_overall )
    self.data_core.ustar( self.scaling_parameters.u_star )
    self.data_core.ksol( self.scaling_parameters.k_sol )
    self.data_core.usol( self.scaling_parameters.u_sol )


  def update_f_atoms(self):
    """Get f calc from the xray structure"""
    if self.miller_set is None:
      completion = xray.twin_completion( self.f_obs_array.indices(),
                                         self.xs.space_group(),
                                         self.f_obs_array.anomalous_flag(),
                                         self.twin_law.as_double_array()[0:9] )
      indices = completion.twin_complete()
      self.miller_set = miller.set(
        crystal_symmetry = self.xs,
        indices =indices,
        anomalous_flag = self.f_obs_array.anomalous_flag() ).map_to_asu()
      assert self.miller_set.is_unique_set_under_symmetry()
    tmp = self.miller_set.structure_factors_from_scatterers(
      xray_structure = self.xray_structure )

    self.f_atoms = tmp.f_calc()
    if self.miller_set is None:
      self.miller_set = tmp.miller_set()


  def update_f_mask(self):
    self.mask = masks.bulk_solvent(
      xray_structure=self.xray_structure,
      gridding_n_real=self.mask_parameters.gridding_n_real,
      grid_step=self.mask_parameters.grid_step,
      solvent_radius=self.mask_parameters.solvent_radius,
      shrink_truncation_radius=self.mask_parameters.shrink_truncation_radius)

    self.f_mask = self.mask.structure_factors( self.miller_set )


  def update_xray_structure(self,
                            xray_structure,
                            update_mask=False):
    self.xray_structure=xray_structure
    self.update_f_atoms()
    if update_mask:
      self.update_f_mask()

  def r_values(self, table=True):
    r_abs_work_f_overall = self.r_work_object.r_amplitude_abs(
      f_obs         = self.f_obs_array_work.data(),
      f_model       = self.data_core.f_model(),
      selection     = None,
      twin_fraction = self.twin_fraction_object.twin_fraction)

    r_abs_free_f_overall = self.r_free_object.r_amplitude_abs(
      self.f_obs_array_free.data(),
      self.data_core.f_model(),
      None,
      self.twin_fraction_object.twin_fraction)

    if table:
      r_abs_work_f_bin = []
      r_abs_free_f_bin = []
      bin_low = []
      bin_high= []
      n_free = []
      n_work = []
      rows = []
      for i_bin in self.f_obs_array_free.binner().range_used():
        selection = flex.bool( self.f_obs_array_work.binner().bin_indices() == i_bin )
        n_work = selection.count(True)
        tmp_work = self.r_work_object.r_amplitude_abs(
          f_obs         = self.f_obs_array_work.data(),
          f_model       = self.data_core.f_model(),
          selection     = selection,
          twin_fraction = self.twin_fraction_object.twin_fraction)
        selection = flex.bool( self.f_obs_array_free.binner().bin_indices() == i_bin )
        n_free = selection.count(True)
        tmp_free = self.r_free_object.r_amplitude_abs(
          f_obs         = self.f_obs_array_free.data(),
          f_model       = self.data_core.f_model(),
          selection     = selection,
          twin_fraction = self.twin_fraction_object.twin_fraction)

        r_abs_work_f_bin.append(tmp_work)
        r_abs_free_f_bin.append(tmp_free)
        d_max,d_min = self.f_obs_array_work.binner().bin_d_range( i_bin )
        tmp = [ str( "%3i"%(i_bin)    ),
                str( "%5.2f"%(d_max)  ),
                str( "%5.2f"%(d_min)  ),
                str( "%5i"%(n_work)   ),
                str( "%3.2f"%(tmp_work) ),
                str( "%5i"%(n_free)   ),
                str( "%3.2f"%(tmp_free) ) ]

        rows.append( tmp )

      header = ("bin", "d_max", "d_min", "n_work", "r_work", "n_free", "r_free")
      comments = """
Overall r values
R Work : %4.3f
R Free : %4.3f

R  = \sum_h( |Ft-Fo| )/ \sum_h(Fo)
Ft = Sqrt(tf*F1^2 + (1-tf)F2^2)
F1,F2 are twin related model amplitudes.
tf is the twin fractrion and Fo is an observed amplitude."""%(r_abs_work_f_overall, r_abs_free_f_overall)

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
    else:
      return r_abs_work_f_overall, r_abs_free_f_overall

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
      self.update_bulk_solvent_parameters( twin_fraction_parameters =  tmp_twin_fraction )
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
    self.update_bulk_solvent_parameters( twin_fraction_parameters =  current_twin_fraction )



  def target(self, print_it=True):
    tmp_w=self.target_evaluator.target( self.data_core.f_model() )
    tmp_f=self.free_target_evaluator.target( self.data_core.f_model() )
    if print_it:
      print >> self.out
      print >> self.out, "----------------- Target values -----------------"
      print >> self.out, " Basic values "
      print >> self.out, "   working set  : %8.6e "%(tmp_w)
      print >> self.out, "   free set     : %8.6e "%(tmp_f)
      print >> self.out
      print >> self.out, " Target values devided by number of contributers"
      print >> self.out, "    working set : %8.6e "%(tmp_w/self.f_obs_array_work.data().size() )
      print >> self.out, "    free set    : %8.6e "%(tmp_f/self.f_obs_array_free.data().size() )
      print >> self.out, "-------------------------------------------------"
    else:
      return(tmp_w,tmp_f)


  def map_coefficients(self, apply_weight=True):
    print >> self.out
    print >> self.out, "--------------- Constructing map coefficients ---------------"
    #first make a detwinned fobs
    tmp_i_obs = self.f_obs_array.f_as_f_sq()
    dt_iobs, dt_isigma = self.full_detwinner.detwin_with_model_data(
      tmp_i_obs.data(),
      tmp_i_obs.sigmas(),
      self.data_core.f_model(),
      self.twin_fraction_object.twin_fraction )
    tmp_i_obs = tmp_i_obs.customized_copy(
      data = dt_iobs,
      sigmas = dt_isigma ).set_observation_type( tmp_i_obs )
    dt_f_obs = tmp_i_obs.f_sq_as_f()

    # get the proper set of fmodel data please
    tmp_f_model = self.f_atoms.customized_copy(
      data = self.data_core.f_model() ).common_set(
        dt_f_obs )
    tmp_abs_f_model = tmp_f_model.customized_copy(
      data = flex.abs( tmp_f_model.data()) ).set_observation_type( dt_f_obs )
    # we now have to perform a scaling of the detwinned data to the calculated data
    # we will use a no-nonsense local scalnig implemented elsewheer
    print >> self.out
    print >> self.out, "Local scaling of detwinned data with model data"
    print >> self.out
    local_scaler = relative_scaling.local_scaling_driver(
      miller_native=tmp_abs_f_model,
      miller_derivative=dt_f_obs,
      use_intensities=False,
      local_scaling_dict={'local_nikonov':True, 'local_moment':False, 'local_lsq':False} )

    dt_f_obs = dt_f_obs.customized_copy(
      data =  dt_f_obs.data()*local_scaler.local_scaler.get_scales()
      ).set_observation_type( tmp_i_obs )

    m = None
    alpha_d = None
    twofofc_map = None
    fofc_map = None

    # get coefficients for a gradient map please
    gradients = self.target_evaluator.d_target_d_fmodel(
      self.data_core.f_model() )

    gradients = self.f_atoms.customized_copy(
      data = -flex.conj(gradients) ).common_set( dt_f_obs )

    sigmaa_object = sigmaa_estimation.sigmaa_estimator(
      miller_obs   = dt_f_obs,
      miller_calc  = tmp_f_model,
      r_free_flags = self.free_array )
    print >> self.out
    print >> self.out, " Sigma estimation on detwinned data using model info. "
    sigmaa_object.show(out=self.out)
    d = sigmaa_object.alpha
    m = sigmaa_object.fom

    wtwofofc_map = dt_f_obs.customized_copy(
      data=2.0*m.data()*dt_f_obs.data()-d.data()*flex.abs(tmp_f_model.data()) )
    twofofc_map = dt_f_obs.customized_copy(
      data=2.0*dt_f_obs.data()-flex.abs(tmp_f_model.data()) )

    # phase transfers
    twofofc_map = twofofc_map.phase_transfer( tmp_f_model )
    wtwofofc_map = wtwofofc_map.phase_transfer( tmp_f_model )
    print >> self.out, "-------------------------------------------------------------"
    return twofofc_map, wtwofofc_map, gradients
