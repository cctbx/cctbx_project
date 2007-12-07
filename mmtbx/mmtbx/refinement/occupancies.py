from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
import copy, math
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time


class manager(object):
  def __init__(self, fmodel,
                     model,
                     max_number_of_iterations    = 25,
                     number_of_macro_cycles      = 5,
                     log                         = None):
    scatterers = fmodel.xray_structure.scatterers()
    scatterers.flags_set_grads(state=False)
    sel = flex.bool(scatterers.size(), True).iselection()
    scatterers.flags_set_grad_occupancy(iselection = sel)
    fmodel_copy = fmodel.deep_copy()
    sc_start = fmodel.xray_structure.scatterers().deep_copy()
    fmodel.xray_structure.adjust_occupancy(occ_max = 1,
                                             occ_min = 0.0)
    par_initial = sc_start.extract_occupancies()
    minimized = None
    for macro_cycle in xrange(3):
      if(minimized is not None): par_initial = minimized.par_min
      minimized = group_minimizer(
                   fmodel                      = fmodel_copy,
                   model                       = model,
                   sc_start                    = sc_start,
                   par_initial                 = par_initial,
                   max_number_of_iterations    = max_number_of_iterations)
      if(minimized is not None): par_initial = minimized.par_min
      apply_transformation(
        xray_structure = fmodel.xray_structure,
        par            = par_initial,
        sc_start       = sc_start)
      #fmodel.xray_structure.adjust_occupancy(occ_max = 1,
      #                                       occ_min = 0.1)
      fmodel_copy.update_xray_structure(
        xray_structure = fmodel.xray_structure,
        update_f_calc  = True,
        out            = log)
      #fmodel.xray_structure.adjust_occupancy(occ_max = 1,
      #                                       occ_min = 0)
      rwork = minimized.fmodel.r_work()
      rfree = minimized.fmodel.r_free()
      assert approx_equal(rwork, fmodel_copy.r_work())
    ###
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    self.fmodel = fmodel


class group_minimizer(object):
  def __init__(self,
               fmodel,
               model,
               sc_start,
               par_initial,
               max_number_of_iterations):
    adopt_init_args(self, locals())
    self.target_functor = fmodel.target_functor()
    self.par_min = copy.deepcopy(self.par_initial)
    self.x = self.pack(self.par_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
               target_evaluator = self,
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_number_of_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True,
                    ignore_line_search_failed_step_at_upper_bound = True,)
                              )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, par):
    #par = flex.double([par[0], 1-par[0]-par[2], 1-par[0]-par[1],    par[3], 1-par[3]])
    #par = flex.double([par[0], 1-par[0]-par[2]-par[3], 1-par[0]-par[1]-par[3], 1-par[0]-par[1]-par[2],
    #                   par[4], 1-par[4]])
    #return flex.double(tuple(par))
    return pack_unpack(x = par, table = self.model.atom_i_seqs_grouped_by_altconf)

  def unpack_x(self):
    #x = self.x
    #self.par_min = flex.double([x[0], 1-x[0]-x[2], 1-x[0]-x[1],    x[3], 1-x[3]])
    #self.par_min = flex.double([x[0], 1-x[0]-x[2]-x[3], 1-x[0]-x[1]-x[3], 1-x[0]-x[1]-x[2],
    #                            x[4], 1-x[4]])
    self.par_min = pack_unpack(x = self.x, table = self.model.atom_i_seqs_grouped_by_altconf)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    apply_transformation(
      xray_structure = self.fmodel.xray_structure,
      par            = self.par_min,
      sc_start       = self.sc_start)
    self.fmodel.update_xray_structure(update_f_calc=True)
    tg_obj = target_and_grads(target_functor = self.target_functor)
    self.f = tg_obj.target()
    self.g = flex.double(tg_obj.gradients_wrt_par())
    ###
    #g = self.g
    #g = flex.double([g[0]-g[1]-g[2],-g[2],-g[1],     g[3]-g[4],-g[4]])
    #g = flex.double([g[0]-g[1]-g[2],-g[2],-g[1],     g[3]-g[4],0])
    #g = flex.double([g[0]-g[2], g[1]-g[2], -g[1]+g[2]-g[3], -g[2]+g[3],
    #                 g[4]-g[5],-g[5]])
    #self.g = g
    ###
    self.g = pack_gradients(x = self.g, table = self.model.atom_i_seqs_grouped_by_altconf)
    #print list(g)
    #print list(self.g)
    #print
    return self.f, self.g

def pack_unpack(x, table):
  result = flex.double(x.size(), 0)
  for indices in table:
    if(len(indices) == 1):
      result[indices[0]] = x[indices[0]]
    elif(len(indices) == 2):
      result[indices[0]] = x[indices[0]]
      result[indices[1]] = 1. - x[indices[0]]
    elif(len(indices) == 3):
      result[indices[0]] = x[indices[0]]
      result[indices[1]] = 1. - x[indices[0]] - x[indices[2]]
      result[indices[2]] = 1. - x[indices[0]] - x[indices[1]]
    elif(len(indices) == 4):
      result[indices[0]] = x[indices[0]]
      result[indices[1]] = 1. - x[indices[0]] - x[indices[2]] - x[indices[3]]
      result[indices[2]] = 1. - x[indices[0]] - x[indices[1]] - x[indices[3]]
      result[indices[3]] = 1. - x[indices[0]] - x[indices[1]] - x[indices[2]]
  return result

def pack_gradients(x, table):
  result = flex.double(x.size(), 0)
  for indices in table:
    if(len(indices) == 1):
      result[indices[0]] = x[indices[0]]
    elif(len(indices) == 2):
      result[indices[0]] = x[indices[0]] - x[indices[1]]
      result[indices[1]] =    -x[indices[1]] # ??? zero or this value?
    elif(len(indices) == 3):
      result[indices[0]] =  x[indices[0]]-x[indices[1]]-x[indices[2]]
      result[indices[1]] = -x[indices[2]]
      result[indices[2]] = -x[indices[1]]
    elif(len(indices) == 4):
      result[indices[0]] = x[indices[0]] - x[indices[2]]
      result[indices[1]] = x[indices[1]] - x[indices[2]]
      result[indices[2]] = x[indices[2]] - x[indices[1]] - x[indices[3]]
      result[indices[3]] = x[indices[3]] - x[indices[2]]
  return result


def apply_transformation(xray_structure,
                         par,
                         sc_start):
  new_sc = sc_start.deep_copy()
  for i, sc in enumerate(new_sc):
    sc.occupancy = par[i]
  xray_structure.replace_scatterers(new_sc)

class target_and_grads(object):
  def __init__(self, target_functor):
    t_r = target_functor(compute_gradients= True)
    self.f = t_r.target_work()
    target_grads_wrt_par = t_r.gradients_wrt_atomic_parameters(
      occupancy = True)
    self.grads_wrt_par = target_grads_wrt_par

  def target(self):
    return self.f

  def gradients_wrt_par(self):
    return self.grads_wrt_par
