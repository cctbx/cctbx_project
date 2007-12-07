from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys, copy
from libtbx.test_utils import approx_equal
from scitbx import lbfgs
from cctbx import xray

class manager(object):
  def __init__(self, fmodel,
                     model,
                     max_number_of_iterations    = 25,
                     number_of_macro_cycles      = 3,
                     occypancy_max               = 1,
                     occupancy_min               = 0,
                     log                         = None):
    scatterers = fmodel.xray_structure.scatterers()
    scatterers.flags_set_grads(state=False)
    sel = flex.bool(scatterers.size(), True).iselection()
    scatterers.flags_set_grad_occupancy(iselection = sel)
    fmodel_copy = fmodel.deep_copy()
    fmodel.xray_structure.adjust_occupancy(occ_max = occypancy_max,
                                           occ_min = occupancy_min)
    par_initial = fmodel.xray_structure.scatterers().extract_occupancies()
    minimized = None
    r_work_start = fmodel.r_work()
    for macro_cycle in xrange(number_of_macro_cycles):
      if(minimized is not None): par_initial = minimized.par_min
      minimized = group_minimizer(
        fmodel                   = fmodel_copy,
        model                    = model,
        par_initial              = par_initial,
        max_number_of_iterations = max_number_of_iterations)
      if(minimized is not None): par_initial = minimized.par_min
      fmodel.xray_structure.set_occupancies(par_initial)
      fmodel.xray_structure.adjust_occupancy(occ_max = occypancy_max,
                                             occ_min = occupancy_min)
      fmodel_copy.update_xray_structure(xray_structure = fmodel.xray_structure,
        update_f_calc  = True)
      r_work = minimized.fmodel.r_work()
      assert approx_equal(r_work, fmodel_copy.r_work())
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    self.fmodel = fmodel


class group_minimizer(object):
  def __init__(self,
               fmodel,
               model,
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
      ignore_line_search_failed_step_at_upper_bound = True))
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, par):
    return pack_unpack(x = par, table =
      self.model.atom_i_seqs_grouped_by_altconf)

  def unpack_x(self):
    self.par_min = pack_unpack(x = self.x, table =
      self.model.atom_i_seqs_grouped_by_altconf)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.fmodel.xray_structure.set_occupancies(self.par_min)
    self.fmodel.update_xray_structure(update_f_calc = True)
    t_r = self.target_functor(compute_gradients= True)
    self.f = t_r.target_work()
    self.g = pack_gradients(x = t_r.gradients_wrt_atomic_parameters(
      occupancy = True), table = self.model.atom_i_seqs_grouped_by_altconf)
    return self.f, self.g

def pack_unpack(x, table):
  result = flex.double(x.size(), 0)
  for indices in table:
    if(len(indices) == 1):
      i0 = indices[0]
      result[i0] = x[i0]
    elif(len(indices) == 2):
      i0,i1 = indices
      result[i0] = x[i0]
      result[i1] = 1. - x[i0]
    elif(len(indices) == 3):
      i0,i1,i2 = indices
      result[i0] = x[i0]
      result[i1] = 1. - x[i0] - x[i2]
      result[i2] = 1. - x[i0] - x[i1]
    elif(len(indices) == 4):
      i0,i1,i2,i3 = indices
      result[i0] = x[i0]
      result[i1] = 1. - x[i0] - x[i2] - x[i3]
      result[i2] = 1. - x[i0] - x[i1] - x[i3]
      result[i3] = 1. - x[i0] - x[i1] - x[i2]
    else:
      raise RuntimeError("Exceed maximum allowable number of conformers (=4).")
  return result

def pack_gradients(x, table):
  result = flex.double(x.size(), 0)
  for indices in table:
    if(len(indices) == 1):
      i0 = indices[0]
      result[i0] = x[i0]
    elif(len(indices) == 2):
      i0,i1 = indices
      result[i0] = x[i0] - x[i1]
      result[i1] =-x[i1] # ??? zero or this value?
    elif(len(indices) == 3):
      i0,i1,i2 = indices
      result[i0] =  x[i0] - x[i1] - x[i2]
      result[i1] = -x[i2]
      result[i2] = -x[i1]
    elif(len(indices) == 4):
      i0,i1,i2,i3 = indices
      result[i0] = x[i0] - x[i2]
      result[i1] = x[i1] - x[i2]
      result[i2] = x[i2] - x[i1] - x[i3]
      result[i3] = x[i3] - x[i2]
    else:
      raise RuntimeError("Exceed maximum allowable number of conformers (=4).")
  return result
