from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys, copy
from libtbx.test_utils import approx_equal
from scitbx import lbfgs
from cctbx import xray
from libtbx.str_utils import format_value

class manager(object):
  def __init__(self, fmodel,
                     model,
                     max_number_of_iterations    = 25,
                     number_of_macro_cycles      = 3,
                     occupancy_max               = 1,
                     occupancy_min               = 0,
                     r_increase_tolerance        = 0.001,
                     log                         = None):
    self.show(fmodel=fmodel, log= log, message =
      "individual occupancy refinement: start")
    selection = model.refinement_flags.occupancies_individual
    i_selection = flex.size_t()
    for s in selection:
      for ss in s:
        i_selection.append(ss)
    b_selection = flex.bool(fmodel.xray_structure.scatterers().size(),
      i_selection)
    fmodel.xray_structure.scatterers().flags_set_grad_occupancy(
      iselection = i_selection)
    has_alt_conf = False
    for i_seq_list in selection:
      if(len(i_seq_list) > 1):
        has_alt_conf = True
        break
    fmodel.xray_structure.adjust_occupancy(occ_max   = occupancy_max,
                                           occ_min   = occupancy_min,
                                           selection = b_selection)
    xray_structure_dc = fmodel.xray_structure.deep_copy_scatterers()
    par_initial = fmodel.xray_structure.scatterers().extract_occupancies()
    minimized = None
    r_work_start = fmodel.r_work()
    for macro_cycle in xrange(number_of_macro_cycles):
      if(minimized is not None): par_initial = minimized.par_min
      minimized = group_minimizer(
        fmodel                   = fmodel,
        model                    = model,
        b_selection              = b_selection,
        selection                = selection,
        par_initial              = par_initial,
        has_alt_conf             = has_alt_conf,
        max_number_of_iterations = max_number_of_iterations)
      if(minimized is not None): par_initial = minimized.par_min
      fmodel.xray_structure.set_occupancies(par_initial)
      fmodel.xray_structure.adjust_occupancy(occ_max   = occupancy_max,
                                             occ_min   = occupancy_min,
                                             selection = b_selection)
      r_work = minimized.fmodel.r_work()
      assert approx_equal(r_work, fmodel.r_work())
    delta = r_work_start - r_work
    if(abs(delta) > abs(r_increase_tolerance) and delta < 0):
      xray_structure_final = xray_structure_dc
    else:
      xray_structure_final = fmodel.xray_structure
    model.xray_structure = xray_structure_final
    fmodel.update_xray_structure(xray_structure = xray_structure_final,
                                 update_f_calc  = True)
    self.fmodel = fmodel
    refined_occ = self.fmodel.xray_structure.scatterers().extract_occupancies(
      ).select(i_selection)
    assert flex.min(refined_occ) >= occupancy_min
    assert flex.max(refined_occ) <= occupancy_max
    self.show(fmodel=fmodel, log= log, message =
      "individual occupancy refinement: end")

  def show(self, fmodel, message, log):
    r_work = format_value("%6.4f", fmodel.r_work())
    r_free = format_value("%6.4f", fmodel.r_free())
    target = format_value("%-13.3f", fmodel.target_w())
    target_name = format_value("%s", fmodel.target_name)
    occupancies = fmodel.xray_structure.scatterers().extract_occupancies()
    occ_max = format_value("%4.2f", flex.max(occupancies))
    occ_min = format_value("%4.2f", flex.min(occupancies))
    number_small = format_value("%8d", (occupancies < 0.1).count(True))
    print >> log, "|-"+message+"-"*(79-len("|-"+message+"|"))+"|"
    p1 = "| r_work = %s r_free = %s" % (r_work, r_free)
    p2 = "target_work(%s) = %s |" % (target_name, target)
    print >> log, p1+" "*(79-len(p1+p2))+p2
    print >> log, \
      "| occupancies: max = %s  min = %s   number of occupancies < 0.1: %s |"%(
      occ_max, occ_min, number_small)
    print >> log, "|"+"-"*77+"|"

class group_minimizer(object):
  def __init__(self,
               fmodel,
               model,
               selection,
               b_selection,
               par_initial,
               has_alt_conf,
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
    if(self.has_alt_conf):
      return pack_unpack(x = par, table = self.selection)
    else:
      return par

  def unpack_x(self):
    if(self.has_alt_conf):
      self.par_min = pack_unpack(x = self.x, table = self.selection)
    else:
      self.par_min = self.x

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.fmodel.xray_structure.set_occupancies(self.par_min, self.b_selection)
    self.fmodel.update_xray_structure(update_f_calc = True)
    t_r = self.target_functor(compute_gradients= True)
    self.f = t_r.target_work()
    g =  t_r.gradients_wrt_atomic_parameters(occupancy = True)
    if(self.has_alt_conf):
      self.g = pack_gradients(x = g, table = self.selection)
    else:
      self.g = flex.double(g.size(), 0).set_selected(self.b_selection, g)
    return self.f, self.g

def pack_unpack(x, table):
  result = x.deep_copy()
  #result = flex.double(x.size(), 0)  # why not this ?
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
