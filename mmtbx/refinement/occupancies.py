from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys, copy
from libtbx.test_utils import approx_equal
from scitbx import lbfgs
from cctbx import xray
from libtbx.str_utils import format_value
from phenix.refinement import weight_xray_chem

class manager(object):
  def __init__(self, fmodels,
                     model,
                     max_number_of_iterations    = 25,
                     number_of_macro_cycles      = 3,
                     occupancy_max               = None,
                     occupancy_min               = None,
                     log                         = None):
    self.show(fmodels=fmodels, log= log, message="occupancy refinement: start")
    fmodels.update_xray_structure(xray_structure = model.xray_structure,
                                  update_f_calc  = True)
    selections = model.refinement_flags.s_occupancies
    i_selection = flex.size_t()
    for s in selections:
      for ss in s:
        i_selection.extend(ss)
    fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
      state=False)
    fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grad_occupancy(
      iselection = i_selection)
    fmodels.fmodel_xray().xray_structure.adjust_occupancy(
      occ_max   = occupancy_max,
      occ_min   = occupancy_min,
      selection = i_selection)
    xray_structure_dc = fmodels.fmodel_xray().xray_structure.\
      deep_copy_scatterers()
    par_initial = flex.double()
    occupancies = xray_structure_dc.scatterers().extract_occupancies()
    constrained_groups_selections = []
    group_counter = 0
    for sel in selections:
      ss = []
      for sel_ in sel:
        ss.append(group_counter)
        group_counter += 1
        val = flex.mean(occupancies.select(sel_))
        par_initial.append(val)
      constrained_groups_selections.append(ss)
    minimized = None
    r_work_start = fmodels.fmodel_xray().r_work()
    for macro_cycle in xrange(number_of_macro_cycles):
      if(minimized is not None): par_initial = minimized.par_min
      minimized = minimizer(
        fmodels                       = fmodels,
        selections                    = selections,
        constrained_groups_selections = constrained_groups_selections,
        par_initial                   = par_initial,
        max_number_of_iterations      = max_number_of_iterations)
      if(minimized is not None): par_initial = minimized.par_min
      set_refinable_parameters(
        xray_structure = fmodels.fmodel_xray().xray_structure,
        parameters     = par_initial,
        selections     = selections)
      fmodels.fmodel_xray().xray_structure.adjust_occupancy(
        occ_max   = occupancy_max,
        occ_min   = occupancy_min,
        selection = i_selection)
    xray_structure_final = fmodels.fmodel_xray().xray_structure
    model.xray_structure = xray_structure_final
    fmodels.update_xray_structure(xray_structure = xray_structure_final,
                                  update_f_calc  = True)
    refined_occ = xray_structure_final.scatterers().extract_occupancies().\
      select(i_selection)
    assert flex.min(refined_occ) >= occupancy_min
    assert flex.max(refined_occ) <= occupancy_max
    self.show(fmodels= fmodels, log = log, message="occupancy refinement: end")

  def show(self, fmodels, message, log):
    if(log is not None):
      print >> log, "|-"+message+"-"*(79-len("|-"+message+"|"))+"|"
      fm_x, fm_n = fmodels.fmodel_xray(), fmodels.fmodel_neutron()
      if(fm_n is not None):
        print >> log, "|"+" "*36+"X-ray"+" "*36+"|"
      self.show_helper(fmodel = fm_x, log = log)
      if(fm_n is not None):
        print >> log, "|"+" "*35+"neutron"+" "*35+"|"
        self.show_helper(fmodel = fm_n, log = log)
      occupancies = fm_x.xray_structure.scatterers().extract_occupancies()
      occ_max = format_value("%4.2f", flex.max(occupancies))
      occ_min = format_value("%4.2f", flex.min(occupancies))
      number_small = format_value("%8d", (occupancies < 0.1).count(True))
      print >> log, \
        "| occupancies: max = %s  min = %s   number of occupancies < 0.1: %s |"%(
        occ_max, occ_min, number_small)
      print >> log, "|"+"-"*77+"|"

  def show_helper(self, fmodel, log):
    r_work = format_value("%6.4f", fmodel.r_work())
    r_free = format_value("%6.4f", fmodel.r_free())
    target = format_value("%-13.3f", fmodel.target_w())
    target_name = format_value("%s", fmodel.target_name)
    p1 = "| r_work = %s r_free = %s" % (r_work, r_free)
    p2 = "target_work(%s) = %s |" % (target_name, target)
    print >> log, p1+" "*(79-len(p1+p2))+p2

class minimizer(object):
  def __init__(self,
               fmodels,
               constrained_groups_selections,
               selections,
               par_initial,
               max_number_of_iterations):
    adopt_init_args(self, locals())
    self.fmodels.create_target_functors()
    self.fmodels.prepare_target_functors_for_minimization()
    self.weights = weight_xray_chem.weights(wx       = 1,
                                            wx_scale = 1,
                                            angle_x  = None,
                                            wn       = 1,
                                            wn_scale = 1,
                                            angle_n  = None,
                                            w        = 0,
                                            wxn      = 1) # XXX
    self.par_min = self.par_initial.deep_copy()
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
    return pack_unpack(x = par, table = self.constrained_groups_selections)

  def unpack_x(self):
    self.par_min = pack_unpack(x = self.x,
      table = self.constrained_groups_selections)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    set_refinable_parameters(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      parameters     = self.par_min,
      selections     = self.selections)
    self.fmodels.update_xray_structure(update_f_calc = True)
    fmodels_target_and_gradients = self.fmodels.target_and_gradients(
      weights           = self.weights,
      compute_gradients = True,
      occupancy         = True)
    self.f = fmodels_target_and_gradients.target()
    g =  fmodels_target_and_gradients.gradients()
    # do group grads first then pack for constraints
    grads_wrt_par = flex.double()
    for sel in self.selections:
      for sel_ in sel:
        grads_wrt_par.append(flex.sum( g.select(sel_) ))
    # now apply constraints
    self.g = pack_gradients(x = grads_wrt_par,
      table = self.constrained_groups_selections)
    return self.f, self.g

def set_refinable_parameters(xray_structure, parameters, selections):
  sz = xray_structure.scatterers().size()
  i = 0
  for sel in selections:
    for sel_ in sel:
      sel__b = flex.bool(sz, flex.size_t(sel_))
      xray_structure.set_occupancies(parameters[i], sel__b)
      i+=1

def pack_unpack(x, table) :
  result = x.deep_copy()
  for indices in table :
    if (len(indices) == 1) :
      i0 = indices[0]
      result[i0] = x[i0]
    else :
      xsum = 0
      for i in indices[0:-1] :
        result[i] = x[i]
        xsum += x[i]
      result[indices[-1]] = 1. - xsum
  return result

def pack_gradients (x, table) :
  result = flex.double(x.size(), 0)
  for indices in table :
    if(len(indices) == 1):
      i0 = indices[0]
      result[i0] = x[i0]
    elif(len(indices) == 2):
      i0,i1 = indices
      result[i0] = x[i0] - x[i1]
      result[i1] =-x[i1] # ??? zero or this value?
    else :
      for i in indices[0:-1] :
        result[i] = x[i] - x[indices[-1]]
      result[indices[-1]] = 0
  return result
