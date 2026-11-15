
from __future__ import absolute_import, division, print_function
from mmtbx.utils import get_atom_selections
import mmtbx.utils
from cctbx.array_family import flex
from libtbx.str_utils import format_value, make_sub_header
from libtbx.utils import Sorry, null_out
from six.moves import range

class manager(object):
  def __init__(self, fmodels,
                     model,
                     max_number_of_iterations    = 25,
                     number_of_macro_cycles      = 3,
                     occupancy_max               = None,
                     occupancy_min               = None,
                     log                         = None,
                     exclude_hd                  = False):
    self.show(fmodels=fmodels, log= log, message="occupancy refinement: start")
    fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
                                  update_f_calc  = True)
    selections = model.refinement_flags.s_occupancies
    # exclude H or D from refinement if requested
    if(exclude_hd):
      hd_sel = model.get_hd_selection()
      tmp_sel = []
      for sel in selections:
        tmp_sel_ = []
        for sel_ in sel:
          tmp_sel__ = flex.size_t()
          for sel__ in sel_:
            if(not hd_sel[sel__]):
              tmp_sel__.append(sel__)
          if(tmp_sel__.size()>0):
            tmp_sel_.append(tmp_sel__)
        if(len(tmp_sel_)>0):
          tmp_sel.append(tmp_sel_)
      selections = tmp_sel
    #
    if(len(selections)>0):
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
      occupancies = fmodels.fmodel_xray().xray_structure.scatterers(
        ).extract_occupancies()
      par_initial = flex.double()
      for sel in selections:
        for sel_ in sel:
          par_initial.append(flex.mean(occupancies.select(sel_)))
      n_iterations_total, n_fun_total = 0, 0
      for macro_cycle in range(number_of_macro_cycles):
        minimized = minimizer(
          fmodels                  = fmodels,
          selections               = selections,
          par_initial              = par_initial,
          max_number_of_iterations = max_number_of_iterations,
          occupancy_min            = occupancy_min,
          occupancy_max            = occupancy_max)
        par_initial = minimized.par_min
        n_iterations_total += minimized.n_iterations
        n_fun_total        += minimized.n_fun
      # Constrained groups are bounded by 0 and 1 in minimization.
      fmodels.fmodel_xray().xray_structure.adjust_occupancy(
        occ_max   = occupancy_max,
        occ_min   = occupancy_min,
        selection = i_selection)
      xray_structure_final = fmodels.fmodel_xray().xray_structure
      model.set_xray_structure(xray_structure_final)
      fmodels.update_xray_structure(xray_structure = xray_structure_final,
                                    update_f_calc  = True)
      refined_occ = xray_structure_final.scatterers().extract_occupancies().\
        select(i_selection)
      assert flex.min(refined_occ) >= occupancy_min
      assert flex.max(refined_occ) <= occupancy_max
      self.show(fmodels= fmodels, log = log, message="occupancy refinement: end")
      import mmtbx.refinement
      mmtbx.refinement.show_number_of_minimizer_iterations(
        n_iterations = n_iterations_total, n_fun = n_fun_total, log = log)

  def show(self, fmodels, message, log):
    if(log is not None):
      print("|-"+message+"-"*(79-len("|-"+message+"|"))+"|", file=log)
      fm_x, fm_n = fmodels.fmodel_xray(), fmodels.fmodel_neutron()
      if(fm_n is not None):
        print("|"+" "*36+"X-ray"+" "*36+"|", file=log)
      self.show_helper(fmodel = fm_x, log = log)
      if(fm_n is not None):
        print("|"+" "*35+"neutron"+" "*35+"|", file=log)
        self.show_helper(fmodel = fm_n, log = log)
      occupancies = fm_x.xray_structure.scatterers().extract_occupancies()
      occ_max = format_value("%4.2f", flex.max(occupancies))
      occ_min = format_value("%4.2f", flex.min(occupancies))
      number_small = format_value("%8d", (occupancies < 0.1).count(True))
      print("| occupancies: max = %s  min = %s   number of occupancies < 0.1: %s |"%(
        occ_max, occ_min, number_small), file=log)
      print("|"+"-"*77+"|", file=log)

  def show_helper(self, fmodel, log):
    r_work = format_value("%6.4f", fmodel.r_work())
    r_free = format_value("%6.4f", fmodel.r_free())
    target = format_value("%-13.3f", fmodel.target_w())
    target_name = format_value("%s", fmodel.target_name)
    p1 = "| r_work = %s r_free = %s" % (r_work, r_free)
    p2 = "target_work(%s) = %s |" % (target_name, target)
    print(p1+" "*(79-len(p1+p2))+p2, file=log)

class minimizer(object):
  """
  Minimize the data target with respect to occupancies using L-BFGS-B.

  Atoms of one selection share one occupancy, q. A constrained group of one
  selection has one parameter, x = q, bounded by occupancy_min and
  occupancy_max. Occupancies of a constrained group of n > 1 selections add
  up to one; the group has n-1 parameters, each bounded by 0 and 1:
  q1 = x1, q2 = (1-x1)*x2, ..., qn = (1-x1)*...*(1-x(n-1)). This way box
  bounds on x keep all q of the group between 0 and 1.

  Structure factors are a sum over atoms, so only the contribution of
  selected atoms is recalculated in minimization (not with twinning or joint
  X-ray and neutron refinement). For this sum to be equal to structure
  factors calculated from all atoms at once, sampling of atoms in FFT
  calculations is fixed by the structure the minimization starts with.

  Parameters
  ----------
  fmodels : mmtbx.fmodels
      Occupancy gradient flags of the scatterers must be set for the
      selected atoms.
  selections : list of list of flex.size_t
      Constrained groups, each a list of atom selections.
  par_initial : flex.double
      Starting occupancy for each selection, in the order of selections.
  max_number_of_iterations : int
      Zero means no minimization.
  occupancy_min, occupancy_max : float
      Bounds for occupancies of constrained groups of one selection.

  Attributes
  ----------
  par_min : flex.double
      Refined occupancy for each selection.
  n_iterations, n_fun : int
      Number of iterations and of target evaluations.
  """
  def __init__(self,
               fmodels,
               selections,
               par_initial,
               max_number_of_iterations,
               occupancy_min = 0,
               occupancy_max = 1):
    self.fmodels = fmodels
    self.fmodels.create_target_functors()
    self.fmodels.prepare_target_functors_for_minimization()
    from mmtbx.refinement import weights
    self.weights = weights.weights(wx = 1, wx_scale = 1, w = 0)
    self.group_sizes = []
    self.i_selection = flex.size_t() # all selected atoms
    self.i_q = flex.size_t()         # index of q for each selected atom
    self.q_ranges = []               # range in i_selection for each q
    lower, upper = flex.double(), flex.double()
    for sel in selections:
      self.group_sizes.append(len(sel))
      for sel_ in sel:
        start = self.i_selection.size()
        self.i_selection.extend(sel_)
        self.i_q.extend(flex.size_t(sel_.size(), len(self.q_ranges)))
        self.q_ranges.append((start, self.i_selection.size()))
      if(len(sel) == 1):
        lower.append(occupancy_min)
        upper.append(occupancy_max)
      else:
        lower.extend(flex.double(len(sel)-1, 0))
        upper.extend(flex.double(len(sel)-1, 1))
    fmodel = self.fmodels.fmodel_xray()
    self.occupancies = fmodel.xray_structure.scatterers(
      ).extract_occupancies()
    self.f_calc_fixed = None # from atoms that are not selected
    if(self.fmodels.fmodel_neutron() is None and fmodel.twin_law is None and
       fmodel.sfg_params.algorithm in ["fft", "direct"]):
      import cctbx.xray.structure_factors
      p = fmodel.sfg_params
      self.sf_managers = []
      for manager_type, miller_set in [
            (cctbx.xray.structure_factors.from_scatterers, fmodel.f_obs()),
            (cctbx.xray.structure_factors.gradients, fmodel.f_obs_work())]:
        manager = manager_type(
          miller_set                   = miller_set,
          cos_sin_table                = p.cos_sin_table,
          grid_resolution_factor       = p.grid_resolution_factor,
          quality_factor               = p.quality_factor,
          u_base                       = p.u_base,
          b_base                       = p.b_base,
          wing_cutoff                  = p.wing_cutoff,
          exp_table_one_over_step_size = p.exp_table_one_over_step_size)
        manager.fix_sampling(xray_structure = fmodel.xray_structure)
        self.sf_managers.append(manager)
      xray_structure_fixed = fmodel.xray_structure.deep_copy_scatterers()
      xray_structure_fixed.set_occupancies(0, self.i_selection)
      self.f_calc_fixed = self.f_calc(xray_structure = xray_structure_fixed)
    self.x = self.q_as_x(q = par_initial)
    for i in range(self.x.size()):
      self.x[i] = min(upper[i], max(lower[i], self.x[i]))
    self.x_current = None # parameters the structure factors correspond to
    self.n_fun = 0
    self.n_iterations = 0
    # L-BFGS-B stops on absolute changes of a target that is less than one.
    self.target_scale = None
    self.occupancy_small = 0.01
    if(max_number_of_iterations > 0):
      from scitbx import lbfgsb
      n = self.x.size()
      m = lbfgsb.run(
        target_evaluator = self,
        bound_flags      = flex.int(n, 2),
        lower_bound      = lower,
        upper_bound      = upper,
        n                = n,
        max_iterations   = max_number_of_iterations)
      self.n_iterations = m.n_iteration()
    self.par_min = self.x_as_q(x = self.x)
    if(self.x_current is None or not self.x_current.all_eq(self.x)):
      self.set_occupancies(q = self.par_min)

  def q_as_x(self, q):
    x = flex.double()
    i = 0
    for n in self.group_sizes:
      if(n == 1):
        x.append(q[i])
      else:
        rest = 1.
        for k in range(n-1):
          if(rest > 1.e-6): x.append(min(1., max(0., q[i+k]/rest)))
          else:             x.append(0.5)
          rest -= q[i+k]
      i += n
    return x

  def x_as_q(self, x):
    q = flex.double()
    j = 0
    for n in self.group_sizes:
      if(n == 1):
        q.append(x[j])
        j += 1
      else:
        rest = 1.
        for k in range(n-1):
          q.append(rest*x[j])
          rest *= (1-x[j])
          j += 1
        q.append(rest)
    return q

  def gradients_wrt_x(self, x, gq):
    result = flex.double()
    i, j = 0, 0
    for n in self.group_sizes:
      if(n == 1):
        result.append(gq[i])
      else:
        # v[k]: derivative wrt the part of the group left for selections k,
        # k+1, ..., n-1, per unit of that part.
        v = [0]*n
        v[n-1] = gq[i+n-1]
        for k in range(n-2, -1, -1):
          v[k] = x[j+k]*gq[i+k] + (1-x[j+k])*v[k+1]
        rest = 1.
        for k in range(n-1):
          result.append(rest*(gq[i+k]-v[k+1]))
          rest *= (1-x[j+k])
      i += n
      j += max(1, n-1)
    return result

  def f_calc(self, xray_structure):
    fmodel = self.fmodels.fmodel_xray()
    return self.sf_managers[0](
      xray_structure = xray_structure,
      miller_set     = fmodel.f_obs(),
      algorithm      = fmodel.sfg_params.algorithm).f_calc()

  def set_occupancies(self, q):
    self.occupancies.set_selected(self.i_selection, q.select(self.i_q))
    fmodel = self.fmodels.fmodel_xray()
    fmodel.xray_structure.scatterers().set_occupancies(self.occupancies)
    if(self.f_calc_fixed is None):
      self.fmodels.update_xray_structure(update_f_calc = True)
    else:
      f_calc = self.f_calc(
        xray_structure = fmodel.xray_structure.select(self.i_selection))
      fmodel.update(f_calc = f_calc.customized_copy(
        data = f_calc.data() + self.f_calc_fixed.data()))

  def compute_functional_and_gradients(self):
    q = self.x_as_q(x = self.x)
    self.set_occupancies(q = q)
    self.x_current = self.x.deep_copy()
    # The derivative of a structure factor wrt occupancy does not depend on
    # occupancy, but the FFT calculation gives zero for an atom with zero
    # occupancy: such an atom would stay at the lower bound forever. So
    # gradients are calculated with occupancies that are not too small.
    scatterers = self.fmodels.fmodel_xray().xray_structure.scatterers()
    small = q < self.occupancy_small
    if(small.count(True) > 0):
      q_g = q.deep_copy().set_selected(small, self.occupancy_small)
      scatterers.set_occupancies(self.occupancies.deep_copy().set_selected(
        self.i_selection, q_g.select(self.i_q)))
    if(self.f_calc_fixed is None):
      t_g = self.fmodels.target_and_gradients(
        weights           = self.weights,
        compute_gradients = True,
        occupancy         = True)
      f = t_g.target()
      g = t_g.gradients().select(self.i_selection)
    else:
      t_g = self.fmodels.target_functor_result_xray(compute_gradients = True)
      f = t_g.target_work()
      fmodel = self.fmodels.fmodel_xray()
      d_target_d_f_calc = t_g.d_target_d_f_calc_work()
      g = self.sf_managers[1](
        xray_structure         = fmodel.xray_structure.select(self.i_selection),
        u_iso_refinable_params = None,
        miller_set             = d_target_d_f_calc,
        d_target_d_f_calc      = d_target_d_f_calc.data(),
        n_parameters           = 0,
        algorithm              = fmodel.sfg_params.algorithm
        ).d_target_d_occupancy()
    scatterers.set_occupancies(self.occupancies)
    gq = flex.double()
    for start, end in self.q_ranges:
      gq.append(flex.sum(g[start:end]))
    self.n_fun += 1
    if(self.target_scale is None):
      self.target_scale = 1
      if(f > 0): self.target_scale = 1/f
    gx = self.gradients_wrt_x(x = self.x, gq = gq)
    return self.x, f*self.target_scale, gx*self.target_scale

#-----------------------------------------------------------------------
# SELECTION HANDLING

def list_3d_as_bool_selection(list_3d, size):
  result = flex.bool(size, False)
  for i in list_3d:
    for j in i:
      for k in j:
        if (result[k]):
          raise Sorry("Duplicate selection for occupancies.")
        result[k] = True
  return result

def add_occupancy_selection(result, size, selection, hd_special=None):
  result_as_1d_array_b = list_3d_as_bool_selection(list_3d=result, size=size)
  sel_b = selection
  if(isinstance(selection, flex.size_t)):
    sel_b = flex.bool(size, selection)
  if(hd_special is not None):
    not_common = ((sel_b != result_as_1d_array_b) & (sel_b == True))
    not_common_ = ((not_common != hd_special) & (not_common == True)).iselection()
    not_common = not_common_
  else:
    not_common = ((sel_b != result_as_1d_array_b) & (sel_b == True)).iselection()
  sel_checked = []
  for i in not_common:
    sel_checked.append([[i]])
  if(len(sel_checked) > 0):
    result.extend(sel_checked)
  return result

def extract_partial_occupancy_selections(hierarchy):
  result = []
  for model in hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        if(not residue_group.have_conformers()):
          assert len(residue_group.atom_groups()) == 1
          occs = flex.double()
          i_seqs = []
          for atom in residue_group.atoms():
            occs.append(atom.occ)
            i_seqs.append(atom.i_seq)
          if(occs[0]<1 and occs[0]!=0 and occs.all_eq(occs[0]) and occs.size()>1):
            result.append([i_seqs])
  return result

def occupancy_selections_wrapper(
      model,
      add_water                          = False,
      other_individual_selection_strings = None,
      other_constrained_groups           = None,
      remove_selection                   = None,
      as_flex_arrays                     = True,
      constrain_correlated_3d_groups     = False,
      log                                = None,
      ungrouped                          = False):
  result = occupancy_selections(
    model                              = model,
    add_water                          = add_water,
    other_individual_selection_strings = other_individual_selection_strings,
    other_constrained_groups           = other_constrained_groups,
    remove_selection                   = remove_selection,
    as_flex_arrays                     = as_flex_arrays,
    constrain_correlated_3d_groups     = constrain_correlated_3d_groups,
    log                                = log)
  if ungrouped:
    result = [[v] for sub in result for v in sub]
  return result

def occupancy_selections(
      model,
      add_water                          = False,
      other_individual_selection_strings = None,
      other_constrained_groups           = None,
      remove_selection                   = None,
      as_flex_arrays                     = True,
      constrain_correlated_3d_groups     = False,
      log                                = None):
  # set up defaults
  if(other_individual_selection_strings is not None and
     len(other_individual_selection_strings) == 0):
    other_individual_selection_strings = None
  if(other_constrained_groups is not None and
     len(other_constrained_groups) == 0):
    other_constrained_groups = None
  if(remove_selection is not None and len(remove_selection) == 0):
    remove_selection = None
  result = model.get_hierarchy().occupancy_groups_simple(
    common_residue_name_class_only = None,
    always_group_adjacent          = False,
    ignore_hydrogens               = False)
  exchangable_hd_pairs = model.get_hierarchy().exchangeable_hd_selections()
  if(len(exchangable_hd_pairs)==0 and result is not None):
    occupancy_regroupping(
      pdb_hierarchy = model.get_hierarchy(),
      cgs           = result)
  result = mmtbx.utils.remove_selections(selection = result,
    other = exchangable_hd_pairs,
    size = model.get_number_of_atoms())
  result.extend(exchangable_hd_pairs)
  # extract group-[0,1]-constrained atoms withing a residue
  pogl = extract_partial_occupancy_selections(hierarchy = model.get_hierarchy())
  rm_duplicate_with_pogl = []
  for t_ in pogl:
    for t__ in t_:
      for t___ in t__:
        rm_duplicate_with_pogl.append(t___)
  result = mmtbx.utils.remove_selections(selection = result, other = pogl,
    size = model.get_number_of_atoms())
  result.extend(pogl)
  # add partial occupancies
  occupancies = model.get_xray_structure().scatterers().extract_occupancies()
  sel = (occupancies != 1.) & (occupancies != 0.)
  result = add_occupancy_selection(
    result     = result,
    size       = model.get_number_of_atoms(),
    selection  = sel,
    hd_special = None)
  # check user's input
  all_sel_strgs = []
  if(other_individual_selection_strings is not None):
    all_sel_strgs = all_sel_strgs + other_individual_selection_strings
  if(other_constrained_groups is not None):
    for other_constrained_group in other_constrained_groups:
      for other_constrained_group in other_constrained_groups:
        if(len(other_constrained_group.selection)>0):
          all_sel_strgs = all_sel_strgs + other_constrained_group.selection
  if(len(all_sel_strgs) > 0):
    for sel_str in all_sel_strgs:
      sel_str_sel = get_atom_selections(
        model               = model,
        selection_strings   = [sel_str],
        iselection          = True,
        one_selection_array = True)
      if(sel_str_sel.size() == 0):
        raise Sorry("Empty selection: %s"%sel_str)
  #
  if([other_individual_selection_strings,
      other_constrained_groups].count(None) == 0):
    sel1 = get_atom_selections(
      model               = model,
      selection_strings   = other_individual_selection_strings,
      iselection          = True,
      one_selection_array = True)
    for other_constrained_group in other_constrained_groups:
      for other_constrained_group in other_constrained_groups:
        for cg_sel_strs in other_constrained_group.selection:
          sel2 = get_atom_selections(
            model               = model,
            selection_strings   = cg_sel_strs,
            iselection          = True,
            one_selection_array = True)
          if(sel1.intersection(sel2).size() > 0):
            raise Sorry("Duplicate selection: same atoms selected for individual and group occupancy refinement.")
  # check user's input and apply remove_selection to default selection
  if(remove_selection is not None):
    sel1 = get_atom_selections(
      model               = model,
      selection_strings   = remove_selection,
      iselection          = True,
      one_selection_array = True)
    if(sel1.size() == 0): # XXX check all and not total.
      raise Sorry("Empty selection: remove_selection.")
    if(other_individual_selection_strings is not None):
      sel2 = get_atom_selections(
        model               = model,
        selection_strings   = other_individual_selection_strings,
        iselection          = True,
        one_selection_array = True)
      if(sel1.intersection(sel2).size() > 0):
        raise Sorry("Duplicate selection: occupancies of same atoms selected to be fixed and to be refined.")
    if(other_constrained_groups is not None):
      for other_constrained_group in other_constrained_groups:
        for cg_sel_strs in other_constrained_group.selection:
          sel2 = get_atom_selections(
            model               = model,
            selection_strings   = cg_sel_strs,
            iselection          = True,
            one_selection_array = True)
          if(sel1.intersection(sel2).size() > 0):
            raise Sorry("Duplicate selection: occupancies of same atoms selected to be fixed and to be refined.")
    result = mmtbx.utils.remove_selections(selection = result, other = sel1,
      size = model.get_number_of_atoms())
  #
  if(other_individual_selection_strings is not None):
    sel = get_atom_selections(
      model               = model,
      selection_strings   = other_individual_selection_strings,
      iselection          = True,
      one_selection_array = True)
    result = mmtbx.utils.remove_selections(selection = result, other = sel,
      size = model.get_number_of_atoms())
    result = add_occupancy_selection(
      result     = result,
      size       = model.get_number_of_atoms(),
      selection  = sel,
      hd_special = None)
  if(other_constrained_groups is not None):
    for other_constrained_group in other_constrained_groups:
      cg_sel = []
      for cg_sel_strs in other_constrained_group.selection:
        sel = get_atom_selections(
          model               = model,
          selection_strings   = cg_sel_strs,
          iselection          = True,
          one_selection_array = True)
        result = mmtbx.utils.remove_selections(selection = result, other = sel,
          size = model.get_number_of_atoms())
        if(sel.size() > 0):
          cg_sel.append(list(sel))
      if(len(cg_sel) > 0):
        result.append(cg_sel)
  # capping H (occupancy 0, not part of the deposited model) are never refined:
  # not as members of their conformer's constrained group, not as partial
  # occupancies, not as members of the user's selections
  result = mmtbx.utils.remove_selections(selection = result,
    other = model.selection('capping_h').iselection(),
    size = model.get_number_of_atoms())

  if (constrain_correlated_3d_groups) and (len(result) > 0):
      result = assemble_constraint_groups_3d(
        xray_structure=model.get_xray_structure(),
        pdb_atoms=model.get_atoms(),
        constraint_groups=result,
        log=log)

  if(add_water):
    water_selection = get_atom_selections(
      model                 = model,
      selection_strings     = ['water'],
      iselection            = True,
      allow_empty_selection = True,
      one_selection_array   = True)
    if(remove_selection is not None):
      sel_rm = get_atom_selections(
        model               = model,
        selection_strings   = remove_selection,
        iselection          = True,
        one_selection_array = True)
      if water_selection is not None:
        water_selection = flex.size_t(
          list(set(water_selection).difference(sel_rm)))
    if water_selection.size()>0:
      def flatten(lst):
        flat_list = []
        for item in lst:
            if isinstance(item, list):
                flat_list.extend(flatten(item))
            else:
                flat_list.append(item)
        return flat_list
      occ_groups_of_more_than_one = []
      for g in result:
        if len(g)>1:
          occ_groups_of_more_than_one.extend( flatten(g) )
      water_selection = list(water_selection)
      wocc = model.get_hierarchy().atoms().extract_occ()
      wsel = model.solvent_selection().iselection()
      for i in wsel:
         # Drop zero-occupancy and multi-conformer solvent atoms.
         # Single IF so no chance to call water_seelction.remove(i)
         # twice on the same water.
         if i in water_selection and (
              wocc[i]<1.e-6 or i in occ_groups_of_more_than_one):
           water_selection.remove(i)
      result = add_occupancy_selection(
        result     = result,
        size       = model.get_number_of_atoms(),
        selection  = flex.size_t(water_selection),
        hd_special = None)

  list_3d_as_bool_selection(
    list_3d=result, size=model.get_number_of_atoms())

  if(len(result) == 0): result = None
  if(as_flex_arrays and result is not None):
    result_ = []
    for gsel in result:
      result__ = []
      for sel in gsel:
        result__.append(flex.size_t(sel))
      result_.append(result__)
    result = result_

  return result

def occupancy_regroupping(pdb_hierarchy, cgs):
  h = pdb_hierarchy
  awl = list(h.atoms_with_labels())
  elements = h.atoms().extract_element()
  rgs = list(h.residue_groups())
  for cg in cgs: # loop over constraint groups
    for c in cg: # loop over conformers of constrained group
      altloc_h = awl[c[0]].altloc
      ind=None
      for ind_ in c:
        if(elements[ind_].strip().upper() in ["H","D"] and
           awl[ind_].name.strip().upper() in ["H","D"] and
           altloc_h.strip() != ""):
          ind = ind_
          break
      if(ind is not None):
        # find "-1" (previous to given constraint group) residue group
        rg_prev = None
        for i_rg, rg in enumerate(rgs):
          if(i_rg > 0 and ind in rg.atoms().extract_i_seq()):
            rg_prev = rgs[i_rg-1]
            break
        if(rg_prev is None): continue
        rg_prev_i_seqs = rg_prev.atoms().extract_i_seq()
        rg_i_seqs = rg.atoms().extract_i_seq()
        # find i_seq of C of previous residue
        i_seqs_c_prev=[]
        for a_prev in rg_prev.atoms():
          if(a_prev.name.strip().upper()=="C"):
            i_seqs_c_prev.append(a_prev.i_seq)
        if(i_seqs_c_prev == []): continue
        # find constarint group corresponding to rg_prev
        cg_prev=None
        for cg2 in cgs:
          for c2 in cg2:
            i_seqs_inter = set(c2).intersection(set(rg_prev_i_seqs))
            if(len(i_seqs_inter)>1 or
               (len(i_seqs_inter)==1 and not
                awl[list(i_seqs_inter)[0]].name.strip().upper() in ["H","D"])):
              for cg2_ in cg2:
                for i_seq_c_prev in i_seqs_c_prev:
                  if(i_seq_c_prev in cg2_):
                    cg_prev = cg2
                    break
          if(cg_prev is not None): break
        if(cg_prev is None): continue
        # identify to which constraint group H belongs and move it there
        for cg2 in cgs:
          for c2 in cg2:
            if(ind in c2):
              c2.remove(ind)
              break
        found = False
        for conformer_prev in rg_prev.conformers():
          if(conformer_prev.altloc == altloc_h):
            conformer_prev_i_seqs = conformer_prev.atoms().extract_i_seq()
            for cg2 in cgs:
              for c2 in cg2:
                i_seqs_inter = set(c2).intersection(set(conformer_prev_i_seqs))
                if(len(i_seqs_inter)>1 or
                   (len(i_seqs_inter)==1 and not
                    awl[list(i_seqs_inter)[0]].name.strip().upper() in ["H","D"])):
                  for i_seq_c_prev in i_seqs_c_prev:
                    if(i_seq_c_prev in c2):
                      c2.append(ind)
                      found = True
  for cg in cgs:
    while cg.count([None])>0:
      cg.remove([None])
    while cg.count([])>0:
      cg.remove([])
  while cgs.count([])>0:
    cgs.remove([])

def assemble_constraint_groups_3d(
    xray_structure,
    pdb_atoms,
    constraint_groups,
    interaction_distance_cutoff=4.0,
    verbose=False,
    log=None):
  """
  Re-sorts occupancy constraint groups so that conformers whose motion is
  correlated (i.e. they interact in 3D, without necessarily being part of
  the same fragment/molecule/ASU) are grouped together.  As input, it expects
  the constraint groups output by mmtbx.utils.occupancy_selections(), which
  will already have connectivity taken into account.  This function will exit
  with an error if the occupancies for the new groups are not consistent.
  """
  if (log is None):
    log = null_out()
  make_sub_header("Correlated occupancy grouping", out=log)
  print("""
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                  WARNING - EXPERIMENTAL FEATURE                        !!
  !!                                                                        !!
  !! Grouping of occupancy constraints in 3D is experimental and not fully  !!
  !! tested.  Use at your own risk!  For bug reports, etc. contact us by    !!
  !! email at help@phenix-online.org.                                       !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
""", file=log)
  occupancies = pdb_atoms.extract_occ()
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=interaction_distance_cutoff)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  k = 0
  n_groups_start = len(constraint_groups)
  while (k < len(constraint_groups)):
    groups = constraint_groups[k]
    print("Constraint group %d: %d conformers" % (k+1, len(groups)), file=log)
    merge_constraints = []
    for i_sel, selection in enumerate(groups):
      occ = occupancies.select(selection)
      altloc = pdb_atoms[selection[0]].fetch_labels().altloc
      print("  conformer '%s': %d atoms" % (altloc, len(selection)), file=log)
      if (not occ.all_eq(occ[0])):
        raise Sorry("At least one occupancy constraint group has "+
          "inconsistent occupancies for atoms in a single conformer.  To use "+
          "the automatic 3D constraints, the starting occupancies must be "+
          "uniform within each selection.")
      for i_seq in selection :
        labels = pdb_atoms[i_seq].fetch_labels()
        if (labels.altloc.strip() == ''):
          continue
        pair_sym_dict = pair_sym_table[i_seq]
        if (verbose):
          print("%s (group %d):" % (pdb_atoms[i_seq].id_str(), k+1))
        for j_seq, sym_ops in pair_sym_dict.items():
          kk = k + 1
          while (kk < len(constraint_groups)):
            combine_group = False
            for other_selection in constraint_groups[kk] :
              if (j_seq in other_selection):
                if (verbose):
                  print("  %s (group %d)" % (pdb_atoms[j_seq].id_str(), kk+1))
                merge_constraints.append(constraint_groups[kk])
                del constraint_groups[kk]
                combine_group = True
                break
            if (not combine_group):
              kk += 1
    if (len(merge_constraints) > 0):
      print("Merging %d constraint groups with group %d" % (
        len(merge_constraints), (k+1)), file=log)
      for selection in groups :
        first_atom = pdb_atoms[selection[0]]
        altloc = first_atom.fetch_labels().altloc
        if (altloc.strip() == ''):
          raise RuntimeError(("Atom '%s' in occupancy constraint group has "+
            "blank altloc ID") % first_atom.id_str())
        for merge_groups in merge_constraints :
          kk = 0
          while (kk < len(merge_groups)):
            other_selection = merge_groups[kk]
            altloc2 = pdb_atoms[other_selection[0]].fetch_labels().altloc
            if (altloc2 == altloc):
              print("  combining %d atoms with altloc %s" % \
                (len(other_selection), altloc), file=log)
              occ1 = occupancies.select(selection)
              occ2 = occupancies.select(other_selection)
              if (not occ1.all_eq(occ2[0])) or (not occ2.all_eq(occ1[0])):
                raise Sorry(
                  ("Inconsistent occupancies in spatially related groups "+
                  "(%.2f versus %.2f).  To use automatic 3D occupancy "+
                  "restraints, the correlated conformers must start with "+
                  "the same initial occupancy.") % (occ1[0], occ2[0]))
              selection.extend(other_selection)
              del merge_groups[kk]
            else :
              kk += 1
      for merge_groups in merge_constraints :
        if (len(merge_groups) > 0):
          for other_selection in merge_groups :
            altloc = pdb_atoms[other_selection[0]].fetch_labels().altloc
            print(("  warning: %d atoms with altloc %s do not "+
              "correspond to an existing group") % (len(other_selection),
              altloc), file=log)
            groups.append(other_selection)
    k += 1
  if (len(constraint_groups) != n_groups_start):
    print("New occupancy constraint groups:", file=log)
    for i_group, constraint_group in enumerate(constraint_groups):
      print("  group %d:" % (i_group+1), file=log)
      for selection in constraint_group :
        resids = []
        altlocs = set()
        for i_seq in selection :
          atom_group = pdb_atoms[i_seq].parent()
          ag_id = atom_group.id_str()
          altlocs.add(atom_group.altloc)
          if (not ag_id in resids):
            resids.append(ag_id)
        assert len(altlocs) == 1
        print("    conformer '%s' (%d atoms):" % (list(altlocs)[0],
          len(selection)), file=log)
        for ag_id in resids :
          print("      atom_group %s" % ag_id, file=log)
  else :
    print("Occupancy constraint groups unmodified.", file=log)
  print("", file=log)
  return constraint_groups
