
from __future__ import absolute_import, division, print_function
from mmtbx.utils import get_atom_selections
import mmtbx.utils
from cctbx.array_family import flex
from scitbx import lbfgs
from libtbx.str_utils import format_value, make_sub_header
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args
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
      for macro_cycle in range(number_of_macro_cycles):
        if(minimized is not None): par_initial = minimized.par_min
        minimized = minimizer(
          fmodels                       = fmodels,
          selections                    = selections,
          constrained_groups_selections = constrained_groups_selections,
          par_initial                   = par_initial,
          max_number_of_iterations      = max_number_of_iterations)
        if(minimized is not None): par_initial = minimized.par_min
        set_refinable_parameters(
          xray_structure     = fmodels.fmodel_xray().xray_structure,
          parameters         = par_initial,
          selections         = selections,
          enforce_positivity = (occupancy_min>=0))
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
  def __init__(self,
               fmodels,
               constrained_groups_selections,
               selections,
               par_initial,
               max_number_of_iterations):
    adopt_init_args(self, locals())
    self.fmodels.create_target_functors()
    self.fmodels.prepare_target_functors_for_minimization()
    from mmtbx.refinement import weights
    self.weights = weights.weights(wx = 1, wx_scale = 1, w = 0)
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

def set_refinable_parameters(xray_structure, parameters, selections,
                             enforce_positivity=False):
  # XXX PVA: Code below is terribly inefficient and MUST be moved into C++
  sz = xray_structure.scatterers().size()
  i = 0
  for sel in selections:
    # pre-check for positivity begin
    # spread negative occupancies across i_seqs having positive ones
    par_all = flex.double()
    par_neg = flex.double()
    i_p = i
    for sel_ in sel:
      p = parameters[i_p]
      par_all.append(p)
      if(p<0): par_neg.append(p)
      i_p += 1
    if(enforce_positivity and par_neg.size()>0):
      par_all = par_all - flex.min(par_all)
      fs = flex.sum(par_all)
      if(fs != 0):
        par_all = par_all / fs
    # pre-check for positivity end
    for j, sel_ in enumerate(sel):
      sel__b = flex.bool(sz, flex.size_t(sel_))
      xray_structure.set_occupancies(par_all[j], sel__b)
      i+=1

def pack_unpack(x, table):
  result = x.deep_copy()
  for indices in table :
    if (len(indices) == 1):
      i0 = indices[0]
      result[i0] = x[i0]
    else :
      xsum = 0
      for i in indices[0:-1] :
        result[i] = x[i]
        xsum += x[i]
      result[indices[-1]] = 1. - xsum
  return result

def pack_gradients(x, table):
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
      wremove = []
      for i in wsel:
         if wocc[i]<1.e-6:
           water_selection.remove(i)
         if i in occ_groups_of_more_than_one:
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
