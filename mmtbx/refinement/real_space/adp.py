from __future__ import division
import mmtbx.refinement.real_space.utils
import mmtbx.refinement.utils
from scitbx.array_family import flex
from libtbx import adopt_init_args
from cctbx import adptbx
from libtbx import easy_mp

from libtbx.test_utils import approx_equal
import mmtbx.secondary_structure
from libtbx.utils import null_out
from itertools import groupby
from operator import itemgetter
import iotbx.phil

def get_chain_selections(pdb_hierarchy):
  pdb_hierarchy.atoms().reset_i_seq()
  params = mmtbx.secondary_structure.sec_str_master_phil.extract()
  params.secondary_structure.protein.search_method="from_ca"
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy                = pdb_hierarchy,
    sec_str_from_pdb_file        = None,
    params                       = params.secondary_structure,
    log                          = null_out())
  n_atoms = pdb_hierarchy.atoms().size()
  helix_sel = [s.iselection() for s in ssm.helix_selections()]
  beta_sel  = [s.iselection() for s in ssm.beta_selections()]
  bp_sel    = [s.iselection() for s in ssm.base_pair_selections()]
  ss_all = helix_sel + beta_sel + bp_sel
  #
  chain_selections = []
  for sel in ss_all:
    ph = pdb_hierarchy.select(sel)
    isels = ph.atoms().extract_i_seq()
    chain_selections.append(isels)
  #
  ss_all_as_one_array = ss_all[0]
  if(len(ss_all)>1):
    for s in ss_all[1:]:
      ss_all_as_one_array.extend(s)
  ss_all_as_one_array_bool = flex.bool(n_atoms, ss_all_as_one_array)
  #
  iseqs = []
  for atom in pdb_hierarchy.atoms():
    if(not ss_all_as_one_array_bool[atom.i_seq]):
      iseqs.append(atom.i_seq)
  for k, g in groupby(enumerate(iseqs), lambda (i, x): i-x):
    isels = list(map(itemgetter(1), g))
    chain_selections.append(flex.size_t(isels))
  # DEBUG
  C1 = flex.size_t()
  for s in chain_selections:
    C1.extend(s)
  s = flex.sort_permutation(C1)
  C1 = C1.select(s)
  C2 = flex.size_t(xrange(n_atoms))
  asert approx_equal(C1, C2)
  assert C1.size()==C2.size()
  #
  return chain_selections

class real_space_group_adp_refinery_via_reciprocal_space(object):
  def __init__(self,
               target_map,
               pdb_hierarchy,
               atom_radius,
               use_adp_restraints,
               nproc,
               log=None):
    adopt_init_args(self, locals())
    self.xray_structure = self.pdb_hierarchy.extract_xray_structure(
      crystal_symmetry = self.target_map.miller_array.crystal_symmetry())
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)

    self.xray_structure = self.xray_structure.set_b_iso(
      value = flex.mean(b_isos))

    #for rg in self.pdb_hierarchy.residue_groups():
    #  sel = rg.atoms().extract_i_seq()
    #  sel = flex.bool(b_isos.size(), sel)
    #  self.xray_structure = self.xray_structure.set_b_iso(
    #    value     = flex.mean(b_isos.select(sel)),
    #    selection = sel)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    #self.chain_selections = []
    #for chain in self.pdb_hierarchy.chains():
    #  self.chain_selections.append(chain.atoms().extract_i_seq())
    self.chain_selections = get_chain_selections(
      pdb_hierarchy = self.pdb_hierarchy)

  def refine_box_with_selected(self, selection=None):
    if(selection is None): selections = self.xray_structure.all_selection()
    ph_box = self.pdb_hierarchy.select(selection)
    ph_box.atoms().reset_i_seq()
    box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure         = self.xray_structure.select(selection),
      map_data               = self.target_map.map_data,
      box_cushion            = self.atom_radius,
      mask_atoms             = True,
      mask_atoms_atom_radius = self.atom_radius)
    ph_box.adopt_xray_structure(box.xray_structure_box)
    group_adp_sel = []
    
    #group_adp_sel.append(ph_box.atoms().extract_i_seq())
    
    for rg in ph_box.residue_groups():
      group_adp_sel.append(rg.atoms().extract_i_seq())
      
    f_obs_box = abs(box.box_map_coefficients(d_min = self.target_map.d_min))
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs_box,
      xray_structure = box.xray_structure_box)

    #XXX
    fmodel.update_all_scales(update_f_part1=False, apply_back_trace=True,
      remove_outliers=False)
    #XXX

    fmodel.update(target_name="ls_wunit_k1")
    if(self.nproc>1): log = None
    else:             log = self.log
    group_b_manager = mmtbx.refinement.group.manager(
      fmodel                   = fmodel,
      selections               = group_adp_sel,
      convergence_test         = False,
      max_number_of_iterations = 50,
      number_of_macro_cycles   = 3,
      run_finite_differences_test = False,
      use_restraints           = self.use_adp_restraints,
      refine_adp               = True,
      log                      = log)
    return fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)

  def refine(self):
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
    if(self.nproc==1):
      for sel in self.chain_selections:
        b_isos_refined = self.refine_box_with_selected(selection=sel)
        b_isos = b_isos.set_selected(sel, b_isos_refined)
    else:
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.refine_box_with_selected,
        args         = self.chain_selections,
        func_wrapper = "buffer_stdout_stderr")
      for i, it in enumerate(stdout_and_results):
        so, b_isos_refined = it
        b_isos = b_isos.set_selected(self.chain_selections[i], b_isos_refined)
        print >> self.log, so
    self.xray_structure = self.xray_structure.set_b_iso(values = b_isos)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
