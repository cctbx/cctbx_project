from __future__ import division
from cctbx.array_family import flex
from scitbx import lbfgs
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args
import sys

class minimizer(object):

  def __init__(self,
        fmodel,
        groups,
        call_back_after_minimizer_cycle=None,
        number_of_minimizer_cycles=3,
        lbfgs_max_iterations=20,
        number_of_finite_difference_tests=0):
    adopt_init_args(self, locals())
    self.x = flex.double()
    for group in groups:
      if (group.refine_f_prime): self.x.append(group.f_prime)
      if (group.refine_f_double_prime): self.x.append(group.f_double_prime)
    fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    for group in groups:
      if (group.refine_f_prime):
        fmodel.xray_structure.scatterers().flags_set_grad_fp(
          iselection=group.iselection)
      if (group.refine_f_double_prime):
        fmodel.xray_structure.scatterers().flags_set_grad_fdp(
          iselection=group.iselection)
    self.target_functor = fmodel.target_functor()
    self.target_functor.prepare_for_minimization()
    for self.i_cycle in xrange(number_of_minimizer_cycles):
      self.lbfgs = lbfgs.run(
        target_evaluator=self,
        termination_params=lbfgs.termination_parameters(
          max_iterations=lbfgs_max_iterations),
        exception_handling_params=lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = True))
      if (call_back_after_minimizer_cycle is not None):
        self.unpack()
        if (not call_back_after_minimizer_cycle(minimizer=self)):
          break
    if (call_back_after_minimizer_cycle is None):
      self.unpack()
    del self.i_cycle
    del self.lbfgs
    del self.x
    del self.target_functor
    del self.fmodel
    del self.groups

  def unpack(self):
    xi = iter(self.x)
    for group in self.groups:
      if (group.refine_f_prime): group.f_prime = xi.next()
      if (group.refine_f_double_prime): group.f_double_prime = xi.next()
    for group in self.groups:
      group.copy_to_scatterers_in_place(
        scatterers=self.fmodel.xray_structure.scatterers())
    self.fmodel.update_xray_structure(update_f_calc=True)

  def compute_functional_and_gradients(self):
    self.unpack()
    t_r = self.target_functor(compute_gradients=True)
    fmodel = self.fmodel
    f = t_r.target_work()
    d_target_d_f_calc = t_r.d_target_d_f_calc_work()
    sfg = fmodel.structure_factor_gradients_w(
      u_iso_refinable_params=None,
      d_target_d_f_calc=d_target_d_f_calc.data(),
      xray_structure=fmodel.xray_structure,
      n_parameters=0,
      miller_set=d_target_d_f_calc,
      algorithm=fmodel.sfg_params.algorithm)
    d_t_d_fp = sfg.d_target_d_fp()
    d_t_d_fdp = sfg.d_target_d_fdp()
    del sfg
    g = flex.double()
    for group in self.groups:
      if (group.refine_f_prime):
        g.append(flex.sum(d_t_d_fp.select(group.iselection)))
      if (group.refine_f_double_prime):
        g.append(flex.sum(d_t_d_fdp.select(group.iselection)))
    if (self.number_of_finite_difference_tests != 0):
      self.number_of_finite_difference_tests -= 1
      g_fin = []
      eps = 1.e-5
      x = self.x
      for i in xrange(x.size()):
        fs = []
        xi0 = x[i]
        for signed_eps in [eps,-eps]:
          x[i] = xi0 + signed_eps
          self.unpack()
          x[i] = xi0
          t_r = self.target_functor(compute_gradients=False)
          fs.append(t_r.target_work())
        g_fin.append((fs[0]-fs[1])/(2*eps))
      self.unpack()
      assert approx_equal(g_fin, g)
    return f, g

def get_single_atom_selection_string (atom) :
  labels = atom.fetch_labels()
  altloc = labels.altloc
  if (altloc == '') : altloc = ' ' # XXX this is gross
  sele = \
    "chain '%s' and resname %s and name '%s' and altloc '%s' and resid %s" % \
      (labels.chain_id, labels.resname, labels.name, altloc, labels.resid())
  return sele

def find_anomalous_scatterer_groups (
    pdb_atoms,
    xray_structure,
    group_same_element=True, # XXX should this be True by default?
    out=None) :
  """
  Automatic setup of anomalously scattering atoms, defined here as anything
  with atomic number 15 (P) or greater.  Not yet accessible from phenix.refine.
  """
  from cctbx.eltbx import sasaki
  from cctbx import xray
  if (out is None) : out = sys.stdout
  element_i_seqs = {}
  groups = []
  if (out is None) : out = null_out()
  for i_seq, scatterer in enumerate(xray_structure.scatterers()) :
    element = scatterer.element_symbol()
    atomic_number = sasaki.table(element).atomic_number()
    if (atomic_number >= 15) :
      if (group_same_element) :
        if (not element in element_i_seqs) :
          element_i_seqs[element] = flex.size_t()
        element_i_seqs[element].append(i_seq)
      else :
        print >> out, "  creating anomalous group for %s" % \
          pdb_atoms[i_seq].id_str()
        asg = xray.anomalous_scatterer_group(
          iselection=flex.size_t([i_seq]),
          f_prime=0,
          f_double_prime=0,
          refine=["f_prime","f_double_prime"],
          selection_string=get_single_atom_selection_string(pdb_atoms[i_seq]))
        groups.append(asg)
  if (group_same_element) :
    for elem in sorted(element_i_seqs.keys()) :
      iselection = element_i_seqs[elem]
      print >> out, \
        "  creating anomalous group for element %s with %d atoms" % \
        (elem, len(iselection))
      asg = xray.anomalous_scatterer_group(
        iselection=iselection,
        f_prime=0,
        f_double_prime=0,
        refine=["f_prime","f_double_prime"],
        selection_string="element %s" % elem)
      groups.append(asg)
  return groups
