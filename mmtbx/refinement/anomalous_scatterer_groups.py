from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from scitbx import lbfgs
from libtbx.str_utils import make_sub_header
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args
import time
import sys
from six.moves import zip
from six.moves import range

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
    for self.i_cycle in range(number_of_minimizer_cycles):
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
      if (group.refine_f_prime): group.f_prime = next(xi)
      if (group.refine_f_double_prime): group.f_double_prime = next(xi)
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
      algorithm=fmodel.sfg_params.algorithm,
      extra_params=fmodel.sfg_params.extra)
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
      for i in range(x.size()):
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

def get_single_atom_selection_string(atom):
  labels = atom.fetch_labels()
  altloc = labels.altloc
  if (altloc == '') : altloc = ' ' # XXX this is gross
  sele = \
    "chain '%s' and resname %s and name '%s' and altloc '%s' and resid %s" % \
      (labels.chain_id, labels.resname, labels.name, altloc, labels.resid())
  return sele

def find_anomalous_scatterer_groups(
    pdb_atoms,
    xray_structure,
    group_same_element=True, # XXX should this be True by default?
    out=None):
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
  hd_selection = xray_structure.hd_selection()
  for i_seq, scatterer in enumerate(xray_structure.scatterers()):
    if (hd_selection[i_seq]):
      continue
    element = scatterer.element_symbol().strip()
    try :
      atomic_number = sasaki.table(element).atomic_number()
    except RuntimeError as e :
      print("Error for %s" % pdb_atoms[i_seq].id_str(), file=out)
      print("  " + str(e), file=out)
      continue
    if (atomic_number >= 15):
      if (group_same_element):
        if (not element in element_i_seqs):
          element_i_seqs[element] = flex.size_t()
        element_i_seqs[element].append(i_seq)
      else :
        print("  creating anomalous group for %s" % \
          pdb_atoms[i_seq].id_str(), file=out)
        asg = xray.anomalous_scatterer_group(
          iselection=flex.size_t([i_seq]),
          f_prime=0,
          f_double_prime=0,
          refine=["f_prime","f_double_prime"],
          selection_string=get_single_atom_selection_string(pdb_atoms[i_seq]))
        groups.append(asg)
  if (group_same_element):
    for elem in sorted(element_i_seqs.keys()):
      iselection = element_i_seqs[elem]
      print("  creating anomalous group for element %s with %d atoms" % \
        (elem, len(iselection)), file=out)
      asg = xray.anomalous_scatterer_group(
        iselection=iselection,
        f_prime=0,
        f_double_prime=0,
        refine=["f_prime","f_double_prime"],
        selection_string="element %s" % elem)
      groups.append(asg)
  return groups

def refine_anomalous_substructure(
    fmodel,
    pdb_hierarchy,
    wavelength=None,
    map_type="anom_residual",
    exclude_waters=False,
    exclude_non_water_light_elements=True,
    n_cycles_max=None,
    map_sigma_min=3.0,
    refine=("f_prime","f_double_prime"),
    reset_water_u_iso=True,
    use_all_anomalous=True,
    verbose=True,
    out=sys.stdout):
  """
  Crude mimic of Phaser's substructure completion, with two essential
  differences: only the existing real scatterers in the input model will be
  used (with the assumption that the model is already more or less complete),
  and the anomalous refinement will be performed in Phenix, yielding both
  f-prime and f-double-prime.  The refined f-prime provides us with an
  orthogonal estimate of the number of electrons missing from an incorrectly
  labeled scatterer.

  :param wavelength: X-ray wavelenth in Angstroms
  :param exclude_waters: Don't refine anomalous scattering for water oxygens
  :param exclude_non_water_light_elements: Don't refine anomalous scattering
    for light atoms other than water (CHNO).
  :param n_cycles_max: Maximum number of refinement cycles
  :param map_sigma_min: Sigma cutoff for identify anomalous scatterers
  :param reset_water_u_iso: Reset B-factors for water atoms prior to f'
    refinement
  :param use_all_anomalous: include any scatterers which are already modeled
    as anomalous in the refinement
  """
  from cctbx import xray
  assert (fmodel.f_obs().anomalous_flag())
  assert (map_type in ["llg", "anom_residual"])
  make_sub_header("Iterative anomalous substructure refinement", out=out)
  fmodel.update(target_name="ls")
  pdb_atoms = pdb_hierarchy.atoms()
  non_water_non_hd_selection = pdb_hierarchy.atom_selection_cache().selection(
    "(not element H and not element D and not resname HOH)")
  sites_frac = fmodel.xray_structure.sites_frac()
  scatterers = fmodel.xray_structure.scatterers()
  u_iso_mean = flex.mean(
    fmodel.xray_structure.extract_u_iso_or_u_equiv().select(
      non_water_non_hd_selection))
  anomalous_iselection = flex.size_t()
  anomalous_groups = []
  t_start = time.time()
  n_cycle = 0
  while ((n_cycles_max is None) or (n_cycle < n_cycles_max)):
    n_cycle += 1
    n_new_groups = 0
    t_start_cycle = time.time()
    print("Cycle %d" % n_cycle, file=out)
    anom_map = fmodel.map_coefficients(map_type=map_type).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
    map_min = abs(flex.min(anom_map.as_1d()))
    map_max = flex.max(anom_map.as_1d())
    print("  map range: -%.2f sigma to %.2f sigma" % (map_min, map_max), file=out)
    reset_u_iso_selection = flex.size_t()
    for i_seq, atom in enumerate(pdb_atoms):
      resname = atom.parent().resname
      elem = atom.element.strip()
      if  ((i_seq in anomalous_iselection) or
           ((exclude_waters) and (resname == "HOH")) or
           ((elem in ["H","D","N","C","O"]) and (resname != "HOH") and
            exclude_non_water_light_elements)):
        continue
      scatterer = scatterers[i_seq]
      site_frac = sites_frac[i_seq]
      anom_map_value = anom_map.tricubic_interpolation(site_frac)
      if ((anom_map_value >= map_sigma_min) or
          ((scatterer.fdp != 0) and use_all_anomalous)):
        if (verbose):
          if (n_new_groups == 0):
            print("", file=out)
            print("  new anomalous scatterers:", file=out)
          print("    %-34s  map height: %6.2f sigma" % (atom.id_str(),
            anom_map_value), file=out)
        anomalous_iselection.append(i_seq)
        selection_string = get_single_atom_selection_string(atom)
        group = xray.anomalous_scatterer_group(
          iselection=flex.size_t([i_seq]),
          f_prime=0,
          f_double_prime=0,
          refine=list(refine),
          selection_string=selection_string)
        anomalous_groups.append(group)
        n_new_groups += 1
        if (resname == "HOH") and (reset_water_u_iso):
          water_u_iso = scatterer.u_iso
          if (water_u_iso < u_iso_mean):
            reset_u_iso_selection.append(i_seq)
    if (n_new_groups == 0):
      print("", file=out)
      print("No new groups - anomalous scatterer search terminated.", file=out)
      break
    elif (not verbose):
      print("  %d new groups" % n_new_groups, file=out)
    for i_seq in anomalous_iselection :
      sc = scatterers[i_seq]
      sc.fp = 0
      sc.fdp = 0
    if (verbose):
      print("", file=out)
      print("Anomalous refinement:", file=out)
      fmodel.info().show_targets(text="before minimization", out=out)
      print("", file=out)
    u_iso = fmodel.xray_structure.extract_u_iso_or_u_equiv()
    u_iso.set_selected(reset_u_iso_selection, u_iso_mean)
    fmodel.xray_structure.set_u_iso(values=u_iso)
    fmodel.update_xray_structure(update_f_calc=True)
    minimizer(fmodel=fmodel, groups=anomalous_groups)
    if (verbose):
      fmodel.info().show_targets(text="after minimization", out=out)
      print("", file=out)
      print("  Refined sites:", file=out)
      for i_seq, group in zip(anomalous_iselection, anomalous_groups):
        print("    %-34s  f' = %6.3f  f'' = %6.3f" % (
          pdb_atoms[i_seq].id_str(), group.f_prime, group.f_double_prime), file=out)
    t_end_cycle = time.time()
    print("", file=out)
    if (verbose):
      print("  time for this cycle: %.1fs" % (t_end_cycle-t_start_cycle), file=out)
  fmodel.update(target_name="ml")
  print("%d anomalous scatterer groups refined" % len(anomalous_groups), file=out)
  t_end = time.time()
  print("overall time: %.1fs" % (t_end - t_start), file=out)
  return anomalous_groups
