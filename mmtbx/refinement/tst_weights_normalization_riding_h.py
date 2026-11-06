import random
import re
import iotbx.pdb
import libtbx
import mmtbx.f_model
import mmtbx.model
import mmtbx.restraints
from cctbx import geometry_restraints
from cctbx.array_family import flex
from io import StringIO
from mmtbx.refinement import refinement_flags, weights
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.refinement.tst_minimization_riding_h import pdb_str

#-----------------------------------------------------------------------------
# The normalized geometry and ADP restraint targets divide by a restraint
# count that excludes riding H atoms, so that a weight means the same thing
# with and without riding H. H atoms the riding manager does not
# parameterize (here the water H) stay in the count; with
# normalization=False nothing changes.
#-----------------------------------------------------------------------------

# A water whose H have no riding parameterization.
water_str = """\
HETATM   24  O   HOH S   1      11.000  13.000   9.500  1.00 10.00           O
HETATM   25  H1  HOH S   1      11.800  13.300   9.500  1.00 10.00           H
HETATM   26  H2  HOH S   1      10.600  13.700   9.800  1.00 10.00           H
"""

def get_models(setup_riding=True, scattering_table="n_gaussian"):
  """
  The ARG residue plus the water, with and without H; the heavy-atom
  coordinates are identical. The nonbonded listing cutoff is fixed: left
  to pdb_interpretation it follows the vdw radii, which are larger when H
  are absent, and the heavy-atom pair lists would differ.
  """
  lines_h = (pdb_str + water_str).splitlines()
  lines_no_h = [l for l in lines_h if not l.endswith(" H")]
  result = []
  for lines in [lines_no_h, lines_h]:
    pdb_inp = iotbx.pdb.input(lines=lines, source_info=None)
    model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.nonbonded_distance_cutoff = 4.0
    model.process(pdb_interpretation_params=params, make_restraints=True,
      grm_normalization=True)
    model.setup_scattering_dictionaries(scattering_table=scattering_table)
    if(model.has_hd() and setup_riding):
      model.setup_riding_h_manager(idealize=True)
    result.append(model)
  model_no_h, model_h = result
  heavy = ~model_h.get_hd_selection()
  assert approx_equal(
    model_h.get_sites_cart().select(heavy), model_no_h.get_sites_cart())
  return model_no_h, model_h

def n_without(index_lists, flagged):
  """
  Entries none of whose indices is flagged; a plain Python count to check
  the C++ against.
  """
  result = 0
  for i_seqs in index_lists:
    if (not any(flagged[i] for i in i_seqs)): result += 1
  return result

def exercise_geometry_count():
  model_no_h, model_h = get_models()
  rm_h = model_h.get_restraints_manager()
  assert rm_h.normalization
  e = model_no_h.restraints_manager_energies_sites(compute_gradients=True)
  e_h = model_h.restraints_manager_energies_sites(compute_gradients=True)
  g, g_h = e.geometry, e_h.geometry
  riding = model_h.get_riding_h_manager().riding_selection()
  hd = model_h.get_hd_selection()
  # The water H do not ride.
  assert riding.count(True) == hd.count(True) - 2
  sites_cart = model_h.get_sites_cart()
  pp = rm_h.geometry.pair_proxies(sites_cart=sites_cart)
  bonds = [p.i_seqs for p in pp.bond_proxies.simple] + \
    [(p.i_seq, p.j_seq) for p in pp.bond_proxies.asu]
  nonbonded = [p.i_seqs for p in pp.nonbonded_proxies.simple] + \
    [(p.i_seq, p.j_seq) for p in pp.nonbonded_proxies.asu]
  planes = [list(p.i_seqs) for p in rm_h.geometry.planarity_proxies]
  # Per type: the C++ count excluding riding H against the Python count,
  # and against the model without H (the water H add two bonds, one angle
  # and their nonbonded pairs; the guanidinium plane holds riding H).
  n_bond = pp.bond_proxies.n_total_excluding(exclude=riding)
  assert n_bond == n_without(bonds, riding) == g.n_bond_proxies + 2
  n_nonbonded = pp.nonbonded_proxies.n_total_excluding(exclude=riding)
  assert n_nonbonded == n_without(nonbonded, riding)
  assert n_without(nonbonded, hd) == g.n_nonbonded_proxies
  n_water_h_nonbonded = n_nonbonded - n_without(nonbonded, hd)
  assert n_water_h_nonbonded > 0
  n_angle = geometry_restraints.n_proxies_excluding(
    proxies=rm_h.geometry.angle_proxies, exclude=riding)
  assert n_angle == g.n_angle_proxies + 1
  n_dihedral = geometry_restraints.n_proxies_excluding(
    proxies=rm_h.geometry.dihedral_proxies, exclude=riding)
  assert n_dihedral == g.n_dihedral_proxies
  n_chirality = geometry_restraints.n_proxies_excluding(
    proxies=rm_h.geometry.chirality_proxies, exclude=riding)
  assert n_chirality == g.n_chirality_proxies
  n_planarity = geometry_restraints.n_proxies_excluding(
    proxies=rm_h.geometry.planarity_proxies, exclude=riding)
  assert n_planarity == n_without(planes, riding) == g.n_planarity_proxies - 1
  # The normalization count is the sum of the per-type counts.
  assert e_h.number_of_restraints == \
    n_bond + n_nonbonded + n_angle + n_dihedral + n_chirality + n_planarity
  assert e_h.number_of_restraints == \
    e.number_of_restraints + 2 + 1 + n_water_h_nonbonded - 1
  # Without a selection every restraint is counted, as before.
  n_counted = rm_h.energies_sites(
    sites_cart=sites_cart, normalization_exclude=None).number_of_restraints
  assert n_counted == g_h.n_bond_proxies + g_h.n_nonbonded_proxies + \
    g_h.n_angle_proxies + g_h.n_dihedral_proxies + g_h.n_chirality_proxies + \
    g_h.n_planarity_proxies
  assert n_counted > 3 * e_h.number_of_restraints
  # normalization=False: the selection is ignored.
  rm_plain = mmtbx.restraints.manager(
    geometry=rm_h.geometry, normalization=False)
  e0 = rm_plain.energies_sites(sites_cart=sites_cart, compute_gradients=True)
  e1 = rm_plain.energies_sites(sites_cart=sites_cart, compute_gradients=True,
    normalization_exclude=riding)
  assert e1.number_of_restraints == e0.number_of_restraints
  assert approx_equal(e0.target, e1.target)
  assert approx_equal(e0.gradients, e1.gradients)
  assert e0.target > e_h.target
  # The normalized heavy-atom gradients now have comparable size with and
  # without riding H (the H forces transfer to the parents).
  g_h = model_h.get_riding_h_manager().gradients_reduced_cpp(
    gradients=e_h.gradients, sites_cart=sites_cart, hd_selection=hd)
  ratio = g_h.norm() / e.gradients.norm()
  assert 0.9 <= ratio <= 1.6, ratio

def get_fmodel(model):
  """
  Least-squares fmodel against the model's own structure, with the heavy
  atoms displaced the same way in both models.
  """
  xrs = model.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=1.2, algorithm="direct").f_calc())
  random.seed(0)
  flex.set_random_seed(0)
  r_free_flags = f_obs.generate_r_free_flags()
  sites_cart = xrs.sites_cart()
  hd = model.get_hd_selection()
  random.seed(1)
  for i in range(sites_cart.size()):
    if(not hd[i]):
      sites_cart[i] = [c + random.uniform(-0.15, 0.15) for c in sites_cart[i]]
  if(model.riding_h_manager is not None):
    model.get_riding_h_manager().idealize_riding_h_positions(
      sites_cart=sites_cart)
  model.set_sites_cart(sites_cart)
  sfg = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = model.get_xray_structure(),
    target_name                  = "ls_wunit_k1",
    sf_and_grads_accuracy_params = sfg)
  fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
  heavy = ~hd
  model.set_refinement_flags(refinement_flags.manager(
    individual_sites=True, sites_individual=heavy,
    individual_adp=True, adp_individual_iso=heavy))
  return fmodel

def parse_counts(log):
  """
  The geometry restraint and ADP pair counts of the Normalization counts
  line.
  """
  for line in log.splitlines():
    if(line.strip().startswith("Normalization counts:")):
      numbers = [int(w) for w in re.findall(r"\b\d+\b", line)]
      return line, numbers[0], numbers[-1]
  raise AssertionError("Normalization counts line missing:\n" + log)

def exercise_weights_and_adp():
  model_no_h, model_h = get_models()
  # phenix.refine defaults except shake_sites and shake_adp, which are
  # off here: shaking needs phenix.refinement.memory_eraser and does not
  # touch the counts.
  twp = libtbx.group_args(
    wxc_scale=0.5, wxu_scale=1.0, wc=1.0, wu=1.0, fix_wxc=None,
    fix_wxu=None, shake_sites=False, shake_adp=None, regularize_ncycles=50,
    rmsd_cutoff_for_gradient_filtering=3.0)
  results = []
  for model in [model_no_h, model_h]:
    fmodel = get_fmodel(model)
    log = StringIO()
    # The ADP weight shakes the (all equal) ADPs: same shake for both.
    random.seed(0)
    flex.set_random_seed(0)
    w = weights.weight(
      fmodel=fmodel, model=model, correct_special_position_tolerance=1.0,
      target_weights_params=twp, macro_cycle=1, show_summary=True, log=log)
    results.append((w, log.getvalue()))
  (w, log), (w_h, log_h) = results
  # Both weights are computed on the H-stripped model (a guard on that
  # stripping, which this change does not touch). H are half the atoms of
  # this fixture and their scattering moves the X-ray ADP gradient by
  # about 10%.
  assert abs(w_h.xyz_weights_result.wx / w.xyz_weights_result.wx - 1) < 0.1
  assert abs(w_h.adp_weights_result.wx / w.adp_weights_result.wx - 1) < 0.2
  line, n_geo, n_adp = parse_counts(log)
  line_h, n_geo_h, n_adp_h = parse_counts(log_h)
  assert "riding" not in line, line
  assert "involving riding H excluded" in line_h, line_h
  assert "ADP" in line_h, line_h
  # The printed geometry counts are the model's normalization counts.
  assert n_geo == \
    model_no_h.restraints_manager_energies_sites().number_of_restraints
  assert n_geo_h == \
    model_h.restraints_manager_energies_sites().number_of_restraints
  assert abs(n_adp_h - n_adp) <= 0.01 * n_adp, (line, line_h)
  # The ADP pair count itself, through the model.
  n = model_no_h.energies_adp(
    iso_restraints=None, compute_gradients=False,
    use_hd=False).number_of_restraints
  n_h = model_h.energies_adp(
    iso_restraints=None, compute_gradients=False,
    use_hd=False).number_of_restraints
  assert n == n_adp, (n, n_adp)
  assert abs(n_h - n) <= 0.01 * n, (n, n_h)

class geometry_without_keyword(object):
  """
  A geometry manager whose energies_sites has the signature of the
  subclasses outside cctbx (amber, QM gradients): no normalization_exclude.
  """
  def __init__(self, geometry):
    self.geometry = geometry

  def energies_sites(self, sites_cart, flags=None,
        custom_nonbonded_function=None, compute_gradients=False,
        gradients=None, disable_asu_cache=False, normalization=False,
        external_energy_function=None, extension_objects=[],
        site_labels=None):
    return self.geometry.energies_sites(sites_cart=sites_cart, flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients, gradients=gradients,
      disable_asu_cache=disable_asu_cache, normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects, site_labels=site_labels)

def exercise_geometry_without_keyword():
  """
  Without a selection the keyword is not passed on, so a geometry manager
  that does not know it keeps working.
  """
  model_no_h, model_h = get_models()
  rm = mmtbx.restraints.manager(
    geometry=geometry_without_keyword(model_no_h.restraints_manager.geometry),
    normalization=True)
  e = rm.energies_sites(sites_cart=model_no_h.get_sites_cart(),
    compute_gradients=True)
  assert e.number_of_restraints == \
    model_no_h.restraints_manager_energies_sites().number_of_restraints

def exercise_neutron_riding_wxc():
  """
  With neutron data the weight is computed on the model with its H; the
  restraint gradient must then be normalized as in the minimizer, with
  the riding H excluded from the count.
  """
  model_no_h, model_h = get_models(scattering_table="neutron")
  assert model_h.is_neutron()
  fmodel = get_fmodel(model_h)
  sg = weights.site_gradients(fmodel=fmodel, model=model_h,
    correct_special_position_tolerance=1.0,
    cartesian_dynamics_parameters=None, shake=False, regularize_ncycles=50,
    rmsd_cutoff_for_gradient_filtering=3.0, gradient_filtering=False)
  e_h = model_h.restraints_manager_energies_sites(compute_gradients=True)
  sel = model_h.refinement_flags.sites_individual
  assert approx_equal(sg.gc, e_h.gradients.select(sel))

def exercise_adp_count_without_riding_manager():
  """
  H present but no riding manager (phenix.refine's first macro-cycle with
  hydrogens.refine=riding): the caller says with use_hd whether the H
  pairs count; without it H are taken as individual.
  """
  model_no_h, model_h = get_models(setup_riding=False)
  assert model_h.riding_h_manager is None
  twp = libtbx.group_args(
    wxc_scale=0.5, wxu_scale=1.0, wc=1.0, wu=1.0, fix_wxc=None,
    fix_wxu=None, shake_sites=False, shake_adp=None, regularize_ncycles=50,
    rmsd_cutoff_for_gradient_filtering=3.0)
  counts = []
  for model, use_hd in [(model_no_h, None), (model_h, False), (model_h, None)]:
    fmodel = get_fmodel(model)
    log = StringIO()
    random.seed(0)
    flex.set_random_seed(0)
    weights.weight(
      fmodel=fmodel, model=model, correct_special_position_tolerance=1.0,
      target_weights_params=twp, macro_cycle=1, show_summary=True, log=log,
      use_hd=use_hd)
    counts.append(parse_counts(log.getvalue())[2])
  assert counts[0] == counts[1], counts
  assert counts[2] > counts[1], counts

def run():
  exercise_geometry_count()
  exercise_weights_and_adp()
  exercise_geometry_without_keyword()
  exercise_neutron_riding_wxc()
  exercise_adp_count_without_riding_manager()
  print("OK")

if (__name__ == "__main__"):
  run()
