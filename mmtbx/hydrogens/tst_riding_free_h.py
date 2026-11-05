from __future__ import absolute_import, division, print_function
import random
import time
import iotbx.pdb
import mmtbx.f_model
import mmtbx.model
import mmtbx.refinement.geometry_minimization
import scitbx.lbfgs
from cctbx.array_family import flex
from cctbx.geometry_restraints.lbfgs import lbfgs
from libtbx.utils import null_out
from mmtbx.refinement import minimization, refinement_flags, weights

#-----------------------------------------------------------------------------
# Minimization with the riding H model when the manager does not parameterize
# an H atom: a truncated side chain (ARG cut at CB: HB2, HB3 have no third
# neighbor), water H and the H of an ammonium ion. The drivers used to freeze
# every H of the model while the parent moved, so the X-H lengths in the
# output were arbitrary. Now an H without a parameterization is an ordinary
# refinable atom and ends at the restraint target, while the parameterized H
# (HA of the ARG) still ride, in every driver: cctbx.geometry_restraints.lbfgs
# (model.geometry_minimization, regularize), mmtbx.refinement.
# geometry_minimization.run2 (phenix.geometry_minimization) and
# mmtbx.refinement.minimization (phenix.refine).
#-----------------------------------------------------------------------------

pdb_arg_truncated = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  50       6.642  11.452   6.691  1.00 10.00           N
ATOM      2  CA  ARG A  50       7.570  10.240   6.682  1.00 10.00           C
ATOM      3  C   ARG A  50       9.005  10.280   6.707  1.00 10.00           C
ATOM      4  O   ARG A  50       9.602  11.108   6.810  1.00 10.00           O
ATOM      5  CB  ARG A  50       6.962   9.371   5.075  1.00 10.00           C
ATOM      6  HA  ARG A  50       6.817   9.405   7.203  1.00 10.00           H
ATOM      7  HB2 ARG A  50       7.451   9.961   4.409  1.00 10.00           H
ATOM      8  HB3 ARG A  50       7.391   8.588   5.297  1.00 10.00           H
"""

pdb_water = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH S   1      10.000  10.000  10.000  1.00 10.00           O
HETATM    2  H1  HOH S   1      10.800  10.300  10.000  1.00 10.00           H
HETATM    3  H2  HOH S   1       9.600  10.700  10.300  1.00 10.00           H
"""

pdb_ammonium = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  N   NH4 B   1       5.000   5.000   5.000  1.00 10.00           N
HETATM    2  HN1 NH4 B   1       5.600   5.600   5.600  1.00 10.00           H
HETATM    3  HN2 NH4 B   1       4.400   4.400   5.600  1.00 10.00           H
HETATM    4  HN3 NH4 B   1       5.600   4.400   4.400  1.00 10.00           H
HETATM    5  HN4 NH4 B   1       4.400   5.600   4.400  1.00 10.00           H
"""

# (pdb string, parent of the free H, free H, parameterized H and its parent)
cases = [
  ("truncated ARG", pdb_arg_truncated, "CB", ["HB2", "HB3"], ("HA", "CA")),
  ("water", pdb_water, "O", ["H1", "H2"], None),
  ("ammonium", pdb_ammonium, "N", ["HN1", "HN2", "HN3", "HN4"], None),
]

def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  model.setup_riding_h_manager(idealize=True)
  return model

def i_seqs(model, names):
  atoms = model.get_hierarchy().atoms()
  by_name = dict((a.name.strip(), a.i_seq) for a in atoms)
  return [by_name[n] for n in names]

def ideal_length(model, i, j):
  bpt = model.get_restraints_manager().geometry.bond_params_table
  return bpt.lookup(i, j).distance_ideal

def regularize_geometry(model):
  """
  Ideal X-H lengths for the free H: the riding manager does not place them,
  so an idealized model has them wherever the input put them.
  """
  sites_cart = model.get_sites_cart().deep_copy()
  lbfgs(
    sites_cart                  = sites_cart,
    geometry_restraints_manager = model.get_restraints_manager().geometry,
    lbfgs_termination_params    = scitbx.lbfgs.termination_parameters(
      max_iterations=500))
  model.set_sites_cart(sites_cart)

def distort(model, i_parent, i_free):
  """
  Move the parent and lengthen every X-H bond of the free H, so that no
  position of the parent alone satisfies the restraints.
  """
  sites_cart = model.get_sites_cart().deep_copy()
  parent = flex.vec3_double([sites_cart[i_parent]])[0]
  for i in i_free:
    v = flex.vec3_double([sites_cart[i]]) - flex.vec3_double([parent])
    v = v * (1 + 0.3/v.norms()[0])
    sites_cart[i] = (flex.vec3_double([parent]) + v)[0]
  sites_cart[i_parent] = (flex.vec3_double([parent]) +
    flex.vec3_double([(0.3, 0.1, -0.2)]))[0]
  # The riding H follow the parent; the free H stay as set. Not through
  # model.idealize_h_riding(): that also rebuilds the water H about the O.
  model.get_riding_h_manager().idealize_riding_h_positions(
    sites_cart=sites_cart)
  model.set_sites_cart(sites_cart)

def check(label, model, sites_start, sites_cart, i_parent, i_free, riding):
  moved = (sites_cart - sites_start).norms()
  assert moved[i_parent] > 0.05, (label, moved[i_parent])
  for i in i_free:
    # The free H moved with the parent, to the restraint target.
    assert moved[i] > 0.05, (label, i, moved[i])
    d = (flex.vec3_double([sites_cart[i]]) -
         flex.vec3_double([sites_cart[i_parent]])).norms()[0]
    ideal = ideal_length(model, i_parent, i)
    assert abs(d - ideal) < 0.01, (label, i, d, ideal)
  if(riding is not None):
    i_h, i_p = riding
    idealized = sites_cart.deep_copy()
    model.get_riding_h_manager().idealize_riding_h_positions(
      sites_cart=idealized)
    assert (idealized - sites_cart).norms()[i_h] < 1.e-6, label
    d = (flex.vec3_double([sites_cart[i_h]]) -
         flex.vec3_double([sites_cart[i_p]])).norms()[0]
    assert abs(d - ideal_length(model, i_p, i_h)) < 0.01, (label, d)

def prepared(name, pdb_str, parent, free, riding):
  model = get_model(pdb_str)
  i_parent = i_seqs(model, [parent])[0]
  i_free = i_seqs(model, free)
  if(riding is not None): riding = tuple(i_seqs(model, riding))
  # The free H are indeed free, the control H rides.
  parameterized = [rc.ih for rc in
    model.get_riding_h_manager().parameterization_cpp]
  for i in i_free: assert i not in parameterized, (name, i)
  if(riding is not None): assert riding[0] in parameterized, name
  regularize_geometry(model)
  distort(model, i_parent, i_free)
  return model, i_parent, i_free, riding

def exercise_lbfgs(name, pdb_str, parent, free, riding):
  model, i_parent, i_free, riding = prepared(name, pdb_str, parent, free,
    riding)
  sites_start = model.get_sites_cart().deep_copy()
  sites_cart = sites_start.deep_copy()
  minimized = lbfgs(
    sites_cart                  = sites_cart,
    geometry_restraints_manager = model.get_restraints_manager().geometry,
    riding_h_manager            = model.get_riding_h_manager(),
    lbfgs_termination_params    = scitbx.lbfgs.termination_parameters(
      max_iterations=500))
  assert minimized.final_target_value < minimized.first_target_value, (
    minimized.first_target_value, minimized.final_target_value)
  check("lbfgs %s" % name, model, sites_start, sites_cart, i_parent, i_free,
    riding)

def exercise_run2(name, pdb_str, parent, free, riding):
  model, i_parent, i_free, riding = prepared(name, pdb_str, parent, free,
    riding)
  sites_start = model.get_sites_cart().deep_copy()
  hierarchy = model.get_hierarchy()
  mmtbx.refinement.geometry_minimization.run2(
    restraints_manager       = model.get_restraints_manager(),
    pdb_hierarchy            = hierarchy,
    riding_h_manager         = model.get_riding_h_manager(),
    correct_special_position_tolerance = 1.0,
    max_number_of_iterations = 500,
    number_of_macro_cycles   = 1,
    selection                = flex.bool(sites_start.size(), True),
    bond                     = True,
    nonbonded                = True,
    angle                    = True,
    dihedral                 = True,
    chirality                = True,
    planarity                = True,
    parallelity              = True,
    fix_rotamer_outliers     = False,
    log                      = null_out())
  sites_cart = hierarchy.atoms().extract_xyz()
  check("run2 %s" % name, model, sites_start, sites_cart, i_parent, i_free,
    riding)

class target_weights(object):
  def __init__(self, w):
    self.xyz_weights_result = weights.weights(wx=1, wx_scale=1, w=w)

def exercise_refinement(name, pdb_str, parent, free, riding):
  """
  The phenix.refine coordinate minimizer: the data are computed from the
  regular model, so the data and the restraints agree on the answer.
  """
  model, i_parent, i_free, riding = prepared(name, pdb_str, parent, free,
    riding)
  xrs = model.get_xray_structure()
  # f_obs of the regular model: minimize the distorted one against them.
  regular = model.deep_copy()
  regularize_geometry(regular)
  f_obs = abs(regular.get_xray_structure().structure_factors(
    d_min=1.0, algorithm="direct").f_calc())
  random.seed(0)
  flex.set_random_seed(0)
  r_free_flags = f_obs.generate_r_free_flags()
  sfg = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs,
    target_name                  = "ls_wunit_k1",
    sf_and_grads_accuracy_params = sfg)
  fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
  fmodels = mmtbx.fmodels(fmodel_xray=fmodel)
  n = model.get_number_of_atoms()
  model.set_refinement_flags(refinement_flags.manager(
    individual_sites=True, sites_individual=flex.bool(n, True)))
  model.set_refine_individual_sites()
  sites_start = xrs.sites_cart()
  minimized = minimization.lbfgs(
    fmodels                  = fmodels,
    restraints_manager       = model.get_restraints_manager(),
    model                    = model,
    refine_xyz               = True,
    target_weights           = target_weights(w=1),
    collect_monitor          = False,
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations=200))
  assert minimized.minimizer.iter() > 5
  sites_cart = xrs.sites_cart()
  check("refinement %s" % name, model, sites_start, sites_cart, i_parent,
    i_free, riding)

def run():
  for name, pdb_str, parent, free, riding in cases:
    exercise_lbfgs(name, pdb_str, parent, free, riding)
    exercise_run2(name, pdb_str, parent, free, riding)
    exercise_refinement(name, pdb_str, parent, free, riding)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
