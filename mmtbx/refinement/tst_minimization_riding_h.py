from __future__ import absolute_import, division, print_function
import random
import time
import iotbx.pdb
import mmtbx
import mmtbx.f_model
import mmtbx.model
import scitbx.lbfgs
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.refinement import minimization, refinement_flags, weights

#-----------------------------------------------------------------------------
# Coordinate refinement with riding H atoms (mmtbx.refinement.minimization):
# the gradient handed to the minimizer must be the derivative of the target
# with the H atoms following their parents, which is checked against finite
# differences; H atoms are never parameters and always ride, with and without
# an atom selection.
#-----------------------------------------------------------------------------

# Distorted ARG residue with riding H atoms.
pdb_str = """\
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  50       6.642  11.452   6.691  1.00 10.00           N
ATOM      2  CA  ARG A  50       7.570  10.240   6.682  1.00 10.00           C
ATOM      3  C   ARG A  50       9.005  10.280   6.707  1.00 10.00           C
ATOM      4  O   ARG A  50       9.602  11.108   6.810  1.00 10.00           O
ATOM      5  CB  ARG A  50       6.962   9.371   5.075  1.00 10.00           C
ATOM      6  CG  ARG A  50       5.382   9.275   4.806  1.00 10.00           C
ATOM      7  CD  ARG A  50       4.830   7.928   6.084  1.00 10.00           C
ATOM      8  NE  ARG A  50       5.434   6.841   5.374  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.450   6.116   6.035  1.00 10.00           C
ATOM     10  NH1 ARG A  50       7.392   6.694   6.803  1.00 10.00           N
ATOM     11  NH2 ARG A  50       6.898   5.026   5.305  1.00 10.00           N
ATOM     12  HA  ARG A  50       6.817   9.405   7.203  1.00 10.00           H
ATOM     13  HB2 ARG A  50       7.451   9.961   4.409  1.00 10.00           H
ATOM     14  HB3 ARG A  50       7.391   8.588   5.297  1.00 10.00           H
ATOM     15  HG2 ARG A  50       5.016  10.089   5.174  1.00 10.00           H
ATOM     16  HG3 ARG A  50       5.385   9.031   4.141  1.00 10.00           H
ATOM     17  HD2 ARG A  50       5.268   8.147   6.643  1.00 10.00           H
ATOM     18  HD3 ARG A  50       3.909   8.226   5.899  1.00 10.00           H
ATOM     19  HE  ARG A  50       5.208   6.478   4.962  1.00 10.00           H
ATOM     20 HH11 ARG A  50       6.970   7.301   7.302  1.00 10.00           H
ATOM     21 HH12 ARG A  50       7.866   6.356   7.069  1.00 10.00           H
ATOM     22 HH21 ARG A  50       6.604   4.626   5.005  1.00 10.00           H
ATOM     23 HH22 ARG A  50       7.522   4.396   5.734  1.00 10.00           H
"""

class target_weights(object):
  def __init__(self, w):
    self.xyz_weights_result = weights.weights(wx=1, wx_scale=1, w=w)

class capture(minimization.lbfgs):
  """
  The minimizer, recording every (f, g) it hands to LBFGS. At the first
  evaluation it also differentiates the target by finite differences, with
  the target functor the minimizer prepared for the run but without its
  shift and gradient machinery: the heavy atoms displaced directly, the H
  atoms placed on them.
  """
  def __init__(self, fd_atoms=None, eps=1.e-4, **kw):
    self.fd_atoms = fd_atoms
    self.eps = eps
    self.evaluations = []
    self.fd = {}
    minimization.lbfgs.__init__(self, **kw)

  def compute_functional_and_gradients(self):
    f, g = minimization.lbfgs.compute_functional_and_gradients(self)
    if(len(self.evaluations) == 0 and self.fd_atoms is not None):
      self.finite_differences(f)
    self.evaluations.append((f, g.deep_copy()))
    return f, g

  def target(self, sites_cart):
    sites_cart = sites_cart.deep_copy()
    self.riding_h_manager.idealize_riding_h_positions(sites_cart=sites_cart)
    self.xray_structure.set_sites_cart(sites_cart)
    self.fmodels.update_xray_structure(
      xray_structure=self.xray_structure, update_f_calc=True)
    result = self.fmodels.target_and_gradients(
      weights=self.weights, compute_gradients=False).target()
    if(self.restraints_manager is not None and self.weights.w > 0):
      result += self.weights.w * \
        self.model.restraints_manager_energies_sites().target
    return result

  def finite_differences(self, f):
    sites_start = self.xray_structure.sites_cart()
    assert approx_equal(self.target(sites_start), f, eps=1.e-8*abs(f))
    for i_seq in self.fd_atoms:
      for k in range(3):
        values = []
        for sign in [-1, 1]:
          sites_cart = sites_start.deep_copy()
          site = list(sites_cart[i_seq])
          site[k] += sign*self.eps
          sites_cart[i_seq] = site
          values.append(self.target(sites_cart))
        self.fd[(i_seq, k)] = (values[1]-values[0])/(2*self.eps)
    self.target(sites_start)

def get_model_and_fmodels():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  model.setup_riding_h_manager(idealize=True)
  xrs = model.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=1.2, algorithm="direct").f_calc())
  random.seed(0)
  flex.set_random_seed(0)
  r_free_flags = f_obs.generate_r_free_flags()
  # Move the heavy atoms away from the "observed" structure so that the X-ray
  # target has a gradient, and let the H atoms follow.
  sites_cart = xrs.sites_cart()
  hd = model.get_hd_selection()
  for i in range(sites_cart.size()):
    if(not hd[i]):
      sites_cart[i] = [c + random.uniform(-0.15, 0.15) for c in sites_cart[i]]
  model.get_riding_h_manager().idealize_riding_h_positions(
    sites_cart=sites_cart)
  model.set_sites_cart(sites_cart)
  sfg = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg.algorithm = "direct"
  # A least-squares target: no prepared state, so the function can be
  # evaluated outside the minimizer exactly as inside it.
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = model.get_xray_structure(),
    target_name                  = "ls_wunit_k1",
    sf_and_grads_accuracy_params = sfg)
  assert fmodel.xray_structure is model.get_xray_structure()
  fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
  return model, mmtbx.fmodels(fmodel_xray=fmodel)

def set_sites_individual(model, selection_string):
  n = model.get_number_of_atoms()
  if(selection_string is None):
    selection = flex.bool(n, True)
  else:
    selection = model.selection(selection_string)
  model.set_refinement_flags(refinement_flags.manager(
    individual_sites=True, sites_individual=selection))
  model.set_refine_individual_sites()
  return selection

def h_displacement_from_ideal(model, sites_cart):
  idealized = sites_cart.deep_copy()
  model.get_riding_h_manager().idealize_riding_h_positions(
    sites_cart=idealized)
  return flex.max((idealized - sites_cart).norms())

def exercise_gradient(selection_string, w, tolerance=1.e-3):
  model, fmodels = get_model_and_fmodels()
  selection = set_sites_individual(model, selection_string)
  hd = model.get_hd_selection()
  # Every H atom carries a site gradient during the minimization (the chain
  # rule needs its X-ray gradient), so the parameters are the selected atoms
  # plus all H atoms, packed in atom order, three per atom.
  in_x = selection | hd
  atoms_in_x = in_x.iselection()
  restraints_manager = model.get_restraints_manager() if w > 0 else None
  minimized = capture(
    fd_atoms                 = (selection & ~hd).iselection(),
    fmodels                  = fmodels,
    restraints_manager       = restraints_manager,
    model                    = model,
    refine_xyz               = True,
    target_weights           = target_weights(w),
    collect_monitor          = False,
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations=0))
  f0, g0 = minimized.evaluations[0]
  assert g0.size() == 3*in_x.count(True), (g0.size(), 3*in_x.count(True))
  max_deviation = 0
  for i, i_seq in enumerate(atoms_in_x):
    for k in range(3):
      if(hd[i_seq]):
        # Not a parameter: the function does not depend on it.
        assert g0[3*i+k] == 0, (i_seq, k, g0[3*i+k])
        continue
      fd = minimized.fd[(i_seq, k)]
      max_deviation = max(max_deviation, abs(g0[3*i+k]-fd))
      assert abs(g0[3*i+k]-fd) < tolerance, (i_seq, k, g0[3*i+k], fd)
  assert len(minimized.fd) == 3*(selection & ~hd).count(True)
  print("  selection=%s w=%s: max |g - fd| = %.2e" % (
    selection_string, w, max_deviation))

def exercise_minimization(selection_string, max_iterations=30):
  model, fmodels = get_model_and_fmodels()
  selection = set_sites_individual(model, selection_string)
  hd = model.get_hd_selection()
  xrs = model.get_xray_structure()
  n_parameters = xrs.n_parameters()
  n_parameterized_h = len(model.get_riding_h_manager().parameterization_cpp)
  sites_start = xrs.sites_cart()
  minimized = capture(
    fmodels                  = fmodels,
    restraints_manager       = model.get_restraints_manager(),
    model                    = model,
    refine_xyz               = True,
    target_weights           = target_weights(w=1),
    collect_monitor          = False,
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations))
  assert minimized.minimizer.iter() > 5
  f_first = minimized.evaluations[0][0]
  f_last = minimized.evaluations[-1][0]
  assert f_last < f_first, (f_first, f_last)
  # The minimizer moves the xray_structure; callers sync the model from it.
  assert fmodels.fmodel_xray().xray_structure is xrs
  sites_cart = xrs.sites_cart()
  # Unselected atoms stayed, selected heavy atoms moved, all H atoms ride.
  moved = (sites_cart-sites_start).norms()
  if((~selection).count(True) > 0):
    assert flex.max(moved.select(~selection & ~hd)) < 1.e-6
  assert flex.max(moved.select(selection & ~hd)) > 0.05
  assert h_displacement_from_ideal(model, sites_cart) < 1.e-6
  # The minimizer leaves the gradient flags and the riding manager as found.
  assert xrs.n_parameters() == n_parameters
  assert len(model.get_riding_h_manager().parameterization_cpp) == \
    n_parameterized_h

def run():
  exercise_gradient(selection_string=None, w=1)
  exercise_gradient(selection_string=None, w=0)
  exercise_gradient(selection_string="not (name CZ or name HA or name HB2)",
    w=1)
  exercise_minimization(selection_string=None)
  exercise_minimization(
    selection_string="not (name CZ or name HA or name HB2)")

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
