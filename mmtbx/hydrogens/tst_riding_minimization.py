from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
import mmtbx.refinement.geometry_minimization
import scitbx.lbfgs
from cctbx.array_family import flex
from cctbx.geometry_restraints.lbfgs import lbfgs
from libtbx.utils import null_out

#-----------------------------------------------------------------------------
# Geometry minimization with the riding H model: the H atoms must follow their
# parents, with and without an atom selection, and the minimizer must converge.
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

def get_model():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  model.setup_riding_h_manager(idealize=True)
  return model

def h_displacement_from_ideal(model, sites_cart):
  """Largest distance of a riding H atom from its idealized position."""
  idealized = sites_cart.deep_copy()
  model.get_riding_h_manager().idealize_riding_h_positions(
    sites_cart=idealized)
  return flex.max((idealized - sites_cart).norms())

def deviations(model, sites_cart):
  es = model.get_restraints_manager().geometry.energies_sites(
    sites_cart=sites_cart, compute_gradients=False)
  return es.bond_deviations()[2], es.angle_deviations()[2]

def minimize(model, sites_cart_selection=None, max_iterations=200):
  sites_cart = model.get_sites_cart().deep_copy()
  minimized = lbfgs(
    sites_cart                  = sites_cart,
    geometry_restraints_manager = model.get_restraints_manager().geometry,
    riding_h_manager            = model.get_riding_h_manager(),
    sites_cart_selection        = sites_cart_selection,
    lbfgs_termination_params    = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations))
  return minimized, sites_cart

def exercise_no_selection():
  model = get_model()
  start = model.get_sites_cart().deep_copy()
  hd = model.get_hd_selection()
  minimized, sites_cart = minimize(model)
  assert minimized.final_target_value < minimized.first_target_value
  # The minimizer converged to good geometry: a wrong shift bookkeeping
  # (shifts accumulating between evaluations) blows the geometry up instead.
  b_rmsd, a_rmsd = deviations(model, sites_cart)
  assert b_rmsd < 0.01, b_rmsd
  assert a_rmsd < 1.5, a_rmsd
  # Heavy atoms moved, and the H atoms rode along with them.
  assert flex.max((sites_cart - start).select(~hd).norms()) > 0.1
  assert flex.max((sites_cart - start).select(hd).norms()) > 0.1
  assert h_displacement_from_ideal(model, sites_cart) < 1.e-6

def exercise_selection():
  model = get_model()
  start = model.get_sites_cart().deep_copy()
  hd = model.get_hd_selection()
  # What phenix.geometry_minimization passes with selection=all: the riding
  # model must still be used.
  all_atoms = flex.bool(start.size(), True)
  minimized_all, sites_all = minimize(model, sites_cart_selection=all_atoms)
  minimized, sites_cart = minimize(model)
  assert minimized_all.final_target_value == minimized.final_target_value
  assert flex.max((sites_all - sites_cart).norms()) < 1.e-6
  # A partial selection: unselected atoms stay, selected ones move, and the
  # H atoms ride on whatever their parents do.
  backbone = model.selection("name N or name CA or name C or name O")
  minimized, sites_cart = minimize(model, sites_cart_selection=~backbone)
  assert minimized.final_target_value < minimized.first_target_value
  assert flex.max((sites_cart - start).select(backbone).norms()) == 0
  assert flex.max((sites_cart - start).select(~backbone & ~hd).norms()) > 0.1
  assert h_displacement_from_ideal(model, sites_cart) < 1.e-6

def exercise_run2():
  # The route phenix.geometry_minimization takes.
  model = get_model()
  hd = model.get_hd_selection()
  mmtbx.refinement.geometry_minimization.run2(
    restraints_manager       = model.get_restraints_manager(),
    pdb_hierarchy            = model.get_hierarchy(),
    riding_h_manager         = model.get_riding_h_manager(),
    correct_special_position_tolerance = 1.0,
    max_number_of_iterations = 200,
    number_of_macro_cycles   = 1,
    selection                = flex.bool(hd.size(), True),
    bond                     = True,
    nonbonded                = True,
    angle                    = True,
    dihedral                 = True,
    chirality                = True,
    planarity                = True,
    parallelity              = True,
    fix_rotamer_outliers     = False,
    log                      = null_out())
  sites_cart = model.get_hierarchy().atoms().extract_xyz()
  b_rmsd, a_rmsd = deviations(model, sites_cart)
  assert b_rmsd < 0.01, b_rmsd
  assert h_displacement_from_ideal(model, sites_cart) < 1.e-6

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_no_selection()
  exercise_selection()
  exercise_run2()
  print("OK. Time: %8.3f"%(time.time()-t0))
