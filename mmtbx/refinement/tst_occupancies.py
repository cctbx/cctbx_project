from __future__ import absolute_import, division, print_function
import time
import random
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
import mmtbx.refinement.refinement_flags
from mmtbx.refinement import occupancies
from cctbx.array_family import flex
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

# (atom name, altloc, residue name, residue number, x, y, z)
# true occupancy and B, starting occupancy and B (None: atom is not in truth)
atoms_info = [
  (" N  ", " ", "ALA", 1,  3.0,  3.0,  3.0, 1.0, 10, 1.0, 10),
  (" CA ", " ", "ALA", 1,  4.4,  3.3,  3.2, 1.0, 10, 1.0, 10),
  (" C  ", " ", "ALA", 1,  5.1,  2.2,  4.0, 1.0, 10, 1.0, 10),
  (" O  ", " ", "ALA", 1,  4.6,  1.1,  4.2, 1.0, 10, 1.0, 10),
  (" CB ", " ", "ALA", 1,  5.0,  4.6,  2.6, 1.0, 10, 1.0, 10),
  # two conformers
  (" N  ", " ", "SER", 2,  6.4,  2.5,  4.5, 1.0, 10, 1.0, 10),
  (" CA ", " ", "SER", 2,  7.2,  1.5,  5.3, 1.0, 10, 1.0, 10),
  (" C  ", " ", "SER", 2,  8.6,  2.0,  5.6, 1.0, 10, 1.0, 10),
  (" O  ", " ", "SER", 2,  9.0,  3.1,  5.2, 1.0, 10, 1.0, 10),
  (" CB ", "A", "SER", 2,  6.5,  1.1,  6.6, 0.7, 10, 0.5, 10),
  (" OG ", "A", "SER", 2,  6.3,  2.2,  7.5, 0.7, 12, 0.5, 12),
  (" CB ", "B", "SER", 2,  7.3,  0.2,  4.5, 0.3, 10, 0.5, 10),
  (" OG ", "B", "SER", 2,  8.3, -0.7,  5.0, 0.3, 12, 0.5, 12),
  # three conformers
  (" N  ", " ", "SER", 3,  9.4,  1.2,  6.3, 1.0, 10, 1.0, 10),
  (" CA ", " ", "SER", 3, 10.8,  1.5,  6.7, 1.0, 10, 1.0, 10),
  (" C  ", " ", "SER", 3, 11.7,  0.3,  6.4, 1.0, 10, 1.0, 10),
  (" O  ", " ", "SER", 3, 11.3, -0.7,  5.8, 1.0, 10, 1.0, 10),
  (" CB ", "A", "SER", 3, 10.9,  1.9,  8.2, 0.5, 10, 0.34, 10),
  (" OG ", "A", "SER", 3, 10.1,  3.0,  8.5, 0.5, 12, 0.34, 12),
  (" CB ", "B", "SER", 3, 11.3,  2.7,  5.9, 0.3, 10, 0.33, 10),
  (" OG ", "B", "SER", 3, 12.7,  2.9,  6.1, 0.3, 12, 0.33, 12),
  (" CB ", "C", "SER", 3, 11.6,  2.3,  7.6, 0.2, 10, 0.33, 10),
  (" OG ", "C", "SER", 3, 12.6,  3.1,  8.3, 0.2, 12, 0.33, 12),
  # individual occupancies: best value above 1 (B is too large in the model),
  # inside (0,1), and below or at 0 (atom is not in truth)
  (" O  ", " ", "HOH", 4,  2.0,  8.0,  8.0, 1.0, 10, 0.7, 20),
  (" O  ", " ", "HOH", 5,  8.0,  8.0,  2.0, 0.4, 15, 0.8, 15),
  (" O  ", " ", "HOH", 6,  2.0,  8.0,  2.0, None, None, 0.5, 15),
]

def pdb_str(truth):
  lines = ["CRYST1   16.000   14.000   12.000  90.00  90.00  90.00 P 1"]
  for i, a in enumerate(atoms_info):
    name, altloc, resname, resseq, x, y, z, occ_t, b_t, occ_s, b_s = a
    if(truth):
      occ, b = occ_t, b_t
    else:
      occ, b = occ_s, b_s
    if(occ is None): continue
    record = "ATOM  "
    if(resname == "HOH"): record = "HETATM"
    lines.append("%s%5d %s%s%s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"%(
      record, i+1, name, altloc, resname, resseq, x, y, z, occ, b,
      name.strip()[0]))
  lines.append("END")
  return "\n".join(lines)

def get_model_and_fmodels(target_name="ls_wunit_k1"):
  xrs_true = iotbx.pdb.input(
    source_info=None, lines=pdb_str(truth=True)).xray_structure_simple()
  f_obs = abs(xrs_true.structure_factors(d_min=1.2).f_calc())
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(
      source_info=None, lines=pdb_str(truth=False)),
    log = null_out())
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = model.get_xray_structure(),
    target_name    = target_name)
  selections = occupancies.occupancy_selections(model = model)
  assert [len(s) for s in selections] == [2, 3, 1, 1, 1]
  model.set_refinement_flags(
    flags = mmtbx.refinement.refinement_flags.manager(
      occupancies = True, s_occupancies = selections))
  return model, mmtbx.fmodels(fmodel_xray = fmodel), selections

def exercise_gradients():
  """
  Gradients given to the minimizer agree with finite differences.
  """
  model, fmodels, selections = get_model_and_fmodels()
  fmodels.fmodel_xray().sfg_params.algorithm = "direct"
  fmodels.update_xray_structure(update_f_calc = True)
  scatterers = fmodels.fmodel_xray().xray_structure.scatterers()
  scatterers.flags_set_grads(state=False)
  i_selection = flex.size_t()
  for s in selections:
    for ss in s: i_selection.extend(ss)
  scatterers.flags_set_grad_occupancy(iselection = i_selection)
  m = occupancies.minimizer(
    fmodels                  = fmodels,
    selections               = selections,
    par_initial = flex.double([0.6, 0.4, 0.5, 0.2, 0.3, 0.7, 0.8, 0.5]),
    max_number_of_iterations = 0)
  # one parameter per group of atoms, less one per constrained group
  assert m.x.size() == 6
  x0 = m.x.deep_copy()
  _, _, g = m.compute_functional_and_gradients()
  g = g.deep_copy()
  eps = 1.e-4
  for i in range(x0.size()):
    fs = []
    for sign in [1, -1]:
      m.x = x0.deep_copy()
      m.x[i] += sign*eps
      fs.append(m.compute_functional_and_gradients()[1])
    assert approx_equal(g[i], (fs[0]-fs[1])/(2*eps), 1.e-5*max(1, abs(g[i])))
  assert flex.max(flex.abs(g)) > 1.e-3

def exercise_refinement():
  """
  Refined occupancies match the truth, constrained groups add up to one, and
  occupancies whose best value is outside [0,1] end up at the bound.
  """
  for target_name in ["ls_wunit_k1", "ml"]:
    model, fmodels, selections = get_model_and_fmodels(target_name=target_name)
    r_start = fmodels.fmodel_xray().r_work()
    occupancies.manager(
      fmodels                  = fmodels,
      model                    = model,
      max_number_of_iterations = 50,
      number_of_macro_cycles   = 3,
      occupancy_max            = 1,
      occupancy_min            = 0,
      log                      = None)
    occ = model.get_xray_structure().scatterers().extract_occupancies()
    assert approx_equal(occ,
      fmodels.fmodel_xray().xray_structure.scatterers().extract_occupancies())
    q = []
    for s in selections:
      for ss in s:
        o = occ.select(ss)
        assert approx_equal(flex.min(o), flex.max(o))
        q.append(o[0])
    tol = 0.02
    assert approx_equal(q[0:2], [0.7, 0.3], tol)
    assert approx_equal(q[0]+q[1], 1, 1.e-6)
    assert approx_equal(q[2:5], [0.5, 0.3, 0.2], tol)
    assert approx_equal(q[2]+q[3]+q[4], 1, 1.e-6)
    assert approx_equal(q[5], 1, 1.e-6)
    assert approx_equal(q[6], 0.4, tol)
    assert q[7] < 0.03
    assert flex.min(occ) >= 0 and flex.max(occ) <= 1
    assert fmodels.fmodel_xray().r_work() < 0.7*r_start

if(__name__ == "__main__"):
  t0 = time.time()
  random.seed(0)
  flex.set_random_seed(0)
  exercise_gradients()
  exercise_refinement()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
