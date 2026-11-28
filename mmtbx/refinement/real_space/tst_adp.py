from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import mmtbx.refinement.real_space.adp as adp
from mmtbx.secondary_structure.build import ss_idealization as ssb
from cctbx import uctbx
from iotbx import map_manager, map_model_manager
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from libtbx import group_args

def get_inputs(d_min = 2.5, b_add = 0, b_start = 50, buffer_layer = 5):
  """
  Helix with B-factors that vary along the chain (plus b_add), map calculated
  from it at resolution d_min, and the same model with all B set to one value
  (refinement start).
  """
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = buffer_layer)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry())
  model.setup_scattering_dictionaries(scattering_table = "electron")
  b_answer = flex.double()
  for i, rg in enumerate(ph.residue_groups()):
    for a in rg.atoms():
      b = 30 + 3*abs(i-15) + b_add
      if(not a.name.strip() in ["N", "CA", "C", "O"]): b += 10
      b_answer.append(b)
  model.set_b_iso(values = b_answer)
  fc = model.get_xray_structure().structure_factors(d_min = d_min).f_calc()
  map_data = fc.fft_map(resolution_factor = 0.25).real_map_unpadded()
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = False)
  model.set_b_iso(values = flex.double(b_answer.size(), b_start))
  mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)
  return mmm, b_answer

class count_evaluations(object):
  """
  Count target evaluations in each restrained minimization of individual B.
  """
  def __init__(self):
    self.per_minimization = []
    self._init, self._compute = adp.tg.__init__, adp.tg._compute
    counter = self
    def init(self, *args, **kwds):
      counter.per_minimization.append(0)
      counter._init(self, *args, **kwds)
    def compute(self, x):
      counter.per_minimization[-1] += 1
      return counter._compute(self, x)
    adp.tg.__init__, adp.tg._compute = init, compute

  def restore(self):
    adp.tg.__init__, adp.tg._compute = self._init, self._compute

def exercise(d_min = 2.5):
  """
  Individual isotropic B refinement against a map: result is close to the
  B-factors the map was calculated with, and each minimization stops once the
  target stops decreasing (it used to run to the limit of 100 iterations in
  numerical noise).
  """
  mmm, b_answer = get_inputs(d_min = d_min)
  log = StringIO()
  counter = count_evaluations()
  try:
    o = adp.ncs_aware_refinement(
      map_model_manager = mmm,
      d_min             = d_min,
      atom_radius       = 3,
      individual        = True,
      restraints_weight = 1,
      group_mode        = "one_per_residue",
      log               = log)
  finally:
    counter.restore()
  b = o.mmm.model().get_b_iso()
  cc = flex.linear_correlation(x = b, y = b_answer).coefficient()
  diff = flex.mean(flex.abs(b - b_answer))
  n = counter.per_minimization
  print("evaluations per minimization:", n)
  print("CC(B, B answer): %6.4f mean |B - B answer|: %6.2f"%(cc, diff))
  lines = [l for l in log.getvalue().splitlines() if l.count("r_work")]
  assert lines[-1].startswith("r_work (final)"), lines
  r_final = float(lines[-1].split()[-1])
  assert len(n) >= 2, n       # at least two restraint weights tried
  assert max(n) < 90, n       # used to be about 100 and more for each
  assert sum(n) < 60*len(n), n
  assert cc > 0.99, cc
  assert diff < 1.5, diff # was 10.1 (X-ray scattering table, mask)
  assert abs(flex.mean(b - b_answer)) < 1
  assert r_final < 0.03, r_final
  # the number of evaluations is reported in the log
  trials = [l for l in lines if l.count("restraints_weight:")]
  assert len(trials) == len(n), [trials, n]
  for l, n_ in zip(trials, n):
    assert l.split()[-2] == "evaluations:" and int(l.split()[-1]) == n_, l

def exercise_tg_evaluates_once():
  """
  Creating the target/gradient calculator evaluates the target once (it used
  to do it twice).
  """
  mmm, b_answer = get_inputs()
  model = mmm.model()
  fmodel = adp.map_and_model_to_fmodel(
    map_data       = mmm.map_data().deep_copy(),
    xray_structure = model.get_xray_structure(),
    atom_radius    = 3,
    d_min          = 2.5)
  counter = count_evaluations()
  try:
    x = fmodel.xray_structure.extract_u_iso_or_u_equiv()*79.
    c = adp.tg(fmodel = fmodel, x = x, restraints_weight = 1,
      bound_flags = flex.int(x.size(), 2),
      lower_bound = flex.double(x.size(), 0),
      upper_bound = flex.double(x.size(), 200))
  finally:
    counter.restore()
  assert counter.per_minimization == [1], counter.per_minimization
  assert c.target() is not None and c.gradients().size() == x.size()

def refine(d_min = 2.5, inputs = None, **kw):
  if(inputs is None): inputs = dict(d_min = d_min)
  mmm, b_answer = get_inputs(**inputs)
  log = StringIO()
  counter = count_evaluations()
  try:
    o = adp.ncs_aware_refinement(
      map_model_manager = mmm,
      d_min             = d_min,
      atom_radius       = 3,
      individual        = True,
      log               = log,
      **kw)
  finally:
    counter.restore()
  b = o.mmm.model().get_b_iso()
  cc = 0
  if(flex.max(b) - flex.min(b) > 1.e-6):
    cc = flex.linear_correlation(x = b, y = b_answer).coefficient()
  return group_args(
    cc     = cc,
    diff   = flex.mean(flex.abs(b - b_answer)),
    offset = flex.mean(b - b_answer),
    n      = counter.per_minimization)

def exercise_scattering_table():
  """
  The model is in a box of its own when it is compared with the map: it
  keeps its scattering table there (electron in this test; it used to be
  replaced with the X-ray one), so the model the map was calculated from fits
  the data made from this map.
  """
  mmm, b_answer = get_inputs()
  xrs = mmm.model().get_xray_structure().deep_copy_scatterers()
  assert xrs.scattering_type_registry().last_table() == "electron"
  xrs.set_b_iso(values = b_answer)
  fmodel = adp.map_and_model_to_fmodel(
    map_data       = mmm.map_data().deep_copy(),
    xray_structure = xrs,
    atom_radius    = 3,
    d_min          = 2.5,
    reset_adp      = False)
  table = fmodel.xray_structure.scattering_type_registry().last_table()
  assert table == "electron", table
  print("r_work, model the map was calculated from: %6.4f"%fmodel.r_work())
  assert fmodel.r_work() < 0.06, fmodel.r_work() # was 0.107

def exercise_uniform_start():
  """
  Start with all B equal and no group refinement before the individual one.
  Restraints have zero gradient then, which used to set the weight of the
  data to zero (nothing refined); without restraints the minimizer used to
  stop after a few evaluations.
  """
  for rw in [1, None]:
    r = refine(restraints_weight = rw, group_mode = None)
    print("uniform start, restraints_weight=%s: CC %6.4f mean |B-B answer| "
      "%5.2f offset %5.2f evaluations %s"%(str(rw), r.cc, r.diff, r.offset,
      str(r.n)))
    assert r.cc > 0.95, r.cc
    assert r.diff < 4, r.diff
    assert abs(r.offset) < 3, r.offset
    assert max(r.n) > 10, r.n

def exercise_low_resolution():
  """
  Lower resolution (4 and 6 A), maps that fade with resolution the way
  experimental ones do (calculated to twice the resolution used, with large
  B). Density of an atom is wide there and the mask around atoms cuts much
  of it off; with the map outside the mask taken from the model
  (adp.fill_outside_mask) B-factors are close to the answer, without it
  they are not.
  """
  for d_min, b_add, limit in [(4., 100, 5), (6., 300, 20)]:
    diff = {}
    for fill in [True, False]:
      adp.fill_outside_mask = fill
      try:
        r = refine(
          d_min             = d_min,
          inputs            = dict(d_min = d_min/2, b_add = b_add,
                                   b_start = 57 + b_add, buffer_layer = 10),
          restraints_weight = 1,
          group_mode        = "one_per_residue")
      finally:
        adp.fill_outside_mask = True
      print("d_min %3.1f fill_outside_mask=%-5s: CC %6.4f mean |B-B answer| "
        "%6.2f offset %6.2f"%(d_min, str(fill), r.cc, r.diff, r.offset))
      diff[fill] = r.diff
      if(fill):
        assert r.cc > 0.9, r.cc
        assert r.diff < limit, r.diff
    assert diff[True] < 0.6*diff[False], diff

def exercise_gradients():
  """
  Gradients given to the minimizer are those of its target with respect to
  its variables (B): finite differences. (They used to be with respect to
  U, 79 times larger.)
  """
  mmm, b_answer = get_inputs()
  fmodel = adp.map_and_model_to_fmodel(
    map_data       = mmm.map_data().deep_copy(),
    xray_structure = mmm.model().get_xray_structure(),
    atom_radius    = 3,
    d_min          = 2.5)
  flex.set_random_seed(0)
  n = fmodel.xray_structure.scatterers().size()
  for rw in [1, None]:
    x = 40 + flex.random_double(n)*30
    c = adp.tg(fmodel = fmodel, x = x, restraints_weight = rw,
      bound_flags = flex.int(n, 2), lower_bound = flex.double(n, 0),
      upper_bound = flex.double(n, 200))
    g = c.gradients().deep_copy()
    direction = g*(1./g.norm())
    eps = 0.5
    c.update(x + direction*eps); tp = c.target()
    c.update(x - direction*eps); tm = c.target()
    fd = (tp - tm)/(2*eps)
    an = flex.sum(g*direction)
    print("restraints_weight=%s: gradient along itself: analytical %g finite "
      "differences %g"%(str(rw), an, fd))
    assert abs(fd - an) < 0.05*abs(an), [fd, an]

def exercise_group_only():
  """
  Group B refinement without individual one (what is done at 3.5 A and
  lower resolution): rounds stop when one no longer lowers the R-factor (all
  five used to be done, the last four changing R in the fourth digit at five
  times the cost), and the result is close to the answer.
  """
  mmm, b_answer = get_inputs(d_min = 2., b_add = 100, b_start = 157,
    buffer_layer = 10)
  log = StringIO()
  o = adp.ncs_aware_refinement(
    map_model_manager = mmm,
    d_min             = 4.,
    atom_radius       = 3,
    individual        = False,
    group_mode        = "one_per_residue",
    log               = log)
  b = o.mmm.model().get_b_iso()
  cc = flex.linear_correlation(x = b, y = b_answer).coefficient()
  diff = flex.mean(flex.abs(b - b_answer))
  rounds = [l for l in log.getvalue().splitlines()
            if l.count("group b-factor refinement (macro cycle =") and
               not l.count("(macro cycle = 0;")]
  print("group only: rounds %d CC %6.4f mean |B - B answer| %5.2f"%(
    len(rounds), cc, diff))
  assert 1 <= len(rounds) <= 3, len(rounds)
  assert cc > 0.85, cc
  assert diff < 8, diff

def exercise_hydrogens(d_min = 2.5):
  """
  Hydrogens ride on the atoms they are bonded to: B(H) = 1.2 * B(parent),
  they are not refined as atoms of their own (which gave B = 0 to many of
  them: H are invisible in the map, and they are not restrained either).
  """
  import mmtbx.hydrogens.reduce_hydrogen
  from libtbx.utils import null_out
  ph = ssb.secondary_structure_from_sequence(
    ssb.alpha_helix_str, "ACDEFGHIKLMN")
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = ph.atoms().extract_xyz(), buffer_layer = 5)
  ph.atoms().set_xyz(box.sites_cart)
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = ph,
    crystal_symmetry = box.crystal_symmetry(),
    log              = null_out())
  o = mmtbx.hydrogens.reduce_hydrogen.place_hydrogens(model = model)
  o.run()
  model = o.get_model()
  model.set_log(null_out())
  model.process(make_restraints = True)
  model.setup_scattering_dictionaries(scattering_table = "electron")
  hd = model.get_hd_selection()
  assert hd.count(True) > 50
  pairs = model.xh_connectivity_table()
  assert len(pairs) == hd.count(True)
  # parents are the same found from restraints and, without restraints, as
  # the nearest non-hydrogen atom of the residue
  model_no_restraints = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = model.get_hierarchy().deep_copy(),
    crystal_symmetry = model.crystal_symmetry(),
    log              = null_out())
  assert not model_no_restraints.restraints_manager_available()
  r1 = adp.hydrogen_riding(model = model)
  r2 = adp.hydrogen_riding(model = model_no_restraints)
  assert r1.scale == 1.2
  assert sorted(zip(r1.i_x, r1.i_h)) == sorted(zip(r2.i_x, r2.i_h))
  assert sorted(zip(r1.i_x, r1.i_h)) == sorted([(p[0],p[1]) for p in pairs])
  assert adp.hydrogen_riding(model = model.select(~hd)) is None
  def ride(b):
    for p in pairs: b[p[1]] = 1.2*b[p[0]]
    return b
  b_answer = flex.double()
  for i, rg in enumerate(model.get_hierarchy().residue_groups()):
    for a in rg.atoms():
      b = 30 + 5*abs(i-6)
      if(not a.name.strip() in ["N", "CA", "C", "O"]): b += 10
      b_answer.append(b)
  b_answer = ride(b_answer)
  model.set_b_iso(values = b_answer)
  fc = model.get_xray_structure().structure_factors(d_min = d_min).f_calc()
  map_data = fc.fft_map(resolution_factor = 0.25).real_map_unpadded()
  mm = map_manager.map_manager(
    map_data                   = map_data,
    unit_cell_grid             = map_data.all(),
    unit_cell_crystal_symmetry = model.crystal_symmetry(),
    wrapping                   = False)
  model.set_b_iso(values = flex.double(b_answer.size(), 50))
  mmm = map_model_manager.map_model_manager(model = model, map_manager = mm)
  o = adp.ncs_aware_refinement(
    map_model_manager = mmm,
    d_min             = d_min,
    atom_radius       = 3,
    individual        = True,
    log               = StringIO())
  b = o.mmm.model().get_b_iso()
  d = flex.double([abs(b[p[1]] - 1.2*b[p[0]]) for p in pairs])
  cc = flex.linear_correlation(
    x = b.select(~hd), y = b_answer.select(~hd)).coefficient()
  print("hydrogens: max |B(H) - 1.2*B(parent)| %6.3f, min B(H) %5.1f, "
    "CC(B, answer) non-H %6.4f"%(flex.max(d), flex.min(b.select(hd)), cc))
  assert flex.max(d) < 0.01, flex.max(d)
  assert flex.min(b.select(hd)) > 10
  assert cc > 0.8, cc

if(__name__ == "__main__"):
  t0 = time.time()
  exercise_hydrogens()
  exercise_group_only()
  exercise_gradients()
  exercise_low_resolution()
  exercise_scattering_table()
  exercise_uniform_start()
  exercise()
  exercise_tg_evaluates_once()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
