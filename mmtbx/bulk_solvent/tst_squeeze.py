from cctbx.array_family import flex
from cctbx import crystal, xray, adptbx
from mmtbx.bulk_solvent import squeeze
from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
import mmtbx.f_model
from libtbx import group_args
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import random, time, math

"""
Squeeze (mmtbx.bulk_solvent.squeeze): simulated data with man-made blobs of
density of different size in the solvent region that the model does not have.
"""

D_MIN = 1.0
# Blobs: number of atoms, radius of the blob (A), occupancy, B (A^2)
BLOBS = [(1, 0.0, 0.8, 12.), (4, 1.3, 0.7, 18.), (12, 2.2, 0.6, 25.)]

def scatterer(label, site_frac, b, occupancy, element):
  return xray.scatterer(label = label, site = site_frac,
    u = adptbx.b_as_u(b), occupancy = occupancy, scattering_type = element)

def distance_to(uc, site_frac, others_frac):
  result = 1.e9
  for o in others_frac:
    result = min(result, uc.mod_short_distance(site_frac, o))
  return result

def synthetic(d_min=D_MIN, blobs=BLOBS, n_water=6):
  """
  A globule of atoms (model), a few water molecules (model), and blobs of
  atoms in the solvent (in the data only). F_obs = |F_model| of everything
  with a flat bulk solvent.

  Returns group_args: f_obs, r_free_flags, xray_structure (model),
  water_selection, blob_sites_frac (list of flex.vec3_double, one per blob),
  f_blobs (flex.complex_double, true contribution of the blob atoms),
  control_sites_frac (points in the solvent away from everything).
  """
  cs = crystal.symmetry((32, 35, 37, 90, 90, 90), "P 21 21 21")
  uc = cs.unit_cell()
  sp = crystal.special_position_settings(cs)
  # the globule: a jittered lattice inside a sphere
  center, radius, spacing = (7., 8., 8.5), 6.0, 1.9
  atoms = []
  n = int(radius/spacing)+1
  for i in range(-n, n+1):
    for j in range(-n, n+1):
      for k in range(-n, n+1):
        d = [i*spacing, j*spacing, k*spacing]
        if(math.sqrt(d[0]**2+d[1]**2+d[2]**2) > radius): continue
        site = [center[m]+d[m]+random.uniform(-0.25, 0.25) for m in range(3)]
        atoms.append(scatterer("C%d"%len(atoms), uc.fractionalize(site),
          b = random.uniform(6, 10), occupancy = 1,
          element = random.choice(["C", "C", "N", "O"])))
  model = xray.structure(sp, flex.xray_scatterer(atoms))
  # centres of all symmetry copies of the globule
  centers = list(xray.structure(sp, flex.xray_scatterer([scatterer("X",
    uc.fractionalize(center), 10, 1, "O")])).expand_to_p1().sites_frac())
  def distance_to_globule(s):
    return distance_to(uc, s, centers) - radius
  def free_site(d_atoms, d_others, others):
    for trial in range(100000):
      s = (random.random(), random.random(), random.random())
      if(distance_to_globule(s) < d_atoms): continue
      if(len(others) > 0 and distance_to(uc, s, others) < d_others): continue
      return s
    raise RuntimeError("no place found")
  # all symmetry copies of what is placed already
  taken = []
  def register(sites_frac):
    x = xray.structure(sp, flex.xray_scatterer(
      [scatterer("X", s, 10, 1, "O") for s in sites_frac]))
    taken.extend(list(x.expand_to_p1().sites_frac()))
  # blobs, well inside the solvent
  blob_sites, blob_atoms = [], []
  for n_atoms, r_blob, occ, b in blobs:
    c = uc.orthogonalize(free_site(4.5+r_blob, 7.0+r_blob, taken))
    sites = flex.vec3_double()
    while sites.size() < n_atoms:
      d = [random.uniform(-r_blob, r_blob) for m in range(3)]
      if(math.sqrt(d[0]**2+d[1]**2+d[2]**2) > r_blob): continue
      s = uc.fractionalize([c[m]+d[m] for m in range(3)])
      if(sites.size() > 0 and distance_to(uc, s, list(sites)) < 1.0 and
         n_atoms < 8): continue
      sites.append(s)
      blob_atoms.append(scatterer("B", s, b, occ, "O"))
    blob_sites.append(sites)
    register(list(sites))
  # water molecules of the model: at the surface of the globule
  waters = []
  while len(waters) < n_water:
    s = free_site(2.5, 5.0, taken)
    if(distance_to_globule(s) > 3.2): continue
    waters.append(scatterer("O", s, 15, 1, "O"))
    register([s])
  # control points: in the solvent, far from everything
  control = flex.vec3_double()
  while control.size() < 5:
    s = free_site(4.5, 6.0, taken)
    control.append(s)
    register([s])
  n_model = len(atoms)
  model = xray.structure(sp, flex.xray_scatterer(atoms + waters))
  full  = xray.structure(sp, flex.xray_scatterer(atoms + waters + blob_atoms))
  water_selection = flex.bool(n_model, False)
  water_selection.extend(flex.bool(len(waters), True))
  f_dummy = abs(full.structure_factors(d_min = d_min).f_calc())
  fm = mmtbx.f_model.manager(f_obs = f_dummy, xray_structure = full,
    k_sol = 0.35, b_sol = 46.)
  f_obs = abs(fm.f_model())
  flags = f_obs.generate_r_free_flags(fraction = 0.05, max_free = None)
  f_blobs = f_obs.structure_factors_from_scatterers(
    xray_structure = xray.structure(sp, flex.xray_scatterer(blob_atoms)),
    algorithm = "direct").f_calc().data()
  return group_args(
    f_obs              = f_obs,
    r_free_flags       = flags,
    xray_structure     = model,
    water_selection    = water_selection,
    blob_sites_frac    = blob_sites,
    f_blobs            = f_blobs,
    control_sites_frac = control)

def write_files(o, pdb_file_name, mtz_file_name):
  """
  Model (water as HOH) and data (FOBS, R-free-flags) of synthetic() as
  files.
  """
  import iotbx.pdb
  xrs = o.xray_structure
  uc = xrs.unit_cell()
  lines = [iotbx.pdb.format_cryst1_record(crystal_symmetry = xrs)]
  for i, sc in enumerate(xrs.scatterers()):
    x, y, z = uc.orthogonalize(sc.site)
    e = sc.scattering_type
    record, name, resname = "ATOM  ", " %-3s"%(e+"1"), "GLY"
    if(o.water_selection[i]): record, name, resname = "HETATM", " O  ", "HOH"
    lines.append("%s%5d %s %s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"%(
      record, i+1, name, resname, i+1, x, y, z, sc.occupancy,
      adptbx.u_as_b(sc.u_iso), e))
  with open(pdb_file_name, "w") as f: f.write("\n".join(lines)+"\nEND\n")
  mtz = o.f_obs.as_mtz_dataset(column_root_label = "FOBS")
  mtz.add_miller_array(miller_array = o.r_free_flags,
    column_root_label = "R-free-flags")
  mtz.mtz_object().write(file_name = mtz_file_name)

def new_fmodel(o, scale=True):
  fm = mmtbx.f_model.manager(
    f_obs          = o.f_obs,
    r_free_flags   = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers())
  if(scale):
    fm.update_all_scales(remove_outliers = False,
      apply_scale_k1_to_f_obs = False)
  return fm

def density_at(map_data, sites_frac):
  """Largest map value at the given sites."""
  return max([map_data.tricubic_interpolation(s) for s in sites_frac])

def correlation(a, b):
  """Correlation of two complex arrays (as vectors)."""
  num = flex.sum(flex.real(a*flex.conj(b)))
  den = math.sqrt(flex.sum(flex.norm(a))*flex.sum(flex.norm(b)))
  return num/den

def r_factor(f_obs, f_model):
  fm = flex.abs(f_model)
  k = flex.sum(f_obs*fm)/flex.sum(fm*fm)
  return flex.sum(flex.abs(f_obs-k*fm))/flex.sum(f_obs)

def exercise_fit_kb():
  """
  k and B are recovered from amplitudes.
  """
  n = 5000
  ss = flex.random_double(n)*0.25
  f0 = flex.complex_double(flex.random_double(n)-0.5,
                           flex.random_double(n)-0.5)*20
  f  = flex.complex_double(flex.random_double(n)-0.5,
                           flex.random_double(n)-0.5)*5
  k_total = 0.5+flex.random_double(n)
  params = squeeze.default_params()
  for k, b in [(3.0, 22.), (0.7, 5.), (1.5, -10.), (10., 80.)]:
    f_obs = flex.abs(k_total*(f0 + f*(k*flex.exp(-b*ss))))
    k_, b_ = squeeze.fit_kb(f_obs = f_obs, k_total = k_total, f0 = f0, f = f,
      ss = ss, params = params)
    assert approx_equal([k_, b_], [k, b], 1.e-3), [k_, b_, k, b]
  # bounds
  f_obs = flex.abs(k_total*(f0 + f*(3*flex.exp(-22*ss))))
  params.k_max = 2
  params.b_min, params.b_max = 30, 40
  params.b_start = 35
  k_, b_ = squeeze.fit_kb(f_obs = f_obs, k_total = k_total, f0 = f0, f = f,
    ss = ss, params = params)
  assert k_ <= 2 and 30 <= b_ <= 40, [k_, b_]

def exercise_run(o):
  fm = new_fmodel(o)
  f_model_before = fm.f_model().data().deep_copy()
  params = squeeze.default_params()
  gridding = squeeze.crystal_gridding(fmodel = fm, params = params)
  step = max([a/n for a, n in zip(fm.f_obs().unit_cell().parameters()[:3],
    gridding.n_real())])
  assert step <= 0.2 + 1.e-6, step
  # the blobs are in the map, and they are in the solvent mask
  m = squeeze.difference_map(fmodel = fm, gridding = gridding)
  sigma = m.sample_standard_deviation()
  mask = squeeze.solvent_mask(fmodel = fm, params = params,
    n_real = gridding.n_real(), water_selection = o.water_selection)
  for sites in o.blob_sites_frac:
    assert density_at(m, sites) > 5*sigma, density_at(m, sites)/sigma
    for s in sites: assert mask.value_at_closest_grid_point(s) == 1
  # water: part of the region with the mask computed without it only
  i_water = o.water_selection.iselection()[0]
  s = o.xray_structure.sites_frac()[i_water]
  assert mask.value_at_closest_grid_point(s) == 1
  mask_w = squeeze.solvent_mask(fmodel = fm, params = params,
    n_real = gridding.n_real(), water_selection = None)
  assert mask_w.value_at_closest_grid_point(s) == 0
  assert flex.mean(mask_w.as_1d()) < flex.mean(mask.as_1d())
  params.exclude_water_from_mask = False
  mask_w2 = squeeze.solvent_mask(fmodel = fm, params = params,
    n_real = gridding.n_real(), water_selection = o.water_selection)
  assert mask_w2.as_1d().all_eq(mask_w.as_1d())
  # run
  sio = StringIO()
  r = squeeze.run(fmodel = fm, water_selection = o.water_selection,
    log = sio)
  assert r is not None, sio.getvalue()
  assert approx_equal(fm.f_model().data(), f_model_before) # not changed
  assert r.f.size() == fm.f_obs().data().size()
  assert r.threshold == 2.5 and r.water_excluded
  assert approx_equal(r.sigma, sigma)
  assert r.k > 1 and r.volume > 0 and r.n_electrons > 0
  assert 0.5 < r.solvent_fraction < 0.9, r.solvent_fraction
  # it is the contribution of the blob atoms that is missing in the model
  f_sq = r.scaled(ss = fm.ss)
  cc = correlation(f_sq, o.f_blobs)
  assert cc > 0.85, cc
  # and F_model with it fits the data better: work and free reflections
  k_total = fm.k_isotropic()*fm.k_anisotropic()
  f_new = f_model_before + f_sq*k_total
  for sel in [~o.r_free_flags.data(), o.r_free_flags.data()]:
    fo = fm.f_obs().data().select(sel)
    r0 = r_factor(fo, f_model_before.select(sel))
    r1 = r_factor(fo, f_new.select(sel))
    assert r1 < r0*0.8, [r0, r1]
  # show, summary
  assert r.summary().startswith("squeeze: density above 2.5 sigma")
  sio = StringIO()
  r.show(log = sio)
  assert "solvent mask computed without water" in sio.getvalue()
  assert "scale k:" in sio.getvalue()
  # select
  sel = flex.random_bool(r.f.size(), 0.5)
  rs = r.select(sel)
  assert approx_equal(rs.f, r.f.select(sel))
  assert rs.k == r.k and rs.b == r.b
  return r

def exercise_settings(o, r_default):
  fm = new_fmodel(o)
  # higher threshold: less density
  params = squeeze.default_params()
  params.threshold = 6
  r = squeeze.run(fmodel = fm, params = params,
    water_selection = o.water_selection)
  assert 0 < r.volume < r_default.volume
  # nothing above the threshold
  params.threshold = 1000
  sio = StringIO()
  r = squeeze.run(fmodel = fm, params = params, log = sio)
  assert r is None
  assert "no density above" in sio.getvalue()
  # without water selection: mask from all atoms
  r = squeeze.run(fmodel = fm)
  assert not r.water_excluded
  assert r.solvent_fraction < r_default.solvent_fraction
  # bad settings
  for name, value in [("threshold", 0), ("number_of_cycles", 0),
                      ("grid_step_max", -1), ("k_max", 0), ("b_min", 1000)]:
    params = squeeze.default_params()
    setattr(params, name, value)
    try: squeeze.run(fmodel = fm, params = params)
    except Exception as e: assert "squeeze" in str(e), str(e)
    else: raise RuntimeError("Exception expected")

def exercise_not_applicable(o):
  fm = new_fmodel(o)
  assert squeeze.not_applicable(fmodel = fm) is None
  # resolution
  fm_low = fm.resolution_filter(d_min = 1.5)
  reason = squeeze.not_applicable(fmodel = fm_low)
  assert reason == "data resolution 1.50 A, 1.10 A or better needed", reason
  sio = StringIO()
  assert squeeze.run(fmodel = fm_low, log = sio) is None
  assert "squeeze: skipped (data resolution" in sio.getvalue()
  params = squeeze.default_params()
  params.d_min_max = None
  assert squeeze.not_applicable(fmodel = fm_low, params = params) is None
  assert squeeze.run(fmodel = fm_low, params = params) is not None
  # no atoms
  fm_no = mmtbx.f_model.manager(
    f_obs        = o.f_obs,
    r_free_flags = o.r_free_flags,
    f_calc       = fm.f_calc(),
    f_mask       = fm.f_masks()[0])
  assert squeeze.not_applicable(fmodel = fm_no) == "no atomic model"
  # neutron
  xrs = o.xray_structure.deep_copy_scatterers()
  xrs.switch_to_neutron_scattering_dictionary()
  fm_n = mmtbx.f_model.manager(f_obs = o.f_obs,
    r_free_flags = o.r_free_flags, xray_structure = xrs)
  assert squeeze.not_applicable(fmodel = fm_n) == "neutron data"
  # twinned
  fm_t = mmtbx.f_model.manager(f_obs = o.f_obs,
    r_free_flags = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers(),
    twin_law = "k,h,-l")
  assert squeeze.not_applicable(fmodel = fm_t) == "twinned data"

def exercise_params():
  p = bss.master_params.extract().squeeze
  assert p.threshold == 2.5
  assert p.exclude_water_from_mask is True
  assert p.d_min_max == 1.1
  assert p.number_of_cycles == 2

def run():
  random.seed(0)
  flex.set_random_seed(0)
  exercise_params()
  exercise_fit_kb()
  o = synthetic()
  r = exercise_run(o)
  exercise_settings(o, r)
  exercise_not_applicable(o)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
