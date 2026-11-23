from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import crystal, xray, adptbx
from scitbx import matrix
from mmtbx.bulk_solvent import mosaic
import mmtbx.masks
from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
import mmtbx.f_model
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.test_utils import Exception_expected
from six.moves import cStringIO as StringIO
import math
from libtbx.test_utils import approx_equal
import boost_adaptbx.boost.python as bp
import random, time
mosaic_ext = bp.import_ext("mmtbx_mosaic_ext")

BIG = 1.e30

def random_complex(n):
  return flex.complex_double(flex.random_double(n)-0.5,
                             flex.random_double(n)-0.5)

# ------------------------------------------------------------------------------
# C++ kernels
# ------------------------------------------------------------------------------

def exercise_normal_matrix():
  G = [random_complex(200) for i in range(4)]
  A = mosaic_ext.normal_matrix(G)
  assert A.focus() == (4,4)
  for n in range(4):
    for m in range(4):
      ref = flex.sum(flex.real(G[n]*flex.conj(G[m])))
      assert approx_equal(A[n*4+m], ref, 1.e-9)
  A = mosaic_ext.normal_matrix([])
  assert A.size() == 0

def gradient(A, b, x, n):
  return [sum([A[i*n+j]*x[j] for j in range(n)]) - b[i] for i in range(n)]

def exercise_solve_bounded():
  n = 6
  # symmetric positive definite matrix
  G = [random_complex(50) for i in range(n)]
  A = mosaic_ext.normal_matrix(G)
  x_true = flex.double([0.3, -0.7, 1.4, 0.0, 2.5, -1.2])
  b = flex.double(gradient(A, [0]*n, x_true, n))
  zero = flex.double(n, 0)
  # (a) bounds not active: solution of A x = b
  x = mosaic_ext.alg4_solve_bounded(a=A, b=b, lower=flex.double(n, -BIG),
    upper=flex.double(n, BIG), x_start=zero)
  assert approx_equal(x, x_true, 1.e-6)
  x_inv = matrix.sqr(list(A)).inverse() * matrix.col(list(b))
  assert approx_equal(x, x_inv, 1.e-6)
  # the start does not matter
  x = mosaic_ext.alg4_solve_bounded(a=A, b=b, lower=flex.double(n, -BIG),
    upper=flex.double(n, BIG), x_start=flex.double(n, 100))
  assert approx_equal(x, x_true, 1.e-6)
  # (b) bounds active: inside the box, and the Karush-Kuhn-Tucker conditions
  for lo, hi in [(0., 1.), (0., BIG), (-BIG, 1.), (-0.5, 0.5)]:
    lower, upper = flex.double(n, lo), flex.double(n, hi)
    x = mosaic_ext.alg4_solve_bounded(a=A, b=b, lower=lower, upper=upper,
      x_start=zero)
    g = gradient(A, b, x, n)
    n_active = 0
    for xi, gi in zip(x, g):
      assert lo <= xi <= hi
      if(xi == lo):   assert gi >= -1.e-6; n_active += 1
      elif(xi == hi): assert gi <=  1.e-6; n_active += 1
      else:           assert abs(gi) < 1.e-6
    if((lo, hi) != (-0.5, 0.5)): assert 0 < n_active < n
    # not the same as clipping the unconstrained solution
    clipped = flex.double([min(hi, max(lo, v)) for v in x_true])
    def target(v):
      return sum([0.5*v[i]*sum([A[i*n+j]*v[j] for j in range(n)]) - b[i]*v[i]
        for i in range(n)])
    assert target(x) <= target(clipped) + 1.e-9
  # (c) empty column (region that does not contribute): no NaN, value 0
  G0 = G[:3] + [flex.complex_double(50, 0)] + G[4:]
  A0 = mosaic_ext.normal_matrix(G0)
  x = mosaic_ext.alg4_solve_bounded(a=A0, b=b, lower=flex.double(n, -1),
    upper=flex.double(n, BIG), x_start=flex.double(n, 0.5))
  for xi in x: assert xi == xi and abs(xi) < 1.e6
  assert x[3] == 0.5 # not determined by the data: stays where it started
  # (d) two identical columns (singular matrix): finite, and reproduces b
  G1 = G[:5] + [G[4]]
  A1 = mosaic_ext.normal_matrix(G1)
  b1 = flex.double(gradient(A1, [0]*n, [1,1,1,1,1,1], n))
  x = mosaic_ext.alg4_solve_bounded(a=A1, b=b1, lower=flex.double(n, 0),
    upper=flex.double(n, BIG), x_start=zero)
  assert approx_equal(gradient(A1, b1, x, n), [0]*n, 1.e-5)
  assert approx_equal(x[4]+x[5], 2, 1.e-4)

def exercise_alg4_kernel():
  """
  F_obs = |base + sum x_n G_n| with known x: recovered from amplitudes alone
  starting from x=0; bounds respected.
  """
  size, n = 2000, 4
  base = random_complex(size)*20
  G = [random_complex(size) for i in range(n)]
  for x_true in [[0.2, 0.9, -0.6, 0.0], [-1.0, 1.5, 0.3, -0.2]]:
    m = base.deep_copy()
    for xn, g in zip(x_true, G): m += g*xn
    f_obs = flex.abs(m)
    x = mosaic_ext.alg4(G=G, f_obs=f_obs, base=base,
      lower=flex.double(n, -BIG), upper=flex.double(n, BIG),
      x_start=flex.double(n, 0), n_cycles=50)
    assert approx_equal(x, x_true, 1.e-3)
    x = mosaic_ext.alg4(G=G, f_obs=f_obs, base=base,
      lower=flex.double(n, -0.5), upper=flex.double(n, 0.5),
      x_start=flex.double(n, 0), n_cycles=20)
    assert flex.min(x) >= -0.5 and flex.max(x) <= 0.5
    for xi, ti in zip(x, x_true):
      if(ti >  0.5): assert xi ==  0.5
      if(ti < -0.5): assert xi == -0.5
  # nothing to fit
  x = mosaic_ext.alg4(G=[], f_obs=f_obs, base=base, lower=flex.double(),
    upper=flex.double(), x_start=flex.double())
  assert x.size() == 0

def exercise_alg2_tg_update():
  """
  Target of the unphased search does not depend on the history of calls.
  """
  size, n = 300, 3
  F = [random_complex(size) for i in range(n)]
  i_obs = flex.random_double(size)+1
  x1 = flex.double([1.0, 0.3, 0.7])
  x2 = flex.double([0.9, 0.1, 0.2])
  o = mosaic_ext.alg2_tg(F=F, i_obs=i_obs)
  o.update(x1)
  t1, g1 = o.target(), o.gradient().deep_copy()
  o.update(x2)
  o.update(x1)
  assert approx_equal(o.target(), t1, 1.e-12)
  assert approx_equal(o.gradient(), g1, 1.e-12)
  fresh = mosaic_ext.alg2_tg(F=F, i_obs=i_obs)
  fresh.update(x1)
  assert approx_equal(fresh.target(), t1, 1.e-12)

# ------------------------------------------------------------------------------
# Synthetic model: sealed spherical shells of atoms, each around an empty
# cavity that the bulk-solvent mask fills.
# ------------------------------------------------------------------------------

# Positions chosen so that no two shells (symmetry copies included) are closer
# than 10 A, surface to surface.
SHELLS = [((18., 9., 56.), 9.0), ((16., 26., 14.), 7.5), ((32., 17., 34.), 6.5)]
CELL = (60., 62., 64., 90, 90, 90)

def shell_sites(center, radius, spacing=1.5):
  n = int(4*math.pi*radius**2/spacing**2*1.3)
  sites = []
  ga = math.pi*(3.-math.sqrt(5.))
  for i in range(n):
    z = 1 - 2*(i+0.5)/n
    r = math.sqrt(max(0, 1-z*z))
    phi = i*ga
    sites.append((center[0]+radius*r*math.cos(phi),
                  center[1]+radius*r*math.sin(phi),
                  center[2]+radius*z))
  return sites

def shell_model(shells=SHELLS, space_group_symbol="P 21 21 21", extra=None):
  """
  xray structure of C atoms: closed shells around the given centers, and
  optional extra sites (list of Cartesian coordinates).
  """
  cs = crystal.symmetry(CELL, space_group_symbol)
  sites = []
  for center, radius in shells: sites.extend(shell_sites(center, radius))
  if(extra is not None): sites.extend(extra)
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer()
  uc = cs.unit_cell()
  for i, site in enumerate(sites):
    scatterers.append(xray.scatterer(label="C%d"%i, site=uc.fractionalize(site),
      u=adptbx.b_as_u(20), scattering_type="C"))
  return xray.structure(sp, scatterers)

def cavity_of(regs, center):
  """Secondary region that contains the given Cartesian point."""
  nx,ny,nz = regs.n_real
  f = regs.unit_cell.fractionalize(center)
  p = ((int(round(f[0]*nx))%nx)*ny + int(round(f[1]*ny))%ny)*nz + \
       int(round(f[2]*nz))%nz
  found = [r for r in regs.secondary if (r.points == p).count(True) == 1]
  assert len(found) == 1, len(found)
  return found[0]

def exercise_regions():
  mp = mmtbx.masks.mask_master_params.extract()
  xrs = shell_model()
  regs = mosaic.regions(xray_structure = xrs, mask_params = mp)
  n_grid = regs.n_real[0]*regs.n_real[1]*regs.n_real[2]
  assert regs.n_grid == n_grid
  assert approx_equal(regs.volume_per_point, xrs.unit_cell().volume()/n_grid)
  # three cavities, each with its four symmetry copies merged into one region
  assert len(regs.secondary) == 3, len(regs.secondary)
  assert regs.main.n_points > 10*regs.secondary[0].n_points
  volumes = [r.volume for r in regs.secondary]
  assert volumes == sorted(volumes, reverse=True)
  # mask boundary: atom centre - (vdW radius + solvent radius - shrink radius)
  margin = 1.775 + mp.solvent_radius - mp.shrink_truncation_radius
  for (center, radius), r in zip(SHELLS, regs.secondary):
    assert cavity_of(regs, center) is r
    expected = 4*4./3*math.pi*(radius-margin)**3
    assert abs(r.volume-expected)/expected < 0.2, [r.volume, expected]
    assert approx_equal(r.volume, r.n_points*regs.volume_per_point)
    assert r.points.size() == r.n_points
    assert r.n_points % 4 == 0
    # the reported centre is the centre of one of the copies
    uc = xrs.unit_cell()
    cf = uc.fractionalize(center)
    d = min([uc.mod_short_distance(op*cf, r.centroid_frac)
      for op in xrs.space_group().all_ops()])
    assert d < 1.0, d
  # every grid point is protein, main solvent or secondary solvent
  n_sec = sum([r.n_points for r in regs.secondary])
  n_solvent = regs.main.n_points + n_sec
  assert approx_equal(regs.solvent_fraction, n_solvent/n_grid)
  assert 0.5 < regs.solvent_fraction < 0.95
  # masks of regions
  m = regs.mask_of(regs.secondary)
  assert m.all() == regs.n_real
  assert approx_equal(flex.sum(m.as_1d()), n_sec)
  assert approx_equal(flex.sum(regs.mask_of([regs.secondary[2]]).as_1d()),
    regs.secondary[2].n_points)
  # volumes do not depend on the grid step (the old code assumed step 0.6)
  mp.step = 0.4
  regs4 = mosaic.regions(xray_structure = xrs, mask_params = mp)
  assert regs4.n_real[0] > regs.n_real[0]
  assert len(regs4.secondary) == 3
  for r6, r4 in zip(regs.secondary, regs4.secondary):
    assert abs(r6.volume-r4.volume)/r6.volume < 0.1, [r6.volume, r4.volume]
  # P1: no merging, 3 regions again
  regs_p1 = mosaic.regions(
    xray_structure = shell_model(space_group_symbol="P 1"),
    mask_params    = mmtbx.masks.mask_master_params.extract())
  assert len(regs_p1.secondary) == 3
  for r, rp in zip(regs.secondary, regs_p1.secondary):
    assert abs(r.volume/4-rp.volume)/rp.volume < 0.05
  # no cavities
  xrs0 = shell_model(shells=[((18., 9., 56.), 3.0)])
  regs0 = mosaic.regions(xray_structure = xrs0,
    mask_params = mmtbx.masks.mask_master_params.extract())
  assert len(regs0.secondary) == 0

def exercise_padded_maps():
  """
  Maps with the padding of the in-place Fourier transform: values land on
  the same grid points as without padding, for even and odd grid sizes.
  """
  from cctbx import maptbx
  for n_real in [(4, 6, 8), (5, 6, 9), (3, 5, 7), (6, 4, 2)]:
    regs = mosaic.regions.__new__(mosaic.regions)
    regs.n_real = n_real
    n = n_real[0]*n_real[1]*n_real[2]
    all_points = list(range(n))
    random.shuffle(all_points)
    p1 = flex.size_t(sorted(all_points[:n//3]))
    p2 = flex.size_t(sorted(all_points[n//3:n//2]))
    for points, values in [([p1, p2], [0.5, -2.]), ([p1], None),
                           ([flex.size_t(), p2], [3., 1.]),
                           ([flex.size_t([0, n-1])], [7.])]:
      m = regs.mask_of_points(points = points, values = values)
      mp = regs.mask_of_points(points = points, values = values,
        padded = True)
      assert not m.is_padded() and m.all() == n_real
      assert mp.is_padded() and mp.focus() == n_real
      assert mp.all()[2] == 2*(n_real[2]//2+1)
      maptbx.unpad_in_place(map = mp)
      assert mp.all() == n_real
      assert approx_equal(mp.as_1d(), m.as_1d())
      expected = flex.double(n, 0)
      for p, v in zip(points, values if values is not None else [1]):
        expected.set_selected(p, v)
      assert approx_equal(m.as_1d(), expected)

def exercise_region_structure_factors():
  """
  Region structure factors: direct summation and FFT agree, and all regions
  together give the F_mask that fmodel uses.
  """
  mp = mmtbx.masks.mask_master_params.extract()
  xrs = shell_model()
  f_obs = abs(xrs.structure_factors(d_min=3.0).f_calc())
  fmodel = mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs,
    mask_params = mp)
  regs = mosaic.regions(xray_structure = xrs, mask_params = mp)
  sf = mosaic.structure_factors(regions = regs, miller_array = f_obs,
    d_min = mp.Fmask_res_high)
  for r in regs.secondary:
    a = sf.of_points(r.points, method="direct")
    b = sf.of_points(r.points, method="fft")
    c = sf.of_points(r.points)
    assert a.size() == f_obs.data().size()
    assert flex.sum(flex.abs(a-b))/flex.sum(flex.abs(b)) < 1.e-9
    assert flex.sum(flex.abs(a-c))/flex.sum(flex.abs(b)) < 1.e-9
  # main region = all of the solvent but the secondary regions
  main = regs.solvent_mask()
  n_solvent = flex.sum(main.as_1d())
  for r in regs.secondary: main.as_1d().set_selected(r.points, 0)
  assert approx_equal(flex.sum(main.as_1d()), regs.main.n_points)
  assert approx_equal(n_solvent/regs.n_grid, regs.solvent_fraction)
  f_all = sf.of_map(map_data = main)
  for r in regs.secondary: f_all += sf.of_points(r.points)
  f_mask = fmodel.f_masks()[0].data()
  assert mosaic._same_mask(fmodel = fmodel, regions = regs, sf = sf)
  assert approx_equal(fmodel.mask_manager.solvent_content_via_mask,
    regs.solvent_fraction, 1.e-9)
  # another mask (different radii) is recognized as such
  mp2 = mmtbx.masks.mask_master_params.extract()
  mp2.solvent_radius = 1.4
  regs2 = mosaic.regions(xray_structure = xrs, mask_params = mp2)
  sf2 = mosaic.structure_factors(regions = regs2, miller_array = f_obs,
    d_min = mp.Fmask_res_high)
  assert not mosaic._same_mask(fmodel = fmodel, regions = regs2, sf = sf2)
  assert sf2.mask_mismatch(f_mask = fmodel.f_masks()[0]) > 1.e-3
  assert flex.sum(flex.abs(f_all-f_mask))/flex.sum(flex.abs(f_mask)) < 1.e-6
  assert sf.mask_mismatch(f_mask = fmodel.f_masks()[0]) < 1.e-6
  # zero beyond the resolution limit of the mask contribution
  d = f_obs.d_spacings().data()
  assert flex.abs(f_all.select(d < mp.Fmask_res_high)).all_eq(0)
  # anomalous set of reflections (Friedel mates are separate)
  f_anom = f_obs.generate_bijvoet_mates()
  sf_a = mosaic.structure_factors(regions = regs, miller_array = f_anom,
    d_min = mp.Fmask_res_high)
  for r in regs.secondary:
    a = sf_a.of_points(r.points, method="direct")
    b = sf_a.of_points(r.points, method="fft")
    assert a.size() == f_anom.data().size()
    assert flex.sum(flex.abs(a-b))/flex.sum(flex.abs(b)) < 1.e-9
  fmodel_a = mmtbx.f_model.manager(f_obs = f_anom, xray_structure = xrs,
    mask_params = mp)
  assert sf_a.mask_mismatch(f_mask = fmodel_a.f_masks()[0]) < 1.e-6
  # subset of reflections, values for the subset only
  subset = flex.random_bool(f_obs.data().size(), 0.7)
  sf_s = mosaic.structure_factors(regions = regs, miller_array = f_obs,
    d_min = mp.Fmask_res_high, selection = subset)
  assert sf_s.selection.all_eq(subset & (d >= mp.Fmask_res_high))
  assert sf_s.n_selected == sf_s.selection.count(True) < sf.n_selected
  r = regs.secondary[0]
  full = sf.of_points(r.points)
  for method in ["direct", "fft", None]:
    part = sf_s.of_points(r.points, method = method)
    assert part.size() == f_obs.data().size()
    assert approx_equal(part.select(sf_s.selection),
      full.select(sf_s.selection), 1.e-6)
    assert flex.abs(part.select(~sf_s.selection)).all_eq(0)
    compact = sf_s.of_points(r.points, method = method, compact = True)
    assert compact.size() == sf_s.n_selected
    assert approx_equal(compact, full.select(sf_s.selection), 1.e-6)
  # weighted sum over regions in one go
  weights = [0.3, -1., 0.]
  expected = flex.complex_double(f_obs.data().size(), 0)
  for w, r in zip(weights, regs.secondary):
    expected += sf.of_points(r.points)*w
  m = regs.mask_of_points(points = [r.points for r in regs.secondary],
    values = weights)
  assert approx_equal(sf.of_map(map_data = m), expected, 1.e-6)
  assert approx_equal(flex.sum(regs.mask_of(regs.secondary).as_1d()),
    sum([r.n_points for r in regs.secondary]))
  # padded maps (ready for in-place Fourier transform): same values, same
  # structure factors, no second copy of the map
  mp_ = regs.mask_of_points(points = [r.points for r in regs.secondary],
    values = weights, padded = True)
  assert mp_.is_padded() and mp_.focus() == regs.n_real
  assert not m.is_padded()
  from cctbx import maptbx
  mu = mp_.deep_copy()
  maptbx.unpad_in_place(map = mu)
  assert approx_equal(mu.as_1d(), m.as_1d())
  assert approx_equal(sf.of_map(map_data = mp_), expected, 1.e-6)
  sm = regs.solvent_mask()
  smp = regs.solvent_mask(padded = True)
  assert smp.is_padded() and not sm.is_padded()
  assert approx_equal(sf.of_map(map_data = smp), sf.of_map(map_data = sm),
    1.e-6)
  assert approx_equal(flex.sum(sm.as_1d()), n_solvent) # sm is not changed

# ------------------------------------------------------------------------------
# Data with known region occupancies
# ------------------------------------------------------------------------------

class synthetic(object):
  """
  Shell model (three cavities) with F_obs computed from

    F_calc + k_mask(s) * (F_whole + sum_n (q_n-1) F_n),

  k_mask = 0.35*exp(-50*s^2/4), with given occupancies q of the three
  cavities (largest first). extra: Cartesian sites of atoms that are in the
  crystal (F_obs) but not in the model.
  """

  def __init__(self, q, d_min=3.5, extra=None, free_fraction=0.05,
               space_group_symbol="P 21 21 21"):
    self.q = q
    self.mask_params = mmtbx.masks.mask_master_params.extract()
    self.xray_structure = shell_model(space_group_symbol=space_group_symbol)
    f_calc = self.xray_structure.structure_factors(d_min=d_min).f_calc()
    fm = mmtbx.f_model.manager(f_obs = abs(f_calc),
      xray_structure = self.xray_structure, mask_params = self.mask_params)
    self.regions = mosaic.regions(xray_structure = self.xray_structure,
      mask_params = self.mask_params)
    self.sf = mosaic.structure_factors(regions = self.regions,
      miller_array = f_calc, d_min = self.mask_params.Fmask_res_high)
    self.f_whole = fm.f_masks()[0].data()
    self.F = [self.sf.of_points(r.points) for r in self.regions.secondary]
    assert len(self.F) == len(q) == 3
    self.f_calc = f_calc
    f_atoms = f_calc.data()
    if(extra is not None):
      full = shell_model(space_group_symbol=space_group_symbol, extra=extra)
      f_atoms = f_calc.structure_factors_from_scatterers(
        xray_structure = full).f_calc().data()
    f_mask = self.f_whole.deep_copy()
    for qn, f in zip(q, self.F): f_mask += f*(qn-1.)
    self.k_mask = 0.35*flex.exp(-50.*fm.ss)
    self.f_obs = f_calc.array(data = flex.abs(f_atoms + self.k_mask*f_mask))
    self.r_free_flags = self.f_obs.generate_r_free_flags(
      fraction = free_fraction, max_free = None) if free_fraction > 0 else None

  def fmodel(self, f_obs=None, r_free_flags=None):
    if(f_obs is None): f_obs = self.f_obs
    if(r_free_flags is None): r_free_flags = self.r_free_flags
    fm = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = self.xray_structure,
      mask_params    = self.mask_params)
    fm.update_all_scales(params = self.scaling_params(),
      remove_outliers=False, apply_scale_k1_to_f_obs=False)
    return fm

  def scaling_params(self):
    # k_sol, b_sol model: the one the data were made with
    params = bss.master_params.extract()
    params.mode = "slow"
    return params

def occupancies(result):
  return [r.occupancy for r in result.regions]

def exercise_recover():
  """
  Known occupancies are recovered by both algorithms; the correction arrays
  are what they are documented to be.
  """
  for q in [(0., 0.5, 0.8), (0.8, 0., 0.3), (0.7, 0.7, 0.7)]:
    o = synthetic(q = q)
    fm = o.fmodel()
    r_start = fm.r_work()
    n0 = mosaic.number_of_runs
    res = mosaic.run(fmodel = fm)
    assert mosaic.number_of_runs == n0+1
    assert res.n_fitted == 3 and res.n_small == 0 and res.n_regions() == 3
    assert res.algorithm == "alg4"
    assert approx_equal(occupancies(res), q, 0.08), [occupancies(res), q]
    volumes = [r.volume for r in o.regions.secondary]
    assert approx_equal([r.volume for r in res.regions], volumes)
    # correction = sum (q_n - 1) F_n, with the fitted q
    c = flex.complex_double(o.f_whole.size(), 0)
    for qn, f in zip(occupancies(res), o.F): c += f*(qn-1.)
    assert approx_equal(res.correction, c, 1.e-6)
    # using it lowers R
    fm2 = fm.deep_copy()
    fm2.update(f_mask = fm.f_masks()[0].array(
      data = fm.f_masks()[0].data() + res.correction))
    fm2.update_all_scales(params = o.scaling_params(),
      remove_outliers=False, apply_scale_k1_to_f_obs=False)
    if(q != (1., 1., 1.)): assert fm2.r_work() < r_start - 0.005
    assert fm2.r_work() < 0.02, fm2.r_work()
    # fmodel is not changed by run
    assert approx_equal(fm.r_work(), r_start, 1.e-12)
    # unphased search
    params = mosaic.default_params()
    params.algorithm = "alg2"
    res2 = mosaic.run(fmodel = fm, params = params)
    assert res2.algorithm == "alg2"
    assert approx_equal(occupancies(res2), q, 0.08), [occupancies(res2), q]

def exercise_bounds():
  """
  0 <= q <= 1: a region cannot have more solvent than the main bulk solvent,
  nor less than nothing.
  """
  assert mosaic.OCCUPANCY_MAX == 1
  # denser than bulk solvent: occupancy is one
  o = synthetic(q = (0.8, 0.3, 2.))
  fm = o.fmodel()
  for algorithm in ["alg4", "alg2"]:
    params = mosaic.default_params()
    params.algorithm = algorithm
    res = mosaic.run(fmodel = fm, params = params)
    q = occupancies(res)
    assert q[2] == 1, q
    assert max(q) <= 1 and min(q) >= 0
    assert approx_equal(q[:2], o.q[:2], 0.15), q
    c = flex.complex_double(o.f_whole.size(), 0)
    for qn, f in zip(q, o.F): c += f*(qn-1.)
    assert approx_equal(res.correction, c, 1.e-6)
  # the fit itself, with and without the limit
  sel = mosaic.work_selection(fm)
  scale = (fm.k_isotropic()*fm.k_anisotropic()*fm.k_masks()[0]).select(sel)
  G = [f.select(sel)*scale for f in o.F]
  whole = o.f_whole.select(sel)*scale
  atoms = fm.f_model().data().select(sel) - whole
  rest = whole - G[0] - G[1] - G[2]
  f_obs = fm.f_obs().data().select(sel)
  for algorithm in ["alg4", "alg2"]:
    for refine_scales in [True, False]:
      kw = dict(f_obs = f_obs, atoms = atoms, rest = rest, G = G,
        algorithm = algorithm, refine_scales = refine_scales)
      q = mosaic.fit_occupancies(q_max = None, **kw)
      assert q[2] > 1.6, list(q)
      q = mosaic.fit_occupancies(**kw)
      assert q[2] == 1 and max(q) <= 1 and min(q) >= 0, list(q)
      q = mosaic.fit_occupancies(q_max = 1.5, **kw)
      assert q[2] == 1.5 and max(q) <= 1.5, list(q)
  # nothing in the output about a limit being hit
  sio = StringIO()
  res.show(log = sio)
  assert not "q > 1" in sio.getvalue()
  assert "0 - empty, 1 - same" in sio.getvalue()
  assert not hasattr(res, "map_correction")
  # less than empty: occupancy cannot be negative
  o = synthetic(q = (-0.5, 1., 0.5))
  fm = o.fmodel()
  for algorithm in ["alg4", "alg2"]:
    params = mosaic.default_params()
    params.algorithm = algorithm
    res = mosaic.run(fmodel = fm, params = params)
    q = occupancies(res)
    assert q[0] == 0, q
    assert min(q) >= 0

def exercise_free_independent():
  """
  Free reflections are not used: changing their F_obs does not change q.
  Without free flags all reflections are used.
  """
  o = synthetic(q = (0., 0.5, 1.))
  fm = o.fmodel()
  q1 = occupancies(mosaic.run(fmodel = fm))
  free = o.r_free_flags.data()
  assert 50 < free.count(True) < free.size()//10
  # same scales, different F_obs of the free reflections
  fm2 = fm.deep_copy()
  data = fm.f_obs().data().deep_copy()
  data.set_selected(free, data.select(free)*
    (flex.random_double(free.count(True))*3+0.5))
  fm2.update(f_obs = fm.f_obs().array(data = data))
  assert fm2.r_free() > fm.r_free() + 0.1
  assert approx_equal(fm2.r_work(), fm.r_work(), 1.e-9)
  for algorithm in ["alg4", "alg2"]:
    params = mosaic.default_params()
    params.algorithm = algorithm
    qa = occupancies(mosaic.run(fmodel = fm,  params = params))
    qb = occupancies(mosaic.run(fmodel = fm2, params = params))
    assert approx_equal(qa, qb, 1.e-6)
  assert mosaic.work_selection(fm).count(False) == free.count(True)
  # no free reflections
  o = synthetic(q = (0., 0.5, 1.), free_fraction = 0)
  fm = o.fmodel()
  assert mosaic.work_selection(fm).all_eq(True)
  assert approx_equal(occupancies(mosaic.run(fmodel = fm)), o.q, 0.08)

def sites_in_sphere(center, radius, spacing):
  result = []
  n = int(radius/spacing)+1
  for i in range(-n, n+1):
    for j in range(-n, n+1):
      for k in range(-n, n+1):
        d = (i*spacing, j*spacing, k*spacing)
        if(d[0]**2+d[1]**2+d[2]**2 <= radius**2):
          result.append((center[0]+d[0], center[1]+d[1], center[2]+d[2]))
  return result

def exercise_unmodelled():
  """
  Atoms in a cavity that are in the crystal but not in the model, and no
  solvent there: the region comes out as full of bulk solvent (occupancy
  one, the upper limit), the others are not disturbed.
  """
  center, radius = SHELLS[2]
  extra = sites_in_sphere(center = center, radius = 3.6, spacing = 1.6)
  assert len(extra) > 40
  o = synthetic(q = (1., 1., 0.), extra = extra)
  fm = o.fmodel()
  res = mosaic.run(fmodel = fm)
  q = occupancies(res)
  assert q[2] == 1, q
  assert approx_equal(q[:2], [1, 1], 0.1), q
  assert max(q) <= 1

def exercise_small():
  o = synthetic(q = (0.5, 0.5, 0.))
  fm = o.fmodel()
  volumes = [r.volume for r in o.regions.secondary]
  params = mosaic.default_params()
  assert params.volume_min == 50 and params.algorithm == "alg4"
  # the smallest cavity is below the limit: emptied, not fitted
  params.volume_min = (volumes[1]+volumes[2])/2
  res = mosaic.run(fmodel = fm, params = params)
  assert res.n_small == 1 and res.n_fitted == 2 and res.n_regions() == 3
  assert approx_equal(res.volume_small, volumes[2])
  assert approx_equal(occupancies(res), [0.5, 0.5], 0.05), occupancies(res)
  c = -1.*o.F[2]
  for qn, f in zip(occupancies(res), o.F[:2]): c += f*(qn-1.)
  assert approx_equal(res.correction, c, 1.e-6)
  # all regions small: nothing to fit, all emptied
  params.volume_min = volumes[0]*2
  res = mosaic.run(fmodel = fm, params = params)
  assert res.n_small == 3 and res.n_fitted == 0 and res.regions == []
  assert res.mean_occupancy() is None
  assert approx_equal(res.correction, -1.*(o.F[0]+o.F[1]+o.F[2]), 1.e-6)
  sio = StringIO()
  res.show(log = sio)
  assert "emptied: 3" in sio.getvalue()
  assert res.summary() == "mosaic: regions 3, small 3, fitted 0"
  # volume_min=0: all fitted
  params.volume_min = 0
  res = mosaic.run(fmodel = fm, params = params)
  assert res.n_small == 0 and res.n_fitted == 3
  # bad value
  params.volume_min = -1
  try: mosaic.run(fmodel = fm, params = params)
  except Sorry as e: assert "volume_min" in str(e)
  else: raise Exception_expected

def exercise_not_applicable():
  o = synthetic(q = (0., 0.5, 1.))
  fm = o.fmodel()
  assert mosaic.not_applicable(fm) is None
  n0 = mosaic.number_of_runs
  def check(fmodel, text):
    sio = StringIO()
    assert mosaic.run(fmodel = fmodel, log = sio) is None
    assert sio.getvalue().startswith("mosaic: skipped (") and \
      text in sio.getvalue(), sio.getvalue()
    assert mosaic.run(fmodel = fmodel) is None # no log: no crash
  # no cavities
  xrs0 = shell_model(shells=[((18., 9., 56.), 3.0)])
  f_obs = abs(xrs0.structure_factors(d_min=3.5).f_calc())
  f_obs = abs(mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs0,
    k_sol = 0.35, b_sol = 50.).f_model()) # with bulk solvent
  fm0 = mmtbx.f_model.manager(f_obs = f_obs, xray_structure = xrs0)
  fm0.update_all_scales(params = o.scaling_params(), remove_outliers=False)
  assert not fm0.k_masks()[0].all_eq(0)
  check(fm0, "no isolated regions")
  # data without bulk solvent: scaling finds none
  fm00 = mmtbx.f_model.manager(xray_structure = o.xray_structure,
    f_obs = abs(o.xray_structure.structure_factors(d_min=3.5).f_calc()))
  params = bss.master_params.extract()
  params.bulk_solvent = False
  fm00.update_all_scales(params = params, remove_outliers=False)
  check(fm00, "no bulk-solvent contribution")
  # the cases below are rejected before any region analysis
  n0 = mosaic.number_of_runs
  # no low-resolution data: no bulk-solvent contribution to fit
  sel = fm.f_obs().d_spacings().data() < 3.9
  hi = fm.f_obs().select(sel)
  mp = mmtbx.masks.mask_master_params.extract()
  mp.Fmask_res_high = 4.0
  fmh = mmtbx.f_model.manager(f_obs = hi, xray_structure = o.xray_structure,
    r_free_flags = o.r_free_flags.select(sel), mask_params = mp)
  check(fmh, "work reflections at d >= 4.0 A")
  # few low-resolution reflections
  mp = mmtbx.masks.mask_master_params.extract()
  mp.Fmask_res_high = 20.0 # about 25 reflections
  fml = mmtbx.f_model.manager(f_obs = fm.f_obs(),
    xray_structure = o.xray_structure, r_free_flags = o.r_free_flags,
    mask_params = mp)
  check(fml, "at least 50 needed")
  # neutron
  xrs_n = o.xray_structure.deep_copy_scatterers()
  xrs_n.switch_to_neutron_scattering_dictionary()
  fmn = mmtbx.f_model.manager(f_obs = fm.f_obs(), xray_structure = xrs_n,
    r_free_flags = o.r_free_flags)
  check(fmn, "neutron")
  assert mosaic.number_of_runs == n0

def exercise_too_few_reflections():
  """
  A minimal number of work reflections per fitted occupancy is required.
  """
  assert mosaic._REFLECTIONS_PER_REGION == 20
  o = synthetic(q = (0., 0.5, 0.8))
  d = o.f_obs.d_spacings().data()
  work = ~o.r_free_flags.data()
  limit = None
  for dd in [20, 18, 16, 14, 13, 12, 11, 10, 9, 8]:
    n = ((d >= dd) & work).count(True)
    if(50 <= n < 60): limit = dd
  if(limit is None): # no such resolution limit for this model: 60..120
    for dd in [20, 18, 16, 14, 13, 12, 11, 10, 9, 8]:
      n = ((d >= dd) & work).count(True)
      if(60 <= n < 120): limit = dd
    mosaic._REFLECTIONS_PER_REGION = 40
  assert limit is not None
  mp = mmtbx.masks.mask_master_params.extract()
  mp.Fmask_res_high = limit
  fm = mmtbx.f_model.manager(f_obs = o.f_obs, r_free_flags = o.r_free_flags,
    xray_structure = o.xray_structure, mask_params = mp)
  fm.update_all_scales(params = o.scaling_params(), remove_outliers=False,
    apply_scale_k1_to_f_obs=False)
  assert not fm.k_masks()[0].all_eq(0)
  assert mosaic.not_applicable(fm) is None
  sio = StringIO()
  assert mosaic.run(fmodel = fm, log = sio) is None
  assert "mosaic: skipped (" in sio.getvalue() and \
    "work reflections for 3 regions" in sio.getvalue(), sio.getvalue()
  # enough when fewer regions are fitted
  params = mosaic.default_params()
  params.volume_min = (o.regions.secondary[0].volume +
                       o.regions.secondary[1].volume)/2
  res = mosaic.run(fmodel = fm, params = params)
  assert res is not None and res.n_fitted == 1 and res.n_small == 2
  mosaic._REFLECTIONS_PER_REGION = 20

def exercise_many_regions():
  """
  Run time and memory grow with the number of fitted regions: only the
  largest ones are fitted, the others stay as in the flat model.
  """
  assert mosaic._MAX_FITTED_REGIONS == 300
  o = synthetic(q = (0., 0.5, 0.2))
  fm = o.fmodel()
  mosaic._MAX_FITTED_REGIONS = 2
  try:
    res = mosaic.run(fmodel = fm)
  finally:
    mosaic._MAX_FITTED_REGIONS = 300
  assert res.n_fitted == 2 and res.n_not_fitted == 1 and res.n_small == 0
  assert res.n_regions() == 3
  volumes = [r.volume for r in o.regions.secondary]
  assert approx_equal([r.volume for r in res.regions], volumes[:2])
  assert approx_equal(res.volume_not_fitted, volumes[2])
  # the third region is not touched: correction has the first two only
  c = flex.complex_double(o.f_whole.size(), 0)
  for qn, f in zip(occupancies(res), o.F[:2]): c += f*(qn-1.)
  assert approx_equal(res.correction, c, 1.e-6)
  sio = StringIO()
  res.show(log = sio)
  assert "the smaller ones are not fitted (occupancy 1): 1 (" in \
    sio.getvalue(), sio.getvalue()
  assert res.summary().startswith(
    "mosaic: regions 3, small 0, fitted 2, not fitted 1")
  # same limit expressed as memory: (regions) x (reflections of the fit)
  assert mosaic._MAX_FIT_ELEMENTS == 60000000
  n_work = (mosaic.work_selection(fm) & (fm.f_obs().d_spacings().data() >=
    fm.mask_params.Fmask_res_high)).count(True)
  mosaic._MAX_FIT_ELEMENTS = n_work*2 + n_work//2
  try:
    res2 = mosaic.run(fmodel = fm)
  finally:
    mosaic._MAX_FIT_ELEMENTS = 60000000
  assert res2.n_fitted == 2 and res2.n_not_fitted == 1
  assert approx_equal(res2.correction, res.correction, 1.e-6)
  # default: nothing left out
  res = mosaic.run(fmodel = fm)
  assert res.n_fitted == 3 and res.n_not_fitted == 0
  sio = StringIO()
  res.show(log = sio)
  assert not "not fitted" in sio.getvalue()
  assert not "not fitted" in res.summary()

def exercise_p1_and_anomalous():
  """
  P1, and anomalous data (Friedel mates are separate reflections).
  """
  o = synthetic(q = (0., 0.5, 1.), space_group_symbol = "P 1", d_min = 4.5)
  fm = o.fmodel()
  assert approx_equal(occupancies(mosaic.run(fmodel = fm)), o.q, 0.08)
  o = synthetic(q = (0., 0.5, 1.))
  f_obs = o.f_obs.generate_bijvoet_mates()
  flags = o.r_free_flags.generate_bijvoet_mates()
  assert f_obs.anomalous_flag()
  fm = o.fmodel(f_obs = f_obs, r_free_flags = flags)
  res = mosaic.run(fmodel = fm)
  assert res.correction.size() == f_obs.data().size()
  assert approx_equal(occupancies(res), o.q, 0.08), occupancies(res)

def exercise_annotate():
  pdb_str = """\
CRYST1   60.000   62.000   64.000  90.00  90.00  90.00 P 21 21 21
ATOM      1  O   HOH A   1      24.000   9.000  56.000  1.00 20.00           O
ATOM      2  CA  GLY B   7      32.000  23.000  34.000  1.00 20.00           C
ATOM      3  CA  ALA B  99      16.000  40.000  14.000  1.00 20.00           C
END
"""
  import iotbx.pdb
  import mmtbx.model
  model = mmtbx.model.manager(model_input = iotbx.pdb.input(
    source_info = None, lines = pdb_str.splitlines()))
  o = synthetic(q = (0., 0.5, 1.))
  res = mosaic.run(fmodel = o.fmodel())
  assert [r.near for r in res.regions] == [None]*3
  sio = StringIO()
  res.show(log = sio)
  assert not "nearest residues" in sio.getvalue()
  res.annotate(model = model)
  # the water is 6 A from the centre of the largest cavity
  assert res.regions[0].near.startswith("A HOH1 (6.0)"), res.regions[0].near
  assert res.regions[2].near.startswith("B GLY7 (6.0)"), res.regions[2].near
  assert len(res.regions[1].near.split(",")) == 3
  sio = StringIO()
  res.show(log = sio, prefix = "  |")
  out = sio.getvalue()
  assert "nearest residues" in out and "A HOH1 (6.0)" in out
  for l in out.splitlines(): assert l.startswith("  |")
  assert res.summary().startswith("mosaic: regions 3, small 0, fitted 3, <q> ")
  # the largest regions only
  sio = StringIO()
  res.show(log = sio, max_regions = 1)
  out = sio.getvalue()
  assert "A HOH1 (6.0)" in out and not "B GLY7 (6.0)" in out
  assert "... and 2 smaller regions, <q> " in out
  sio = StringIO()
  res.show(log = sio, max_regions = 10)
  assert not "smaller regions" in sio.getvalue()

# ------------------------------------------------------------------------------
# Algorithms of the paper kept as library functions
# ------------------------------------------------------------------------------

def unphased_inputs():
  """
  I_obs = |F_calc + sum k_n F_n|^2 with known k_n, F_n being the main region
  and the three cavities.
  """
  o = synthetic(q = (1., 1., 1.))
  f_main = o.f_whole - o.F[0] - o.F[1] - o.F[2]
  k = [0.35, 0.1, 0.25, 0.5]
  masks = [o.f_calc.array(data = d) for d in [f_main]+o.F]
  data = o.f_calc.data().deep_copy()
  for kn, m in zip(k, masks): data += m.data()*kn
  i_obs = o.f_calc.array(data = flex.norm(data))
  return group_args(i_obs = i_obs, fc = o.f_calc, masks = masks, k = k)

def exercise_algorithm_2():
  inp = unphased_inputs()
  F = [inp.fc]+inp.masks
  answer = [1.]+inp.k
  x = flex.double([1, .01, .01, .01, .01])
  # gradients against finite differences
  g_anal = mosaic.tg(i_obs = inp.i_obs, F = F, x = x,
    use_curvatures = False).gradients()
  e = 1.e-6
  g_fd = flex.double()
  for i in range(x.size()):
    xp, xm = x.deep_copy(), x.deep_copy()
    xp[i] += e
    xm[i] -= e
    tp = mosaic.tg(i_obs = inp.i_obs, F = F, x = xp,
      use_curvatures = False).target()
    tm = mosaic.tg(i_obs = inp.i_obs, F = F, x = xm,
      use_curvatures = False).target()
    g_fd.append((tp-tm)/(2*e))
  assert flex.sum(flex.abs(g_anal-g_fd))*2/flex.sum(flex.abs(g_anal+g_fd)) \
    < 1.e-6
  # target is zero at the answer
  assert mosaic.tg(i_obs = inp.i_obs, F = F, x = flex.double(answer),
    use_curvatures = False).target() < 1.e-12
  # called as the old driver did (no use_lbfgsb): used to raise TypeError
  r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
    use_curvatures = False)
  assert r.size() == 5
  # the answer is found
  r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
    use_curvatures = True, use_lbfgsb = False)
  assert approx_equal(r, answer, 1.e-4)
  for use_lbfgsb in [False, True]:
    r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
      use_curvatures = False, use_lbfgsb = use_lbfgsb, macro_cycles = 100)
    assert approx_equal(r, answer, 1.e-2), list(r)
    r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
      use_curvatures = None, use_lbfgsb = use_lbfgsb)
    assert approx_equal(r, answer, 1.e-2), list(r)
  # Bounds. Without phases |F_calc + k*F_mask| has a second, false minimum
  # at large k (mask stronger than atoms); with F_calc fixed by the bounds
  # the search has to start on the right side of it.
  x = flex.double([1, .3, .3, .3, .3])
  # bounds not active: same answer
  r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
    lower = flex.double([1, 0, 0, 0, 0]),
    upper = flex.double([1]+[mosaic.NO_LIMIT]*4), macro_cycles = 100)
  assert approx_equal(r, answer, 1.e-2), list(r)
  # bounds active
  lower = flex.double([1, 0, 0.15, 0, 0])
  upper = flex.double([1, 0.3, 1, 1, 0.4])
  r = mosaic.algorithm_2(i_obs = inp.i_obs, F = F, x = x.deep_copy(),
    lower = lower, upper = upper, macro_cycles = 100)
  for ri, lo, hi in zip(r, lower, upper): assert lo <= ri <= hi
  assert r[0] == 1
  # at least one bound is active, and the result differs from the answer
  n_active = len([1 for ri, lo, hi in zip(r[1:], lower[1:], upper[1:])
    if abs(ri-lo) < 1.e-6 or abs(ri-hi) < 1.e-6])
  assert n_active >= 1, list(r)
  assert approx_equal(r[1], 0.3, 1.e-6), list(r)
  t_bounded = mosaic.tg(i_obs = inp.i_obs, F = F, x = r,
    use_curvatures = False).target()
  clipped = flex.double([min(hi, max(lo, v))
    for v, lo, hi in zip(answer, lower, upper)])
  t_clipped = mosaic.tg(i_obs = inp.i_obs, F = F, x = clipped,
    use_curvatures = False).target()
  assert t_bounded <= t_clipped*(1+1.e-6)

def exercise_algorithm_2_macro_cycles():
  """
  With bounds: macro-cycles go on until the result stops changing. A few
  iterations per macro-cycle then reach the same answer as many.
  """
  inp = unphased_inputs()
  F = [inp.fc]+inp.masks
  answer = [1.]+inp.k
  x = flex.double([1, .3, .3, .3, .3])
  kw = dict(i_obs = inp.i_obs, F = F, lower = flex.double([1, 0, 0, 0, 0]),
    upper = flex.double([1]+[mosaic.NO_LIMIT]*4))
  r_one = mosaic.algorithm_2(x = x.deep_copy(), max_iterations = 3,
    macro_cycles = 1, **kw)
  assert not approx_equal(r_one, answer, 1.e-2, out = None) # too short
  r_many = mosaic.algorithm_2(x = x.deep_copy(), max_iterations = 3,
    macro_cycles = 200, **kw)
  assert approx_equal(r_many, answer, 1.e-2), list(r_many)
  # the start is not changed
  x0 = x.deep_copy()
  mosaic.algorithm_2(x = x, macro_cycles = 5, **kw)
  assert approx_equal(x, x0, 1.e-12)

def exercise_algorithm_3():
  inp = unphased_inputs()
  r = mosaic.algorithm_3(i_obs = inp.i_obs, fc = inp.fc, f_masks = inp.masks)
  assert approx_equal(r, [1.]+inp.k, 1.e-6)

def exercise_alg4():
  """
  Python interface of the phased search: scales of F_calc and of all masks
  found from amplitudes, starting from F_calc alone (as in the paper).
  """
  inp = unphased_inputs()
  f_obs = flex.sqrt(inp.i_obs.data())
  zero = flex.complex_double(f_obs.size(), 0)
  x_start = flex.double([1, 0, 0, 0, 0])
  r = mosaic.alg4(f_obs = f_obs, base = zero,
    G = [inp.fc.data()]+[m.data() for m in inp.masks],
    x_start = x_start, n_cycles = 200)
  assert approx_equal(r, [1.]+inp.k, 1.e-3), list(r)
  # F_calc fixed, masks fitted, with bounds
  r = mosaic.alg4(f_obs = f_obs, base = inp.fc.data(),
    G = [m.data() for m in inp.masks], n_cycles = 100,
    lower = flex.double(4, 0), upper = flex.double([1, 1, 0.2, 1]))
  assert r[2] == 0.2
  assert flex.min(r) >= 0

def run():
  random.seed(0)
  flex.set_random_seed(0)
  exercise_normal_matrix()
  exercise_solve_bounded()
  exercise_alg4_kernel()
  exercise_alg2_tg_update()
  exercise_regions()
  exercise_padded_maps()
  exercise_region_structure_factors()
  exercise_alg4()
  exercise_algorithm_2()
  exercise_algorithm_2_macro_cycles()
  exercise_algorithm_3()
  exercise_recover()
  exercise_bounds()
  exercise_free_independent()
  exercise_unmodelled()
  exercise_small()
  exercise_not_applicable()
  exercise_p1_and_anomalous()
  exercise_many_regions()
  exercise_too_few_reflections()
  exercise_annotate()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %.2f"%(time.time()-t0))
  print("OK")
