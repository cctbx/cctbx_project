from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import maptbx
from cctbx.eltbx.caasf import wk1995
from cctbx import xray
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_conversions():
  d = flex.double((10,-1))
  s = flex.double((1,2))
  r = xray.array_f_sq_as_f(d, s)
  r = xray.array_f_sq_as_f(d, s, 1.e-6)
  assert approx_equal(tuple(r.f), (3.1622777, 0))
  assert approx_equal(tuple(r.sigma_f), (0.1543471, 1.4142136))
  r = xray.array_f_sq_as_f(d)
  assert approx_equal(tuple(r.f), (3.1622777, 0))
  assert r.sigma_f.size() == 0
  r = xray.array_f_as_f_sq(d, s)
  assert approx_equal(tuple(r.f_sq), (100, 1))
  assert approx_equal(tuple(r.sigma_f_sq), (20, -4))
  r = xray.array_f_as_f_sq(d)
  assert approx_equal(tuple(r.f_sq), (100, 1))
  assert r.sigma_f_sq.size() == 0

def exercise_gradient_flags():
  f = xray.ext.gradient_flags(00000, 0001, 00000, 0001, 00000, 0001)
  assert not f.site
  assert f.u_iso
  assert not f.u_aniso
  assert f.occupancy
  assert not f.fp
  assert f.fdp
  f.site = 0001
  f.u_iso = 00000
  f.u_aniso = 0001
  f.occupancy = 00000
  f.fp = 0001
  f.fdp = 00000
  assert f.site
  assert not f.u_iso
  assert f.u_aniso
  assert not f.occupancy
  assert f.fp
  assert not f.fdp
  c = xray.ext.gradient_flags(f)
  assert c.site
  assert not c.u_iso
  assert c.u_aniso
  assert not c.occupancy
  assert c.fp
  assert not c.fdp
  assert not f.all_false()
  assert xray.ext.gradient_flags(
    00000, 00000, 00000, 00000, 00000, 00000).all_false()
  f.u_iso = 0001
  assert f.adjust(00000).u_iso == 0001
  assert f.adjust(0001).u_iso == 00000
  assert f.adjust(00000).u_aniso == 00000
  assert f.adjust(0001).u_aniso == 0001

def exercise_xray_scatterer():
  caasf = wk1995("const")
  x = xray.scatterer("a", (0.1,0.2,0.3), 0.25, 0.9, caasf, 0j)
  assert x.label == "a"
  x.label = "b"
  assert x.label == "b"
  assert x.caasf.label() == "const"
  x.caasf = wk1995("SI")
  assert x.caasf.label() == "Si"
  assert x.fp_fdp == 0j
  x.fp_fdp = 1+2j
  assert x.fp_fdp == 1+2j
  assert approx_equal(x.site, (0.1,0.2,0.3))
  x.site = (0.3,-0.4,0.5)
  assert approx_equal(x.site, (0.3,-0.4,0.5))
  assert approx_equal(x.occupancy, 0.9)
  x.occupancy = 0.3
  assert approx_equal(x.occupancy, 0.3)
  assert not x.anisotropic_flag
  assert approx_equal(x.u_iso, 0.25)
  x.u_iso = 0.52
  assert approx_equal(x.u_iso, 0.52)
  x = xray.scatterer("a", (0.1,0.2,0.3), (1,2,3,4,5,6), 0.9, caasf, 0j)
  assert x.anisotropic_flag
  assert approx_equal(x.u_star, (1,2,3,4,5,6))
  x.u_star = (3,2,1,6,5,4)
  assert approx_equal(x.u_star, (3,2,1,6,5,4))
  x.anisotropic_flag = 0
  assert not x.anisotropic_flag
  x = xray.scatterer(
    "si1", site=(0.01,0.02,0.3), occupancy=0.9, u=(0.3, 0.3, 0.2, 0,0,0))
  assert x.caasf.label() == "Si"
  uc = uctbx.unit_cell((10, 10, 13))
  sg = sgtbx.space_group_info("P 4")
  ss = x.apply_symmetry(uc, sg.group())
  assert x.multiplicity() == 1
  assert approx_equal(x.weight_without_occupancy(), 1/4.)
  assert approx_equal(x.weight(), 0.9/4.)
  assert approx_equal(x.site, (0,0,0.3))
  assert ss.multiplicity() == x.multiplicity()
  x.occupancy = 0.8
  x.update_weight(sg.group().order_z())
  assert approx_equal(x.weight_without_occupancy(), 1/4.)
  assert approx_equal(x.weight(), 0.8/4.)
  u_cart = (0.3354, 0.3771, 0.4874, -0.05161, 0.026763, -0.02116)
  x.u_star = adptbx.u_cart_as_u_star(uc, u_cart)
  x.anisotropic_flag = 1
  try: x.apply_symmetry(uc, sg.group(), u_star_tolerance=0.1)
  except: pass
  else: raise AssertionError, "Exception expected."
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0, 0)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0, 0, 0)
  assert ss.is_compatible_u_star(x.u_star)
  assert approx_equal(x.u_star, (0.0035625, 0.0035625, 0.002884, 0, 0, 0))

def exercise_rotate():
  uc = uctbx.unit_cell((10, 10, 13))
  s = flex.xray_scatterer((xray.scatterer("Si1", site=(0.01,0.02,0.3)),))
  r = xray.rotate(uc, ((1,0,0, 0,1,0, 0,0,1)), s)
  assert r.size() == 1
  assert approx_equal(s[0].site, r[0].site)
  r = xray.rotate(uc, ((0,-1,0, -1,0,0, 0,0,-1)), s)
  assert approx_equal(r[0].site, (-0.02,-0.01,-0.3))

def exercise_pack_parameters():
  unit_cell = uctbx.unit_cell((10, 10, 10))
  sg = sgtbx.space_group_info(symbol="P 1")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), u=0.2),
    xray.scatterer("O1", site=(0.3,0.4,0.5), u=(0.4,0.5,0.6,-.05,0.2,-0.02))))
  for scatterer in scatterers:
    scatterer.apply_symmetry(unit_cell, sg.group())
  for uc in (unit_cell, None):
    x = flex.double()
    assert xray.pack_parameters(uc, scatterers, x, 00000, 00000, 00000) == 0
  x = flex.double()
  assert xray.pack_parameters(unit_cell, scatterers, x, 0001, 00000, 00000)==6
  assert approx_equal(tuple(x), (0.1, 0.2, 3, 3, 4, 5))
  x = flex.double()
  assert xray.pack_parameters(None, scatterers, x, 0001, 00000, 00000)==6
  assert approx_equal(tuple(x), (0.01, 0.02, 0.3, 0.3, 0.4, 0.5))
  x = flex.double()
  assert xray.pack_parameters(None, scatterers, x, 00000, 0001, 00000)==1
  assert approx_equal(tuple(x), (0.2,))
  x = flex.double()
  assert xray.pack_parameters(None, scatterers, x, 00000, 00000, 0001)==2
  assert approx_equal(tuple(x), (1,1))
  for start in (0, 10):
    x = flex.double(start)
    assert xray.pack_parameters(None, scatterers, x, 0001, 0001, 0001)==9+start
    x *= 2
    sc = scatterers.deep_copy()
    assert xray.unpack_parameters(
      None, 0, x, start, sc, 00000, 00000, 00000)==0+start
    assert xray.unpack_parameters(
      None, sg.group().order_z(), x, start, sc, 0001, 0001, 0001)==9+start
    assert approx_equal(sc[0].site, (0.02, 0.04, 0.6))
    assert approx_equal(sc[0].u_iso, 0.4)
    assert approx_equal(sc[0].occupancy, 2)
    assert approx_equal(sc[0].weight_without_occupancy(), 1)
    assert approx_equal(sc[0].weight(), 2)
    assert approx_equal(sc[1].site, (0.6, 0.8, 1))
    assert approx_equal(sc[1].occupancy, 2)

def exercise_structure_factors():
  uc = uctbx.unit_cell((10, 10, 13))
  sg = sgtbx.space_group_info("P 4")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3)),
    xray.scatterer("O1", site=(0.3,0.4,0.5), u=(0.4,0.5,0.6,-.05,0.2,-0.02))))
  for s in scatterers:
    assert s.multiplicity() == 0
  assert tuple(xray.apply_symmetry(uc, sg.group(), scatterers)) == (0,)
  for s in scatterers:
    assert s.multiplicity() != 0
  mi = flex.miller_index(((1,2,3), (2,3,4)))
  fc = xray.structure_factors_direct_with_first_derivatives(
    uc, sg.group(), mi, scatterers,
    flex.complex_double(),
    xray.ext.gradient_flags(00000, 00000, 00000, 00000, 00000, 00000)).f_calc()
  a = flex.abs(fc)
  p = flex.arg(fc, 1)
  assert approx_equal(tuple(a), (10.50871, 9.049631))
  assert approx_equal(tuple(p), (-36, 72))

def exercise_targets():
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, 0001)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, 00000, 3)
  assert approx_equal(ls.scale_factor(), 3)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual(f_obs, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double((10,20,30,40,50))
  ls = xray.targets_least_squares_residual(f_obs, f_calc, 0001)
  assert approx_equal(ls.scale_factor(), 1/10.)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual(f_obs, f_calc, 00000, 3/10.)
  assert approx_equal(ls.scale_factor(), 3/10.)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j,5-4j,-5+6j))
  w = flex.double((1,2,3,2,4))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, 0001)
  assert approx_equal(ls.scale_factor(), 0.6307845)
  assert approx_equal(ls.target(), 0.06211855)
  assert approx_equal(tuple(ls.derivatives()), (
    (0.0013784963+0.002756992j), (0.0103982354+0.013864313j),
    (0.0160141831+0.032028366j), (0.0004572786-0.000365822j),
    (0.0014117387-0.001694086j)))
  f_obs = flex.double((1,2,3,4,5))
  w = flex.int((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  ic = xray.targets_intensity_correlation(f_obs, w, f_calc)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert ic.derivatives().size() == 0
  ic = xray.targets_intensity_correlation(f_obs, w, f_calc, 0001)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert approx_equal(tuple(ic.derivatives()), (0j,0j,0j,0j,0j))
  ic = xray.targets_intensity_correlation(f_obs, f_calc)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert ic.derivatives().size() == 0
  f_calc = flex.complex_double((10,20,30,40,50))
  ic = xray.targets_intensity_correlation(f_obs, f_calc, 0001)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert approx_equal(tuple(ic.derivatives()), (0j,0j,0j,0j,0j))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j,5-4j,-5+6j))
  w = flex.int((1,2,3,2,4))
  ic = xray.targets_intensity_correlation(f_obs, w, f_calc, 0001)
  assert approx_equal(ic.correlation(), 0.8932460)
  assert approx_equal(ic.target(), 1-ic.correlation())
  assert approx_equal(tuple(ic.derivatives()), (
    (0.002855645+0.005711291j), (0.035410006+0.047213342j),
    (0.010851453+0.021702907j), (0.005711291-0.004569033j),
    (0.024748929-0.029698715j)))

def exercise_sampled_model_density():
  assert approx_equal(xray.calc_u_extra(2, 1./3), 0.1350949)
  uc = uctbx.unit_cell((20, 20, 23))
  sg = sgtbx.space_group_info("P 4")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp_fdp=(-1+2j)),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  for scatterer in scatterers:
    scatterer.apply_symmetry(uc, sg.group())
  d = xray.sampled_model_density(uc, scatterers, (20,20,22), (20,20,23))
  assert d.unit_cell().is_similar_to(uc)
  assert approx_equal(d.u_extra(), 0.25)
  assert approx_equal(d.wing_cutoff(), 1.e-3)
  assert approx_equal(d.exp_table_one_over_step_size(), -100)
  assert d.n_scatterers_passed() == 2
  assert d.n_contributing_scatterers() == 2
  assert d.n_anomalous_scatterers() == 1
  assert d.anomalous_flag()
  assert d.real_map().size() == 0
  assert d.complex_map().size() == (20*20*22)
  assert d.exp_table_size() == 1968
  assert d.max_shell_radii() == (2,2,2)
  assert approx_equal(d.max_shell_radii_frac(), (1/10.,1/10.,1/11.))
  t = maptbx.grid_tags((20,20,22))
  f = maptbx.symmetry_flags(0001)
  t.build(sg.type(), f)
  d.apply_symmetry(t)
  i = flex.miller_index(((1,2,3), (2,3,4)))
  f = flex.complex_double((1+2j, 2+3j))
  d.eliminate_u_extra_and_normalize(i, f)
  f_orig = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), 0.2, i, f)
  f_elim = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), -0.2, i, f, 1)
  assert approx_equal(f, f_orig)
  m = flex.double((1,2))
  xray.apply_u_extra(d.unit_cell(), 0.2, i, f, m)
  assert approx_equal(f[0], f_elim[0])
  assert approx_equal(f[1], f_elim[1]*2)
  m = flex.double((1,1/2.))
  xray.apply_u_extra(d.unit_cell(), -0.2, i, f, m)
  assert approx_equal(f, f_orig)

def run():
  exercise_conversions()
  exercise_gradient_flags()
  exercise_xray_scatterer()
  exercise_rotate()
  exercise_pack_parameters()
  exercise_structure_factors()
  exercise_targets()
  exercise_sampled_model_density()
  print "OK"

if (__name__ == "__main__"):
  run()
