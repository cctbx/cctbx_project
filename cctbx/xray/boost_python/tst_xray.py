from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import maptbx
from cctbx import eltbx
from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import xray
from cctbx import math_module
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import pickle

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
  x = xray.scatterer("a", (0.1,0.2,0.3), 0.25, 0.9, "const", 0, 0)
  assert x.label == "a"
  x.label = "b"
  assert x.label == "b"
  assert x.scattering_type == "const"
  x.scattering_type = "Si"
  assert x.scattering_type == "Si"
  assert x.fp == 0
  assert x.fdp == 0
  x.fp = 1
  assert x.fp == 1
  x.fdp = 2
  assert x.fdp == 2
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
  x = xray.scatterer("a", (0.1,0.2,0.3), (1,2,3,4,5,6), 0.9, "const", 0, 0)
  assert x.anisotropic_flag
  assert approx_equal(x.u_star, (1,2,3,4,5,6))
  x.u_star = (3,2,1,6,5,4)
  assert approx_equal(x.u_star, (3,2,1,6,5,4))
  x.anisotropic_flag = 0
  assert not x.anisotropic_flag
  x = xray.scatterer(
    "si1", site=(0.01,0.02,0.3), occupancy=0.9, u=(0.3, 0.3, 0.2, 0,0,0))
  assert x.scattering_type == "Si"
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

def exercise_scattering_dictionary():
  sd = xray.scattering_dictionary()
  assert sd.n_scatterers() == 0
  assert sd.dict_size() == 0
  assert len(sd.dict()) == 0
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1"),
    xray.scatterer("Si2"),
    xray.scatterer("O1"),
    xray.scatterer("O2"),
    xray.scatterer("Al1"),
    xray.scatterer("O3"),
    xray.scatterer("Al2"),
    xray.scatterer("const", scattering_type="const"),
    xray.scatterer("custom", scattering_type="custom")))
  sd = xray.scattering_dictionary(scatterers)
  assert sd.n_scatterers() == 9
  assert sd.dict_size() == 5
  sd_dict = sd.dict()
  assert len(sd_dict) == 5
  all_keys = sd_dict.keys()
  all_keys.sort()
  assert all_keys == ["Al", "O", "Si", "const", "custom"]
  for k,v in sd_dict.items():
    if   (k == "Si"): assert tuple(v.member_indices) == (0,1)
    elif (k == "O"): assert tuple(v.member_indices) == (2,3,5)
    elif (k == "Al"): assert tuple(v.member_indices) == (4,6)
    elif (k == "const"): assert tuple(v.member_indices) == (7,)
    elif (k == "custom"): assert tuple(v.member_indices) == (8,)
    assert v.gaussian.n_terms() == 0
    assert v.gaussian.c() == 0
    assert tuple(sd.lookup(k).member_indices) == tuple(v.member_indices)
  z = list(sd.find_undefined())
  z.sort()
  assert z == all_keys
  p = list(sd.scatterer_permutation())
  p.sort()
  assert p == range(9)
  for table,n_terms in (("IT1992",4), ("WK1995",5)):
    sd = xray.scattering_dictionary(scatterers)
    sd.assign("const", eltbx.xray_scattering.gaussian(10))
    sd.assign("custom", eltbx.xray_scattering.gaussian((1,2),(3,4),5))
    sd.assign_from_table(table)
    for k,v in sd.dict().items():
      if (k in ("Al", "O", "Si")):
        assert v.gaussian.n_terms() == n_terms
      elif (k == "const"):
        assert v.gaussian.n_terms() == 0
        assert approx_equal(v.gaussian.c(), 10)
      else:
        assert v.gaussian.n_terms() == 2
        assert approx_equal(v.gaussian.c(), 5)
    sd.assign("Al", eltbx.xray_scattering.gaussian(20))
    assert approx_equal(sd.lookup("Al").gaussian.c(), 20)
  assert sd.find_undefined().size() == 0
  g = sd.dict()["custom"]
  c = g.gaussian
  assert c.n_terms() == 2
  assert approx_equal(c.array_of_a(), (1,2))
  assert approx_equal(c.array_of_b(), (3,4))
  assert approx_equal(c.c(), 5)
  assert tuple(g.member_indices) == (8,)
  s = pickle.dumps(g)
  l = pickle.loads(s)
  c = l.gaussian
  assert c.n_terms() == 2
  assert approx_equal(c.array_of_a(), (1,2))
  assert approx_equal(c.array_of_b(), (3,4))
  assert approx_equal(c.c(), 5)
  assert tuple(l.member_indices) == (8,)
  s = pickle.dumps(sd)
  l = pickle.loads(s)
  l_dict = l.dict()
  for k,v in sd.dict().items():
    w = l_dict[k]
    assert tuple(v.member_indices) == tuple(w.member_indices)
    vc = v.gaussian
    wc = w.gaussian
    assert vc.array_of_a() == wc.array_of_a()
    assert vc.array_of_b() == wc.array_of_b()
    assert vc.c() == wc.c()
  try:
    sd.lookup("undef")
  except RuntimeError, e:
    assert str(e).startswith(
      "cctbx Error: Label not in scattering dictionary: ")
  else:
    raise RuntimeError("Exception expected.")

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
  scattering_dict = xray.ext.scattering_dictionary(scatterers)
  scattering_dict.assign_from_table("WK1995")
  for sf in (xray.ext.structure_factors_simple,
             xray.ext.structure_factors_direct):
    fc = sf(uc, sg.group(), mi, scatterers, scattering_dict).f_calc()
    a = flex.abs(fc)
    p = flex.arg(fc, 1)
    assert approx_equal(tuple(a), (10.50871, 9.049631))
    assert approx_equal(tuple(p), (-36, 72))
  xray.ext.structure_factors_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, scattering_dict).f_calc()
  xray.ext.structure_factors_gradients_direct(
    uc, sg.group(), mi, scatterers, scattering_dict,
    flex.complex_double(mi.size()),
    xray.ext.gradient_flags(00001, 00001, 00001, 00001, 00001, 00001),
    0)
  xray.ext.structure_factors_gradients_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, scattering_dict,
    flex.complex_double(mi.size()),
    xray.ext.gradient_flags(00001, 00001, 00001, 00001, 00001, 00001),
    0)

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
    (0.0013784963-0.002756992j), (0.0103982354-0.013864313j),
    (0.0160141831-0.032028366j), (0.0004572786+0.000365822j),
    (0.0014117387+0.001694086j)))
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
    (0.002855645-0.005711291j), (0.035410006-0.047213342j),
    (0.010851453-0.021702907j), (0.005711291+0.004569033j),
    (0.024748929+0.029698715j)))

def exercise_sampled_model_density():
  assert approx_equal(xray.calc_u_base(2, 1./3), 0.1350949)
  uc = uctbx.unit_cell((20, 20, 23))
  sg = sgtbx.space_group_info("P 4")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  for scatterer in scatterers:
    scatterer.apply_symmetry(uc, sg.group())
  scattering_dict = xray.ext.scattering_dictionary(scatterers)
  scattering_dict.assign_from_table("WK1995")
  d = xray.sampled_model_density(
    unit_cell=uc,
    scatterers=scatterers,
    scattering_dict=scattering_dict,
    fft_n_real=(20,20,22),
    fft_m_real=(20,20,23))
  assert d.unit_cell().is_similar_to(uc)
  assert approx_equal(d.u_base(), 0.25)
  assert approx_equal(d.u_extra(), 0.25)
  assert approx_equal(d.u_min(), 0)
  assert approx_equal(d.ave_u_iso_or_equiv(), 0.025)
  assert approx_equal(d.wing_cutoff(), 1.e-3)
  assert approx_equal(d.exp_table_one_over_step_size(), -100)
  assert approx_equal(d.tolerance_positive_definite(), 1.e-5)
  assert d.n_scatterers_passed() == 2
  assert d.n_contributing_scatterers() == 2
  assert d.n_anomalous_scatterers() == 1
  assert d.anomalous_flag()
  assert d.real_map().size() == 0
  assert d.complex_map().size() == (20*20*22)
  assert d.exp_table_size() == 2834
  assert d.max_sampling_box_n_points() == 180
  assert d.sum_sampling_box_n_points() == 305
  assert approx_equal(d.ave_sampling_box_n_points(), 305/2.)
  assert d.max_sampling_box_edges() == (5,6,6)
  assert approx_equal(d.max_sampling_box_edges_frac(), (5/20.,6/20.,6/22.))
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

def exercise_minimization_apply_shifts():
  uc = uctbx.unit_cell((20, 20, 23))
  sg = sgtbx.space_group_info("P 4")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  scattering_dict = xray.ext.scattering_dictionary(scatterers)
  scattering_dict.assign_from_table("WK1995")
  f = xray.ext.gradient_flags(0001, 0001, 0001, 0001, 0001, 0001)
  shifts = flex.double(19, 0.001)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  for a,b in zip(scatterers, shifted_scatterers):
    assert a.scattering_type == b.scattering_type
    assert a.site != b.site
    if (not a.anisotropic_flag):
      assert a.u_iso != b.u_iso
      assert approx_equal(a.u_star, b.u_star)
    else:
      assert a.u_iso == b.u_iso
      assert not approx_equal(a.u_star, b.u_star)
    assert a.occupancy != b.occupancy
    assert a.fp != b.fp
    assert a.fdp != b.fdp
  shifts = flex.double(19, -0.001)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), shifted_scatterers, scattering_dict, f, shifts, 0)
  for a,b in zip(scatterers, shifted_scatterers):
    assert a.scattering_type == b.scattering_type
    assert approx_equal(a.site, b.site)
    assert approx_equal(a.u_iso, b.u_iso)
    assert approx_equal(a.u_star, b.u_star)
    assert approx_equal(a.occupancy, b.occupancy)
    assert approx_equal(a.fp, b.fp)
    assert approx_equal(a.fdp, b.fdp)
  f = xray.ext.gradient_flags(0001, 00000, 00000, 00000, 00000, 00000)
  shifts = flex.double((-1,2,-3,4,-5,-6))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  assert approx_equal(
    shifted_scatterers[0].site,
    (0.01-1/20.,0.02+2/20.,0.3-3/23.))
  assert approx_equal(
    shifted_scatterers[1].site,
    (0.3+4/20.,0.4-5/20.,0.5-6/23.))
  f = xray.ext.gradient_flags(00000, 0001, 00000, 00000, 00000, 00000)
  shifts = flex.double(1, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  assert shifted_scatterers[0].u_iso == 0
  f = xray.ext.gradient_flags(00000, 00000, 0001, 00000, 00000, 00000)
  shifts = flex.double(6, -100)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  assert not approx_equal(scatterers[1].u_star, shifted_scatterers[1].u_star)
  shifts = flex.double(6, -1000)
  shifted_scatterers_2 = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  # due to eigenvalue filtering two extreme shifts
  # should produce the same result
  assert approx_equal(shifted_scatterers[1].u_star,
                      shifted_scatterers_2[1].u_star)
  f = xray.ext.gradient_flags(00000, 00000, 00000, 0001, 00000, 00000)
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  for i in xrange(2):
    assert shifted_scatterers[i].occupancy == 0
  f = xray.ext.gradient_flags(00000, 00000, 00000, 00000, 0001, 00000)
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  assert approx_equal(shifted_scatterers[0].fp, -11)
  assert shifted_scatterers[1].fp == -10
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 3)
  for i in xrange(2):
    assert approx_equal(
      shifted_scatterers[i].fp,
      -(scattering_dict
          .lookup(shifted_scatterers[i].scattering_type)
          .gaussian.at_d_star_sq(1/9.)))
    assert shifted_scatterers[i].fdp == scatterers[i].fdp
  f = xray.ext.gradient_flags(00000, 00000, 00000, 00000, 00000, 0001)
  shifts = flex.double((2,3))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, sg.type(), scatterers, scattering_dict, f, shifts, 3)
  assert shifted_scatterers[0].fp == -1
  assert approx_equal(shifted_scatterers[0].fdp, 4)
  assert shifted_scatterers[1].fp == 0
  assert shifted_scatterers[1].fdp == 3
  shifts = flex.double(1)
  try:
    xray.ext.minimization_apply_shifts(
      uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  except Exception, e:
    assert str(e) == "scitbx Error: Array of shifts is too small."
  else:
    raise RuntimeError("Exception expected.")
  shifts = flex.double(3)
  try:
    xray.ext.minimization_apply_shifts(
      uc, sg.type(), scatterers, scattering_dict, f, shifts, 0)
  except Exception, e:
    assert str(e) == "cctbx Error: Array of shifts is too large."
  else:
    raise RuntimeError("Exception expected.")

def exercise_asu_mappings():
  from cctbx.development import random_structure
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["O"]*10)
  asu_mappings = crystal.direct_space_asu.asu_mappings(
    space_group=structure.space_group(),
    asu=structure.direct_space_asu().as_float_asu(),
    buffer_thickness=3)
  xray.asu_mappings_process(
    asu_mappings=asu_mappings,
    scatterers=structure.scatterers())
  assert asu_mappings.mappings().size() == structure.scatterers().size()

def run():
  exercise_conversions()
  exercise_gradient_flags()
  exercise_xray_scatterer()
  exercise_scattering_dictionary()
  exercise_rotate()
  exercise_structure_factors()
  exercise_targets()
  exercise_sampled_model_density()
  exercise_minimization_apply_shifts()
  exercise_asu_mappings()
  print "OK"

if (__name__ == "__main__"):
  run()
