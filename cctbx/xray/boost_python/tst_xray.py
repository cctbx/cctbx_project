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
from libtbx.test_utils import approx_equal, not_approx_equal, show_diff
from cStringIO import StringIO
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
  f = xray.ext.gradient_flags(
    False, True, False, True, False, True, False, True)
  assert not f.site
  assert f.u_iso
  assert not f.u_aniso
  assert f.occupancy
  assert not f.fp
  assert f.fdp
  assert not f.sqrt_u_iso
  f.site = True
  f.u_iso = False
  f.u_aniso = True
  f.occupancy = False
  f.fp = True
  f.fdp = False
  f.sqrt_u_iso = True
  assert f.site
  assert not f.u_iso
  assert f.u_aniso
  assert not f.occupancy
  assert f.fp
  assert not f.fdp
  assert f.sqrt_u_iso
  c = xray.ext.gradient_flags(f)
  assert c.site
  assert not c.u_iso
  assert c.u_aniso
  assert not c.occupancy
  assert c.fp
  assert not c.fdp
  assert not f.all_false()
  assert xray.ext.gradient_flags(
    False, False, False, False, False, False, False, False).all_false()
  f.u_iso = True
  assert f.adjust(False).u_iso == True
  assert f.adjust(True).u_iso == False
  assert f.adjust(False).u_aniso == False
  assert f.adjust(True).u_aniso == True

def exercise_refinement_flags():
  f = xray.refinement_flags(
    site=False, u_iso=True, u_aniso=False, occupancy=True,
    fp=False, fdp=True, tan_u_iso=False, param=42)
  assert not f.site()
  assert f.u_iso()
  assert not f.u_aniso()
  assert f.occupancy()
  assert not f.fp()
  assert f.fdp()
  assert not f.tan_u_iso()
  assert f.param == 42
  #
  f = xray.refinement_flags(
    site=True, u_iso=False, u_aniso=True, occupancy=False,
    fp=True, fdp=False, tan_u_iso=True, param=-24)
  assert f.site()
  assert not f.u_iso()
  assert f.u_aniso()
  assert not f.occupancy()
  assert f.fp()
  assert not f.fdp()
  assert f.tan_u_iso()
  assert f.param == -24
  #
  for state in [False, True]:
    f.set_site(state=state)
    assert f.site() == state
    f.set_u_iso(state=state)
    assert f.u_iso() == state
    f.set_u_aniso(state=state)
    assert f.u_aniso() == state
    f.set_occupancy(state=state)
    assert f.occupancy() == state
    f.set_fp(state=state)
    assert f.fp() == state
    f.set_fdp(state=state)
    assert f.fdp() == state
  f.param = 35
  assert f.param == 35

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
  assert not x.refinement_flags.site()
  x.refinement_flags.set_site(state=True)
  assert x.refinement_flags.site()
  #
  x = xray.scatterer(
    "si1", site=(0.01,0.02,0.3), occupancy=0.9, u=(0.3, 0.3, 0.2, 0,0,0))
  assert x.scattering_type == "Si"
  uc = uctbx.unit_cell((10, 10, 13))
  sg = sgtbx.space_group_info("P 4")
  ss = x.apply_symmetry(unit_cell=uc, space_group=sg.group())
  assert x.multiplicity() == 1
  assert approx_equal(x.weight_without_occupancy(), 1/4.)
  assert approx_equal(x.weight(), 0.9/4.)
  assert approx_equal(x.site, (0,0,0.3))
  assert ss.multiplicity() == x.multiplicity()
  x.occupancy = 0.8
  assert approx_equal(x.weight(), 0.8/4.)
  u_cart = (0.3354, 0.3771, 0.4874, -0.05161, 0.026763, -0.02116)
  x.u_star = adptbx.u_cart_as_u_star(uc, u_cart)
  x.anisotropic_flag = 1
  try:
    x.apply_symmetry(uc, sg.group(), u_star_tolerance=0.1)
  except RuntimeError, e:
    assert str(e).find("is_compatible_u_star") > 0
  else:
    raise AssertionError("Exception expected.")
  x.apply_symmetry(site_symmetry_ops=ss)
  x.apply_symmetry(site_symmetry_ops=ss, u_star_tolerance=0.5)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0, 0)
  ss = x.apply_symmetry(
    unit_cell=uc,
    space_group=sg.group(),
    min_distance_sym_equiv=0.5,
    u_star_tolerance=0,
    assert_min_distance_sym_equiv=False)
  assert ss.is_compatible_u_star(x.u_star)
  assert approx_equal(x.u_star, (0.0035625, 0.0035625, 0.002884, 0, 0, 0))
  site = (0.2,0.5,0.4)
  x.apply_symmetry_site(ss)
  assert approx_equal(x.site, (0,0,0.3))
  x.u_star = (1,2,3,4,5,6)
  x.apply_symmetry_u_star(
    site_symmetry_ops=ss,
    u_star_tolerance=0)
  assert approx_equal(x.u_star, (1.5,1.5,3.0,0,0,0))
  x.site = (0.2,0.5,0.4)
  ss = x.apply_symmetry(uc, sg.group(), 1.e-10, 0)
  assert ss.is_point_group_1()
  assert x.anisotropic_flag
  x.convert_to_isotropic(unit_cell=uc)
  assert not x.anisotropic_flag
  assert approx_equal(x.u_iso, 269)
  assert approx_equal(x.u_star, (-1,-1,-1,-1,-1,-1))
  x.convert_to_anisotropic(unit_cell=uc)
  assert x.anisotropic_flag
  assert approx_equal(x.u_iso, -1)
  assert approx_equal(x.u_star, (2.69, 2.69, 1.59171598, 0, 0, 0))
  x.u_star = (1,2,3,4,5,6)
  assert not x.is_positive_definite_u(unit_cell=uc)
  assert not x.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=1.e2)
  assert x.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=1.e3)
  x.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=0)
  assert approx_equal(x.u_star,
    (3.3379643647809192, 4.5640522609325131, 4.4690204772593507,
     3.9031581835726965, 3.8623090371651934, 4.5162864184404032))
  x.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=1)
  assert approx_equal(x.u_star,
    (3.3458045216665266, 4.5710990727698393, 4.4720459395534728,
     3.9006326295505751, 3.8598099147456764, 4.5133641373560351))
  assert x.is_positive_definite_u(unit_cell=uc)
  y = x.customized_copy(u=-1)
  assert not y.anisotropic_flag
  assert approx_equal(y.u_iso, -1)
  assert not y.is_positive_definite_u(unit_cell=uc)
  assert not y.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=0.5)
  assert y.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=2)
  a = flex.xray_scatterer([x,y])
  assert list(xray.is_positive_definite_u(
    scatterers=a, unit_cell=uc)) == [True, False]
  a = flex.xray_scatterer([y,x])
  assert list(xray.is_positive_definite_u(
    scatterers=a, unit_cell=uc)) == [False, True]
  assert list(xray.is_positive_definite_u(
    scatterers=a, unit_cell=uc, u_cart_tolerance=2)) == [True, True]
  y.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=1)
  assert approx_equal(y.u_iso, 1)
  assert y.is_positive_definite_u(unit_cell=uc)
  x_u_star_orig = x.u_star
  x.shift_u(unit_cell=uc, u_shift=10)
  assert approx_equal(x.u_star,
    (3.4458045216665267, 4.6710990727698389, 4.5312175371866088,
     3.9006326295505751, 3.8598099147456764, 4.5133641373560351))
  y.shift_u(unit_cell=uc, u_shift=10)
  assert approx_equal(y.u_iso, 11)
  a = flex.xray_scatterer([x,y])
  xray.shift_us(scatterers=a, unit_cell=uc, u_shift=-10)
  assert approx_equal(a[0].u_star, x_u_star_orig)
  assert approx_equal(a[1].u_iso, 1)
  a[0].fp = 3;
  a[1].fdp = 4;
  assert not show_diff(a[0].report_details(unit_cell=uc, prefix="&%"), """\
&%scatterer label: si1
&%scattering type: Si
&%fractional coordinates: 0.200000 0.500000 0.400000
&%cartesian coordinates: 2.000000 5.000000 5.200000
&%u_star: 3.3458 4.5711 4.47205 3.90063 3.85981 4.51336
&%u_cart: 334.58 457.11 755.776 390.063 501.775 586.737
&%occupancy: 0.8
&%f-prime: 3
&%f-double-prime: 0""")
  assert not show_diff(a[1].report_details(unit_cell=uc, prefix="=#"), """\
=#scatterer label: si1
=#scattering type: Si
=#fractional coordinates: 0.200000 0.500000 0.400000
=#cartesian coordinates: 2.000000 5.000000 5.200000
=#u_iso: 1
=#b_iso: 78.9568
=#occupancy: 0.8
=#f-prime: 0
=#f-double-prime: 4""")

def exercise_rotate():
  uc = uctbx.unit_cell((10, 10, 13))
  s = flex.xray_scatterer((xray.scatterer("Si1", site=(0.01,0.02,0.3)),))
  r = xray.rotate(
    unit_cell=uc,
    rotation_matrix=((1,0,0, 0,1,0, 0,0,1)),
    scatterers=s)
  assert r.size() == 1
  assert approx_equal(s[0].site, r[0].site)
  r = xray.rotate(
    unit_cell=uc,
    rotation_matrix=((0,-1,0, -1,0,0, 0,0,-1)),
    scatterers=s)
  assert approx_equal(r[0].site, (-0.02,-0.01,-0.3))

def exercise_scattering_type_registry():
  reg = xray.scattering_type_registry()
  assert len(reg.type_index_pairs_as_dict()) == 0
  assert len(reg.unique_gaussians_as_list()) == 0
  assert reg.unique_counts.size() == 0
  assert reg.size() == 0
  assert reg.process(scattering_type="const") == 0
  assert reg.size() == 1
  assert reg.type_index_pairs_as_dict() == {"const": 0}
  assert list(reg.unique_counts) == [1]
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
  unique_indices = reg.process(scatterers=scatterers)
  assert list(unique_indices) == [1,1,2,2,3,2,3,0,4]
  assert reg.unique_indices(scatterers=scatterers).all_eq(unique_indices)
  assert reg.type_index_pairs_as_dict() \
      == {"const": 0, "Al": 3, "O": 2, "custom": 4, "Si": 1}
  assert reg.unique_gaussians_as_list().count(None) == 5
  assert list(reg.unique_counts) == [2,2,3,2,1]
  for t,i in reg.type_index_pairs_as_dict().items():
    assert reg.unique_index(scattering_type=t) == i
    assert reg.gaussian(scattering_type=t) is None
  assert list(reg.unassigned_types()) == ['Al', 'O', 'Si', 'const', 'custom']
  assert reg.assign(
    scattering_type="const",
    gaussian=eltbx.xray_scattering.gaussian(10))
  assert reg.unique_gaussians_as_list().count(None) == 4
  assert reg.unique_gaussians_as_list()[0].n_parameters() == 1
  assert not reg.assign("const", eltbx.xray_scattering.gaussian(10))
  assert reg.gaussian(scattering_type="const").n_parameters() == 1
  assert reg.gaussian_not_optional(scattering_type="const").n_parameters() == 1
  assert reg.assign("custom", eltbx.xray_scattering.gaussian((1,2),(3,4),5))
  assert reg.unique_gaussians_as_list().count(None) == 3
  assert reg.unique_gaussians_as_list()[4].n_parameters() == 5
  assert list(reg.unassigned_types()) == ['Al', 'O', 'Si']
  for table,n_terms in (("IT1992",4), ("WK1995",5)):
    reg = xray.scattering_type_registry()
    reg.process(scatterers=scatterers)
    reg.assign("const", eltbx.xray_scattering.gaussian(10))
    reg.assign("custom", eltbx.xray_scattering.gaussian((1,2),(3,4),5))
    reg.assign_from_table(table=table)
    ugs = reg.unique_gaussians_as_list()
    for t,i in reg.type_index_pairs_as_dict().items():
      if (t in ["Al", "O", "Si"]):
        assert ugs[i].n_terms() == n_terms
      elif (t == "const"):
        assert ugs[i].n_terms() == 0
        assert approx_equal(ugs[i].c(), 10)
      else:
        assert ugs[i].n_terms() == 2
        assert approx_equal(ugs[i].c(), 5)
    reg.assign("Al", eltbx.xray_scattering.gaussian(20))
    ugs = reg.unique_gaussians_as_list()
    assert approx_equal(ugs[reg.unique_index("Al")].c(), 20)
    assert reg.unassigned_types().size() == 0
    ff = reg.unique_form_factors_at_d_star_sq(d_star_sq=0)
    if (n_terms == 4):
      assert approx_equal(ff, [13.9976, 7.9994, 20.0, 10.0, 8.0])
    else:
      assert approx_equal(ff, [13.998917, 7.999706, 20.0, 10.0, 8.0])
    ff = reg.unique_form_factors_at_d_star_sq(d_star_sq=0.123)
    if (n_terms == 4):
      assert approx_equal(ff, [10.173729, 6.042655, 20.0, 10.0, 7.680404])
    else:
      assert approx_equal(ff, [10.174771, 6.042745, 20.0, 10.0, 7.680404])
    reg.assign(scattering_type="custom", gaussian=None)
    assert reg.gaussian(scattering_type="custom") is None
    s = StringIO()
    reg.show_summary(out=s, prefix="=-")
    if (n_terms == 4):
      assert not show_diff(s.getvalue(),
        "=-Al:0+c*2 Si:4+c*2 const:0+c*1 O:4+c*3 custom:None*1\n")
    else:
      assert not show_diff(s.getvalue(),
        "=-Al:0+c*2 Si:5+c*2 const:0+c*1 O:5+c*3 custom:None*1\n")
    s = StringIO()
    for show_weights in [False, True]:
      for show_gaussians in [False, True]:
        reg.show(
          show_weights=show_weights,
          show_gaussians=show_gaussians,
          out=s,
          prefix=":#")
    assert not show_diff(s.getvalue(), """\
:#Number of scattering types: 5
:#  Type    Number
:#   Al         2
:#   Si         2
:#   const      1
:#   O          3
:#   custom     1
:#Number of scattering types: 5
:#  Type    Number   Gaussians
:#   Al         2        0+c
:#   Si         2        %(n_terms)d+c
:#   const      1        0+c
:#   O          3        %(n_terms)d+c
:#   custom     1       None
:#Number of scattering types: 5
:#  Type    Number   Weight
:#   Al         2     20.00
:#   Si         2     14.00
:#   const      1     10.00
:#   O          3      8.00
:#   custom     1      None
:#Number of scattering types: 5
:#  Type    Number   Weight   Gaussians
:#   Al         2     20.00       0+c
:#   Si         2     14.00       %(n_terms)d+c
:#   const      1     10.00       0+c
:#   O          3      8.00       %(n_terms)d+c
:#   custom     1      None      None
""" % vars())
    assert reg.wilson_dict() \
        == {'Si': 2, 'const': 1, 'Al': 2, 'O': 3, 'custom': 1}
    type_index_pairs = reg.type_index_pairs_as_dict()
    unique_gaussians = reg.unique_gaussians_as_list()
    unique_counts = reg.unique_counts
    reg = xray.scattering_type_registry(
      type_index_pairs, unique_gaussians, unique_counts)
    assert reg.type_index_pairs_as_dict() == type_index_pairs
    for orig,restored in zip(unique_gaussians, reg.unique_gaussians_as_list()):
      if (restored is None):
        assert orig is None
        continue
      assert restored is not orig
      assert restored.n_parameters() == orig.n_parameters()
      for d_star_sq in [0, 0.1, 0.234]:
        assert approx_equal(
          restored.at_d_star_sq(d_star_sq), orig.at_d_star_sq(d_star_sq))
    assert reg.unique_counts.all_eq(unique_counts)
    s = pickle.dumps(reg)
    l = pickle.loads(s)
    orig = StringIO()
    restored = StringIO()
    reg.show(out=orig)
    l.show(out=restored)
    assert not show_diff(restored.getvalue(), orig.getvalue())
  try: reg.unique_index("foo")
  except RuntimeError, e:
    assert str(e) == 'scattering_type "foo" not in scattering_type_registry.'
  else: raise RuntimeError("Exception expected.")
  try: reg.gaussian_not_optional(scattering_type="custom")
  except RuntimeError, e:
    assert str(e) == 'gaussian not defined for scattering_type "custom".'
  else: raise RuntimeError("Exception expected.")
  try: reg.unique_form_factors_at_d_star_sq(d_star_sq=0)
  except RuntimeError, e:
    assert str(e) == 'gaussian not defined for scattering_type "custom".'
  else: raise RuntimeError("Exception expected.")

def exercise_structure_factors():
  uc = uctbx.unit_cell((10, 10, 13))
  sg = sgtbx.space_group_info("P 4")
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3)),
    xray.scatterer("O1", site=(0.3,0.4,0.5), u=(0.4,0.5,0.6,-.05,0.2,-0.02))))
  for s in scatterers:
    assert s.multiplicity() == 0
  assert xray.n_undefined_multiplicities(scatterers) == 2
  site_symmetry_table = sgtbx.site_symmetry_table()
  xray.add_scatterers_ext(
    unit_cell=uc,
    space_group=sg.group(),
    scatterers=scatterers,
    site_symmetry_table=site_symmetry_table,
    site_symmetry_table_for_new=sgtbx.site_symmetry_table(),
    min_distance_sym_equiv=0.5,
    u_star_tolerance=0,
    assert_min_distance_sym_equiv=True)
  assert list(site_symmetry_table.special_position_indices()) == [0]
  xray.tidy_us(
    scatterers=scatterers,
    unit_cell=uc,
    site_symmetry_table=site_symmetry_table,
    u_min=0)
  assert approx_equal(scatterers[0].u_iso, 0)
  assert approx_equal(scatterers[1].u_star, (0.4,0.5,0.6,-.05,0.2,-0.02))
  for s in scatterers:
    assert s.multiplicity() != 0
  assert xray.n_undefined_multiplicities(scatterers) == 0
  mi = flex.miller_index(((1,2,3), (2,3,4)))
  scattering_type_registry = xray.scattering_type_registry()
  scattering_type_registry.process(scatterers=scatterers)
  scattering_type_registry.assign_from_table("WK1995")
  fc = xray.ext.structure_factors_simple(
    uc, sg.group(), mi, scatterers, scattering_type_registry).f_calc()
  assert approx_equal(flex.abs(fc), (10.50871, 9.049631))
  assert approx_equal(flex.arg(fc, 1), (-36, 72))
  assert approx_equal(flex.abs(fc), (10.50871, 9.049631))
  assert approx_equal(flex.arg(fc, 1), (-36, 72))
  fc = xray.ext.structure_factors_direct(
    uc, sg.group(), mi, scatterers, scattering_type_registry).f_calc()
  assert approx_equal(flex.abs(fc), (10.50871, 9.049631))
  assert approx_equal(flex.arg(fc, 1), (-36, 72))
  xray.tidy_us(
    scatterers=scatterers,
    unit_cell=uc,
    site_symmetry_table=site_symmetry_table,
    u_min=100)
  assert approx_equal(scatterers[0].u_iso, 100)
  assert approx_equal(scatterers[1].u_star,
    (1.0134539945616343, 1.0005190241807682, 0.64980451464405997,
     -0.0026425269166861672, 0.027955730692513142, -0.0054908429234285239))
  xray.ext.structure_factors_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, scattering_type_registry).f_calc()
  xray.ext.structure_factors_gradients_direct(
    uc, sg.group(), mi, scatterers, None,
    scattering_type_registry, site_symmetry_table,
    flex.complex_double(mi.size()),
    xray.ext.gradient_flags(
      True, True, True, True, True, True, False, False),
    0)
  xray.ext.structure_factors_gradients_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, None,
    scattering_type_registry, site_symmetry_table,
    flex.complex_double(mi.size()),
    xray.ext.gradient_flags(
      True, True, True, True, True, True, False, False),
    0)

def exercise_targets():
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, True)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, False, 3)
  assert approx_equal(ls.scale_factor(), 3)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual(f_obs, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double((10,20,30,40,50))
  ls = xray.targets_least_squares_residual(f_obs, f_calc, True)
  assert approx_equal(ls.scale_factor(), 1/10.)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual(f_obs, f_calc, False, 3/10.)
  assert approx_equal(ls.scale_factor(), 3/10.)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j,5-4j,-5+6j))
  w = flex.double((1,2,3,2,4))
  ls = xray.targets_least_squares_residual(f_obs, w, f_calc, True)
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
  ic = xray.targets_intensity_correlation(f_obs, w, f_calc, True)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert approx_equal(tuple(ic.derivatives()), (0j,0j,0j,0j,0j))
  ic = xray.targets_intensity_correlation(f_obs, f_calc)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert ic.derivatives().size() == 0
  f_calc = flex.complex_double((10,20,30,40,50))
  ic = xray.targets_intensity_correlation(f_obs, f_calc, True)
  assert approx_equal(ic.correlation(), 1)
  assert approx_equal(ic.target(), 0)
  assert approx_equal(tuple(ic.derivatives()), (0j,0j,0j,0j,0j))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j,5-4j,-5+6j))
  w = flex.int((1,2,3,2,4))
  ic = xray.targets_intensity_correlation(f_obs, w, f_calc, True)
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
  scattering_type_registry = xray.ext.scattering_type_registry()
  scattering_type_registry.process(scatterers)
  scattering_type_registry.assign_from_table("WK1995")
  d = xray.sampled_model_density(
    unit_cell=uc,
    scatterers=scatterers,
    scattering_type_registry=scattering_type_registry,
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
  assert d.excessive_sampling_radius_i_seqs().size() == 0
  i = flex.miller_index(((1,2,3), (2,3,4)))
  f = flex.complex_double((1+2j, 2+3j))
  d.eliminate_u_extra_and_normalize(i, f)
  f_orig = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), 0.2, i, f)
  f_elim = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), -0.2, i, f, 1)
  assert approx_equal(f, f_orig)

def exercise_minimization_apply_shifts():
  uc = uctbx.unit_cell((20, 20, 23))
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  f = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  s = xray.ext.minimization_shift_scales(
    scatterers=scatterers,
    gradient_flags=f,
    n_parameters=19,
    site_cart=1,
    u_iso=2,
    u_cart=3,
    occupancy=4,
    fp=5,
    fdp=6)
  assert [int(x) for x in s] \
      == [1, 1, 1, 2, 4, 5, 6, 1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 5, 6]
  shifts = flex.double(19, 0.001)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  for a,b in zip(scatterers, shifted_scatterers):
    assert a.scattering_type == b.scattering_type
    assert a.site != b.site
    if (not a.anisotropic_flag):
      assert a.u_iso != b.u_iso
      assert approx_equal(a.u_star, b.u_star)
    else:
      assert a.u_iso == b.u_iso
      assert not_approx_equal(a.u_star, b.u_star)
    assert a.occupancy != b.occupancy
    assert a.fp != b.fp
    assert a.fdp != b.fdp
  shifts = flex.double(19, -0.001)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    unit_cell=uc,
    scatterers=shifted_scatterers,
    gradient_flags=f,
    shifts=shifts).shifted_scatterers
  for a,b in zip(scatterers, shifted_scatterers):
    assert a.scattering_type == b.scattering_type
    assert approx_equal(a.site, b.site)
    assert approx_equal(a.u_iso, b.u_iso)
    assert approx_equal(a.u_star, b.u_star)
    assert approx_equal(a.occupancy, b.occupancy)
    assert approx_equal(a.fp, b.fp)
    assert approx_equal(a.fdp, b.fdp)
  f = xray.ext.gradient_flags(
    True, False, False, False, False, False, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 6, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [1]*6
  shifts = flex.double((-1,2,-3,4,-5,-6))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  assert approx_equal(
    shifted_scatterers[0].site,
    (0.01-1/20.,0.02+2/20.,0.3-3/23.))
  assert approx_equal(
    shifted_scatterers[1].site,
    (0.3+4/20.,0.4-5/20.,0.5-6/23.))
  f = xray.ext.gradient_flags(
    False, True, False, False, False, False, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 1, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [2]
  shifts = flex.double(1, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  assert approx_equal(shifted_scatterers[0].u_iso, -10)
  f = xray.ext.gradient_flags(
    False, False, True, False, False, False, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 6, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [3]*6
  shifts = flex.double(6, -100)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  assert not_approx_equal(shifted_scatterers[1].u_star,
    [u_ij-100 for u_ij in scatterers[1].u_star])
  f = xray.ext.gradient_flags(
    False, False, False, True, False, False, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [4]*2
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  for i in xrange(2):
    assert approx_equal(shifted_scatterers[i].occupancy, -9)
  f = xray.ext.gradient_flags(
    False, False, False, False, True, False, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [5]*2
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  assert approx_equal(shifted_scatterers[0].fp, -11)
  assert shifted_scatterers[1].fp == -10
  for i in xrange(2):
    assert shifted_scatterers[i].fdp == scatterers[i].fdp
  f = xray.ext.gradient_flags(
    False, False, False, False, False, True, False, False)
  s = xray.ext.minimization_shift_scales(scatterers, f, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [6]*2
  shifts = flex.double((2,3))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, f, shifts).shifted_scatterers
  assert shifted_scatterers[0].fp == -1
  assert approx_equal(shifted_scatterers[0].fdp, 4)
  assert shifted_scatterers[1].fp == 0
  assert shifted_scatterers[1].fdp == 3
  shifts = flex.double(1)
  try:
    xray.ext.minimization_apply_shifts(uc, scatterers, f, shifts)
  except Exception, e:
    assert str(e) == "scitbx Error: Array of shifts is too small."
  else:
    raise RuntimeError("Exception expected.")
  shifts = flex.double(3)
  try:
    xray.ext.minimization_apply_shifts(uc, scatterers, f, shifts)
  except Exception, e:
    assert str(e) == "cctbx Error: Array of shifts is too large."
  else:
    raise RuntimeError("Exception expected.")

def exercise_minimization_add_gradients():
  uc = uctbx.unit_cell((20, 20, 23))
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  gradient_flags = xray.ext.gradient_flags(
    True, False, False, False, False, False, False, False)
  xray_gradients = flex.double(xrange(6))
  geometry_restraints_site_gradients = flex.vec3_double([(1,-2,3),(-4,-5,6)])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,-1,-1,11])
  gradient_flags = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  xray_gradients = flex.double(xrange(19))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,5,6,3,3,15,10,11,12,13,14,15,16,17,18])
  gradient_flags = xray.ext.gradient_flags(
    True, True, False, True, True, True, False, False)
  xray_gradients = flex.double(xrange(13))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,5,6,3,3,15,10,11,12])
  gradient_flags = xray.ext.gradient_flags(
    True, False, True, True, False, True, False, False)
  xray_gradients = flex.double(xrange(16))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,1,1,13,8,9,10,11,12,13,14,15])
  site_gradients = xray.ext.minimization_extract_site_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients)
  assert approx_equal(site_gradients, [(1,-1,5), (1,1,13)])
  #
  gradient_flags = xray.ext.gradient_flags(
    False, True, False, False, False, False, False, False)
  xray_gradients = flex.double([1])
  u_iso_gradients = flex.double([3,0])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients, [4])
  gradient_flags = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  xray_gradients = flex.double(xrange(19))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [0,1,2,6,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
  gradient_flags = xray.ext.gradient_flags(
    True, True, False, True, True, True, False, False)
  xray_gradients = flex.double(xrange(13))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [0,1,2,6,4,5,6,7,8,9,10,11,12])
  gradient_flags = xray.ext.gradient_flags(
    False, True, True, True, False, True, False, False)
  xray_gradients = flex.double(xrange(11))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [3,1,2,3,4,5,6,7,8,9,10])
  #
  gradient_flags = xray.ext.gradient_flags(
    False, False, False, True, False, False, False, False)
  xray_gradients = flex.double([1,2])
  occupancy_gradients = flex.double([3,4])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=None,
    occupancy_gradients=occupancy_gradients)
  assert approx_equal(xray_gradients, [4,6])
  xray_gradients = flex.double([2,1])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    gradient_flags=gradient_flags,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients, [2,1])

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
    scatterers=structure.scatterers(),
    site_symmetry_table=structure.site_symmetry_table())
  assert asu_mappings.mappings().size() == structure.scatterers().size()

def exercise_ls_target_with_scale_k1():
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = False,
                                    fix_scale           = False)
  assert approx_equal(ls.scale(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = False,
                                    fix_scale           = True,
                                    scale               = 2.0)
  assert approx_equal(ls.scale(), 2.0)
  assert approx_equal(ls.target(),1.0)
  assert ls.derivatives().size() == 0
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((4,5,6))
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = False)
  assert approx_equal(ls.scale(), 50./134.)
  assert approx_equal(ls.target(), 0.0671641791)
  assert approx_equal( tuple(ls.derivatives()),
                     ((0.0551347738+0j),(-0.0100245043+0j),(-0.0284027623+0j)) )
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = True,
                                    scale               = 2.0)
  assert approx_equal(ls.scale(), 2.0)
  assert approx_equal(ls.target(),17.8)
  assert approx_equal( tuple(ls.derivatives()), ((4.2+0j),(3.2+0j),(1.8+0j)) )
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j))
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = False)
  assert approx_equal(ls.scale(), 0.4773772552)
  assert approx_equal(ls.target(), 0.2023883467)
  assert approx_equal( tuple(ls.derivatives()),
                       ((0.0043198335244903152-0.0086396670489806305j),
                        (0.022162885026120613-0.029550513368160818j),
                        (0.041257975234691303-0.082515950469382607j)) )

def exercise_ls_target_with_scale_k2():
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = False,
                                    fix_scale           = False)
  assert approx_equal(ls.scale(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  #
  ls = xray.ls_target_with_scale_k1(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = False,
                                    fix_scale           = True,
                                    scale               = 2.0)
  assert approx_equal(ls.scale(), 2.0)
  assert approx_equal(ls.target(),1.0)
  assert ls.derivatives().size() == 0
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((4,5,6))
  #
  ls = xray.ls_target_with_scale_k2(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = False)
  assert approx_equal(ls.scale(), 50./20.)
  assert approx_equal(ls.target(), 0.45)
  assert approx_equal( tuple(ls.derivatives()),
                     ((0.45000000000000001+0j), 0j, (-0.15000000000000002+0j)) )
  #
  ls = xray.ls_target_with_scale_k2(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = True,
                                    scale               = 2.0)
  assert approx_equal(ls.scale(), 2.0)
  assert approx_equal(ls.target(),0.7)
  assert approx_equal( tuple(ls.derivatives()), ((0.6+0j),(0.2+0j),(0.0+0j)) )
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j))
  #
  ls = xray.ls_target_with_scale_k2(f_obs               = f_obs,
                                    weights             = w,
                                    f_calc              = f_calc,
                                    compute_derivatives = True,
                                    fix_scale           = False)
  scale = flex.sum(w*flex.abs(f_calc)*f_obs)/flex.sum(w*f_obs*f_obs)
  assert approx_equal(ls.scale(), 1.6708203932)
  assert approx_equal(ls.target(), 0.7083592135)
  assert approx_equal( tuple(ls.derivatives()),
                       ((0.075835921350012631-0.15167184270002526j),
                        (0.19900310562001516-0.26533747416002024j),
                        (0.12416407864998737-0.24832815729997473j)) )

def run():
  exercise_conversions()
  exercise_gradient_flags()
  exercise_refinement_flags()
  exercise_xray_scatterer()
  exercise_scattering_type_registry()
  exercise_rotate()
  exercise_structure_factors()
  exercise_targets()
  exercise_sampled_model_density()
  exercise_minimization_apply_shifts()
  exercise_minimization_add_gradients()
  exercise_asu_mappings()
  exercise_ls_target_with_scale_k1()
  exercise_ls_target_with_scale_k2()
  print "OK"

if (__name__ == "__main__"):
  run()
