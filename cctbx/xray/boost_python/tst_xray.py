from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import eltbx
from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import xray
from cctbx import math_module
from cctbx.array_family import flex
from scitbx.array_family import shared
from libtbx.test_utils import Exception_expected, approx_equal, \
  not_approx_equal, show_diff
from cStringIO import StringIO
import pickle

def exercise_scatterer_flags():
  f = xray.scatterer_flags()
  assert f.use()                      == True
  assert f.use_u_iso()                == False
  assert f.use_u_aniso()              == False
  assert f.grad_site()                == False
  assert f.grad_u_iso()               == False
  assert f.grad_u_aniso()             == False
  assert f.grad_occupancy()           == False
  assert f.grad_fp()                  == False
  assert f.grad_fdp()                 == False
  assert f.curv_site_site()           == False
  assert f.curv_site_u_iso()          == False
  assert f.curv_site_u_aniso()        == False
  assert f.curv_site_occupancy()      == False
  assert f.curv_site_fp()             == False
  assert f.curv_site_fdp()            == False
  assert f.curv_u_iso_u_iso()         == False
  assert f.curv_u_iso_u_aniso()       == False
  assert f.curv_u_iso_occupancy()     == False
  assert f.curv_u_iso_fp()            == False
  assert f.curv_u_iso_fdp()           == False
  assert f.curv_u_aniso_u_aniso()     == False
  assert f.curv_u_aniso_occupancy()   == False
  assert f.curv_u_aniso_fp()          == False
  assert f.curv_u_aniso_fdp()         == False
  assert f.curv_occupancy_occupancy() == False
  assert f.curv_occupancy_fp()        == False
  assert f.curv_occupancy_fdp()       == False
  assert f.curv_fp_fp()               == False
  assert f.curv_fp_fdp()              == False
  assert f.curv_fdp_fdp()             == False
  assert f.tan_u_iso()                == False
  assert f.param                      == 0
  for state in [True, False]:
    f.set_use(state)
    f.set_use_u_iso(state)
    f.set_use_u_aniso(state)
    f.set_grad_site(state)
    f.set_grad_u_iso(state)
    f.set_grad_u_aniso(state)
    f.set_grad_occupancy(state)
    f.set_grad_fp(state)
    f.set_grad_fdp(state)
    f.set_curv_site_site(state)
    f.set_curv_site_u_iso(state)
    f.set_curv_site_u_aniso(state)
    f.set_curv_site_occupancy(state)
    f.set_curv_site_fp(state)
    f.set_curv_site_fdp(state)
    f.set_curv_u_iso_u_iso(state)
    f.set_curv_u_iso_u_aniso(state)
    f.set_curv_u_iso_occupancy(state)
    f.set_curv_u_iso_fp(state)
    f.set_curv_u_iso_fdp(state)
    f.set_curv_u_aniso_u_aniso(state)
    f.set_curv_u_aniso_occupancy(state)
    f.set_curv_u_aniso_fp(state)
    f.set_curv_u_aniso_fdp(state)
    f.set_curv_occupancy_occupancy(state)
    f.set_curv_occupancy_fp(state)
    f.set_curv_occupancy_fdp(state)
    f.set_curv_fp_fp(state)
    f.set_curv_fp_fdp(state)
    f.set_curv_fdp_fdp(state)
    f.set_tan_u_iso(state)
    f.param = 42
    assert f.use()                      == state
    assert f.use_u_iso()                == state
    assert f.use_u_aniso()              == state
    assert f.grad_site()                == state
    assert f.grad_u_iso()               == state
    assert f.grad_u_aniso()             == state
    assert f.grad_occupancy()           == state
    assert f.grad_fp()                  == state
    assert f.grad_fdp()                 == state
    assert f.curv_site_site()           == state
    assert f.curv_site_u_iso()          == state
    assert f.curv_site_u_aniso()        == state
    assert f.curv_site_occupancy()      == state
    assert f.curv_site_fp()             == state
    assert f.curv_site_fdp()            == state
    assert f.curv_u_iso_u_iso()         == state
    assert f.curv_u_iso_u_aniso()       == state
    assert f.curv_u_iso_occupancy()     == state
    assert f.curv_u_iso_fp()            == state
    assert f.curv_u_iso_fdp()           == state
    assert f.curv_u_aniso_u_aniso()     == state
    assert f.curv_u_aniso_occupancy()   == state
    assert f.curv_u_aniso_fp()          == state
    assert f.curv_u_aniso_fdp()         == state
    assert f.curv_occupancy_occupancy() == state
    assert f.curv_occupancy_fp()        == state
    assert f.curv_occupancy_fdp()       == state
    assert f.curv_fp_fp()               == state
    assert f.curv_fp_fdp()              == state
    assert f.curv_fdp_fdp()             == state
    assert f.tan_u_iso()                == state
    assert f.param                      == 42

  f.set_use_u_iso(state=False)
  f.set_use_u_aniso(state=True)
  f.set_use_u_iso_only()
  assert f.use_u_iso()
  assert not f.use_u_aniso()
  assert f.use_u_iso_only()
  assert not f.use_u_aniso_only()
  f.set_use_u_aniso_only()
  assert not f.use_u_iso()
  assert f.use_u_aniso()
  assert not f.use_u_iso_only()
  assert f.use_u_aniso_only()
  f.set_use_u(iso=True, aniso=True)
  try: f.use_u_iso_only()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "scatterer.flags.u_iso_only(): u_iso and u_aniso both true.")
  else: raise Exception_expected
  try: f.use_u_aniso_only()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "scatterer.flags.u_aniso_only(): u_iso and u_aniso both true.")
  else: raise Exception_expected
  f.set_use_u(iso=False, aniso=False)
  try: f.use_u_iso_only()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "scatterer.flags.u_iso_only(): u_iso and u_aniso both false.")
  else: raise Exception_expected
  try: f.use_u_aniso_only()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "scatterer.flags.u_aniso_only(): u_iso and u_aniso both false.")
  else: raise Exception_expected

  f1 = xray.scatterer_flags()
  f1.set_grad_site(True)
  f1.set_grad_u_aniso(True)
  f1.set_grad_fp(True)
  f2 = xray.scatterer_flags()
  assert f2.implies(f1)
  f2.set_grad_site(True)
  f2.set_grad_u_aniso(True)
  f2.set_grad_fp(True)
  assert f2.implies(f1)
  f2.set_grad_u_iso(True)
  assert not f2.implies(f1)

  flags = xray.shared_scatterer_flags()
  f1 = xray.scatterer_flags()
  f1.set_use_u_iso(True)
  f1.set_grad_u_iso(True)
  f1.set_grad_fp(True)
  flags.append(f1)
  f2 = xray.scatterer_flags()
  f2.set_grad_site(True)
  f2.set_grad_u_aniso(True)
  flags.append(f2)
  assert flags[0].bits == f1.bits
  assert flags[1].bits == f2.bits
  for f in flags:
    f.set_grad_occupancy(True)
  assert flags[0].grad_occupancy()
  assert flags[1].grad_occupancy()
  assert flags.n_parameters() == 7

  x = xray.scatterer("C", site=(0.4, 0.5, 0.6))
  scatterers = flex.xray_scatterer(10, x)
  for sc in scatterers:
    sc.flags.set_grad_site(True)
    sc.flags.set_use_u_aniso(True)
    sc.flags.set_grad_u_aniso(True)
  grad_flags = xray.shared_scatterer_flags(scatterers)
  assert [ f.bits for f in grad_flags ] \
      == [ sc.flags.bits for sc in scatterers ]

def exercise_set_scatterer_grad_flags():
  x = xray.scatterer("c", site=(0.1,0.2,0.3), occupancy=0.0, u=(0,0,0,0,0,0))
  y = xray.scatterer("c", site=(0.1,0.2,0.3), occupancy=0.0, u=1.0)
  scatterers = flex.xray_scatterer(10, x)
  scatterers.extend(flex.xray_scatterer(10, y))
  for scatterer in scatterers:
      assert scatterer.flags.grad_site()      == False
      assert scatterer.flags.grad_u_iso()     == False
      assert scatterer.flags.grad_u_aniso()   == False
      assert scatterer.flags.grad_occupancy() == False
      assert scatterer.flags.grad_fp()        == False
      assert scatterer.flags.grad_fdp()       == False
      assert scatterer.flags.tan_u_iso()      == False
      assert scatterer.flags.param            == 0
  xray.set_scatterer_grad_flags(scatterers = scatterers)
  for scatterer in scatterers:
      assert scatterer.flags.grad_site()      == False
      assert scatterer.flags.grad_u_iso()     == False
      assert scatterer.flags.grad_u_aniso()   == False
      assert scatterer.flags.grad_occupancy() == False
      assert scatterer.flags.grad_fp()        == False
      assert scatterer.flags.grad_fdp()       == False
      assert scatterer.flags.tan_u_iso()      == False
      assert scatterer.flags.param            == 0
  cn1 = 0
  cn2 = 0
  for site in [True, False]:
    for u_iso in [True, False]:
      for u_aniso in [True, False]:
        for occupancy in [True, False]:
          for fp in [True, False]:
            for fdp in [True, False]:
              for tan_u_iso in [True, False]:
                for param in [0, 1, 2]:
                    xray.set_scatterer_grad_flags(scatterers = scatterers,
                                                  site       = site,
                                                  u_iso      = u_iso,
                                                  u_aniso    = u_aniso,
                                                  occupancy  = occupancy,
                                                  fp         = fp,
                                                  fdp        = fdp,
                                                  tan_u_iso  = tan_u_iso,
                                                  param      = param)
                    for scatterer in scatterers:
                        assert scatterer.flags.grad_site()      == site
                        if(scatterer.flags.use_u_iso()):
                           assert scatterer.flags.grad_u_iso()     == u_iso
                           assert scatterer.flags.grad_u_aniso()   == False
                           cn1+=1
                        if(scatterer.flags.use_u_aniso()):
                           assert scatterer.flags.grad_u_iso()     == False
                           assert scatterer.flags.grad_u_aniso()   == u_aniso
                           cn2+=1
                        assert scatterer.flags.grad_occupancy() == occupancy
                        assert scatterer.flags.grad_fp()        == fp
                        assert scatterer.flags.grad_fdp()       == fdp
                        assert scatterer.flags.param            == param
  assert cn1 != 0
  assert cn2 != 0
  xray.set_scatterer_grad_flags(scatterers = scatterers)
  for scatterer in scatterers:
      assert scatterer.flags.grad_site()      == False
      assert scatterer.flags.grad_u_iso()     == False
      assert scatterer.flags.grad_u_aniso()   == False
      assert scatterer.flags.grad_occupancy() == False
      assert scatterer.flags.grad_fp()        == False
      assert scatterer.flags.grad_fdp()       == False
      assert scatterer.flags.tan_u_iso()      == False
      assert scatterer.flags.param            == 0
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                fp         = True,
                                param      = -1)
  for scatterer in scatterers:
      assert scatterer.flags.grad_site()      == False
      assert scatterer.flags.grad_u_iso()     == False
      assert scatterer.flags.grad_u_aniso()   == False
      assert scatterer.flags.grad_occupancy() == False
      assert scatterer.flags.grad_fp()        == True
      assert scatterer.flags.grad_fdp()       == False
      assert scatterer.flags.tan_u_iso()      == False
      assert scatterer.flags.param            == -1
  sc0 = scatterers[0]
  for iso in [False, True]:
      for aniso in [True, False]:
          sc0.flags.set_use_u_iso(iso)
          sc0.flags.set_use_u_aniso(aniso)
          assert sc0.flags.use_u_iso() == iso
          assert sc0.flags.use_u_aniso() == aniso
          sc0.flags.set_use_u(iso = aniso, aniso = iso)
          assert sc0.flags.use_u_iso() == aniso
          assert sc0.flags.use_u_aniso() == iso
          if (iso): sc0.flags.set_use_u_iso_only()
          else:     sc0.flags.set_use_u_aniso_only()
          assert sc0.flags.use_u_iso() == iso
          assert sc0.flags.use_u_aniso() != iso
  for iso in [False, True]:
      for aniso in [True, False]:
          sc0.flags.set_use_u_iso(iso)
          sc0.flags.set_use_u_aniso(aniso)
          assert sc0.flags.use_u_iso() == iso
          assert sc0.flags.use_u_aniso() == aniso
          sc0.u_iso = 0.5
          sc0.u_star = [0.5,0.5,0.5,0.5,0.5,0.5]
          sc0.set_use_u(iso = aniso, aniso = iso)
          assert sc0.flags.use_u_iso() == aniso
          assert sc0.flags.use_u_aniso() == iso
          if(aniso == False): assert sc0.u_iso == -1.0
          if(iso == False): assert sc0.u_star == (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
          if(aniso): assert sc0.u_iso == 0.5
          if(iso): assert sc0.u_star == (0.5,0.5,0.5,0.5,0.5,0.5)
          sc0.u_iso = 0.5
          sc0.u_star = [0.5,0.5,0.5,0.5,0.5,0.5]
          if (iso): sc0.set_use_u_iso_only()
          else:     sc0.set_use_u_aniso_only()
          assert sc0.flags.use_u_iso() == iso
          assert sc0.flags.use_u_aniso() != iso
          if(iso):
             assert sc0.u_iso == 0.5
             assert sc0.u_star == (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
          else:
             assert sc0.u_iso == -1.0
             assert sc0.u_star == (0.5,0.5,0.5,0.5,0.5,0.5)

def exercise_set_selected_scatterer_grad_flags():
  x = xray.scatterer("c", site=(0.1,0.2,0.3), occupancy=0.0, u=(0,0,0,0,0,0))
  y = xray.scatterer("c", site=(0.1,0.2,0.3), occupancy=0.0, u=1.0)
  scatterers = flex.xray_scatterer(50, x)
  scatterers.extend(flex.xray_scatterer(50, y))
  mt = flex.mersenne_twister(seed=0)
  scatterers = scatterers.select(mt.random_permutation(size=100))
  def all_clear():
    for scatterer in scatterers:
      f = scatterer.flags
      if (f.grad_site()): return False
      if (f.grad_u_iso()): return False
      if (f.grad_u_aniso()): return False
      if (f.grad_occupancy()): return False
      if (f.grad_fp()): return False
      if (f.grad_fdp()): return False
    return True
  assert all_clear()
  xray.set_scatterer_grad_flags(scatterers=scatterers)
  assert all_clear()
  sels = {}
  for attr in "site u_iso u_aniso occupancy fp fdp".split():
    sels[attr] = mt.random_bool(size=100, threshold=0.5)
  for i_seq, scatterer in enumerate(scatterers):
    if (not scatterer.flags.use_u_iso()): sels["u_iso"][i_seq] = False
    if (not scatterer.flags.use_u_aniso()): sels["u_aniso"][i_seq] = False
  scatterers.flags_set_grads(state=False)
  assert all_clear()
  scatterers.flags_set_grads(state=True)
  assert not all_clear()
  for attr,sel in sels.items():
    scatterers.flags_set_grads(state=False)
    assert all_clear()
    getattr(scatterers, "flags_set_grad_"+attr)(iselection=sel.iselection())
    assert not all_clear()
    for scatterer,flag in zip(scatterers, sel):
      assert getattr(scatterer.flags, "grad_"+attr)() == flag

def exercise_scatterer_flags_counts():
  x = xray.scatterer("c", site=(0.1,0.2,0.3), occupancy=0.0, u=(0,0,0,0,0,0))
  scatterers = flex.xray_scatterer(10, x)
  manager = xray.scatterer_grad_flags_counts(scatterers)
  assert manager.site      == 0
  assert manager.u_iso     == 0
  assert manager.u_aniso   == 0
  assert manager.occupancy == 0
  assert manager.fp        == 0
  assert manager.fdp       == 0
  scatterers[1].flags.set_grad_site      (True)
  scatterers[1].flags.set_grad_u_iso     (True)
  scatterers[1].flags.set_grad_u_aniso   (True)
  scatterers[1].flags.set_grad_occupancy (True)
  scatterers[1].flags.set_grad_fp        (True)
  scatterers[1].flags.set_grad_fdp       (True)
  manager = xray.scatterer_grad_flags_counts(scatterers)
  #assert manager.site      == 1     #XXX
  #assert manager.u_iso     == 1     #XXX
  #assert manager.u_aniso   == 1     #XXX
  assert manager.occupancy == 1
  assert manager.fp        == 1
  assert manager.fdp       == 1

def exercise_conversions():
  d = flex.double((10,-1))
  s = flex.double((1,2))
  r = xray.array_f_sq_as_f_xtal_3_7(d, s)
  r = xray.array_f_sq_as_f_xtal_3_7(d, s, 1.e-6)
  assert approx_equal(r.f, (3.1622777, 0))
  assert approx_equal(r.sigma_f, (0.1543471, 1.4142136))
  r = xray.array_f_sq_as_f_xtal_3_7(d)
  assert approx_equal(r.f, (3.1622777, 0))
  assert r.sigma_f.size() == 0
  r = xray.array_f_sq_as_f_crystals(d, s)
  r = xray.array_f_sq_as_f_crystals(d, s, 1.e-6)
  assert approx_equal(r.f, (3.1622777, -1))
  assert approx_equal(r.sigma_f, (0.1581139, 2))
  r = xray.array_f_sq_as_f_crystals(d)
  assert approx_equal(r.f, (3.1622777, -1))
  assert r.sigma_f.size() == 0
  r = xray.array_f_as_f_sq(d, s)
  assert approx_equal(r.f_sq, (100, 1))
  assert approx_equal(r.sigma_f_sq, (20, -4))
  r = xray.array_f_as_f_sq(d)
  assert approx_equal(r.f_sq, (100, 1))
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
  c = xray.ext.scatterer(other=x)
  assert c.fp == 1
  c.fp = 3
  assert x.fp == 1
  assert c.fp == 3
  f = xray.scatterer_flags()
  f.set_use_u_iso_only()
  c.flags = f
  assert c.flags.use_u_iso_only()
  c.flags.set_use_u_aniso_only()
  assert f.use_u_iso_only()
  assert c.flags.use_u_aniso_only()
  for flag0 in [True, False]:
    x.flags.set_use(flag0)
    for flag1 in [True, False]:
      x.flags.set_use_u_iso(flag1)
      for flag2 in [True, False]:
        x.flags.set_use_u_aniso(flag2)
        for flag3 in [True, False]:
          x.flags.set_grad_u_iso(flag3)
          for flag4 in [True, False]:
            x.flags.set_grad_u_aniso(flag4)
            assert x.flags.use()          == flag0
            assert x.flags.use_u_iso()    == flag1
            assert x.flags.use_u_aniso()  == flag2
            assert x.flags.grad_u_iso()   == flag3
            assert x.flags.grad_u_aniso() == flag4
  assert approx_equal(x.u_iso, 0.25)
  x.u_iso = 0.52
  assert approx_equal(x.u_iso, 0.52)
  x = xray.scatterer("a", (0.1,0.2,0.3), (1,2,3,4,5,6), 0.9, "const", 0, 0)
  assert x.flags.use_u_aniso()
  assert approx_equal(x.u_star, (1,2,3,4,5,6))
  x.u_star = (3,2,1,6,5,4)
  assert approx_equal(x.u_star, (3,2,1,6,5,4))
  assert x.flags.use() == True
  x.flags.set_grad_site(state=True)
  assert x.flags.grad_site()
  #
  x.u_iso = 4
  x.set_use_u(iso=True, aniso=True)
  assert approx_equal(x.u_iso, 4)
  assert approx_equal(x.u_star, (3,2,1,6,5,4))
  x.set_use_u_iso_only()
  assert approx_equal(x.u_iso, 4)
  assert approx_equal(x.u_star, (-1,-1,-1,-1,-1,-1))
  x.u_star = (9,3,6,0,1,5)
  x.set_use_u_aniso_only()
  assert approx_equal(x.u_iso, -1)
  assert approx_equal(x.u_star, (9,3,6,0,1,5))
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
  x.flags.set_use_u_aniso(True)
  try:
    x.apply_symmetry(uc, sg.group(), u_star_tolerance=0.1)
  except RuntimeError, e:
    assert str(e).find("is_compatible_u_star") > 0
  else:
    raise Exception_expected
  x.apply_symmetry(site_symmetry_ops=ss)
  x.apply_symmetry(site_symmetry_ops=ss, u_star_tolerance=0.5)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0)
  ss = x.apply_symmetry(uc, sg.group(), 0.5, 0, False)
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
  assert x.flags.use_u_aniso_only()
  x.convert_to_isotropic(unit_cell=uc)
  assert not x.flags.use_u_aniso()
  assert approx_equal(x.u_iso, 269)
  assert approx_equal(x.u_star, (-1,-1,-1,-1,-1,-1))
  x.convert_to_anisotropic(unit_cell=uc)
  assert x.flags.use_u_aniso()
  assert approx_equal(x.u_iso, -1)
  assert approx_equal(x.u_star, (2.69, 2.69, 1.59171598, 0, 0, 0))
  x.u_star = (1,2,3,4,5,6)
  assert not x.is_positive_definite_u(unit_cell=uc)
  assert not x.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=1.e2)
  assert x.is_positive_definite_u(unit_cell=uc, u_cart_tolerance=1.e3)
  x.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=0, u_max=9999,
    anisotropy_min=0)
  assert approx_equal(x.u_star,
    (3.3379643647809192, 4.5640522609325131, 4.4690204772593507,
     3.9031581835726965, 3.8623090371651934, 4.5162864184404032))
  x.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=1, u_max=9999,
    anisotropy_min=0)
  assert approx_equal(x.u_star,
    (3.3458045216665266, 4.5710990727698393, 4.4720459395534728,
     3.9006326295505751, 3.8598099147456764, 4.5133641373560351))
  assert x.is_positive_definite_u(unit_cell=uc)
  f = xray.scatterer_flags()
  f.set_use_u_iso_only()
  y = x.customized_copy(u=-1, flags=f)
  assert not y.flags.use_u_aniso()
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
  y.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=1, u_max=1.0,
    anisotropy_min=0)
  assert approx_equal(y.u_iso, 1.0)
  yy = xray.scatterer("si1", site=(0.01,0.02,0.3), occupancy=0.8, u=158.0)
  assert approx_equal(yy.u_iso, 158.0)
  yy.tidy_u(unit_cell=uc, site_symmetry_ops=ss, u_min=1, u_max=0.9,
    anisotropy_min=0)
  assert approx_equal(yy.u_iso, 0.9)
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
  #
  f = xray.scatterer_flags()
  f.set_use_u(iso=True, aniso=True)
  y = x.customized_copy(u_iso=0.3, u_star=(1,3,5,-2,6,4), flags=f)
  assert y.flags.use_u_iso()
  assert y.flags.use_u_aniso()
  assert approx_equal(y.u_iso, 0.3)
  assert approx_equal(y.u_star, (1,3,5,-2,6,4))

  scs = flex.xray_scatterer((
    xray.scatterer("o", site=(0.01,0.02,0.1), u=0.1),
    xray.scatterer("o", site=(0.02,0.02,0.2), u=0.2),
    xray.scatterer("o", site=(0.03,0.02,0.3), u=0.3),
    xray.scatterer("o", site=(0.04,0.02,0.4), u=0.4),
    xray.scatterer("o", site=(0.05,0.02,0.4), u=0.5),
    xray.scatterer("o", site=(0.1,0.2,0.3),   u=(0.1,0.2,0.3,-.04,0.5,-0.06)),
    xray.scatterer("o", site=(0.3,0.4,0.5),   u=(0.4,0.5,0.6,-.05,0.2,-0.02))))
  uisos = scs.extract_u_iso()
  xray.shift_us(scatterers=scs, unit_cell=uc, u_shift=1,
                selection = flex.size_t([1,3,5]))
  xray.shift_occupancies(scatterers=scs, q_shift=2.5,
                         selection = flex.size_t([1,3,5]))
  assert approx_equal(scs.extract_occupancies(), [1.0, 3.5, 1.0, 3.5, 1.0, 3.5, 1.0])
  xray.shift_occupancies(scatterers=scs, q_shift=1.0)
  assert approx_equal(scs.extract_occupancies(), [2.0, 4.5, 2.0, 4.5, 2.0, 4.5, 2.0])
  assert approx_equal(scs[0].u_iso,  0.1) and approx_equal(scs[0].u_star, (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  assert approx_equal(scs[1].u_iso,  1.2) and approx_equal(scs[1].u_star, (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  assert approx_equal(scs[2].u_iso,  0.3) and approx_equal(scs[2].u_star, (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  assert approx_equal(scs[3].u_iso,  1.4) and approx_equal(scs[3].u_star, (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  assert approx_equal(scs[4].u_iso,  0.5) and approx_equal(scs[4].u_star, (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  assert approx_equal(scs[5].u_iso, -1.0)
  assert approx_equal(scs[6].u_iso, -1.0) and approx_equal(scs[6].u_star, (0.40, 0.50, 0.59999999999999998, -0.05, 0.2, -0.02))
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
  #
  pc = a[0].as_py_code(indent="%&")
  assert not show_diff(pc, """\
xray.scatterer(
%&  label="si1",
%&  site=(0.200000, 0.500000, 0.400000),
%&  u=(3.345805, 4.571099, 4.472046,
%&     3.900633, 3.859810, 4.513364),
%&  occupancy=0.800000,
%&  fp=3.000000,
%&  fdp=0.000000)""")
  pc = a[1].as_py_code(indent="Q")
  assert not show_diff(pc, """\
xray.scatterer(
Q  label="si1",
Q  site=(0.200000, 0.500000, 0.400000),
Q  u=1.000000,
Q  occupancy=0.800000,
Q  fp=0.000000,
Q  fdp=4.000000)""")
  for xs in a:
    xs2 = eval(xs.as_py_code())
    assert approx_equal(xs2.site, xs.site)
    assert approx_equal(xs2.fp, xs.fp)
    assert approx_equal(xs2.fdp, xs.fdp)
  #
  uc = uctbx.unit_cell((1,2,3, 89, 96, 107))
  results_u_cart_plus_u_iso = []
  for u_iso, u_aniso in [(0, 0), (0.04, 0), (0, 0.04), (0.04, 0.04)]:
    xs = xray.scatterer(label="C")
    xs.flags.set_use_u_iso(False)
    xs.flags.set_use_u_aniso(False)
    if (u_iso != 0):
      xs.u_iso = u_iso
      xs.flags.set_use_u_iso(True)
    if (u_aniso != 0):
      xs.u_star = adptbx.u_cart_as_u_star(
        uc, (u_aniso/2, u_aniso, 2*u_aniso, 0, 0, 0))
      xs.flags.set_use_u_aniso(True)
    expected = u_iso + 3.5*u_aniso/3
    assert approx_equal(xs.u_iso_or_equiv(unit_cell=uc), expected)
    if (not xs.flags.use_u_aniso()):
      assert approx_equal(xs.u_iso_or_equiv(unit_cell=None), expected)
    else:
      try: xs.u_iso_or_equiv(unit_cell=None)
      except RuntimeError, e:
        assert str(e).startswith("cctbx Internal Error: ")
      else: raise Exception_expected
    results_u_cart_plus_u_iso.append(xs.u_cart_plus_u_iso(unit_cell=uc))
    if (not xs.flags.use_u_aniso()):
      assert approx_equal(
        xs.u_cart_plus_u_iso(unit_cell=None), (u_iso,u_iso,u_iso,0,0,0))
    else:
      try: xs.u_cart_plus_u_iso(unit_cell=None)
      except RuntimeError, e:
        assert str(e).startswith("cctbx Internal Error: ")
      else: raise Exception_expected
  assert approx_equal(results_u_cart_plus_u_iso, [
    (0,0,0,0,0,0),
    (0.04,0.04,0.04,0,0,0),
    (0.02,0.04,0.08,0,0,0),
    (0.06,0.08,0.12,0,0,0)])

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
  assert not reg.has_key("const")
  assert not reg.has_key("unknown")
  assert reg.process(scattering_type="const") == 0
  assert reg.size() == 1
  assert reg.has_key("const")
  assert not reg.has_key("unknown")
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
  assert approx_equal(reg.occupancy_sums(scatterers=scatterers), [1,2,3,2,1])
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
    all_dilated_ffs = reg.dilated_form_factors_at_d_star_sq(
      d_star_sq=0.1,
      dilation_coeffs=flex.double_range(1,scatterers.size()+1),
      unique_indices=unique_indices)
    for i, dilated_ff in enumerate(all_dilated_ffs):
      j = unique_indices[i]
      ref_ff = reg.unique_gaussians_as_list()[j].at_d_star_sq(0.1/(i+1))
      assert approx_equal(dilated_ff, ref_ff)
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
    for show_sf0 in [False, True]:
      for show_gaussians in [False, True]:
        reg.show(
          show_sf0=show_sf0,
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
:#  Type    Number    sf(0)
:#   Al         2     20.00
:#   Si         2     14.00
:#   const      1     10.00
:#   O          3      8.00
:#   custom     1      None
:#  sf(0) = scattering factor at diffraction angle 0.
:#Number of scattering types: 5
:#  Type    Number    sf(0)   Gaussians
:#   Al         2     20.00       0+c
:#   Si         2     14.00       %(n_terms)d+c
:#   const      1     10.00       0+c
:#   O          3      8.00       %(n_terms)d+c
:#   custom     1      None      None
:#  sf(0) = scattering factor at diffraction angle 0.
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
  else: raise Exception_expected
  try: reg.gaussian_not_optional(scattering_type="custom")
  except RuntimeError, e:
    assert str(e) == 'gaussian not defined for scattering_type "custom".'
  else: raise Exception_expected
  try: reg.unique_form_factors_at_d_star_sq(d_star_sq=0)
  except RuntimeError, e:
    assert str(e) == 'gaussian not defined for scattering_type "custom".'
  else: raise Exception_expected

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
    assert_min_distance_sym_equiv=True,
    non_unit_occupancy_implies_min_distance_sym_equiv_zero=False)
  assert list(site_symmetry_table.special_position_indices()) == [0]
  xray.tidy_us(
    scatterers=scatterers,
    unit_cell=uc,
    site_symmetry_table=site_symmetry_table,
    u_min=0,
    u_max = 99999,
    anisotropy_min=0)
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
  assert approx_equal(flex.arg(fc, True), (-36, 72))
  assert approx_equal(flex.abs(fc), (10.50871, 9.049631))
  assert approx_equal(flex.arg(fc, True), (-36, 72))
  fc = xray.ext.structure_factors_direct(
    uc, sg.group(), mi, scatterers, scattering_type_registry).f_calc()
  assert approx_equal(flex.abs(fc), (10.50871, 9.049631))
  assert approx_equal(flex.arg(fc, True), (-36, 72))
  xray.tidy_us(
    scatterers=scatterers,
    unit_cell=uc,
    site_symmetry_table=site_symmetry_table,
    u_min=100,
    u_max = 99999,
    anisotropy_min=0)
  assert approx_equal(scatterers[0].u_iso, 100)
  assert approx_equal(scatterers[1].u_star,
    (1.0134539945616343, 1.0005190241807682, 0.64980451464405997,
     -0.0026425269166861672, 0.027955730692513142, -0.0054908429234285239))
  xray.ext.structure_factors_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, scattering_type_registry).f_calc()
  gradient_flags = xray.ext.gradient_flags(
                              True, True, True, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray.ext.structure_factors_gradients_direct(
    uc, sg.group(), mi, scatterers, None,
    scattering_type_registry, site_symmetry_table,
    flex.complex_double(mi.size()),
    0)
  gradient_flags = xray.ext.gradient_flags(
                              True, True, True, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray.ext.structure_factors_gradients_direct(
    math_module.cos_sin_table(12),
    uc, sg.group(), mi, scatterers, None,
    scattering_type_registry, site_symmetry_table,
    flex.complex_double(mi.size()),
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
    (0.0013784963+0.002756992j), (0.0103982354+0.013864313j),
    (0.0160141831+0.032028366j), (0.0004572786-0.000365822j),
    (0.0014117387-0.001694086j)))

  import cmath
  f_obs_sqr = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double([ cmath.sqrt(x) for x in f_obs_sqr ])
  ls = xray.targets_least_squares_residual_for_intensity(f_obs_sqr, w, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, w, f_calc, True)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, w, f_calc, False, 3)
  assert approx_equal(ls.scale_factor(), 3)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, f_calc)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target(), 0)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double([ cmath.sqrt(x) for x in (10,20,30,40,50) ])
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, f_calc, True)
  assert approx_equal(ls.scale_factor(), 1/10.)
  assert approx_equal(ls.target(), 0)
  assert approx_equal(tuple(ls.derivatives()), (0j,0j,0j,0j,0j))
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, f_calc, False, 3/10.)
  assert approx_equal(ls.scale_factor(), 3/10.)
  assert approx_equal(ls.target(), 4)
  assert ls.derivatives().size() == 0
  f_calc = flex.complex_double([ cmath.sqrt(x)
                                 for x in (1+2j,3+4j,-1-2j,5-4j,-5+6j) ])
  w = flex.double((1,2,3,2,4))
  ls = xray.targets_least_squares_residual_for_intensity(
    f_obs_sqr, w, f_calc, True)
  assert approx_equal(ls.scale_factor(), 0.6307845)
  assert approx_equal(ls.target(), 0.06211855)
  assert approx_equal(tuple(ls.derivatives()), (
    (0.00784178+0.00484648j), (0.0693216+0.0346608j),
    (-0.0563023+0.091099j), (0.0027966-0.000980993j),
    (-0.00522801-0.011162j)))

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
  assert approx_equal(d.wing_cutoff(), 1.e-3)
  assert approx_equal(d.exp_table_one_over_step_size(), -100)
  assert approx_equal(d.tolerance_positive_definite(), 1.e-5)
  assert d.n_scatterers_passed() == 2
  assert d.n_contributing_scatterers() == 2
  assert d.n_anomalous_scatterers() == 1
  assert d.anomalous_flag()
  assert d.real_map().size() == 0
  assert d.complex_map().size() == (20*20*22)
  assert d.exp_table_size() == 2889
  assert d.max_sampling_box_n_points() == 216
  assert d.sum_sampling_box_n_points() == 341
  assert approx_equal(d.ave_sampling_box_n_points(), 341/2.)
  assert d.max_sampling_box_edges() == (6,6,6)
  assert approx_equal(d.max_sampling_box_edges_frac(), (6/20.,6/20.,6/22.))
  assert d.excessive_sampling_radius_i_seqs().size() == 0
  assert d.grid_indices_for_each_scatterer().size() == 0
  i = flex.miller_index(((1,2,3), (2,3,4)))
  f = flex.complex_double((1+2j, 2+3j))
  d.eliminate_u_extra_and_normalize(i, f)
  f_orig = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), 0.2, i, f)
  f_elim = f.deep_copy()
  xray.apply_u_extra(d.unit_cell(), -0.2, i, f, 1)
  assert approx_equal(f, f_orig)
  #
  scatterers[0].fdp = 0
  for sgifes in [1, -1]:
    d = xray.sampled_model_density(
      unit_cell=uc,
      scatterers=scatterers,
      scattering_type_registry=scattering_type_registry,
      fft_n_real=(20,20,22),
      fft_m_real=(20,20,23),
      store_grid_indices_for_each_scatterer=sgifes)
    if (sgifes > 0):
      assert d.real_map().size() == 20*20*23
      assert d.real_map().focus() == (20,20,22)
      assert d.real_map().all() == (20,20,23)
    else:
      assert d.real_map().size() == 0
    assert d.complex_map().size() == 0
    assert d.grid_indices_for_each_scatterer().size() == scatterers.size()
    map = d.real_map()
    n_non_zero = 0
    for gi,expected_sizes in zip(d.grid_indices_for_each_scatterer(),
                                 [[88], [33,31]]):
      assert gi.size() in expected_sizes
      if (sgifes < 0): continue
      for i_map in gi:
        assert map[i_map] > 0
      n_non_zero += gi.size()
    if (sgifes > 0):
      assert map.count(0) + n_non_zero == 20*20*23
    comb_sel = shared.stl_set_unsigned()
    comb_sel.append_union_of_selected_arrays(
      arrays=d.grid_indices_for_each_scatterer(),
      selection=flex.size_t([0,1]))
    assert comb_sel[0].size() in [121,119]
    if (sgifes > 0):
      assert map.select(comb_sel[0]).all_gt(0)
      assert comb_sel[0].size() == n_non_zero

def exercise_minimization_apply_shifts():
  uc = uctbx.unit_cell((20, 20, 23))
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  f = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(
    scatterers=scatterers,
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
    uc, scatterers, shifts).shifted_scatterers
  for a,b in zip(scatterers, shifted_scatterers):
    assert a.scattering_type == b.scattering_type
    assert a.site != b.site
    if (not a.flags.use_u_aniso()):
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
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 6, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [1]*6
  shifts = flex.double((-1,2,-3,4,-5,-6))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  assert approx_equal(
    shifted_scatterers[0].site,
    (0.01-1/20.,0.02+2/20.,0.3-3/23.))
  assert approx_equal(
    shifted_scatterers[1].site,
    (0.3+4/20.,0.4-5/20.,0.5-6/23.))
  f = xray.ext.gradient_flags(
    False, True, False, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 1, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [2]
  shifts = flex.double(1, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  assert approx_equal(shifted_scatterers[0].u_iso, -10)
  f = xray.ext.gradient_flags(
    False, False, True, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 6, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [3]*6
  shifts = flex.double(6, -100)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  assert not_approx_equal(shifted_scatterers[1].u_star,
    [u_ij-100 for u_ij in scatterers[1].u_star])
  f = xray.ext.gradient_flags(
    False, False, False, True, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [4]*2
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  for i in xrange(2):
    assert approx_equal(shifted_scatterers[i].occupancy, -9)
  f = xray.ext.gradient_flags(
    False, False, False, False, True, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [5]*2
  shifts = flex.double(2, -10)
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  assert approx_equal(shifted_scatterers[0].fp, -11)
  assert shifted_scatterers[1].fp == -10
  for i in xrange(2):
    assert shifted_scatterers[i].fdp == scatterers[i].fdp
  f = xray.ext.gradient_flags(
    False, False, False, False, False, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = f.site,
                                u_iso      = f.u_iso,
                                u_aniso    = f.u_aniso,
                                occupancy  = f.occupancy,
                                fp         = f.fp,
                                fdp        = f.fdp)
  s = xray.ext.minimization_shift_scales(scatterers, 2, 1,2,3,4,5,6)
  assert [int(x) for x in s] == [6]*2
  shifts = flex.double((2,3))
  shifted_scatterers = xray.ext.minimization_apply_shifts(
    uc, scatterers, shifts).shifted_scatterers
  assert shifted_scatterers[0].fp == -1
  assert approx_equal(shifted_scatterers[0].fdp, 4)
  assert shifted_scatterers[1].fp == 0
  assert shifted_scatterers[1].fdp == 3
  shifts = flex.double(1)
  try:
    xray.ext.minimization_apply_shifts(uc, scatterers, shifts)
  except Exception, e:
    assert str(e) == "scitbx Error: Array of shifts is too small."
  else:
    raise Exception_expected
  shifts = flex.double(3)
  try:
    xray.ext.minimization_apply_shifts(uc, scatterers, shifts)
  except Exception, e:
    assert str(e) == "cctbx Error: Array of shifts is too large."
  else:
    raise Exception_expected

def exercise_minimization_add_gradients():
  uc = uctbx.unit_cell((20, 20, 23))
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", site=(0.01,0.02,0.3), fp=-1, fdp=2),
    xray.scatterer("O1", site=(0.3,0.4,0.5),
                   u=adptbx.u_cart_as_u_star(uc,
                     (0.04,0.05,0.06,-.005,0.02,-0.002)))))
  gradient_flags = xray.ext.gradient_flags(
    True, False, False, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(6))
  geometry_restraints_site_gradients = flex.vec3_double([(1,-2,3),(-4,-5,6)])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,-1,-1,11])
  gradient_flags = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(19))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,5,6,3,3,15,10,11,12,13,14,15,16,17,18])
  gradient_flags = xray.ext.gradient_flags(
    True, True, False, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(13))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,5,6,3,3,15,10,11,12])
  gradient_flags = xray.ext.gradient_flags(
    True, False, True, True, False, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(16))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=geometry_restraints_site_gradients,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,3,4,1,1,13,8,9,10,11,12,13,14,15])
  site_gradients = xray.ext.minimization_extract_site_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients)
  assert approx_equal(site_gradients, [(1,-1,5), (1,1,13)])
  #
  gradient_flags = xray.ext.gradient_flags(
    False, True, False, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double([1])
  u_iso_gradients = flex.double([3,0])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients, [4])
  gradient_flags = xray.ext.gradient_flags(
    True, True, True, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(19))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [0,1,2,6,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
  gradient_flags = xray.ext.gradient_flags(
    True, True, False, True, True, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(13))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [0,1,2,6,4,5,6,7,8,9,10,11,12])
  gradient_flags = xray.ext.gradient_flags(
    False, True, True, True, False, True, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(11))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=u_iso_gradients,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [3,1,2,3,4,5,6,7,8,9,10])
  #
  gradient_flags = xray.ext.gradient_flags(
    False, False, False, True, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double([1,2])
  occupancy_gradients = flex.double([3,4])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=occupancy_gradients)
  assert approx_equal(xray_gradients, [4,6])
  xray_gradients = flex.double([2,1])
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=None,
    u_aniso_gradients=None,
    occupancy_gradients=None)
  assert approx_equal(xray_gradients, [2,1])
  #
  gradient_flags = xray.ext.gradient_flags(
    False, False, True, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(6))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=None,
    u_iso_gradients=None,
    u_aniso_gradients=flex.sym_mat3_double([(0,0,0,0,0,0),(1,-1,-4,-7,12,-3)]),
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1, 0, -2, -4, 16, 2])
  #
  gradient_flags = xray.ext.gradient_flags(
    True, False, True, False, False, False, False, False)
  xray.set_scatterer_grad_flags(scatterers = scatterers,
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray_gradients = flex.double(xrange(12))
  xray.ext.minimization_add_gradients(
    scatterers=scatterers,
    xray_gradients=xray_gradients,
    site_gradients=flex.vec3_double([(1,-2,3),(-4,-5,6)]),
    u_iso_gradients=None,
    u_aniso_gradients=flex.sym_mat3_double([(0,0,0,0,0,0),(1,-1,-4,-7,12,-3)]),
    occupancy_gradients=None)
  assert approx_equal(xray_gradients,
    [1,-1,5,-1,-1,11,7, 6, 4, 2, 22, 8])

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

def exercise_targets_common_results():
  r = xray.targets_common_results(
    target_per_reflection=flex.double([1,2,3]),
    target_work=13,
    target_test=None,
    gradients_work=flex.complex_double([1+2j,2+3j,3-4j]))
  assert approx_equal(r.target_per_reflection(), [1,2,3])
  assert approx_equal(r.target_work(), 13)
  assert r.target_test() is None
  assert approx_equal(r.gradients_work(), [1+2j,2+3j,3-4j])
  r = xray.targets_common_results(
    target_per_reflection=flex.double(),
    target_work=42,
    target_test=53,
    gradients_work=flex.complex_double())
  assert r.target_per_reflection().size() == 0
  assert approx_equal(r.target_work(), 42)
  assert approx_equal(r.target_test(), 53)
  assert r.gradients_work().size() == 0

def exercise_targets_least_squares():
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=True,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=0,
    scale_factor=0)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target_per_reflection(), [0]*5)
  assert approx_equal(ls.target_work(), 0)
  assert ls.target_test() is None
  assert ls.gradients_work().size() == 0
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=True,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=0,
    scale_factor=2.0)
  assert approx_equal(ls.scale_factor(), 2.0)
  assert approx_equal(ls.target_per_reflection(), [1,4,9,16,25])
  assert approx_equal(ls.target_work(),1.0)
  assert ls.gradients_work().size() == 0
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((4,5,6))
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=True,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=1,
    scale_factor=0)
  assert approx_equal(ls.scale_factor(), 50./134.)
  assert approx_equal(ls.target_work(), 0.0671641791)
  assert approx_equal(ls.gradients_work(),
                     ((0.0551347738+0j),(-0.0100245043+0j),(-0.0284027623+0j)) )
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=True,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=1,
    scale_factor=2.0)
  assert approx_equal(ls.scale_factor(), 2.0)
  assert approx_equal(ls.target_work(),17.8)
  assert approx_equal(ls.gradients_work(), ((4.2+0j),(3.2+0j),(1.8+0j)) )
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j))
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=True,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=1,
    scale_factor=0)
  assert approx_equal(ls.scale_factor(), 0.4773772552)
  assert approx_equal(ls.target_work(), 0.2023883467)
  assert approx_equal(ls.gradients_work(),
                       ((0.0043198335244903152+0.0086396670489806305j),
                        (0.022162885026120613+0.029550513368160818j),
                        (0.041257975234691303+0.082515950469382607j)) )
  #
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  f_calc = flex.complex_double((1,2,3,4,5))
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((4,5,6))
  #
  f_obs = flex.double((1,2,3))
  w = flex.double((3,2,1))
  f_calc = flex.complex_double((1+2j,3+4j,-1-2j))
  #
  f_obs = flex.double((1,2,3,4,5))
  w = flex.double((1,1,1,1,1))
  r_free_flags = flex.bool([False,True,False,True,False])
  f_calc = flex.complex_double((2,3,4,5,6))
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=False,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=r_free_flags,
    f_calc=f_calc,
    derivatives_depth=1,
    scale_factor=0)
  assert not ls.compute_scale_using_all_data()
  assert approx_equal(ls.scale_factor(), 0.785714285714)
  assert approx_equal(ls.target_work(), 0.0122448979592)
  assert approx_equal(ls.target_test(), 0.00663265306122)
  assert approx_equal(ls.gradients_work(), [
    0.025655976676384838+0j,
    0.0064139941690962068+0j,
    -0.012827988338192414-0j])
  #
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=False,
    obs_type="F",
    obs=f_obs,
    weights=w,
    r_free_flags=r_free_flags,
    f_calc=flex.complex_double(5), # all f_calc 0
    derivatives_depth=1,
    scale_factor=1)
  assert approx_equal(ls.scale_factor(), 1)
  assert approx_equal(ls.target_work(), 1)
  assert approx_equal(ls.target_test(), 1)
  assert approx_equal(ls.gradients_work(), [0j,0j,0j])

def exercise_maximum_likelihood_targets():
  f_obs = flex.double([
    35.6965, 308.042, 128.949, 35.4385, 29.259, 108.11, 35.3133, 67.2968,
    109.585, 130.959, 78.7692, 138.602])
  r_free_flags = flex.bool([
    False, False, True, True, False, False, False, False, False, False,
    False, False])
  experimental_phases = flex.hendrickson_lattman([
    [0.457488,0,0,0], [-0.769364,0,0,0], [1.19222,0,0,0], [-0.0931637,0,0,0],
    [0.890576,0,0,0], [-1.66137,1.69479,-0.107808,1.41357],
    [-0.00674594,0.847689,-1.4854,0.547142],
    [-0.0295139,-0.00293545,0.130557,0.0780442], [-0.409673,0,0,0],
    [-0.416943,0,0,0], [0.101886,2.88433,2.87327,2.71556], [0.901622,0,0,0]])
  f_calc = flex.complex_double([
    -19.5412, -171.232, 69.435, -32.7787, -41.0725, 32.8168+39.3201j,
    36.8631-10.3995j, -7.23216-79.5072j, 71.4597, 69.1851,
    41.9767-14.4786j, 73.0327])
  alpha = flex.double([
    0.211893, 0.368057, 0.71595, 0.6732, 0.858985, 0.987573, 0.973331,
    0.161511, 0.439952, 0.530318, 0.272504, 0.95601])
  beta = flex.double([
    0.963688, 0.243077, 0.755248, 0.98828, 0.668677, 0.0767036, 0.745345,
    0.55896, 0.415472, 0.0690112, 0.563668, 0.329134])
  epsilons = flex.int([2, 4, 1, 8, 1, 1, 2, 1, 4, 1, 8, 1])
  centric_flags = flex.bool([
    True, True, True, True, True, False, False, False, True, True,
    False, True])
  for rff in [None, r_free_flags]:
    for compute_gradients in [False, True]:
      ml = xray.targets_maximum_likelihood_criterion(
        f_obs=f_obs,
        r_free_flags=rff,
        f_calc=f_calc,
        alpha=alpha,
        beta=beta,
        scale_factor=1,
        epsilons=epsilons,
        centric_flags=centric_flags,
        compute_gradients=compute_gradients)
      tpr = ml.target_per_reflection()
      tw = ml.target_work()
      tt = ml.target_test()
      gw = ml.gradients_work()
      assert approx_equal(tpr, [
        259.57027227308151, 30872.934956654299, 4157.3629318162803,
        13.260741195206734, 27.831255194364076, 43149.780572387972,
        3.3949285154155859, 5294.345206664615, 1838.4876581905569,
        64384.962967657091, 986.06706016214844, 7187.3520116172294])
      if (not compute_gradients):
        assert gw.size() == 0
      if (rff is None):
        assert approx_equal(tw, 13181.2792135)
        assert tt is None
        if (compute_gradients):
          assert approx_equal(gw, [
            (0.28910053111650247+0j), (7.7291101812048604+0j),
            (-6.2595044456051188+0j), (0.094882323440219324+0j),
            (-0.64462074968137506+0j), (-79.103928224377-94.779941011168859j),
            (0.20707720782494976-0.058418836798195622j),
            (0.23728742278013587+2.6086340153515435j),
            (-1.7239709842603137+0j), (-60.367607555268741+0j),
            (-0.63389304157017523+0.21864233709838887j),
            (-16.648813735508114+0j)])
      else:
        assert approx_equal(tw, 15400.4726889)
        assert approx_equal(tt, 2085.31183651)
        if (compute_gradients):
          assert approx_equal(gw, [
            (0.34692063733980305+0j), (9.2749322174458335+0j),
            (-0.77354489961765016+0j),
            (-94.924713869252415-113.73592921340264j),
            (0.24849264938993973-0.070102604157834744j),
            (0.28474490733616309+3.1303608184218525j),
            (-2.0687651811123766+0j), (-72.441129066322489+0j),
            (-0.76067164988421032+0.26237080451806666j),
            (-19.978576482609739+0j)])
        mlw = xray.targets_maximum_likelihood_criterion(
          f_obs=f_obs.select(~rff),
          r_free_flags=None,
          f_calc=f_calc.select(~rff),
          alpha=alpha.select(~rff),
          beta=beta.select(~rff),
          scale_factor=1,
          epsilons=epsilons.select(~rff),
          centric_flags=centric_flags.select(~rff),
          compute_gradients=compute_gradients)
        assert approx_equal(mlw.target_per_reflection(), tpr.select(~rff))
        assert approx_equal(mlw.target_work(), tw)
        assert mlw.target_test() is None
        assert approx_equal(mlw.gradients_work(), gw)
      ml = xray.targets_maximum_likelihood_criterion_hl(
        f_obs=f_obs,
        r_free_flags=rff,
        experimental_phases=experimental_phases,
        f_calc=f_calc,
        alpha=alpha,
        beta=beta,
        epsilons=epsilons,
        centric_flags=centric_flags,
        integration_step_size=5,
        compute_gradients=compute_gradients)
      tpr = ml.target_per_reflection()
      tw = ml.target_work()
      tt = ml.target_test()
      gw = ml.gradients_work()
      assert approx_equal(tpr, [
        259.4738891742619, 30871.95384262779, 4156.0852750170579,
        11.907959982009571, 28.697266914523933, 43148.119151934428,
        2.9161647472283221, 5293.4083862783209, 1838.4175626851754,
        64386.490862538834, 983.86578354410517, 7186.7802434227397])
      if (not compute_gradients):
        assert gw.size() == 0
      if (rff is None):
        assert approx_equal(tw, 13180.6763657)
        assert tt is None
        if (compute_gradients):
          assert approx_equal(gw, [
            (0.28910053111650247-4.778461360683447e-19j),
            (7.7291101812048604+9.1707854448857317e-20j),
            (-6.2595044456051188+0j),
            (0.094882323440219324+5.8011567980489493e-20j),
            (-0.64462074968137495-4.4256743068048503e-19j),
            (-79.574611695327491-94.387109668225506j),
            (0.19351597841962068-0.10406336738042934j),
            (0.22621461424524394+2.6101454947879947j),
            (-1.7239709842603133+0j), (-60.367607555268734+0j),
            (-0.64001526423037869+0.20000512034506801j),
            (-16.648813735508117+0j)])
      else:
        assert approx_equal(tw, 15400.0123154)
        assert approx_equal(tt, 2083.9966175)
        if (compute_gradients):
          assert approx_equal(gw, [
            (0.34692063733980305-5.7341536328201376e-19j),
            (9.2749322174458335+1.100494253386288e-19j),
            (-0.77354489961765005-5.3108091681658213e-19j),
            (-95.489534034392989-113.26453160187062j),
            (0.23221917410354484-0.12487604085651523j),
            (0.27145753709429277+3.132174593745594j),
            (-2.0687651811123762+0j), (-72.441129066322489+0j),
            (-0.76801831707645452+0.24000614441408163j),
            (-19.978576482609743+0j)])
        mlw = xray.targets_maximum_likelihood_criterion_hl(
          f_obs=f_obs.select(~rff),
          r_free_flags=None,
          experimental_phases=experimental_phases.select(~rff),
          f_calc=f_calc.select(~rff),
          alpha=alpha.select(~rff),
          beta=beta.select(~rff),
          epsilons=epsilons.select(~rff),
          centric_flags=centric_flags.select(~rff),
          integration_step_size=5,
          compute_gradients=compute_gradients)
        assert approx_equal(mlw.target_per_reflection(), tpr.select(~rff))
        assert approx_equal(mlw.target_work(), tw)
        assert mlw.target_test() is None
        assert approx_equal(mlw.gradients_work(), gw)

def exercise_twin_components():
  twin_a = xray.twin_component(sgtbx.rot_mx((0,1,0,1,0,0,0,0,-1)),
                               value=0.5,
                               grad=True)
  assert twin_a.value == 0.5
  assert twin_a.twin_law == sgtbx.rot_mx((0,1,0,1,0,0,0,0,-1))
  assert twin_a.grad == True
  twin_a.grad = False
  assert twin_a.grad == False
  import cPickle
  pickled_twin_a = cPickle.dumps(twin_a, cPickle.HIGHEST_PROTOCOL)
  unpickled_twin_a = cPickle.loads(pickled_twin_a)
  assert twin_a.value == unpickled_twin_a.value
  assert twin_a.twin_law == unpickled_twin_a.twin_law
  assert twin_a.grad == unpickled_twin_a.grad

  twin_b = xray.twin_component(sgtbx.rot_mx((-1,0,0,0,-1,0,0,0,-1)),
                               value=0.2,
                               grad = False)
  assert twin_b.grad == False
  twins = (twin_a, twin_b)
  xray.set_grad_twin_fraction(twins, True)
  assert twin_a.grad == True
  assert twin_b.grad == True
  assert xray.sum_twin_fractions(twins) == 0.7

def exercise_extinction_correction():
  uc = uctbx.unit_cell((13,15,17,85,95,105))
  wavelength = 1.3
  extinction_val = 0.134
  h = (1,2,3)
  fc_sq = 3.456
  ec = xray.shelx_extinction_correction(uc, wavelength, value=extinction_val)
  ec.grad = True
  c = ec.compute(h, fc_sq, True)
  # just check that the formulae are the same one...
  p = 0.001*pow(wavelength,3)*fc_sq/uc.sin_two_theta(h, wavelength)
  c_v = pow(1+extinction_val*p,-0.5)
  c_g = -0.5*fc_sq*p*pow(1+p*extinction_val,-1.5)
  assert approx_equal(c[0], c_v)
  assert approx_equal(c[1], c_g)
  #now try numerical differentiation and compare
  step = 1e-8
  c_g_n = -ec.compute(h, fc_sq, False)[0]
  ec.value += step
  c_g_n += ec.compute(h, fc_sq, False)[0]
  c_g_n = c_g_n*fc_sq/step
  assert approx_equal(c[1], c_g_n)


def run():
  exercise_extinction_correction()
  exercise_twin_components()
  exercise_scatterer_flags()
  exercise_set_selected_scatterer_grad_flags()
  exercise_set_scatterer_grad_flags()
  exercise_scatterer_flags_counts()
  exercise_xray_scatterer()
  exercise_conversions()
  exercise_gradient_flags()
  exercise_scattering_type_registry()
  exercise_rotate()
  exercise_structure_factors()
  exercise_targets()
  exercise_sampled_model_density()
  exercise_minimization_apply_shifts()
  exercise_minimization_add_gradients()
  exercise_asu_mappings()
  exercise_targets_common_results()
  exercise_targets_least_squares()
  exercise_maximum_likelihood_targets()
  print "OK"

if (__name__ == "__main__"):
  run()
