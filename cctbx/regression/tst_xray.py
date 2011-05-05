from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal, \
  is_above_limit, show_diff
from libtbx.test_utils import not_approx_equal
import random
import pickle
from cStringIO import StringIO
import sys, random, math
from itertools import count

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def exercise_scatterer():
  assert xray.scatterer(scattering_type="si4+").element_and_charge_symbols() \
      == ("Si", "4+")
  assert xray.scatterer(scattering_type="si+4").element_and_charge_symbols() \
      == ("Si", "")
  assert xray.scatterer(scattering_type="x").element_and_charge_symbols() \
      == ("", "")
  assert xray.scatterer(scattering_type="Cval").element_symbol() == "C"
  assert xray.scatterer(scattering_type="si+4").element_symbol() == "Si"
  assert xray.scatterer(scattering_type="x").element_symbol() is None

def exercise_anomalous_scatterer_group():
  scatterers = flex.xray_scatterer(
    [xray.scatterer(scattering_type=scattering_type)
     for scattering_type in ["S", "AU", "C", "S"]])
  groups = [
    xray.anomalous_scatterer_group(
      iselection=flex.size_t([1]),
      f_prime=-3,
      f_double_prime=4,
      refine=["f_double_prime", "f_prime"]),
    xray.anomalous_scatterer_group(
      iselection=flex.size_t([0,3]),
      f_prime=-1,
      f_double_prime=2,
      refine=[],
      selection_string="name S")
  ]
  assert groups[0].labels_refine() == ["f_prime", "f_double_prime"]
  assert groups[1].labels_refine() == []
  for i_pass in [0,1]:
    s = StringIO()
    for group in groups:
      group.show_summary(out=s, prefix="{.")
    assert not show_diff(s.getvalue(), """\
{.Anomalous scatterer group:
{.  Number of selected scatterers: 1
{.  f_prime:        -3
{.  f_double_prime: 4
{.  refine: f_prime f_double_prime
{.Anomalous scatterer group:
{.  Selection: "name S"
{.  Number of selected scatterers: 2
{.  f_prime:        -1
{.  f_double_prime: 2
{.  refine: None
""")
    for group in groups:
      group.copy_to_scatterers_in_place(scatterers=scatterers)
      group.extract_from_scatterers_in_place(scatterers=scatterers)
  scatterers[1].fp = -4
  scatterers[1].fdp = 5
  group = groups[0]
  group.extract_from_scatterers_in_place(scatterers=scatterers)
  assert approx_equal(group.f_prime, -4)
  assert approx_equal(group.f_double_prime, 5)
  scatterers[0].fdp = 3
  try: groups[1].extract_from_scatterers_in_place(scatterers=scatterers)
  except RuntimeError, e:
    assert str(e) == """\
Anomalous scatterer group with significantly different f_double_prime:
  Selection: "name S"
  Number of selected scatterers: 2
  f_double_prime min:  2
  f_double_prime max:  3
  f_double_prime mean: 2.5
  tolerance: 0.0001"""
  else: raise Exception_expected

def exercise_structure():
  cs1 = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp1 = crystal.special_position_settings(cs1)
  scatterers1 = flex.xray_scatterer((
    xray.scatterer("o", (0.5, 0, 0)),
    xray.scatterer("c", (0, 0, 0))))
  xs1 = xray.structure(sp1, scatterers1)
  cs2 = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp2 = crystal.special_position_settings(cs2)
  scatterers1 = flex.xray_scatterer((
    xray.scatterer("o", (0, 0, 0)),
    xray.scatterer("c", (0,   0, .5))))
  xs2 = xray.structure(sp1, scatterers1)
  assert approx_equal(list(xs1.distances(other = xs2)), [5,15])
  cs = crystal.symmetry((5.01, 5.01, 5.47, 90, 90, 120), "P 62 2 2")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", (1./2, 1./2, 0.3)),
    xray.scatterer("O1", (0.18700, -0.20700, 0.83333))))
  xs = xray.structure(sp, scatterers)
  assert xs.scatterers().size() == 2
  assert not xs.non_unit_occupancy_implies_min_distance_sym_equiv_zero()
  assert xs.n_undefined_multiplicities() == 0
  assert tuple(xs.special_position_indices()) == (0, 1)
  s = StringIO()
  xs.show_special_position_shifts(
    sites_cart_original=cs.unit_cell().orthogonalize(
      sites_frac=scatterers.extract_sites()),
    out=s,
    prefix="%^")
  assert s.getvalue() == """\
%^Number of sites at special positions: 2
%^  Minimum distance between symmetrically equivalent sites: 0.5
%^  Label   Mult   Shift    Fractional coordinates
%^  Si1       3    0.182 (  0.5000   0.5000   0.3000) original
%^        site sym 222   (  0.5000   0.5000   0.3333) exact
%^                               1/2,1/2,1/3
%^  O1        6    0.050 (  0.1870  -0.2070   0.8333) original
%^        site sym 2     (  0.1970  -0.1970   0.8333) exact
%^                        1/2*x-1/2*y,-1/2*x+1/2*y,5/6
"""
  ys = xs.deep_copy_scatterers()
  ys.add_scatterers(ys.scatterers())
  assert ys.scatterers().size() == 4
  assert xs.scatterers().size() == 2
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3)
  ys.add_scatterer(ys.scatterers()[0])
  assert ys.scatterers().size() == 5
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3, 4)
  sx = xs.primitive_setting()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 62 2 2"
  sx = xs.change_hand()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 64 2 2"
  assert approx_equal(sx.scatterers()[0].site, (-1./2, -1./2, -1./3))
  assert approx_equal(sx.scatterers()[1].site, (-0.19700, 0.19700, -0.833333))
  p1 = xs.asymmetric_unit_in_p1()
  assert p1.scatterers().size() == 2
  for i in xrange(2):
    assert p1.scatterers()[i].weight() == xs.scatterers()[i].weight()
  assert str(p1.space_group_info()) == "P 1"
  p1 = xs.expand_to_p1()
  assert p1.scatterers().size() == 9
  for i in xrange(2):
    assert p1.scatterers()[i].weight() != xs.scatterers()[i].weight()
  sh = p1.apply_shift((0.2,0.3,-1/6.))
  assert approx_equal(sh.scatterers()[0].site, (0.7,0.8,1/6.))
  assert approx_equal(sh.scatterers()[3].site, (0.3970,0.1030,2/3.))
  sl = sh[:1]
  assert sl.scatterers().size() == 1
  assert sl.scatterers()[0].label == sh.scatterers()[0].label
  sl = sh[1:4]
  assert sl.scatterers().size() == 3
  for i in xrange(3):
    assert sl.scatterers()[i].label == sh.scatterers()[i+1].label
  xs.scatterers().set_occupancies(flex.double((0.5,0.2)))
  s = xs.sort(by_value="occupancy")
  assert approx_equal(s.scatterers().extract_occupancies(), (0.2,0.5))
  assert s.scatterers()[0].label == "O1"
  assert s.scatterers()[1].label == "Si1"
  aw = xs.atomic_weights()
  assert approx_equal(aw, (28.086, 15.999))
  center_of_mass = xs.center_of_mass(atomic_weights=aw)
  assert approx_equal(center_of_mass, (1.335228, 1.071897, 2.815899))
  center_of_mass = xs.center_of_mass()
  assert approx_equal(center_of_mass, (1.335228, 1.071897, 2.815899))
  ys = xs.apply_shift(
    shift=xs.unit_cell().fractionalize([-e for e in center_of_mass]),
    recompute_site_symmetries=True)
  assert approx_equal(ys.center_of_mass(), (0,0,0))
  ys = xray.structure(xs)
  assert ys.atomic_weights().size() == 0
  pa1 = xs.principal_axes_of_inertia()
  pa2 = xs.principal_axes_of_inertia(atomic_weights=aw)
  assert approx_equal(pa1.center_of_mass(), pa2.center_of_mass())
  assert approx_equal(pa1.inertia_tensor(), pa2.inertia_tensor())
  assert approx_equal(pa1.center_of_mass(), (1.335228, 1.071897, 2.815899))
  assert approx_equal(pa1.inertia_tensor(),
                      (169.46094427548525, 76.773803949363497, 93.746443127965918,
                       7.0265499590972862, -6.3547480051846588, 84.304420337921911))
  ys = xray.structure(sp, scatterers)
  ys.scatterers()[1].occupancy = 0.5
  assert approx_equal(ys.scatterers()[1].weight(),0.25)
  ys.scatterers()[1].occupancy -= 0.1
  assert approx_equal(ys.scatterers()[1].weight(),0.2)
  gradient_flags = xray.structure_factors.gradient_flags(default=True)
  xray.set_scatterer_grad_flags(scatterers = xs.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  assert xs.n_parameters() == 14
  assert xs.n_parameters(considering_site_symmetry_constraints=True) == 9
  g = flex.vec3_double(((0.1,0.2,0.3),(0.2,0.3,0.4)))
  xs.apply_special_position_ops_d_target_d_site(g)
  assert approx_equal(g[0], (0,0,0))
  assert approx_equal(g[1], (-0.05,0.05,0))
  xs.replace_scatterers(xs.scatterers()[:1], None)
  assert xs.scatterers().size() == 1
  assert tuple(xs.special_position_indices()) == (0,)
  reg = ys.scattering_type_registry(table="electron")
  assert reg.gaussian("Si").n_terms() == 5
  s = StringIO()
  reg.show(out=s)
  assert not show_diff(s.getvalue(), """\
Number of scattering types: 2
  Type Number    sf(0)   Gaussians
   Si      1      5.83       5
   O       1      1.98       5
  sf(0) = scattering factor at diffraction angle 0.
""")
  reg = ys.scattering_type_registry(table="wk1995")
  assert reg.gaussian("Si").n_terms() == 5
  reg = ys.scattering_type_registry(table="it1992")
  assert reg.gaussian("Si").n_terms() == 4
  reg = ys.scattering_type_registry(
    custom_dict={"Si":eltbx.xray_scattering.gaussian(1)})
  assert reg.gaussian("Si").n_terms() == 0
  s = StringIO()
  reg.show_summary(out=s)
  assert s.getvalue().strip() == "O:4+c*1 Si:0+c*1"
  s = StringIO()
  reg.show(out=s)
  assert not show_diff(s.getvalue(), """\
Number of scattering types: 2
  Type Number    sf(0)   Gaussians
   O       1      8.00       4+c
   Si      1      1.00       0+c
  sf(0) = scattering factor at diffraction angle 0.
""")
  s = StringIO()
  reg.show(out=s, show_sf0=False)
  assert s.getvalue() == """\
Number of scattering types: 2
  Type Number   Gaussians
   O       1        4+c
   Si      1        0+c
"""
  s = StringIO()
  reg.show(out=s, show_gaussians=False)
  assert s.getvalue() == """\
Number of scattering types: 2
  Type Number    sf(0)
   O       1      8.00
   Si      1      1.00
  sf(0) = scattering factor at diffraction angle 0.
"""
  s = StringIO()
  reg.show(out=s, show_sf0=False, show_gaussians=False)
  assert s.getvalue() == """\
Number of scattering types: 2
  Type Number
   O       1
   Si      1
"""
  wd = reg.wilson_dict()
  assert len(wd) == 2
  assert wd["O"] == 1
  assert wd["Si"] == 1
  tgd = reg.as_type_gaussian_dict()
  assert len(tgd) == 2
  assert tgd["O"].n_terms() == 4
  assert tgd["Si"].n_terms() == 0
  #
  reg = ys.scattering_type_registry()
  assert reg.type_index_pairs_as_dict() == {"Si": 0, "O": 1}
  assert list(reg.unique_indices(scatterers=ys.scatterers())) == [0,1]
  assert reg.type_count_dict() == {"Si": 1, "O": 1}
  #
  am = xs.asu_mappings(buffer_thickness=1)
  assert am.mappings().size() == xs.scatterers().size()
  rs = p1.random_shift_sites(max_shift_cart=0.2)
  assert flex.max(flex.abs(p1.difference_vectors_cart(rs).as_double())) <= 0.2
  assert approx_equal(p1.rms_difference(p1), 0)
  assert approx_equal(rs.rms_difference(rs), 0)
  assert p1.rms_difference(rs) > 0
  for s in [xs, ys, p1, rs]:
    p = pickle.dumps(s)
    l = pickle.loads(p)
    assert l.scatterers().size() == s.scatterers().size()
    assert l.special_position_indices().all_eq(s.special_position_indices())
  xs0 = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10,10,10,90,90,90),
      space_group_symbol="P 2 2 2"))
  xs0.add_scatterer(xray.scatterer(label="C", site=(0.1,0.1,0.1)))
  assert xs0.site_symmetry_table().get(0).is_point_group_1()
  xs1 = xs0.apply_shift(shift=(-0.08,-0.1,0), recompute_site_symmetries=False)
  assert xs1.site_symmetry_table().get(0).is_point_group_1()
  xs2 = xs0.apply_shift(shift=(-0.08,-0.1,0), recompute_site_symmetries=True)
  assert str(xs2.site_symmetry_table().get(0).special_op()) == "0,0,z"
  assert approx_equal(xs1.scatterers()[0].site, (0.02, 0, 0.1))
  assert approx_equal(xs2.scatterers()[0].site, (0, 0, 0.1))
  assert list(xs1.deep_copy_scatterers().special_position_indices()) == []
  assert list(xs2.deep_copy_scatterers().special_position_indices()) == [0]
  assert list(xs1[:].special_position_indices()) == []
  assert list(xs2[:].special_position_indices()) == [0]
  xs1.add_scatterer(xs1.scatterers()[0])
  assert list(xs1[:].special_position_indices()) == [1]
  xs1.add_scatterer(xs1.scatterers()[0], xs1.site_symmetry_table().get(0))
  assert list(xs1[:].special_position_indices()) == [1]
  xs1.add_scatterers(xs1.scatterers())
  assert list(xs1[:].special_position_indices()) == [1,3,4,5]
  xs1.add_scatterers(xs1.scatterers(), xs1.site_symmetry_table())
  assert list(xs1[:].special_position_indices()) == [1,3,4,5,7,9,10,11]
  for selection in [flex.size_t([1,4,6]),
                    flex.bool([False,True,False,False,True,False,
                               True,False,False,False,False,False])]:
    xs2 = xs1.select(selection=selection)
    assert xs2.scatterers().size() == 3
    assert list(xs2.special_position_indices()) == [0,1]
    if (isinstance(selection, flex.bool)):
      xs2 = xs1.select(selection=selection, negate=True)
      assert xs2.scatterers().size() == 9
  xs2 = xs1[2::2]
  assert xs2.scatterers().size() == 5
  assert list(xs2[:].special_position_indices()) == [1,4]
  xs1.replace_scatterers(xs2.scatterers(), None)
  assert list(xs1[:].special_position_indices()) == [0,1,2,3,4]
  xs1.replace_scatterers(xs2.scatterers(), xs2.site_symmetry_table())
  assert list(xs1[:].special_position_indices()) == [1,4]
  xs2 = xs1.asymmetric_unit_in_p1()
  assert xs2.unit_cell().is_similar_to(xs1.unit_cell())
  assert xs2.space_group().order_z() == 1
  assert list(xs2[:].special_position_indices()) == [1,4]
  for append_number_to_labels in [False, True]:
    xs2 = xs1.expand_to_p1(append_number_to_labels=append_number_to_labels)
    assert xs2.scatterers().size() == 16
    xs2 = xray.structure(xs1, xs1.scatterers())
    xs2 = xs2.expand_to_p1(append_number_to_labels=append_number_to_labels)
    assert xs2.scatterers().size() == 10
  cb_op = sgtbx.change_of_basis_op("z,x,y")
  xs2 = xs1.change_basis(cb_op)
  for i,sc in enumerate(xs2.scatterers()):
    assert sc.multiplicity() > 0
    ss = xs2.site_symmetry_table().get(i)
    assert ss.multiplicity() == sc.multiplicity()
    assert ss.multiplicity() * ss.n_matrices() == xs2.space_group().order_z()
  xs2 = xs2.expand_to_p1()
  for sc in xs2.scatterers():
    assert sc.multiplicity() == 1
  assert xs2.scatterers().size() == 16
  assert approx_equal(xs2.scatterers()[0].site, (0.1, 0.02, 0))
  assert approx_equal(xs2.scatterers()[4].site, (0.1, 0, 0))
  sx = ys.sites_mod_positive()
  assert approx_equal(sx.scatterers()[1].site, (0.197,0.803,0.8333333))
  sx = ys.sites_mod_short()
  assert approx_equal(sx.scatterers()[1].site, (0.197,-0.197,-0.1666667))
  xs1.scatterers().set_occupancies(flex.random_double(size=5))
  xs2 = xs1.sort(by_value="occupancy")
  assert xs2.special_position_indices().size() == 2
  xs2.set_sites_frac(xs2.sites_frac()+(0.1,0.2,0.3))
  xs2.set_sites_cart(xs2.sites_cart()+(1,2,3))
  assert approx_equal(xs2.extract_u_iso_or_u_equiv(), [0,0,0,0,0])
  i = xs2.special_position_indices()[0]
  assert approx_equal(xs2.scatterers()[i].site, (0.2, 0.4, 0.7))
  assert list(xs2.is_positive_definite_u()) == [False]*5
  assert list(xs2.is_positive_definite_u(u_cart_tolerance=0)) == [True]*5
  xs2.tidy_us(u_min=1,u_max=888)
  assert approx_equal(xs2.scatterers().extract_u_iso(), [1]*5)
  xs2.shift_us(u_shift=2)
  assert approx_equal(xs2.scatterers().extract_u_iso(), [3]*5)
  xs2.shift_us(b_shift=-adptbx.u_as_b(2))
  assert approx_equal(xs2.scatterers().extract_u_iso(), [1]*5)
  xs2.scatterers().set_occupancies(flex.double(xs2.scatterers().size(), 1.0))
  xs2.shift_occupancies(q_shift = 1.0)
  assert approx_equal(xs2.scatterers().extract_occupancies(),
                      [2.0, 2.0, 2.0, 2.0, 2.0])
  xs2.shift_occupancies(q_shift = -1.0, selection = flex.size_t([0,3]))
  assert approx_equal(xs2.scatterers().extract_occupancies(),
                      [1.0, 2.0, 2.0, 1.0, 2.0])
  xs2.apply_symmetry_sites()
  assert approx_equal(xs2.scatterers()[i].site, (0, 0, 0.7))
  xs2.apply_symmetry_u_stars()
  s = StringIO()
  xs1.show_distances(distance_cutoff=0.1, out=s)
  assert s.getvalue() == """\
C  pair count:   2       <<  0.0200,  0.0000,  0.1000>>
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
C  pair count:   1       <<  0.0000,  0.0000,  0.1000>>
  C:   0.0000             (  0.0000,  0.0000,  0.1000)
C  pair count:   2       <<  0.0200,  0.0000,  0.1000>>
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
C  pair count:   2       <<  0.0200,  0.0000,  0.1000>>
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
  C:   0.0000             (  0.0200,  0.0000,  0.1000)
C  pair count:   1       <<  0.0000,  0.0000,  0.1000>>
  C:   0.0000             (  0.0000,  0.0000,  0.1000)
"""
  quartz = xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))
  assert approx_equal(quartz.mean_scattering_density(), 0.184934969936)
  s = StringIO()
  quartz.show_angles(distance_cutoff=2, out=s)
  assert not show_diff(s.getvalue(), """\
O*1   Si    O*2    101.31
O*3   Si    O*2    111.31
O*3   Si    O*1    116.13
O*4   Si    O*2    116.13
O*4   Si    O*1    111.31
O*4   Si    O*3    101.31
Si*5  O     Si*6   146.93
*1 x-y,x,z-2/3
*2 -y,x-y,z-1/3
*3 y+1,-x+y+1,z-1/3
*4 -x+y+1,-x+1,z-2/3
*5 y,-x+y,z+2/3
*6 -x+y,-x,z+1/3
""")
### shake_adp()
  cs = crystal.symmetry((5.01, 6.01, 5.47, 60, 80, 120), "P1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer([xray.scatterer("o")]*100)
  selection = flex.bool()
  rd = flex.mersenne_twister(seed=0).random_double
  for i_sc, sc in enumerate(scatterers):
    scale = scatterers.size()/(1.+i_sc)
    sc.u_iso = 1.0 * scale
    sc.u_star = (1.0 * scale, 2.0 * scale, 3.0 * scale,
                 4.0 * scale, 5.0 * scale, 6.0 * scale)
    sc.flags.set_use_u_iso(rd() > 0.5)
    sc.flags.set_use_u_aniso(rd() > 0.5)
    selection.append(rd() > 0.5)
  xs = xray.structure(sp, scatterers)
  for sel in [None, selection]:
    for b_min, b_max, spread in zip([None,10.0], [None,20.0], [10.0,0.0]):
      for keep_anisotropic in [True, False]:
        xs_mod = xs.deep_copy_scatterers()
        xs_mod.shake_adp(keep_anisotropic = keep_anisotropic,
                         b_max = b_max, b_min = b_min, selection=sel, spread = spread)
        if(sel is None): sel = flex.bool(xs.scatterers().size(), True)
        for sc,sc_mod,s in zip(xs.scatterers(),xs_mod.scatterers(),sel):
          assert sc.flags.use_u_iso()   == sc_mod.flags.use_u_iso()
          assert sc.flags.use_u_aniso() == sc_mod.flags.use_u_aniso()
          if(sc.flags.use_u_iso() and s):
            assert not_approx_equal(abs(sc.u_iso - sc_mod.u_iso), 0)
          else:
            assert approx_equal(sc.u_iso, sc_mod.u_iso)
          if(keep_anisotropic):
            assert approx_equal(sc.u_star, sc_mod.u_star)
          else:
            if(sc.flags.use_u_aniso() and s):
              a = flex.double(sc.u_star)
              b = flex.double(sc_mod.u_star)
              # quick-and-dirty test, likely to fail if random seed
              # is changed
              assert is_above_limit(
                value=flex.max(flex.abs((a-b))), limit=0.1)
              assert is_above_limit(
                value=flex.min(flex.abs((a-b))), limit=0.001)
            else:
              assert approx_equal(sc.u_star, sc_mod.u_star)
### shake_occupancies()
  cs = crystal.symmetry((5.01, 6.01, 5.47, 60, 80, 120), "P1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer([xray.scatterer("o")]*100)
  selection = flex.bool()
  for sc in scatterers:
    sc.occupancy = random.random()
    selection.append(random.choice((False,True)))
  xs = xray.structure(sp, scatterers)
  for sel in [None, selection]:
    xs_mod = xs.deep_copy_scatterers()
    xs_mod.shake_occupancies(selection = sel)
    if(sel is None): sel = flex.bool(xs.scatterers().size(), True)
    for sc, sc_mod, s in zip(xs.scatterers(),xs_mod.scatterers(), sel):
      if(s):
        assert abs(sc.occupancy - sc_mod.occupancy) > 1.e-4
      else:
        assert approx_equal(sc.occupancy, sc_mod.occupancy)
### shake_adp_if_all_equal()
  cs = crystal.symmetry((5.01, 6.01, 5.47, 60, 80, 120), "P1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer([xray.scatterer("o", u=0.5)]*100)
  xs = xray.structure(sp, scatterers)
  xs_mod = xs.deep_copy_scatterers()
  assert xs_mod.shake_adp_if_all_equal()
  assert not xs.scatterers().extract_u_iso().all_eq(
    xs_mod.scatterers().extract_u_iso())
  xs.scatterers()[3].u_iso = 0.1
  xs_mod = xs.deep_copy_scatterers()
  assert not xs_mod.shake_adp_if_all_equal()
  assert xs.scatterers().extract_u_iso().all_eq(
    xs_mod.scatterers().extract_u_iso())
# exercise set_b_iso()
  cs = crystal.symmetry((5.01, 6.01, 5.47, 60, 80, 120), "P1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("C", (1./2, 1./2, 0.3)),
    xray.scatterer("C", (0.18700, -0.20700, 0.83333))))
  xs = xray.structure(sp, scatterers)
  uc = xs.unit_cell()
  b_iso_value = 25.0
  xs.set_b_iso(value = b_iso_value)
  result = flex.double([25,25])
  assert approx_equal(xs.scatterers().extract_u_iso()/adptbx.b_as_u(1), result)
  b_iso_values = flex.double([7,9])
  xs.set_b_iso(values = b_iso_values)
  assert approx_equal(xs.scatterers().extract_u_iso()/adptbx.b_as_u(1), b_iso_values)
  #
  xs.scatterers().set_u_iso(flex.double([0.1,0.2]),
                            flex.bool(xs.scatterers().size(), True), uc)
  assert xs.scatterers().count_anisotropic() == 0
  xs.convert_to_anisotropic()
  assert xs.scatterers().count_anisotropic() == 2
  assert approx_equal(xs.scatterers().extract_u_cart(unit_cell=xs.unit_cell()),
                      [[0.1]*3+[0]*3, [0.2]*3+[0]*3])
  xs.convert_to_isotropic()
  assert xs.scatterers().count_anisotropic() == 0
  assert approx_equal(xs.scatterers().extract_u_iso(), [0.1,0.2])
  #
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  uc = cs.unit_cell()
  scatterers = flex.xray_scatterer((
    xray.scatterer("o", site=(0.5, 0, 0),u=1.0),
    xray.scatterer("o", site=(0.5, 1.0, 0),u=0.1),
    xray.scatterer("o", site=(0.5, 1.0, 10),u=0.7),
    xray.scatterer("n", site=(0.5,-1.0, 0),u=adptbx.u_cart_as_u_star(uc,(1,2,3,-.3,-2,1))),
    xray.scatterer("c", site=(0, 0, 0),u=adptbx.u_cart_as_u_star(uc,(6,7,9,2,3,-.7)))))
  xs = xray.structure(sp, scatterers)
  b_isos = [adptbx.u_as_b(i) for i in xs.extract_u_iso_or_u_equiv()]
  assert not_approx_equal(b_isos, [20.0, 20.0, 20.0, 20.0, 20.0])
  xs.set_b_iso(value=20)
  b_isos = [adptbx.u_as_b(i) for i in xs.extract_u_iso_or_u_equiv()]
  assert approx_equal(b_isos, [20.0, 20.0, 20.0, 20.0, 20.0])
### scale_adp:
  cs = crystal.symmetry((5.01, 6.01, 5.47, 60, 80, 120), "P1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer([xray.scatterer("o")]*100)
  selection = flex.bool()
  rd = flex.mersenne_twister(seed=0).random_double
  for sc in scatterers:
    sc.u_iso = 1.0 * rd()
    sc.u_star = (1.0 * rd(), 2.0 * rd(), 3.0 * rd(),
                 4.0 * rd(), 5.0 * rd(), 6.0 * rd())
    sc.flags.set_use_u_iso(rd() > 0.5)
    sc.flags.set_use_u_aniso(rd() > 0.5)
    selection.append(rd() > 0.5)
  xs = xray.structure(sp, scatterers)
  xs_dc = xs.deep_copy_scatterers()
  xs_dc.scale_adp(factor=2, selection=selection)
  for i_seq,sc,sc_dc,sel in zip(count(), xs.scatterers(), xs_dc.scatterers(),
                                selection):
    if(sel and sc.flags.use()):
      if(sc.flags.use_u_iso()):
        if(sc.u_iso != 0.0):
          assert approx_equal(sc_dc.u_iso / sc.u_iso, 2.0)
      if(sc.flags.use_u_aniso()):
        for i in xrange(6):
          if(sc.u_star[i] != 0.0):
            assert approx_equal(sc_dc.u_star[i] / sc.u_star[i], 2.0)

### shake_sites
  selection_ = flex.bool([random.choice((False,True)) for i in xrange(500)])
  xs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["N"]*500,
    unit_cell        = (10, 20, 30, 70, 80, 120))
  mt = flex.mersenne_twister(seed=0)
  errors = [0.0, 0.01, 0.1, 0.5, 1.5, 3.0, 10.0]
  for selection in [None, selection_]:
    for error in errors:
      xs_shaked = xs.deep_copy_scatterers()
      xs_shaked.shake_sites_in_place(
        mean_distance=error,
        selection=selection,
        random_double=[None, mt.random_double][int(error<0.2)])
      sites_cart_xs        = xs.sites_cart()
      sites_cart_xs_shaked = xs_shaked.sites_cart()
      if(selection is None): selection=flex.bool(xs.scatterers().size(),True)
      mean_err = flex.mean(
        flex.sqrt((sites_cart_xs.select(selection) -
                   sites_cart_xs_shaked.select(selection)).dot()))
      assert approx_equal(error, mean_err, 0.001)
      dummy = ~selection
      if(dummy.count(True) > 0):
        mean_err_fixed = flex.mean(
          flex.sqrt((sites_cart_xs.select(~selection) -
                     sites_cart_xs_shaked.select(~selection)).dot()))
        assert approx_equal(mean_err_fixed, 0.0)
  ### random remove sites
  for fraction in xrange(1, 99+1, 10):
    fraction /= 100.
    selection = xs.random_remove_sites_selection(fraction = fraction)
    deleted = selection.count(False) / float(selection.size())
    retained= selection.count(True)  / float(selection.size())
    assert approx_equal(fraction, deleted)
    assert approx_equal(1-deleted, retained)
  #
  xs.scatterers()[0] = xs.scatterers()[0].customized_copy()
  assert approx_equal(xs.scatterers()[0].weight(), 1.0)
  xs.re_apply_symmetry(i_scatterer=0)
  assert approx_equal(xs.scatterers()[0].weight(), 1.0)
# apply_rigid_body_shift
  selection_=flex.bool([random.choice((False,True))
    for i in xrange(100)]).iselection()
  xs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["N"]*100,
    unit_cell        = (10, 20, 30, 70, 80, 120))
  for selection in [None, selection_]:
    for r in [[1.,2.,3.], [0.,0.,0.]]:
      xs_mod = xs.deep_copy_scatterers()
      xs_mod.apply_rigid_body_shift(
        rot       = (1,0,0,0,1,0,0,0,1),
        trans     = [r[0],r[1],r[2]],
        selection = selection)
      d = math.sqrt(r[0]**2+r[1]**2+r[2]**2)
      assert approx_equal(
        d, xs.mean_distance(other = xs_mod, selection = selection))
  selection_=flex.bool([random.choice((False,True)) for i in xrange(100)])
  xs_mod = xs.deep_copy_scatterers()
  xs_mod.apply_rigid_body_shift(
    rot=(0.999,-0.017,0.035,0.019,0.998,-0.052,-0.033,0.052,0.998),
    trans=[1., 2., 3.],
    selection=selection_.iselection())
  assert xs.mean_distance(other = xs_mod, selection = selection_) > 1.0
  assert approx_equal(
    xs.mean_distance(other = xs_mod, selection = ~selection_), 0.0)
  #
  assert xs.scatterers().size() == xs.all_selection().size() == \
    xs.all_selection().count(True)
  #
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(20,30,40,90,90,90),
      space_group_symbol="P222"))
  bs = xs.orthorhombic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(7, 7, 7, 90, 90, 90)"
  bs = xs.cubic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(7, 7, 7, 90, 90, 90)"
  xs.add_scatterer(xray.scatterer(label="S1", site=[0.1,0.2,-0.3]))
  bs = xs.orthorhombic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(7, 7, 7, 90, 90, 90)"
  bs = xs.cubic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(7, 7, 7, 90, 90, 90)"
  xs.add_scatterer(xray.scatterer(label="S1", site=[-0.1,-0.2,0.3]))
  bs = xs.orthorhombic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(11, 19, 31, 90, 90, 90)"
  bs = xs.cubic_unit_cell_around_centered_scatterers(buffer_size=3.5)
  assert str(bs.unit_cell()) == "(31, 31, 31, 90, 90, 90)"
  #
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(3,4,5,90,90,90),
      space_group_symbol="Pmmm"))
  xs.add_scatterer(xray.scatterer("C1", site=(0.5, 0.5, 0.5)))
  xs.add_scatterer(xray.scatterer("C2", site=(0.7, 0.7, 0.7)))
  xs.add_scatterer(xray.scatterer("C3", site=(0.3, 0.3, 0.3)))
  c1, c2, c3 = xs.scatterers()
  c1.flags.set_grad_site(True)
  c2.flags.set_use_u_iso(True)
  c3.flags.set_grad_u_aniso(True)
  grad_flags = xs.scatterer_flags()
  assert ([ f.bits for f in grad_flags ]
          ==
          [ sc.flags.bits for sc in xs.scatterers() ])
  #
  from cctbx.eltbx import henke, sasaki, wavelengths
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(3,4,5,90,90,90),
      space_group_symbol="Pmmm"))
  xs.add_scatterer(xray.scatterer("C1", site=(0.5, 0.5, 0.5)))
  xs.add_scatterer(xray.scatterer("C2", site=(0.7, 0.7, 0.7)))
  xs.add_scatterer(xray.scatterer("S3", site=(0.3, 0.3, 0.3)))
  xs1 = xs.deep_copy_scatterers()
  xs1.set_inelastic_form_factors(wavelengths.characteristic("Mo"), "sasaki")
  for sc in xs1.scatterers():
    assert sc.flags.use_fp_fdp() == True
    assert sc.fp != 0
    assert sc.fdp != 0
    sc_sasaki = sasaki.table(sc.element_symbol())
    sc_fp_fdp_sasaki = sc_sasaki.at_angstrom(
      wavelengths.characteristic('Mo').as_angstrom())
    assert approx_equal(sc.fp, sc_fp_fdp_sasaki.fp())
    assert approx_equal(sc.fdp, sc_fp_fdp_sasaki.fdp())
  xs2 = xs.deep_copy_scatterers()
  xs2.set_inelastic_form_factors(0.71073, "henke") # angstrom
  for sc in xs2.scatterers():
    assert sc.flags.use_fp_fdp() == True
    assert sc.fp != 0
    assert sc.fdp != 0
    sc_henke = henke.table(sc.element_symbol())
    sc_fp_fdp_henke = sc_henke.at_angstrom(0.71073)
    assert approx_equal(sc.fp, sc_fp_fdp_henke.fp())
    assert approx_equal(sc.fdp, sc_fp_fdp_henke.fdp())
  #
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      (5.01,5.01,5.47,90,90,120), "P6222"),
    scatterers=flex.xray_scatterer([
      xray.scatterer("Si", (1/2.,1/2.,1/3.)),
      xray.scatterer("C", (0.1234, 0.5432, 0.4321)),
      xray.scatterer("O", (0.197,-0.197,0.83333))]))
  s = StringIO()
  xs.show_scatterers(f=s)
  assert not show_diff(s.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
Si   Si     3 ( 0.5000  0.5000  0.3333) 1.00 0.0000 [ - ]
C    C     12 ( 0.1234  0.5432  0.4321) 1.00 0.0000 [ - ]
O    O      6 ( 0.1970 -0.1970  0.8333) 1.00 0.0000 [ - ]
""")
  s = StringIO()
  xs.show_scatterers(f=s, special_positions_only=True)
  assert not show_diff(s.getvalue(), """\
Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
Si   Si     3 ( 0.5000  0.5000  0.3333) 1.00 0.0000 [ - ]
O    O      6 ( 0.1970 -0.1970  0.8333) 1.00 0.0000 [ - ]
""")
  #
  xs_p1 = xs.customized_copy(
    space_group_info=sgtbx.space_group_info(symbol="P1"),
    non_unit_occupancy_implies_min_distance_sym_equiv_zero=True)
  assert xs_p1.space_group_info().type().number() == 1
  assert xs_p1.non_unit_occupancy_implies_min_distance_sym_equiv_zero()

def exercise_closest_distances():
  xs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["N"]*3,
    unit_cell        = (10, 20, 30, 70, 80, 120))
  xs_other = xs.deep_copy_scatterers()
  result = xs.closest_distances(sites_frac = xs_other.sites_frac(),
    distance_cutoff = 6)
  assert approx_equal(result.smallest_distances, [0.0, 0.0, 0.0])
  xs_other = xs_other.translate(x=1,y=2,z=3)
  result = xs.closest_distances(sites_frac = xs_other.sites_frac(),
    distance_cutoff = 6)
  assert not_approx_equal(result.smallest_distances, [0.0, 0.0, 0.0])

def exercise_set_occupancies():
  xs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["N"]*5,
    unit_cell        = (10, 20, 30, 70, 80, 120))
  occ = xs.scatterers().extract_occupancies()
  assert occ.all_eq(1.0)
  xs.set_occupancies(value = 2)
  occ = xs.scatterers().extract_occupancies()
  assert occ.all_eq(2.0)
  xs.set_occupancies(
    value = -1, selection = flex.bool([True,True,False,True,False]))
  occ = xs.scatterers().extract_occupancies()
  assert approx_equal(occ, [-1.0, -1.0, 2.0, -1.0, 2.0])

def exercise_u_base():
  d_min = 9
  grid_resolution_factor = 1/3.
  for quality_factor in (1,2,4,8,10,100,200,1000):
    u_base = xray.calc_u_base(d_min, grid_resolution_factor, quality_factor)
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        d_min=d_min,
        grid_resolution_factor=grid_resolution_factor,
        u_base=u_base))
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        d_min=d_min,
        grid_resolution_factor=grid_resolution_factor,
        b_base=adptbx.u_as_b(u_base)))
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        quality_factor=quality_factor))
  #
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  uc = cs.unit_cell()
  scatterers = flex.xray_scatterer((
    xray.scatterer("o", site=(0.5, 0, 0),u=1.0),
    xray.scatterer("o", site=(0.5, 1.0, 0),u=0.1),
    xray.scatterer("o", site=(0.5, 1.0, 10),u=0.7),
    xray.scatterer("n", site=(0.5,-1.0, 0),u=adptbx.u_cart_as_u_star(uc,(1,2,3,0,0,0))),
    xray.scatterer("c", site=(0, 0, 0),u=adptbx.u_cart_as_u_star(uc,(6,7,9,0,0,0)))))
  xs = xray.structure(sp, scatterers)
  assert xs.n_grad_u_iso()==0
  assert xs.n_grad_u_aniso()==0
  xray.set_scatterer_grad_flags(scatterers = xs.scatterers(),
                                u_iso      = True,
                                u_aniso    = True)
  assert xs.n_grad_u_iso()==3
  assert xs.n_grad_u_aniso()==2
  assert xs.use_u_iso().count(True) == 3
  assert xs.use_u_aniso().count(True) == 2
  answer = [(1.0, 1.0, 1.0), (0.1, 0.1, 0.1), (0.7, 0.7, 0.7), (3., 2., 1.),
            (9., 7., 6.)]
  assert approx_equal(answer, list(xs.scatterers().u_cart_eigenvalues(uc)))
  assert approx_equal(list(xs.scatterers().anisotropy(uc)),
                      [1.0, 1.0, 1.0, 1./3, 6./9])

def exercise_from_scatterers_direct(space_group_info,
                                    element_type,
                                    allow_mix,
                                    n_elements=5,
                                    volume_per_atom=1000,
                                    d_min=3,
                                    anomalous_flag=0,
                                    use_u_iso =0,
                                    use_u_aniso =0,
                                    verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=[element_type]*n_elements,
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=True,
    random_f_prime_d_min=d_min-1,
    random_f_prime_scale=0.6,
    random_f_double_prime=anomalous_flag,
    use_u_iso = True,
    use_u_aniso = True,
    random_u_iso = True,
    random_u_iso_scale=.3,
    random_u_cart_scale=.3,
    random_u_iso_min = 0.0,
    random_occupancy=True)
  random_structure.random_modify_adp_and_adp_flags_2(
    scatterers         = structure.scatterers(),
    use_u_iso          = use_u_iso,
    use_u_aniso        = use_u_aniso,
    allow_mix          = allow_mix)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_obs_exact = structure.structure_factors(
    d_min=d_min, algorithm="direct",
    cos_sin_table=False).f_calc()
  assert f_obs_exact.anomalous_flag() == anomalous_flag
  f_obs_simple = xray.ext.structure_factors_simple(
    f_obs_exact.unit_cell(),
    f_obs_exact.space_group(),
    f_obs_exact.indices(),
    structure.scatterers(),
    structure.scattering_type_registry()).f_calc()
  if (0 or verbose):
    for i,h in enumerate(f_obs_exact.indices()):
      print h
      print f_obs_simple[i]
      print f_obs_exact.data()[i]
      if (abs(f_obs_simple[i]-f_obs_exact.data()[i]) >= 1.e-10):
        print "MISMATCH"
      print
  mismatch = flex.max(flex.abs(f_obs_exact.data() - f_obs_simple))
  assert mismatch < 1.e-10, mismatch
  f_obs_table = f_obs_exact.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct",
    cos_sin_table=True).f_calc()
  ls = xray.targets_least_squares_residual(
    abs(f_obs_exact).data(), f_obs_table.data(), False, 1)
  if (0 or verbose):
    print "r-factor:", ls.target()
  assert ls.target() < 1.e-4

def exercise_f_obs_minus_xray_structure_f_calc(
  space_group_info,
  d_min=3,
  verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["C"]*3,
    volume_per_atom=1000,
    min_distance=5,
    general_positions_only=True,
    random_u_iso=False)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_obs_exact = structure.structure_factors(
    d_min=d_min, algorithm="direct",
    cos_sin_table=False).f_calc()
  two_f_obs_minus_f_calc=abs(f_obs_exact).f_obs_minus_xray_structure_f_calc(
    f_obs_factor=2,
    xray_structure=structure,
    structure_factor_algorithm="direct",
    cos_sin_table=False)
  phase_error = two_f_obs_minus_f_calc.mean_weighted_phase_error(
    phase_source=f_obs_exact)
  if (0 or verbose):
    print "%.2f" % phase_error
  assert approx_equal(phase_error, 0)
  two_f_obs_minus_f_calc=abs(f_obs_exact).f_obs_minus_xray_structure_f_calc(
    f_obs_factor=2,
    xray_structure=structure[:-1],
    structure_factor_algorithm="direct",
    cos_sin_table=False)
  fft_map = two_f_obs_minus_f_calc.fft_map()
  fft_map.apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  density_at_sites = [real_map.eight_point_interpolation(scatterer.site)
                      for scatterer in structure.scatterers()]
  try:
    assert min(density_at_sites[:-1]) > 6.9
    assert density_at_sites[-1] > 2.5
  except AssertionError:
    print "density_at_sites:", density_at_sites
    raise

def exercise_n_gaussian(space_group_info, verbose=0):
  structure_5g = random_structure.xray_structure(
    space_group_info,
    elements=["H", "C", "N", "O", "S"]*3)
  if (0 or verbose):
    structure_5g.show_summary().show_scatterers()
  structure_4g = structure_5g.deep_copy_scatterers()
  structure_2g = structure_5g.deep_copy_scatterers()
  structure_5g.scattering_type_registry(table="wk1995")
  structure_4g.scattering_type_registry(table="it1992")
  structure_2g.scattering_type_registry(
    custom_dict=eltbx.xray_scattering.two_gaussian_agarwal_isaacs.table)
  for gaussian in \
      structure_5g.scattering_type_registry().unique_gaussians_as_list():
    assert gaussian.n_terms() == 5
  for gaussian in \
      structure_4g.scattering_type_registry().unique_gaussians_as_list():
    assert gaussian.n_terms() == 4
  for gaussian in \
      structure_2g.scattering_type_registry().unique_gaussians_as_list():
    assert gaussian.n_terms() == 2
  d_min = 1
  f_calc_5g = structure_5g.structure_factors(
    d_min=d_min,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  f_calc_4g = f_calc_5g.structure_factors_from_scatterers(
    xray_structure=structure_4g,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  f_calc_2g = f_calc_5g.structure_factors_from_scatterers(
    xray_structure=structure_2g,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  for n,f_calc_ng in ((4,f_calc_4g), (2,f_calc_2g)):
    ls = xray.targets_least_squares_residual(
      abs(f_calc_5g).data(), f_calc_ng.data(), False, 1)
    if (0 or verbose):
      print "%d-gaussian r-factor:" % n, ls.target()
    if (n == 2):
      assert ls.target() < 0.002
    else:
      assert ls.target() < 0.0002
  #
  for element in ["H", "D", "T"]:
    structure = random_structure.xray_structure(
      space_group_info, elements=[element])
    ugs = structure.scattering_type_registry(table="n_gaussian") \
        .unique_gaussians_as_list()
    assert len(ugs) == 1
    assert ugs[0].n_terms() == 6
    s = StringIO()
    ugs[0].show(f=s)
    assert not show_diff(s.getvalue(), """\
a: -1.0938988 0.76752101 0.44291771 0.42681501 0.35006501 0.10647464
b: 1.7298482 2.0196679 1.4769121 9.3088777 20.966682 44.631255
c: 0
""")
    ugs = structure.scattering_type_registry(table="it1992") \
        .unique_gaussians_as_list()
    assert len(ugs) == 1
    assert ugs[0].n_terms() == 4
    s = StringIO()
    ugs[0].show(f=s)
    assert not show_diff(s.getvalue(), """\
a: 0.493002 0.32291201 0.140191 0.04081
b: 10.5109 26.1257 3.14236 57.799702
c: 0.0030380001
""")
    ugs = structure.scattering_type_registry(table="wk1995") \
        .unique_gaussians_as_list()
    assert len(ugs) == 1
    assert ugs[0].n_terms() == 5
    s = StringIO()
    ugs[0].show(f=s)
    assert not show_diff(s.getvalue().replace("e-005","e-05"), """\
a: -0.11710366 0.0093485946 0.27006859 0.28434139 0.5528717
b: 3.0598466 0.74655777 3.2917862 32.645653 11.546356
c: 0
""")

def run_call_back(flags, space_group_info):
  if (1):
    for element_type in ("Se", "const"):
      for anomalous_flag in [0,1]:
        for (use_u_iso,use_u_aniso) in [(True,True),(False,True),
                                        (True,False),(False,False)]:
          for with_shift in [0,1]:
            if (with_shift):
              sgi = debug_utils.random_origin_shift(space_group_info)
            else:
              sgi = space_group_info
            for allow_mix in [False, True]:
              exercise_from_scatterers_direct(
                space_group_info=sgi,
                element_type=element_type,
                anomalous_flag=anomalous_flag,
                use_u_iso = use_u_iso,
                use_u_aniso = use_u_aniso,
                verbose=flags.Verbose,
                allow_mix = allow_mix)
  if (1):
    exercise_n_gaussian(
      space_group_info=space_group_info)
  if (1):
    exercise_f_obs_minus_xray_structure_f_calc(
      space_group_info=space_group_info)

def exercise_concatenate_inplace():
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("o", (0.5, 0, 0)),
    xray.scatterer("c", (0, 0, 0))))
  xs = xray.structure(sp, scatterers)
  #
  custom_gaussians = {
    "X1": eltbx.xray_scattering.gaussian(
      [1], [2], 0),
    "Z1": eltbx.xray_scattering.gaussian(
      (1,2), (3,5), 0),
  }
  new_scatterers = flex.xray_scatterer()
  new_scatterers.append(xray.scatterer(
    label = "X1", scattering_type = "X1", site=(1,2,3), u=1.1, occupancy=0.5))
  new_scatterers.append(xray.scatterer(
    label = "Z1", scattering_type = "Z1", site=(4,5,6), u=9.1, occupancy=1.5))
  xs1 = xray.structure(sp, new_scatterers)
  xs1.scattering_type_registry(custom_dict=custom_gaussians)
  ##
  out = StringIO()
  xs.concatenate_inplace(other = xs1)
  xs.scattering_type_registry().show(out=out)
  expected_result = """\
Number of scattering types: 4
  Type Number    sf(0)   Gaussians
   O       1      8.00       6
   C       1      6.00       6
   Z1      1      3.00       2
   X1      1      1.00       1
  sf(0) = scattering factor at diffraction angle 0.
"""
  assert out.getvalue() == expected_result
  #
  out = sys.stdout
  sys.stdout = StringIO()
  try:
    custom_gaussians = {"C": eltbx.xray_scattering.gaussian([1],[2], 0)}
    new_scatterers = flex.xray_scatterer()
    new_scatterers.append(xray.scatterer(
      label = "C", scattering_type = "C", site=(7,8,9), u=0.1, occupancy=1.1))
    xs1 = xray.structure(sp, new_scatterers)
    xs1.scattering_type_registry(custom_dict = custom_gaussians)
    xs.concatenate_inplace(other = xs1)
    xs.scattering_type_registry().show()
  except Exception, e: pass
  assert str(e) == "Cannot concatenate: conflicting scatterers"
  sys.stdout = out
  #
  assert [(r.scattering_type, r.count, "%.1f" % r.occupancy_sum)
    for r in xs.scattering_types_counts_and_occupancy_sums()] \
      == [('C', 1, "1.0"), ('X1', 1, "0.5"), ('Z1', 1, "1.5"), ('O', 1, "1.0")]

def exercise_min_u_cart_eigenvalue():
  cs = crystal.symmetry((1, 1, 1, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  a = flex.xray_scatterer()
  assert a.size() == 0
  s1 = xray.scatterer(label = "C", u = -0.0278)
  s2 = xray.scatterer(label = "C", u = -10.0)
  s2.flags.set_use_u_iso(False)
  s3 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s4 = xray.scatterer(label = "C", u = (-91,1,1,1,1,1))
  s4.flags.set_use_u_aniso(False)
  s5 = xray.scatterer(label = "C", u = 0.1)
  s5.u_star=(1,1,1,1,1,1)
  s5.flags.set_use_u_aniso(True)
  s6 = xray.scatterer(label = "C", u = 0.1)
  s6.u_star=(1,1,1,1,1,1)
  s7 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s7.u_iso=0.1
  s8 = xray.scatterer(label = "C", u = (1,1,1,1,1,1))
  s8.u_iso=0.1
  s8.flags.set_use_u_iso(True)
  s9 = xray.scatterer(label = "C")
  s10 = xray.scatterer(label = "C")
  s10.flags.set_use_u_iso(False)
  scatterers = flex.xray_scatterer((s1,s2,s3,s4,s5,s6,s7,s8,s9,s10))
  xs = xray.structure(sp, scatterers)
  assert approx_equal(xs.min_u_cart_eigenvalue(), -0.0278)

def exercise_replace_sites():
  cs = crystal.symmetry((10, 10, 10, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer(
    [xray.scatterer("c", (-1, -1, -1))])
  xrs = xray.structure(sp, scatterers)
  #
  sites_cart = flex.vec3_double([(1,2,3)])
  xrs_ = xrs.replace_sites_cart(sites_cart)
  assert approx_equal(flex.mean(
    xrs_.sites_cart().as_double()-sites_cart.as_double()), 0)
  assert approx_equal(flex.mean(
    xrs_.sites_frac().as_double()-flex.double([0.1,0.2,0.3])), 0)
  #
  sites_frac = flex.vec3_double([(0.1,0.2,0.3)])
  xrs_ = xrs.replace_sites_frac(sites_frac)
  assert approx_equal(flex.mean(
    xrs_.sites_frac().as_double()-sites_frac.as_double()), 0)
  assert approx_equal(flex.mean(
    xrs_.sites_cart().as_double()-flex.double([1,2,3])), 0)

def exercise_add_scatterer_insert():
  cs = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 2")
  sp = crystal.special_position_settings(cs)
  xs = xray.structure(sp, scatterers=flex.xray_scatterer((
    xray.scatterer("c", site=(0.5,0,0.1)),
    xray.scatterer("o", site=(0.5,0,0.2)),
    xray.scatterer("n", site=(0.0,2,1.0)),
    xray.scatterer("p", site=(0.5,0,0.3)),
    xray.scatterer("s", site=(0.5,0,0.5)))))
  n_sites = xs.scatterers().size()
  scatterer = xray.scatterer("zn", site=(0,0,0))
  ls = [sc.label for sc in xs.scatterers()]
  ss = [str(xs.site_symmetry_table().get(i_seq).special_op())
    for i_seq in xrange(n_sites)]
  for i_seq in xrange(n_sites+1):
    xsw = xs.deep_copy_scatterers()
    xsw.add_scatterer(scatterer=scatterer, insert_at_index=i_seq)
    lsw = list(ls)
    lsw.insert(i_seq, "zn")
    assert lsw == [sc.label for sc in xsw.scatterers()]
    ssw = list(ss)
    ssw.insert(i_seq, "0,y,0")
    assert ssw == [str(xsw.site_symmetry_table().get(i_seq).special_op())
      for i_seq in xrange(n_sites+1)]

def exercise_select_on_name_or_chemical_element():
  cs = crystal.symmetry((5,7,9, 90, 120, 90), 'P2')
  xs = xray.structure(crystal.special_position_settings(cs),
                      scatterers=flex.xray_scatterer((
    xray.scatterer("C1", site=(0,0,0)),
    xray.scatterer("C2", site=(0.5,0,0)),
    xray.scatterer("C3", site=(0,0.5,0)),
    xray.scatterer("O1", site=(0,0,0.5)),
    xray.scatterer("C4", site=(0.2,0,0)),
    xray.scatterer("N1", site=(0,0.2,0)),
    )))
  sel = xs.element_selection('O', 'N')
  assert tuple(sel) == (0, 0, 0, 1, 0, 1)
  sel = xs.label_selection('C4')
  assert tuple(sel) == (0, 0, 0, 0, 1, 0)
  sel = xs.label_regex_selection("^(O|C)(1|2)$")
  assert tuple(sel) == (1, 1, 0, 1, 0, 0)

def exercise_chemical_formula():
  cs = crystal.symmetry((10,10,10, 90,90,90), 'hall: P 2 2 3')
  xs = xray.structure(crystal.special_position_settings(cs),
                      scatterers=flex.xray_scatterer((
    xray.scatterer('C1', site=(0,0,0)),
    xray.scatterer("C2", site=(0.5,0,0), occupancy=0.2),
    xray.scatterer("C2'", site=(0.5,0,0), occupancy=0.8),
    xray.scatterer("C3", site=(0,0.5,0)),
    xray.scatterer("O1", site=(0,0,0.5)),
    xray.scatterer("C4", site=(0.2,0,0)),
    xray.scatterer("N1", site=(0,0.2,0)),
    )))
  unit_cell_content = {'C':13, 'O':3, 'N':6}
  assert xs.unit_cell_content() == unit_cell_content
  assert approx_equal(xs.f_000(), 144, eps=1e-2)
  xs.set_inelastic_form_factors(0.71073, "henke")
  assert approx_equal(xs.f_000(), 144, eps=1e-2)
  assert approx_equal(xs.f_000(include_inelastic_part=True), 144.116285145)
  del unit_cell_content['C']
  assert xs.unit_cell_content(omit=set('C')) == unit_cell_content
  assert approx_equal(xs.crystal_density(), 0.47850314720502857)

def exercise_parameter_map():
  cs = crystal.symmetry((8,9,10, 85, 95, 105), "P1")
  xs = xray.structure(cs.special_position_settings())
  for i in xrange(5): xs.add_scatterer(xray.scatterer("C%i" % i))
  grad_site      = (True , False, False, True , False)
  grad_u_iso     = (False, True , True , False, True )
  grad_u_aniso   = (True , False, False, True , True )
  grad_occupancy = (False, True , True , False, False)
  grad_fp        = (True , True , False, False, False)
  grad_fdp       = (True , False, True , False, False)
  for sc, site, u_iso, u_aniso, occ, fp, fdp in zip(xs.scatterers(),
    grad_site, grad_u_iso, grad_u_aniso, grad_occupancy, grad_fp, grad_fdp):
    f = sc.flags
    f.set_grad_site(site)
    f.set_use_u_iso(u_iso)
    f.set_grad_u_iso(u_iso)
    f.set_use_u_aniso(u_aniso)
    f.set_grad_u_aniso(u_aniso)
    f.set_grad_occupancy(occ)
    f.set_grad_fp(fp)
    f.set_grad_fdp(fdp)
  m1 = xs.parameter_map()
  twins = (xray.twin_component(sgtbx.rot_mx((1,0,0,0,1,0,0,0,-1)),0.5, True),
           xray.twin_component(sgtbx.rot_mx((-1,0,0,0,-1,0,0,0,-1)), 0.2, False))
  m2 = xray.parameter_map(xs.scatterers())
  for t in twins:
    if t.grad:
      m2.add_independent_scalar()
  assert m1.n_parameters == xs.n_parameters()
  assert m2.n_parameters == xs.n_parameters()+1

  for m in (m1, m2):
    indices = m[0]
    assert indices.site == 0
    assert indices.u_iso == xray.parameter_indices.invariable
    assert indices.u_aniso == 3
    assert indices.occupancy == xray.parameter_indices.invariable
    assert indices.fp == 9
    assert indices.fdp == 10

    indices = m[1]
    assert indices.site == xray.parameter_indices.invariable
    assert indices.u_iso == 11
    assert indices.u_aniso == xray.parameter_indices.invariable
    assert indices.occupancy == 12
    assert indices.fp == 13
    assert indices.fdp == xray.parameter_indices.invariable

    for i, indices in enumerate(m):
      assert indices.site == m[i].site
      assert indices.u_iso == m[i].u_iso
      assert indices.u_aniso == m[i].u_aniso
      assert indices.occupancy == m[i].occupancy
      assert indices.fp == m[i].fp
      assert indices.fdp == m[i].fdp


def exercise_xray_structure_as_py_code():
  import itertools
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry((2, 2, 3, 90, 90, 80), "hall: P 2z"),
    scatterers=flex.xray_scatterer((
      xray.scatterer('C1', site=(0.5, 0.5, 0.5), u=0.1),
      xray.scatterer('O1', site=(0.1, 0.2, 0.3), u=(0.1, 0.2, 0.3,
                                                    0.4, 0.5, 0.6)),
      xray.scatterer('Fe', site=(-0.8, 0.2, 0), u=0.2,
                     scattering_type="Fe3+")
      )))
  pc = xs.as_py_code(indent="V")
  assert not show_diff(pc, """\
Vxray.structure(
V  crystal_symmetry=crystal.symmetry(
V    unit_cell=(2, 2, 3, 90, 90, 80),
V    space_group_symbol="P 1 1 2"),
V  scatterers=flex.xray_scatterer([
V    xray.scatterer( #0
V      label="C1",
V      site=(0.500000, 0.500000, 0.500000),
V      u=0.100000),
V    xray.scatterer( #1
V      label="O1",
V      site=(0.100000, 0.200000, 0.300000),
V      u=(0.100000, 0.200000, 0.300000,
V         0.400000, 0.500000, 0.600000)),
V    xray.scatterer( #2
V      label="Fe",
V      scattering_type="Fe3+",
V      site=(-0.800000, 0.200000, 0.000000),
V      u=0.200000)]))""")
  pc = xs.as_py_code()
  xs1 = eval(pc)
  assert xs.crystal_symmetry().is_similar_symmetry(
    xs1.crystal_symmetry(),
    relative_length_tolerance=0,
    absolute_angle_tolerance=0)
  for sc, sc1 in itertools.izip(xs.scatterers(), xs1.scatterers()):
    assert sc.flags.bits == sc1.flags.bits
    assert sc.site  == sc1.site
    if sc.flags.use_u_iso():
      assert sc.u_iso == sc1.u_iso
    if sc.flags.use_u_aniso():
      assert sc.u_star == sc1.u_star
    assert sc.occupancy == sc1.occupancy
    assert sc.fp  == sc1.fp
    assert sc.fdp == sc1.fdp

def exercise_delta_sites_cart_measure():
  rnd_delta = 1e-3
  for hall_symbol, continuously_shift in [
    ('P 2yb', lambda x,y,z: (x    , y+0.7, z)    ),
    ('P -2x', lambda x,y,z: (x    , y+0.9, z-0.7)),
    ('P 1'  , lambda x,y,z: (x+0.1, y+0.2, z+0.9)),
    ]:
    xs0 = random_structure.xray_structure(
      sgtbx.space_group_info('hall: %s' % hall_symbol),
      use_u_iso=False,
      use_u_aniso=True,
      n_scatterers=5,
      elements="random")
    xs1 = xs0.deep_copy_scatterers()
    delta = flex.vec3_double()
    for sc in xs1.scatterers():
      x, y, z = continuously_shift(*sc.site)
      dx, dy, dz = [ random.uniform(-rnd_delta, rnd_delta)
                     for i in xrange(3) ]
      delta.append((dx, dy, dz))
      sc.site = (x + dx, y + dy, z + dz)
    delta = flex.vec3_double([ xs0.unit_cell().orthogonalize(d)
                               for d in delta ])
    diff_ref = flex.max_absolute(delta.as_double())
    mscd = xray.meaningful_site_cart_differences(xs1=xs1, xs2=xs0)
    assert approx_equal(mscd.max_absolute(), diff_ref, eps=rnd_delta)

def run():
  exercise_delta_sites_cart_measure()
  exercise_xray_structure_as_py_code()
  exercise_parameter_map()
  exercise_chemical_formula()
  exercise_select_on_name_or_chemical_element()
  exercise_add_scatterer_insert()
  exercise_replace_sites()
  exercise_min_u_cart_eigenvalue()
  exercise_set_occupancies()
  exercise_closest_distances()
  exercise_concatenate_inplace()
  exercise_scatterer()
  exercise_anomalous_scatterer_group()
  exercise_structure()
  exercise_u_base()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
