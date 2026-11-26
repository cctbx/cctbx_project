from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx import matrix
import iotbx.phil
import mmtbx.refinement.refinement_flags
import mmtbx.solvent.ordered_solvent as ordered_solvent
from mmtbx.regression import tst_ordered_solvent_altlocs as t

"""
Options of mmtbx.solvent.ordered_solvent: parameters, keep_input_water,
ADP type of new water (new_solvent, anisotropic_d_min).
"""

def exercise_params():
  p = ordered_solvent.master_params().extract()
  assert p.new_solvent == "Auto", p.new_solvent
  assert approx_equal(p.anisotropic_d_min, 1.2)
  assert p.keep_input_water is False
  assert p.coordination_restraints.enabled is True
  assert approx_equal(p.coordination_restraints.sigma, 0.5)
  assert approx_equal(p.dist_min_altloc, 0.5)
  # parameters that are gone
  master = ordered_solvent.master_params()
  for s in ["include_altlocs=True", "correct_drifted_waters=True",
            "output_chain_id=S"]:
    w, unused = master.fetch(sources=[iotbx.phil.parse(s)],
      track_unused_definitions=True)
    assert len(unused) == 1, s
  w, unused = master.fetch(sources=[iotbx.phil.parse(
    "keep_input_water=True\nnew_solvent=anisotropic")],
    track_unused_definitions=True)
  assert len(unused) == 0
  # model.add_solvent puts water into the chain it is asked to
  m = t.get_model(t.ser_lines())
  m.add_solvent(sites_frac=flex.vec3_double([(0.8, 0.8, 0.1)]), chain_id="W")
  assert [c.id for c in m.get_hierarchy().chains()] == ["A", "W"]
  m.add_solvent(sites_frac=flex.vec3_double([(0.7, 0.8, 0.1)]))
  assert [c.id for c in m.get_hierarchy().chains()] == ["A", "W", "S"]
  # ensemble refinement has its own
  import mmtbx.solvent.ensemble_ordered_solvent as eos
  assert eos.master_params().extract().output_chain_id == "S"
  assert "output_chain_id" in ordered_solvent.output_params_str

W_O = t.O + matrix.col((0, 2.8, 0))
W_N = t.N + matrix.col((0, 0, -2.8))

def exercise_keep_input_water():
  """
  keep_input_water=True: water of input model stays as it is whatever it is;
  new water is still added, but not next to it.
  """
  w_far   = t.FREE                                # nothing there, no partner
  w_clash = t.CB + matrix.col((0, 0, 1.0))        # 1 A from CB
  w_b     = W_O                                   # real, but B=150
  w_new   = W_N                                   # to be found
  w_alt   = W_O + matrix.col((1.5, 0, 0))         # in data, next to w_b
  lines = t.ser_lines() + [
    t.hoh(1, w_far, b=20), t.hoh(2, w_clash, occ=0.3, b=30),
    t.hoh(3, w_b, b=150)]
  lines_true = t.ser_lines() + [t.hoh(3, w_b, occ=0.5), t.hoh(4, w_new),
    t.hoh(5, w_alt, occ=0.5)]
  for keep in [True, False]:
    params = ordered_solvent.master_params().extract()
    params.add_hydrogens.d_min = None
    params.keep_input_water = keep
    model = t.get_model(lines)
    ws = ordered_solvent.water_state(model=model, keep_input_water=keep)
    assert len(ws.input_water_keys) == (3 if keep else 0)
    o = t.run_manager(lines_true, lines, d_min=1.0, params=params,
      water_state=ws)
    m = o.model
    if(keep):
      for xyz, occ, b in [(w_far, 1, 20), (w_clash, 0.3, 30), (w_b, 1, 150)]:
        a = t.find_water(m, xyz, tolerance=1.e-3)
        assert approx_equal(a.occ, occ) and approx_equal(a.b, b, 1.e-2)
        assert a.parent().altloc == ""
      # no new water next to the one that is kept
      for a in t.waters(m):
        d = (matrix.col(a.xyz)-w_b).length()
        assert d < 1.e-3 or d > 2.2, d
      assert ws.fixed_selection(m).count(True) == 3
      # second update: still there
      o2 = ordered_solvent.manager(fmodel=o.fmodel, fmodels=None, model=m,
        is_neutron_scat_table=False, params=params, water_state=ws,
        log=null_out())
      for xyz in [w_far, w_clash, w_b]:
        t.find_water(o2.model, xyz, tolerance=1.e-3)
    else:
      for xyz in [w_far, w_clash]:
        assert len([a for a in t.waters(m)
          if (matrix.col(a.xyz)-xyz).length() < 0.5]) == 0
      for a in t.waters(m): assert a.b < 90, a.b

def aniso_ser_lines():
  """ SER with anisotropic ADP """
  result = []
  for l in t.ser_lines():
    result.append(l)
    u = int(round(adptbx.b_as_u(float(l[60:66]))*1.e4))
    result.append("ANISOU"+l[6:28]+"%7d%7d%7d%7d%7d%7d"%(u,u,u,0,0,0)+
      "      "+l[76:78])
  return result

def exercise_resolve():
  p = ordered_solvent.master_params().extract()
  m_iso = t.get_model(t.ser_lines())
  m_aniso = t.get_model(aniso_ser_lines())
  assert m_aniso.get_xray_structure().use_u_aniso().count(True) == 7
  r = ordered_solvent.resolve_new_solvent
  assert p.new_solvent == "Auto"
  assert r(params=p, d_min=1.0, model=m_aniso) == "anisotropic"
  assert r(params=p, d_min=1.5, model=m_aniso) == "isotropic"
  assert r(params=p, d_min=1.0, model=m_iso)   == "isotropic"
  assert r(params=p, d_min=1.5, model=m_iso)   == "isotropic"
  assert r(params=p, d_min=1.2, model=m_aniso) == "isotropic" # not better
  p.anisotropic_d_min = 2.0
  assert r(params=p, d_min=1.5, model=m_aniso) == "anisotropic"
  for choice in ["isotropic", "anisotropic"]:
    p.new_solvent = choice
    for m in [m_iso, m_aniso]:
      for d_min in [1.0, 3.0]:
        assert r(params=p, d_min=d_min, model=m) == choice
  # water does not count
  p.new_solvent = "Auto"
  p.anisotropic_d_min = 1.2
  lines = t.ser_lines() + [t.hoh(i, t.FREE+matrix.col((i*3,0,0)))
    for i in range(10)]
  m = t.get_model(lines)
  m.get_xray_structure().convert_to_anisotropic(
    selection = m.solvent_selection())
  assert r(params=p, d_min=1.0, model=m) == "isotropic"

def exercise_adp_type():
  """
  ADP of new water as asked for; refined in water update as isotropic, then
  converted; water of the model that is anisotropic stays so.
  """
  aniso_w = "ANISOU  101  O   HOH A 101     1500   2500   3500    100    200" \
            "    300       O"
  lines = aniso_ser_lines() + [t.hoh(1, W_O, b=20), aniso_w]
  lines_true = lines + [t.hoh(2, W_N, b=12)]
  for choice, d_min, expected in [
        ("anisotropic", 1.0, True),
        ("Auto",        1.0, True),
        ("isotropic",   1.0, False),
        ("Auto",        1.5, False),
        ("anisotropic", 1.5, True)]:
    params = ordered_solvent.master_params().extract()
    params.add_hydrogens.d_min = None
    params.new_solvent = choice
    model_true = t.get_model(lines_true)
    model = t.get_model(lines)
    for m in [model_true, model]:
      m.setup_scattering_dictionaries(scattering_table="n_gaussian")
    model.process(make_restraints=True)
    n = model.size()
    aniso = model.get_xray_structure().use_u_aniso()
    model.set_refinement_flags(mmtbx.refinement.refinement_flags.manager(
      individual_sites     = True,
      individual_adp       = True,
      sites_individual     = flex.bool(n, True),
      adp_individual_iso   = ~aniso,
      adp_individual_aniso = aniso))
    fmodel = t.get_fmodel(model_true, model, d_min)
    o = t.manager_no_map_filter(fmodel=fmodel, fmodels=None, model=model,
      is_neutron_scat_table=False, params=params, log=null_out())
    m = o.model
    xrs = m.get_xray_structure()
    sc = xrs.scatterers()
    a_old = t.find_water(m, W_O)
    assert sc[a_old.i_seq].flags.use_u_aniso()
    assert approx_equal(a_old.uij, (0.15,0.25,0.35,0.01,0.02,0.03), 1.e-4)
    a_new = t.find_water(m, W_N)
    assert sc[a_new.i_seq].flags.use_u_aniso() == expected, (choice, d_min)
    b = adptbx.u_as_b(sc[a_new.i_seq].u_iso_or_equiv(xrs.unit_cell()))
    assert 5 < b < 30, b
    rf = m.refinement_flags
    aniso = xrs.use_u_aniso()
    assert rf.adp_individual_aniso.all_eq(aniso)
    assert rf.adp_individual_iso.all_eq(~aniso)
    if(expected):
      assert a_new.uij != (-1,-1,-1,-1,-1,-1)
      # as written out
      pdb = m.model_as_pdb()
      assert len([l for l in pdb.splitlines()
        if l.startswith("ANISOU") and l[17:20] == "HOH"]) >= 2

def exercise_b_at_limit():
  """
  Water with B at the limit b_iso_min is kept: B refined with that lower
  bound comes back as the limit give or take rounding errors (U to B
  conversion), for example -4.e-15 for b_iso_min=0. This is what a strong
  peak (an ion picked as water) gets.
  """
  lines = t.ser_lines() + [t.hoh(1, W_O, b=20), t.hoh(2, W_N, b=20)]
  for b_min, b_max, b_value, n_expected in [
      (0, 80, -4.4e-15, 2), (0, 80, 0., 2), (0, 80, 1.e-9, 2),
      (0, 80, -0.01, 1),
      (1, 80, 1-1.e-12, 2), (1, 80, 0.99, 1),
      (0, 80, 80+1.e-12, 2), (0, 80, 80.01, 1)]:
    params = ordered_solvent.master_params().extract()
    params.add_hydrogens.d_min = None
    params.mode = "filter_only" # nothing added or refined
    params.b_iso_min = b_min
    params.b_iso_max = b_max
    o = t.run_manager(lines, lines, d_min=1.0, params=params,
      map_filter=False)
    assert len(t.waters(o.model)) == 2
    a = t.find_water(o.model, W_N, tolerance=1.e-2)
    xrs = o.model.get_xray_structure()
    u = xrs.extract_u_iso_or_u_equiv()
    u[a.i_seq] = adptbx.b_as_u(b_value)
    xrs.set_u_iso(values = u)
    o.model.set_xray_structure(xrs)
    o._filter_q_b()
    assert len(t.waters(o.model)) == n_expected, (b_value,
      len(t.waters(o.model)))

def run():
  exercise_params()
  exercise_resolve()
  exercise_keep_input_water()
  exercise_adp_type()
  exercise_b_at_limit()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
