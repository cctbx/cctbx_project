"""
Tests for mmtbx.hydrogens.xh_mapper.

The model is ARG 3 (truncated at CB, so it has the cap HCG at occupancy 0)
and SER 4 with OXT, plus a water, hydrogenated by mmtbx.reduce2
qmready=true capping=true. SER 4 OG and HG are split into altlocs A and B
(occupancy 0.5, B displaced by 0.3 A) and the backbone amide H of SER 4 is
renamed D (element D).
"""
import math
import iotbx.pdb
import iotbx.phil
import mmtbx.model
import mmtbx.hydrogens
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from mmtbx.hydrogens import xh_mapper

pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A   3      27.911  23.777  26.214  1.00  5.99           N
ATOM      2  CA  ARG A   3      27.558  25.104  25.723  1.00  6.90           C
ATOM      3  C   ARG A   3      26.047  25.238  25.567  1.00  6.28           C
ATOM      4  O   ARG A   3      25.450  26.196  26.055  1.00  7.24           O
ATOM      5  CB  ARG A   3      28.246  25.375  24.384  1.00  6.79           C
ATOM      6  H1  ARG A   3      27.335  23.552  27.050  1.00  5.99           H
ATOM      7  H2  ARG A   3      27.731  23.073  25.470  1.00  5.99           H
ATOM      8  H3  ARG A   3      28.918  23.759  26.472  1.00  5.99           H
ATOM      9  HA  ARG A   3      27.896  25.854  26.438  1.00  6.90           H
ATOM     10  HB2 ARG A   3      29.327  25.310  24.511  1.00  6.79           H
ATOM     11  HB3 ARG A   3      27.922  24.635  23.652  1.00  6.79           H
ATOM     12  HCG ARG A   3      27.981  26.373  24.034  0.00  6.79           H
ATOM     13  N   SER A   4      25.443  24.270  24.884  1.00  5.99           N
ATOM     14  CA  SER A   4      24.000  24.258  24.649  1.00  6.90           C
ATOM     15  C   SER A   4      23.547  25.461  23.816  1.00  6.28           C
ATOM     16  O   SER A   4      23.472  25.372  22.590  1.00  7.24           O
ATOM     17  CB  SER A   4      23.234  24.186  25.975  1.00  6.79           C
ATOM     18  OG ASER A   4      21.834  24.229  25.757  0.50  8.89           O
ATOM     18  OG BSER A   4      22.134  24.229  25.757  0.50  8.89           O
ATOM     19  OXT SER A   4      23.253  26.536  24.382  1.00  6.28           O
ATOM     20  D   SER A   4      25.929  23.471  24.476  1.00  5.99           D
ATOM     21  HA  SER A   4      23.749  23.365  24.076  1.00  6.90           H
ATOM     22  HB2 SER A   4      23.485  23.253  26.479  1.00  6.79           H
ATOM     23  HB3 SER A   4      23.522  25.034  26.596  1.00  6.79           H
ATOM     24  HG ASER A   4      21.570  25.122  25.452  0.50  8.89           H
ATOM     24  HG BSER A   4      21.870  25.122  25.452  0.50  8.89           H
TER
HETATM   25  O   HOH S   1      21.834  24.229  29.757  1.00  8.00           O
HETATM   26  H1  HOH S   1      22.814  24.229  29.757  1.00  8.00           H
HETATM   27  H2  HOH S   1      21.598  25.180  29.757  1.00  8.00           H
END
"""

pdb_str_no_h = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A   3      27.911  23.777  26.214  1.00  5.99           N
ATOM      2  CA  ARG A   3      27.558  25.104  25.723  1.00  6.90           C
ATOM      3  C   ARG A   3      26.047  25.238  25.567  1.00  6.28           C
ATOM      4  O   ARG A   3      25.450  26.196  26.055  1.00  7.24           O
ATOM      5  CB  ARG A   3      28.246  25.375  24.384  1.00  6.79           C
ATOM     13  N   SER A   4      25.443  24.270  24.884  1.00  5.99           N
ATOM     14  CA  SER A   4      24.000  24.258  24.649  1.00  6.90           C
ATOM     15  C   SER A   4      23.547  25.461  23.816  1.00  6.28           C
ATOM     16  O   SER A   4      23.472  25.372  22.590  1.00  7.24           O
ATOM     17  CB  SER A   4      23.234  24.186  25.975  1.00  6.79           C
ATOM     18  OG ASER A   4      21.834  24.229  25.757  0.50  8.89           O
ATOM     18  OG BSER A   4      22.134  24.229  25.757  0.50  8.89           O
ATOM     19  OXT SER A   4      23.253  26.536  24.382  1.00  6.28           O
TER
HETATM   25  O   HOH S   1      21.834  24.229  29.757  1.00  8.00           O
END
"""

edit_water_oh = """\
geometry_restraints.edits {
  bond {
    action = add *change delete
    atom_selection_1 = name O and resname HOH
    atom_selection_2 = name H1 and resname HOH
    distance_ideal = 1.5
    sigma = 0.02
  }
}
"""

def get_model(pdb_str, scattering_table, use_neutron_distances, edits=None):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  pip = mmtbx.model.manager.get_default_pdb_interpretation_params()
  if edits is not None:
    scope = mmtbx.model.manager.get_default_pdb_interpretation_scope()
    pip = scope.fetch(iotbx.phil.parse(edits)).extract()
  pip.pdb_interpretation.use_neutron_distances = use_neutron_distances
  model.process(make_restraints=True, pdb_interpretation_params=pip)
  if scattering_table is not None:
    model.setup_scattering_dictionaries(scattering_table=scattering_table)
  return model

def xh_lengths(model, sites_cart):
  """{i_h: (i_x, ideal X-H length)} as the model's restraints give them."""
  xrs = model.get_xray_structure().deep_copy_scatterers()
  xrs.set_sites_cart(sites_cart)
  table = mmtbx.model.xh_connectivity_table(
    geometry=model.get_restraints_manager(), xray_structure=xrs).table
  return dict((e[1], (e[0], e[3])) for e in table)

def dist(sites, i, j):
  return flex.vec3_double([sites[i]]).__sub__(
    flex.vec3_double([sites[j]])).norms()[0]

def unit(sites, i, j):
  v = flex.vec3_double([sites[i]]) - flex.vec3_double([sites[j]])
  return tuple(v[0][k] / v.norms()[0] for k in range(3))

def test_lengths():
  m = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  mp = xh_mapper.mapper(m)
  assert not mp.is_identity
  sites = m.get_sites_cart()
  x = mp.shrink(sites)
  n = mp.expand(sites)
  nuclear = xh_lengths(m, sites)
  mx = get_model(pdb_str, "n_gaussian", use_neutron_distances=False)
  xray = xh_lengths(mx, sites)
  assert mp.selection().count(True) == len(nuclear) == mp.i_h.size()
  for i_h, (i_x, d_n) in nuclear.items():
    d_x = xray[i_h][1]
    assert d_n > d_x
    assert approx_equal(dist(x, i_h, i_x), d_x, 1e-6)
    assert approx_equal(dist(n, i_h, i_x), d_n, 1e-6)
    assert approx_equal(unit(x, i_h, i_x), unit(sites, i_h, i_x), 1e-6)
  heavy = ~m.get_hd_selection()
  assert approx_equal(x.select(heavy), sites.select(heavy), 1e-12)
  assert approx_equal(mp.shrink(x), x, 1e-9)
  assert approx_equal(mp.expand(n), n, 1e-9)
  assert approx_equal(mp.expand(mp.shrink(sites)), n, 1e-9)
  assert sites.rms_difference(m.get_sites_cart()) == 0
  names = dict((a.i_seq, a.name.strip()) for a in m.get_hierarchy().atoms())
  mapped = set(names[i] for i in mp.i_h)
  assert set(["HCG", "H1", "H2", "D"]) <= mapped, mapped

def test_identity():
  for table in ["neutron", "electron"]:
    m = get_model(pdb_str, table, use_neutron_distances=True)
    mp = xh_mapper.mapper(m)
    assert mp.is_identity and mp.selection().count(True) == 0
    sites = m.get_sites_cart()
    for f in [mp.shrink, mp.expand]:
      out = f(sites)
      assert out is not sites and out.rms_difference(sites) == 0

def test_entry_from_either_convention():
  mn = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  mx = get_model(pdb_str, "n_gaussian", use_neutron_distances=False)
  mx.set_hydrogen_bond_length(use_neutron_distances=False)
  a = xh_mapper.mapper(mn).expand(mn.get_sites_cart())
  b = xh_mapper.mapper(mx).expand(mx.get_sites_cart())
  assert approx_equal(a, b, 1e-6)

def test_non_default_params():
  """
  The X-ray counterpart of the model is processed with the model's own
  settings, not with defaults: a restraint edit setting the water O-H1 length
  to 1.5 must be seen in both the nuclear and the X-ray length.
  """
  for neutron in [True, False]:
    m = get_model(pdb_str, "n_gaussian", neutron, edits=edit_water_oh)
    mp = xh_mapper.mapper(m)
    water_h1 = [a.i_seq for a in m.get_hierarchy().atoms()
      if a.name.strip() == "H1" and a.parent().resname == "HOH"][0]
    k = list(mp.i_h).index(water_h1)
    assert approx_equal(mp.d_xray[k], 1.5)
    assert approx_equal(mp.d_nuclear[k], 1.5)
    other = get_model(pdb_str, "n_gaussian", not neutron, edits=edit_water_oh)
    mo = xh_mapper.mapper(other)
    assert approx_equal(mo.d_xray, mp.d_xray)
    assert approx_equal(mo.d_nuclear, mp.d_nuclear)

def test_no_h():
  m = get_model(pdb_str_no_h, "n_gaussian", use_neutron_distances=True)
  mp = xh_mapper.mapper(m)
  assert not mp.is_identity and mp.selection().count(True) == 0
  sites = m.get_sites_cart()
  assert mp.shrink(sites).rms_difference(sites) == 0

def test_no_scattering_table():
  m = get_model(pdb_str, None, use_neutron_distances=True)
  # process() already registers a default table, so make the model report
  # that none is set.
  m.get_scattering_table = lambda: None
  try:
    xh_mapper.mapper(m)
  except Sorry:
    pass
  else:
    raise AssertionError("Sorry expected")

def smooth_target(sites):
  """
  Sum over all pairs i<j of w_ij * exp(-|ri-rj|^2/4), the weights coming from
  a seeded random.Random(0). Returns the value and its analytic gradient.
  """
  import random
  rng = random.Random(0)
  n = sites.size()
  f = 0
  g = [[0., 0., 0.] for i in range(n)]
  for i in range(n):
    for j in range(i + 1, n):
      w = 0.5 + rng.random()
      d = [sites[i][k] - sites[j][k] for k in range(3)]
      e = w * math.exp(-(d[0]**2 + d[1]**2 + d[2]**2) / 4)
      f += e
      for k in range(3):
        g[i][k] -= 0.5 * e * d[k]
        g[j][k] += 0.5 * e * d[k]
  return f, flex.vec3_double([tuple(x) for x in g])

def fd(func, sites, i, k, eps):
  """Central difference of func(sites) w.r.t. coordinate k of atom i."""
  values = []
  for sign in [1, -1]:
    s = sites.deep_copy()
    x = list(s[i])
    x[k] += sign * eps
    s[i] = tuple(x)
    values.append(func(s))
  return (values[0] - values[1]) / (2 * eps)

def test_smooth_target_self_check():
  m = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  sites = m.get_sites_cart()
  f, g = smooth_target(sites)
  for k in range(3):
    gp = fd(lambda s: smooth_target(s)[0], sites, 5, k, 1.e-5)
    assert approx_equal(g[5][k], gp, 1.e-6), (k, g[5][k], gp)

def test_gradients_fd():
  m = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  mp = xh_mapper.mapper(m)
  sites = m.get_sites_cart()
  f0, gx = smooth_target(mp.shrink(sites))
  g = mp.gradients(sites, gx)
  eps = 1.e-5
  for i in range(sites.size()):
    for k in range(3):
      gp = fd(lambda s: smooth_target(mp.shrink(s))[0], sites, i, k, eps)
      assert approx_equal(g[i][k], gp, 1.e-6), (i, k, g[i][k], gp)
  assert g.size() == sites.size()

def test_gradients_fd_fmodel():
  import random
  import mmtbx.f_model
  from cctbx import xray
  m = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  mp = xh_mapper.mapper(m)
  f_obs = abs(m.get_xray_structure().structure_factors(
    d_min=1.5, algorithm="direct").f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  xrs = m.get_xray_structure().deep_copy_scatterers()
  random.seed(0)
  flex.set_random_seed(0)
  xrs.shake_sites_in_place(mean_distance=0.3)
  # The shaken sites are the nuclear evaluation point.
  nuclear = xrs.sites_cart().deep_copy()
  # The default fft algorithm makes the gradients only approximate.
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=xrs, target_name="ls_wunit_k1",
    sf_and_grads_accuracy_params=params)
  def target(sites):
    fmodel.xray_structure.set_sites_cart(mp.shrink(sites))
    fmodel.update_xray_structure(update_f_calc=True)
    return fmodel.target_w()
  target(nuclear)
  xray.set_scatterer_grad_flags(
    scatterers=fmodel.xray_structure.scatterers(), site=True)
  dg = flex.vec3_double(fmodel.target_functor()(
    compute_gradients=True).gradients_wrt_atomic_parameters(
      site=True).packed())
  g = mp.gradients(nuclear, dg)
  floor = 1.e-3 * (flex.sum(g.dot(g)) / g.size())**0.5
  atoms = sorted(set(mp.i_h).union(set(mp.i_x)))
  for i in atoms:
    for k in range(3):
      gp = fd(target, nuclear, i, k, 1.e-4)
      assert abs(g[i][k] - gp) <= 1.e-3 * max(abs(gp), floor), (
        i, k, g[i][k], gp)

def test_gradients_identity():
  m = get_model(pdb_str, "neutron", use_neutron_distances=True)
  mp = xh_mapper.mapper(m)
  g = flex.vec3_double(m.size(), (1, 2, 3))
  out = mp.gradients(m.get_sites_cart(), g)
  assert out is not g and out.rms_difference(g) == 0

def test_fit_rotatable2_keeps_nuclear():
  import mmtbx.f_model
  from scitbx.matrix import rotate_point_around_axis
  m = get_model(pdb_str, "n_gaussian", use_neutron_distances=True)
  assert m.riding_h_manager is None
  mp = xh_mapper.mapper(m)
  sites = m.get_sites_cart()
  lengths_before = [dist(sites, i, j) for i, j in zip(mp.i_h, mp.i_x)]
  f_obs = abs(m.get_xray_structure().structure_factors(
    d_min=1.5, algorithm="direct").f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  # fmodel is built on a copy with the SER HG rotated 120 degrees about
  # CB-OG, so the fit has something to do.
  rotated = sites.deep_copy()
  atoms = m.get_hierarchy().atoms()
  n_rotated = 0
  for a in atoms:
    if(a.parent().resname == "SER" and a.name.strip() == "HG"):
      rg_atoms = a.parent().parent().atoms()
      cb = [b.i_seq for b in rg_atoms if b.name.strip() == "CB"]
      og = [b.i_seq for b in a.parent().atoms() if b.name.strip() == "OG"]
      assert len(cb) == 1 and len(og) == 1
      rotated[a.i_seq] = rotate_point_around_axis(
        axis_point_1=sites[cb[0]], axis_point_2=sites[og[0]],
        point=sites[a.i_seq], angle=120, deg=True)
      n_rotated += 1
  assert n_rotated > 0
  xrs = m.get_xray_structure().deep_copy_scatterers()
  xrs.set_sites_cart(mp.shrink(rotated))
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=xrs, target_name="ls_wunit_k1",
    sf_and_grads_accuracy_params=params)
  mmtbx.hydrogens.fit_rotatable2(model=m, fmodel=fmodel, reset_h=False,
    mapper=mp)
  sites = m.get_sites_cart()
  for k, (i, j) in enumerate(zip(mp.i_h, mp.i_x)):
    assert approx_equal(dist(sites, i, j), lengths_before[k], 1.e-6)
  assert m.riding_h_manager is None
  assert approx_equal(list(fmodel.xray_structure.sites_cart().as_double()),
    list(mp.shrink(sites).as_double()), 1.e-9)

def run():
  test_lengths()
  test_identity()
  test_entry_from_either_convention()
  test_non_default_params()
  test_no_h()
  test_no_scattering_table()
  test_smooth_target_self_check()
  test_gradients_fd()
  test_gradients_fd_fmodel()
  test_gradients_identity()
  test_fit_rotatable2_keeps_nuclear()
  print("OK")

if __name__ == "__main__":
  run()
