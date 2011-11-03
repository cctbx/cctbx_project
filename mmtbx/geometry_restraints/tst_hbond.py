
import libtbx.load_env
from libtbx.test_utils import approx_equal
import os
import sys

def simple_pdb () :
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   ALA A   1       4.288  62.025  -3.678  1.00 27.96           N
ATOM      2  CA  ALA A   1       3.119  62.736  -4.184  1.00 28.04           C
ATOM      3  C   ALA A   1       2.683  63.793  -3.199  1.00 28.16           C
ATOM      4  O   ALA A   1       2.311  64.910  -3.605  1.00 28.25           O
ATOM      5  CB  ALA A   1       1.962  61.780  -4.504  1.00 27.64           C
ATOM      6  H   ALA A   1       4.150  61.050  -3.504  1.00 27.96           H
ATOM     25  N   ALA A   5       2.586  67.940  -3.614  1.00 29.64           N
ATOM     26  CA  ALA A   5       1.415  68.721  -3.231  1.00 30.71           C
ATOM     27  C   ALA A   5       1.799  69.877  -2.302  1.00 29.77           C
ATOM     28  O   ALA A   5       1.213  70.964  -2.365  1.00 30.52           O
ATOM     29  CB  ALA A   5       0.345  67.845  -2.558  1.00 30.05           C
ATOM     30  H   ALA A   5       2.547  66.981  -3.333  1.00 29.64           H
""")
  return pdb_in

def exercise_simple () :
  from mmtbx.geometry_restraints import hbond
  import cctbx.geometry_restraints
  from scitbx.array_family import flex
  import iotbx.pdb
  pdb_in = simple_pdb()
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  sites_cart = hierarchy.atoms().extract_xyz()
  cache = hierarchy.atom_selection_cache()
  i_seq_1 = cache.selection("resseq 1 and name O").iselection()
  i_seq_2 = cache.selection("resseq 5 and name H").iselection()
  i_seq_3 = cache.selection("resseq 5 and name N").iselection()
  assert (i_seq_1.size() == i_seq_2.size() == i_seq_3.size() == 1)
  build_proxies = hbond.build_simple_hbond_proxies()
  # XXX note: actual bond length is approximately 2.1021
  build_proxies.add_proxy(
    i_seqs=[i_seq_1[0], i_seq_2[0]],
    distance_ideal=1.975,
    distance_cut=-1,
    weight=1/(0.05**2))
  cctbx_bond = cctbx.geometry_restraints.bond(
    sites=[sites_cart[i_seq_1[0]], sites_cart[i_seq_2[0]]],
    distance_ideal=1.975,
    weight=1/(0.05**2))
  grads = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads)
  assert approx_equal(residual, cctbx_bond.residual())
  assert approx_equal(residual, 6.45926322152, eps=0.00001)
  # as before, but with N-O restrained
  build_proxies = hbond.build_simple_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[i_seq_1[0], i_seq_3[0]],
    distance_ideal=2.9,
    distance_cut=3.5,
    weight=0.5/(0.05**2))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads)
  cctbx_bond = cctbx.geometry_restraints.bond(
    sites=[sites_cart[i_seq_1[0]], sites_cart[i_seq_3[0]]],
    distance_ideal=2.9,
    weight=0.5/(0.05**2))
  assert approx_equal(residual, cctbx_bond.residual())
  simple_bonds = hbond.get_simple_bonds(build_proxies.proxies)
  assert (simple_bonds.size() == 1)
  assert (list(simple_bonds[0]) == [3,6])
  # test proxy_select
  build_proxies = hbond.build_simple_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[0,1],
    distance_ideal=2.9,
    distance_cut=3.5,
    weight=100)
  build_proxies.add_proxy(
    i_seqs=[4,5],
    distance_ideal=2.9,
    distance_cut=3.5,
    weight=100)
  build_proxies.add_proxy(
    i_seqs=[8,9],
    distance_ideal=2.9,
    distance_cut=3.5,
    weight=100)
  proxies = build_proxies.proxies
  isel = flex.size_t([0,1,2,3,7,8,9])
  selection = proxies.proxy_select(n_seq=10, iselection=isel)
  assert (selection.size() == 2)
  proxy = selection[1]
  #print proxy.i_seqs

def compare_analytical_and_fd (proxies) :
  from mmtbx.geometry_restraints import hbond
  from scitbx.array_family import flex
  distances = []
  for i in range(250) :
    distances.append(2.0 + (i * 0.01))
  for dist in distances :
    sites_cart = flex.vec3_double([(0.0,dist,0.0), (0.0,0.0,0.0),
      (0.5774,-1.0,0.0)])
    g1 = flex.vec3_double(3, (0.0,0.0,0.0))
    r1 = hbond.target_and_gradients(
      proxies=proxies,
      sites_cart=sites_cart,
      gradient_array=g1,
      use_finite_differences=True)
    g2 = flex.vec3_double(3, (0.0,0.0,0.0))
    r2 = hbond.target_and_gradients(
      proxies=proxies,
      sites_cart=sites_cart,
      gradient_array=g2,
      use_finite_differences=False)
    assert (r2 == r1)
    if (dist == 3.5) or (dist == 3.55) :
      continue
    for (dx1,dy1,dz1), (dx2,dy2,dz2) in zip(g1, g2) :
      if (dist >= 3.5) and (dist <= 3.55) :
        assert approx_equal(dx1, dx2, eps=0.1)
        assert approx_equal(dy1, dy2, eps=0.1)
        assert approx_equal(dz1, dz2, eps=0.1)
      else :
        assert approx_equal(dx1, dx2, eps=0.00001)
        assert approx_equal(dy1, dy2, eps=0.00001)
        assert approx_equal(dz1, dz2, eps=0.00001)

def exercise_lennard_jones () :
  from mmtbx.geometry_restraints import hbond
  import cctbx.geometry_restraints
  from scitbx.array_family import flex
  build_proxies = hbond.build_lennard_jones_proxies()
  build_proxies.add_proxy(
    i_seqs=[0,1],
    distance_ideal=2.9,
    distance_cut=3.5)
  sites_cart = flex.vec3_double([(0.0,2.9,0.0), (0.0,0.0,0.0),
    (0.5774,-1.0,0.0)])
  grads = flex.vec3_double(3, (0,0,0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads,
    use_finite_differences=False)
  assert approx_equal(residual, -59.259259, eps=0.0001)
  for (dx,dy,dz) in grads :
    assert approx_equal(dx, 0, eps=0.000000001)
    assert approx_equal(dy, 0, eps=0.000000001)
    assert approx_equal(dz, 0, eps=0.000000001)
  sites_cart[0] = (0.0, 3.0, 0.0)
  grads = flex.vec3_double(3, (0,0,0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads,
    use_finite_differences=False)
  assert approx_equal(residual, -58.5286436, eps=0.0001)
  compare_analytical_and_fd(build_proxies.proxies)
  grads = flex.vec3_double(3, (0,0,0))
  residual = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads,
    use_finite_differences=False,
    lennard_jones_potential="6_12")
  assert approx_equal(residual, -60.38271638, eps=0.00001)
  simple_bonds = hbond.get_simple_bonds(build_proxies.proxies)
  assert (simple_bonds.size() == 1)
  assert (list(simple_bonds[0]) == [0,1])

def exercise_implicit () :
  from mmtbx.geometry_restraints import hbond
  import cctbx.geometry_restraints
  from scitbx.array_family import flex
  pdb_in = simple_pdb()
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  sites_cart = hierarchy.atoms().extract_xyz()
  cache = hierarchy.atom_selection_cache()
  xrs = pdb_in.xray_structure_simple(unit_cube_pseudo_crystal=True)
  params = hbond.master_phil.extract()
  params.implicit.theta_low = 110.0
  params.implicit.theta_high = 160.0
  build_proxies = hbond.find_implicit_hydrogen_bonds(
    pdb_hierarchy=hierarchy,
    xray_structure=xrs,
    params=params,
    log=None)
  assert (len(build_proxies.proxies) == 1)
  grads_an = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual_an = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=grads_an)
  #assert approx_equal(residual_fd, -0.11357803, eps=0.000001)
  assert approx_equal(residual_an, -54.517427, eps=0.00001)
  (21.36455693329907, 13.332353726790171, 17.10498505822002)
  assert approx_equal(grads_an[2][0], 21.36455693, eps=0.00001)
  assert approx_equal(grads_an[2][2], 17.104985058, eps=0.00001)
  assert approx_equal(grads_an[3][1], -29.0568888, eps=0.00001)
  # only three sites, modified directly
  build_proxies = hbond.build_implicit_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[0,1,2],
    distance_ideal=2.9,
    distance_cut=3.5,
    theta_low=110.0,
    theta_high=150.0,
    weight=1.0)
  sites_cart = flex.vec3_double([(0.0,2.9,0.0), (0.0,0.0,0.0),
    (0.5774,-1.0,0.0)])
  sites_cart[0] = (0.0, 3.0, 0.0)
  g1 = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  r1 = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=g1)
  assert approx_equal(r1, -58.5286434789, eps=0.00001)
  sites_cart[0] = (0.0, 2.8, 0.0)
  g2 = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  r2 = hbond.target_and_gradients(
    proxies=build_proxies.proxies,
    sites_cart=sites_cart,
    gradient_array=g2)
  assert approx_equal(r2, -58.273999, eps=0.000001)
  assert approx_equal(g2[0][1], -21.2470231559, eps=0.00001)
  assert approx_equal(g2[0][0], 0.00310495, eps=0.00001)
  assert approx_equal(g2[2][0], 0.00652011, eps=0.00001)
  simple_bonds = hbond.get_simple_bonds(build_proxies.proxies)
  assert (simple_bonds.size() == 1)
  assert (list(simple_bonds[0]) == [0,1])
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is not None) :
    from iotbx import file_reader
    pdb_in = file_reader.any_file(pdb_file, force_type="pdb").file_object
    hierarchy = pdb_in.construct_hierarchy()
    xrs = pdb_in.xray_structure_simple()
    build_proxies = hbond.find_implicit_hydrogen_bonds(
      pdb_hierarchy=hierarchy,
      xray_structure=xrs,
      params=params,
      log=None)
    assert (len(build_proxies.proxies) == 232)
  # test proxy_select
  build_proxies = hbond.build_implicit_hbond_proxies()
  build_proxies.add_proxy(
    i_seqs=[0,1,2],
    distance_ideal=2.9,
    distance_cut=3.5,
    theta_low=110.0,
    theta_high=150.0,
    weight=1.0)
  build_proxies.add_proxy(
    i_seqs=[7,8,9],
    distance_ideal=2.9,
    distance_cut=3.5,
    theta_low=110.0,
    theta_high=150.0,
    weight=1.0)
  build_proxies.add_proxy(
    i_seqs=[14,15,16],
    distance_ideal=2.9,
    distance_cut=3.5,
    theta_low=110.0,
    theta_high=150.0,
    weight=1.0)
  iselection = flex.size_t([0,1,2,3,4,5,12,13,14,15,16,17,18])
  proxies = build_proxies.proxies
  selected = proxies.proxy_select(n_seq=20, iselection=iselection)
  assert (selected.size() == 2)
  assert (selected[1].i_seqs == (8,9,10))

def exercise_switching_function () :
  from mmtbx.geometry_restraints import hbond
  from cctbx.array_family import flex
  from cctbx.geometry import distance
  import boost.python
  ext = boost.python.import_ext("mmtbx_hbond_restraints_ext")
  sites = [[0.,0.,0.], [3.49, 0.,0.]]
  x = []
  y0 = []
  y1 = []
  y2 = []
  y3 = []
  for i in range(1000) :
    sites[1][0] += 0.0001
    d = distance(sites)
    sw = ext.switch_fn(d.distance_model, 3.5, 3.55)
    dsw_dr = ext.d_switch_d_distance(d.distance_model, 3.5, 3.55)
    d2sw_dr2 = ext.d2_switch_d_distance2(d.distance_model, 3.5, 3.55)
    x.append(d.distance_model)
    y0.append(sw)
    y1.append(dsw_dr)
    y2.append(d2sw_dr2)
    sites[1][0] -= 0.0000001
    d = distance(sites)
    dsw_dr_1 = ext.d_switch_d_distance(d.distance_model, 3.5, 3.55)
    sites[1][0] += 0.0000002
    d = distance(sites)
    dsw_dr_2 = ext.d_switch_d_distance(d.distance_model, 3.5, 3.55)
    d2sw_dr2_fd = (dsw_dr_2 - dsw_dr_1) / 0.0000002
    sites[1][0] -= 0.0000001
    y3.append(d2sw_dr2_fd)
    # avoid corners!
    if ((sites[1][0] < 3.5) or ((sites[1][0] > 3.500001) and
        (sites[1][0] < 3.55)) or (sites[1][0] > 3.550001)) :
      assert approx_equal(d2sw_dr2_fd, d2sw_dr2)
  #from matplotlib import pyplot as plt
  #fig = plt.figure()
  #p = fig.add_subplot(111)
  #p.plot(x, y0, color='r')
  #p.plot(x, y1, color='g')
  #p.plot(x, y2, color='b')
  #p.plot(x, y3, color='k')
  #plt.show()

def plot_potentials () :
  from mmtbx.geometry_restraints import hbond
  from scitbx.array_family import flex
  lj_proxies = hbond.build_lennard_jones_proxies()
  lj_proxies.add_proxy(
    i_seqs=[0,1],
    distance_ideal=2.9,
    distance_cut=3.5)
  distances = []
  for i in range(200) :
    distances.append(2.5 + (i * 0.005))
  y_vals = []
  y2_vals = []
  y3_vals = []
  for dist in distances :
    sites_cart = flex.vec3_double([(0.0,dist,0.0), (0.0,0.0,0.0),
      (0.5774,-1.0,0.0)])
    g = flex.vec3_double(3, (0.0,0.0,0.0))
    r = hbond.target_and_gradients(
      proxies=lj_proxies.proxies,
      sites_cart=sites_cart,
      gradient_array=g)
    y_vals.append(r)
    y2_vals.append(g[0][1])
    g2 = flex.vec3_double(3, (0.0,0.0,0.0))
    r = hbond.target_and_gradients(
      proxies=lj_proxies.proxies,
      sites_cart=sites_cart,
      gradient_array=g2,
      use_finite_differences=True)
    y3_vals.append(g2[0][1])
  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(distances, y_vals, '-')
  ax.plot(distances, y2_vals, '-')
  ax.plot(distances, y3_vals, '-')
  ax.set_title("Lennard-Jones H-bond potential")
  plt.show()

if (__name__ == "__main__") :
  if ("--plot" in sys.argv) :
    plot_potentials()
  else :
    exercise_simple()
    exercise_lennard_jones()
    exercise_implicit()
    exercise_switching_function()
    print "OK"
