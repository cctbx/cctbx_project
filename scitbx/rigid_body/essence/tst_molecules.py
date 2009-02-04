from scitbx.rigid_body.essence import featherstone
from scitbx.rigid_body.essence import joint_lib
from scitbx.rigid_body.essence import utils
from scitbx.graph import tardy_tree
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
from libtbx.test_utils import approx_equal
import random
import sys

def random_wells(sites):
  result = []
  for site in sites:
    result.append(site+matrix.col.random(n=3, a=-1, b=1))
  return result

class potential_object(object):

  def __init__(O,
        sites,
        wells,
        restraint_edges,
        restraint_edge_weight=1/0.1**2,
        epsilon=1.e-100):
    O.wells = wells
    O.restraints = []
    for edge in restraint_edges:
      s = [sites[i] for i in edge]
      O.restraints.append((edge, abs(s[0]-s[1]), restraint_edge_weight))
    O.epsilon = epsilon

  def e_pot(O, sites_moved):
    result = 0
    for s, w in zip(sites_moved, O.wells):
      result += (s - w).dot()
    for edge,d_ideal,w in O.restraints:
      s = [sites_moved[i] for i in edge]
      d_model = abs(s[0]-s[1])
      if (d_model < O.epsilon): continue
      delta = d_ideal - d_model
      result += w * delta**2
    return result

  def d_e_pot_d_sites(O, sites_moved):
    result = []
    for s, w in zip(sites_moved, O.wells):
      result.append(2 * (s - w))
    for edge,d_ideal,w in O.restraints:
      s = [sites_moved[i] for i in edge]
      d_model = abs(s[0]-s[1])
      if (d_model < O.epsilon): continue
      delta = d_ideal - d_model
      g0 = -w * 2 * delta / d_model * (s[0] - s[1])
      result[edge[0]] += g0
      result[edge[1]] -= g0
    return result

class simulation(object):

  def __init__(O,
        labels, sites, bonds, cluster_manager, potential_obj, bodies):
    O.labels = labels
    O.sites = sites
    O.bonds = bonds
    O.cluster_manager = cluster_manager
    O.potential_obj = potential_obj
    O.bodies = bodies
    O.energies_and_accelerations_update()

  def AJA_update(O):
    O.AJA = []
    for B in O.bodies:
      AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
      if (B.parent != -1):
        AJA = O.AJA[B.parent] * AJA
      O.AJA.append(AJA)

  def JAr_update(O):
    O.JAr = []
    for B in O.bodies:
      JAr = B.J.Tps.r * B.A.T0b.r
      if (B.parent != -1):
        JAr *= O.AJA[B.parent].r.transpose()
      O.JAr.append(JAr)

  def sites_moved_update(O):
    O.sites_moved = [None] * len(O.sites)
    n_done = 0
    for iB,B in enumerate(O.bodies):
      AJA = O.AJA[iB]
      for i_seq in O.cluster_manager.clusters[iB]:
        assert O.sites_moved[i_seq] is None
        O.sites_moved[i_seq] = AJA * O.sites[i_seq]
        n_done += 1
    assert n_done == len(O.sites)

  def f_ext_bf_update(O, d_e_pot_d_sites):
    O.f_ext_bf = []
    for iB,B in enumerate(O.bodies):
      f = matrix.col((0,0,0))
      nc = matrix.col((0,0,0))
      for i_seq in O.cluster_manager.clusters[iB]:
        s = O.sites[i_seq]
        force_bf = -(O.JAr[iB] * d_e_pot_d_sites[i_seq])
        f += force_bf
        nc += (B.A.T0b * s).cross(force_bf)
      O.f_ext_bf.append(matrix.col((nc, f)).resolve_partitions())

  def energies_and_accelerations_update(O):
    model = featherstone.system_model(bodies=O.bodies)
    O.e_kin = model.e_kin()
    O.e_pot_and_f_ext_update()
    O.qdd = model.FDab(tau=None, f_ext=O.f_ext_bf)

  def e_pot_and_f_ext_update(O):
    O.AJA_update()
    O.JAr_update()
    O.sites_moved_update()
    O.e_pot = O.potential_obj.e_pot(sites_moved=O.sites_moved)
    O.f_ext_bf_update(
      d_e_pot_d_sites=O.potential_obj.d_e_pot_d_sites(
        sites_moved=O.sites_moved))
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def d_pot_d_q(O):
    return featherstone.system_model(bodies=O.bodies).d_pot_d_q(
      f_ext=O.f_ext_bf)

  def d_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for B in O.bodies:
      gs = []
      J_orig = B.J
      q_orig = list(J_orig.get_q())
      for iq in xrange(J_orig.q_size):
        fs = []
        for signed_eps in [eps, -eps]:
          q_eps = list(q_orig)
          q_eps[iq] += signed_eps
          B.J = J_orig.new_q(q=q_eps)
          O.e_pot_and_f_ext_update()
          fs.append(O.e_pot)
        gs.append((fs[0]-fs[1])/(2*eps))
      B.J = J_orig
      result.append(matrix.col(gs))
    O.energies_and_accelerations_update()
    return result

  def check_d_pot_d_q(O, verbose=0):
    qdd_orig = O.qdd
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    if (verbose):
      for a,f in zip(ana, fin):
        print "fin:", f.elems
        print "ana:", a.elems
      print
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd, qdd_orig)

  def minimization(O, max_iterations=None, callback_after_step=None):
    refinery(
      sim=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)

class refinery(object):

  def __init__(O, sim, max_iterations=None, callback_after_step=None):
    O.sim = sim
    O.callback_after_step = callback_after_step
    O.x = flex.double()
    for B in sim.bodies:
      O.x.extend(flex.double(B.J.get_q()))
    import scitbx.lbfgs
    scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=scitbx.lbfgs.termination_parameters(
       max_iterations=max_iterations))
    O.sim.energies_and_accelerations_update()

  def unpack_x(O):
    x = O.x
    i = 0
    for B in O.sim.bodies:
      n = B.J.q_size
      B.J = B.J.new_q(q=x[i:i+n])
      i += n
    assert i == x.size()
    O.sim.e_pot_and_f_ext_update()

  def compute_functional_and_gradients(O):
    O.unpack_x()
    f = O.sim.e_pot
    g = flex.double()
    for d in O.sim.d_pot_d_q():
      g.extend(flex.double(d))
    return f, g

class six_dof_body(object):

  def __init__(O, sites):
    O.A = joint_lib.six_dof_alignment(
      center_of_mass=utils.center_of_mass_from_sites(sites=sites))
    O.I = utils.spatial_inertia_from_sites(sites=sites, alignment_T=O.A.T0b)
    #
    qE = matrix.col((1,0,0,0))
    qr = matrix.col((0,0,0))
    O.J = joint_lib.six_dof(qE=qE, qr=qr)
    O.qd = O.J.qd_zero

class revolute_body(object):

  def __init__(O, sites, pivot, normal):
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = utils.spatial_inertia_from_sites(sites=sites, alignment_T=O.A.T0b)
    #
    O.J = joint_lib.revolute(qE=matrix.col([0]))
    O.qd = O.J.qd_zero

def pdb_extract(pdb):
  labels, sites = [], []
  for line in pdb.splitlines():
    labels.append(line[22:26].strip()+"."+line[12:16].strip())
    sites.append(matrix.col([float(line[30+i*8:38+i*8]) for i in [0,1,2]]))
  return labels, sites

def construct_bodies(sites, cluster_manager):
  result = []
  cm = cluster_manager
  for ic,cluster in enumerate(cm.clusters):
    he = cm.hinge_edges[ic]
    if (he[0] == -1):
      body = six_dof_body(
        sites=[matrix.col(sites[i]) for i in cluster])
      body.parent = -1
    else:
      normal_sites = [matrix.col(sites[i]) for i in he]
      body = revolute_body(
        sites=[matrix.col(sites[i]) for i in cluster],
        pivot=normal_sites[1],
        normal=(normal_sites[1]-normal_sites[0]).normalize())
      body.parent = cm.cluster_indices[he[0]]
    body.i_seqs = cluster
    result.append(body)
  return result

def construct_simulation(labels, sites, bonds, size_max=8):
  tt = tardy_tree.construct(
    n_vertices=len(sites), edge_list=bonds, size_max=size_max)
  cm = tt.cluster_manager
  cm.merge_clusters_with_multiple_connections(edge_sets=tt.edge_sets)
  cm.construct_spanning_trees(edge_sets=tt.edge_sets)
  cm.find_loop_edge_bendings(edge_sets=tt.edge_sets)
  return simulation(
    labels=labels,
    sites=sites,
    bonds=bonds,
    cluster_manager=cm,
    potential_obj=potential_object(
      sites=sites,
      wells=random_wells(sites),
      restraint_edges=cm.loop_edges+cm.loop_edge_bendings),
    bodies=construct_bodies(sites=sites, cluster_manager=cm))

def construct_sim(pdb, bonds):
  labels, sites = pdb_extract(pdb=pdb)
  return construct_simulation(labels=labels, sites=sites, bonds=bonds)

def simulation_gly_no_h():
  pdb = """\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
"""
  bonds = [(0,1),(1,2),(2,3)]
  return construct_sim(pdb=pdb, bonds=bonds)

def simulation_gly_with_nh():
  pdb = """\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
ATOM      4  H   GLY A   1      11.792  12.691  15.311  0.00  0.00           H
"""
  bonds = [(0,1),(0,4),(1,2),(2,3)]
  return construct_sim(pdb=pdb, bonds=bonds)

def simulation_ala_no_h():
  pdb = """\
ATOM      0  N   ALA A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  ALA A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   ALA A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  CB  ALA A   1      10.908  13.950  17.351  0.00  0.00           C
ATOM      4  O   ALA A   1       9.916  16.090  14.936  0.00  0.00           O
"""
  bonds = [(0,1),(1,2),(1,3),(2,4)]
  return construct_sim(pdb=pdb, bonds=bonds)

def simulation_ala_with_h():
  pdb = """\
ATOM      0  N   ALA A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  ALA A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   ALA A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  HA  ALA A   1       9.428  13.887  15.936  0.00  0.00           H
ATOM      4  O   ALA A   1       9.916  16.090  14.936  0.00  0.00           O
ATOM      5  H   ALA A   1      11.792  12.691  15.311  0.00  0.00           H
ATOM      6  CB  ALA A   1      10.908  13.950  17.351  0.00  0.00           C
ATOM      7  HB1 ALA A   1      10.627  13.138  17.778  0.00  0.00           H
ATOM      8  HB2 ALA A   1      10.540  14.707  17.813  0.00  0.00           H
ATOM      9  HB3 ALA A   1      11.867  14.004  17.346  0.00  0.00           H
"""
  bonds = [(0,1),(0,5),(1,2),(1,3),(1,6),(2,4),(6,7),(6,8),(6,9)]
  return construct_sim(pdb=pdb, bonds=bonds)

def simulation_tyr_with_h():
  pdb = """\
ATOM      0  CG  TYR A   1      11.007   9.417   9.446  1.00  0.79           C
ATOM      1  CD1 TYR A   1       9.923  10.155   8.940  1.00  1.42           C
ATOM      2  CD2 TYR A   1      10.765   8.288  10.238  1.00  1.41           C
ATOM      3  CE1 TYR A   1       8.612   9.760   9.229  1.00  1.61           C
ATOM      4  CE2 TYR A   1       9.453   7.895  10.525  1.00  1.42           C
ATOM      5  CZ  TYR A   1       8.377   8.631  10.021  1.00  1.11           C
ATOM      6  HD1 TYR A   1      10.092  11.024   8.328  1.00  2.14           H
ATOM      7  HD2 TYR A   1      11.596   7.718  10.630  1.00  2.21           H
ATOM      8  HE1 TYR A   1       7.780  10.329   8.841  1.00  2.44           H
ATOM      9  HE2 TYR A   1       9.270   7.023  11.135  1.00  2.13           H
ATOM     10  OH  TYR A   1       7.083   8.244  10.304  1.00  1.32           O
ATOM     11  HH  TYR A   1       6.494   8.723   9.717  1.00  2.00           H
ATOM     12  CB  TYR A   1      12.440   9.818   9.148  1.00  0.74           C
ATOM     13  HB2 TYR A   1      12.827   9.193   8.358  1.00  0.78           H
ATOM     14  HB3 TYR A   1      13.036   9.677  10.037  1.00  0.78           H
ATOM     15  N   TYR A   1      11.593  12.101   9.550  1.00  0.82           N
ATOM     16  CA  TYR A   1      12.527  11.286   8.721  1.00  0.75           C
ATOM     17  C   TYR A   1      12.160  11.413   7.239  1.00  0.76           C
ATOM     18  HA  TYR A   1      13.536  11.638   8.870  1.00  0.85           H
ATOM     19  O   TYR A   1      12.298  12.462   6.643  1.00  0.83           O
ATOM     20  H   TYR A   1      10.948  12.701   9.122  1.00  0.88           H
"""
  bonds=[
    (0, 1), (0, 2), (0, 12), (1, 3), (1, 6), (2, 4), (2, 7), (3, 5),
    (3, 8), (4, 5), (4, 9), (5, 10), (10, 11), (12, 13), (12, 14),
    (12, 16), (15, 16), (15, 20), (16, 17), (16, 18), (17, 19)]
  return construct_sim(pdb=pdb, bonds=bonds)

def simulation_van_fragment():
  # PDB code 1qd8
  pdb = """\
HETATM   47  C44 VAN     1       0.718   5.411   2.269  1.00  3.66           C
HETATM   48  C47 VAN     1       0.913   3.899   2.010  1.00  3.90           C
HETATM   49  C48 VAN     1       1.937   3.222   2.700  1.00  3.87           C
HETATM   50  C50 VAN     1       2.013   1.800   2.726  1.00  4.31           C
HETATM   51  C52 VAN     1       1.044   1.071   2.004  1.00  4.71           C
HETATM   52  O53 VAN     1       1.077  -0.297   2.083  1.00  6.04           O
HETATM   53  C51 VAN     1       0.066   1.767   1.246  1.00  5.04           C
HETATM   54  C49 VAN     1       0.006   3.153   1.237  1.00  4.22           C
HETATM   55  C45 VAN     1      -0.704   5.507   2.917  1.00  3.80           C
HETATM   56  O46 VAN     1      -1.666   5.899   2.216  1.00  4.55           O
HETATM   57  N54 VAN     1      -0.869   5.098   4.194  1.00  4.02           N
HETATM   58  C55 VAN     1       0.141   4.679   5.156  1.00  4.23           C
HETATM   69  C56 VAN     1       0.045   3.151   5.360  1.00  4.51           C
HETATM   70  O57 VAN     1      -1.022   2.542   5.259  1.00  5.51           O
HETATM   71  N68 VAN     1       1.226   2.525   5.641  1.00  4.67           N
HETATM   72  C69 VAN     1       1.340   1.069   5.484  1.00  5.14           C
HETATM   73  C72 VAN     1       2.754   0.737   4.924  1.00  4.58           C
HETATM   74  C73 VAN     1       3.049   1.131   3.577  1.00  4.34           C
HETATM   75  C75 VAN     1       4.354   0.899   3.099  1.00  4.29           C
HETATM   76  O79 VAN     1       4.736   1.286   1.834  1.00  4.73           O
HETATM   77  C76 VAN     1       5.342   0.269   3.883  1.00  4.57           C
HETATM   78  C77 VAN     1       5.023  -0.107   5.195  1.00  4.56           C
HETATM   79  O78 VAN     1       5.912  -0.707   6.053  1.00  5.51           O
HETATM   80  C74 VAN     1       3.724   0.123   5.713  1.00  4.95           C
HETATM   81  C70 VAN     1       1.069   0.287   6.838  1.00  5.82           C
HETATM   82  O71 VAN     1       1.149   0.912   7.924  1.00  7.21           O
HETATM   83  O80 VAN     1       0.816  -0.957   6.693  1.00  8.10           O
"""
  bonds=[
    (0, 1), (0, 8), (1, 2), (1, 7), (2, 3), (3, 4), (3, 17), (4, 5), (4, 6),
    (6, 7), (8, 9), (8, 10), (10, 11), (11, 12), (12, 13), (12, 14), (14, 15),
    (15, 16), (15, 24), (16, 17), (16, 23), (17, 18), (18, 19), (18, 20),
    (20, 21), (21, 22), (21, 23), (24, 25), (24, 26)]
  return construct_sim(pdb=pdb, bonds=bonds)

simulation_factories = [
  simulation_gly_no_h,
  simulation_gly_with_nh,
  simulation_ala_no_h,
  simulation_ala_with_h,
  simulation_tyr_with_h,
  simulation_van_fragment]

def exercise_sim(out, n_dynamics_steps, delta_t, sim):
  sim.check_d_pot_d_q()
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  sim.check_d_pot_d_q()
  print >> out, "energy samples:", e_tots.size()
  print >> out, "e_pot min, max:", min(e_pots), max(e_pots)
  print >> out, "e_kin min, max:", min(e_kins), max(e_kins)
  print >> out, "e_tot min, max:", min(e_tots), max(e_tots)
  print >> out, "start e_tot:", e_tots[0]
  print >> out, "final e_tot:", e_tots[-1]
  ave = flex.sum(e_tots) / e_tots.size()
  range = flex.max(e_tots) - flex.min(e_tots)
  if (ave == 0): relative_range = 0
  else:          relative_range = range / ave
  print >> out, "ave:", ave
  print >> out, "range:", range
  print >> out, "relative range:", relative_range
  print >> out
  out.flush()
  return relative_range

def exercise_dynamics_quick(out, sim, n_dynamics_steps, delta_t=0.001):
  relative_range = exercise_sim(
    out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
  if (out is not sys.stdout):
    assert relative_range < 1.e-4
  print >> out

def exercise_minimization_quick(out, sim, max_iterations=3):
  print >> out, "Minimization:"
  print >> out, "  start e_pot:", sim.e_pot
  e_pot_start = sim.e_pot
  sim.minimization(max_iterations=max_iterations)
  print >> out, "  final e_pot:", sim.e_pot
  e_pot_final = sim.e_pot
  if (out is not sys.stdout):
    assert e_pot_final < e_pot_start * 0.98
  print >> out

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 30
    out = null_out()
  else:
    n_dynamics_steps = max(1, int(args[0]))
    out = sys.stdout
  show_times_at_exit()
  random.seed(0)
  for sim_factory in simulation_factories:
    sim = sim_factory()
    exercise_dynamics_quick(
      out=out, sim=sim, n_dynamics_steps=n_dynamics_steps)
    exercise_minimization_quick(out=out, sim=sim)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
