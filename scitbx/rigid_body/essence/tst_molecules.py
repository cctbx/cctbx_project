from scitbx.rigid_body.essence import featherstone
from scitbx.rigid_body.essence import joint_lib
from scitbx.rigid_body.essence import utils
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
from libtbx.test_utils import approx_equal
import sys

def potential_energy(sites, wells, A, J, AJA_tree=None):
  result = 0
  AJA = A.Tb0 * J.Tsp * A.T0b
  if (AJA_tree is not None): AJA = AJA_tree * AJA
  for s, w in zip(sites, wells):
    result += (AJA * s - w).dot()
  return result

def potential_f_ext_bf(sites, wells, A, J, AJA_tree=None):
  AJA = A.Tb0 * J.Tsp * A.T0b
  if (AJA_tree is not None): AJA = AJA_tree * AJA
  f_cart_ff = [-2 * (AJA * s - w) for s, w in zip(sites, wells)]
  JAr = J.Tps.r * A.T0b.r
  if (AJA_tree is not None): JAr *= AJA_tree.r.transpose()
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force_ff in zip(sites, f_cart_ff):
    force_bf = JAr * force_ff
    f += force_bf
    nc += (A.T0b * s).cross(force_bf)
  return matrix.col((nc, f)).resolve_partitions()

class simulation(object):

  def __init__(O, bodies):
    O.bodies = bodies
    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):
    model = featherstone.system_model(bodies=O.bodies)
    O.e_kin = model.e_kin()
    O.e_pot_and_f_ext_update()
    O.qdd = model.FDab(tau=None, f_ext=O.f_ext_bf)

  def e_pot_and_f_ext_update(O):
    O.AJA_accu = []
    O.e_pot = 0
    O.f_ext_bf = []
    for B in O.bodies:
      AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
      if (B.parent == -1):
        AJA_tree = None
      else:
        AJA_tree = O.AJA_accu[B.parent]
        AJA = AJA_tree * AJA
      O.AJA_accu.append(AJA)
      e_pot_bf = potential_energy(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      f_ext_using_bf = potential_f_ext_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      O.f_ext_bf.append(f_ext_using_bf)
      O.e_pot += e_pot_bf
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

def create_wells(sites, mersenne_twister, r=None):
  "overall random rotation and translation + noise"
  if (r is None):
    r = matrix.sqr(mersenne_twister.random_double_r3_rotation_matrix())
  t = matrix.col(mersenne_twister.random_double(size=3)-0.5)
  wells = []
  for site in sites:
    t_noise = t + matrix.col(mersenne_twister.random_double(size=3)-0.5)*0.2
    wells.append(r * site + t_noise)
  return wells

def shift_gently(sites, mersenne_twister, angle=5):
  axis = mersenne_twister.random_double_point_on_sphere()
  r = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    axis=axis, angle=angle, deg=True))
  return create_wells(sites=sites, mersenne_twister=mersenne_twister, r=r)

class six_dof_body(object):

  def __init__(O, labels, sites, bonds, mersenne_twister):
    O.labels = labels
    O.sites = sites
    O.bonds = bonds
    O.A = joint_lib.six_dof_alignment(
      center_of_mass=utils.center_of_mass_from_sites(sites=sites))
    O.I = utils.spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = shift_gently(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col((1,0,0,0))
    qr = matrix.col((0,0,0))
    O.J = joint_lib.six_dof(qE=qE, qr=qr)
    O.qd = O.J.qd_zero

class revolute_body(object):

  def __init__(O, labels, sites, bonds, pivot, normal, mersenne_twister):
    O.labels = labels
    O.sites = sites
    O.bonds = bonds
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = utils.spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = shift_gently(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    O.J = joint_lib.revolute(qE=matrix.col([0]))
    O.qd = O.J.qd_zero

def simulation_zigzag(n_bodies=5):
  mersenne_twister = flex.mersenne_twister(seed=0)
  body = six_dof_body(
    labels=["00", "01", "02"],
    sites=matrix.col_list([
      (0.3,-0.5,0),
      (0.4,0.5,0),
      (0,0,0)]),
    bonds=[(0,2),(1,2)],
    mersenne_twister=mersenne_twister)
  body.parent = -1
  bodies = [body]
  vu = matrix.col((0,1,0)).rotate(axis=matrix.col((1,0,0)), angle=75, deg=True)
  vr = matrix.col((0,1,0))
  v = vu
  pivot = matrix.col((0,0,0))
  for ib in xrange(1,n_bodies):
    body = revolute_body(
      labels=[str(ib)],
      sites=[pivot + v*0.5],
      bonds=[(-1,0)],
      pivot=pivot,
      normal=matrix.col((1,0,0)),
      mersenne_twister=mersenne_twister)
    body.parent = ib-1
    bodies.append(body)
    pivot += v
    if (v is vu): v = vr
    else:         v = vu
  return simulation(bodies=bodies)

def pdb_extract(pdb):
  labels, sites = [], []
  for line in pdb.splitlines():
    labels.append(line[22:26].strip()+"."+line[12:16].strip())
    sites.append(matrix.col([float(line[30+i*8:38+i*8]) for i in [0,1,2]]))
  return labels, sites

def simulation_gly_no_h():
  pdb = """\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
"""
  labels, sites = pdb_extract(pdb=pdb)
  mersenne_twister = flex.mersenne_twister(seed=0)
  body0 = six_dof_body(
    labels=labels[:3],
    sites=sites[:3],
    bonds=[(0,1),(1,2)],
    mersenne_twister=mersenne_twister)
  body0.parent = -1
  body1 = revolute_body(
    labels=labels[3:],
    sites=sites[3:],
    bonds=[(-1,0)],
    pivot=sites[2],
    normal=(sites[2]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  return simulation(bodies=[body0, body1])

def simulation_gly_with_nh():
  pdb = """\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
ATOM      4  H   GLY A   1      11.792  12.691  15.311  0.00  0.00           H
"""
  labels, sites = pdb_extract(pdb=pdb)
  mersenne_twister = flex.mersenne_twister(seed=0)
  body0 = six_dof_body(
    labels=labels[:3],
    sites=sites[:3],
    bonds=[(0,1),(1,2)],
    mersenne_twister=mersenne_twister)
  body0.parent = -1
  body1 = revolute_body(
    labels=labels[3:4],
    sites=sites[3:4],
    bonds=[(-1,0)],
    pivot=sites[2],
    normal=(sites[2]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  body2 = revolute_body(
    labels=labels[4:],
    sites=sites[4:],
    bonds=[(-3,0)],
    pivot=sites[0],
    normal=(sites[0]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body2.parent = 0
  return simulation(bodies=[body0, body1, body2])

def simulation_ala_no_h():
  pdb = """\
ATOM      0  N   ALA A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  ALA A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   ALA A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  CB  ALA A   1      10.908  13.950  17.351  0.00  0.00           C
ATOM      4  O   ALA A   1       9.916  16.090  14.936  0.00  0.00           O
"""
  labels, sites = pdb_extract(pdb=pdb)
  mersenne_twister = flex.mersenne_twister(seed=0)
  body0 = six_dof_body(
    labels=labels[:4],
    sites=sites[:4],
    bonds=[(0,1),(1,2),(1,3)],
    mersenne_twister=mersenne_twister)
  body0.parent = -1
  body1 = revolute_body(
    labels=labels[4:],
    sites=sites[4:],
    bonds=[(-2,0)],
    pivot=sites[2],
    normal=(sites[2]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  return simulation(bodies=[body0, body1])

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
  labels, sites = pdb_extract(pdb=pdb)
  mersenne_twister = flex.mersenne_twister(seed=0)
  body0 = six_dof_body(
    labels=labels[:4],
    sites=sites[:4],
    bonds=[(0,1),(1,2),(1,3)],
    mersenne_twister=mersenne_twister)
  body0.parent = -1
  body1 = revolute_body(
    labels=labels[4:5],
    sites=sites[4:5],
    bonds=[(-2,0)],
    pivot=sites[2],
    normal=(sites[2]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  body2 = revolute_body(
    labels=labels[5:6],
    sites=sites[5:6],
    bonds=[(-4,0)],
    pivot=sites[0],
    normal=(sites[0]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body2.parent = 0
  body3 = revolute_body(
    labels=labels[6:],
    sites=sites[6:],
    bonds=[(-3,0),(0,1),(0,2),(0,3)],
    pivot=sites[6],
    normal=(sites[6]-sites[1]).normalize(),
    mersenne_twister=mersenne_twister)
  body3.parent = 0
  return simulation(bodies=[body0, body1, body2, body3])

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
  labels, sites = pdb_extract(pdb=pdb)
  mersenne_twister = flex.mersenne_twister(seed=0)
  body0 = six_dof_body(
    labels=labels[:11],
    sites=sites[:11],
    bonds=[(0,1),(0,2),(1,3),(2,4),(3,5),(4,5),(5,10),(1,6),(2,7),(3,8),(4,9)],
    mersenne_twister=mersenne_twister)
  body0.parent = -1
  body1 = revolute_body(
    labels=labels[11:12],
    sites=sites[11:12],
    bonds=[(-1,0)],
    pivot=sites[10],
    normal=(sites[10]-sites[5]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  body2 = revolute_body(
    labels=labels[12:15],
    sites=sites[12:15],
    bonds=[(-11,0),(0,1),(0,2)],
    pivot=sites[12],
    normal=(sites[12]-sites[0]).normalize(),
    mersenne_twister=mersenne_twister)
  body2.parent = 0
  body3 = revolute_body(
    labels=labels[15:19],
    sites=sites[15:19],
    bonds=[(-3,1),(0,1),(1,2),(1,3)],
    pivot=sites[16],
    normal=(sites[16]-sites[12]).normalize(),
    mersenne_twister=mersenne_twister)
  body3.parent = 2
  body4 = revolute_body(
    labels=labels[19:20],
    sites=sites[19:20],
    bonds=[(-2,0)],
    pivot=sites[17],
    normal=(sites[17]-sites[16]).normalize(),
    mersenne_twister=mersenne_twister)
  body4.parent = 3
  body5 = revolute_body(
    labels=labels[20:21],
    sites=sites[20:21],
    bonds=[(-4,0)],
    pivot=sites[15],
    normal=(sites[15]-sites[16]).normalize(),
    mersenne_twister=mersenne_twister)
  body5.parent = 3
  return simulation(bodies=[body0, body1, body2, body3, body4, body5])

simulation_factories = [
  simulation_zigzag,
  simulation_gly_no_h,
  simulation_gly_with_nh,
  simulation_ala_no_h,
  simulation_ala_with_h,
  simulation_tyr_with_h]

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
    assert e_pot_final < e_pot_start * 0.9
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
  for sim_factory in simulation_factories:
    sim = sim_factory()
    exercise_dynamics_quick(
      out=out, sim=sim, n_dynamics_steps=n_dynamics_steps)
    exercise_minimization_quick(out=out, sim=sim)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
