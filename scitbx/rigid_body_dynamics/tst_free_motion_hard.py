from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.tst_free_motion import \
  featherstone_system_model
from scitbx.rigid_body_dynamics.free_motion_reference_impl import \
  body_inertia, \
  create_triangle_with_center_of_mass_at_origin
from scitbx.array_family import flex
from scitbx import matrix
import sys

def create_wells(sites, mersenne_twister):
  "overall random rotation and translation + noise"
  r = matrix.sqr(mersenne_twister.random_double_r3_rotation_matrix())
  t = matrix.col(mersenne_twister.random_double(size=3)-0.5)
  wells = []
  for site in sites:
    t_noise = t + matrix.col(mersenne_twister.random_double(size=3)-0.5)*0.2
    wells.append(r * site + t_noise)
  return wells

def potential_energy(sites, wells, T_inv):
  result = 0
  for xyz_s, xyz_w in zip(sites, wells):
    result += (xyz_s - T_inv * xyz_w).dot()
  return result

def potential_f_ext(sites, wells, T_inv):
  f_cart = [2 * (xyz_s - T_inv * xyz_w) for xyz_s, xyz_w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for xyz,force in zip(sites, f_cart):
    f += force
    nc += xyz.cross(force)
  return matrix.col((nc, f)).resolve_partitions()

class kinetic_energy(object):

  def __init__(O, m, I, v_spatial):
    "Shabana (2005) p. 148 eqs. 3.125, 3.126"
    w, v = matrix.col_list([v_spatial.elems[:3], v_spatial.elems[3:]])
    O.ang = 0.5 * w.dot(I * w)
    O.lin = 0.5 * m * v.dot()
    O.tot = O.ang + O.lin

class simulation(object):

  def __init__(O, mersenne_twister):
    O.sites = create_triangle_with_center_of_mass_at_origin()
    O.m = len(O.sites) # unit masses
    O.I = body_inertia(sites_cart=O.sites)
    #
    O.wells = create_wells(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    O.J = joint_lib.six_dof_joint_euler_params(qE=qE, qr=qr)
    O.v_spatial = matrix.col(mersenne_twister.random_double(size=6)*2-1)
    #
    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):
    O.e_kin = kinetic_energy(m=O.m, I=O.I, v_spatial=O.v_spatial)
    O.e_pot = potential_energy(sites=O.sites, wells=O.wells, T_inv=O.J.T_inv)
    O.f_ext = potential_f_ext(sites=O.sites, wells=O.wells, T_inv=O.J.T_inv)
    O.e_tot = O.e_kin.tot + O.e_pot
    #
    model = featherstone_system_model(m=O.m, I=O.I, J=O.J)
    q = [None] # already stored in joint as qE and qr
    qd = [O.v_spatial]
    tau = [matrix.col((0,0,0,0,0,0))]
    grav_accn = [0,0,0]
    qdd = featherstone.FDab(model, q, qd, tau, [O.f_ext], grav_accn)
    O.a_spatial = qdd[0]

  def dynamics_step(O, delta_t):
    O.v_spatial = O.J.time_step_velocity(
      v_spatial=O.v_spatial, a_spatial=O.a_spatial, delta_t=delta_t)
    O.J = O.J.time_step_position(v_spatial=O.v_spatial, delta_t=delta_t)
    O.energies_and_accelerations_update()

def run_simulation(mersenne_twister, n_time_steps, delta_t):
  sim = simulation(mersenne_twister=mersenne_twister)
  e_pots = flex.double()
  e_kins = flex.double()
  for i_step in xrange(n_time_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin.tot)
  e_tots = e_pots + e_kins
  print "start e_tot:", e_tots[0]
  print "final e_tot:", e_tots[-1]
  print "min:", min(e_tots)
  print "max:", max(e_tots)
  print "ave:", sum(e_tots)/len(e_tots)
  print "range, max, min:", max(e_tots)-min(e_tots), max(e_tots), min(e_tots)
  rr = (max(e_tots)-min(e_tots))/(sum(e_tots)/len(e_tots))
  print "relative range:", rr
  print "e_pot min, max:", min(e_pots), max(e_pots)
  print "e_kin min, max:", min(e_kins), max(e_kins)
  print
  return rr

def exercise(n_trials=10, n_time_steps=1000, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  relative_ranges = flex.double()
  for i in xrange(n_trials):
    relative_ranges.append(run_simulation(
      mersenne_twister=mersenne_twister,
      n_time_steps=n_time_steps,
      delta_t=delta_t))
  print "relative ranges:"
  relative_ranges.min_max_mean().show(prefix="  ")
  print

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
