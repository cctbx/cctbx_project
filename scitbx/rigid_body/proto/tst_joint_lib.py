from __future__ import absolute_import, division, print_function
from scitbx.rigid_body.proto import featherstone
from scitbx.rigid_body.proto import joint_lib
from scitbx.rigid_body.proto.utils import \
  spatial_inertia_from_sites, \
  T_as_X, \
  featherstone_system_model, \
  e_kin_from_model
from scitbx.rigid_body.proto import test_utils
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
import math
import sys
from six.moves import range
from six.moves import zip

class random_revolute(object):

  def __init__(O, mersenne_twister):
    def random_vector():
      return matrix.col(mersenne_twister.random_double(size=3)*2-1)
    def random_angle():
      return (mersenne_twister.random_double()*2-1)*math.pi
    #
    O.sites = [random_vector()]
    for i_trial in range(100): # guard against unlikely singularity
      O.A = joint_lib.revolute_alignment(
        pivot=random_vector(),
        normal=random_vector().normalize())
      if (abs(O.A.normal.cos_angle(O.sites[0] - O.A.pivot)) > 1.e-3):
        break
    else:
      raise RuntimeError
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = [random_vector()]
    #
    O.J = joint_lib.revolute(qE=matrix.col([random_angle()]))
    O.qd = matrix.col([random_angle()])

class revolute_z(object):

  def __init__(O, pivot, x=1/3.):
    O.sites = [matrix.col((x,0,2/3.))+pivot]
    O.A = joint_lib.revolute_alignment(
      pivot=pivot,
      normal=matrix.col((0,0,1)))
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    O.wells = O.sites
    O.J = joint_lib.revolute(qE=matrix.col([math.pi/8.]))
    O.qd = matrix.col([0])

class revolute_x(object):

  def __init__(O, pivot):
    O.sites = [matrix.col((2/3.,0,1/3.))+pivot]
    O.A = joint_lib.revolute_alignment(
      pivot=pivot,
      normal=matrix.col((1,0,0)))
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    O.wells = O.sites
    O.J = joint_lib.revolute(qE=matrix.col([math.pi/6.]))
    O.qd = matrix.col([0])

class revolute_simulation(object):

  def __init__(O, mersenne_twister, NB, config):
    O.bodies = []
    if (config == "random"):
      for ib in range(NB):
        B = random_revolute(mersenne_twister=mersenne_twister)
        B.parent = -1+ib
        O.bodies.append(B)
    elif (config == "zigzag"):
      assert NB <= 3
      B = revolute_z(pivot=matrix.col((0,0,0)))
      B.parent = -1
      O.bodies.append(B)
      if (NB > 1):
        B = revolute_x(pivot=matrix.col((0,0,1)))
        B.parent = 0
        O.bodies.append(B)
      if (NB > 2):
        B = revolute_z(pivot=matrix.col((1,0,1)))
        B.parent = 1
        O.bodies.append(B)
    elif (config == "singular"):
      assert NB == 1
      B = revolute_z(pivot=matrix.col((0,0,0)), x=0)
      B.parent = -1
      O.bodies.append(B)
    else:
      raise RuntimeError
    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.qd for B in O.bodies]
    #
    O.e_kin = e_kin_from_model(model, q, qd)
    O.e_pot_and_f_ext_update()
    #
    tau = None
    grav_accn = [0,0,0]
    qdd_using_f_ext_ff = featherstone.FDab(
      model, q, qd, tau, O.f_ext_ff, grav_accn, f_ext_in_ff=True)
    qdd_using_f_ext_bf = featherstone.FDab(
      model, q, qd, tau, O.f_ext_bf, grav_accn, f_ext_in_ff=False)
    assert approx_equal(qdd_using_f_ext_bf, qdd_using_f_ext_ff)
    O.qdd = qdd_using_f_ext_ff
    #
    X0s = FDab_X0(model, q, qd)
    e_pot_vfy = check_transformations(O.bodies, model.Ttree, X0s)
    assert approx_equal(e_pot_vfy, O.e_pot)

  def e_pot_and_f_ext_update(O):
    O.AJA_accu = []
    O.e_pot = 0
    O.f_ext_ff = []
    O.f_ext_bf = []
    for B in O.bodies:
      AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
      if (B.parent == -1):
        AJA_tree = None
      else:
        AJA_tree = O.AJA_accu[B.parent]
        AJA = AJA_tree * AJA
      O.AJA_accu.append(AJA)
      e_pot_ff = test_utils.potential_energy(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      e_pot_bf = test_utils.potential_energy_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      assert approx_equal(e_pot_bf, e_pot_ff)
      f_ext_using_ff = test_utils.potential_f_ext_ff(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      f_ext_using_bf = test_utils.potential_f_ext_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      O.f_ext_ff.append(f_ext_using_ff)
      O.f_ext_bf.append(f_ext_using_bf)
      O.e_pot += e_pot_ff
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def d_pot_d_q(O):
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [matrix.col((0,)) for B in O.bodies]
    qdd = [matrix.col((0,)) for B in O.bodies]
    grav_accn = [0,0,0]
    return featherstone.ID(model, q, qd, qdd, O.f_ext_bf, grav_accn)

  def d_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for B in O.bodies:
      fs = []
      J_orig = B.J
      for signed_eps in [eps, -eps]:
        B.J = joint_lib.revolute(qE=matrix.col((J_orig.qE[0]+signed_eps,)))
        O.e_pot_and_f_ext_update()
        fs.append(O.e_pot)
      B.J = J_orig
      result.append(matrix.col(((fs[0]-fs[1])/(2*eps),)))
    O.energies_and_accelerations_update()
    return result

  def check_d_pot_d_q(O):
    qdd_orig = O.qdd
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd, qdd_orig)

def FDab_X0(model, q, qd):
  Xup = [None] * model.NB
  X0 = [None] * model.NB
  for i in range(model.NB):
    XJ, S = featherstone.jcalc( model.pitch[i], q[i] )
    Xup[i] = XJ * model.Xtree[i]
    if model.parent[i] == -1:
      X0[i] = Xup[i]
    else:
      X0[i] = Xup[i] * X0[model.parent[i]]
  return X0

def check_transformations(bodies, Ttree, X0s):
  T0s = []
  for B,Tt,X0 in zip(bodies, Ttree, X0s):
    Tj = B.J.Tps
    Tup = Tj * Tt
    if (B.parent == -1):
      T0s.append(Tup)
    else:
      T0s.append(Tup * T0s[B.parent])
    X0_from_T0 = T_as_X(T0s[-1])
    assert approx_equal(X0_from_T0, X0)
  e_pot = 0
  AJA_accu = []
  for B,T0 in zip(bodies, T0s):
    AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
    if (B.parent == -1):
      AJA_tree = None
      AJA_accu.append(AJA)
    else:
      AJA_tree = AJA_accu[B.parent]
      AJA = AJA_tree * AJA
      AJA_accu.append(AJA)
    for s,w in zip(B.sites, B.wells):
      s_bf = B.A.T0b * s
      s_mv1 = T0.inverse_assuming_orthogonal_r() * s_bf
      s_mv2 = AJA * s
      assert approx_equal(s_mv1, s_mv2)
      e_pot += (s_mv1 - w).dot()
  return e_pot

plot_number = [0]

def exercise_sim(out, n_dynamics_steps, delta_t, sim):
  sim.check_d_pot_d_q()
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in range(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  sim.check_d_pot_d_q()
  print("energy samples:", e_tots.size(), file=out)
  print("e_pot min, max:", min(e_pots), max(e_pots), file=out)
  print("e_kin min, max:", min(e_kins), max(e_kins), file=out)
  print("e_tot min, max:", min(e_tots), max(e_tots), file=out)
  print("start e_tot:", e_tots[0], file=out)
  print("final e_tot:", e_tots[-1], file=out)
  ave = flex.sum(e_tots) / e_tots.size()
  range_ = flex.max(e_tots) - flex.min(e_tots)
  if (ave == 0): relative_range = 0
  else:          relative_range = range_ / ave
  print("ave:", ave, file=out)
  print("range:", range_, file=out)
  print("relative range:", relative_range, file=out)
  print(file=out)
  out.flush()
  if (out is sys.stdout):
    f = open("tmp%02d.xy" % plot_number[0], "w")
    for es in [e_pots, e_kins, e_tots]:
      for e in es:
        print(e, file=f)
      print("&", file=f)
    f.close()
    plot_number[0] += 1
  return relative_range

def exercise_revolute_sim(
      out,
      mersenne_twister,
      n_dynamics_steps,
      delta_t,
      NB,
      config):
  print("config:", config, file=out)
  sim = revolute_simulation(
    mersenne_twister=mersenne_twister,
    NB=NB,
    config=config)
  return exercise_sim(
    out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)

def exercise_revolute(out, n_trials, n_dynamics_steps, delta_t=0.001, NB=3):
  mersenne_twister = flex.mersenne_twister(seed=0)
  relative_ranges = flex.double()
  for i_trial in range(n_trials):
    relative_ranges.append(exercise_revolute_sim(
      out=out,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t,
      NB=[1, NB][min(i_trial, 1)],
      config=["singular", "zigzag", "random"][min(i_trial, 2)]))
  print("relative ranges:", file=out)
  relative_ranges.min_max_mean().show(out=out, prefix="  ")
  if (out is not sys.stdout):
    assert flex.max(relative_ranges) < 0.0006
  print(file=out)

def run(args):
  assert len(args) in [0,2]
  if (len(args) == 0):
    n_trials = 3
    n_dynamics_steps = 30
    out = null_out()
  else:
    n_trials = max(1, int(args[0]))
    n_dynamics_steps = max(1, int(args[1]))
    out = sys.stdout
  show_times_at_exit()
  exercise_revolute(
    out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
