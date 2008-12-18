from scitbx.rigid_body.proto import featherstone
from scitbx.rigid_body.proto import test_utils
from scitbx.rigid_body.proto.utils import \
  e_kin_from_model, \
  featherstone_system_model
import scitbx.lbfgs
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal

class simulation(object):

  def __init__(O, bodies):
    O.bodies = bodies
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
    O.qdd = featherstone.FDab(model, q, qd, tau, O.f_ext_bf, grav_accn)

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
      e_pot_bf = test_utils.potential_energy_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      f_ext_using_bf = test_utils.potential_f_ext_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      O.f_ext_bf.append(f_ext_using_bf)
      O.e_pot += e_pot_bf
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def sensitivity_test(O, n_significant_digits):
    "RBDA section 10.2, p. 199-201"
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.qd for B in O.bodies]
    qdd = [matrix.col([1]*len(B.qd)) for B in O.bodies]
    grav_accn = [0,0,0]
    tau = featherstone.ID(model, q, qd, qdd, O.f_ext_bf, grav_accn)
    if (n_significant_digits is not None):
      assert n_significant_digits > 0
      fmt = "%%.%dg" % n_significant_digits
      tau_trunc = []
      for v in tau:
        tau_trunc.append(matrix.col([float(fmt%e) for e in v]))
      tau = tau_trunc
    qdd = featherstone.FDab(model, q, qd, tau, O.f_ext_bf, grav_accn)
    result = []
    for v in qdd: result.extend(v.elems)
    return result

  def d_pot_d_q(O):
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.J.qd_zero for B in O.bodies]
    qdd = [B.J.qdd_zero for B in O.bodies]
    grav_accn = [0,0,0]
    taus = featherstone.ID(model, q, qd, qdd, O.f_ext_bf, grav_accn)
    result = []
    for B,tau in zip(O.bodies, taus):
      tau_as_d_pot_d_q = getattr(B.J, "tau_as_d_pot_d_q", None)
      if (tau_as_d_pot_d_q is None):
        result.append(tau)
      else:
        result.append(tau_as_d_pot_d_q(tau=tau))
    return result

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
