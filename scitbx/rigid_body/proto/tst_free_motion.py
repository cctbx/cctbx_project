from __future__ import absolute_import, division, print_function
from scitbx.rigid_body.proto import free_motion_reference_impl as fmri
from scitbx.rigid_body.proto import featherstone
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times, null_out
import sys
from six.moves import range

def exercise_reference_impl_quick():
  sites_cart = fmri.create_triangle_with_center_of_mass_at_origin()
  assert approx_equal(flex.vec3_double(sites_cart).mean(), (0,0,0))
  inertia1 = fmri.body_inertia(sites_cart=sites_cart)
  inertia2 = matrix.sym(sym_mat3=scitbx.math.inertia_tensor(
    points=flex.vec3_double(sites_cart), pivot=(0,0,0)))
  assert approx_equal(inertia1, inertia2)
  #
  for use_classical_accel in [False, True]:
    sim = fmri.simulation()
    assert approx_equal(
      [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot],
      [0.64030878777041611,
       0.012310594130384761, 0.02835, 0.04066059413038476,
       0.68096938190080092])
    for i in range(100):
      sim.dynamics_step(delta_t=0.01, use_classical_accel=use_classical_accel)
    expected = [
      [0.028505221929112364,
       0.091503230553568404, 0.56329655444242244, 0.65479978499599079,
       0.6833050069251031],
      [0.053276067541032097,
       0.091503230553568404, 0.53805622991666513, 0.62955946047023348,
       0.68283552801126557]][int(use_classical_accel)]
    assert approx_equal(
      [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot],
      expected)

def exercise_reference_impl_long(n_dynamics_steps, out):
  sim = fmri.simulation()
  e_tots = flex.double([sim.e_tot])
  print("i_step, [e_pot, e_kin_ang, e_kin_lin, e_kin, e_tot]", file=out)
  def show(i_step):
    print(i_step, [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot], file=out)
    out.flush()
  n_show = max(1, n_dynamics_steps // 10)
  for i_step in range(n_dynamics_steps):
    sim.dynamics_step(delta_t=0.001)
    e_tots.append(sim.e_tot)
    if (i_step % n_show == 0):
      show(i_step)
  show(n_dynamics_steps)
  print(file=out)
  print("number of dynamics steps:", n_dynamics_steps, file=out)
  print("e_tot start:", e_tots[0], file=out)
  print("      final:", e_tots[-1], file=out)
  print("        min:", flex.min(e_tots), file=out)
  print("        max:", flex.max(e_tots), file=out)
  print("    max-min:", flex.max(e_tots) - flex.min(e_tots), file=out)
  print(file=out)
  out.flush()

class featherstone_system_model(object):

  def __init__(model, m, I, J):
    model.NB = 1
    model.pitch = [J]
    model.parent =[-1]
    model.Xtree = [matrix.identity(n=6)]
    model.I = [featherstone.mcI(m, (0,0,0), I)]

class six_dof_joint_euler_params_featherstone(fmri.six_dof_joint_euler_params):

  def Xj_S(O, q):
    assert q is None
    Xj = featherstone.Xrot(O.E) \
       * featherstone.Xtrans(O.r) # RBDA Tab. 4.1 footnote
    S = None
    return Xj, S

def exercise_featherstone_FDab(out):
  def check():
    model = featherstone_system_model(
      m=sim.m,
      I=sim.I_F1,
      J=six_dof_joint_euler_params_featherstone(qE=sim.J1.qE, qr=sim.J1.qr))
    q = [None] # already stored in J1 as qE and qr
    qd = [sim.qd]
    tau = None
    f_ext = [matrix.col((sim.nc_F1, sim.f_F1)).resolve_partitions()]
    grav_accn = [0,0,0]
    qdd = featherstone.FDab(model, q, qd, tau, f_ext, grav_accn)
    if (i_step % 10 == 0):
      print("ang acc 3D:", sim.wd_F1.elems, file=out)
      print("        6D:", qdd[0].elems[:3], file=out)
      print(file=out)
      print("lin acc 3D:", sim.as_F1.elems, file=out)
      print("        6D:", qdd[0].elems[3:], file=out)
      print(file=out)
    assert approx_equal(qdd[0].elems[:3], sim.wd_F1)
    assert approx_equal(qdd[0].elems[3:], sim.as_F1)
  sim = fmri.simulation()
  for i_step in range(100):
    check()
    sim.dynamics_step(delta_t=0.1) # large time step to sample
  check()                          # diverse configurations

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 1
    out = null_out()
  else:
    n_dynamics_steps = max(1, int(args[0]))
    out = sys.stdout
  #
  exercise_reference_impl_quick()
  exercise_featherstone_FDab(out=out)
  exercise_reference_impl_long(n_dynamics_steps=n_dynamics_steps, out=out)
  #
  print(format_cpu_times())

if (__name__ == "__main__"):
  run(sys.argv[1:])
