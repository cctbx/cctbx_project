from scitbx.rigid_body_dynamics import free_motion_reference_impl as fmri
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times, null_out
import sys

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 1
    out = null_out()
  else:
    n_dynamics_steps = max(1, int(args[0]))
    out = sys.stdout
  #
  sites_cart = fmri.create_triangle_with_center_of_mass_at_origin()
  assert approx_equal(flex.vec3_double(sites_cart).mean(), (0,0,0))
  inertia1 = fmri.body_inertia(sites_cart=sites_cart)
  inertia2 = matrix.sym(sym_mat3=scitbx.math.inertia_tensor(
    points=flex.vec3_double(sites_cart), pivot=(0,0,0)))
  assert approx_equal(inertia1, inertia2)
  #
  sim = fmri.simulation()
  assert approx_equal(
    [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot],
    [0.62637925394862359,
     0.012310594130384761, 0.02835, 0.04066059413038476,
     0.6670398480790084])
  for i in xrange(100):
    sim.dynamics_step(delta_t=0.01)
  assert approx_equal(
    [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot],
    [0.018277821171901298,
     0.093940194296481649, 0.55483639080176617, 0.64877658509824787,
     0.66705440627014911])
  #
  sim = fmri.simulation()
  e_tots = flex.double([sim.e_tot])
  print >> out, "i_step, [e_pot, e_kin_ang, e_kin_lin, e_kin, e_tot]"
  def show(i_step):
    print >> out, \
      i_step, [sim.e_pot, sim.e_kin_ang, sim.e_kin_lin, sim.e_kin, sim.e_tot]
    out.flush()
  n_show = max(1, n_dynamics_steps // 10)
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=0.001)
    e_tots.append(sim.e_tot)
    if (i_step % n_show == 0):
      show(i_step)
  show(n_dynamics_steps)
  print >> out
  print >> out, "number of dynamics steps:", n_dynamics_steps
  print >> out, "e_tot start:", e_tots[0]
  print >> out, "      final:", e_tots[-1]
  print >> out, "        min:", flex.min(e_tots)
  print >> out, "        max:", flex.max(e_tots)
  print >> out, "    max-min:", flex.max(e_tots) - flex.min(e_tots)
  print >> out
  out.flush()
  #
  print format_cpu_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
