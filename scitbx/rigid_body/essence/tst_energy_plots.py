from scitbx.rigid_body.essence import tst_tardy
from scitbx.array_family import flex
import sys

def run(args):
  if (len(args) == 0):
    model_index = 1
    n_dynamics_steps = 100
    delta_t = 0.01
  else:
    assert len(args) == 3, "model_index, n_dynamics_steps, delta_t"
    model_index = int(eval(args[0]))
    n_dynamics_steps = int(eval(args[1]))
    delta_t = float(eval(args[2]))
  tardy_model = tst_tardy.get_test_model_by_index(i=model_index)
  tardy_model.assign_random_velocities(e_kin_target=1)
  e_pots = flex.double([tardy_model.e_pot()])
  e_kins = flex.double([tardy_model.e_kin()])
  def show_e_tot():
    print "e_tot: %.6g" % (e_pots[-1]+e_kins[-1])
    sys.stdout.flush()
  for i_step in xrange(n_dynamics_steps):
    tardy_model.dynamics_step(delta_t=delta_t)
    e_pots.append(tardy_model.e_pot())
    e_kins.append(tardy_model.e_kin())
    if (i_step % 1000 == 0):
      print "i_step:", i_step,
      show_e_tot()
  print "i_step:", n_dynamics_steps,
  show_e_tot()
  e_tots = e_pots + e_kins
  f = open("e_pot_kin_tot_i=%02d_n=%d_d=%.5f.xy" % (
    model_index, n_dynamics_steps, delta_t), "w")
  for es in [e_pots, e_kins, e_tots]:
    for e in es:
      print >> f, e
    print >> f, "&"
  f.close()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
