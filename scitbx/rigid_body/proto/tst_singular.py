from scitbx.rigid_body.proto.tst_joint_lib import exercise_sim
from scitbx.rigid_body.proto import joint_lib
from scitbx.rigid_body.proto.test_simulation import simulation
from scitbx.rigid_body.proto.test_utils import create_wells
from scitbx.rigid_body.proto.utils import spatial_inertia_from_sites
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
import math
import sys

class six_dof_body(object):

  def __init__(O, mersenne_twister, n_sites):
    if (n_sites > 0):
      O.sites = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    while (len(O.sites) != n_sites):
      O.sites.append(O.sites[0]
        + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.six_dof_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    O.J = joint_lib.six_dof(type="euler_params", qE=qE, qr=qr, r_is_qr=True)
    O.qd = matrix.col(mersenne_twister.random_double(size=6)*2-1)

class spherical_body(object):

  def __init__(O, mersenne_twister, n_sites):
    if (n_sites > 0):
      O.sites = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    while (len(O.sites) != n_sites):
      O.sites.append(O.sites[0]
        + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.spherical_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    O.J = joint_lib.spherical(type="euler_params", qE=qE)
    O.qd = matrix.col(mersenne_twister.random_double(size=3)*2-1)

class revolute_body(object):

  def __init__(O, mersenne_twister, prev=None):
    if (prev is None):
      normal = matrix.col(mersenne_twister.random_double_point_on_sphere())
      pivot = matrix.col(mersenne_twister.random_double_point_on_sphere()) * 0.6
    else:
      normal = prev.A.normal
      pivot = prev.A.pivot + normal * 0.3
    O.sites = [pivot + normal * (-0.8)]
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    O.J = joint_lib.revolute(qE=matrix.col([math.pi/8]))
    O.qd = matrix.col([math.pi/12])

def exercise_six_dof(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,4):
    for i_trial in xrange(n_trials):
      body = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body.parent = -1
      sim = simulation(bodies=[body])
      print >> out, "six_dof number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is not sys.stdout):
        assert relative_range < 1.e-4

def exercise_six_dof2(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,4):
    for i_trial in xrange(n_trials):
      body1 = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body1.parent = -1
      body2 = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body2.parent = 0
      sim = simulation(bodies=[body1, body2])
      print >> out, "six_dof2 number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is not sys.stdout):
        assert relative_range < 1.e-4

def exercise_revolute(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    body = revolute_body(mersenne_twister=mersenne_twister)
    body.parent = -1
    sim = simulation(bodies=[body])
    print >> out, "revolute:"
    relative_range = exercise_sim(
      out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
    if (out is not sys.stdout):
      assert relative_range < 1.e-4

def exercise_revolute2(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    body1 = revolute_body(mersenne_twister=mersenne_twister)
    body1.parent = -1
    body2 = revolute_body(mersenne_twister=mersenne_twister, prev=body1)
    body2.parent = 0
    sim = simulation(bodies=[body1, body2])
    print >> out, "revolute2:"
    relative_range = exercise_sim(
      out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
    if (out is not sys.stdout):
      assert relative_range < 1.e-4

def exercise_spherical(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,3):
    for i_trial in xrange(n_trials):
      body = spherical_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body.parent = -1
      sim = simulation(bodies=[body])
      print >> out, "spherical number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is not sys.stdout):
        assert relative_range < 1.e-4

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
  if (1):
    exercise_six_dof(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_six_dof2(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_spherical(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_revolute(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_revolute2(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
