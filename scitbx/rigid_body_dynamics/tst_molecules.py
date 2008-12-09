from scitbx.rigid_body_dynamics.tst_joint_lib import exercise_sim
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.test_simulation import simulation
from scitbx.rigid_body_dynamics.test_utils import create_wells
from scitbx.rigid_body_dynamics.utils import spatial_inertia_from_sites
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
import sys

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
    O.A = joint_lib.six_dof_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = shift_gently(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col((0,0,0))
    qr = matrix.col((0,0,0))
    O.J = joint_lib.six_dof(type="euler_params", qE=qE, qr=qr, r_is_qr=True)
    O.qd = O.J.qd_zero

class revolute_body(object):

  def __init__(O, labels, sites, bonds, pivot, normal, mersenne_twister):
    O.labels = labels
    O.sites = sites
    O.bonds = bonds
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = shift_gently(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    O.J = joint_lib.revolute(qE=matrix.col([0]))
    O.qd = O.J.qd_zero

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
    pivot=sites[3],
    normal=(sites[3]-sites[2]).normalize(),
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
    pivot=sites[3],
    normal=(sites[3]-sites[2]).normalize(),
    mersenne_twister=mersenne_twister)
  body1.parent = 0
  body2 = revolute_body(
    labels=labels[4:],
    sites=sites[4:],
    bonds=[(-3,0)],
    pivot=sites[4],
    normal=(sites[4]-sites[0]).normalize(),
    mersenne_twister=mersenne_twister)
  body2.parent = 0
  return simulation(bodies=[body0, body1, body2])

simulation_factories = [
  simulation_gly_no_h,
  simulation_gly_with_nh]

def exercise_dynamics_quick(out, sim, n_dynamics_steps, delta_t=0.001):
  relative_range = exercise_sim(
    out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
  if (out is not sys.stdout):
    assert relative_range < 1.e-4

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
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
