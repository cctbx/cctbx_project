"""\
6-DoF joint (free motion) reference implementation using the notation of
Featherstone (2007).

Dynamics simulation of the motion of a rigid triangle attached with
springs to an arbitrary position in space.

RBDA:
  Rigid Body Dynamics Algorithms.
  Roy Featherstone,
  Springer, New York, 2007.
  ISBN-10: 0387743146

Shabana (2005):
  Dynamics of Multibody Systems
  Ahmed A. Shabana
  Cambridge University Press, 3 edition, 2005
  ISBN-10: 0521850118
"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip

try:
  from scitbx import matrix
except ImportError:
  import scitbx_matrix as matrix

class simulation(object):

  def __init__(O):

    sites_cart_F0 = create_triangle_with_center_of_mass_at_origin()
    O.sites_cart_F1 = sites_cart_F0 # Xtree1, XJ1 = identity
    # F0, F01, Xtree1, J1: see RBDA Fig. 4.7

    # body inertia in F01, assuming unit masses, and total mass
    #
    O.I_F1 = body_inertia(O.sites_cart_F1)
    O.m = len(O.sites_cart_F1)

    # Euler parameters for some arbitrary rotation and translation
    #
    qE_J1 = matrix.col((0.9, -0.14, 0.1, 0.06)).normalize() # RBDA Eq. 4.11
    qr_J1 = matrix.col((-0.15, 0.1, 0.05))
    O.J1 = six_dof_joint_euler_params(qE_J1, qr_J1)

    # arbitrary spatial velocity in F01
    #
    w_F1 = matrix.col((0.04, -0.15, 0.10))
    v_F1 = matrix.col((0.1, 0.08, -0.05))
    O.qd = matrix.col((w_F1, v_F1)).resolve_partitions() # RBDA Eq. 2.4

    # arbitrary rotation and translation of potential wells
    #
    qE_Jp = matrix.col((0.99, -0.05, 0.04, -0.11)).normalize()
    qr_Jp = matrix.col((-0.15, 0.14, -0.37))
    Jp = six_dof_joint_euler_params(qE_Jp, qr_Jp)
    O.sites_cart_wells_F01 = [Jp.Tsp * xyz for xyz in O.sites_cart_F1]

    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):

    # positional coordinates of moved triangle in F01 (not actually used here)
    #
    O.sites_cart_moved_F01 = [O.J1.Tsp * xyz for xyz in O.sites_cart_F1]

    # potential pulling triangle to wells (spring-like forces)
    #
    sites_cart_wells_moved_F1 = [O.J1.Tps * xyz
      for xyz in O.sites_cart_wells_F01]
    O.e_pot = 0
    for xyz, xyz_wells_moved in zip(O.sites_cart_F1,
                                    sites_cart_wells_moved_F1):
      O.e_pot += (xyz - xyz_wells_moved).dot()

    # corresponding forces
    #
    f_cart_F1 = [-2 * (xyz - xyz_wells_moved)
      for xyz, xyz_wells_moved in zip(O.sites_cart_F1,
                                      sites_cart_wells_moved_F1)]

    # f and nc in 3D, for RBDA Eq. 1.4
    #
    O.f_F1 = matrix.col((0,0,0))
    O.nc_F1 = matrix.col((0,0,0))
    for xyz,force in zip(O.sites_cart_F1, f_cart_F1):
      O.f_F1 += force
      O.nc_F1 += xyz.cross(force)

    # solution of RBDA Eq. 1.4 for ac and wd
    #
    O.ac_F1 = O.f_F1 / O.m
    w_F1, v_F1 = matrix.col(O.qd.elems[:3]), matrix.col(O.qd.elems[3:])
    O.wd_F1 = O.I_F1.inverse() * (O.nc_F1 - w_F1.cross(O.I_F1 * w_F1))

    # spatial acceleration
    O.as_F1 = O.ac_F1 - w_F1.cross(v_F1) # RBDA Eq. 2.47

    # Kinetic energy, angular velocity (Shabana (2005) p. 148 eq. 3.126)
    #
    O.e_kin_ang = 0.5 * w_F1.dot(O.I_F1 * w_F1)

    # Kinetic energy, linear velocity (Shabana (2005) p. 148 eq. 3.125)
    #
    O.e_kin_lin = 0.5 * O.m * v_F1.dot()

    O.e_kin = O.e_kin_ang + O.e_kin_lin
    O.e_tot = O.e_pot + O.e_kin

  def dynamics_step(O, delta_t, use_classical_accel=False):
    if (use_classical_accel):
      qdd = matrix.col((O.wd_F1, O.ac_F1)).resolve_partitions()
    else:
      qdd = matrix.col((O.wd_F1, O.as_F1)).resolve_partitions()
    O.qd = O.J1.time_step_velocity(O.qd, qdd, delta_t)
    O.J1 = O.J1.time_step_position(O.qd, delta_t)
    O.energies_and_accelerations_update()

class six_dof_joint_euler_params(object):

  def __init__(O, qE, qr):
    O.qE = qE
    O.qr = qr
    #
    O.E = RBDA_Eq_4_12(qE)
    O.r = O.E.transpose() * qr # RBDA Tab. 4.1
    #
    O.Tps = matrix.rt((O.E, -O.E * O.r)) # RBDA Eq. 2.28
    O.Tsp = matrix.rt((O.E.transpose(), O.r))
    # RBDA p. 69, Sec. 4.1.3:
    #   s = frame fixed in successor
    #   p = frame fixed in predecessor

  def time_step_position(O, qd, delta_t):
    w_F1, v_F1 = matrix.col_list([qd.elems[:3], qd.elems[3:]])
    qEd = RBDA_Eq_4_13(O.qE.elems) * w_F1
    qrd = v_F1 - w_F1.cross(O.qr) # RBDA Eq. 2.38 p. 27
    new_qE = (O.qE + qEd * delta_t).normalize() # RBDA, bottom of p. 86
    new_qr = O.qr + qrd * delta_t
    return six_dof_joint_euler_params(new_qE, new_qr)

  def time_step_velocity(O, qd, qdd, delta_t):
    return qd + qdd * delta_t

def RBDA_Eq_4_12(q):
  p0, p1, p2, p3 = q
  return matrix.sqr((
    p0**2+p1**2-0.5,   p1*p2+p0*p3,     p1*p3-p0*p2,
      p1*p2-p0*p3,   p0**2+p2**2-0.5,   p2*p3+p0*p1,
      p1*p3+p0*p2,     p2*p3-p0*p1,   p0**2+p3**2-0.5)) * 2

def RBDA_Eq_4_13(q):
  p0, p1, p2, p3 = q
  return matrix.rec((
    -p1, -p2, -p3,
    p0, -p3, p2,
    p3, p0, -p1,
    -p2, p1, p0), n=(4,3)) * 0.5

def body_inertia(sites_cart):
  m = [0] * 9
  for x,y,z in sites_cart:
    m[0] += y*y+z*z
    m[4] += x*x+z*z
    m[8] += x*x+y*y
    m[1] -= x*y
    m[2] -= x*z
    m[5] -= y*z
  m[3] = m[1]
  m[6] = m[2]
  m[7] = m[5]
  return matrix.sqr(m)

def create_triangle_with_center_of_mass_at_origin():
  sites_cart = matrix.col_list([
    (0,0,0),
    (1,0,0),
    (0.5,0.5*3**0.5,0)])
  com = matrix.col((0,0,0))
  for xyz in sites_cart:
    com += xyz
  com /= len(sites_cart)
  sites_cart = [xyz-com for xyz in sites_cart]
  E = RBDA_Eq_4_12(matrix.col((0.81, -0.32, -0.42, -0.26)).normalize())
  sites_cart = [E*xyz for xyz in sites_cart]
  return sites_cart

def run():
  O = simulation()
  for i_time_step in range(10):
    print("e_kin tot ang lin:", O.e_kin, O.e_kin_ang, O.e_kin_lin)
    print("            e_pot:", O.e_pot)
    print("            e_tot:", O.e_tot)
    print("ang acc 3D:", O.wd_F1.elems)
    print("lin acc 3D:", O.as_F1.elems)
    print()
    O.dynamics_step(delta_t=0.01)
  print("OK")

if (__name__ == "__main__"):
  run()
