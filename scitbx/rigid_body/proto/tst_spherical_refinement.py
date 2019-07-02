from __future__ import absolute_import, division, print_function
from scitbx.rigid_body.proto import joint_lib
from scitbx.rigid_body.proto import test_utils
from scitbx.rigid_body.proto.utils import center_of_mass_from_sites
import scitbx.lbfgs
import scitbx.math
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
from libtbx.str_utils import show_sorted_by_counts
import math
import sys
from six.moves import range

class setT_mixin(object):

  def setT(O):
    O.Tps = matrix.rt((O.E, (0,0,0)))
    O.Tsp = matrix.rt((O.E.transpose(), (0,0,0)))

class euler_params(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((1,0,0,0))
    O.qE = qE
    O.unit_quaternion = qE.normalize()
    O.E = joint_lib.RBDA_Eq_4_12(q=O.unit_quaternion)
    O.setT()

  def tau_as_d_pot_d_q(O, tau):
    d = joint_lib.d_unit_quaternion_d_qE_matrix(q=O.qE)
    c = d * 4 * joint_lib.RBDA_Eq_4_13(q=O.unit_quaternion)
    n = tau
    return c * n

class euler_angles_xyz(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0))
    O.qE = qE
    O.E = joint_lib.RBDA_Eq_4_7(q=qE)
    O.setT()

  def tau_as_d_pot_d_q(O, tau):
    c = joint_lib.RBDA_Eq_4_8(q=O.qE).transpose()
    n = tau
    return c * n

class euler_angles_zxz(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0))
    O.qE = qE
    O.E = euler_angles_zxz_matrix(q=qE)
    O.setT()

  def tau_as_d_pot_d_q(O, tau):
    c = euler_angles_zxz_S(q=O.qE).transpose()
    n = tau
    return c * n

class euler_angles_yxyz(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0,0))
    O.qE = qE
    O.E = euler_angles_yxyz_matrix(q=qE)
    O.setT()

class euler_angles_xyzy(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0,0))
    O.qE = qE
    O.E = euler_angles_xyzy_matrix(q=qE)
    O.setT()

class inf_euler_params(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0))
    O.qE = qE
    th = abs(qE)
    if (th == 0):
      p = matrix.col((1,0,0,0))
    else:
      p0 = math.cos(th)
      p1, p2, p3 = math.sin(th) / th * qE
      p = matrix.col((p0,p1,p2,p3))
      assert abs(abs(p)-1) < 1.e-12
    O.E = joint_lib.RBDA_Eq_4_12(q=p)
    O.setT()

class inf_axis_angle(setT_mixin):

  def __init__(O, qE):
    if (qE is None): qE = matrix.col((0,0,0))
    O.qE = qE
    angle = abs(qE)
    if (angle == 0):
      O.E = matrix.sqr((1,0,0,0,1,0,0,0,1))
    else:
      axis = qE / angle
      O.E = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis=axis, angle=angle, deg=False))
    O.setT()

def euler_angles_zxz_matrix(q):
  q1,q2,q3 = q
  def cs(a): return math.cos(a), math.sin(a)
  c1,s1 = cs(q1)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  """
  Mathematica code:
    rz1 = {{c1, s1, 0}, {-s1, c1, 0}, {0, 0, 1}}
    rx2 = {{1, 0, 0}, {0, c2, s2}, {0, -s2, c2}}
    rz3 = {{c3, s3, 0}, {-s3, c3, 0}, {0, 0, 1}}
    rz3.rx2.rz1
    Same as Goldstein Eq. 4.46
  """
  return matrix.sqr((
    c1*c3-c2*s1*s3, c3*s1+c1*c2*s3, s2*s3,
    -c2*c3*s1-c1*s3, c1*c2*c3-s1*s3, c3*s2,
    s1*s2, -c1*s2, c2))

def euler_angles_zxz_S(q):
  q1,q2,q3 = q
  def cs(a): return math.cos(a), math.sin(a)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  # see RBDA p. 85, text after Eq. 4.8
  return matrix.sqr((
    s2*s3,  c3, 0,
    c3*s2, -s3, 0,
       c2,   0, 1))

def euler_angles_yxyz_matrix(q):
  q1,q2,q3,q4 = q
  def cs(a): return math.cos(a), math.sin(a)
  c1,s1 = cs(q1)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  c4,s4 = cs(q4)
  """
  Mathematica code:
    rz1 = {{c1, s1, 0}, {-s1, c1, 0}, {0, 0, 1}}
    ry2 = {{c2, 0, -s2}, {0, 1, 0}, {s2, 0, c2}}
    rx3 = {{1, 0, 0}, {0, c3, s3}, {0, -s3, c3}}
    ry4 = {{c4, 0, -s4}, {0, 1, 0}, {s4, 0, c4}}
    ry4.rx3.ry2.rz1
  """
  return matrix.sqr((
    -(s1*s3*s4)+c1*(c2*c4-c3*s2*s4),
    c1*s3*s4+s1*(c2*c4-c3*s2*s4),
    -(c4*s2)-c2*c3*s4,
    -(c3*s1)+c1*s2*s3,
    c1*c3+s1*s2*s3,
    c2*s3,
    c4*s1*s3+c1*(c3*c4*s2+c2*s4),
    -(c1*c4*s3)+s1*(c3*c4*s2+c2*s4),
    c2*c3*c4-s2*s4))

def euler_angles_xyzy_matrix(q):
  q1,q2,q3,q4 = q
  def cs(a): return math.cos(a), math.sin(a)
  c1,s1 = cs(q1)
  c2,s2 = cs(q2)
  c3,s3 = cs(q3)
  c4,s4 = cs(q4)
  """
  Mathematica code:
    ry1 = {{c1, 0, -s1}, {0, 1, 0}, {s1, 0, c1}}
    rz2 = {{c2, s2, 0}, {-s2, c2, 0}, {0, 0, 1}}
    ry3 = {{c3, 0, -s3}, {0, 1, 0}, {s3, 0, c3}}
    rx4 = {{1, 0, 0}, {0, c4, s4}, {0, -s4, c4}}
    rx4.ry3.rz2.ry1
  """
  return matrix.sqr((
    c1*c2*c3-s1*s3,
    c3*s2,
    -(c2*c3*s1)-c1*s3,
    c3*s1*s4+c1*(-(c4*s2)+c2*s3*s4),
    c2*c4+s2*s3*s4,
    c1*c3*s4-s1*(-(c4*s2)+c2*s3*s4),
    c3*c4*s1+c1*(c2*c4*s3+s2*s4),
    c4*s2*s3-c2*s4,
    c1*c3*c4-s1*(c2*c4*s3+s2*s4)))

class body_with_spherical_joint(object):

  def __init__(O, sites, wells, spherical_type):
    O.sites = sites
    O.A = joint_lib.spherical_alignment(sites=O.sites)
    O.wells = wells
    O.spherical_type = spherical_type
    O.J = spherical_type(qE=None)

  def d_pot_d_q_via_finite_differences(O, J=None, eps=1.e-6):
    if (J is None): J = O.J
    result = []
    q = J.qE
    for i in range(len(q)):
      fs = []
      for signed_eps in [eps, -eps]:
        q_eps = list(q)
        q_eps[i] += signed_eps
        q_eps = matrix.col(q_eps)
        J = O.spherical_type(qE=q_eps)
        e_pot = test_utils.potential_energy(
          sites=O.sites, wells=O.wells, A=O.A, J=J)
        fs.append(e_pot)
      result.append((fs[0]-fs[1])/(2*eps))
    return matrix.col(result)

  def d_pot_d_q(O, J=None):
    if (J is None): J = O.J
    f_ext_bf = test_utils.potential_f_ext_bf(
      sites=O.sites, wells=O.wells, A=O.A, J=J)
    tau = matrix.col((-f_ext_bf).elems[:3])
    return J.tau_as_d_pot_d_q(tau=tau)

  def compute_functional_and_gradients(O, qE, use_analytical_gradients):
    J = O.spherical_type(qE=qE)
    f = test_utils.potential_energy(
      sites=O.sites, wells=O.wells, A=O.A, J=J)
    fin = O.d_pot_d_q_via_finite_differences(J=J)
    if (use_analytical_gradients):
      ana = O.d_pot_d_q(J=J)
      assert approx_equal(ana, fin, eps=1.e-4)
      g = ana
    else:
      g = fin
    return f, g

class refinery(object):

  def __init__(O, spherical_type_info, sites, wells, out):
    O.spherical_type_info = spherical_type_info
    O.body = body_with_spherical_joint(
      spherical_type=spherical_type_info.type, sites=sites, wells=wells)
    O.x = flex.double(O.body.J.qE)
    minimizer = scitbx.lbfgs.run(target_evaluator=O)
    f, g = O.compute_functional_and_gradients()
    label = str(O.spherical_type_info)
    print(label, \
      "functional: %12.6g" % f, \
      "gradient norm: %12.6g" % g.norm(), \
      "iter: %3d" % minimizer.iter(), \
      "nfun: %3d" % minimizer.nfun(), file=out)
    O.failed = (f > 1.e-3)
    if (O.failed): print("               FAILED", label, f, file=out)
    sys.stdout.flush()
    O.nfun = minimizer.nfun()

  def compute_functional_and_gradients(O):
    f, g = O.body.compute_functional_and_gradients(
      qE=matrix.col(O.x),
      use_analytical_gradients=O.spherical_type_info.use_analytical_gradients)
    return f, flex.double(g)

def run(args):
  assert len(args) in [0,2], "n_sites, n_trials"
  if (len(args) == 0):
    n_sites, n_trials = 3, 2
    out = null_out()
  else:
    n_sites, n_trials = [int(arg) for arg in args]
    out = sys.stdout
  #
  show_times_at_exit()
  class type_info(object):
    def __init__(O, type, use_analytical_gradients):
      O.type = type
      O.use_analytical_gradients = use_analytical_gradients
    def __str__(O):
      return "%s(use_analytical_gradients=%s)" % (
        O.type.__name__, str(O.use_analytical_gradients))
  spherical_types = [
    type_info(euler_params, False),
    type_info(euler_params, True),
    type_info(euler_angles_xyz, False),
    type_info(euler_angles_xyz, True),
    type_info(euler_angles_zxz, False),
    type_info(euler_angles_zxz, True),
    type_info(euler_angles_yxyz, False),
    type_info(euler_angles_xyzy, False),
    type_info(inf_euler_params, False),
    type_info(inf_axis_angle, False)]
  nfun_accu = {}
  n_failed = {}
  for ti in spherical_types:
    nfun_accu[str(ti)] = flex.size_t()
    n_failed[str(ti)] = 0
  mersenne_twister = flex.mersenne_twister(seed=0)
  for i_trial in range(n_trials):
    sites = [matrix.col(s) for s in flex.vec3_double(
      mersenne_twister.random_double(size=n_sites*3)*2-1)]
    c = center_of_mass_from_sites(sites)
    r = matrix.sqr(mersenne_twister.random_double_r3_rotation_matrix())
    wells = [r*(s-c)+c for s in sites]
    for ti in spherical_types:
      r = refinery(spherical_type_info=ti, sites=sites, wells=wells, out=out)
      nfun_accu[str(ti)].append(r.nfun)
      if (r.failed):
        n_failed[str(ti)] += 1
  nfun_sums = []
  annotations = []
  for ti in spherical_types:
    print(ti, file=out)
    nfuns = nfun_accu[str(ti)]
    stats = nfuns.as_double().min_max_mean()
    stats.show(out=out, prefix="  ")
    nfun_sums.append((str(ti), flex.sum(nfuns)))
    if (n_failed[str(ti)] == 0):
      annotations.append(None)
    else:
      annotations.append("failed: %d" % n_failed[str(ti)])
  print(file=out)
  show_sorted_by_counts(
    label_count_pairs=nfun_sums,
    reverse=False,
    out=out,
    annotations=annotations)
  print(file=out)
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
