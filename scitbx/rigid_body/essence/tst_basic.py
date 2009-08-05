import featherstone
import body_lib
import joint_lib
import spatial_lib

matrix = featherstone.matrix
if (featherstone.scitbx is not None):
  from libtbx.test_utils import approx_equal
else:
  def approx_equal(a1, a2): return True
  print "libtbx.test_utils not available: approx_equal() disabled"
  def sum(l):
    result = 0
    for e in l: result += e
    return result

import sys

def exercise_basic():
  assert approx_equal(sum(spatial_lib.xrot((1,2,3,4,5,6,7,8,9))), 90)
  assert approx_equal(sum(spatial_lib.xtrans((1,2,3))), 6)
  assert approx_equal(
    sum(spatial_lib.cb_as_spatial_transform(
      cb=matrix.rt(((1,2,3,4,5,6,7,8,9), (1,2,3))))), 90)
  assert approx_equal(sum(spatial_lib.crm((1,2,3,4,5,6))), 0)
  assert approx_equal(sum(spatial_lib.crf((1,2,3,4,5,6))), 0)
  i_spatial = spatial_lib.mci(
    m=1.234,
    c=matrix.col((1,2,3)),
    i=matrix.sym(sym_mat3=(2,3,4,0.1,0.2,0.3)))
  assert approx_equal(sum(i_spatial), 21.306)
  assert approx_equal(spatial_lib.kinetic_energy(
    i_spatial=i_spatial, v_spatial=matrix.col((1,2,3,4,5,6))), 75.109)
  #
  mass_points = body_lib.mass_points(
    sites=matrix.col_list([
      (0.949, 2.815, 5.189),
      (0.405, 3.954, 5.917),
      (0.779, 5.262, 5.227)]),
    masses=[2.34, 3.56, 1.58])
  assert approx_equal(mass_points.sum_of_masses(), 7.48)
  assert approx_equal(mass_points.center_of_mass(),
    [0.654181818182, 3.87397058824, 5.54350802139])
  assert approx_equal(mass_points._sum_of_masses, 7.48)
  assert approx_equal(
    mass_points.inertia(pivot=matrix.col((0.9,-1.3,0.4))),
    [404.7677928, 10.04129606, 10.09577652,
     10.04129606, 199.7384559, -199.3511949,
     10.09577652, -199.3511949, 206.8314171])

class zero_dof_body(object):

  def __init__(O):
    O.alignment = joint_lib.zero_dof_alignment()
    O.i_spatial = matrix.sqr([0]*36)
    O.joint = joint_lib.zero_dof()
    O.qd = O.joint.qd_zero
    assert O.joint.get_linear_velocity(qd=O.qd) is None
    assert O.joint.new_linear_velocity(qd=None, value=None) is None
    O.parent = -1

class six_dof_body(object):

  def __init__(O):
    sites = matrix.col_list([
      (0.949, 2.815, 5.189),
      (0.405, 3.954, 5.917),
      (0.779, 5.262, 5.227)])
    mass_points = body_lib.mass_points(sites=sites, masses=[1.0, 1.0, 1.0])
    O.alignment = joint_lib.six_dof_alignment(
      center_of_mass=mass_points.center_of_mass())
    O.i_spatial = mass_points.spatial_inertia(
      alignment_cb_0b=O.alignment.cb_0b)
    qe = matrix.col((0.18, 0.36, 0.54, -0.73)).normalize()
    qr = matrix.col((-0.1,0.3,0.2))
    O.joint = joint_lib.six_dof(qe=qe, qr=qr)
    O.qd = matrix.col((0.18,-0.02,-0.16,0.05,0.19,-0.29))
    assert O.joint.get_linear_velocity(qd=O.qd).elems == (0.05,0.19,-0.29)
    O.qd = O.joint.new_linear_velocity(
      qd=O.qd, value=matrix.col((-0.05,-0.19,0.29)))
    assert O.joint.get_linear_velocity(qd=O.qd).elems == (-0.05,-0.19,0.29)
    O.parent = -1

class spherical_body(object):

  def __init__(O):
    sites = matrix.col_list([
      (0.04, -0.16, 0.19),
      (0.10, -0.15, 0.18)])
    mass_points = body_lib.mass_points(sites=sites, masses=[1.0, 1.0])
    O.alignment = joint_lib.spherical_alignment(
      pivot=mass_points.center_of_mass())
    O.i_spatial = mass_points.spatial_inertia(
      alignment_cb_0b=O.alignment.cb_0b)
    qe = matrix.col((-0.50, -0.33, 0.67, -0.42)).normalize()
    O.joint = joint_lib.spherical(qe=qe)
    O.qd = matrix.col((0.12, -0.08, 0.11))
    assert O.joint.get_linear_velocity(qd=O.qd) is None
    assert O.joint.new_linear_velocity(qd=None, value=None) is None
    O.parent = 0

class revolute_body(object):

  def __init__(O, parent):
    pivot = matrix.col((0.779, 5.262, 5.227))
    normal = matrix.col((0.25, 0.86, -0.45)).normalize()
    sites = matrix.col_list([(-0.084, 6.09, 4.936)])
    mass_points = body_lib.mass_points(sites=sites, masses=[1.0])
    O.alignment = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.i_spatial = mass_points.spatial_inertia(
      alignment_cb_0b=O.alignment.cb_0b)
    O.joint = joint_lib.revolute(qe=matrix.col([0.26]))
    O.qd = matrix.col([-0.19])
    assert O.joint.get_linear_velocity(qd=O.qd) is None
    assert O.joint.new_linear_velocity(qd=None, value=None) is None
    O.parent = parent

class translational_body(object):

  def __init__(O):
    sites = [matrix.col((0.949, 2.815, 5.189))]
    mass_points = body_lib.mass_points(sites=sites, masses=[1.0])
    O.alignment = joint_lib.translational_alignment(
      center_of_mass=mass_points.center_of_mass())
    O.i_spatial = mass_points.spatial_inertia(
      alignment_cb_0b=O.alignment.cb_0b)
    qr = matrix.col((-0.1,0.3,0.2))
    O.joint = joint_lib.translational(qr=qr)
    O.qd = matrix.col((0.05,0.19,-0.29))
    assert O.joint.get_linear_velocity(qd=O.qd).elems == O.qd.elems
    O.qd = O.joint.new_linear_velocity(qd=O.qd, value=-O.qd)
    assert O.joint.get_linear_velocity(qd=O.qd).elems == O.qd.elems
    O.parent = -1

def exercise_system_model():
  bodies = [
    six_dof_body(),
    spherical_body(),
    revolute_body(parent=1),
    translational_body()]
  body_lib.set_cb_tree(bodies=bodies)
  model = featherstone.system_model(bodies=bodies)
  assert approx_equal(model.e_kin(), 5.10688665235)
  assert approx_equal(model.qd_e_kin_scales(), [
    0.1036643, 0.1054236, 0.1187526, 0.5773503, 0.5773503, 0.5773503,
    0.1749883, 0.2830828, 0.2225619,
    1.334309,
    1.414214, 1.414214, 1.414214])
  #
  qdd = matrix.col_list([
    (-0.04,0.05,0.23,-0.01,-0.08,0.04),
    (0.08,-0.08,-0.01),
    (0.14,),
    (-0.01, -0.34, 0.28)])
  f_ext = matrix.col_list([
    (-0.10, 0.30, -0.01, -0.01, 0.01, 0.06),
    (0.28, 0.09, 0.14, 0.23, 0.00, 0.07),
    (-0.11, 0.03, -0.07, -0.11, 0.06, 0.08),
    (-0.16, 0.14, -0.33, 0.35, -0.02, -0.20)])
  grav_accn = matrix.col((0.02, -0.13, 0.15, 0.26, -0.16, 0.14))
  #
  qdd0 = [
    (0,0,0,0,0,0),
    (0,0,0),
    (0,),
    (0,0,0)]
  for qdd_array in [None, qdd0]:
    tau = model.inverse_dynamics(qdd_array=qdd_array)
    assert approx_equal(tau, [
      (3.0067673409496019, -1.9747070164103167, -0.96510418705493095,
       0.62119145987365987, 0.79528549692226591, 0.50582706679253908),
      (-3.5089966946778439, 0.85280077414986188, -1.1466929846982585),
      (-0.68035016279655447,),
      (0.0, 0.0, 0.0)])
  qdd2 = model.forward_dynamics_ab(tau_array=tau)
  assert approx_equal(qdd2, qdd0)
  #
  tau = model.inverse_dynamics(qdd_array=qdd)
  assert approx_equal(tau, [
    (-28.4935967396, -13.9449610757, 37.119813341,
     3.09036984758, -3.29209848977, 1.51871803584),
    (22.0723956539, 1.85204188959, -1.96550741514),
    (0.526685445884,),
    (-0.01,-0.34,0.28)])
  qdd2 = model.forward_dynamics_ab(tau_array=tau)
  assert approx_equal(qdd2, qdd)
  #
  tau = model.inverse_dynamics(qdd_array=qdd, f_ext_array=f_ext)
  assert approx_equal(tau, [
    (-28.2069898504, -14.1650076325, 38.2656278316,
     3.24402886492, -3.30833610224, 1.36344785723),
    (22.6135703212, 2.25806832722, -2.77922603881),
    (0.596685445884,),
    (-0.36,-0.32,0.48)])
  qdd2 = model.forward_dynamics_ab(tau_array=tau, f_ext_array=f_ext)
  assert approx_equal(qdd2, qdd)
  #
  tau = model.inverse_dynamics(
    qdd_array=qdd, f_ext_array=f_ext, grav_accn=grav_accn)
  assert approx_equal(tau, [
    (29.0716177639, 1.548329665, -9.90799285557,
     -2.51132634591, 5.78686348626, -3.77518591503),
    (-4.70390288163, -1.2432778965, 1.1909533225),
    (-0.354638347024,),
    (0.54782, -0.17957, 0.16733)])
  qdd2 = model.forward_dynamics_ab(
    tau_array=tau, f_ext_array=f_ext, grav_accn=grav_accn)
  assert approx_equal(qdd2, qdd)
  #
  new_q = [
    body.joint.time_step_position(qd=body.qd, delta_t=0.01).get_q()
      for body in model.bodies]
  assert approx_equal(new_q, [
    (0.18036749, 0.36210928, 0.54329229, -0.7356480,
     -0.10189282, 0.29788946, 0.2020574),
    (-0.50329486, -0.33273860, 0.67548709, -0.4239072),
    (0.2581,),
    (-0.1005, 0.2981, 0.2029)])
  new_qd = [
    body.joint.time_step_velocity(qd=body.qd, qdd=qdd_i, delta_t=0.01)
      for body,qdd_i in zip(model.bodies, qdd)]
  assert approx_equal(new_qd, [
    (0.1796, -0.0195, -0.1577, -0.0501, -0.1908, 0.2904),
    (0.1208, -0.0808, 0.1099),
    (-0.1886,),
    (-0.0501, -0.1934, 0.2928)])
  for body,q in zip(model.bodies, [(1,2,3,4,5,6,7),(8,9,10,11),(12,)]):
    assert approx_equal(body.joint.new_q(q=q).get_q(), q)
  #
  qdd = []
  for body in model.bodies:
    body.qd = body.joint.qd_zero
    qdd.append(body.joint.qdd_zero)
  model.flag_velocities_as_changed()
  tau = model.inverse_dynamics(qdd_array=qdd, f_ext_array=f_ext)
  assert approx_equal(tau, [
    (0.286606889188, -0.220046556736, 1.14581449056,
     0.153659017347, -0.0162376124727, -0.155270178613),
    (0.541174667239, 0.406026437633, -0.813718623664),
    (0.07,),
    (-0.35,0.02,0.2)])
  tau0 = model.f_ext_as_tau(f_ext_array=f_ext)
  assert approx_equal(tau0, tau)

def exercise_system_model_with_zero_dof_body():
  bodies = [
    zero_dof_body(),
    revolute_body(parent=0)]
  body_lib.set_cb_tree(bodies=bodies)
  model = featherstone.system_model(bodies=bodies)
  assert approx_equal(model.e_kin(), 0.0202765671829)
  assert approx_equal(model.qd_e_kin_scales(), [1.334309])
  #
  qdd = matrix.col_list([
    (),
    (0.14,)])
  f_ext = matrix.col_list([
    (-0.10, 0.30, -0.01, -0.01, 0.01, 0.06),
    (-0.11, 0.03, -0.07, -0.11, 0.06, 0.08)])
  grav_accn = matrix.col((0.02, -0.13, 0.15, 0.26, -0.16, 0.14))
  #
  tau = model.inverse_dynamics(qdd_array=qdd)
  assert approx_equal(tau, [
    (),
    (0.15726977316344815,)])
  qdd2 = model.forward_dynamics_ab(tau_array=tau)
  assert approx_equal(qdd2, qdd)
  #
  tau = model.inverse_dynamics(qdd_array=qdd, f_ext_array=f_ext)
  assert approx_equal(tau, [
    (),
    (0.22726977316344815,)])
  qdd2 = model.forward_dynamics_ab(tau_array=tau, f_ext_array=f_ext)
  assert approx_equal(qdd2, qdd)
  #
  tau = model.inverse_dynamics(
    qdd_array=qdd, f_ext_array=f_ext, grav_accn=grav_accn)
  assert approx_equal(tau, [
    (),
    (0.59601742875022201,)])
  qdd2 = model.forward_dynamics_ab(
    tau_array=tau, f_ext_array=f_ext, grav_accn=grav_accn)
  assert approx_equal(qdd2, qdd)
  #
  new_q = [
    body.joint.time_step_position(qd=body.qd, delta_t=0.01).get_q()
      for body in model.bodies]
  assert approx_equal(new_q, [
    (), (0.2581,)])
  new_qd = [
    body.joint.time_step_velocity(qd=body.qd, qdd=qdd_i, delta_t=0.01)
      for body,qdd_i in zip(model.bodies, qdd)]
  assert approx_equal(new_qd, [
    (),
    (-0.1886,)])
  for body,q in zip(model.bodies, [(),(13,)]):
    assert approx_equal(body.joint.new_q(q=q).get_q(), q)
  #
  qdd = []
  for body in model.bodies:
    body.qd = body.joint.qd_zero
    qdd.append(body.joint.qdd_zero)
  model.flag_velocities_as_changed()
  tau = model.inverse_dynamics(qdd_array=qdd, f_ext_array=f_ext)
  assert approx_equal(tau, [
    (),
    (0.07,)])
  tau0 = model.f_ext_as_tau(f_ext_array=f_ext)
  assert approx_equal(tau0, tau)

def run(args):
  assert len(args) == 0
  exercise_basic()
  exercise_system_model()
  exercise_system_model_with_zero_dof_body()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
