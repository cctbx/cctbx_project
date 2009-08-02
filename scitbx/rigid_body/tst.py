import scitbx.rigid_body.essence
from scitbx.rigid_body.essence import tst_tardy
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import \
  approx_equal, is_above_limit, is_below_limit, Exception_expected
from libtbx.utils import format_cpu_times
import random
import sys

def compare_essence_and_fast_tardy_models(etm, have_singularity=False):
  etm = scitbx.rigid_body.essence.tardy.model( # new instance to reset q, qd
    labels=etm.labels,
    sites=etm.sites,
    masses=etm.masses,
    tardy_tree=etm.tardy_tree,
    potential_obj=etm.potential_obj,
    near_singular_hinges_angular_tolerance_deg=
      etm.near_singular_hinges_angular_tolerance_deg)
  ftm = scitbx.rigid_body.tardy_model(
    labels=etm.labels,
    sites=flex.vec3_double(etm.sites),
    masses=flex.double(etm.masses),
    tardy_tree=etm.tardy_tree,
    potential_obj=etm.potential_obj,
    near_singular_hinges_angular_tolerance_deg=
      etm.near_singular_hinges_angular_tolerance_deg)
  assert ftm.bodies_size() == len(etm.bodies)
  assert ftm.number_of_trees == etm.number_of_trees
  assert ftm.degrees_of_freedom == etm.degrees_of_freedom
  assert ftm.q_packed_size == etm.q_packed_size
  ftm.flag_positions_as_changed()
  ftm.flag_velocities_as_changed()
  #
  assert list(ftm.root_indices()) == etm.root_indices()
  #
  q_packed_zero = etm.pack_q()
  qd_packed_zero = etm.pack_qd()
  #
  def check_packed():
    e = etm.pack_q()
    f = ftm.pack_q()
    assert approx_equal(e, f)
    e = etm.pack_qd()
    f = ftm.pack_qd()
    assert approx_equal(e, f)
  check_packed()
  mt = flex.mersenne_twister(seed=0)
  q_packed_rand = mt.random_double(size=q_packed_zero.size())*2-1
  qd_packed_rand = mt.random_double(size=qd_packed_zero.size())*2-1
  for tm in [etm, ftm]: tm.unpack_q(q_packed=q_packed_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_q(q_packed=q_packed_zero)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(qd_packed=qd_packed_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(qd_packed=qd_packed_zero)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_q(q_packed=q_packed_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(qd_packed=qd_packed_rand)
  check_packed()
  #
  fnosiet = ftm.number_of_sites_in_each_tree()
  assert list(fnosiet) == etm.number_of_sites_in_each_tree()
  e = etm.sum_of_masses_in_each_tree()
  f = ftm.sum_of_masses_in_each_tree()
  assert approx_equal(e, f)
  #
  def xxx_check_spatial_inertia():
    e = [body.i_spatial for body in etm.bodies]
    #print "e:", [sum(ei) for ei in e]
    f = ftm.xxx_spatial_inertia()
    #print "f:", [sum(fi) for fi in f]
    for ei,fi in zip(e, f):
      assert approx_equal(ei, fi)
  xxx_check_spatial_inertia()
  #
  def check_mean_linear_velocity():
    e = etm.mean_linear_velocity(number_of_sites_in_each_tree=None)
    for f in [ftm.mean_linear_velocity(number_of_sites_in_each_tree=None),
              ftm.mean_linear_velocity(number_of_sites_in_each_tree=fnosiet)]:
      if (e is None): assert f is None
      else:           assert approx_equal(e, f)
  check_mean_linear_velocity()
  value = matrix.col(mt.random_double(size=3)*2-1)
  for tm in [etm, ftm]:
    tm.subtract_from_linear_velocities(
      number_of_sites_in_each_tree=None,
      value=value)
  check_mean_linear_velocity()
  for tm in [etm, ftm]:
    tm.subtract_from_linear_velocities(
      number_of_sites_in_each_tree=fnosiet,
      value=value)
  check_mean_linear_velocity()
  #
  e = etm.sites_moved()
  f = ftm.sites_moved()
  assert approx_equal(e, f)
  e = etm.e_pot()
  f = ftm.e_pot()
  assert approx_equal(e, f)
  e = etm.d_e_pot_d_sites()
  f = ftm.d_e_pot_d_sites()
  assert approx_equal(e, f)
  e = etm.d_e_pot_d_q_packed()
  f = ftm.d_e_pot_d_q_packed()
  if (e.size() == 0):
    assert f.size() == 0
  else:
    e_max_abs = flex.max(flex.abs(e))
    if (etm.potential_obj is None):
      assert is_below_limit(value=e_max_abs, limit=1e-10)
    else:
      # note: e_max_abs is random generator dependent
      # by change much smaller e_max_abs are possible
      assert is_above_limit(value=e_max_abs, limit=0.1)
    assert approx_equal(e, f)
  #
  e = etm.e_kin()
  f = ftm.e_kin()
  assert approx_equal(e, f)
  e = etm.e_tot()
  f = ftm.e_tot()
  assert approx_equal(e, f)
  etm.reset_e_kin(e_kin_target=1.234)
  ftm.reset_e_kin(e_kin_target=1.234)
  e = etm.e_kin()
  if (etm.degrees_of_freedom == 0):
    assert approx_equal(e, 0)
  else:
    assert approx_equal(e, 1.234)
  f = ftm.e_kin()
  assert approx_equal(e, f)
  try:
    ftm.reset_e_kin(e_kin_target=1, e_kin_epsilon=0)
  except RuntimeError, e:
    assert str(e).find("e_kin_epsilon > 0") > 0
  else: raise Exception_expected
  etm.assign_zero_velocities()
  ftm.assign_zero_velocities()
  assert etm.e_kin() == 0
  assert ftm.e_kin() == 0
  random.seed(0)
  etm.assign_random_velocities()
  random.seed(0)
  ftm.assign_random_velocities()
  e = etm.featherstone_system_model().spatial_velocities()
  f = ftm.xxx_spatial_velocities()
  assert approx_equal(e, f)
  e = etm.e_kin()
  f = ftm.e_kin()
  xxx_check_spatial_inertia()
  assert approx_equal(e, f)
  #
  mt = flex.mersenne_twister(seed=0)
  qdd_rand_array = []
  for body in etm.bodies:
    qdd_rand_array.append(matrix.col(
      mt.random_double(size=body.joint.degrees_of_freedom)*2-1))
  tau_rand_array = []
  for body in etm.bodies:
    tau_rand_array.append(matrix.col(
      mt.random_double(size=body.joint.degrees_of_freedom)*2-1))
  f_ext_rand_array = []
  for ib in xrange(len(etm.bodies)):
    f_ext_rand_array.append(matrix.col(
      mt.random_double(size=6)*2-1))
  grav_accn_rand = matrix.col(mt.random_double(size=6)*2-1)
  def pack_array(array, packed_size=None):
    result = flex.double()
    if (packed_size is not None):
      result.reserve(packed_size)
    for sub in array:
      result.extend(flex.double(sub))
    if (packed_size is not None):
      assert result.size() == packed_size
    return result
  qdd_rand_packed = pack_array(
    array=qdd_rand_array, packed_size=etm.degrees_of_freedom)
  tau_rand_packed = pack_array(
    array=tau_rand_array, packed_size=etm.degrees_of_freedom)
  f_ext_rand_packed = pack_array(
    array=f_ext_rand_array, packed_size=len(etm.bodies)*6)
  def etm_inverse_dynamics_packed(
        qdd_array=None, f_ext_array=None, grav_accn=None):
    return pack_array(
      array=etm.featherstone_system_model().inverse_dynamics(
        qdd_array=qdd_array,
        f_ext_array=f_ext_array,
        grav_accn=grav_accn),
      packed_size=etm.degrees_of_freedom)
  def etm_forward_dynamics_ab_packed(
        tau_array=None, f_ext_array=None, grav_accn=None):
    return pack_array(
      array=etm.featherstone_system_model().forward_dynamics_ab(
        tau_array=tau_array,
        f_ext_array=f_ext_array,
        grav_accn=grav_accn),
      packed_size=etm.degrees_of_freedom)
  #
  e = etm_inverse_dynamics_packed(
    qdd_array=qdd_rand_array)
  f = ftm.inverse_dynamics_packed(
    qdd_packed=qdd_rand_packed)
  assert approx_equal(e, f)
  e = etm_inverse_dynamics_packed(
    f_ext_array=f_ext_rand_array)
  f = ftm.inverse_dynamics_packed(
    f_ext_packed=f_ext_rand_packed)
  assert approx_equal(e, f)
  e = etm_inverse_dynamics_packed(
    grav_accn=grav_accn_rand)
  f = ftm.inverse_dynamics_packed(
    grav_accn=flex.double(grav_accn_rand))
  assert approx_equal(e, f)
  for qdd_array,qdd_packed in [(None,None),
                               (qdd_rand_array,qdd_rand_packed)]:
    for f_ext_array,f_ext_packed in [(None,None),
                                     (f_ext_rand_array,f_ext_rand_packed)]:
      for grav_accn in [None, grav_accn_rand]:
        if ([qdd_array,f_ext_array,grav_accn].count(None) >= 2):
          continue # exercised above already
        e = etm_inverse_dynamics_packed(
          qdd_array=qdd_rand_array,
          f_ext_array=f_ext_rand_array,
          grav_accn=grav_accn_rand)
        f = ftm.inverse_dynamics_packed(
          qdd_packed=qdd_rand_packed,
          f_ext_packed=f_ext_rand_packed,
          grav_accn=flex.double(grav_accn_rand))
        assert approx_equal(e, f)
  #
  e = etm_forward_dynamics_ab_packed()
  f = ftm.forward_dynamics_ab_packed()
  assert approx_equal(e, f)
  e = etm_forward_dynamics_ab_packed(
    tau_array=tau_rand_array)
  f = ftm.forward_dynamics_ab_packed(
    tau_packed=tau_rand_packed)
  assert approx_equal(e, f)
  e = etm_forward_dynamics_ab_packed(
    f_ext_array=f_ext_rand_array)
  f = ftm.forward_dynamics_ab_packed(
    f_ext_packed=f_ext_rand_packed)
  assert approx_equal(e, f)
  e = etm_forward_dynamics_ab_packed(
    grav_accn=grav_accn_rand)
  f = ftm.forward_dynamics_ab_packed(
    grav_accn=flex.double(grav_accn_rand))
  assert approx_equal(e, f)
  for tau_array,tau_packed in [(None,None),
                               (tau_rand_array,tau_rand_packed)]:
    for f_ext_array,f_ext_packed in [(None,None),
                                     (f_ext_rand_array,f_ext_rand_packed)]:
      for grav_accn in [None, grav_accn_rand]:
        if (grav_accn is None):
          grav_accn_f = None
        else:
          grav_accn_f = flex.double(grav_accn)
        qdd_array_e = etm.featherstone_system_model().forward_dynamics_ab(
          tau_array=tau_array,
          f_ext_array=f_ext_array,
          grav_accn=grav_accn)
        qdd_packed_e = pack_array(
          array=qdd_array_e,
          packed_size=etm.degrees_of_freedom)
        qdd_packed_f = ftm.forward_dynamics_ab_packed(
          tau_packed=tau_packed,
          f_ext_packed=f_ext_packed,
          grav_accn=grav_accn_f)
        assert approx_equal(qdd_packed_e, qdd_packed_f)
        tau2_array_e = etm.featherstone_system_model().inverse_dynamics(
          qdd_array=qdd_array_e,
          f_ext_array=f_ext_array,
          grav_accn=grav_accn)
        tau2_packed_e = pack_array(
          array=tau2_array_e,
          packed_size=etm.degrees_of_freedom)
        tau2_packed_f = ftm.inverse_dynamics_packed(
          qdd_packed=qdd_packed_f,
          f_ext_packed=f_ext_packed,
          grav_accn=grav_accn_f)
        assert approx_equal(tau2_packed_e, tau2_packed_f)
        if (etm.degrees_of_freedom == 0):
          assert tau2_packed_e.size() == 0
        elif (not have_singularity):
          if (tau_packed is None):
            assert is_below_limit(
              flex.max(flex.abs(tau2_packed_e)), 0, eps=1e-5)
          else:
            assert approx_equal(tau2_packed_e, tau_packed, eps=1e-5)
        qdd2_array_e = etm.featherstone_system_model().forward_dynamics_ab(
          tau_array=tau2_array_e,
          f_ext_array=f_ext_array,
          grav_accn=grav_accn)
        qdd2_packed_e = pack_array(
          array=qdd2_array_e,
          packed_size=etm.degrees_of_freedom)
        qdd2_packed_f = ftm.forward_dynamics_ab_packed(
          tau_packed=tau2_packed_f,
          f_ext_packed=f_ext_packed,
          grav_accn=grav_accn_f)
        assert approx_equal(qdd2_packed_e, qdd2_packed_f)
        assert approx_equal(qdd2_packed_e, qdd_packed_e, eps=1e-4)
  #
  delta_t = 0.1234
  etm.dynamics_step(delta_t=delta_t)
  ftm.dynamics_step(delta_t=delta_t)
  e = etm.pack_q()
  f = etm.pack_q()
  assert approx_equal(e, f)
  e = etm.pack_qd()
  f = etm.pack_qd()
  assert approx_equal(e, f)
  #
  etm.assign_zero_velocities()
  ftm.assign_zero_velocities()
  e = etm_inverse_dynamics_packed(
    f_ext_array=f_ext_rand_array)
  f = ftm.inverse_dynamics_packed(
    f_ext_packed=f_ext_rand_packed)
  assert approx_equal(e, f)
  e_ref = e
  e = pack_array(
    array=etm.featherstone_system_model().f_ext_as_tau(
      f_ext_array=f_ext_array),
    packed_size=etm.degrees_of_freedom)
  assert approx_equal(e_ref, e)
  f = ftm.f_ext_as_tau_packed(f_ext_packed=f_ext_packed)
  assert approx_equal(e_ref, f)

def exercise_fixed_vertices():
  etm = tst_tardy.get_test_model_by_index(
    i=0, fixed_vertex_lists=[[0]])
  assert etm.degrees_of_freedom == 0
  compare_essence_and_fast_tardy_models(etm=etm)
  #
  for i_case,etm in enumerate(
                      tst_tardy.exercise_fixed_vertices_special_cases()):
    assert etm.potential_obj is not None
    compare_essence_and_fast_tardy_models(
      etm=etm,
      have_singularity=(i_case < 4))
  #
  for fixed_vertices,expected_dof in \
        tst_tardy.test_case_5_fixed_vertices_expected_dof:
    etm = tst_tardy.get_test_model_by_index(
      i=5, fixed_vertex_lists=[fixed_vertices])
    assert etm.degrees_of_freedom == expected_dof
    compare_essence_and_fast_tardy_models(etm=etm)

def run(args):
  assert len(args) == 0
  n_tested = 0
  for tc in tst_tardy_pdb.test_cases:
    tt = tc.tardy_tree_construct()
    masses = [1.0]*len(tc.sites)
    if (n_tested < 2):
      etm = scitbx.rigid_body.essence.tardy.model(
        labels=tc.labels,
        sites=tc.sites,
        masses=masses,
        tardy_tree=tt,
        potential_obj=None)
      compare_essence_and_fast_tardy_models(etm=etm)
    etm = tst_tardy.construct_tardy_model(
      labels=tc.labels,
      sites=tc.sites,
      masses=masses,
      tardy_tree=tt)
    assert etm.potential_obj is not None
    compare_essence_and_fast_tardy_models(etm=etm)
    n_tested += 1
  #
  exercise_fixed_vertices()
  #
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
