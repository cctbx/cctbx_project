#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <scitbx/rigid_body/featherstone.h>

namespace scitbx { namespace rigid_body { namespace ext {

  void init_module()
  {
    using namespace boost::python;
    joint_lib::zero_dof_alignment<>();
    joint_lib::zero_dof<> zero_dof;
    joint_lib::six_dof_alignment<> six_dof_alignment(vec3<double>(0,0,0));
    six_dof_alignment.cb_0b * six_dof_alignment.cb_b0;
    joint_lib::six_dof<> six_dof(
      af::tiny<double, 4>(1,0,0,0), vec3<double>(0,0,0));
    six_dof.qd_zero();
    six_dof.qdd_zero();
    six_dof.motion_subspace();
    six_dof.get_linear_velocity(af::const_ref<double>(0, 6));
    six_dof.new_linear_velocity(
      af::const_ref<double>(0, 6),
      af::const_ref<double>(0, 3));
    six_dof.time_step_position(af::const_ref<double>(0, 6), 1e-6);
    six_dof.time_step_velocity(
      af::const_ref<double>(0, 6),
      af::const_ref<double>(0, 6), 1e-6);
    six_dof.tau_as_d_pot_d_q(af::small<double, 6>(6));
    six_dof.get_q();
    six_dof.new_q(af::const_ref<double>(0, 7));
    joint_lib::spherical_alignment<> spherical_alignment(vec3<double>(0,0,0));
    spherical_alignment.cb_0b * spherical_alignment.cb_b0;
    joint_lib::spherical<> spherical(af::tiny<double, 4>(1,0,0,0));
    spherical.qd_zero();
    spherical.qdd_zero();
    spherical.motion_subspace();
    spherical.get_linear_velocity(af::const_ref<double>(0, 3));
    spherical.new_linear_velocity(
      af::const_ref<double>(0, 3),
      af::const_ref<double>(0, 3));
    spherical.time_step_position(af::const_ref<double>(0, 3), 1e-6);
    spherical.time_step_velocity(
      af::const_ref<double>(0, 3),
      af::const_ref<double>(0, 3), 1e-6);
    spherical.tau_as_d_pot_d_q(af::small<double, 6>(3));
    spherical.get_q();
    spherical.new_q(af::const_ref<double>(0, 4));
    joint_lib::revolute_alignment<> revolute_alignment(
      vec3<double>(0,0,0), vec3<double>(1,0,0));
    revolute_alignment.cb_0b * revolute_alignment.cb_b0;
    joint_lib::revolute<> revolute(af::tiny<double, 1>(0));
    revolute.qd_zero();
    revolute.qdd_zero();
    revolute.motion_subspace();
    revolute.get_linear_velocity(af::const_ref<double>(0, 1));
    revolute.new_linear_velocity(
      af::const_ref<double>(0, 1),
      af::const_ref<double>(0, 1));
    revolute.time_step_position(af::const_ref<double>(0, 1), 1e-6);
    revolute.time_step_velocity(
      af::const_ref<double>(0, 1),
      af::const_ref<double>(0, 1), 1e-6);
    revolute.tau_as_d_pot_d_q(af::small<double, 6>(1));
    revolute.get_q();
    revolute.new_q(af::const_ref<double>(0, 1));
    joint_lib::translational_alignment<> translational_alignment(
      vec3<double>(0,0,0));
    translational_alignment.cb_0b * translational_alignment.cb_b0;
    joint_lib::translational<> translational(vec3<double>(0,0,0));
    translational.qd_zero();
    translational.qdd_zero();
    translational.motion_subspace();
    translational.get_linear_velocity(af::const_ref<double>(0, 3));
    translational.new_linear_velocity(
      af::const_ref<double>(0, 3),
      af::const_ref<double>(0, 3));
    translational.time_step_position(af::const_ref<double>(0, 3), 1e-6);
    translational.time_step_velocity(
      af::const_ref<double>(0, 3),
      af::const_ref<double>(0, 3), 1e-6);
    translational.tau_as_d_pot_d_q(af::small<double, 6>(3));
    translational.get_q();
    translational.new_q(af::const_ref<double>(0, 3));
    af::shared<boost::shared_ptr<featherstone::body_t<double> > > bodies;
    featherstone::xrot(mat3<double>());
    featherstone::xtrans(vec3<double>());
    featherstone::cb_as_spatial_transform(rotr3<double>());
    featherstone::system_model<> system_model(bodies);
    system_model.bodies_size();
    system_model.cb_up_array();
    system_model.xup_array();
    system_model.spatial_velocities();
    system_model.e_kin();
    system_model.accumulated_spatial_inertia();
    system_model.qd_e_kin_scales();
    system_model.inverse_dynamics(
      af::const_ref<af::small<double, 6> >(0, 0),
      af::const_ref<af::tiny<double, 6> >(0, 0),
      af::const_ref<double>(0,0));
    system_model.f_ext_as_tau(
      af::const_ref<af::tiny<double, 6> >(0, 0));
    system_model.d_pot_d_q(
      af::const_ref<af::tiny<double, 6> >(0, 0));
    system_model.forward_dynamics_ab(
      af::const_ref<af::small<double, 6> >(0, 0),
      af::const_ref<af::tiny<double, 6> >(0, 0),
      af::const_ref<double>(0, 0));
  }

}}} // namespace scitbx::rigid_body::ext

BOOST_PYTHON_MODULE(scitbx_rigid_body_ext)
{
  scitbx::rigid_body::ext::init_module();
}
