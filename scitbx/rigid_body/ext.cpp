#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <scitbx/rigid_body/joint_lib.h>

namespace scitbx { namespace rigid_body { namespace ext {

  void init_module()
  {
    using namespace boost::python;
    joint_lib::zero_dof_alignment<> zero_dof_alignment;
    joint_lib::zero_dof<> zero_dof;
    joint_lib::six_dof_alignment<> six_dof_alignment(vec3<double>(0,0,0));
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
    six_dof.tau_as_d_pot_d_q(af::const_ref<double>(0, 6));
    six_dof.get_q();
    six_dof.new_q(af::const_ref<double>(0, 7));
  }

}}} // namespace scitbx::lbfgs::ext

BOOST_PYTHON_MODULE(scitbx_rigid_body_ext)
{
  scitbx::rigid_body::ext::init_module();
}
