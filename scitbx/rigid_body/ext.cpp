#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/boost_python/array_as_list.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/import.hpp>
#include <boost/python/tuple.hpp>

#include <scitbx/rigid_body/tardy.h>

namespace scitbx { namespace rigid_body { namespace ext {

  boost::python::tuple
  joint_lib_six_dof_aja_simplified_wrapper(
    vec3<double> const& center_of_mass,
    af::const_ref<double> const& q)
  {
    rotr3<double>
      result = joint_lib::six_dof_aja_simplified(center_of_mass, q);
    return boost::python::make_tuple(result.r, result.t);
  }

  struct featherstone_system_model_wrappers
  {
    typedef featherstone::system_model<> wt;
    typedef wt::ft ft;

    static
    boost::python::object
    sum_of_masses_in_each_tree(
      wt const& O)
    {
      af::shared<std::pair<int, double> >
        somiet = O.sum_of_masses_in_each_tree();
      return boost_python::array_as_list(somiet.begin(), somiet.size());
    }

    template <typename FloatType>
    struct random_gauss_adaptor_python
      : featherstone::random_gauss_adaptor<FloatType>
    {
      boost::python::object callable;

      random_gauss_adaptor_python(
        boost::python::object const& callable_)
      :
        callable(callable_)
      {
        namespace bp = boost::python;
        bp::object none;
        if (callable.ptr() == none.ptr()) {
          callable = bp::import("random").attr("gauss");
        }
      }

      virtual
      FloatType
      operator()(
        FloatType const& mu,
        FloatType const& sigma)
      {
        return boost::python::extract<FloatType>(callable(mu, sigma))();
      }
    };

    static
    boost::optional<af::shared<ft> >
    assign_random_velocities(
      wt& O,
      boost::optional<ft> const& e_kin_target,
      ft const& e_kin_epsilon,
      boost::python::object const& random_gauss)
    {
      random_gauss_adaptor_python<ft> rga(random_gauss);
      return O.assign_random_velocities(rga, e_kin_target, e_kin_epsilon);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      object none;
      class_<wt, boost::noncopyable>("featherstone_system_model", no_init)
        .def_readonly("number_of_trees", &wt::number_of_trees)
        .def_readonly("degrees_of_freedom", &wt::degrees_of_freedom)
        .def_readonly("q_packed_size", &wt::q_packed_size)
        .def("bodies_size", &wt::bodies_size)
        .def("degrees_of_freedom_each_joint",
          &wt::degrees_of_freedom_each_joint)
        .def("q_size_each_joint", &wt::q_size_each_joint)
        .def("root_indices", &wt::root_indices)
        .def("pack_q", &wt::pack_q)
        .def("unpack_q", &wt::unpack_q, (arg("q_packed")))
        .def("pack_qd", &wt::pack_qd)
        .def("unpack_qd", &wt::unpack_qd, (arg("qd_packed")))
        .def("number_of_sites_in_each_tree",
          &wt::number_of_sites_in_each_tree)
        .def("sum_of_masses_in_each_tree", sum_of_masses_in_each_tree)
        .def("mean_linear_velocity", &wt::mean_linear_velocity, (
          arg("number_of_sites_in_each_tree")))
        .def("subtract_from_linear_velocities",
          &wt::subtract_from_linear_velocities, (
            arg("number_of_sites_in_each_tree"),
            arg("value")))
        .def("e_kin", &wt::e_kin, ccr())
        .def("reset_e_kin", &wt::reset_e_kin, (
           arg("e_kin_target"),
           arg("e_kin_epsilon")=1e-12))
        .def("assign_zero_velocities", &wt::assign_zero_velocities)
        .def("assign_random_velocities", assign_random_velocities, (
           arg("e_kin_target")=none,
           arg("e_kin_epsilon")=1e-12,
           arg("random_gauss")=none))
        .def("inverse_dynamics_packed", &wt::inverse_dynamics_packed, (
          arg("qdd_packed")=none,
          arg("f_ext_packed")=none,
          arg("grav_accn")=none))
        .def("f_ext_as_tau_packed", &wt::f_ext_as_tau_packed, (
          arg("f_ext_packed")))
        .def("forward_dynamics_ab_packed", &wt::forward_dynamics_ab_packed, (
          arg("tau_packed")=none,
          arg("f_ext_packed")=none,
          arg("grav_accn")=none))
      ;
    }
  };

  struct tardy_model_wrappers
  {
    typedef tardy::model<> wt;
    typedef wt::ft ft;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_value_policy<return_by_value> rbv;
      object none;
      class_<wt,
             bases<featherstone::system_model<ft> >,
             boost::noncopyable
            >("tardy_model", no_init)
        .enable_pickling()
        .def(init<
          object const&,
          af::shared<vec3<ft> > const&,
          af::shared<ft> const&,
          object const&,
          object const&,
          optional<ft const&> >((
            arg("labels"),
            arg("sites"),
            arg("masses"),
            arg("tardy_tree"),
            arg("potential_obj"),
            arg("near_singular_hinges_angular_tolerance_deg")=5)))
        .def_readonly("labels", &wt::labels)
        .add_property("sites", make_getter(&wt::sites, rbv()))
        .add_property("masses", make_getter(&wt::masses, rbv()))
        .def_readonly("tardy_tree", &wt::tardy_tree)
        .def_readonly("potential_obj", &wt::potential_obj)
        .def_readonly("near_singular_hinges_angular_tolerance_deg",
          &wt::near_singular_hinges_angular_tolerance_deg)
        .def("flag_positions_as_changed", &wt::flag_positions_as_changed)
        .def("flag_velocities_as_changed", &wt::flag_velocities_as_changed)
        .def("sites_moved_is_cached", &wt::sites_moved_is_cached)
        .def("qdd_array_is_cached", &wt::qdd_array_is_cached)
        .def("sites_moved", &wt::sites_moved, ccr())
        .def("e_pot", &wt::e_pot, ccr())
        .def("d_e_pot_d_sites", &wt::d_e_pot_d_sites, ccr())
        .def("d_e_pot_d_q_packed", &wt::d_e_pot_d_q_packed)
        .def("e_tot", &wt::e_tot)
        .def("qdd_packed", &wt::qdd_packed)
        .def("dynamics_step", &wt::dynamics_step, (arg("delta_t")))
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    def("joint_lib_six_dof_aja_simplified",
      joint_lib_six_dof_aja_simplified_wrapper, (
        arg("center_of_mass"),
        arg("q")));

    featherstone_system_model_wrappers::wrap();
    tardy_model_wrappers::wrap();
  }

}}} // namespace scitbx::rigid_body::ext

BOOST_PYTHON_MODULE(scitbx_rigid_body_ext)
{
  scitbx::rigid_body::ext::init_module();
}
