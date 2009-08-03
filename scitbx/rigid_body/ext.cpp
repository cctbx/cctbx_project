#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/boost_python/array_as_list.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/import.hpp>

#include <scitbx/rigid_body/tardy.h>

namespace scitbx { namespace rigid_body { namespace ext {

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
        .def("root_indices", &wt::root_indices)
        .def("pack_q", &wt::pack_q)
        .def("unpack_q", &wt::unpack_q, (arg_("q_packed")))
        .def("pack_qd", &wt::pack_qd)
        .def("unpack_qd", &wt::unpack_qd, (arg_("qd_packed")))
        .def("number_of_sites_in_each_tree",
          &wt::number_of_sites_in_each_tree)
        .def("sum_of_masses_in_each_tree", sum_of_masses_in_each_tree)
        .def("mean_linear_velocity", &wt::mean_linear_velocity, (
          arg_("number_of_sites_in_each_tree")))
        .def("subtract_from_linear_velocities",
          &wt::subtract_from_linear_velocities, (
            arg_("number_of_sites_in_each_tree"),
            arg_("value")))
        .def("e_kin", &wt::e_kin, ccr())
        .def("reset_e_kin", &wt::reset_e_kin, (
           arg_("e_kin_target"),
           arg_("e_kin_epsilon")=1e-12))
        .def("assign_zero_velocities", &wt::assign_zero_velocities)
        .def("assign_random_velocities", assign_random_velocities, (
           arg_("e_kin_target")=none,
           arg_("e_kin_epsilon")=1e-12,
           arg_("random_gauss")=none))
        .def("inverse_dynamics_packed", &wt::inverse_dynamics_packed, (
          arg_("qdd_packed")=none,
          arg_("f_ext_packed")=none,
          arg_("grav_accn")=none))
        .def("f_ext_as_tau_packed", &wt::f_ext_as_tau_packed, (
          arg_("f_ext_packed")))
        .def("forward_dynamics_ab_packed", &wt::forward_dynamics_ab_packed, (
          arg_("tau_packed")=none,
          arg_("f_ext_packed")=none,
          arg_("grav_accn")=none))
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
      object none;
      class_<wt,
             bases<featherstone::system_model<ft> >,
             boost::noncopyable
            >("tardy_model", no_init)
        .def(init<
          object const&,
          af::shared<vec3<ft> > const&,
          af::shared<ft> const&,
          object const&,
          object const&,
          optional<ft const&> >((
            arg_("labels"),
            arg_("sites"),
            arg_("masses"),
            arg_("tardy_tree"),
            arg_("potential_obj"),
            arg_("near_singular_hinges_angular_tolerance_deg")=5)))
        .def("flag_positions_as_changed", &wt::flag_positions_as_changed)
        .def("flag_velocities_as_changed", &wt::flag_velocities_as_changed)
        .def("sites_moved", &wt::sites_moved, ccr())
        .def("e_pot", &wt::e_pot, ccr())
        .def("d_e_pot_d_sites", &wt::d_e_pot_d_sites, ccr())
        .def("d_e_pot_d_q_packed", &wt::d_e_pot_d_q_packed)
        .def("e_tot", &wt::e_tot)
        .def("dynamics_step", &wt::dynamics_step, (arg_("delta_t")))
      ;
    }
  };

  void init_module()
  {
    featherstone_system_model_wrappers::wrap();
    tardy_model_wrappers::wrap();
  }

}}} // namespace scitbx::rigid_body::ext

BOOST_PYTHON_MODULE(scitbx_rigid_body_ext)
{
  scitbx::rigid_body::ext::init_module();
}
