#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/boost_python/array_as_list.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/rigid_body/tardy.h>

namespace scitbx { namespace rigid_body { namespace ext {

  struct tardy_model_wrappers
  {
    typedef tardy::model<> wt;
    typedef wt::ft ft;

    static
    boost::python::list
    xxx_spatial_inertia(
      wt const& O)
    {
      boost::python::list result;
      unsigned nb = O.bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        result.append(O.bodies[ib]->i_spatial);
      }
      return result;
    }

    static
    boost::python::list
    xxx_spatial_velocities(
      wt& O)
    {
      boost::python::list result;
      af::shared<af::tiny<ft, 6> >
        sv = O.featherstone_system_model().spatial_velocities();
      SCITBX_ASSERT(sv.size() == O.bodies_size());
      unsigned nb = O.bodies_size();
      for(unsigned ib=0;ib<nb;ib++) {
        result.append(boost_python::array_as_list(
          sv[ib].begin(), sv[ib].size()));
      }
      return result;
    }

    static
    boost::python::object
    sum_of_masses_in_each_tree(
      wt const& O)
    {
      af::shared<std::pair<int, double> >
        somiet = O.sum_of_masses_in_each_tree();
      return boost_python::array_as_list(somiet.begin(), somiet.size());
    }

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      reset_e_kin_overloads, reset_e_kin, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      assign_random_velocities_overloads, assign_random_velocities, 0, 3)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<wt, boost::noncopyable>("tardy_model", no_init)
        .def("bodies_size", &wt::bodies_size)
        .def_readonly("number_of_trees", &wt::number_of_trees)
        .def_readonly("degrees_of_freedom", &wt::degrees_of_freedom)
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
        .def("root_indices", &wt::root_indices)
        .def("number_of_sites_in_each_tree",
          &wt::number_of_sites_in_each_tree)
        .def("sum_of_masses_in_each_tree", sum_of_masses_in_each_tree)
        .def("mean_linear_velocity", &wt::mean_linear_velocity, (
          arg_("number_of_sites_in_each_tree")))
        .def("subtract_from_linear_velocities",
          &wt::subtract_from_linear_velocities, (
            arg_("number_of_sites_in_each_tree"),
            arg_("value")))
        .def("sites_moved", &wt::sites_moved, ccr())
        .def("e_pot", &wt::e_pot, ccr())
        .def("d_e_pot_d_sites", &wt::d_e_pot_d_sites, ccr())
        .def("e_kin", &wt::e_kin, ccr())
        .def("e_tot", &wt::e_tot)
        .def("reset_e_kin", &wt::reset_e_kin, reset_e_kin_overloads((
           arg_("e_kin_target"),
           arg_("e_kin_epsilon")=1e-12)))
        .def("assign_zero_velocities", &wt::assign_zero_velocities)
        .def("assign_random_velocities", &wt::assign_random_velocities,
          assign_random_velocities_overloads((
             arg_("e_kin_target")=object(),
             arg_("e_kin_epsilon")=1e-12,
             arg_("random_gauss")=object())))
        .def("dynamics_step", &wt::dynamics_step, (arg_("delta_t")))
        .def("pack_q", &wt::pack_q)
        .def("unpack_q", &wt::unpack_q, (arg_("packed_q")))
        .def("pack_qd", &wt::pack_qd)
        .def("unpack_qd", &wt::unpack_qd, (arg_("packed_qd")))
        .def("xxx_spatial_inertia", xxx_spatial_inertia)
        .def("xxx_spatial_velocities", xxx_spatial_velocities)
      ;
    };
  };

  void init_module()
  {
    tardy_model_wrappers::wrap();
  }

}}} // namespace scitbx::rigid_body::ext

BOOST_PYTHON_MODULE(scitbx_rigid_body_ext)
{
  scitbx::rigid_body::ext::init_module();
}
