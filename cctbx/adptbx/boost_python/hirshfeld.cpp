#include <cctbx/adptbx/hirshfeld.h>
#include <scitbx/boost_python/container_conversions.h>
#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/make_function.hpp>

namespace cctbx { namespace adptbx { namespace boost_python {

  template <typename T>
  struct mean_square_displacement_wrapper
  {
    typedef mean_square_displacement<T> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<wt>(name, no_init)
        .def(init<uctbx::unit_cell const &,
                  scitbx::vec3<T> const &>
             ((arg("unit_cell"), arg("z")))
             [with_custodian_and_ward<1, 2>()])
        .def("__call__", &wt::operator(), return_self<>())
        .add_property("value", &wt::value)
        .add_property("grad_u",
                      make_function(&wt::grad_u, rbv))
        .add_property("grad_z",
                      make_function(&wt::grad_z, rbv))
        .add_property("grad_g",
                      make_function(&wt::grad_g, rbv))
        .add_property("grad_unit_cell_params",
                      make_function(&wt::grad_unit_cell_params, rbv))
        .add_property("well_defined", &wt::well_defined)
        ;
    }
  };

  template <typename T>
  struct relative_hirshfeld_difference_wrapper
  {
    typedef relative_hirshfeld_difference<T> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      using namespace scitbx::boost_python::container_conversions;
      return_value_policy<return_by_value> rbv;
      tuple_mapping_fixed_size<af::tiny<T, 2*(3+6)> >();
      class_<wt>(name, no_init)
      .def(init<uctbx::unit_cell const &,
                scitbx::vec3<T> const &,
                scitbx::sym_mat3<T> const &,
                scitbx::vec3<T> const &,
                scitbx::sym_mat3<T> const &,
                sgtbx::rt_mx const &>
           ((arg("unit_cell"),
             arg("site_1"), arg("u_star_1"),
             arg("site_2"), arg("u_star_2"), arg("rt_mx_2"))))
      .add_property("value", &wt::value)
      .add_property("grad_sites_adps",
                    make_function(&wt::grad_sites_adps, rbv))
      .add_property("grad_unit_cell_params",
                    make_function(&wt::grad_unit_cell_params, rbv))
      ;
    }
  };


  void wrap_hirshfeld() {
    mean_square_displacement_wrapper<double>
    ::wrap("mean_square_displacement");
    relative_hirshfeld_difference_wrapper<double>
    ::wrap("relative_hirshfeld_difference");
  }

}}}
