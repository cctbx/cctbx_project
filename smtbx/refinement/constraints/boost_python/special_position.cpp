#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/special_position.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  struct special_position_site_wrapper
  {
    typedef special_position_site wt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt, bases<site_parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("special_position_site", no_init)
        .def(init<sgtbx::site_symmetry_ops const &, wt::scatterer_type *>
             ((arg("site_symmetry"), arg("scatterer"))))
        .add_property("independent_params",
                      make_function(&wt::independent_params, rir))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct special_position_u_star_parameter_wrapper
  {
    typedef special_position_u_star_parameter wt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt, bases<u_star_parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("special_position_u_star_parameter", no_init)
        .def(init<sgtbx::site_symmetry_ops const &,
                  wt::scatterer_type *>
             ((arg("site_symmetry"), arg("scatterer"))))
        .add_property("independent_params",
                      make_function(&wt::independent_params, rir))
      ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  void wrap_special_position() {
    special_position_site_wrapper::wrap();
    special_position_u_star_parameter_wrapper::wrap();
  }

}}}}
