#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/special_position.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  template <class wt, class ancestor>
  struct special_position_wrapper
  {
    static void wrap(char const *name) {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt,
             bases<ancestor>,
             std::auto_ptr<wt> >(name, no_init)
        .def(init<sgtbx::site_symmetry_ops const &,
                  typename wt::scatterer_type *>
             ((arg("site_symmetry"), arg("scatterer"))))
        .add_property("independent_params",
                      make_function(&wt::independent_params, rir))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  void wrap_special_position() {
    using namespace boost::python;
    special_position_wrapper<special_position_site_parameter,
                             asu_site_parameter>
    ::wrap("special_position_site_parameter");

    special_position_wrapper<special_position_u_star_parameter,
                             asu_u_star_parameter>
    ::wrap("special_position_u_star_parameter");
  }

}}}}
