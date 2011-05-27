#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <smtbx/refinement/constraints/direction.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  struct direction_base_wrapper {
    typedef direction_base wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("direction_base", no_init)
        .def("direction", &wt::direction)
        ;
    }
  };

  struct static_direction_wrapper {
    typedef static_direction wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<direction_base>,
          std::auto_ptr<wt> >("static_direction", no_init)
        .def(init<cart_t const &>
          ((arg("direction"))))
        .def("best_line", &wt::best_line,
          (arg("points")))
          .staticmethod("best_line")
        .def("calc_best_line", &wt::calc_best_line,
          (arg("unit_cell"), arg("sites")))
          .staticmethod("calc_best_line")
        .def("best_plane_normal", &wt::best_plane_normal,
          (arg("points")))
          .staticmethod("best_plane_normal")
        .def("calc_best_plane_normal", &wt::calc_best_plane_normal,
          (arg("unit_cell"), arg("sites")))
          .staticmethod("calc_best_plane_normal")
        ;
    }
  };

  struct vector_direction_wrapper {
    typedef vector_direction wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<direction_base>,
          std::auto_ptr<wt> >("vector_direction", no_init)
        .def(init<site_parameter *,
          site_parameter *>
            ((arg("from"),
              arg("to"))))
        .def(init<af::shared<site_parameter *> const &>
          ((arg("sites"))))
        ;
    }
  };

  struct normal_direction_wrapper {
    typedef normal_direction wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<direction_base>,
          std::auto_ptr<wt> >("normal_direction", no_init)
        .def(init<af::shared<site_parameter *> const &>
          ((arg("sites"))))
        ;
    }
  };

  void wrap_direction() {
    using namespace boost::python;
    direction_base_wrapper::wrap();
    static_direction_wrapper::wrap();
    vector_direction_wrapper::wrap();
    normal_direction_wrapper::wrap();
  }

}}}}
