#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/symmetry_equivalent.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace boost_python {

  struct symmetry_equivalent_site_parameter_wrapper
  {
    typedef symmetry_equivalent_site_parameter wt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt, bases<site_parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("symmetry_equivalent_site_parameter", no_init)
        .def(init<site_parameter *, sgtbx::rt_mx const &>
             ((arg("site"), arg("motion"))))
        .add_property("original", make_function(&wt::original, rir))
        .add_property("motion", make_function(&wt::motion, rir))

        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  void wrap_symmetry_equivalent_site_parameter() {
    symmetry_equivalent_site_parameter_wrapper::wrap();
  }


}}}}
