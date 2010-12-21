#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/occupancy.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct dependent_occu_wrapper  {
      typedef dependent_occupancy wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_occupancy_parameter>,
               std::auto_ptr<wt> >("dependent_occupancy", no_init)
          .def(init<scalar_parameter *,
                    double,
                    bool,
                    wt::scatterer_type *>
               ((arg("occupancy"),
                 arg("multiplier"),
                 arg("as_one"),
                 arg("scatterer"))))
          .add_property("occupancy", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    void wrap_occupancy() {
      dependent_occu_wrapper::wrap();
    }


  }}}}
