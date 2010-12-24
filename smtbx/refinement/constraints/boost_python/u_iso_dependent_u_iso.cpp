#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/u_iso_dependent_u_iso.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct u_iso_proportional_to_pivot_u_iso_wrapper  {
      typedef u_iso_proportional_to_pivot_u_iso wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_iso_parameter>,
               std::auto_ptr<wt> >("u_iso_proportional_to_pivot_u_iso", no_init)
          .def(init<scalar_parameter *,
                    double,
                    wt::scatterer_type *>
               ((arg("pivot_u_iso"),
                 arg("multiplier"),
                 arg("scatterer"))))
          .add_property("pivot_u_iso", make_function(&wt::pivot_u_iso, rir))
          .def_readwrite("multiplier", &wt::multiplier)
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    void wrap_u_iso_dependent_u_iso() {
      u_iso_proportional_to_pivot_u_iso_wrapper::wrap();
    }


  }}}}
