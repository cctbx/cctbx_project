#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/scaled_adp.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct scalar_scaled_u_star_parameter_wrapper {
      typedef scalar_scaled_u_star_parameter wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_star_parameter>,
               boost::shared_ptr<wt> >("scalar_scaled_u_star", no_init)
          .def(init<independent_scalar_parameter *,
                    wt::scatterer_type *>
               ((arg("scalar"),
                 arg("scatterer"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct scalar_scaled_u_iso_parameter_wrapper {
      typedef scalar_scaled_u_iso_parameter wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_iso_parameter>,
               boost::shared_ptr<wt> >("scalar_scaled_u_iso", no_init)
          .def(init<independent_scalar_parameter *,
                    wt::scatterer_type *>
               ((arg("scalar"),
                 arg("scatterer"))))
          .add_property("reference", &wt::reference)
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    void wrap_scaled_adp() {
      scalar_scaled_u_star_parameter_wrapper::wrap();
      scalar_scaled_u_iso_parameter_wrapper::wrap();
    }


}}}}

