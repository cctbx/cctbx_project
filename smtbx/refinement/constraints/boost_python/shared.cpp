#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/shared.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct shared_u_star_wrapper  {
      typedef shared_u_star wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_star_parameter>,
               std::auto_ptr<wt> >("shared_u_star", no_init)
          .def(init<u_star_parameter *,
                    wt::scatterer_type *>
               ((arg("u_star"),
                 arg("scatterer"))))
          .add_property("u_star", make_function(&wt::original, rir))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct shared_u_iso_wrapper  {
      typedef shared_u_iso wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_iso_parameter>,
               std::auto_ptr<wt> >("shared_u_iso", no_init)
          .def(init<scalar_parameter *,
                    wt::scatterer_type *>
               ((arg("u_iso"),
                 arg("scatterer"))))
          .add_property("u_iso", make_function(&wt::original, rir))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct shared_site_wrapper  {
      typedef shared_site wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_site_parameter>,
               std::auto_ptr<wt> >("shared_site", no_init)
          .def(init<site_parameter *,
                    wt::scatterer_type *>
               ((arg("site"),
                 arg("scatterer"))))
          .add_property("site", make_function(&wt::original, rir))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    void wrap_shared() {
      shared_u_star_wrapper::wrap();
      shared_u_iso_wrapper::wrap();
      shared_site_wrapper::wrap();
    }


  }}}}
