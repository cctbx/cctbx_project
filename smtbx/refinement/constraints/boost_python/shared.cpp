#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/shared.h>
#include <smtbx/refinement/constraints/scaled_adp.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct shared_u_star_wrapper  {
      typedef shared_u_star wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_star_parameter>,
               boost::shared_ptr<wt> >("shared_u_star", no_init)
          .def(init<wt::scatterer_type *,
                    u_star_parameter *>
               ((arg("scatterer"),
                 arg("reference"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct shared_rotated_u_star_wrapper {
      typedef shared_rotated_u_star wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_star_parameter>,
               boost::shared_ptr<wt> >("shared_rotated_u_star", no_init)
          .def(init<wt::scatterer_type *,
                    u_star_parameter *,
                    direction_base *,
                    independent_scalar_parameter *>
               ((arg("scatterer"),
                 arg("reference"),
                 arg("direction"),
                 arg("angle"))))
          .add_property("reference", make_function(&wt::reference, rir))
          .add_property("angle", make_function(&wt::angle, rir))
          .add_property("direction", make_function(&wt::direction, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct shared_rotating_u_star_wrapper {
      typedef shared_rotating_u_star wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
          bases<asu_u_star_parameter>,
          boost::shared_ptr<wt> >("shared_rotating_u_star", no_init)
          .def(init<wt::scatterer_type*,
            u_star_parameter*,
            independent_scalar_parameter*,
            independent_scalar_parameter*,
            independent_scalar_parameter*,
            independent_scalar_parameter*>
            ((arg("scatterer"),
              arg("reference"),
              arg("scale"),
              arg("alpha"),
              arg("beta"),
              arg("gamma"))))
          .add_property("reference", make_function(&wt::reference, rir))
          .add_property("scale", make_function(&wt::scale, rir))
          .add_property("alpha", make_function(&wt::alpha, rir))
          .add_property("beta", make_function(&wt::beta, rir))
          .add_property("gamma", make_function(&wt::gamma, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };
    struct shared_u_iso_wrapper  {
      typedef shared_u_iso wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_u_iso_parameter>,
               boost::shared_ptr<wt> >("shared_u_iso", no_init)
          .def(init<wt::scatterer_type *,
                    scalar_parameter *>
               ((arg("scatterer"),
                 arg("reference"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct shared_site_wrapper  {
      typedef shared_site wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_site_parameter>,
               boost::shared_ptr<wt> >("shared_site", no_init)
          .def(init<wt::scatterer_type *,
                    site_parameter *>
               ((arg("scatterer"),
                 arg("reference"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct shared_fp_wrapper {
      typedef shared_fp wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
          bases<asu_fp_parameter>,
          boost::shared_ptr<wt> >("shared_fp", no_init)
          .def(init<wt::scatterer_type *,
            scalar_parameter *>
            ((arg("scatterer"),
              arg("reference"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    struct shared_fdp_wrapper {
      typedef shared_fdp wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
          bases<asu_fdp_parameter>,
          boost::shared_ptr<wt> >("shared_fdp", no_init)
          .def(init<wt::scatterer_type *,
            scalar_parameter *>
            ((arg("scatterer"),
              arg("reference"))))
          .add_property("reference", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<boost::shared_ptr<wt>, boost::shared_ptr<parameter> >();
      }
    };

    void wrap_shared() {
      shared_u_star_wrapper::wrap();
      shared_rotated_u_star_wrapper::wrap();
      shared_rotating_u_star_wrapper::wrap();
      shared_u_iso_wrapper::wrap();
      shared_site_wrapper::wrap();
      shared_fp_wrapper::wrap();
      shared_fdp_wrapper::wrap();
    }

}}}}
