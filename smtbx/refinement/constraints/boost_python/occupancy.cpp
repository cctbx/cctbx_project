#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/refinement/constraints/occupancy.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct affine_asu_occupancy_parameter_wrapper {
      typedef affine_asu_occupancy_parameter wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<affine_scalar_parameter, asu_parameter>,
               std::auto_ptr<wt> >("affine_asu_occupancy_parameter", no_init)
          .def(init<scalar_parameter *, double, double,
                    wt::scatterer_type *>
               ((arg("dependee"), arg("a"), arg("b"), arg("scatterer"))))
          .def(init<scalar_parameter *, double,
                    scalar_parameter *, double,
                    double,
                    wt::scatterer_type *>
             ((arg("dependee_0"), arg("a_0"),
               arg("dependee_1"), arg("a_1"),
               arg("b"),
               arg("scatterer"))))
          .def(init<af::shared<scalar_parameter *> const &,
                    af::shared<double> const &, double,
                    wt::scatterer_type *>
               ((arg("dependees"), arg("a"), arg("b"), arg("scatterer"))))
            ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

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
                    double,
                    bool,
                    wt::scatterer_type *>
               ((arg("occupancy"),
                 arg("original_multiplier"),
                 arg("multiplier"),
                 arg("as_one"),
                 arg("scatterer"))))
          .add_property("occupancy", make_function(&wt::reference, rir))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    void wrap_occupancy() {
      affine_asu_occupancy_parameter_wrapper::wrap();
      dependent_occu_wrapper::wrap();
    }


  }}}}
