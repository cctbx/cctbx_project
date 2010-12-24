#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <scitbx/boost_python/container_conversions.h>

#include <smtbx/refinement/constraints/rigid.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct pivoted_rotable_group_wrapper  {
      typedef pivoted_rotable_group wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("rigid_pivoted_rotable_group", no_init)
          .def(init<site_parameter *,
                    site_parameter *,
                    independent_scalar_parameter *,
                    const af::shared<wt::scatterer_type *>&>
               ((arg("pivot"),
                 arg("pivot_neighbour"),
                 arg("azimuth"),
                 arg("scatterers"))))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct rigid_site_proxy_wrapper {
      typedef rigid_site_proxy wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<site_parameter>,
               std::auto_ptr<wt> >("rigid_site_proxy", no_init)
          .def(init<pivoted_rotable_group *,
                    int>
                ((arg("parent"),
                  arg("index"))))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    void wrap_rigid() {
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping_variable_capacity<
          af::shared<asu_parameter::scatterer_type *> >();
      }
      pivoted_rotable_group_wrapper::wrap();
      rigid_site_proxy_wrapper::wrap();
    }


  }}}}
