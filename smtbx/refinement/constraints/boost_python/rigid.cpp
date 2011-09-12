#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <scitbx/boost_python/container_conversions.h>

#include <smtbx/refinement/constraints/rigid.h>
#include <smtbx/refinement/constraints/proxy.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct pivoted_rotatable_group_wrapper  {
      typedef pivoted_rotatable_group wt;

      static void wrap() {
        using namespace boost::python;
        return_internal_reference<> rir;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("rigid_pivoted_rotatable_group", no_init)
          .def(init<site_parameter *,
                    site_parameter *,
                    independent_scalar_parameter *,
                    independent_scalar_parameter *,
                    const af::shared<wt::scatterer_type *>&>
               ((arg("pivot"),
                 arg("pivot_neighbour"),
                 arg("azimuth"),
                 arg("size"),
                 arg("scatterers"))))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct rotatable_expandable_group_wrapper  {
      typedef rotatable_expandable_group wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("rigid_rotatable_expandable_group", no_init)
          .def(init<site_parameter *,
                    independent_scalar_parameter *,
                    independent_scalar_parameter *,
                    independent_scalar_parameter *,
                    independent_scalar_parameter *,
                    const af::shared<wt::scatterer_type *>&>
               ((arg("pivot"),
                 arg("size"),
                 arg("alpha"),
                 arg("beta"),
                 arg("gamma"),
                 arg("scatterers"))))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };


    struct riding_expandable_group_wrapper  {
      typedef riding_expandable_group wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("rigid_riding_expandable_group", no_init)
          .def(init<site_parameter *,
                    independent_scalar_parameter *,
                    const af::shared<wt::scatterer_type *>&>
               ((arg("pivot"),
                 arg("size"),
                 arg("scatterers"))))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct rigid_site_proxy_wrapper {
      typedef site_proxy<rigid_group_base> wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<site_parameter>,
               std::auto_ptr<wt> >("rigid_site_proxy", no_init)
          .def(init<pivoted_rotatable_group *,
                    int>
                ((arg("parent"),
                  arg("index"))))
          .def(init<riding_expandable_group *,
                    int>
                ((arg("parent"),
                  arg("index"))))
          .def(init<rotatable_expandable_group *,
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
      pivoted_rotatable_group_wrapper::wrap();
      rotatable_expandable_group_wrapper::wrap();
      riding_expandable_group_wrapper::wrap();
      rigid_site_proxy_wrapper::wrap();
    }


  }}}}
