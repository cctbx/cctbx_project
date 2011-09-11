#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <scitbx/boost_python/container_conversions.h>

#include <smtbx/refinement/constraints/same_group.h>

namespace smtbx { namespace refinement { namespace constraints {
  namespace boost_python {

    struct same_group_xyz_wrapper  {
      typedef same_group_xyz wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("same_group_xyz", no_init)
          .def(init<af::shared<wt::scatterer_type *> const &,
                    af::shared<site_parameter *> const &,
                    scitbx::mat3<double> const &,
                    independent_small_vector_parameter<6> *>
               ((arg("scatterers"),
                 arg("sites"),
                 arg("alignment_matrix"),
                 arg("shifts_and_angles")
                 )))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct same_group_u_iso_wrapper  {
      typedef same_group_u_iso wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("same_group_u_iso", no_init)
          .def(init<af::shared<wt::scatterer_type *> const &,
                    af::shared<scalar_parameter *> const &,
                    independent_small_vector_parameter<6> *>
               ((arg("scatterers"),
                 arg("u_isos"),
                 arg("shifts_and_angles")
                 )))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };

    struct same_group_u_star_wrapper  {
      typedef same_group_u_star wt;

      static void wrap() {
        using namespace boost::python;
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >("same_group_u_star", no_init)
          .def(init<af::shared<wt::scatterer_type *> const &,
                    af::shared<u_star_parameter *> const &,
                    scitbx::mat3<double> const &,
                    independent_small_vector_parameter<6> *>
               ((arg("scatterers"),
                 arg("u_stars"),
                 arg("alignment_matrix"),
                 arg("shifts_and_angles")
                 )))
          ;
        implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
      }
    };


    void wrap_same_group() {
      {
        using namespace scitbx::boost_python::container_conversions;
        tuple_mapping_variable_capacity<
          af::shared<u_star_parameter *> >();
        tuple_mapping_variable_capacity<
          af::shared<scalar_parameter *> >();
      }
      same_group_xyz_wrapper::wrap();
      same_group_u_iso_wrapper::wrap();
      same_group_u_star_wrapper::wrap();
    }

}}}}
