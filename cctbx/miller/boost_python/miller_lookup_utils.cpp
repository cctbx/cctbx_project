#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <cctbx/miller/miller_lookup_utils.h>

namespace cctbx { namespace miller {



namespace{

  struct miller_index_lookup_tensor_wrapper
  {
    typedef miller_lookup_utils::miller_index_lookup_tensor<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("miller_index_lookup_tensor", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> >  const&,
                   cctbx::sgtbx::space_group const&,
                   bool const&
                 >((arg_("miller_indices"),
                    arg_("space_group"),
                    arg_("anomalous_flag") )))
        .def("find_hkl", (long(w_t::*)(cctbx::miller::index<> const&))
             &w_t::find_hkl)
        .def("find_hkl", (scitbx::af::shared< long >(w_t::*)
                          (scitbx::af::const_ref< cctbx::miller::index<> >
                           const&))
             &w_t::find_hkl)
        ;
    }

  };



    struct local_neighbourhood_wrapper
  {
    typedef miller_lookup_utils::local_neighbourhood<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("local_neighbourhood", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   cctbx::sgtbx::space_group const&,
             bool const&,
                   long const&
                 >
             ((
                arg_("miller_indices"),
                arg_("space_group"),
                arg_("anomalous_flag"),
                arg_("radius")
                )))
        .def("construct_neighbourhood",
             ( std::vector<unsigned>(w_t::*) ( cctbx::miller::index<> const&) )
             &w_t::construct_neighbourhood)
        .def("construct_neighbourhood",
             ( scitbx::af::shared< std::vector<unsigned> >(w_t::*)
               ( scitbx::af::shared<cctbx::miller::index<> > const&) )
             &w_t::construct_neighbourhood)


        .def("construct_neighbourhood",
             ( std::vector<unsigned>(w_t::*) ( unsigned const&) )
             &w_t::construct_neighbourhood)
        .def("construct_neighbourhood",
             ( scitbx::af::shared< std::vector<unsigned> >(w_t::*)
               ( scitbx::af::shared< long > const&) )
             &w_t::construct_neighbourhood)
        .def("construct_neighbourhood",
             ( scitbx::af::shared< std::vector<unsigned> >(w_t::*)
               ( ) )
             &w_t::construct_neighbourhood)


        ;

        }
  };


  struct local_area_with_property_wrapper
  {
    typedef miller_lookup_utils::local_area_with_property<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("local_area_with_property", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< bool > const&,
                   cctbx::sgtbx::space_group const&,
                   bool const&,
                   std::size_t const&,
                   std::size_t const&,
                   std::size_t const&
                 >
             ((
                arg_("miller_indices"),
                arg_("property"),
                arg_("space_group"),
                arg_("anomalous_flag"),
                arg_("radius"),
                arg_("depth"),
                arg_("target_ref")
                )))
        .def("area", &w_t::area)
        ;

        }
  };











}  // namespace <anonymous>

namespace boost_python {

  void wrap_miller_index_lookup_tensor()
  {
    miller_index_lookup_tensor_wrapper::wrap();
  }

  void wrap_local_neighbourhood()
  {
    local_neighbourhood_wrapper::wrap();
  }

  void wrap_local_area_with_property()
  {
    local_area_with_property_wrapper::wrap();
  }

}}}
