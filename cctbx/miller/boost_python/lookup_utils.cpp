#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <cctbx/miller/lookup_utils.h>

namespace cctbx { namespace miller {



namespace{

  struct lookup_tensor_wrapper
  {
    typedef lookup_utils::lookup_tensor<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("lookup_tensor", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> >  const&,
                   cctbx::sgtbx::space_group const&,
                   bool const&
                 >((arg("miller_indices"),
                    arg("space_group"),
                    arg("anomalous_flag") )))
        .def("n_duplicates", &w_t::n_duplicates)
        .def("find_hkl", (long(w_t::*)(cctbx::miller::index<> const&) const)
             &w_t::find_hkl)
        .def("find_hkl", (scitbx::af::shared< long >(w_t::*)
                          (scitbx::af::const_ref< cctbx::miller::index<> >
                           const&) const)
             &w_t::find_hkl)
        ;
    }

  };



  struct local_neighbourhood_wrapper
  {
    typedef lookup_utils::local_neighbourhood<> w_t;

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
                arg("miller_indices"),
                arg("space_group"),
                arg("anomalous_flag"),
                arg("radius")
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


  struct local_area
  {
    typedef lookup_utils::local_area<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("local_area", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< bool > const&,
                   cctbx::sgtbx::space_group const&,
                   bool const&,
                   std::size_t const&,
                   std::size_t const&,
                   std::size_t const&
                 >
             ((
                arg("miller_indices"),
                arg("property"),
                arg("space_group"),
                arg("anomalous_flag"),
                arg("radius"),
                arg("depth"),
                arg("target_ref")
                )))
        .def("area", &w_t::area)
        ;

        }
  };











}  // namespace <anonymous>

namespace boost_python {

  void wrap_lookup_tensor()
  {
    lookup_tensor_wrapper::wrap();
  }

  void wrap_local_neighbourhood()
  {
    local_neighbourhood_wrapper::wrap();
  }

  void wrap_local_area()
  {
    local_area::wrap();
  }

}}}
