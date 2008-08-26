#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/grouped_data.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/boost_python/is_polymorphic_workaround.h>


namespace cctbx { namespace xray { namespace grouped_data { namespace boost_python {
  namespace {

    struct unmerged_data_wrappers
    {
      typedef unmerged_data<double> w_t;
      static void
      wrap()
      {
        using namespace boost::python;

        class_<w_t>("unmerged_data", no_init)
          .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                     scitbx::af::const_ref< cctbx::miller::index<> > const&,
                     sgtbx::space_group const&,
                     bool const&
                   >
               (( arg_("hkl_obs"),
                  arg_("asu_hkl"),
                  arg_("space_group"),
                  arg_("anomalous_flag") )))
          ;
      }


    };




    struct merger_wrappers
    {
      typedef merger<double> w_t;
      static void
      wrap()
      {
        using namespace boost::python;

        class_<w_t>("merger", no_init)
          .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                     scitbx::af::const_ref< double >                 const&,
                     scitbx::af::const_ref< double >                 const&,
                     sgtbx::space_group                              const&,
                     bool                                            const&,
                     cctbx::uctbx::unit_cell                         const&
                   >
               (( arg_("hkl_obs"),
                  arg_("i_obs"),
                  arg_("s_obs"),
                  arg_("space_group"),
                  arg_("anomalous_flag"),
                  arg_("unit_cell")
               )))
          .def("bic", &w_t::bic)
          .def("r_abs", &w_t::r_abs)
          ;
      }


    };



  }}} // namespace cctbx::xray

  namespace boost_python {

    void wrap_grouped_data()
    {
      grouped_data::boost_python::unmerged_data_wrappers::wrap();
      grouped_data::boost_python::merger_wrappers::wrap();
    }
  }

}}
