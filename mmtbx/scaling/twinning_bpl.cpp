#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <mmtbx/scaling/twinning.h>

namespace mmtbx { namespace scaling {



namespace{


  struct twin_r_wrapper
  {
    typedef twinning::twin_r<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("twin_r", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const& >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("symop") )))
        .def("r_abs_value", &w_t::r_abs_value)
        .def("r_sq_value", &w_t::r_sq_value)
        ;
    }

  };




  struct l_test_wrapper
  {
    typedef twinning::l_test<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("l_test", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             long const&,
             long const&,
             long const&,
             std::size_t const& >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("parity_h"),
               arg_("parity_k"),
               arg_("parity_l"),
               arg_("max_delta_h") )))
        .def("mean_l", &w_t::mean_l)
        .def("mean_l2", &w_t::mean_l2)
        .def("ml_alpha", &w_t::ml_alpha)
        .def("cumul", &w_t::cumul)
        ;
    }

  };


  struct detwin_wrapper
  {
    typedef twinning::detwin<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("detwin", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&
             >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("sigma"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("twin_law") )))
        .def("detwin_with_alpha", &w_t::detwin_with_alpha)
        .def("detwinned_i", &w_t::detwinned_i)
        .def("detwinned_sigi", &w_t::detwinned_sigi)
        .def("completeness", &w_t::completeness)
        .def("location", &w_t::location)
        ;
    }

  };


  struct h_test_wrapper
  {
    typedef twinning::h_test<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("h_test", no_init)
        .def(init<
             scitbx::af::const_ref< cctbx::miller::index<> >  const&,
             scitbx::af::const_ref< double > const&,
             scitbx::af::const_ref< double > const&,
             cctbx::sgtbx::space_group const&,
             bool const&,
             scitbx::mat3<double> const&,
             double const&
             >
             ((arg_("miller_indices"),
               arg_("intensity"),
               arg_("sigma"),
               arg_("space_group"),
               arg_("anomalous_flag"),
               arg_("twin_law"),
               arg_("fraction")) ))
        .def("h_array", &w_t::h_array)
        .def("h_values", &w_t::h_values)
        .def("h_cumul_obs", &w_t::h_cumul_obs)
        .def("h_cumul_fit", &w_t::h_cumul_fit)
        .def("distance", &w_t::distance)
        .def("alpha", &w_t::alpha)
        .def("mean_h", &w_t::mean_h)
        .def("mean_h2", &w_t::mean_h2)
        ;
    }

  };







}  // namespace <anonymous>

namespace boost_python {

  void wrap_twin_r()
  {
    twin_r_wrapper::wrap();
  }

  void wrap_l_test()
  {
    l_test_wrapper::wrap();
  }

  void wrap_detwin()
  {
    detwin_wrapper::wrap();
  }


  void wrap_h_test()
  {
    h_test_wrapper::wrap();
  }


}}}
