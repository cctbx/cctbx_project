#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>


#include <mmtbx/scaling/relative_scaling.h>

namespace mmtbx { namespace scaling {


namespace{


  struct least_squares_on_i_wrapper
  {
    typedef scaling::relative_scaling::least_squares_on_i<> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("least_squares_on_i", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&
                 >
             ((
               arg_("hkl"),
               arg_("i_nat"),
               arg_("sig_nat"),
               arg_("i_der"),
               arg_("sig_nat"),
               arg_("p_scale"),
               arg_("unit_cell"),
               arg_("u_rwgk")
             )))
        .def("get_function",( double(w_t::*)() ) &w_t::get_function)
        .def("get_function",( double(w_t::*)(unsigned) ) &w_t::get_function)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)() ) &w_t::get_gradient)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)(unsigned) ) &w_t::get_gradient)
        .def("set_p_scale", &w_t::set_p_scale)
        .def("set_u_rwgk", &w_t::set_u_rwgk)
        .def("set_params", &w_t::set_params)
        ;
    }
  };

  struct least_squares_on_f_wrapper
  {
    typedef scaling::relative_scaling::least_squares_on_f<> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("least_squares_on_f", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&
                 >
             ((
               arg_("hkl"),
               arg_("f_nat"),
               arg_("sig_nat"),
               arg_("f_der"),
               arg_("sig_nat"),
               arg_("p_scale"),
               arg_("unit_cell"),
               arg_("u_rwgk")
             )))

        .def("get_function",( double(w_t::*)() ) &w_t::get_function)
        .def("get_function",( double(w_t::*)(unsigned) ) &w_t::get_function)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)() ) &w_t::get_gradient)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)(unsigned) ) &w_t::get_gradient)
        .def("set_p_scale", &w_t::set_p_scale)
        .def("set_u_rwgk", &w_t::set_u_rwgk)
        .def("set_params", &w_t::set_params)
        ;
    }
  };


  struct least_squares_on_i_wt_wrapper
  {
    typedef scaling::relative_scaling::least_squares_on_i_wt<> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("least_squares_on_i_wt", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&
                 >
             ((
               arg_("hkl"),
               arg_("i_nat"),
               arg_("sig_nat"),
               arg_("i_der"),
               arg_("sig_nat"),
               arg_("p_scale"),
               arg_("unit_cell"),
               arg_("u_rwgk")
             )))
        .def("get_function",( double(w_t::*)() ) &w_t::get_function)
        .def("get_function",( double(w_t::*)(unsigned) ) &w_t::get_function)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)() ) &w_t::get_gradient)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)(unsigned) ) &w_t::get_gradient)
        .def("set_p_scale", &w_t::set_p_scale)
        .def("set_u_rwgk", &w_t::set_u_rwgk)
        .def("set_params", &w_t::set_params)
        ;
    }
  };




  struct least_squares_on_f_wt_wrapper
  {
    typedef scaling::relative_scaling::least_squares_on_f_wt<> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("least_squares_on_f_wt", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   double const&,
                   cctbx::uctbx::unit_cell const&,
                   scitbx::sym_mat3<double> const&
                 >
             ((
               arg_("hkl"),
               arg_("f_nat"),
               arg_("sig_nat"),
               arg_("f_der"),
               arg_("sig_nat"),
               arg_("p_scale"),
               arg_("unit_cell"),
               arg_("u_rwgk")
             )))
        .def("get_function",( double(w_t::*)() ) &w_t::get_function)
        .def("get_function",( double(w_t::*)(unsigned) ) &w_t::get_function)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)() ) &w_t::get_gradient)
        .def("get_gradient",( scitbx::af::shared<double>(w_t::*)(unsigned) ) &w_t::get_gradient)
        .def("set_p_scale", &w_t::set_p_scale)
        .def("set_u_rwgk", &w_t::set_u_rwgk)
        .def("set_params", &w_t::set_params)
        ;
    }
  };
















  struct local_scaling_moment_based_wrapper
  {
    typedef scaling::relative_scaling::local_scaling_moment_based<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("local_scaling_moment_based", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< cctbx::miller::index<> > const&,

                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,

                   cctbx::sgtbx::space_group const&,
                   bool const&,
                   std::size_t const&,
                   std::size_t const&,
                   std::size_t const&,

                   bool const&
                 >
             ((
                arg_("hkl_master"),
                arg_("hkl_sets"),
                arg_("data_set_a"),
                arg_("sigma_set_a"),
                arg_("data_set_b"),
                arg_("sigma_set_b"),
                arg_("space_group"),
                arg_("anomalous_flag"),
                arg_("radius"),
                arg_("depth"),
                arg_("target_ref"),
                arg_("use_experimental_sigmas")
                )))
        .def("get_scales", &w_t::get_scales)
        ;

        }
  };


  struct local_scaling_ls_based_wrapper
  {
    typedef scaling::relative_scaling::local_scaling_ls_based<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("local_scaling_ls_based", no_init)
        .def(init< scitbx::af::const_ref< cctbx::miller::index<> > const&,
                   scitbx::af::const_ref< cctbx::miller::index<> > const&,

                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,
                   scitbx::af::const_ref< double > const&,

                   cctbx::sgtbx::space_group const&,
                   bool const&,
                   std::size_t const&,
                   std::size_t const&,
                   std::size_t const&,

                   bool const&
                 >
             ((
                arg_("hkl_master"),
                arg_("hkl_sets"),
                arg_("data_set_a"),
                arg_("sigma_set_a"),
                arg_("data_set_b"),
                arg_("sigma_set_b"),
                arg_("space_group"),
                arg_("anomalous_flag"),
                arg_("radius"),
                arg_("depth"),
                arg_("target_ref"),
                arg_("use_experimental_sigmas")
                )))
        .def("get_scales", &w_t::get_scales)
        ;

        }
  };















}  // namespace <anonymous>

namespace boost_python {

  void wrap_local_scaling_moment_based()
  {
    local_scaling_moment_based_wrapper::wrap();
  }

  void wrap_local_scaling_ls_based()
  {
    local_scaling_ls_based_wrapper::wrap();
  }

  void wrap_least_squares_on_i()
  {
    least_squares_on_i_wrapper::wrap();
  }

  void wrap_least_squares_on_f()
  {
    least_squares_on_f_wrapper::wrap();
  }

  void wrap_least_squares_on_i_wt()
  {
    least_squares_on_i_wt_wrapper::wrap();
  }

  void wrap_least_squares_on_f_wt()
  {
    least_squares_on_f_wt_wrapper::wrap();
  }

}}} //namespace mmtbx::scaling::relative_scaling
