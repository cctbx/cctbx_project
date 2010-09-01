#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/math/zernike.h>

namespace scitbx { namespace math {

namespace {


  struct nlm_array_wrapper
  {
    typedef zernike::nlm_array<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("nlm_array", no_init)
        .def(init<int const&>
                  ((arg("n_max") )))
        .def("set_coef", &w_t::set_coef )
        .def("get_coef", &w_t::get_coef )
        .def("load_coefs", &w_t::load_coefs)
        .def("select_on_nl", &w_t::select_on_nl)
        .def("nlm", &w_t::nlm)
        .def("coefs", &w_t::coefs)
       ;
    }
  };

  struct nl_array_wrapper
  {
    typedef zernike::nl_array<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("nl_array", no_init)
        .def(init<int const&>
                  ((arg("n_max") )))
        .def("set_coef", &w_t::set_coef )
        .def("get_coef", &w_t::get_coef )
        .def("load_coefs", &w_t::load_coefs)
        .def("nl", &w_t::nl)
        .def("coefs", &w_t::coefs)
       ;
    }
  };


  struct log_factorial_generator_wrapper
  {
    typedef zernike::log_factorial_generator<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("log_factorial_generator", no_init)
        .def(init<int const&>
                  ((arg("n_max") )))
        .def("log_fact", &w_t::log_fact )
        .def("fact", &w_t::fact )
       ;
    }
  };



  struct nss_spherical_harmonics_wrapper
  {
    typedef zernike::nss_spherical_harmonics<> w_t;
    static void
    wrap()
    {
      using namespace boost::python;
      // int const& max_l, int const& mangle, log_factorial_generator<FloatType> const& lgf
      class_<w_t>("nss_spherical_harmonics", no_init)
        .def(init<int const&, int const&, zernike::log_factorial_generator<double> const&>
                  ((arg("l_max"), arg("mangle"), arg("lgf") )))
        .def("legendre_lm", &w_t::legendre_lm)
        .def("legendre_lm_pc", &w_t::legendre_lm_pc)
        .def("spherical_harmonic_pc", &w_t::spherical_harmonic_pc)
        .def("spherical_harmonic", &w_t::spherical_harmonic_direct)
       ;
    }
  };




  struct zernike_radial_wrapper
  {
    typedef zernike::zernike_radial<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_radial", no_init)
        .def(init<int const& , int const&, zernike::log_factorial_generator<double> const&>
                  ((arg("n"),
                    arg("l"),
                    arg("log_factorial_generator")
                   )) )
        .def("Nnlk", &w_t::Nnlk)
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
       ;
    }
  };



  struct zernike_2d_radial_wrapper
  {
    typedef zernike::zernike_2d_radial<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_2d_radial", no_init)
        .def(init<int const& , int const&, zernike::log_factorial_generator<double> const&>
                  ((arg("n"),
                    arg("l"),
                    arg("log_factorial_generator")
                   )) )
        .def("Nnlk", &w_t::Nnlk)
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
       ;
    }
  };









  struct zernike_polynome_wrapper
  {
    typedef zernike::zernike_polynome<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_polynome", no_init)
        .def(init<int const&, int const&, int const&, zernike::zernike_radial<double> const& >
                  ((arg("n"),
                    arg("l"),
                    arg("m"),
                    arg("Rnl")
                   )))
         .def("f", &w_t::f)
       ;
    }
  };

  struct zernike_2d_polynome_wrapper
  {
    typedef zernike::zernike_2d_polynome<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_2d_polynome", no_init)
        .def(init<int const&, int const&, zernike::zernike_2d_radial<double> const& >
                  ((arg("n"),
                    arg("l"),
                    arg("Rnl")
                   )))
         .def("f", &w_t::f)
       ;
    }
  };


  struct zernike_grid_wrapper
  {
    typedef zernike::zernike_grid<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_grid", no_init)
        .def(init<int const&, int const&, bool const&>
                  ((  arg("m"),
                      arg("n_max"),
                      arg("hex")
                  )))
        .def("xyz", &w_t::xyz )
        .def("rtp", &w_t::rtp )
        //.def("get_coef", &w_t::get_coef )
        .def("load_coefs", &w_t::load_coefs)
        .def("f", &w_t::f)
        .def("f_real", &w_t::f_real)
        .def("nlm", &w_t::nlm)
        .def("coefs", &w_t::coefs)
        .def("slow_moments", &w_t::slow_moments)
       ;
    }
  };



} // namespace <anonymous>

namespace boost_python {


  void wrap_zernike()
  {
    nlm_array_wrapper::wrap();
    nl_array_wrapper::wrap();
    zernike_polynome_wrapper::wrap();
    zernike_radial_wrapper::wrap();
    zernike_2d_radial_wrapper::wrap();
    zernike_2d_polynome_wrapper::wrap();
    log_factorial_generator_wrapper::wrap();
    zernike_grid_wrapper::wrap();
    nss_spherical_harmonics_wrapper::wrap();
  }






}}} // namespace scitbx::math::boost_python
