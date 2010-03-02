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


        //.def("get_coef", &w_t::get_coef )
        //.def("load_coefs", &w_t::load_coefs)
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

  struct zernike_grid_wrapper
  {
    typedef zernike::zernike_grid<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("zernike_grid", no_init)
        .def(init<int const&, int const&>
                  ((  arg("m"),
                      arg("n_max")
                  )))
        .def("xyz", &w_t::xyz )
        //.def("get_coef", &w_t::get_coef )
        //.def("load_coefs", &w_t::load_coefs)
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
    log_factorial_generator_wrapper::wrap();
    zernike_grid_wrapper::wrap();
  }






}}} // namespace scitbx::math::boost_python
