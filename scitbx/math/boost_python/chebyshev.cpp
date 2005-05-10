#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/math/chebyshev.h>

namespace scitbx { namespace math {

namespace {


  struct chebyshev_base_wrappers
  {
    typedef chebyshev::chebyshev_base<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_base", no_init)
        .def(init<std::size_t const&,
                  double const&,
                  double const&>((arg_("n_terms"),
                                  arg_("low_limit"),
                                  arg_("high_limit") )))
        .def(init<std::size_t const& ,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&>
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit"),
               arg_("cheb_coefs") )))
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
        .def("coefs", &w_t::coefs )
        .def("replace", &w_t::replace )
        ;

    }


  };


  struct chebyshev_polynome_wrappers
  {
    typedef chebyshev::chebyshev_polynome<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_polynome", no_init)
        .def(init<std::size_t const& ,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&>
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit"),
               arg_("cheb_coefs") )))
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
        .def("coefs", &w_t::coefs )

        .def("dfdx", (double(w_t::*)(double const&)) &w_t::dfdx)
        .def("dfdx", (scitbx::af::shared<double>(w_t::*)
                      (scitbx::af::const_ref<double> const&))
             &w_t::dfdx )
        .def("dfdx_coefs", &w_t::dfdx_coefs)
        ;

    }


  };


  struct chebyshev_fitter_wrappers
  {
    typedef chebyshev::chebyshev_fitter<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_fitter", no_init)
        .def(init<std::size_t const& ,
                  double const&,
                  double const& >
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit") )))
        .def(init<std::size_t const& ,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&>
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit"),
               arg_("cheb_coefs") )))


        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
        .def("coefs", &w_t::coefs )
        .def("replace", &w_t::replace)
        .def("dfdcoefs", &w_t::dfdcoefs)
        ;

    }


  };





  struct chebyshev_smooth_wrappers
  {
    typedef chebyshev::chebyshev_smooth<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_smooth", no_init)
        .def(init<std::size_t const&,
                  double const&,
                  double const&>((arg_("n_terms"),
                                  arg_("low_limit"),
                                  arg_("high_limit") )))
        .def(init<std::size_t const& ,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&>
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit"),
               arg_("cheb_coefs") )))
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
        .def("coefs", &w_t::smooth_coefs )
        .def("replace", &w_t::replace_and_smooth )
        ;

    }


  };



  struct chebyshev_smooth_fitter_wrappers
  {
    typedef chebyshev::chebyshev_smooth_fitter<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_smooth_fitter", no_init)
        .def(init<std::size_t const&,
                  double const&,
                  double const&>((arg_("n_terms"),
                                  arg_("low_limit"),
                                  arg_("high_limit") )))
        .def(init<std::size_t const& ,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&>
             ((arg_("n_terms"),
               arg_("low_limit"),
               arg_("high_limit"),
               arg_("cheb_coefs") )))
        .def("f", (double(w_t::*)(double const&)) &w_t::f)
        .def("f", (scitbx::af::shared<double>(w_t::*)
                   (scitbx::af::const_ref<double> const&))
                    &w_t::f )
        .def("coefs", &w_t::smooth_coefs )
        .def("replace", &w_t::replace_and_smooth )
        .def("dfdcoefs", &w_t::dfdcoefs)
        ;

    }


  };




  struct chebyshev_lsq_wrappers
  {
    typedef chebyshev::chebyshev_lsq<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<w_t>("chebyshev_lsq", no_init)
        .def(init<std::size_t const&,
                  double const&,
                  double const&,
                  scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<double> const&,
                  scitbx::af::const_ref<bool> const&
                  >((arg_("n_terms"),
                     arg_("low_limit"),
                     arg_("high_limit"),
                     arg_("x_obs"),
                     arg_("y_obs"),
                     arg_("w_obs"),
                     arg_("free_flags") )))

        .def("residual", &w_t::residual)
        .def("free_residual", &w_t::free_residual)
        .def("gradient", &w_t::gradient )
        .def("replace", &w_t::replace )
        .def("coefs", &w_t::coefs)
        ;

    }


  };








} // namespace <anonymous>

namespace boost_python {


  void wrap_chebyshev_base()
  {
    chebyshev_base_wrappers::wrap();
  }
  void wrap_chebyshev_polynome()
  {
    chebyshev_polynome_wrappers::wrap();
  }
  void wrap_chebyshev_fitter()
  {
    chebyshev_fitter_wrappers::wrap();
  }

  void wrap_chebyshev_smooth()
  {
    chebyshev_smooth_wrappers::wrap();
  }
  void wrap_chebyshev_smooth_fitter()
  {
    chebyshev_smooth_fitter_wrappers::wrap();
  }

  void wrap_chebyshev_lsq()
  {
    chebyshev_lsq_wrappers::wrap();
  }



}}} // namespace scitbx::math::boost_python
